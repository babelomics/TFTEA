#####################################################################################
#######################    Pan-cancer survival analisys    ########################## 
#####################################################################################
#By: Matías Marín Falco                                                             #
#                                                                                   #
#The required data format for the analysis is explained in                          #
#the README file in my Git Hub repo: https://github.com/1mati1/pan-cancer-analysis/ #
#or in the babelomics repo:https://github.com/babelomics/TTFEA                      #
#                                                                                   #
#The methods are also described in the paper:                                       #
#The pan-cancer pathological regulatory landscape                                   #
#####################################################################################

rm (list = ls ())# clean the global environment

#define cancers to use in the analysis
cancers<-c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA","UCEC")

#####perform differential expression analysis
for (cancer in cancers){
  normales<- read.table(paste0(cancer,"/adjacent.tsv"), header=T, sep="\t",stringsAsFactors=F)#get normal samples
  tumors<- read.table(paste0(cancer,"/tumor.tsv"), header=T, sep="\t",stringsAsFactors=F)#get tumor samples
  genes<- read.table(paste0(cancer,"/genes.txt"), header=F, sep="\t",stringsAsFactors=F)$V1#get gene names
  
  ##remove duplicated donors
  normals<- normals[!duplicated(normals[ , "icgc_donor_id"]),]
  tumors<-tumors[!duplicated(tumors[ , "icgc_donor_id"]),]
  
  ##remove "gene_id" identifier in genes.txt
  genes<-genes[genes!="gene_id"]
  
  ##create a matrix expression for each cancer
  rawexp<-mat.or.vec(nr=length(c(normals$icgc_donor_id,tumors$icgc_donor_id)),nc=length(genes) )
  rownames(rawexp) <- c(normals$icgc_specimen_id,tumors$icgc_specimen_id)
  colnames(rawexp) <- genes
  
  ##read expression file
  library(ff)#as R loads in memory all object it may be hard for some computers to load the expression file, so this package helps to patch this file and read it
  x<-dir(path=cancer)[grep("raw_exp_seq.",dir(path=cancer))]#find the expression file
  expresion<-read.table.ffdf(file=paste0(cancer,"/",x), header=T, sep="\t")
  
  ##create a column in the expression file indicating for each expression data (row) if it is a tumour or normal specimen (in order to have a faster and easier acces to data later)
  specimen_type <- c( rep("TUMOR",nrow(tumors)), rep("NORMAL",nrow(normals)) )
  names(specimen_type) <- c( tumors$icgc_specimen_id, normals$icgc_specimen_id)
  
  library("ffbase")
  library("ETLUtils")
  expresion$sp_type <- with(expresion[c('icgc_specimen_id')], 
                            recoder(as.character(icgc_specimen_id), from = names(specimen_type), to = specimen_type))
  
  ##fill the expression matrix(this may take a while)
  library("pbapply") #use pbapply to observe function progress
  for(i in chunk(expression)){
    y <- expression[i, ]
    pbapply(y,1,function(x){
      
      if( x["sp_type"]=="TUMOR"){
        rawexp[x["icgc_specimen_id"],x["gene_id"]] <<- as.numeric(x["raw_read_count"])
      }
      if( x["sp_type"]=="NORMAL"){
        rawexp[x["icgc_specimen_id"],x["gene_id"]] <<- as.numeric(x["raw_read_count"])
        
      }
    })
  }
  
  ## remove empty rows
  rawexp<- data.matrix(rawexp[which(apply(rawexp,1,function(x){sum(as.numeric(x))>0})),])
  ##transpose matrix(required by normalization function)
  rawexp<- t(rawexp)
  
  ###TMM normalization
  vec<-specimen_type[colnames(rawexp)]#make a vector indicating specimen(colnames) condition
  library(edgeR)#package for normalization
  vec<-specimen_type[colnames(rawexp)]
  tmmexp<- DGEList(counts=as.matrix(rawexp), group=vec)
  tmmexp<- calcNormFactors(tmmexp)
  
  ###LIMMA diferential expression
  library(limma)
  vec<- factor(vec,levels=c("NORMAL","TUMOR"))
  design<- model.matrix(~vec)
  limma <- voom(tmmexp,design,plot=TRUE)
  fit <- lmFit(limma,design)
  fit <- eBayes(fit)
  tts <- topTable(fit, coef=ncol(design),number=nrow(rawexp))
  ##order by estadistic t value
  tts <- tts[order(tts$t,decreasing=T),]
  ##write ranked gene list that will be used later
  save.image(file=paste0(cancer,"/",stage,".RData"))
  t<-as.data.frame(tts[,'t'], row.names=rownames(tts))
  save.image(file=paste0(cancer,"/survival_donors.RData"))
  print (cancer)
  rm (list = setdiff(ls(),"cancers"))
}

#####extract individual gene expression values for each donor
rm (list = ls ()) # clean the global environment

#define cancers to use in the analysis
cancers<-c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA","UCEC")

###get clinical data and merge it
tcga<-list()
icgc<-list()
exp_table<-list()
normals<-list()
tumors<-list()
library(httr)
for (cancer in cancers){
  #get TGCA clinical data
  x<-content(GET(paste("https://tcga-data.nci.nih.gov/tcgafiles/TFp_auth/distro_TFpusers/anonymous/tumor/",tolower(cancer),"/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_",tolower(cancer),".txt",sep="")))
  tcga[[cancer]]<-read.table(text=x, sep="\t",header=T,quote="",stringsAsFactors=F)
  rownames(tcga[[cancer]])<-tcga[[cancer]]$bcr_patient_barcode
  #get ICGC clinical data
  x<-dir(path=cancer)[grep("clinical",dir(path=cancer))]
  icgc[[cancer]]<- read.table(paste(cancer,"/",x,sep=""), sep="\t",header=T,quote="",stringsAsFactors=F)
  #get normalized expression values
  load(paste0(cancer,"/survival_donors.RData"))
  exp_table[[cancer]]<-2^limma$E      
  normals[[cancer]]<- exp_table[[cancer]][,names(vec[vec=="NORMAL"])]
  tumors[[cancer]]<-exp_table[[cancer]][,names(vec[vec=="TUMOR"])]
  cat(c(cancer,"\n normals:\t",ncol(normals[[cancer]]),"\n tumors: \t", ncol(tumors[[cancer]]),"\n"))
  
}
gene_names<- rownames( normals[[cancer]])

rm(list=setdiff(ls(),c("tcga","icgc","normals","tumors","gene_names")))
cancers<-c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA","UCEC")

### create lists containing normal mean and standard deviation expression values for each gene in each cancer 
mean_list<-list()
sd_list<-list()
for (cancer in cancers){
  mean_list[[cancer]] <- list()
  sd_list[[cancer]] <- list()
  for(gen in gene_names){
    mean_list[[cancer]] [[gen]] <- mean(normals[[cancer]][gen,])
    sd_list[[cancer]] [[gen]] <- sd(normals[[cancer]][gen,])
    ##obtain the individual differential expresion value for each gene comparing the expression value with the mean: (exp-value-mean)/std.dev
    tumors[[cancer]][gen,]<-  (tumors[[cancer]][gen,] - as.numeric(mean_list[[cancer]][[gen]]))/as.numeric(sd_list[[cancer]][[gen]])
  }
  print(cancer)
}

#####calculate individual TTFEA values

##load available TFs and modifie necessary names
TF<-as.matrix (read.table ("anotation_file.tsv", header = FALSE, sep = "\t", quote = "", colClasses = "character"))
ntargets<-as.character(table(TF[,2]))
names(ntargets)<-toupper(names(table(TF[,2])))
TF<-toupper(unique(TF[,2]))
olds<-c("CMYC","CFOS","PU1","GABP","AP2ALPHA","NRSF","AP2GAMMA","CJUN","EBF","NFKB","P300","BRG1","INI1","BAF155","TR4","BAF170")
news<-c("MYC","FOS","SPI1","GABPA?GABPB1?2?","TFAP2A","REST","TFAP2C","JUN","EBF1","NFKB1?NFKB2","EP300","SMARCA4","SMARCB1","SMARCC1","NR2C2","SMARCC2")
for(i in 1:length(olds)){
  oldgene<-olds[i]
  newgene<-news[i]
  TF[grep(oldgene, TF)]=newgene
}
TF<-TF[-grep("\\?",TF)]

rm(list=setdiff(ls(),c("tcga","icgc","normals","tumors","cancers","gene_names","TF","mean_list", "sd_list","ntargets")))

################################################################################
#Transcription Factor Target Gene Analysis (TTFEA) FUNCTION

uvGsa <- function (rankstat, annotation, p.adjust.method = "BH") {
  
  t0 <- proc.time ()
  
  ##variable names (colnames)
  statsnames <- colnames (rankstat)[1]
  if (is.null (statsnames)) statsnames <- "X"
  
  ##transform vector, matrix or data.frame keepsing rownames
  rankstat <- as.matrix (rankstat)[,1] 
  
  
  ##LOOP over GOs
  GOs <- unique (annotation[,2])
  ncolRES <- 7 #3 + 4
  RES <- matrix (nrow = 0, ncol = ncolRES)
  separador <- "." 
  
  colnames (RES) <- c ("size", "conv", "error",                     #COLS: 1:3
                       paste ("LOR", statsnames, sep = separador),  #COLS: 4
                       paste ("sd",  statsnames, sep = separador),  #COLS: 5
                       paste ("z",   statsnames, sep = separador),  #COLS: 6
                       paste ("p",   statsnames, sep = separador))  #COLS: 7
  
  for (go in GOs) {
    
    geneswithgo <- annotation[annotation[,2] == go, 1]
    GObinario   <- as.numeric (names (rankstat) %in% geneswithgo)
    blocksize   <- sum (GObinario)
    
    G  <- try (glm (GObinario ~ rankstat, family = binomial))
    SG <- try (summary (G))
    
    try.error <- "try-error" %in% class (G) | "try-error" %in% class (SG)
    
    if (try.error) {
      RES <- rbind (RES, c(blocksize, NA, try.error, rep (NA, times = ncolRES - 3)))
    } else {
      RES <- rbind (RES, c(blocksize, G$converged, try.error, as.vector (SG$coefficients[-1,])))
    }
  }
  
  rownames (RES) <- GOs
  
  ##adjust p-values
  if (!is.null (p.adjust.method)) {
    adj.p.matrix <- as.matrix (p.adjust (RES[,7], method = p.adjust.method))
    colnames (adj.p.matrix) <- paste ("adj", statsnames, sep = separador)  #COLS: 8
    RES <- cbind (RES, adj.p.matrix)
  }
  
  RES <- as.data.frame (RES, stringsAsFactors = FALSE)
  
  t1 <- proc.time ()
  print ("timing in seconds")
  print (t1-t0)
  
  ##results
  return (RES)
}
################################################################################

#load annotation file
anot <- as.matrix (read.table ("anotation_file.tsv ", header = FALSE, sep = "\t", quote = "", colClasses = "character"))


##perform TTFEA for each cancer
TTFEA<-list()
for (cancer in cancers){
  TTFEA[[cancer]]<-list()
  for(muestra in colnames(tumors[[cancer]])){
    stat <-sort(tumors[[cancer]][,muestra],decreasing = T)
    res <- uvGsa (rankstat = stat, annotation = anot)
    res[,"converged"] <- res[,"conv"] - res[,"error"]            
    res[1:10,]
    table (res[,"conv"], res[,"error"])
    TTFEA[[cancer]][[muestra]] <-res
    
  }
  print(cancer)
}
save.image("sd_mean_survival.RData")

#####prepare survival data
merged_list<- list()
tumor_surv<-list()
p<-list()
for(cancer in cancers){
  ##get pvalues and z score values from TTFEA
  merged_list[[cancer]]<-as.data.frame(mat.or.vec(nr=length(names(TTFEA[[cancer]])),nc=length( rownames(TTFEA[[cancer]][[1]]))),stringsAsFactors = F)
  rownames(merged_list[[cancer]]) <-  names(TTFEA[[cancer]])
  colnames(merged_list[[cancer]]) <- rownames(TTFEA[[cancer]][[1]])
  for(specimen in names(TTFEA[[cancer]])){
    merged_list[[cancer]][specimen,]<-TTFEA[[cancer]][[specimen]]$z.X
  }
  p[[cancer]]<-as.data.frame(mat.or.vec(nr=length(names(TTFEA[[cancer]])),nc=length( rownames(TTFEA[[cancer]][[1]]))),stringsAsFactors = F)
  rownames(p[[cancer]]) <-  names(TTFEA[[cancer]])
  colnames(p[[cancer]]) <- rownames(TTFEA[[cancer]][[1]])
  for(specimen in names(TTFEA[[cancer]])){
    p[[cancer]][specimen,]<-as.numeric(TTFEA[[cancer]][[specimen]]$p.X)
  }
  ##correct TF names
  colnames(merged_list[[cancer]])<-toupper(colnames(merged_list[[cancer]]))
  colnames(p[[cancer]])<-toupper(colnames(p[[cancer]]))
  olds<-c("CMYC","CFOS","PU1","GABP","AP2ALPHA","NRSF","AP2GAMMA","CJUN","EBF","NFKB","P300","BRG1","INI1","BAF155","TR4","BAF170")
  news<-c("MYC","FOS","SPI1","GABPA?GABPB1?2?","TFAP2A","REST","TFAP2C","JUN","EBF1","NFKB1?NFKB2","EP300","SMARCA4","SMARCB1","SMARCC1","NR2C2","SMARCC2")
  for(i in 1:length(olds)){
    oldgene<-olds[i]
    newgene<-news[i]
    colnames(merged_list[[cancer]])[grep(oldgene, colnames(merged_list[[cancer]]))]=newgene
    colnames(p[[cancer]])[grep(oldgene, colnames(p[[cancer]]))]=newgene
  }
  colnames(p[[cancer]])<-paste0(colnames(p[[cancer]]),".p")
  
  #remove "GABPA?GABPB1?2?","NFKB1?NFKB2" from lists
  merged_list[[cancer]]<-merged_list[[cancer]][,-grep("\\?",colnames(merged_list[[cancer]]))]
  p[[cancer]]<-p[[cancer]][,-grep("\\?",colnames(p[[cancer]]))]
  
  #keep TF names
  available_TF<-c(colnames(p[[cancer]]))
  available_TF<-gsub("::","_",available_TF)#available TF tiene los TF y acaba en ".p"
  
  ###merge clinical data, zscore and pvalue on a same list 
  merged_list[[cancer]][,"specimen"]<-rownames(merged_list[[cancer]])
  merged_list[[cancer]]<- merge(p[[cancer]],merged_list[[cancer]], by.x= "row.names",by.y="row.names", all.y=T)
  
  #add icgc info
  merged_list[[cancer]]<- merge(icgc[[cancer]],merged_list[[cancer]], by.x= "icgc_specimen_id",by.y="specimen", all.y=T)
  
  #add tcga info
  merged_list[[cancer]]<-merge(tcga[[cancer]],merged_list[[cancer]], by.x= "bcr_patient_barcode",by.y="submitted_donor_id", all.y=T)
  
  ###keep to-be used columns and prepare data format
  #create days column
  merged_list[[cancer]]<-as.data.frame(append(merged_list[[cancer]], list(days= NA), aTFer = 1))
  merged_list[[cancer]]$last_contact_days_to<-gsub("[Not Available]",as.numeric(0),merged_list[[cancer]]$last_contact_days_to,fixed = T)
  merged_list[[cancer]]$death_days_to<-gsub("[Not Applicable]",as.numeric(0),merged_list[[cancer]]$death_days_to,fixed=T)
  merged_list[[cancer]]$death_days_to<-gsub("[Not Available]",as.numeric(0),merged_list[[cancer]]$death_days_to,fixed=T)
  merged_list[[cancer]]$days<-as.numeric(merged_list[[cancer]]$last_contact_days_to )+ as.numeric(merged_list[[cancer]]$death_days_to)
  
  #pass status to numeric value 
  merged_list[[cancer]]$vital_status<-gsub("Alive",as.numeric(0),merged_list[[cancer]]$vital_status,fixed=T)
  merged_list[[cancer]]$vital_status<-gsub("Dead",as.numeric(1),merged_list[[cancer]]$vital_status,fixed=T)  
  merged_list[[cancer]]$tumor_status<-gsub("[Not Available]",NA,merged_list[[cancer]]$tumor_status,fixed=T)
  merged_list[[cancer]]$tumor_status<-gsub("[Unknown]",NA,merged_list[[cancer]]$tumor_status,fixed=T)
  
  #merge changes some names, restore them
  colnames(merged_list[[cancer]])<-gsub("..","_",colnames(merged_list[[cancer]]),fixed = T)
  
  #keep to-be-used columns
  columns<-c("icgc_specimen_id","days","vital_status","gender","bcr_patient_barcode","donor_age_at_diagnosis","donor_age_at_last_followup","tumor_status","disease_status_last_followup","icgc_donor_id", "specimen_type", "submitted_donor_id",available_TF, gsub(".p","",available_TF,fixed = T) )
  columns<-columns[columns%in%colnames(merged_list[[cancer]])]
  merged_list[[cancer]]<-subset(merged_list[[cancer]],select= columns)    
  
  ###make a list only with tumor samples
  tumor_surv[[cancer]]<-merged_list[[cancer]][merged_list[[cancer]]$specimen_type=="Primary tumour - solid tissue",]
  rownames(tumor_surv[[cancer]])<-tumor_surv[[cancer]]$icgc_specimen_id
  
  #remove discrepancies from vital status
  if( "[Discrepancy]" %in% tumor_surv[[cancer]]$vital_status ){
    tumor_surv[[cancer]] <- tumor_surv[[cancer]][-which(tumor_surv[[cancer]]$vital_status=="[Discrepancy]"),]
  }
  
}


#####perform survival analysis
library(survival)


#### first make cox analysis

##create argumento for cox funciton
x=1
for (factor in gsub(".p","",available_TF,fixed = T)){
  if (x!=0 & factor %in% colnames(tumor_surv[[cancer]])){
    arg<-paste0("survival ~ tumor_surv[[cancer]]$",factor)
    x=0
  }
  else if(x==0 &  factor %in% colnames(tumor_surv[[cancer]])) {
    arg<-paste(arg,paste0("tumor_surv[[cancer]]$",factor),sep=" + ")
  }
}


#create a cox model for each cancer
p_cox_list<-list()
for(cancer in cancers){
  survival<- Surv(tumor_surv[[cancer]]$days, as.numeric(tumor_surv[[cancer]]$vital_status))
  p_cox_list[[cancer]] <- coxph (as.formula(arg))
  print(cancer)
}

#improve cox model using stepwise algorithm
ss<-list()
for(cancer in cancers){
  print (cancer)
  survival<- Surv(tumor_surv[[cancer]]$days, as.numeric(tumor_surv[[cancer]]$vital_status))
  tryCatch({
    ss[[cancer]] <- step(p_cox_list[[cancer]])
  }, error = function(e) {
    print("super error")
    print(e)
  })
  
}

####Kaplan-Meier(KM) analysis

library(ggplot2)
library(gridBase)
library(grid)

##create a dataframe where KM pvalues will be added, with a row per TF and cancer
KM_sig<-c()
for(cancer in canceres){
  KM_sig<-rbind(KM_sig,data.frame(rep(cancer,length(ntargets)),names(ntargets)))     
}
KM_sig<-cbind(KM_sig,rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)),rep(0,nrow(KM_sig)))
colnames(KM_sig)<-c("Cancer","Transcription_Factor","UADA_Pvalue" ,"UD_Pvalue","adj.UADA_Pvalue" ,"adj.UD_Pvalue","UAD","DAD","DD","UD","ND")
contingency_table<-list()
##fill table
for(line in 1:nrow(KM_sig)){ 
  x<-KM_sig[line,]
  #make also a contingency table
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]<-mat.or.vec(5,2)
  colnames(contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]])<-c("Alive","Deceased")
  rownames(contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]])<-c("Up_Altered","Down_Altered","Up","Down","Not_Altered")
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up_Altered","Alive"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]<0.05 & tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]>0 & tumor_surv[[x$Cancer]][["vital_status"]]==0))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up_Altered","Deceased"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]<0.05 & tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]>0 & tumor_surv[[x$Cancer]][["vital_status"]]==1))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down_Altered","Alive"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]<0.05 & tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]<0 & tumor_surv[[x$Cancer]][["vital_status"]]==0))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down_Altered","Deceased"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]<0.05 & tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]<0 & tumor_surv[[x$Cancer]][["vital_status"]]==1))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up","Alive"]<-length(which(tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]>0 & tumor_surv[[x$Cancer]][["vital_status"]]==0))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up","Deceased"]<-length(which(tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]>0 & tumor_surv[[x$Cancer]][["vital_status"]]==1))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down","Alive"]<-length(which(tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]<0 & tumor_surv[[x$Cancer]][["vital_status"]]==0))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down","Deceased"]<-length(which(tumor_surv[[x$Cancer]][[as.character(x$Transcription_Factor)]]<0 & tumor_surv[[x$Cancer]][["vital_status"]]==1))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Not_Altered","Alive"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]>0.05 & tumor_surv[[x$Cancer]][["vital_status"]]==0))
  contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Not_Altered","Deceased"]<-length(which(tumor_surv[[x$Cancer]][[paste0(x$Transcription_Factor,".p")]]>0.05 & tumor_surv[[x$Cancer]][["vital_status"]]==1))
  
  KM_sig[line,"ND"]<-contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Not_Altered","Deceased"]
  KM_sig[line,"DAD"]<-contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down_Altered","Deceased"]
  KM_sig[line,"UAD"]<-contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up_Altered","Deceased"]
  KM_sig[line,"DD"]<-contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Down","Deceased"]
  KM_sig[line,"UD"]<-contingency_table[[paste0(x$Transcription_Factor,"_",x$Cancer)]]["Up","Deceased"]
  
  ##define survival object
  survival<- Surv(tumor_surv[[x$Cancer]]$days, as.numeric(tumor_surv[[x$Cancer]]$vital_status))
  
  ###calculate KM curves
  #dicotomize zscore and pvalue and multiply them (0=not altered, 1=down-regulated, 2=up-regulated) 
  y<-(tumor_surv[[x$Cancer]][,as.character(x$Transcription_Factor)]>0) + 1
  y<-y*((tumor_surv[[x$Cancer]][,paste0(x$Transcription_Factor,".p")]<0.05) + 0)
  zy<-(tumor_surv[[x$Cancer]][,as.character(x$Transcription_Factor)]>0) + 0
  
  ##fill table
  if(length(unique(y))<2){
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UADA_Pvalue"]<-1
  }else if(length(unique(y))==2 & 2%in%y & 0%in%y){
    groups <- factor(y, levels=c(0,2),labels=c("TF Not Altered","TF Up-regulated"))
    table(tumor_surv[[x$Cancer]]$vital_status,groups)
    surv.KM <- survfit (survival ~ groups)
    dif<-survdiff (survival ~ groups)
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UADA_Pvalue"]<-signif(as.numeric(1 - pchisq(dif$chisq, length(dif$n) - 1), digits=3))
    
  }else if(length(unique(y))==2 & 1%in%y & 0%in%y){
    groups <- factor(y, levels=c(0,1),labels=c("TF Not Altered","TF Down-regulated"))
    table(tumor_surv[[x$Cancer]]$vital_status,groups)
    surv.KM <- survfit (survival ~ groups)
    dif<-survdiff (survival ~ groups)
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UADA_Pvalue"]<-signif(as.numeric(1 - pchisq(dif$chisq, length(dif$n) - 1), digits=3))
    
  }else{
    groups <- factor(y, levels=sort(unique(y)),labels=c("TF Not Altered","TF Down-regulated","TF Up-regulated"))
    table(tumor_surv[[x$Cancer]]$vital_status,groups)
    surv.KM <- survfit (survival ~ groups)
    dif<-survdiff (survival ~ groups)
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UADA_Pvalue"]<-signif(as.numeric(1 - pchisq(dif$chisq, length(dif$n) - 1), digits=3))
    
  }
  if(length(unique(zy))<2){
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UD_Pvalue"]<-1
  }
  else{
    #indicar los grupos
    groups <- factor(zy, levels=sort(unique(zy)),labels=c("TF Down-regulated","TF Up-regulated"))
    table(tumor_surv[[x$Cancer]]$vital_status,groups)
    surv.KM <- survfit (survival ~ groups)
    dif<-survdiff (survival ~ groups)
    KM_sig[which(KM_sig[,1]==x$Cancer & KM_sig[,2]==x$Transcription_Factor),"UD_Pvalue"]<-signif(as.numeric(1 - pchisq(dif$chisq, length(dif$n) - 1), digits=3))
  }
}

KM_sig<-KM_sig[as.numeric(KM_sig$UD_Pvalue)<0.05,]
##calculate adjusted pvalue
KM_sig$adj.UADA_Pvalue<-p.adjust(as.numeric(KM_sig$UADA_Pvalue))
KM_sig$adj.UD_Pvalue<-p.adjust(as.numeric(KM_sig$UD_Pvalue))
KM_sig[,2]<-as.character(KM_sig[,2])

##order acording pvalue
KM_sig<-KM_sig[order(as.numeric(KM_sig[,"UD_Pvalue"])),]
rownames(KM_sig)<-1:nrow(KM_sig)

###plot the result of significant KM plots for all cancers
for(cancer in unique(KM_sig$Cancer)){
  contador<-0
  print(cancer)
  tiff(filename = paste0("results/",cancer,"_resumen_KM.tiff"),width = 1450*2.8,height =  (ceiling((length(rownames(KM_sig[KM_sig$Cancer==cancer,]))+1+ceiling((length(rownames(KM_sig[KM_sig$Cancer==cancer,]))+1)/9))/9)*150*2.8)) 
  par(mfrow = c(ceiling((length(rownames(KM_sig[KM_sig$Cancer==cancer,]))+1+ceiling((length(rownames(KM_sig[KM_sig$Cancer==cancer,]))+1)/9))/9),9))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, cancer,cex = 4.6*2.8, col = "black")
  
  contador<-contador+1
  for(line in rownames(KM_sig[KM_sig$Cancer==cancer,])){ 
    
    if(contador!=0 & contador%%9==0){
      print("HOLA")
      plot.new()
      contador<-contador+1
    }
    contador<-contador+1
    x<-KM_sig[line,]
    
    #define survival object
    survival<- Surv(tumor_surv[[x$Cancer]]$days, as.numeric(tumor_surv[[x$Cancer]]$vital_status))
    
    #indicate groups
    zy<-(tumor_surv[[x$Cancer]][,as.character(x$Transcription_Factor)]>0) + 0
    zgroups <- factor(zy, levels=sort(unique(zy)),labels=c("TF Down-regulated","TF Up-regulated"))
    zsurv.KM <- survfit (survival ~ zgroups)
    zdif<-survdiff (survival ~ zgroups)
    
    #bold the title for thos TF wit significant adj.pvalue
    if(as.numeric(x$adj.UD_Pvalue)<0.05){
      plot(zsurv.KM,col=c("darkgreen","darkred"), lty=c(1,2), yaxt='n', xaxt='n', lwd=2.5, mark.time=TRUE,xlab= paste0("p-value: ", signif(as.numeric(1 - pchisq(zdif$chisq, length(zdif$n) - 1), digits=3))),cex.lab=1.5*2.8,font.lab=2)
      title(TF[which(gsub("::","_",TF)==x$Transcription_Factor)], cex.main=2.2*2.8, font.main=2)
      
    }else{
      plot(zsurv.KM,col=c("darkgreen","darkred"), lty=c(1,2), yaxt='n', xaxt='n', lwd=2.5, mark.time=TRUE,xlab= paste0("p-value: ", signif(as.numeric(1 - pchisq(zdif$chisq, length(zdif$n) - 1), digits=3))),cex.lab=1.5*2.8,font.lab=3)
      title(TF[which(gsub("::","_",TF)==x$Transcription_Factor)], cex.main=2.2*2.8, font.main=3)      
    }
    
  }
  dev.off()
}
