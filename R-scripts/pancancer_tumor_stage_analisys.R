#####################################################################################
#######################    Pan-cancer tumour stage analisys    ###################### 
#####################################################################################
#By: Matías Marín Falco                                                             #
#                                                                                   #
#The required data format for the analysis is explained in                          #
#the README file in my Git Hub repo: https://github.com/1mati1/pan-cancer-analysis/ #
#or in the babelomics repo:https://github.com/babelomics/TFTEA                      #     
#                                                                                   #
#The methods are also described in the paper:                                       #
#The pan-cancer pathological regulatory landscape                                   #
#####################################################################################

rm (list = ls ()) # clean the global environment

#define cancers to use in the analysis
cancers<-c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA","UCEC")


###get clinical data and merge it (must run this program on the main folder)

##get data for every cancer
library(curl)
for (cancer in cancers){
  #get TCGA clinical data (since all cancers selected are from TCGA project and ICGC clinical file has no info about tumor stage)
  x<-rawToChar(curl_fetch_memory(paste("https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/",tolower(cancer),"/bcr/biotab/clin/nationwidechildrens.org_clinical_patient_",tolower(cancer),".txt",sep=""))$content)
  assign(paste(cancer,"_tcga_clinical",sep=""),read.table(text=x, sep="\t",header=T,quote="",stringsAsFactors=F))
  
  #get ICGC clinical data
  x<-dir(path=cancer)[grep("clinical",dir(path=cancer))]#search the clinical file from each cancer cancer
  assign(paste(cancer,"_icgc_clinical",sep=""),read.table(paste(cancer,"/",x,sep=""), sep="\t",header=T,quote="",stringsAsFactors=F))
  
  #get donors name from a file(donors.txt is a one-column file which has the unique name of the donors of icgc clinical file. Each cancer folder has his own donors.txt file )
  assign(paste(cancer,"_donors",sep=""),read.table(paste(cancer,"/donors.txt",sep=""), sep="\t",header=F,quote="",stringsAsFactors=F))
  
}

##clean dataframes
for (cancer in cancers){
  assign(paste(cancer,"_icgc_clinical",sep=""),get(paste(cancer,"_icgc_clinical",sep=""))[which(get(paste(cancer,"_icgc_clinical",sep=""))$submitted_donor_id%in%get(paste(cancer,"_tcga_clinical",sep=""))$bcr_patient_barcode),])
  assign(paste(cancer,"_icgc_clinical",sep=""),get(paste(cancer,"_icgc_clinical",sep=""))[which(get(paste(cancer,"_icgc_clinical",sep=""))$icgc_donor_id%in%get(paste(cancer,"_donors",sep=""))[,1]),])
  assign(paste(cancer,"_tcga_clinical",sep=""),get(paste(cancer,"_tcga_clinical",sep=""))[which(get(paste(cancer,"_tcga_clinical",sep=""))$bcr_patient_barcode%in%get(paste(cancer,"_icgc_clinical",sep=""))$submitted_donor_id),])
  
}

##use donors names as rownames identifier
for (cancer in cancers){
  assign(paste(cancer,"_icgc_clinical",sep=""),within(get(paste(cancer,"_icgc_clinical",sep="")),{Stage<-0}))
  a<-get(paste(cancer,"_tcga_clinical",sep=""))
  rownames(a)<-a$bcr_patient_barcode
  assign(paste(cancer,"_tcga_clinical",sep=""),a)
  rm(a)
}

##gather all dataframes in lists
all_tcga <- list()
for(cancer in cancers){
  all_tcga[[cancer]] <- get(paste(cancer,"_tcga_clinical",sep=""))
}

all_icgc <- list()
for(cancer in cancers){
  all_icgc[[cancer]] <- get(paste(cancer,"_icgc_clinical",sep=""))
}



### add stage info to the ICGC dataframes and only keep those specimen with stage info(1,2,3 o 4)
types<-c()
for(cancer in cancers){
  #get stage info from TCGA 
  colnames(all_tcga[[cancer]])<-gsub("clinical_stage","ajcc_pathologic_tumor_stage",colnames(all_tcga[[cancer]]))
  all_icgc[[cancer]]$Stage<-all_tcga[[cancer]][all_icgc[[cancer]]$submitted_donor_id,"ajcc_pathologic_tumor_stage"]
  
  # get the different names used to describe the stages and simplify, if possible, to stage I, II, III or IV
  types<-c(types,unique(all_icgc[[cancer]]$Stage))
  for (i in c("Stage IA","Stage IB","Stage IC")){all_icgc[[cancer]]$Stage<-gsub(i,"Stage I",all_icgc[[cancer]]$Stage)}
  for (i in c("Stage IIA","Stage IIB","Stage IIC")){all_icgc[[cancer]]$Stage<-gsub(i,"Stage II",all_icgc[[cancer]]$Stage)}
  for (i in c("Stage IIIA","Stage IIIB","Stage IIIC","Stage IIIC1","Stage IIIC2","Stage III1","Stage III2")){ all_icgc[[cancer]]$Stage<-gsub(i,"Stage III",all_icgc[[cancer]]$Stage)}
  for (i in c("Stage IVA","Stage IVB","Stage IVC")){all_icgc[[cancer]]$Stage<-gsub(i,"Stage IV",all_icgc[[cancer]]$Stage)}
  
  # once transformed add it to icgc dataframes
  if(length(unique(all_tcga[[cancer]]$ajcc_pathologic_tumor_stage))>1){all_icgc[[cancer]]<-all_icgc[[cancer]][grep("Stage I",all_icgc[[cancer]]$Stage),]}
}

###remove all from global environment but the lists 
rm(list=setdiff(ls(),c("all_icgc","all_tcga","cancers","ptm")))

### make a list for each stage which indicates in each cancer which are the donors tumour samples(specimen) that belongs to the stage

#Create a list for the normal samples, which will be the reference
normals<-list()
for (cancer in cancers){
  normals[[cancer]]<-all_icgc[[cancer]][which(all_icgc[[cancer]]$specimen_type=="Normal - tissue adjacent to primary"),] #it's important to only use the tissue adjacent normal samples, so we try to avoid expression differences due to tissue differences
}

#create a list with stage I samples
stage_I<-list()
for (cancer in cancers){
  stage_I[[cancer]]<-all_icgc[[cancer]][which(all_icgc[[cancer]]$Stage=="Stage I"&all_icgc[[cancer]]$specimen_type=="Primary tumour - solid tissue"),] #in stages we will use solid tissue tumour samples 
  stage_I[[cancer]]<-stage_I[[cancer]][which(!duplicated(stage_I[[cancer]]$submitted_donor_id)),]
}
# create a list with stage II samples
stage_II<-list()
for (cancer in cancers){
  stage_II[[cancer]]<-all_icgc[[cancer]][which(all_icgc[[cancer]]$Stage=="Stage II"&all_icgc[[cancer]]$specimen_type=="Primary tumour - solid tissue"),]
  stage_II[[cancer]]<-stage_II[[cancer]][which(!duplicated(stage_II[[cancer]]$submitted_donor_id)),]
}
#create a list with stage III samples
stage_III<-list()
for (cancer in cancers){
  stage_III[[cancer]]<-all_icgc[[cancer]][which(all_icgc[[cancer]]$Stage=="Stage III"&all_icgc[[cancer]]$specimen_type=="Primary tumour - solid tissue"),]
  stage_III[[cancer]]<-stage_III[[cancer]][which(!duplicated(stage_III[[cancer]]$submitted_donor_id)),]
}
#create a list with stage IV samples
stage_IV<-list()
for (cancer in cancers){
  stage_IV[[cancer]]<-all_icgc[[cancer]][which(all_icgc[[cancer]]$Stage=="Stage IV"&all_icgc[[cancer]]$specimen_type=="Primary tumour - solid tissue"),]
  stage_IV[[cancer]]<-stage_IV[[cancer]][which(!duplicated(stage_IV[[cancer]]$submitted_donor_id)),]
}

stages<-c("stage_I","stage_II","stage_III","stage_IV")

######### apply differential expression (DE) between normal samples and different stages

for (stage in stages){
  print(stage)
  for (cancer in cancers){#apply DE in all cancers separately
    print (cancer)
    ##read necesary files
    reference<- normals[[cancer]]
    versus<- get(stage)[[cancer]]
    genes<- read.table(paste(cancer,"/genes.txt",sep=""), header=F, sep="\t",stringsAsFactors=F)$V1 # genes.txt has in one column all the unique genes from the icgc expression data file.
    ##remove "gene_id" identifier in genes.txt
    genes<-genes[genes!="gene_id"]
    
    ###create a matrix expression for each cancer
    rawexp<-mat.or.vec(nr=length(c(reference$icgc_specimen_id,versus$icgc_specimen_id)),nc=length(genes) )
    rownames(rawexp) <- c(reference$icgc_specimen_id,versus$icgc_specimen_id) #join as row names reference and case specimen names
    colnames(rawexp) <- genes #genes are in the columns
    
    ##read expression file
    library(ff)#as R loads in memory all object it may be hard for some computers to load the expression file, so this package helps to patch this file and read it
    raw<-dir(path=cancer)[grep("raw_exp_seq",dir(path=cancer))]#find the expression file
    expression<-read.table.ffdf(file=paste(cancer,"/",raw,sep=""), header=T, sep="\t")
    
    ##create a column in the expression file indicating for each expression data (row) if it is a tumour or normal specimen (in order to have a faster and easier acces to data later)
    specimen_type <- c( rep("TUMOR",nrow(versus)), rep("NORMAL",nrow(reference)))#indicate in an object which specimen are tumour and normal
    names(specimen_type) <- c( versus$icgc_specimen_id, reference$icgc_specimen_id)
    
    library("ffbase")
    library("ETLUtils")
    
    expression$sp_type <- with(expression[c('icgc_specimen_id')],
                              recoder(as.character(icgc_specimen_id), from = names(specimen_type), to = specimen_type))#create the column
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
    write.table(t,file=paste0(cancer,"/",stage,"_DE.tsv"),quote=F,sep="\t",row.names=T,col.names=F)
    
  }
}


rm (list = ls ())# once DE is performed we clean the environment again

stages<-c("stage_I","stage_II","stage_III","stage_IV")
tables<-c()

######### apply GSEA for each stage
for (stage in stages){
  
  ##load in environment the results of DE
  cancers<-c("BLCA","BRCA","COAD","HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","THCA","UCEC")
  ################################################################################
  #Transcription Factor Target Gene Analysis (TFTEA) FUNCTION
  
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
  
  ##load anotation file for TFTEA (in first column must be the gene targets  and in the second one the TF)
  anot <- as.matrix (read.table ("anotation_file.tsv", header = FALSE, sep = "\t", quote = "", colClasses = "character"))
  
  ###apply TFTEA function for every cancer
  for (cancer in cancers){
    ##Reading DE data
    stat <- read.table (paste(cancer,"/",stage,"_DE.tsv",sep=""), header = FALSE, sep = "\t", quote = "", as.is = TRUE)
    
    ##create, after DE, a vector for TFTEA
    names <- stat[,1]
    stat <- stat[,2]
    names (stat) <- names
    
    ## apply function
    res <- uvGsa (rankstat = stat, annotation = anot)
    
    ##assign results to an object and save them
    assign(paste(cancer,"_TFTEA",sep=""),res)
    outfile<-paste(cancer,"/result/babelomics.tsv",sep="")
    write.table (res, file = outfile, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
    
  }
  rm("stat","outfile","names","cancer","anot","res","uvGsa")
  
  ##create a matrix with adj.p.value from TFTEA
  pvalues<-mat.or.vec(nr=length(rownames(BLCA_TFTEA)),nc=length(cancers) )
  rownames(pvalues) <- rownames(BLCA_TFTEA) #rows are TFs
  colnames(pvalues) <- cancers #columns are cancers
  #fill matrix 
  for (cancer in cancers){
    pvalues[,cancer]<-get(paste(cancer,"_TFTEA",sep=""))$adj.X
  }
  #change some TFs names with confusing names(correct names searched manually), in order to query later in biomart
  rownames(pvalues)<-toupper(rownames(pvalues))
  olds<-c("CMYC","CFOS","PU1","GABP","AP2ALPHA","NRSF","AP2GAMMA","CJUN","EBF","NFKB","P300","BRG1","INI1","BAF155","TR4","BAF170")
  news<-c("MYC","FOS","SPI1","GABPA?GABPB1?2?","TFAP2A","REST","TFAP2C","JUN","EBF1","NFKB1?NFKB2","EP300","SMARCA4","SMARCB1","SMARCC1","NR2C2","SMARCC2")
  for(i in 1:length(olds)){
    oldgene<-olds[i]
    newgene<-news[i]
    rownames(pvalues)[grep(oldgene, rownames(pvalues))]=newgene
  }
  
  ###obtain the normal levels of expression fo TFs
  ##search Ensembl_ID of TFs
  library("biomaRt")
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="www.ensembl.org")#this function may cause some problems
  filters = listFilters(ensembl)
  TFs<-getBM(attributes=c("hgnc_symbol",'ensembl_gene_id' ), filters = "hgnc_symbol", values=TFs, mart=ensembl)
  #delete repeated TFs
  TFs<-TFs[!duplicated(TFs$hgnc_symbol),]
  rm("ensembl","filters")
  ##Look for the normal expression values for each TF in its corresponding cancer tissue, using the csv file downloaded for The Human Protein Atlas
  tissues<-read.table("normal_tissue.csv", header=T, sep=",",stringsAsFactors=F, comment.char="")
  #select available TFs and specify in each cancer which tissue to look at
  tissues<-tissues[ tissues$Gene %in% TFs$ensembl_gene_id,]
  Available_tissues<-unique(tissues$Tissue)
  assigned_tissue<-c("urinary bladder","breast","colon","NS","kidney","kidney","liver","lung","lung","thyroid gland","cervix, uterine")
  names(assigned_tissue)<-cancers
  assigned_tissue
  tissues<-tissues[ tissues$Tissue %in% assigned_tissue,]
  ##select the hicher expression values for each TFs in each cancer tissue
  values<-c(1:4)
  names(values)<-c("Not detected","Low","Medium","High")
  for (i in unique(tissues$Gene)){
    gene<-tissues[ tissues$Gene==i,]
    for (t in unique(assigned_tissue)){
      tis<-gene[gene$Tissue==t,]
      if(nrow(tis)>1){
        tissues<-tissues[! rownames(tissues) %in% rownames(tis[- which.max(values[tis$Level]),]),]
      }
    }
  }  
  rm("gene","tis","t","i","a","name","newgene","oldgene")
  
  ###create a matrix to indicate statistic value(z.X) from TFTEA
  odds<-mat.or.vec(nr=length(rownames(BLCA_TFTEA)),nc=length(cancers) )
  rownames(odds) <- rownames(BRCA_TFTEA)
  colnames(odds) <- cancers
  #fill table
  for (cancer in cancers){
    odds[,cancer]<-get(paste(cancer,"_TFTEA",sep=""))$z.X
  }
  rownames(odds)<-toupper(rownames(odds))
  #change some TFs names with confusing names(correct names searched manually), in order to query later in biomart
  for(i in 1:length(olds)){
    oldgene<-olds[i]
    newgene<-news[i]
    rownames(odds)[grep(oldgene, rownames(odds))]=newgene
  }
  rm("oldgene","i","newgene","cancer")
  
  ##plot a heatmap with the tissue expression levels for each TF in each cancer and its alteration pattern(blue=down-regulated, red=up-regulated, grey=not_differentialy_altered)
  rownames(TFs)<-TFs$ensembl_gene_id
  tissues<-cbind(name=TFs[tissues$Gene,'hgnc_symbol'],tissues)
  labels<-as.data.frame(mat.or.vec(nr=nrow(pvalues),nc=ncol(pvalues)))
  rownames(labels)<-rownames(pvalues)
  colnames(labels)<-colnames(pvalues)
  for (cancer in cancers){
    if(assigned_tissue[cancer]!="NS"){
      for(name in rownames(labels)){
        if(name %in% tissues$name){
          labels[name,cancer]<-tissues[which(tissues$name==name&tissues$Tissue==assigned_tissue[cancer]),'Level']
        }
      }
    }  
  }
  
  library(gplots)
  png(file=paste0("results/up_or_down_",stage,".png"),width = 1800, height = 900)
  up_or_down<-(pvalues<0.01) +(pvalues<0.01&odds>0)+0;up_or_down
  heatmap.2(up_or_down,trace="none",Rowv=T,Colv=T,tracecol=F,cellnote=labels,notecol="white",notecex=0.8,col=colorRampPalette(c("gray","blue","red"))(3))
  dev.off()
  
  #keep
  assign(paste0("pvalues_",stage), pvalues)
  assign(paste0("odds_", stage), odds)
  tables<-c(tables,paste0("pvalues_",stage),paste0("odds_", stage))
  
}

######### create a plot with the stage analysis results

## load necessary packages 
library(grid)
library(gplots)
library(ggplot2)
library(reshape2)
factors<-list() #list where we'll gather all individual TF plots

## load al available FTs from anotation file and correct names
FT<-as.matrix (read.table ("anotation_file.tsv", header = FALSE, sep = "\t", quote = "", colClasses = "character"))
ntargets<-as.character(table(FT[,2]))
names(ntargets)<-toupper(names(table(FT[,2])))
FT<-toupper(unique(FT[,2]))
olds<-c("CMYC","CFOS","PU1","GABP","AP2ALPHA","NRSF","AP2GAMMA","CJUN","EBF","NFKB","P300","BRG1","INI1","BAF155","TR4","BAF170")
news<-c("MYC","FOS","SPI1","GABPA?GABPB1?2?","TFAP2A","REST","TFAP2C","JUN","EBF1","NFKB1?NFKB2","EP300","SMARCA4","SMARCB1","SMARCC1","NR2C2","SMARCC2")
for(i in 1:length(olds)){
  oldgene<-olds[i]
  newgene<-news[i]
  FT[grep(oldgene, FT)]=newgene
}
FT<-setdiff(FT,c("NFKB1?NFKB2","SIN3AK20","GABPA?GABPB1?2?"))
ntargets[setdiff(names(ntargets),c("NFKB","GABP"))]

### create individual plots for each TF taking into account if they are in the edges plot also the axes (the plot has a proportion visual problem that I solved manually)

counter<-0#counter to control the axis

for (factor in FT){
  counter=counter+1
  #create the matrix with the TF info
  pvalue<-mat.or.vec(nr=length(stages),nc=length(canceres) )
  rownames(pvalue) <- stages
  colnames(pvalue) <- canceres
  odds<-mat.or.vec(nr=length(stages),nc=length(canceres) )
  rownames(odds) <- stages
  colnames(odds) <- canceres
  realpvalue<-mat.or.vec(nr=length(stages),nc=length(canceres) )
  rownames(realpvalue) <- stages
  colnames(realpvalue) <- canceres
  for (stage in stages){
    pvalue[stage,]<-(get(paste0("pvalues_",stage))[factor,])
    odds[stage,]<-(get(paste0("odds_",stage))[factor,])
  }
  matrix<-pmin(pvalue)      
  label<-round(matrix, 3)
  
  #when the analysis was made there wasn't enough infor for BLCA in stageI and LIHC and LUSC in stageIV
  label['stage_I','BLCA'] <- "NA"
  label['stage_IV','LIHC'] <- "NA"
  label['stage_IV','LUSC'] <- "NA"
  
  # combine pvalue matrix with odds matrix
  matrix<-(matrix<0.05) +(matrix<0.05&odds>0)+0
  
  u<-unique(as.numeric(matrix))
  
  melted_graf <- melt(t(matrix))
  melted_graf<-cbind(melted_graf,melt(t(label))[,3])
  colnames(melted_graf)<-c("Var1","Var2","value","labels")
  melted_graf$Var2 <- factor(melted_graf$Var2, levels = rev(levels(melted_graf$Var2)))#para ordenar como salen los stages en la grafica
  if (counter<13){
    factors[[factor]]<-ggplot(data = melted_graf, aes(x=Var1, y=Var2, fill=value))+ geom_tile(colour = "white")+ ##estos comandos para hacer un heatmap
      scale_fill_gradient2(low = "grey", high = "red", mid = "blue", midpoint = 1, limit = c(0,2)) +#para ajustar los colores
      theme_minimal()+ # minimal theme
      ggtitle(factor)+#poner el titulo
      theme(axis.ticks = element_blank(),plot.title = element_text(lineheight=.8, face="bold",size = rel(1)),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+ guides(fill=FALSE)+##quitar nombre de ejes y leyenda y modificar el titulo del grafico
      theme(plot.margin = unit(c(0,0,0,0), "cm"))+ylab(NULL)+xlab(NULL)+#reducir los margenes del plot
      geom_text(aes(label=labels),colour="white",size=4)#poner las etiquetas
  }
  else if (counter==13){
    factors[[factor]]<-ggplot(data = melted_graf, aes(x=Var1, y=Var2, fill=value))+ geom_tile(colour = "white")+ ##estos comandos para hacer un heatmap
      scale_fill_gradient2(low = "grey", high = "red", mid = "blue", midpoint = 1, limit = c(0,2)) +#para ajustar los colores
      theme_minimal()+ # minimal theme
      ggtitle(factor)+#poner el titulo
      theme(axis.ticks = element_blank(),plot.title = element_text(lineheight=.8, face="bold",size = rel(1)),axis.title.x = element_blank(),axis.title.y = element_blank())+ guides(fill=FALSE)+##quitar nombre de ejes y leyenda y modificar el titulo del grafico
      theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1),plot.margin = unit(c(0,0,0,0), "cm"))+ylab(NULL)+xlab(NULL)+#reducir los margenes del plot
      geom_text(aes(label=labels),colour="white",size=4)#poner las etiquetas
  }
  else if (counter%%13==0){
    factors[[factor]]<-ggplot(data = melted_graf, aes(x=Var1, y=Var2, fill=value))+ geom_tile(colour = "white")+ ##estos comandos para hacer un heatmap
      scale_fill_gradient2(low = "grey", high = "red", mid = "blue", midpoint = 1, limit = c(0,2)) +#para ajustar los colores
      theme_minimal()+ # minimal theme
      ggtitle(factor)+ # poner el titulo
      theme(axis.ticks = element_blank(),plot.title = element_text(lineheight=.8, face="bold",size = rel(1)),axis.text.y = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+ guides(fill=FALSE)+##quitar nombre de ejes y leyenda y modificar el titulo del grafico
      theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1),plot.margin = unit(c(0,0,0,0), "cm"))+ylab(NULL)+xlab(NULL)+#reducir los margenes del plot
      geom_text(aes(label=labels),colour="white",size=4)#poner las etiquetas
  }
  else{
    factors[[factor]]<-ggplot(data = melted_graf, aes(x=Var1, y=Var2, fill=value))+ geom_tile(colour = "white")+ ##estos comandos para hacer un heatmap
      scale_fill_gradient2(low = "grey", high = "red", mid = "blue", midpoint = 1, limit = c(0,2)) +#para ajustar los colores
      theme_minimal()+ # minimal theme
      ggtitle(factor)+#poner el titulo
      theme(axis.ticks = element_blank(),plot.title = element_text(lineheight=.8, face="bold",size = rel(1)),axis.text.y = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())+ guides(fill=FALSE)+##quitar nombre de ejes y leyenda y modificar el titulo del grafico
      theme(plot.margin = unit(c(0,0,0,0), "cm"))+ylab(NULL)+xlab(NULL)+#reducir los margenes del plot
      geom_text(aes(label=labels),colour="white",size=4)#poner las etiquetas
  }
  
  
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

tiff(filename = "results/summary_stages.tiff",width = 1375,height =  1208)
multiplot(plotlist = factors,cols=4)
dev.off()
