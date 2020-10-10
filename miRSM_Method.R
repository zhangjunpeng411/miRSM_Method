#################################################################################### 
## miRSM: Inferring miRNA sponge modules in human breast cancer
## November 14th, 2016, written by Junpeng Zhang
####################################################################################

## Get data from file ##
Read<-function(dataset){
		
  data<-read.csv(dataset, header=TRUE, sep=",")
  return(data)
	}

## Get data header from dataset ##
readHeader<-function(dataset){
  
  data<-read.csv(dataset, header=F)
  header<-character()
  for (i in 1:ncol(data)){
  header[i]=toString(data[1,i])
  }
  return(header)
}

## Constructing miRNA-target binding matrix using putative miRNA-target interactions ##
Bindingmatrix<-function(miRNA,mRNA,file){
  
  #the column name of file should be tagged as "mir" and "gene"
  data<-Read(file)
  mir=as.character(data$mir)
  #   mir<-paste("hsa-",sep="",mir);mir<-sub('r','R',mir)
  gene=as.character(data$gene)
  #   symbol<-geneSymbol(gene)   

  rep<-replicate(length(miRNA),mir)
  edge=matrix(F,length(miRNA)+length(mRNA),length(miRNA)+length(mRNA))
  for (i in 1:length(mir)){
    
    if (length(which(rep[i,]==miRNA)>0)){
      match1<-which(rep[i,]==miRNA,arr.ind=T)
      #gene: gene[i] #mirna: miRNA[match]
      rep2<-replicate(length(mRNA),gene[i])
      match2<-which(rep2==mRNA,arr.ind=T)
      edge[match1,match2+length(miRNA)]=T
    }
  }
  return(edge)
}


## Extract miRNA-target interactions by combining expression data and putative miRNA-target interactions ##
queryTargetFile<-function(miRNA,mRNA,file){
  
  #the column name of file should be tagged as "mir" and "gene"
  data<-Read(file)
  mir=as.character(data$mir)
  #   mir<-paste("hsa-",sep="",mir);mir<-sub('r','R',mir)
  gene=as.character(data$gene)
  #   symbol<-geneSymbol(gene)
  score=as.numeric(data$score)  

  rep<-replicate(length(miRNA),mir)
  edge=matrix(F,length(miRNA)+length(mRNA),length(miRNA)+length(mRNA))
  for (i in 1:length(mir)){
    
    if (length(which(rep[i,]==miRNA)>0)){
      match1<-which(rep[i,]==miRNA,arr.ind=T)
      #gene: gene[i] #mirna: miRNA[match]
      rep2<-replicate(length(mRNA),gene[i])
      match2<-which(rep2==mRNA,arr.ind=T)
      edge[match1,match2+length(miRNA)]=score[i]
    }
  }
  return(edge)
}

## Calculate colloration scores between miRNAs and mRNAs ##
ColScore<-function(data, header, cause, effect, targetbinding=NA){

num_miRNA=length(cause)
miR=header[1:num_miRNA]
mR=header[-(1:num_miRNA)]

miRNA=data[,cause]
mRNA=data[,effect]
corMat=cor(mRNA, miRNA, method="pearson")# miRNAs in columns and mRNAs in rows

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan=queryTargetFile(miR,mR,targetbinding)
edgeTargetScan=edgeTargetScan+t(edgeTargetScan)
weightededge=edgeTargetScan[effect,cause]

## binaryedge=edgeTargetScan!=0
## edge=binaryedge[effect,cause]

## Regarding the collaboration scores of miRNA-miRNA pairs without targetbinding as zero
## colScore=(0.5*corMat+0.5*weightededge)*edge 

# colScore=a*corMat+b*weightededge
colScore=0.5*corMat+0.5*weightededge # Default value of tuning parameters a and b, for example, a=0.5 and b=0.5

}

return(colScore)
}

## Identifying bicluster modules between miRNAs and mRNAs ##
BCMod<-function(colscore, BCmethod){

## Examples of biclust package ##
## test <- matrix(rbinom(400, 50, 0.4), 20, 20), res1 <- biclust(test, method=BCPlaid()) ## The method is from the publication "Turner H, Bailey T, Krzanowski W. Improved biclustering of microarray data demonstrated through systematic performance tests. Computational statistics & data analysis, 2005, 48(2): 235-254."
## Extract biocluster matrix: bicluster(test, res1), Extract miRNA and mRNA number: biclusternumber(res1)

## install.packages("biclust")
library("biclust")

## BCmethod contains 6 methods including BCBimax()(using binarize(colscore)), BCCC(), BCPlaid(), BCQuest(), BCSpectral() and BCXmotifs()(using discretize(colscore))
if (BCmethod=="BCBimax"){
    colscore=binarize(colscore)
    BCres=biclust(colscore, BCBimax(), minr=2, minc=3)
} else if (BCmethod=="BCCC") {
    BCres=biclust(colscore, BCCC())
} else if (BCmethod=="BCPlaid") {
    BCres=biclust(colscore, BCPlaid())
} else if (BCmethod=="BCQuest") {
    BCres=biclust(colscore, BCQuest(),ns=10, nd=10, sd=5, alpha=0.05)
} else if (BCmethod=="BCSpectral") {
    BCres=biclust(colscore, BCSpectral(), numberOfEigenvalues=3, minr=2, minc=3)
} else if (BCmethod=="BCXmotifs") {
    colscore=discretize(colscore)
    BCres=biclust(colscore, BCXmotifs(), ns=10, nd=10, sd=5, alpha=0.05)
} 

BCresnum=biclusternumber(BCres)

return(BCresnum)
}

## Inferring miRNA sponge interactions with positive correlation of p-value<0.01 ##
PCSponge<-function(data,colscore,BCresnum){

n=dim(data)[2]
ceRInt=matrix(NA,n*(n-1)/2,2)
PC=matrix(NA,n*(n-1)/2,4)

for (i in 1:(n-1)){
      for (j in (i+1):n){
           ceRInt[(i-1)*n+j-sum(1:i),1]=gsub("\\.","-",colnames(data)[i])
           ceRInt[(i-1)*n+j-sum(1:i),2]=gsub("\\.","-",colnames(data)[j])
           PC[(i-1)*n+j-sum(1:i),1]=cor.test(data[,i],data[,j])$estimate
           PC[(i-1)*n+j-sum(1:i),2]=cor.test(data[,i],data[,j])$p.value
           R1=length(which(colscore[BCresnum[i],]<=-0.3))
           R2=length(which(colscore[BCresnum[j],]<=-0.3))
           R3=length(intersect(which(colscore[BCresnum[i],]<=-0.3),which(colscore[BCresnum[j],]<=-0.3)))
           R4=dim(colscore)[2]
           PC[(i-1)*n+j-sum(1:i),3]=R3
           PC[(i-1)*n+j-sum(1:i),4]=1-phyper(R3,R2,R4-R2,R1)
      }
}

## Keep those miRNA sponge interactions with positive correlation of p-value<0.01
ceRInt=ceRInt[which((PC[,1]>0 & PC[,2]<0.01 & PC[,3]>2 & PC[,4]<0.01)=='TRUE'),]
PC=PC[which((PC[,1]>0 & PC[,2]<0.01 & PC[,3]>2 & PC[,4]<0.01)=='TRUE'),]

PCceRInt=cbind(ceRInt,PC)
return(PCceRInt)

}


