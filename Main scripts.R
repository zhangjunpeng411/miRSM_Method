#################################################################################### 
## Main running scripts of miRSM: Inferring miRNA sponge modules in human breast cancer
## November 14th, 2016, written by Junpeng Zhang
####################################################################################
source("miRSM_Method.R")
datacsv="miRNA(453)-mRNA(11157).csv"
data=Read(datacsv)
data=scale(data)
header=readHeader(datacsv)
cause=1:453
effect=454:11610
targetbinding="TargetScan_7.1.csv"
colscore=ColScore(data, header, cause, effect, targetbinding)

## BCmethod contains BCBimax(), BCCC(), BCPlaid(), BCQuest(), BCSpectral() and BCXmotifs() 
BCresnum=BCMod(colscore,BCmethod="BCPlaid")

num_bicluster=length(BCresnum)
num_miRNA=length(cause)
ModExp=lapply(1:num_bicluster, function(i) data[,BCresnum[[i]]$Rows+num_miRNA]) ## Extract mRNA expression data for each bicluster module
ModCol=lapply(1:num_bicluster, function(i) colscore[BCresnum[[i]]$Rows,BCresnum[[i]]$Cols]) ## Extract miRNA-mRNA collarboration score for each bicluster module

miRSponge=lapply(1:num_bicluster, function(i) PCSponge(ModExp[[i]],colscore,BCresnum[[i]]$Rows))

## install.packages("igraph")
library("igraph")

## Graph of PPIs with HPRD_9
PPI_HPRD=Read("HPRD.csv")
PPI_HPRD_graph=make_graph(c(t(PPI_HPRD)),directed = FALSE)

## Graph of TF-target interactions with HTRI
TFtarget_HTRI=Read("HTRI.csv")
TFtarget_HTRI_graph=make_graph(c(t(TFtarget_HTRI)),directed = FALSE)

## Graph of miRNA sponge interactions
miRSponge_graph=lapply(1:num_bicluster, function(i) make_graph(c(t(miRSponge[[i]][,1:2])),directed = FALSE))

## Remove miRNA sponge interactions existing in PPIs and TF-target interactions
miRSponge_PPI_graph=lapply(1:num_bicluster, function(i) miRSponge_graph[[i]] %m% PPI_HPRD_graph)
miRSponge_PPI_TFtarget_graph=lapply(1:num_bicluster, function(i) miRSponge_PPI_graph[[i]] %m% TFtarget_HTRI_graph)
miRSponge_PPI_TFtarget=lapply(1:num_bicluster, function(i) as_data_frame(miRSponge_PPI_TFtarget_graph[[i]], what="edges"))

## Replace "." of miRNA and gene names using "-"
miRSpongeModCol=lapply(1:num_bicluster, function(i) ModCol[[i]][which(gsub("\\.","-",rownames(ModCol[[i]])) %in% unique(c(miRSponge_PPI_TFtarget[[i]][,1],miRSponge_PPI_TFtarget[[i]][,2])) == "TRUE"),])
miRSpongeModCol_rownames=lapply(1:num_bicluster, function(i) gsub("\\.","-",rownames(miRSpongeModCol[[i]])))
miRSpongeModCol_colnames=lapply(1:num_bicluster, function(i) gsub("\\.","-",colnames(miRSpongeModCol[[i]])))

## Extract expression data of miRNA sponges and miRNAs for each module
miRSpongeExp=lapply(1:num_bicluster, function(i) data[,which(gsub("\\.","-",rownames(colscore)) %in% miRSpongeModCol_rownames[[i]] )+num_miRNA])
miRExp=lapply(1:num_bicluster, function(i) data[,which(gsub("\\.","-",colnames(colscore)) %in% miRSpongeModCol_colnames[[i]] )])

## Differential network analysis of four miRNA sponge modules
Diff_14=as_data_frame(miRSponge_PPI_TFtarget_graph[[1]] %m% miRSponge_PPI_TFtarget_graph[[4]])
Diff_41=as_data_frame(miRSponge_PPI_TFtarget_graph[[4]] %m% miRSponge_PPI_TFtarget_graph[[1]])
Diff_23=as_data_frame(miRSponge_PPI_TFtarget_graph[[2]] %m% miRSponge_PPI_TFtarget_graph[[3]])
Diff_32=as_data_frame(miRSponge_PPI_TFtarget_graph[[3]] %m% miRSponge_PPI_TFtarget_graph[[2]])
Diff=rbind(Diff_14,Diff_41,Diff_23,Diff_32)
write.csv(Diff, "Differential miRSponge network.csv")

UnimiRSponge=miRSponge[[1]]

for (i in 1:num_bicluster){
    if ((i+1)<=num_bicluster){
         UnimiRSponge=rbind(UnimiRSponge,miRSponge[[i+1]])
}
}


## Assemble all miRNA sponge interactions of each bicluster module
UnimiRSponge=unique(data.frame(UnimiRSponge))

## Network analysis using R package "igraph"
# %m% is difference of two graph, %u% is the union of two graphs
UnimiRSponge_graph=make_graph(c(t(UnimiRSponge[,1:2])),directed = FALSE)
UnimiRSponge_PPI_graph=UnimiRSponge_graph %m% PPI_HPRD_graph
UnimiRSponge_PPI_TFtarget_graph=UnimiRSponge_PPI_graph %m% TFtarget_HTRI_graph

## Graph of validated miRNA sponge interactions with miRSponge
Validated_miRSponge=Read("Validated_miRSponge.csv")
Validated_miRSponge_graph=make_graph(c(t(Validated_miRSponge)),directed = FALSE)

## Validated miRNA sponge interactions using the third-party database miRSponge 
Intersection_Validated=UnimiRSponge_PPI_TFtarget_graph %s% Validated_miRSponge_graph

## 1228533 interactions (-0.3, a=0.55, b=0.45)
## [1] PTEN --SERINC1 LRCH1--PTEN    ZEB2 --PTEN    MBNL1--PTEN   
## [5] KLF6 --PTEN

## 577544 interactions (-0.3, a=0.5, b=0.5)
## [1] PTEN --SERINC1 LRCH1--PTEN    KLF6 --PTEN    ZEB2 --PTEN   
## [5] MBNL1--PTEN 

## 169617 interactions (-0.3, a=0.45, b=0.55)
## [1] PTEN --ZEB2    PTEN --SERINC1 LRCH1--PTEN    KLF6 --PTEN   
## [5] MBNL1--PTEN

## 17972 interactions (-0.3, a=0.4, b=0.6)
## [1] PTEN--LRCH1 KLF6--PTEN  

## 2029912 interactions (-0.3, a=0.6, b=0.4)
## [1] PTEN --AFF1    PTEN --ZEB2    PTEN --SERINC1 MBNL1--PTEN   
## [5] LRCH1--PTEN    KLF6 --PTEN    FN1  --VCAN

## 2033124 interactions (-0.25, a=0.5, b=0.5)
## [1] PTEN --AFF1    PTEN --ZEB2    PTEN --SERINC1 MBNL1--PTEN   
## [5] LRCH1--PTEN    KLF6 --PTEN    FN1  --VCAN

## 61161 interactions (-0.35, a=0.5, b=0.5)
## /

## 1041 interactions (-0.4, a=0.5, b=0.5)
## /

## 84 interactions (-0.45, a=0.5, b=0.5)
## /

## Save the results of miRNA sponge modules and validated miRNA sponge interactions
save.image("miRSponge_ColScore_p.value=0.01_BCPlaid.RData")

## BRCA-related miRNAs and genes, and hallmark genes
BRCA_GENE=read.csv("BRCA_GENE.csv", header=FALSE, sep=",")
BRCA_miRNA=read.csv("BRCA_miRNA.csv", header=FALSE, sep=",")
Hallmark_GENE=read.csv("Hallmark_GENE.csv", header=FALSE, sep=",")

## Convert the graph of miRNA sponge interactions into data frame
UnimiRSponge_PPI_TFtarget=as_data_frame(UnimiRSponge_PPI_TFtarget_graph, what="edges")

## Extract BRCA-related and Hallmark-related miRNA sponge interactions
BRCA_Interactions=UnimiRSponge_PPI_TFtarget[intersect(which(UnimiRSponge_PPI_TFtarget[,1]%in%BRCA_GENE[,1]),which(UnimiRSponge_PPI_TFtarget[,2]%in%BRCA_GENE[,1])),]
Hallmark_Interactions=UnimiRSponge_PPI_TFtarget[intersect(which(UnimiRSponge_PPI_TFtarget[,1]%in%Hallmark_GENE[,1]),which(UnimiRSponge_PPI_TFtarget[,2]%in%Hallmark_GENE[,1])),]

## Extract PTEN-related miRNA sponge interactions
PTEN_Interactions=UnimiRSponge_PPI_TFtarget[which(UnimiRSponge_PPI_TFtarget[,1]=="PTEN")%u%which(UnimiRSponge_PPI_TFtarget[,2]=="PTEN"),]
write.csv(PTEN_Interactions, "PTEN_Interactions.csv")

## Calculate the percentage of BRCA miRNAs, genes and hallmark genes for each module
length(intersect(miRSpongeModCol_rownames[[1]],as.matrix(BRCA_GENE)))/length(miRSpongeModCol_rownames[[1]])
length(intersect(miRSpongeModCol_rownames[[2]],as.matrix(BRCA_GENE)))/length(miRSpongeModCol_rownames[[2]])
length(intersect(miRSpongeModCol_rownames[[3]],as.matrix(BRCA_GENE)))/length(miRSpongeModCol_rownames[[3]])
length(intersect(miRSpongeModCol_rownames[[4]],as.matrix(BRCA_GENE)))/length(miRSpongeModCol_rownames[[4]])

length(intersect(miRSpongeModCol_colnames[[1]],as.matrix(BRCA_miRNA)))/length(miRSpongeModCol_rownames[[1]])
length(intersect(miRSpongeModCol_colnames[[2]],as.matrix(BRCA_miRNA)))/length(miRSpongeModCol_rownames[[2]])
length(intersect(miRSpongeModCol_colnames[[3]],as.matrix(BRCA_miRNA)))/length(miRSpongeModCol_rownames[[3]])
length(intersect(miRSpongeModCol_colnames[[4]],as.matrix(BRCA_miRNA)))/length(miRSpongeModCol_rownames[[4]])

length(intersect(miRSpongeModCol_rownames[[1]],as.matrix(Hallmark_GENE)))/length(miRSpongeModCol_rownames[[1]])
length(intersect(miRSpongeModCol_rownames[[2]],as.matrix(Hallmark_GENE)))/length(miRSpongeModCol_rownames[[2]])
length(intersect(miRSpongeModCol_rownames[[3]],as.matrix(Hallmark_GENE)))/length(miRSpongeModCol_rownames[[3]])
length(intersect(miRSpongeModCol_rownames[[4]],as.matrix(Hallmark_GENE)))/length(miRSpongeModCol_rownames[[4]])

## Extract miRNA-target interactions for each module
ModmiRNA_target1=cbind(miRSpongeModCol_colnames[[1]][which((miRSpongeModCol[[1]]<(-0.3))==1,arr.ind=TRUE)[,2]],miRSpongeModCol_rownames[[1]][which((miRSpongeModCol[[1]]<(-0.3))==1,arr.ind=TRUE)[,1]])
ModmiRNA_target2=cbind(miRSpongeModCol_colnames[[2]][which((miRSpongeModCol[[2]]<(-0.3))==1,arr.ind=TRUE)[,2]],miRSpongeModCol_rownames[[2]][which((miRSpongeModCol[[2]]<(-0.3))==1,arr.ind=TRUE)[,1]])
ModmiRNA_target3=cbind(miRSpongeModCol_colnames[[3]][which((miRSpongeModCol[[3]]<(-0.3))==1,arr.ind=TRUE)[,2]],miRSpongeModCol_rownames[[3]][which((miRSpongeModCol[[3]]<(-0.3))==1,arr.ind=TRUE)[,1]])
ModmiRNA_target4=cbind(miRSpongeModCol_colnames[[4]][which((miRSpongeModCol[[4]]<(-0.3))==1,arr.ind=TRUE)[,2]],miRSpongeModCol_rownames[[4]][which((miRSpongeModCol[[4]]<(-0.3))==1,arr.ind=TRUE)[,1]])

## Validated miRNA-target interactions with strong evidence from miRTarBase
targetbinding="miRTarBase_Strong.csv"

## Validation of miRNA-target interactions for each module
miR=miRSpongeModCol_colnames[[1]]
mR=miRSpongeModCol_rownames[[1]]
edgeTargetScan=Bindingmatrix(miR,mR,targetbinding)
edgeTargetScan=edgeTargetScan+t(edgeTargetScan)
binaryedge=edgeTargetScan!=0
edge=binaryedge[(length(miR)+1):(length(miR)+length(mR)),1:length(miR)]
which(edge*(miRSpongeModCol[[1]]<(-0.3))==1)
Validated_miRSponge_col1= miRSpongeModCol_colnames[[1]][which(edge*(miRSpongeModCol[[1]]<(-0.3))==1,arr.ind=TRUE)[,2]]
Validated_miRSponge_row1= miRSpongeModCol_rownames[[1]][which(edge*(miRSpongeModCol[[1]]<(-0.3))==1,arr.ind=TRUE)[,1]]

miR=miRSpongeModCol_colnames[[2]]
mR=miRSpongeModCol_rownames[[2]]
edgeTargetScan=Bindingmatrix(miR,mR,targetbinding)
edgeTargetScan=edgeTargetScan+t(edgeTargetScan)
binaryedge=edgeTargetScan!=0
edge=binaryedge[(length(miR)+1):(length(miR)+length(mR)),1:length(miR)]
which(edge*(miRSpongeModCol[[2]]<(-0.3))==1)
Validated_miRSponge_col2= miRSpongeModCol_colnames[[2]][which(edge*(miRSpongeModCol[[2]]<(-0.3))==1,arr.ind=TRUE)[,2]]
Validated_miRSponge_row2= miRSpongeModCol_rownames[[2]][which(edge*(miRSpongeModCol[[2]]<(-0.3))==1,arr.ind=TRUE)[,1]]

miR=miRSpongeModCol_colnames[[3]]
mR=miRSpongeModCol_rownames[[3]]
edgeTargetScan=Bindingmatrix(miR,mR,targetbinding)
edgeTargetScan=edgeTargetScan+t(edgeTargetScan)
binaryedge=edgeTargetScan!=0
edge=binaryedge[(length(miR)+1):(length(miR)+length(mR)),1:length(miR)]
which(edge*(miRSpongeModCol[[3]]<(-0.3))==1)
Validated_miRSponge_col3= miRSpongeModCol_colnames[[3]][which(edge*(miRSpongeModCol[[3]]<(-0.3))==1,arr.ind=TRUE)[,2]]
Validated_miRSponge_row3= miRSpongeModCol_rownames[[3]][which(edge*(miRSpongeModCol[[3]]<(-0.3))==1,arr.ind=TRUE)[,1]]

miR=miRSpongeModCol_colnames[[4]]
mR=miRSpongeModCol_rownames[[4]]
edgeTargetScan=Bindingmatrix(miR,mR,targetbinding)
edgeTargetScan=edgeTargetScan+t(edgeTargetScan)
binaryedge=edgeTargetScan!=0
edge=binaryedge[(length(miR)+1):(length(miR)+length(mR)),1:length(miR)]
which(edge*(miRSpongeModCol[[4]]<(-0.3))==1)
Validated_miRSponge_col4= miRSpongeModCol_colnames[[4]][which(edge*(miRSpongeModCol[[4]]<(-0.3))==1,arr.ind=TRUE)[,2]]
Validated_miRSponge_row4= miRSpongeModCol_rownames[[4]][which(edge*(miRSpongeModCol[[4]]<(-0.3))==1,arr.ind=TRUE)[,1]]

## Differentially expressed miRNAs and genes using limma package
DEmiRNA=read.csv("DEmiRNA_1E-02.csv",header=FALSE,sep=",")
DEG=read.csv("DEG_1E-04.csv",header=FALSE,sep=",")

## Calculate the percentage of differentially expressed miRNAs and genes for each module
length(intersect(miRSpongeModCol_rownames[[1]],as.matrix(DEG)))/length((miRSpongeModCol_rownames[[1]]))
length(intersect(miRSpongeModCol_rownames[[2]],as.matrix(DEG)))/length((miRSpongeModCol_rownames[[2]]))
length(intersect(miRSpongeModCol_rownames[[3]],as.matrix(DEG)))/length((miRSpongeModCol_rownames[[3]]))
length(intersect(miRSpongeModCol_rownames[[4]],as.matrix(DEG)))/length((miRSpongeModCol_rownames[[4]]))

length(intersect(miRSpongeModCol_colnames[[1]],as.matrix(DEmiRNA)))/length((miRSpongeModCol_colnames[[1]]))
length(intersect(miRSpongeModCol_colnames[[2]],as.matrix(DEmiRNA)))/length((miRSpongeModCol_colnames[[2]]))
length(intersect(miRSpongeModCol_colnames[[3]],as.matrix(DEmiRNA)))/length((miRSpongeModCol_colnames[[3]]))
length(intersect(miRSpongeModCol_colnames[[4]],as.matrix(DEmiRNA)))/length((miRSpongeModCol_colnames[[4]]))

