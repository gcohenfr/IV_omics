#January 4th, 2011
#Analysis of continuous cardiac chronic rejection data using RIV

#Search for proteomics markers using microarray gene expression 
#	data as instruments

#--------------------------------------------------------------
#Proteomics Data Set
#--------------------------------------------------------------

library(MASS)
library(riv)

#--------------------------------------------------------------------------
setwd("D:/Documents/iFolder - GC/Stat UBC 2011/RIV in Proteomic/Analysis")

#Proteomics (log2 data)
#READ IN data.riv, gene symbols for these protein groups and Affy information

#data below contains onlty PGCs detected in at least 75% of all proteomic cohort
#Only patients in common with genomics were retained

data.riv=read.csv("proteomics data for RIV.csv",header=T,row.names=1)

colnames(data.riv)=paste(substr(colnames(data.riv),1,1),
		substr(colnames(data.riv),3,100),sep="-")

#data.riv[data.riv==0] <- NA
#data.riv <- log(data.riv,2)

#Genomics (log2 data)
#Impute missing values (not allowed for RIV)
source("imputeFn.txt")
library(SeqKnn)

data.riv <- my.impute(data.riv, "seqKNN", 3)


#RESPONSE
#Order of patients in genomics, proteomics and analysis table are the same

analysis.table=read.table("analysis table RIV.txt",header=T)

Y=as.vector(unlist(subset(analysis.table,subset=pID%in%common.patients,select=CRvsS)))

#Saved in "Data for RIV application.RData"
#-----------------------------------------------------------------------------

####################################################################
# OPTION 4: RIV using subsets of genes
####################################################################

#READ in gene symbols for these protein groups and Affy information

gene.pgc=read.table("Correlated features.txt",header=T)

#Retain the subset of PGCs with high-correlated genes

pgc_genes.withCorr=subset(gene.pgc,subset=Abs.Spearman>0)
pgc.withCorr=unique(as.vector(unlist(subset(gene.pgc,subset=Abs.Spearman>0,select=PGC))))

length(pgc.withCorr)
#[1] 64

data.rivCorr=subset(data.riv,subset=rownames(data.riv)%in%pgc.withCorr)

pgc.list=rownames(data.rivCorr)


####################################################################
#Output matrix: correlations, p-values both with robust and non-robust methods

#TRY GENES ONE BY ONE
i=i+1
gene.i=unique(as.vector(unlist(subset(pgc_genes.withCorr,
		subset=as.numeric(PGC)==as.numeric(pgc.list[i]),
		select=GeneSymbol))))

affy.i=unique(as.vector(unlist(subset(pgc_genes.withCorr,
		subset=GeneSymbol%in%gene.i,select=ProbeSetID))))

gene.i
affy.i

res.i=c()
for (k in 1:length(affy.i)){
W.inst=subset(data.H4,subset=rownames(data.H4)%in%affy.i)[c(k),]
X.end=subset(data.rivCorr,subset=rownames(data.rivCorr)==as.numeric(pgc.list[i]))

rob.riv=riv(Y,Xend=t(X.end),W=W.inst,method="robust")
#rob.riv=riv(Y,Xend=t(X.end),W=t(W.inst),method="robust")
res.i=c(res.i,rob.riv$Summary.Table[2,4])
}
res.i


rob.riv=riv(Y,Xend=t(X.end),W=t(W.inst),method="robust")
rob.riv$Summary.Table[2,4]

classical.oiv=riv(Y,Xend=t(X.end),W=t(W.inst),method="classical")
classical.oiv$Summary.Table[2,4]

#--------------------------------------------------
#CREATE A LIST FOR 8 significant proteins
#--------------------------------------------------

source("instrument list.R")

#CREATE A TABLE WITH UNIVARIATE RESULTS

output.frame=c()

for(j in 1:8){

W.inst=Winst.list[[j]]
X.end=Xend.list[[j]]

n.inst=nrow(W.inst)

W.inst=data.frame(SampleID=colnames(W.inst),t(W.inst))
X.end=data.frame(SampleID=colnames(X.end),t(X.end))
Y.resp=data.frame(SampleID=analysis.table$pID,Y=analysis.table$CRvsS)

Z.all=merge(merge(Y.resp,W.inst,by="SampleID",sort=F),X.end,by="SampleID",sort=F)

W.inst=Z.all[,3:(n.inst+2)]

X.end=Z.all[,(n.inst+3):ncol(Z.all)]

rob.riv=riv(Z.all$Y,Xend=X.end,W=W.inst,method="robust")
classical.oiv=riv(Z.all$Y,Xend=X.end,W=W.inst,method="classical")

#Without instruments
ttest.res=summary(rlm(Z.all$Y~X.end,na.omit=T,method="MM"))
MM.res=2*(1-pt(abs(ttest.res$coeff[2,3]),length(Y)-2))

output.frame=rbind(output.frame,cbind(substr(colnames(Z.all)[ncol(Z.all)],2,100),rob.riv$Summary.Table[2,4],
	classical.oiv$Summary.Table[2,4],
	MM.res,summary(lm(Z.all$Y~X.end))$coefficients[2,4],
	substr(colnames(Z.all)[3:(n.inst+2)],2,100)))
print(paste("j",j,sep="="))
}

colnames(output.frame)=c("PGC","RIV","OIV","MM-est","OLS","ProbeSetID")

#Merge with CORRELATION

cor.frame=c()
for(j in 1:8){
cor.frame=data.frame(rbind(cor.frame,cor.list[[j]]))
}


output.frame=merge(output.frame,cor.frame,all=T,sort=F,by="ProbeSetID")


write.table(output.frame,"RIV Results for Talk.csv",sep=",")


#---------------------------------------------------------------------

#---------------------------------------------------------------------------
#USE THE LIST TO PERFORMED A STEPWISE SELECTION
#	At this stage the strongest IV is selected, using the S-est of the correlation 
#	(Examine others later??)

set.seed(68743422)
W.instM=data.frame(SampleID=colnames(data.H4))
X.endM=data.frame(SampleID=colnames(data.H4))
r2.riv=c()
list.index=1:8
res.summary=c()
count.run=0


while(length(list.index)>0){
#change initial values

W.instI=W.instM
X.endI=X.endM
r2.riv=c()

for(j in list.index){

max.cor=as.vector(unlist(subset(cor.list[[j]],subset=cor.list[[j]][,"Abs.S"]==max(cor.list[[j]][,"Abs.S"]),
	select=ProbeSetID)))[1]
W.inst=subset(Winst.list[[j]],subset=rownames(Winst.list[[j]])==max.cor)
X.end=Xend.list[[j]]

W.inst=data.frame(SampleID=colnames(W.inst),t(W.inst))
X.end=data.frame(SampleID=colnames(X.end),t(X.end))

W.instM=merge(W.instI,W.inst,by="SampleID")
X.endM=merge(X.endI,X.end,by="SampleID")
Y.resp=data.frame(SampleID=analysis.table$pID,Y=analysis.table$CRvsS)

Z.all=merge(merge(Y.resp,W.instM,by="SampleID",sort=F),X.endM,by="SampleID",sort=F)

rob.rivM=riv(Z.all$Y,Xend=Z.all[,(count.run+4):ncol(Z.all)],W=Z.all[,3:(count.run+3)],method="robust")

b0.riv=rob.rivM$Summary.Table[1,1]
b1.riv=rob.rivM$Summary.Table[2:ncol(X.endM),1]

fitted.riv=b0.riv + as.matrix(Z.all[,(count.run+4):ncol(Z.all)]) %*% as.vector(b1.riv)
res=R2.fun(Z.all$Y,rob.rivM$weight,as.vector(fitted.riv))

r2.riv=rbind(r2.riv,c(colnames(X.endM)[ncol(X.endM)],res))

print(paste("j",j,sep="="))
}
colnames(r2.riv)=c("PGC","rob.R2")

max.r2=which(r2.riv[,"rob.R2"]==max(r2.riv[,"rob.R2"]))

#store result
res.summary=rbind(res.summary,r2.riv[max.r2,])

#Selected variable
max.cor=as.vector(unlist(subset(cor.list[[list.index[max.r2]]],
	subset=cor.list[[list.index[max.r2]]][,"Abs.S"]==max(cor.list[[list.index[max.r2]]][,"Abs.S"]),
	select=ProbeSetID)))[1]

W.inst=subset(Winst.list[[list.index[max.r2]]],subset=rownames(Winst.list[[list.index[max.r2]]])==max.cor)
X.end=Xend.list[[list.index[max.r2]]]

W.inst=data.frame(SampleID=colnames(W.inst),t(W.inst))
X.end=data.frame(SampleID=colnames(X.end),t(X.end))

W.instM=merge(W.instI,W.inst,by="SampleID")
X.endM=merge(X.endI,X.end,by="SampleID")

#Refine search
list.index=list.index[-max.r2]
count.run=count.run+1

}

> res.summary
     PGC    rob.R2             
[1,] "X83"  "0.224165680623089"
[2,] "X139" "0.342998229965339"
[3,] "X98"  "0.385546679416441"
[4,] "X7"   "0.176844354580437"
[5,] "X150" "0.171811792089277"

#Using spearman-corr gives a max R^2 of .46

#Try multi-IV partitioning the data to get an S-est



#PLOT OF PROTEIN RATIOS

W139=as.vector(Winst.list[[6]][3,])
X139=as.vector(Xend.list[[6]])


windows(h=7,w=7,pointsize=16)

rr=riv(Y,Xend=X139,W=W139,method="robust")


plot(X139,Y,xlab="Protein X", ylab="Stenosis",pch=19,main="CAV Data")
#identify(X.endI[,(j+1)],Y)
pos2=c( 9, 35 ,37)
points(X139[pos2],Y[pos2],pch=19,cex=1.2,col=4)
abline(lm(Y~X139),col=2,lty=2)
#a=rr$Summary[1,1]+90
#b=rr$Summary[2,1]-120


#abline(a,b,col=3,lty=3,lwd=2)

windows(h=7,w=7,pointsize=14)

plot(W139,X139,xlab="217232_x_at", ylab="PGC 139",pch=19)
#identify(W139,X139)
pos3=c( 9 ,37)
points(W139[pos3],X139[pos3],pch=19,cex=1.2,col=4)

windows(h=7,w=7,pointsize=14)

plot(W116,X116,xlab="204561_x_at", ylab="PGC 139",pch=19)
#identify(W116,X116)
pos3=c( 23,10 ,37)
points(W116[pos3],X116[pos3],pch=19,cex=1.2,col=4)





