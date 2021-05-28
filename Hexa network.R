## AC infection  Angiostrongylus cantonensis infection for brain tissue of mice
#######  Select DEGs from RNA seq Data
setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
Gene_TPM=read.table("F:/parasite RNAseq Data Analysis/孙希-质谱数据/gene_TPM_genename.txt")  # , fill=T,header=T,sep=',', row.names =F, col.names = FALSE,fileEncoding = "UTF-8")
Gene_TPM=as.matrix(Gene_TPM)
colnames(Gene_TPM)=Gene_TPM[1,]

rownames(Gene_TPM)=NULL
rownames(Gene_TPM)=Gene_TPM[,14]
Gene_TPM=Gene_TPM[-1,-c(1,14)]
Gene_TPM=Gene_TPM[,c('ctrl_6h','ctrl_12h','ctrl_24h',"IL4_6h","IL4_12h","IL4_24h", "free_6h", "free_12h", "free_24h", "com_6h","com_12h","com_24h")]
head(Gene_TPM)

Gene_TPM['Chil3',]
Gene_TPM['Hexa',]

Gene_TPM['Arg1',]
Gene_TPM['Nos2',]

Gene=matrix(as.numeric(Gene_TPM),dim(Gene_TPM))
rownames(Gene)=rownames(Gene_TPM)
colnames(Gene)=colnames(Gene_TPM)

Gene_Ctr=Gene[,1:3]
Gene_IL4=Gene[,4:6]
Gene_Free=Gene[,7:9]
Gene_Com=Gene[,10:12]

###### DEG compared with 0 hour
DEG_0=matrix(0,1,3)
for (i in 1:dim(Gene)[1]) 
  {for (j in 1:3)

  {
    if ((max(Gene_IL4[i,])>=10) & (Gene_IL4[i,j]/mean(Gene_Ctr[i,]+1e-3)>5 || Gene_IL4[i,j]/mean(Gene_Ctr[i,]+1e-3)<0.2))
     {
       DEG_0[1]=DEG_0[1]+1
     }
  }
}

for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
  
{
  if ((max(Gene_Free[i,])>=10) & (Gene_Free[i,j]/mean(Gene_Ctr[i,]+1e-3)>5 || Gene_Free[i,j]/mean(Gene_Ctr[i,]+1e-3)<0.2))
  {
    DEG_0[2]=DEG_0[2]+1
  }
}
}

for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
  
{
  if ((max(Gene_Com[i,])>=10) & (Gene_Com[i,j]/mean(Gene_Ctr[i,]+1e-3)>5 || Gene_Com[i,j]/mean(Gene_Ctr[i,]+1e-3)<0.2))
  {
    DEG_0[3]=DEG_0[3]+1
  }
}
}


DEG_0


###### Normalization ###############
# Gene_IL4_Normalized=Gene_IL4
# Gene_Free_Normalized=Gene_Free
# Gene_Com_Normalized=Gene_Com
# 
# for (i in 1:dim(Gene)[1])
# {
#   for (j in 1:3)
#   Gene_IL4_Normalized[i,j]=Gene_IL4[i,j]/max(Gene_Ctr[i,j],1)
#   Gene_Free_Normalized[i,j]=Gene_IL4[i,j]/max(Gene_Ctr[i,j],1)
#   Gene_Com_Normalized[i,j]=Gene_IL4[i,j]/max(Gene_Ctr[i,j],1)
# }
# 
# 
# Gene_IL4=Gene_IL4_Normalized
# Gene_Free=Gene_Free_Normalized
# Gene_Com=Gene_Com_Normalized

###### DEG compared between adjencent time points  under each condition
DEG_num=matrix(0,1,3)

DEG_IL4=NULL
rname=NULL
for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
{
  for (k in j:3)
    
  {
    if (max(Gene_IL4[i,])>=10 & j!=k & (Gene_IL4[i,j]/(Gene_IL4[i,k]+1e-3)>2 || Gene_IL4[i,j]/(Gene_IL4[i,k]+1e-3)<0.5))
    {
      DEG_IL4=rbind(DEG_IL4,Gene_IL4[i,])
      rname=rbind(rname,rownames(Gene_IL4)[i])
    }
  }
}
}

rownames(DEG_IL4)=rname
DEG_IL4=unique(DEG_IL4)
DEG_num[1]=dim(DEG_IL4)[1]


DEG_Free=NULL
rname2=NULL
for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
{
  for (k in j:3)
    
  {
    if (max(Gene_Free[i,])>=10 & j!=k & (Gene_Free[i,j]/(Gene_Free[i,k]+1e-3)>2 || Gene_Free[i,j]/(Gene_Free[i,k]+1e-3)<0.5))
    {
      DEG_Free=rbind(DEG_Free,Gene_Free[i,])
      rname2=rbind(rname2,rownames(Gene_Free)[i])
    }
  }
}
}

rownames(DEG_Free)=rname2
DEG_Free=unique(DEG_Free)
DEG_num[2]=dim(DEG_Free)[1]


DEG_Com=NULL
rname=NULL
for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
{
  for (k in j:3)
    
  {
    if (max(Gene_Com[i,])>=10 & j!=k & (Gene_Com[i,j]/(Gene_Com[i,k]+1e-3)>2 || Gene_Com[i,j]/(Gene_Com[i,k]+1e-3)<0.5))
    {
      DEG_Com=rbind(DEG_Com,Gene_Com[i,])
      rname=rbind(rname,rownames(Gene_Com)[i])
    }
  }
}
}

rownames(DEG_Com)=rname
DEG_Com=unique(DEG_Com)
DEG_num[3]=dim(DEG_Com)[1]

########## DEG under AC-treated condition
DEG=NULL
rname=NULL
for (i in 1:dim(Gene)[1]) 
{for (j in 1:3)
{
  for (k in j:3)
    
  {
    if (max(Gene_Free[i,])>=10 & j!=k & (Gene_Free[i,j]/(Gene_Free[i,k]+1e-3)>5 || Gene_Free[i,j]/(Gene_Free[i,k]+1e-3)<0.2))
    {
      DEG=rbind(DEG,Gene_Free[i,])
      rname=rbind(rname,rownames(Gene_Free)[i])
    }
  }
}
}

rownames(DEG)=rname
DEG=unique(DEG)
dim(DEG)
DEG['Chil3',]



#####  Save DEGs 
setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
write.table(DEG, file="DEG_AC Condition.txt") 


TC_gene_Ctr=Gene_Ctr[rownames(DEG),]   # Temporally changing genes
TC_gene_IL4=Gene_IL4[rownames(DEG),] 
TC_gene_Free=Gene_Free[rownames(DEG),] 
TC_gene_Com=Gene_Com[rownames(DEG),] 


# setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
# write.table(TC_gene_Normalized, file="TC_gene_Normalized.txt",col.names = NA,sep = "\t") 


x1=(as.matrix(TC_gene_Free))
x2=(as.matrix(TC_gene_IL4))
x3=(as.matrix(TC_gene_Com))

x=cbind(x1,x2,x3)

setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
write.table(x, file="Expression_DEGs_ all conditions.txt",col.names = NA,sep = "\t") 
 
 #### LASSO regression gene network Model
 library(ggplot2)
 library(glmnet)
 # library(reshpae)
 
 A=matrix(0,dim(x)[1],dim(x)[1]+1)
 rownames(A)=rownames(x)
 colnames(A)=rownames(x)
 for (i in 1:dim(x)[1])
 {
 cvfit=cv.glmnet(t(x),x[i,],family="gaussian",alpha = 0.5, exclude=i)  # range of lambda: 10^(-5) to 10^(-1)
 Coef = coef(cvfit)
 Coef_min = coef(cvfit,s="lambda.min")
 A[i,]=as.numeric(Coef_min) #AA

 }

summary(A)
sum(A!=0)
dim(A)
########  AIC ti select significant edges
Theta=NULL
BIC=NULL
RSS=NULL

ind=1

for (theta in 10^seq(-5,0,0.1))  # seq(.)里的(to - from)/by 
{
  
  # Y_predict=matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2])
  X_predict=x
  
  for (i in dim(x)[1])
  {    
    D=A[,-1]
    D[abs(D)<=theta]=0
    x_new=x
    
    x_new[which(D[i,]==0),]=0
    
      # cvfit=cv.glmnet(t(x),x[i,],family="gaussian",alpha = 0.5, exclude=i)  # range of lambda: 10^(-5) to 10^(-1): lambda=10^seq(-1,-5,-0.1),
      
      # X_predict[i,]=predict(cvfit,t(x_new),s="lambda.min")  #for predict of glmnet: type="response",exact=TRUE)[,1]
      X_predict=D%*%x_new+A[,1]
  }
  
  
  p=sum(D!=0)
  E=x-X_predict   #Y_predict=D%*%x+A[,1]   #(t(t(x)%*%t(A[,-1]))+A[,1])
  N=dim(E)[1]*dim(E)[2]
  R=sum(E^2)  #/(dim(E)[1]*dim(E)[2])
  RSS[ind]=R
  BIC[ind]=log(1/N*R+1+log(2*pi))+2*p/N
  # BIC[ind]=N*log(R/N)+p*log(N)
  Theta[ind]=theta
  ind=ind+1
}

MSE=(RSS/N)^0.5

par(mfrow=c(2,2))
plot(log(Theta,10),MSE,type="b",main="MSE",sub='',xlim=c(-5,0), ylim=c(min(MSE), max(MSE)),xlab="theta",ylab='MSE')
plot(log(Theta,10),BIC,type="b",main="BIC",sub='',xlim=c(-5,0), ylim=c(min(BIC), max(BIC)),xlab="theta",ylab='BIC',pch=21,col="gray",cex=2,bg="blue",lty=4,lwd=2)

A=A[,-1]
A[A<1e-4]=0
########## 
setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
write.table(A, file="Model coefficient.txt",col.names = NA,sep = "\t") 


# #### Correlation
# # install.packages("ppcor")
# # library(ppcor)
# 
# 
# PCoef=matrix(0,dim(x)[1],dim(x)[1])
# Pvalue=matrix(0,dim(x)[1],dim(x)[1])
# 
# Data=as.data.frame(x)
# Gene_C=NULL
# Netij=matrix(0,dim(Data)[1],dim(Data)[1])
# for (i in 1:dim(Data)[1])
# {
#   for (j in 1:dim(Data)[1])
#   {
#     PCC=cor.test(x[i,],x[j,])
#     PCoef[i,j]=PCC$estimate
#     Pvalue[i,j]=PCC$p.value
#     if (abs(PCoef[i,j])>0.8 & Pvalue[i,j]<0.05)
#     {
#       Gene_C=rbind(Gene_C,rownames(A)[i])
#       Netij[i,j]=PCoef[i,j]
#     }
# 
#   }
# }
# 
# A=Netij
# rownames(A)=rownames(x)
# colnames(A)=rownames(x)
# setwd("F:/parasite RNAseq Data Analysis/Ym1 analysis") 
# write(Gene_C,file="Genes correlated with Chi3l3.xls",sep = "\t")

#### Integrate PPI information
B=matrix(0,dim(A)[1],dim(A)[2])
rownames(B)=rownames(A)
colnames(B)=colnames(A)

for ( i in 1: dim(A)[1])
{
  for (j in 1:dim(A)[2])

  {
    if (A[i,j]!=0 & (length((intersect(c(rownames(A)[i],colnames(A)[j]),as.matrix(PPI[,c(1,2)]))))>0  || length((intersect(c(colnames(A)[j],rownames(A)[i]),as.matrix(PPI[,c(1,2)]))))>0))
    {
      B[i,j]=A[i,j]
    }
    else
      B[i,j]=0

  }

}


### Reform txt data for cytoscape 
B=A
C=matrix(NA,nrow=sum(sum(B!=0)),ncol=3)
row=1
for (i in 1:dim(B)[1])
{
  for (j in 1:dim(B)[1])
  {
    if (B[i,j]!=0)
    {
      C[row,1]=rownames(B)[i]
      C[row,3]=colnames(B)[j]
      C[row,2]=B[i,j]
      row=row+1
    }
  }
}

write.table(C, file="Cytoscape_DEG Network.txt",quote=F,row.names=F,col.names = F,sep = "\t")

########################### Network topology analysis ###############################
install.packages("igraph")
library(igraph) # Load the igraph package

# setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
# C=read.table("Cytoscape_Differential Network.txt")


Net_DEG=as.data.frame(cbind(as.character(C[,1]),as.character(C[,3]),C[,2]))  # Input data for R analysis, data.frame
net=graph_from_data_frame(Net_DEG, directed = T,vertices = NULL)  # 

V(net)  #
E(net)

#Diameter
diam <- get_diameter(net, directed=T)

# Degree
deg <- degree(net, mode="all")
# plot(net, vertex.size=deg*3)
hist(deg, breaks=1:max(deg), main="Histogram of node degree")
deg.dist <- degree_distribution(net, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      xlab="Degree", ylab="Cumulative Frequency")

## neighbor

NB=neighbors(net, "Chil3", mode="all")

############### Entropy of each node
S=matrix(0,length(V(net)),3)

for (i in 1:length(V(net)))
{
  for (t in 1:3)
  {
    
    neigh <- neighbors(net, V(net)[i], mode="all")
    aij=matrix(0,1,length(neigh))
    for (j in 1:length(neigh))
    {
      aij[j]= x1[i,t]*x1[neigh[j],t]#*abs(sign(A[i,neigh[j]]))
    }
    
    S[i,t]=-sum(abs(aij)/sum(1e-10+abs(aij))*log(1e-10+abs(aij/sum(1e-10+abs(aij)))))/(1e-1+log(length(neigh)))
    
  }
  
}
rownames(S)=names(V(net))
colnames(S)=c("Entropy_6h","Entropy_12h","Entropy_24h")

S['Chil3',]

# deltaS=(S[,3]-S[,1])/(S[,1]+1e-10)  # changinf rate 

deltaS=2*S[,1]+S[,3]-3*S[,2]  # second order derivative rate

# deltaS=(S[,3]+S[,2])-2*S[,1]  # average changing rate
 deltaS=S[,3]-S[,1]  # average changing rate

 
 deltaS=(S[,3]-S[,2])*(S[,2]-S[,1])
names(deltaS)=rownames(S)

deltaS_sort=sort(deltaS,decreasing=T)

rank_deltaS=rank(deltaS)
order_deltaS=length(deltaS)-rank(deltaS)+1

DNE=cbind(x1,S,deltaS,rank_deltaS, order_deltaS)

DNE['Chil3',]
DNE[,-c(1:3)]

setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
write.table(DNE, file="Dynamic Network Entropy.txt",quote=F,row.names=T,col.names = T,sep = "\t")

############### Extract ACProtein-DEG interaction ##########################

setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
PPI=read.csv("F:/parasite RNAseq Data Analysis/孙希-质谱数据/String_DEGs_Proteins-MW-50.csv",sep=',')  # , fill=T,header=T,sep=',', row.names =F, col.names = FALSE,fileEncoding = "UTF-8")

head(PPI)
dim(PPI)
PPI=PPI[,c(1,2,15)]

# B=matrix(0,dim(A)[1],dim(A)[2])
# rownames(B)=rownames(A)
# colnames(B)=colnames(A)
# 
# for ( i in 1: dim(A)[1])
# {
#   for (j in 1:dim(A)[2])
#     
#   {
#     if (A[i,j]!=0 & (length((intersect(c(rownames(A)[i],colnames(A)[j]),as.matrix(PPI[,c(1,2)]))))>0  || length((intersect(c(colnames(A)[j],rownames(A)[i]),as.matrix(PPI[,c(1,2)]))))>0))
#     {
#       B[i,j]=A[i,j]
#     }
#     else
#       B[i,j]=0
#     
#   }
#   
# }



setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
ACProtein=read.csv("F:/parasite RNAseq Data Analysis/孙希-质谱数据/AC_Protein_MW_50.csv",sep=',')  # , fill=T,header=T,sep=',', row.names =F, col.names = FALSE,fileEncoding = "UTF-8")

PPI=as.matrix(PPI)
colnames(PPI)=NULL
ACProtein=as.matrix(ACProtein)

P_G=NULL
for (i in 1:dim(PPI)[1])
{
  if ( length(intersect(PPI[i,1], ACProtein))!=0 & length(intersect(PPI[i,2], ACProtein)) == 0 )
  {
    P_G=rbind(P_G,cbind(PPI[i,1],PPI[i,2]))
  }
  
  if ( length(intersect(PPI[i,2],ACProtein))!=0 & length(intersect(PPI[i,1],ACProtein))== 0)
    
  {
    P_G=rbind(P_G,cbind(PPI[i,2],PPI[i,1]))
  }
  
}

rownames(P_G)=P_G[,1]
P_G['Hexa',]
dim(P_G)

### Reform Protein_AC-DEG interaction for cytoscape 
P_G=cbind(P_G,matrix(5,dim(P_G)[1],1))

write.table(P_G, file="Cytoscape_Protein_AC-DEG interaction.txt",quote=F,row.names=F,col.names = F,sep = "\t")

##### 
####  Heatmap

install.packages("gplots") #下载gplots程序包
library(gtools) #加载gplots程序

install.packages("caTools") #下载gplots程序包
library(caTools)

library(gplots) #加载gplots程序


#### Normalization

for (i in 1:dim(TC_gene_Free)[1])
{
  TC_gene_Free[i,]=TC_gene_Free[i,]/(1e-10+max(TC_gene_Free[i,]))
  TC_gene_IL4[i,]=TC_gene_IL4[i,]/(1e-10+max(TC_gene_IL4[i,]))
  TC_gene_Com[i,]=TC_gene_Com[i,]/(1e-10+max(TC_gene_Com[i,]))
}

x1=(as.matrix(TC_gene_Free))
x2=(as.matrix(TC_gene_IL4))
x3=(as.matrix(TC_gene_Com))

dev.new()
heatmap.2(x1,col=redblue(100),distfun = dist,hclustfun = hclust,scale="none", Rowv=T,Colv=NA,dendrogram="row",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1,na.color=par("bg"))

dev.new()
heatmap.2(x2,col=greenred,distfun = dist,hclustfun = hclust,scale="none", Rowv=T,Colv=NA,dendrogram="row",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1,na.color=par("bg"))

dev.new()
heatmap.2(x3,col=greenred,distfun = dist,hclustfun = hclust,scale="none", Rowv=T,Colv=NA,dendrogram="row",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1,na.color=par("bg"))

dev.new()
heatmap.2(cbind(x1,x2,x3),col=greenred,distfun = dist,hclustfun = hclust,scale="row", Rowv=T,Colv=NA,dendrogram="row",key=T,keysize=1.5,trace="none",cexCol=0.5,cexRow=1,na.color=par("bg"))

install.packages("pheatmap")
library(pheatmap)

dev.new()
pheatmap(x1,cluster_row=T, cluster_col=F, cellwidth = 20, cellheight = 2.5,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize=9, fontsize_row=3,labRow=T) #自定义颜色

dev.new()
pheatmap(x2,cluster_row=T, cluster_col=F, cellwidth = 20, cellheight = 2.5,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize=9, fontsize_row=3,labRow=T) #自定义颜色


dev.new()
pheatmap(x3,cluster_row=T, cluster_col=F, cellwidth = 20, cellheight = 2.5,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize=9, fontsize_row=3,labRow=T) #自定义颜色




dev.new()


x_group=cbind(x1[,1],x2[,1],x3[,1],x1[,2],x2[,2],x3[,2],x1[,3],x2[,3],x3[,3])

# x_group=x_group[!is.nan(rowSums(x_group)),]


colnames(x_group)=c('Free_6h','IL4_6h','Com_6h','Free_12h','IL4_12h','Com_12h','Free_24','IL4_24h','Com_24h')
pheatmap(x_group,cluster_row=T, cluster_col=F, cellwidth = 20, cellheight = 2.5,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize=9, fontsize_row=3,labRow=T) #自定义颜色

####################  Protein Molecular Weight

setwd("F:/parasite RNAseq Data Analysis/孙希-质谱数据")  
Protein_MW=read.csv("F:/parasite RNAseq Data Analysis/孙希-质谱数据/Molecular Weight_All.csv")  # , ,fill=T, header=T,sep=' ', row.names =F, col.names = T,fileEncoding = "UTF-8")
head(Protein_MW)

Protein_MW=as.matrix(Protein_MW)
colnames(Protein_MW)=Protein_MW[1,]
rownames(Protein_MW)=Protein_MW[,1]

dim(Protein_MW)  # 1057, 3
  
MW=Protein_MW[,3]

hist(MW)
