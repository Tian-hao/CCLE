#Rcode

library(RColorBrewer) 
hmcol <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
hmcol <- c('#00B050', '#7FCD6A', '#FFEB84', '#FF7642', '#FF0000')
hmcol <- colorRampPalette(hmcol)(100)
library(gplots)
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#a1 <- read.table('RNAseq_sum.FPKM.txt',header=T,row.names=1,stringsAsFactors=F)
#a1 <- read.table('immune_heatmap.txt',header=T,row.names=1,stringsAsFactors=F)
for (gene in c('7gene','5gene')){
  for (cell in c('noSCC','Adeno')){
a1 <- read.table(paste('CCLE_combine_',gene,'_',cell,'.txt',sep=''),header=T,row.names=1,stringsAsFactors=F,comment.char='@')

b1 <- a1[1:(nrow(a1)-2),]
for (i in 1:ncol(b1)){
  b1[,i] <- as.numeric(as.character(b1[,i]))
}
b1 <- log10(b1+0.01)
bo <- kmeans(t(b1),2)
bo2 <- cbind(t(b1),bo$cluster)
bo3 <- order(bo2[,ncol(bo2)])
b1 <- b1[,bo3]
c1 <- a1[(nrow(a1)-1):nrow(a1),]
for (i in 1:ncol(c1)){
  c1[,i] <- as.character(c1[,i])
}
c1 <- c1[,bo3]

dmethodlist = c("euclidean", "maximum", "manhattan", "canberra" ,"minkowski")
hmethodlist = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2", "centroid" ,"median")
for (dmethod in dmethodlist){
  for (hmethod in hmethodlist){
#hmethod = 'mcquitty'
#dmethod = 'manhattan'
png(paste('CCLE_cluster_noA549_',gene,'_',cell,'_',dmethod,'_',hmethod,'.png',sep=''),res=600,width=7200,height=4800)
#png(paste('CCLE_cluster_noA549_',gene,'_',cell,'.png',sep=''),res=600,width=7200,height=4800)
#pdf('CCLE_cluster.pdf',width=12,height=8)
#heatmap.2(as.matrix(log10(b1+0.01)),dendrogram='column',trace='none',scale="row",Rowv=F,col=hmcol,margins=c(8,8),key.par=list(mar=c(3,3,3,0)),ColSideColors=as.character(c1[1,]))
#abc <- heatmap.3(as.matrix(log10(b1+0.01)),dendrogram='column',trace='none',scale="row",Rowv=F,col=hmcol,margins=c(8,8),ColSideColors=t(c1),ColSideColorsSize=5)
#abc <- heatmap.3(as.matrix(log10(b1+0.001)),dendrogram='row',trace='none',scale="row",Colv=F,col=hmcol,margins=c(8,8),ColSideColors=t(c1),ColSideColorsSize=5)
abc <- heatmap.3(as.matrix(b1),dendrogram='column',trace='none',scale="row",Rowv=F,col=hmcol,margins=c(1,8),ColSideColors=t(c1),ColSideColorsSize=5,lmat=rbind(c(5,4),c(3,2),c(0,1),c(0,0)),lhei=c(2,4,1,1),labCol = FALSE,
  hclustfun = function(x) hclust(x, method=hmethod),
   distfun = function(x) dist(x, method=dmethod))
#abc <- heatmap.3(as.matrix(b1),dendrogram='none',trace='none',scale="row",Rowv=F,Colv=F,col=hmcol,margins=c(1,8),ColSideColors=t(c1),ColSideColorsSize=5,lmat=rbind(c(5,4),c(3,2),c(0,1),c(0,0)),lhei=c(2,4,1,1),labCol = FALSE)
par(xpd=T)
text(x=seq(0.231,0.939,length.out=ncol(b1)),y=rep(-0.01,ncol(b1)),colnames(b1)[abc$colInd],adj=1,srt=90,cex=0.4)
write.table(labels(abc$colDendrogram[[1]]),paste('CCLE_CXCL_low_',gene,'_',cell,'_',dmethod,'_',hmethod,'.txt',sep=''),quote=F,sep='\t')
write.table(labels(abc$colDendrogram[[2]]),paste('CCLE_CXCL_high_',gene,'_',cell,'_',dmethod,'_',hmethod,'.txt',sep=''),quote=F,sep='\t')
dev.off()
}}}}
