# Merge cells into X meta cells per transcriptional cluster

exprdf=read.table("data/sc.10x.counts.matrix", as.is=TRUE, header=TRUE,stringsAsFactors=F)
metadata=read.table("data/sc.10x.metadata.tsv", header=T, sep="\t", check.names=F,stringsAsFactors=F)
metadata=data.frame(cbind(metadata[,c('cell_ids')],paste(metadata[,c('dataset')],metadata[,c('reassigned_cluster')],sep="/")))
colnames(metadata)=c("cell","cluster")
metadata$cluster=gsub("PBMMC.*","control",metadata$cluster)

splitdf <- function(df, n) {
	indx <- matrix(seq_len(ncol(df)), ncol = n)
	lapply(seq_len(n), function(x) df[, indx[, x]])
}

n=30
mainDF=data.frame()
metaDF=data.frame()
for(clusterN in unique(metadata$cluster)) {
	tmpDF=exprdf[,as.numeric(rownames(subset(metadata,metadata$cluster==clusterN)))]
	testDF=splitdf(tmpDF,n)
	
	for(i in 1:n) {
		tmpDF2=data.frame(rowSums(testDF[[i]]))
		colnames(tmpDF2)=paste(clusterN,".",i,sep="")
		if(nrow(mainDF)==0) {
			mainDF=tmpDF2
			metaDF=data.frame(cbind(paste(clusterN,".",i,sep=""),clusterN))
			colnames(metaDF)=c("cell","cluster")
		} else {
			mainDF=cbind(mainDF,tmpDF2)
			if(clusterN=="control") {
				metatmp=data.frame(cbind(paste("control.",i,sep=""),"control"))
				colnames(metatmp)=c("cell","cluster")
				metaDF=rbind(metaDF,metatmp)
			} else {
				metatmp=data.frame(cbind(paste(clusterN,".",i,sep=""),clusterN))
				colnames(metatmp)=c("cell","cluster")
				metaDF=rbind(metaDF,metatmp)
		        }

		}
	}
}
write.table(mainDF,"data/metacell.rawcounts.tsv",quote=F,sep="\t")
write.table(metaDF,"data/metacell.metadata.tsv",row.names=F,col.names=F,quote=F,sep="\t")
