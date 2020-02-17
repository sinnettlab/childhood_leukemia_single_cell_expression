set.seed(123456789)

# Clustering solutions

startRes=0.1
endRes=3
step=0.1
for(i in seq(startRes,endRes,step)) {
resCluster=i
cancer.integrated <- FindClusters(cancer.integrated, dims.use = 1:nmbDims, force.recalc=TRUE, save.SNN = T, resolution = resCluster)
}


# Adjusted rand index

ARImatrix=matrix(nrow=length(seq(startRes,endRes,step)),ncol=length(seq(startRes,endRes,step)))
for(i in 1:length(seq(startRes,endRes,step))) {
	for(j in 1:length(seq(startRes,endRes,step))) {
		if(i!=j) {
			vecA=cancer.integrated@meta.data[,paste("res.",seq(startRes,endRes,step)[i],sep="")]
			vecB=cancer.integrated@meta.data[,paste("res.",seq(startRes,endRes,step)[j],sep="")]
			ARImatrix[i,j]=adjustedRandIndex(vecA,vecB)
		}
	}
}

ARImean=NULL
ARIsd=NULL
ARInmbClusters=NULL
for(i in 1:length(seq(startRes,endRes,step))) {
	ARImean=c(ARImean,mean(ARImatrix[i,-i]))
	ARIsd=c(ARIsd,sd(ARImatrix[i,-i]))
	ARInmbClusters=c(ARInmbClusters,length(unique(cancer.integrated@meta.data[,paste("res.",seq(startRes,endRes,step)[i],sep="")])))
}
optimalRes=seq(startRes,endRes,step)[which(ARImean==max(ARImean))[1]]
ARIdf=data.frame(resolution=seq(startRes,endRes,step), ARI_mean=ARImean, ARI_sd=ARIsd, number_clusters=ARInmbClusters)
ARIdf=melt(ARIdf,id='resolution')

# Plot ARI values

print(ggplot(ARIdf, aes(x = resolution, y = value, fill = variable)) + geom_line(aes(group=variable,color=variable), size = 1.2) + labs(y=NULL) + labs(x = "resolution") + ggtitle("Adjusted Rand Index") + facet_grid(variable ~ . , scales = "free") + theme_bw() + geom_vline(xintercept=optimalRes,color="orange",linetype="dotted"))

