set.seed(123456789)

# Get the UMAP 'pseudotime'

B_cell_names=c("CD34+","B cells","CD20+ B cells")
T_cell_names=c("CD34+","Immature Hematopoietic", "T cells")
tmp.seurat=SubsetData(bmmc.integrated,cells.use=rownames(subset(bmmc.integrated@meta.data,bmmc.integrated@ident %in% B_cell_names)))
print(DimPlot(object = tmp.seurat, cols.use=c(cols[4],cols[6],cols[8]),reduction.use = 'umap', do.label = T, pt.size = 0.5, do.return=T) + ggtitle(paste("CD34 + B cells, n=",nrow(tmp.seurat@meta.data),sep="")))
print(FeaturePlot(tmp.seurat, "CD34", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))
print(FeaturePlot(tmp.seurat, "VPREB1", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))
print(FeaturePlot(tmp.seurat, "MS4A1", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))

# Plot horizontaly

umapDF=data.frame(GetDimReduction(tmp.seurat,reduction.type="umap",slot="cell.embeddings"))
umapDF=subset(umapDF,umapDF$UMAP2>-5)
umapDF=cbind(umapDF,rescale(-umapDF$UMAP1),rescale(umapDF$UMAP2))
colnames(umapDF)[3:4]=c("norm_UMAP1","norm_UMAP2")
lo <- loess(umapDF$norm_UMAP1~umapDF$norm_UMAP2)
colors <- c(cols[4],cols[6],cols[8])
names(colors)=c("B cells","CD20+ B cells","CD34+")
plot(umapDF$norm_UMAP2,umapDF$norm_UMAP1,col=colors[tmp.seurat@ident[rownames(umapDF)]], xlab="Scaled UMAP2", ylab="Scaled UMAP1", pch=18)
tmpDF=cbind(umapDF$norm_UMAP2,predict(lo))
tmpDF=tmpDF[order(tmpDF[,1]),]
lines(tmpDF[,2]~tmpDF[,1], col="blue", lwd=2)
legend("topleft",c("CD34+","B cells","CD20+ B cells"),col=c(cols[8],cols[4],cols[6]), pch=15)

# Calculate distance along line fit to get the 'pseudotime'

colnames(tmpDF)=c("x","y")
pts <- st_as_sf(umapDF,coords=c("norm_UMAP2","norm_UMAP1")) # points
newline <- st_linestring(tmpDF[,1:2])
distances=st_distance(pts, newline)
plot_sf(pts)
plot(pts, add = TRUE)
plot(newline, add = TRUE)

# Remove points that are > 3*sd distance to line

pts <- st_as_sf(umapDF[distances<3*sd(distances),],coords=c("norm_UMAP2","norm_UMAP1"))
tmpDF=tmpDF[distances<3*sd(distances),]
tmpCol=colors[tmp.seurat@ident[rownames(umapDF)]][distances<3*sd(distances)]
newline=tmpDF[,1:2]
spline <- as(st_as_sfc(st_as_text(st_linestring(newline))), "Spatial")
position <- gProject(spline, as(pts, "Spatial"))
position <-  coordinates(gInterpolate(spline, position))
colnames(position) <- c("X2", "Y2")
segments <- data.frame(st_coordinates(pts), position)

print(ggplot() + geom_point(data = segments, aes(X,Y),col=tmpCol) + geom_line(data = as.data.frame(newline), aes(x,y)) + geom_segment(data = segments, aes(x = X, y = Y, xend = X2, yend = Y2), color = "orangered", alpha = 0.5) + ylab("Scaled UMAP1") + xlab("Scaled UMAP2") + ylim(0,1))

# Reassign new 'y' coordinates (UMAP1) on the fit line. Then get 'x' distance on the line (reassigned UMAP2)

print(ggplot() + geom_point(data = segments, aes(X,Y),col=tmpCol) + geom_line(data = as.data.frame(newline), aes(x,y)) + geom_segment(data = segments, aes(x = X, y = Y, xend = X2, yend = Y2), color = "orangered", alpha = 0.5) + geom_point(data=as.data.frame(position), aes(x=X2,y=Y2),color=tmpCol) + ylab("Scaled UMAP1") + xlab("Scaled UMAP2") + ylim(0,1)) 
pseudotimeDF=cbind(umapDF[distances<3*sd(distances),],reassigned_X=position[,1],reassigned_Y=position[,2])
spatialPointsObject <- SpatialPoints(coords = cbind(pseudotimeDF$reassigned_X,pseudotimeDF$reassigned_Y))
spatialLinesObject <- SpatialLines(LinesList = list(Lines(slinelist = list(Line(coords = cbind(tmpDF[,1],tmpDF[,2]))), ID = "1")))

# Get length to point and rescale to 0-1

d <- gProject(spatialLinesObject, spatialPointsObject, normalized=FALSE)
pseudotimeDF=cbind(pseudotimeDF,d)
colnames(pseudotimeDF)[ncol(pseudotimeDF)]="pseudotime"
pseudotimeDF$pseudotime=rescale(pseudotimeDF$pseudotime)
pseudotimeDF=merge(pseudotimeDF,tmp.seurat@ident,by='row.names')
colnames(pseudotimeDF)[1]="cell_names"
colnames(pseudotimeDF)[ncol(pseudotimeDF)]="cell_type"
rownames(pseudotimeDF)=pseudotimeDF$cell_names
pseudotimeDF=pseudotimeDF[order(pseudotimeDF$norm_UMAP2),]
write.table(pseudotimeDF, "pseudotime_bcells.tsv", quote=F, sep="\t", row.names=F)
plot(pseudotimeDF$norm_UMAP2,pseudotimeDF$pseudotime, xlab="Scaled UMAP2",ylab="Pseudotime", pch=20, main="Scaled UMAP2 vs pseudotime")
tmpCol=c(cols[4],cols[8],cols[6])
for(i in 1:length(unique(pseudotimeDF$cell_type))) {
	celltype=unique(pseudotimeDF$cell_type)[i]
	if(i==1) {
	plot(density(subset(pseudotimeDF,pseudotimeDF$cell_type==celltype)$pseudotime),col=tmpCol[i],xlim=c(0,1),ylim=c(0,10),xlab="Pseudotime",main="Pseudotime per cell type",lwd=2)
	}
	else {
		lines(density(subset(pseudotimeDF,pseudotimeDF$cell_type==celltype)$pseudotime),col=tmpCol[i],lwd=2)
	}
}
legend("top", c("CD34+","B cells","CD20+ B cells"), col=c(cols[4],cols[8],cols[6]),pch=15)

# Same classifier code for T cells

tmp.seurat=SubsetData(bmmc.integrated,cells.use=rownames(subset(bmmc.integrated@meta.data,bmmc.integrated@ident %in% T_cell_names)))
print(DimPlot(object = tmp.seurat, cols.use=c(cols[8],cols[7],cols[1]),reduction.use = 'umap', do.label = T, pt.size = 0.5, do.return=T) + ggtitle(paste("CD34 + T cells, n=",nrow(tmp.seurat@meta.data),sep="")))
print(FeaturePlot(tmp.seurat, "CD34", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))
print(FeaturePlot(tmp.seurat, "CD3D", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))
print(FeaturePlot(tmp.seurat, "CD4", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))
print(FeaturePlot(tmp.seurat, "CD8A", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap"))

# Plot horizontaly

umapDF=data.frame(GetDimReduction(tmp.seurat,reduction.type="umap",slot="cell.embeddings"))
umapDF=subset(umapDF,umapDF$UMAP2>-5 & umapDF$UMAP2<10)
umapDF=cbind(umapDF,rescale(umapDF$UMAP1),rescale(umapDF$UMAP2))
colnames(umapDF)[3:4]=c("norm_UMAP1","norm_UMAP2")
lo <- loess(umapDF$norm_UMAP1~umapDF$norm_UMAP2)
colors <- rev(c(cols[1],cols[7],cols[8]))
names(colors)=c("CD34","Immature Hematopoietic","T cells")
plot(umapDF$norm_UMAP2,umapDF$norm_UMAP1,col=colors[tmp.seurat@ident[rownames(umapDF)]], xlab="Scaled UMAP2", ylab="Scaled UMAP1", pch=18)
tmpDF=cbind(umapDF$norm_UMAP2,predict(lo))
tmpDF=tmpDF[order(tmpDF[,1]),]
lines(tmpDF[,2]~tmpDF[,1], col="blue", lwd=2)
legend("bottomright",c(names(colors)),col=c(cols[8],cols[7],cols[1]), pch=15)

# Calculate distance along line fit to get the 'pseudotime'

colnames(tmpDF)=c("x","y")
pts <- st_as_sf(umapDF,coords=c("norm_UMAP2","norm_UMAP1")) # points
newline <- st_linestring(tmpDF[,1:2])
distances=st_distance(pts, newline)
plot_sf(pts)
plot(pts, add = TRUE)
plot(newline, add = TRUE)

# Remove points that are > 3*sd distance to line

pts <- st_as_sf(umapDF[distances<3*sd(distances),],coords=c("norm_UMAP2","norm_UMAP1"))
tmpDF=tmpDF[distances<3*sd(distances),]
tmpCol=colors[tmp.seurat@ident[rownames(umapDF)]][distances<3*sd(distances)]

newline=tmpDF[,1:2]
spline <- as(st_as_sfc(st_as_text(st_linestring(newline))), "Spatial")
position <- gProject(spline, as(pts, "Spatial"))
position <-  coordinates(gInterpolate(spline, position))
colnames(position) <- c("X2", "Y2")
segments <- data.frame(st_coordinates(pts), position)

print(ggplot() + geom_point(data = segments, aes(X,Y),col=tmpCol) + geom_line(data = as.data.frame(newline), aes(x,y)) + geom_segment(data = segments, aes(x = X, y = Y, xend = X2, yend = Y2), color = "orangered", alpha = 0.5) + ylab("Scaled UMAP1") + xlab("Scaled UMAP2") + ylim(0,1))

# Reassign new 'y' coordinates (UMAP1) on the fit line. Then get 'x' distance on the line (reassigned UMAP2)

print(ggplot() + geom_point(data = segments, aes(X,Y),col=tmpCol) + geom_line(data = as.data.frame(newline), aes(x,y)) + geom_segment(data = segments, aes(x = X, y = Y, xend = X2, yend = Y2), color = "orangered", alpha = 0.5) + geom_point(data=as.data.frame(position), aes(x=X2,y=Y2),color=tmpCol) + ylab("Scaled UMAP1") + xlab("Scaled UMAP2") + ylim(0,1))

pseudotimeDF=cbind(umapDF[distances<3*sd(distances),],reassigned_X=position[,1],reassigned_Y=position[,2])
spatialPointsObject <- SpatialPoints(coords = cbind(pseudotimeDF$reassigned_X,pseudotimeDF$reassigned_Y))
spatialLinesObject <- SpatialLines(LinesList = list(Lines(slinelist = list(Line(coords = cbind(tmpDF[,1],tmpDF[,2]))), ID = "1")))

# Get length to point

d <- gProject(spatialLinesObject, spatialPointsObject, normalized=FALSE)
pseudotimeDF=cbind(pseudotimeDF,d)
colnames(pseudotimeDF)[ncol(pseudotimeDF)]="pseudotime"
pseudotimeDF$pseudotime=rescale(pseudotimeDF$pseudotime)
pseudotimeDF=merge(pseudotimeDF,tmp.seurat@ident,by='row.names')
colnames(pseudotimeDF)[1]="cell_names"
colnames(pseudotimeDF)[ncol(pseudotimeDF)]="cell_type"
rownames(pseudotimeDF)=pseudotimeDF$cell_names
pseudotimeDF=pseudotimeDF[order(pseudotimeDF$norm_UMAP2),]
write.table(pseudotimeDF, "pseudotime_tcells.tsv", quote=F, sep="\t", row.names=F)
plot(pseudotimeDF$norm_UMAP2,pseudotimeDF$pseudotime, xlab="Scaled UMAP2",ylab="Pseudotime", pch=20, main="Scaled UMAP2 vs pseudotime")
tmpCol=c(cols[1],cols[8],cols[7])
for(i in 1:length(unique(pseudotimeDF$cell_type))) {
        celltype=unique(pseudotimeDF$cell_type)[i]
        if(i==1) {
        plot(density(subset(pseudotimeDF,pseudotimeDF$cell_type==celltype)$pseudotime),col=tmpCol[i],xlim=c(0,1),ylim=c(0,20),xlab="Pseudotime",main="Pseudotime per cell type",lwd=2)
        }
        else {
                lines(density(subset(pseudotimeDF,pseudotimeDF$cell_type==celltype)$pseudotime),col=tmpCol[i],lwd=2)
        }
}
legend("topright",c(names(colors)),col=c(cols[8],cols[7],cols[1]), pch=15)
