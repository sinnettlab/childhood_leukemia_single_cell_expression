set.seed(123456789)

# Read 10X data

ROOT_DIR="."
ABMMC_1=Read10X(paste(ROOT_DIR,"/v2.1/frozen_bmmc_healthy_donor1/grch38/ABMMC_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
ABMMC_2=Read10X(paste(ROOT_DIR,"/v2.1/frozen_bmmc_healthy_donor2/grch38/ABMMC_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PBMMC_1=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_1/grch38/PBMMC_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PBMMC_2=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_2/grch38/PBMMC_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PBMMC_3=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_3/grch38/PBMMC_3/outs/filtered_gene_bc_matrices/GRCh38",sep=""))

objectList=list(ABMMC_1,ABMMC_2,PBMMC_1,PBMMC_2,PBMMC_3)
objectNames=c("ABMMC.1","ABMMC.2","PBMMC.1","PBMMC.2","PBMMC.3")

# For cell cycle estimation

cc.genes <- readLines(con = "regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]


# Keep all genes expressed in min cells, keep all cells with min genes

for(i in 1:length(objectList)) {
	colnames(x = objectList[[i]]) <- paste(objectNames[i], colnames(x = objectList[[i]]), sep = '_')
        objectList[[i]] <- CreateSeuratObject(raw.data = objectList[[i]], min.cells = 5, min.genes = 200, project = objectNames[i])
	mito.genes <- grep("^MT-", rownames(objectList[[i]]@data), value = T)
	objectList[[i]] <- AddMetaData(objectList[[i]], Matrix::colSums(expm1(objectList[[i]]@data[mito.genes, ]))/Matrix::colSums(expm1(objectList[[i]]@data)), "percent.mito")
	print(VlnPlot(object = objectList[[i]], features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3))
	objectList[[i]] <- FilterCells(objectList[[i]], subset.names = "nGene", low.thresholds = 200, high.thresholds=Inf)
	objectList[[i]] <- FilterCells(objectList[[i]], subset.names = "percent.mito", high.thresholds= 0.08)
	objectList[[i]] <- NormalizeData(objectList[[i]])
	objectList[[i]] <- FindVariableGenes(objectList[[i]], do.plot = T, display.progress = F)
	
	# Regress out cc, nUMI, mito.pt

	objectList[[i]] <- CellCycleScoring(object = objectList[[i]], s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
	objectList[[i]]@meta.data$CC.Difference <- objectList[[i]]@meta.data$S.Score - objectList[[i]]@meta.data$G2M.Score
	objectList[[i]] <- ScaleData(objectList[[i]], vars.to.regress = c("nUMI","percent.mito","S.Score", "G2M.Score"))
	print(VlnPlot(object = objectList[[i]], features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by="Phase"))
}	

# Genes to use for the analysis

genes.use <- c()
for (i in 1:length(objectList)) {
	  genes.use <- c(genes.use, head(rownames(objectList[[i]]@hvg.info), 1000))
}
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(objectList)) {
	  genes.use <- genes.use[genes.use %in% rownames(objectList[[i]]@scale.data)]
}

write.table(genes.use,"variable_genes_control_bmmcs.tsv", col.names=F, row.names=F,quote=F,sep="\t")

# Run multi-set CCA

nmbDims=10
bmmc.integrated <- RunMultiCCA(objectList, genes.use = genes.use, num.ccs = nmbDims)
DimPlot(object = bmmc.integrated , reduction.use = "cca", group.by = "orig.ident", pt.size = 0.5, do.return = FALSE)
VlnPlot(object = bmmc.integrated, features.plot = "CC1", group.by = "orig.ident", do.return = FALSE)
MetageneBicorPlot(bmmc.integrated, grouping.var = "orig.ident", dims.eval = 1:nmbDims, display.progress = TRUE)
DimHeatmap(object = bmmc.integrated, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)

# Run rare non-overlapping filtering

bmmc.integrated <- CalcVarExpRatio(object = bmmc.integrated, reduction.type = "pca", grouping.var = "orig.ident", dims.use = 1:nmbDims)
bmmc.integrated <- SubsetData(bmmc.integrated, subset.name = "var.ratio.pca", accept.low = 0.5)

# Alignment

bmmc.integrated <- AlignSubspace(bmmc.integrated, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:nmbDims)
metrics=CalcAlignmentMetric(bmmc.integrated,dims.use=1:nmbDims,grouping.var="orig.ident")
print(metrics)
write.table(metrics, "alignment_score.txt", quote=F, row.names=F, col.names=F)

# t-SNE and Clustering

resCluster=0.4
bmmc.integrated <- FindClusters(bmmc.integrated, reduction.type = "cca.aligned", dims.use = 1:nmbDims, force.recalc=TRUE, save.SNN = T, resolution = resCluster) 
bmmc.integrated <- RunTSNE(bmmc.integrated, reduction.use = "cca.aligned", dims.use = 1:nmbDims)
bmmc.integrated <- RunUMAP(bmmc.integrated, reduction.use = "cca.aligned", dims.use = 1:nmbDims)

# Visualization

bmmc.integrated <- SetAllIdent(bmmc.integrated, id = "orig.ident")
print(TSNEPlot(bmmc.integrated, do.return = TRUE, pt.size = 0.5))
DimPlot(object = bmmc.integrated, reduction.use = 'umap', do.label = F, pt.size = 0.5)
bmmc.integrated <- SetAllIdent(bmmc.integrated, id = "Phase")
print(TSNEPlot(bmmc.integrated,do.label=F, do.return=TRUE, pt.size = 0.5))
DimPlot(object = bmmc.integrated, reduction.use = 'umap', do.label = F, pt.size = 0.5)
bmmc.integrated <- SetAllIdent(bmmc.integrated, id = paste("res.",resCluster,sep=""))
print(TSNEPlot(bmmc.integrated,do.label=T, do.return=TRUE, pt.size = 0.5))
DimPlot(object = bmmc.integrated, reduction.use = 'umap', do.label = T, pt.size = 0.5)
FeaturePlot(bmmc.integrated, "nUMI", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
FeaturePlot(bmmc.integrated, "percent.mito", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
markerSet=c("MS4A1","CD19","MME","VPREB1","VPREB3","DNTT","CD79A","MZB1","NKG7","CD3D","CD8A","CD34","CST3","LYZ","HBA1","FCGR3A","GATA1","GATA2")

for(i in markerSet) {
		FeaturePlot(bmmc.integrated, i, cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
}

