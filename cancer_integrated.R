set.seed(123456789)

# Read 10X data

ROOT_DIR="."
PBMMC_1=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_1/grch38/PBMMC_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PBMMC_2=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_2/grch38/PBMMC_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PBMMC_3=Read10X(paste(ROOT_DIR,"/v2.1/PBMMC_3/grch38/PBMMC_3/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
ETV6_RUNX1_1=Read10X(paste(ROOT_DIR,"/v2.1/ETV6-RUNX1_1/grch38/ETV6_RUNX1_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
ETV6_RUNX1_2=Read10X(paste(ROOT_DIR,"/v2.1/ETV6-RUNX1_2/grch38/ETV6_RUNX1_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
ETV6_RUNX1_3=Read10X(paste(ROOT_DIR,"/v2.1/ETV6-RUNX1_3/grch38/ETV6_RUNX1_3/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
ETV6_RUNX1_4=Read10X(paste(ROOT_DIR,"/v2.1/ETV6-RUNX1_4/grch38/ETV6_RUNX1_4/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
HHD_1=Read10X(paste(ROOT_DIR,"/v2.1/HHD_1/grch38/HHD_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
HHD_2=Read10X(paste(ROOT_DIR,"/v2.1/HHD_2/grch38/HHD_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PRE_T_1=Read10X(paste(ROOT_DIR,"/v2.1/PRE-T_1/grch38/PRE-T_1/outs/filtered_gene_bc_matrices/GRCh38",sep=""))
PRE_T_2=Read10X(paste(ROOT_DIR,"/v2.1/PRE-T_2/grch38/PRE-T_2/outs/filtered_gene_bc_matrices/GRCh38",sep=""))

objectList=list(PBMMC_1,PBMMC_2,PBMMC_3,ETV6_RUNX1_1,ETV6_RUNX1_2,ETV6_RUNX1_3,ETV6_RUNX1_4,HHD_1,HHD_2,PRE_T_1,PRE_T_2)
objectNames=c("PBMMC.1","PBMMC.2","PBMMC.3","ETV6.RUNX1.1","ETV6.RUNX1.2","ETV6.RUNX1.3","ETV6.RUNX1.4","HHD.1","HHD.2","PRE-T.1","PRE-T.2")

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
}	

# Merge datasets 

cancer.integrated=objectList[[1]]
for(i in 2:(length(objectList))) {
		cancer.integrated=MergeSeurat(cancer.integrated,objectList[[i]], do.normalize=F, do.scale=F)
}
cancer.integrated <- NormalizeData(object = cancer.integrated)
cancer.integrated <- ScaleData(object = cancer.integrated)
cancer.integrated <- CellCycleScoring(object = cancer.integrated, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
cancer.integrated@meta.data$CC.Difference <- cancer.integrated@meta.data$S.Score - cancer.integrated@meta.data$G2M.Score
cancer.integrated <- FindVariableGenes(object = cancer.integrated, do.plot = FALSE)
cancer.integrated <- RunPCA(object = cancer.integrated, pc.genes = cancer.integrated@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

print(PCElbowPlot(cancer.integrated, num.pc = 20))
VizPCA(cancer.integrated, 1:9)
PCAPlot(cancer.integrated, 1, 2)

# tSNE / UMAP

nmbDims=20 # based on screeplot
resCluster=0.1
cancer.integrated <- FindClusters(cancer.integrated, dims.use = 1:nmbDims, force.recalc=TRUE, save.SNN = T, resolution = resCluster)
cancer.integrated <- RunTSNE(cancer.integrated, dims.use = 1:nmbDims)
cancer.integrated <- RunUMAP(cancer.integrated, dims.use = 1:nmbDims)

# Visualization

cancer.integrated <- SetAllIdent(cancer.integrated, id = "orig.ident")
print(TSNEPlot(cancer.integrated, do.return = TRUE, pt.size = 0.5))
DimPlot(object = cancer.integrated, reduction.use = 'umap', do.label = F, pt.size = 0.5)
cancer.integrated <- SetAllIdent(cancer.integrated, id = "Phase")
print(TSNEPlot(cancer.integrated,do.label=F, do.return=TRUE, pt.size = 0.5))
DimPlot(object = cancer.integrated, reduction.use = 'umap', do.label = F, pt.size = 0.5)
cancer.integrated <- SetAllIdent(cancer.integrated, id = paste("res.",resCluster,sep=""))
print(TSNEPlot(cancer.integrated,do.label=T, do.return=TRUE, pt.size = 0.5))
DimPlot(object = cancer.integrated, reduction.use = 'umap', do.label = T, pt.size = 0.5)
FeaturePlot(cancer.integrated, "nUMI", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
FeaturePlot(cancer.integrated, "percent.mito", cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
markerSet=c("MS4A1","CD19","MME","VPREB1","VPREB3","DNTT","CD79A","MZB1","NKG7","CD3D","CD8A","CD34","CST3","LYZ","HBA1","FCGR3A","EPCAM","GATA1")
for(i in markerSet) {
		FeaturePlot(cancer.integrated, i, cols.use = c("grey","blue"), no.legend=F, reduction.use="umap")
}


