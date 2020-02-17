library(infercnv)
set.seed(123456789)

# Create object

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="data/metacell.rawcounts.tsv", annotations_file="data/metacell.metadata.tsv", delim="\t", gene_order_file="data/genomic_position_nodups_chr.txt", ref_group_names=c("control"))

# Run

infercnv_obj = infercnv::run(infercnv_obj, cutoff=1, out_dir="output_dir_metacell", cluster_by_groups=T, denoise=T,HMM=F,window_length = 151)
