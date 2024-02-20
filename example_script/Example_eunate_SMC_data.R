library(Seurat)
require(scales)
require(ggplot2)
require(viridis)
require(dplyr)

##### Source the file of predefined functions
source("/data/rq25/DA_project_boris_method/FUNs/FUNs_v0.R", echo = F)

##### Import the example SMC data and take a preliminary visualization
data_S <- readRDS("/data/rq25/DA_project_boris_method/data/data_S_SMC_eunate.rds")
DimPlot(data_S, group.by = "cluster", shuffle = T)
DimPlot(data_S, group.by = "time", shuffle = T)

##### Partition the data based on predefined categories, e.g., cluster labels, cell types...
partition_object_list <- get_partition(data_S, 
                                       partition.by = "cluster",
                                       reduction = "pca", 
                                       dims = 1:30,
                                       K = 25)

##### Perform local two-sample tests based on Boris's method (temporary acronym: RWSS)
RWSS_res_list <- perform_RWSS_test(data_S,
                                   partition_object_list,
                                   q = 0.01, 
                                   condition = "time")                                      


##### Define DA scores - Find the largest local density (normalized) deviation per cell
data_S <- get_DA_score(data_S, RWSS_res_list)

FeaturePlot(data_S, features = "DA_score") & 
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10)), limits = c(-1,1) * max(abs(data_S$DA_score)))


##### Define DA cells
data_S <- get_DA_cell(data_S,
                      RWSS_res_list,
                      cutoff = 0.05)


##### Visualize DA cells
DimPlot(data_S, group.by = "DA", cols = c(hue_pal()(2), "lightgray", "black"))


##### DA clustering
data_S <- get_DA_cluster(data_S,
                         resolution = 0.1,
                         partition.by = "cluster")     

DimPlot(data_S, group.by = "DA_cluster", order = T, reduction = "tsne", 
        label = T, label.size = 5, cols = hue_pal()(length(unique(data_S$DA_cluster))), na.value = "lightgray")


##### Find DA_cluster-specific markers
DefaultAssay(data_S) <- "RNA"
data_S$DA_cluster <- as.factor(data_S$DA_cluster)
Idents(data_S) <- "DA_cluster"
cluster_markers <- FindAllMarkers(
  data_S, only.pos = T, min.diff.pct = 0.1, verbose = T
)
cluster_markers <- cluster_markers[which(cluster_markers$p_val_adj < 0.05), ]

##### Visualization by heatmaps
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(data_S[, !is.na(data_S$DA_cluster)], features = top_markers$gene, slot = "scale.data") + 
  theme(axis.text.y = element_text(size = 8)) + scale_fill_viridis()
