rm(list=ls())
setwd("/Users/lihuazhang/Documents/Vitiligo_revise/revision")
library(devtools)
#devtools::install_github("amsszlh/scMC")
library(scMC)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
library(leiden)
# use data of original, lesional and nonlsional, current, healthy, active
# load current data
load("combined_integration_healthy_nonlesional_activeVitiligo.RData")
combined$conditions <- as.character(combined$Phenotype)
combined$patient <- combined$Patient

object_list1 <- SplitObject(object = combined, split.by = "Phenotype")
#object_list1a <- SplitObject(object = object_list1$Healthy, split.by = "Patient")
#object_list1b <- SplitObject(object = object_list1$Lesional, split.by = "Patient")
#object_list1c <- SplitObject(object = object_list1$Nonlesional, split.by = "Patient")
#names1a <- paste0("Healthy_",names(object_list1a))
#names1b <- paste0("Lesional_",names(object_list1b))
#names1c <- paste0("Nonlesional_",names(object_list1c))
rm(combined)
data <- readRDS(file = "Vitiligo_remove_G_normal_newApril.rds")
# generate split datasets
object_list2 <- SplitObject(object = data, split.by = "conditions")
#object_list2a <- SplitObject(object = object_list2$normal, split.by = "patient")
#object_list2b <- SplitObject(object = object_list2$vitiligo, split.by = "patient")
#names2a <- paste0("normal_",names(object_list2a))
#names2b <- paste0("vitiligo_",names(object_list2b))
#rm(data)

for (i in 1:length(object_list1)) {
  x <- object_list1[[i]]
  if (ncol(x@assays$RNA@counts) > 4000) {
    # randomly selected
    index <- sample(ncol(x@assays$RNA@counts),size = 4000)
    object_list1[[i]] <-subset(x,cells = colnames(x)[index])
  }
}

for (i in 1:length(object_list2)) {
  x <- object_list2[[i]]
  if (ncol(x@assays$RNA@counts) > 4000) {
    # randomly selected
    index <- sample(ncol(x@assays$RNA@counts),size = 4000)
    object_list2[[i]] <-subset(x,cells = colnames(x)[index])
  }
}

object.list <- c(object_list1,object_list2)

for (i in 1:length(object.list)) {
  #x <- NormalizeData(object.list[[i]], verbose = FALSE)
  x <- object.list[[i]]
  x <- FindVariableFeatures(x)
  # perform scaling on the previously identified variable features
  x <- ScaleData(x, verbose = FALSE)
  object.list[[i]] <- x
}

########### Part II: Perform an integrated analysis using scMC ###########
# step2. identify clusters with different resolution for each condition
# compute SNN
object.list <- identifyNeighbors(object.list)
# identify clusters
object.list <- identifyClusters(object.list, mode = "separate",resolution = 1,algorithm = 4)

## step3. detect cluster-specific cells with high confident
features.integration = identifyIntegrationFeatures(object.list, nfeatures = 5000)
object.list <- identifyConfidentCells(object.list, features.integration,quantile.cutoff = 0.75)

## step4. Identify marker genes associated with the putative cell clusters in each dataset
object.list <- identifyMarkers(object.list, test.use = "bimod")
structured_mat <- learnTechnicalVariation(object.list, features.integration, similarity.cutoff = 0.6)

## step 7. Learn a shared embedding of cells across all datasets after removing technical variation
combined <- merge(x = object.list[[1]],y = object.list[2:length(x = object.list)])
VariableFeatures(combined) <- features.integration
combined <- integrateData(combined, structured_mat, lambda = 10)

nPC = 40
combined  <- FindNeighbors(combined , reduction = "scMC")
combined  <- FindClusters(combined, algorithm = 1, resolution = 0.5)
levels(Idents(combined))
#combined.new <- BuildClusterTree(combined.new, reorder = T, reorder.numeric = T, verbose = F)
combined <- RunUMAP(combined, reduction='scMC', dims = 1:nPC)
gg <- DimPlot(combined, reduction = "umap", group.by = "ident")
gg

gg <- DimPlot(combined, reduction = "umap", split.by = "conditions",label = T, pt.size = 0.00000001)
gg

##### include other cells
load("object_list1o.RData")
load("object_list2o.RData")
data1_new <- list()
meta1_new <- list()
centerX = TRUE; scaleX = FALSE;
for (i in 1:length(object_list1)) {
  objo <- object_list1o[[i]]
  obji <- object_list1[[i]]
  index1 <- which(colnames(objo@assays$RNA@counts) %in% colnames(obji@assays$RNA@counts))
  index1 <- setdiff(1:ncol(objo@assays$RNA@counts),index1)
  data.use <- objo@assays$RNA@data
  data.use <- data.use[, index1, drop = FALSE]
  data.use <- scale(data.use, center = centerX, scale = scaleX)
  data1_new[[i]] <- data.use
  meta_data <- objo@meta.data
  meta1_new[[i]] <- meta_data[index1,]
}
data1_all <- do.call(cbind,data1_new)
meta1_all <- do.call(rbind,meta1_new)
# 
index <- match(features.integration,rownames(data1_all))
identical(features.integration,rownames(data1_all)[index])
data1_all <- data1_all[index,]

data2_new <- list()
meta2_new <- list()
centerX = TRUE; scaleX = FALSE;
for (i in 1:length(object_list2)) {
  objo <- object_list2o[[i]]
  obji <- object_list2[[i]]
  index1 <- which(colnames(objo@assays$RNA@counts) %in% colnames(obji@assays$RNA@counts))
  index1 <- setdiff(1:ncol(objo@assays$RNA@counts),index1)
  data.use <- objo@assays$RNA@data
  data.use <- data.use[, index1, drop = FALSE]
  data.use <- scale(data.use, center = centerX, scale = scaleX)
  data2_new[[i]] <- data.use
  meta_data <- objo@meta.data
  meta2_new[[i]] <- meta_data[index1,]
}

data2_all <- do.call(cbind,data2_new)
meta2_all <- do.call(rbind,meta2_new)
# 
index <- match(features.integration,rownames(data2_all))
identical(features.integration,rownames(data2_all)[index])
data2_all <- data2_all[index,]


data_new <- cbind(data1_all,data2_all)
colnames(meta1_all)
colnames(meta2_all)
colnames(combined@meta.data)
meta1_all <- meta1_all[,c(1,2,3,16,17,15)]
colnames(meta1_all) <- c("orig.ident","nCount_RNA","nFeature_RNA","conditions","patient","cell_type_prior")
meta2_all <- meta2_all[,c(1,2,3,5,6,18)]
colnames(meta2_all) <- c("orig.ident","nCount_RNA","nFeature_RNA","conditions","patient","cell_type_prior")

meta_new <- rbind(meta1_all,meta2_all)
source("projectData.R")
combined.new <- projectData(combined, data_new, meta.new = meta_new)


nPC = 40
combined.new  <- FindNeighbors(combined.new , reduction = "scMC")
combined.new  <- FindClusters(combined.new , algorithm = 1, resolution = 0.5)
levels(Idents(combined.new))

combined.8 <- subset(combined.new, idents = "8")
combined.8 <- FindNeighbors(combined.8, reduction = "scMC", dims = 1:nPC)
combined.8 <- FindClusters(combined.8, algorithm = 4, resolution = 0.2)
levels(Idents(combined.8))
DimPlot(combined.8, reduction = "umap", label = TRUE)
features <- c('CXCL9','CXCL10')
gg <- StackedVlnPlot(combined.8, features = features)
gg
new.cluster.ids1 <- c("8","13","14","15")
names(new.cluster.ids1) <- levels(combined.8)
combined.8 <- RenameIdents(combined.8, new.cluster.ids1)
combined.new$cluster <- as.character(Idents(combined.new))
combined.new$cluster[Cells(combined.8)] <- as.character(Idents(combined.8))
Idents(combined.new) <- combined.new$cluster

combined_o <- combined.new
new.order <- c("0","1","2","3","4", "5","6", "7","8","9","10","11", "12","13","14","15")
Idents(combined_o) <- factor(Idents(combined_o), levels = new.order)
combined.new <- combined_o
#combined.new <- BuildClusterTree(combined.new, reorder = T, reorder.numeric = T, verbose = F)
combined.new <- RunUMAP(combined.new, reduction='scMC', dims = 1:nPC)
gg <- DimPlot(combined.new, reduction = "umap", group.by = "ident")
gg
cowplot::save_plot(filename=paste0("integration_", "_umap_all_0.6.pdf"), plot=gg, base_width =4, base_height = 3)

gg <- DimPlot(combined.new, reduction = "umap", split.by = "conditions",label = T, pt.size = 0.00000001)
gg
cowplot::save_plot(filename=paste0("integration_", "_umap_split_all_0.6.pdf"), plot=gg, base_width =15, base_height = 3)

DimPlot(combined, reduction = "umap", group.by = "conditions")

DefaultAssay(combined) <- "RNA"

features = toupper(c("KRT6A","KRT16","KRT6B","CXCL9","CXCL10","KRT6C","S100A9","S100A8"))
features = toupper(c("KRT15","KRT5","KRT1","KRT2","FLG","LOR"))
gg <- FeaturePlot(combined.new, features = features[6],split.by = "conditions")
gg
cowplot::save_plot(filename=paste0("LOR_", "_umap_split_all_0.6.pdf"), plot=gg, base_width =15, base_height = 3)

##### according to original labels to annotated clusters
# build clusters of original cell types, not split to each condition

data_used <- combined.new@assays$RNA@data
cluster_orig <- rep(NA,ncol(data_used))
for (i in 1:length(object_list1o)) {
  objo <- object_list1o[[i]]
  metai <- objo@meta.data
  indexi <- match(rownames(metai),colnames(data_used))
  identical(rownames(metai),colnames(data_used)[indexi])
  cluster_orig[indexi] <- metai$Cluster_Refined
}
for (i in 1:length(object_list2o)) {
  objo <- object_list2o[[i]]
  metai <- objo@meta.data
  indexi <- match(rownames(metai),colnames(data_used))
  identical(rownames(metai),colnames(data_used)[indexi])
  cluster_orig[indexi] <- as.character(metai$clusters.final)
}
table(cluster_orig)
combined.new$cluster_orig <- cluster_orig
DimPlot(combined.new, reduction = "umap", group.by = "cluster_orig")
gg <- DimPlot(combined.new, reduction = "umap", group.by = "cluster_orig")
gg
cowplot::save_plot(filename=paste0("integration_", "_umap_all_0.6_type.pdf"), plot=gg, base_width =8, base_height = 3)

gg <- DimPlot(combined.new, reduction = "umap", group.by = "cluster_orig", split.by = "conditions",label = T, pt.size = 0.00000001)
gg
cowplot::save_plot(filename=paste0("integration_", "_umap_split_all_0.6_type.pdf"), plot=gg, base_width =15, base_height = 3)

############## compare markers between stress and acurate populations
save.image(file = "integrate_scMC_0.6_all.RData")
### compute the ratios of each cluster in identified new identies, plot the heatmap
# build the heatmap of each cluster
cluster <- Idents(combined.new)
c0 <- sort(unique(cluster))
id1 <- which(combined.new$conditions %in% c("Healthy","Lesional","Nonlesional"))
id2 <- which(combined.new$conditions %in% c("normal","vitiligo"))
type1 <- cluster_orig[id1]
c1 <- c("KRT-B1","KRT-B2","KRT-SP","KRT-GR","KRT-ECR", "MEL","DC","CD8","TConv","Treg","MAC","GD","NK")
type2 <- cluster_orig[id2]
c2 <- c("Basal 1" ,"Basal 2" ,"B2S 1","B2S 2","Spinous","S2G 1", "S2G 2", "Granular","Stress 1","Stress 2","Melanocytes","Cycling","DC","TC")
H1 <- matrix(0,nrow = length(unique(cluster)),ncol = length(c1))
H2 <- matrix(0,nrow = length(unique(cluster)),ncol = length(c2))
for (i in 1:length(unique(cluster))) {
  index <- which(cluster == c0[i])
  subc1_id <- intersect(index,id1)
  subc2_id <- intersect(index,id2)
  a1 <- cluster_orig[subc1_id]
  a2 <- cluster_orig[subc2_id]
  for (j in 1:length(c1)) {
    H1[i,j] <- length(which(a1 == c1[j]))
  }
  for (j in 1:length(c2)) {
    H2[i,j] <- length(which(a2 == c2[j]))
  }
}
colnames(H1) <- c1
colnames(H2) <- c2
write.table(H1,file = "cluster_group_1_0.6_all.txt",sep = "\t")
write.table(H2,file = "cluster_group_2_0.6_all.txt",sep = "\t")

write.table(c1,file = "cluster_name_1_0.6_all.txt",sep = "\t",row.names = FALSE)
write.table(c2,file = "cluster_name_2_0.6_all.txt",sep = "\t",row.names = FALSE)

features = toupper(c("KRT2","KRT1","KRT5","CYR61","KRT15","FLG"))
gg <- FeaturePlot(combined.new, features = features[2],split.by = "conditions")
gg

### compare heatmap of cluster 6,8,12
combined.new <- ScaleData(combined.new, features = rownames(combined.new))
combined.part <- subset(combined.new, idents = c("6","12","15"))
markers <- FindAllMarkers(object = combined.part, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, min.diff.pct = 0.1)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
features <- top10$gene
features[2] <- "CXCL10"
features <- c("KRT6A","KRT16","KRT6B","S100A8","S100A9","KRT17","KRT6C","CXCL9","CXCL10","BST2","IFI6","IFI27","MALAT1","ISG15")
colors1 <- c("#2DB567","#A781BA","#F069A0")
colors2 <- c("#CB972B","#01BBDA","#e14891","#328C58","#6A52A3")
gg <- doHeatmap(object = combined.part, features = features,additional.group.by = "conditions", raster = FALSE,colors.use = list(idents = colors1,conditions = colors2)) + NoLegend()
gg
cowplot::save_plot(filename = "heatmap_markers_integrate_all_0.6.pdf", plot=gg, base_width =5, base_height = 5)



