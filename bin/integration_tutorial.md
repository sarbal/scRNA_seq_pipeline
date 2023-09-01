# Single-cell pipeline - integration and case versus control

## Loading libraries and packages

```{r setup, include=FALSE}
all_times <- list()  # store the time for each chunk
knitr::knit_hooks$set(time_it = local({
  now <- NULL
  function(before, options) {
    if (before) {
      now <<- Sys.time()
    } else {
      res <- difftime(Sys.time(), now, units = "secs")
      all_times[[options$label]] <<- res
    }
  }
}))
knitr::opts_chunk$set(
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 95),
  message = FALSE,
  warning = FALSE,
  time_it = TRUE
)

#!/usr/bin/env Rscript

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
  stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}
 
```

```{r}
library(Seurat)
library(dplyr)
library(patchwork) 
library(ggplot2)
library(speckle)
```

 

## Setup the Seurat Objects

```{r init}
# Load the patient datasets
pbmc.data <- Read10X(data.dir = "data/patient6/baseline/")
pbmc_base_p6 <- CreateSeuratObject(counts = pbmc.data, project = "baseline-p6", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir =  "data/patient6/D7/") 
pbmc_d7_p6 <- CreateSeuratObject(counts = pbmc.data, project = "D7-p6", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir =  "data/patient8/baseline/") 
pbmc_base_p8 <- CreateSeuratObject(counts = pbmc.data, project = "baseline-p8", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir =  "data/patient8/D7/") 
pbmc_d7_p8 <- CreateSeuratObject(counts = pbmc.data, project = "D7-p8", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir =  "data/patient9/baseline/") 
pbmc_base_p9 <- CreateSeuratObject(counts = pbmc.data, project = "baseline-p9", min.cells = 3, min.features = 200)

pbmc.data <- Read10X(data.dir =  "data/patient9/D7/") 
pbmc_d7_p9 <- CreateSeuratObject(counts = pbmc.data, project = "D7-p9", min.cells = 3, min.features = 200)

rm(pbmc.data)
```  

```{r mito, fig.height=7, fig.width=13}
pbmc.list = list(pbmc_base_p6, pbmc_d7_p6, pbmc_base_p8, pbmc_d7_p8, pbmc_base_p9, pbmc_d7_p9)
names(pbmc.list) = c("Baseline-p6", "D7-p6","Baseline-p8", "D7-p8","Baseline-p9", "D7-p9")
for( i in 1:length(pbmc.list)){ 
  pbmc.list[[i]][["percent.mt"]] <- PercentageFeatureSet(pbmc.list[[i]], pattern = "^MT-")
}
```
```{r }
#Visualize QC metrics as a violin plot
VlnPlot( pbmc.list[[1]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot( pbmc.list[["Baseline-p9"]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
  
```{r subsetting} 
for( i in 1:length(pbmc.list)){ 
  pbmc.list[[i]] <- subset(pbmc.list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 7)
}
```


```{r prepare for integration}
pbmc.list <- lapply(X = pbmc.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = features)
``` 
```{r }
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", anchor.features = features)
```

```{r integrate data}
pbmc.combined.sct <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")
```

```{r }
pbmc.combined.sct <- RunPCA(pbmc.combined.sct, verbose = FALSE)
pbmc.combined.sct <- RunUMAP(pbmc.combined.sct, reduction = "pca", dims = 1:30)
#pbmc.combined.sct <- RunTSNE(pbmc.combined.sct, reduction = "pca", dims = 1:30)
pbmc.combined.sct <- FindNeighbors(pbmc.combined.sct, dims = 1:10)
pbmc.combined.sct <- FindClusters(pbmc.combined.sct, resolution = 0.5)
```
```{r }
DimPlot(pbmc.combined.sct, reduction = "umap")
DimPlot(pbmc.combined.sct, reduction = "umap", group.by = "orig.ident")


pbmc.combined.sct[["patients"]] = factor(gsub("baseline-|D7-", "", pbmc.combined.sct$orig.ident ))
pbmc.combined.sct[["type"]] = factor(gsub("-p\\d", "", pbmc.combined.sct$orig.ident ))


DimPlot(object = pbmc.combined.sct, reduction = "umap", group.by = "patients")
DimPlot(object = pbmc.combined.sct, reduction = "umap", group.by = "type")

```
 

```{r }
#imputed.assay <- readRDS('data/patient9/azimuth/azimuth_impADT.Rds')
#pbmc.combined.sct <- pbmc.combined.sct[, Cells(imputed.assay)]
#pbmc.combined.sct[['impADT']] <- imputed.assay


predictions <- read.delim('data/patient9/azimuth/azimuth_pred.tsv', row.names = 1)
pbmc.combined.sct <- AddMetaData(object = pbmc.combined.sct, metadata = predictions)


#projected.umap <- readRDS('data/patient9/azimuth/azimuth_umap.Rds')
#pbmc.combined.sct <- pbmc.combined.sct[, Cells(projected.umap)]
#pbmc.combined.sct[['umap.proj']] <- projected.umap
```


```{r fig.width=10, fig.height=7}
DimPlot(pbmc.combined.sct, reduction = "umap", group.by = "predicted.celltype.l2")
DimPlot(pbmc.combined.sct, reduction = "umap", group.by = "predicted.celltype.l1")
```

```{r calculate proportions}
freq = plyr::count(cbind(pbmc.combined.sct$predicted.celltype.l1, as.character(pbmc.combined.sct$type) ))
freq_table <- tidyr::spread(freq, key=1, value=3) 
freq_table[is.na(freq_table)] = 0
frac_table = t((freq_table[,-1])/rowSums(freq_table[,-1] , na.rm=T)) * 100 
colnames(frac_table) = c("baseline", "D7") 
my_color_palette <- scales::hue_pal()(8)
barplot( (frac_table) , col=my_color_palette )

```

```{r proportions per sample}
freq = plyr::count(cbind(as.character(pbmc.combined.sct$predicted.celltype.l1), as.character(pbmc.combined.sct$orig.ident) ))
freq_table <- tidyr::spread(freq, key=1, value=3) 
freq_table[is.na(freq_table)] = 0

frac_table = t((freq_table[,-1])/rowSums(freq_table[,-1] , na.rm=T)) * 100

colnames(frac_table) = freq_table[,1]
my_color_palette <- scales::hue_pal()(8)
barplot( frac_table , col=my_color_palette )

```

## Run either t-test or ANOVA to calculate significance 
```{r test for signficance}
pbmc.combined.sct$predicted.celltype.l1 <- factor(pbmc.combined.sct$predicted.celltype.l1)
props =  propeller(clusters=pbmc.combined.sct$predicted.celltype.l1, sample=pbmc.combined.sct$orig.ident, group=pbmc.combined.sct$type)
props
```
## Prepare data
Since we ran SCTransform, we need to recorrect our data before running the FindMarkers() function. Here we are using the default wilcox-test for differential expression but there are many other approaches that can be used. 
```{r find markers/DE, fig.height=8, fig.width=15}
pbmc.combined.sct = PrepSCTFindMarkers(pbmc.combined.sct)
filt = pbmc.combined.sct@meta.data$predicted.celltype.l1 == "CD8 T"
 
b_v_d7.markers <- FindMarkers(pbmc.combined.sct[,filt], assay="SCT", ident.1 = "D7",  ident.2 = "baseline", group.by = "type", min.pct = 0.25, logfc.threshold = 0, recorrect_umi = FALSE)
 
```

```{r volcano plot}
plot(b_v_d7.markers$avg_log2FC, -log10(b_v_d7.markers$p_val) ,xlab="Average log2FC",  ylab="-log10(p-value)", pch=19, main = "baseline vs D7") 
abline(h=-log10(0.05/dim(b_v_d7.markers)[1]), lwd=2, col=4, lty=2)
abline(v=c(1,-1), lty=3, lwd=2, col=4)
abline(v=c(0.5,-0.5), lty=2, col=4)

head( b_v_d7.markers[  abs(b_v_d7.markers$avg_log2FC) > 1 & b_v_d7.markers$p_val_adj < 0.05, ] ) 
```


## No integration/batch correction 
If you are interested to see what the data would have looked like had we not integrated the data properly, we can just merge the Seurat objects and run the same analyses.   

```{r no batch correction, eval = FALSE}
pbmc.merged <- merge(pbmc_base_p6, y = c(pbmc_d7_p6, pbmc_base_p8, pbmc_d7_p8, pbmc_base_p9, pbmc_d7_p9), add.cell.ids = c("Baseline-p6", "D7-p6","Baseline-p8", "D7-p8","Baseline-p9", "D7-p9"), project = "PBMCs")

pbmc.merged <- NormalizeData(pbmc.merged)
pbmc.merged <- FindVariableFeatures(pbmc.merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc.merged)
pbmc.merged <- ScaleData(pbmc.merged, features = all.genes)
pbmc.merged <- RunPCA(pbmc.merged, verbose = FALSE)
pbmc.merged <- RunUMAP(pbmc.merged, reduction = "pca", dims = 1:30)
pbmc.merged <- RunTSNE(pbmc.merged, reduction = "pca", dims = 1:30)
```

```{r dim plot comparisons}
DimPlot(pbmc.merged, reduction = "umap")
DimPlot(object = pbmc.combined.sct, reduction = "umap", group.by = "orig.ident")
```
