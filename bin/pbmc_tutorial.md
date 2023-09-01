# Single-cell pipeline

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

```{r load libraries}
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
```


## Setup the Seurat Object

For this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are ~10,000 single cells that were sequenced on the Illumina NextSeq 500. 

The raw data can be found [here](https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.tar.gz").

```{r eval = FALSE} 
#download.file("https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k/sc5p_v2_hs_PBMC_10k_filtered_feature_bc_matrix.tar.gz", "pbmcs10k.tar.gz") 

#download.file("https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_melanoma_10k/sc5p_v2_hs_melanoma_10k_filtered_feature_bc_matrix.tar.gz", "melanoma.tar.gz") 

#system("tar -xzvf pbmcs10k.tar.gz")
#system("tar -xzvf melanoma.tar.gz")
``` 

### Reading in the data 

The `Read10X()` function reads in the output of the [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).


```{r init}
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "data/pbmcs10k/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "pbmc10k", min.cells = 3, min.features = 200)
pbmc
```


### Checking the count matrix

What does data in a count matrix look like?

```{r}
# Lets examine a few genes in the first thirty cells
pbmc.data$`Gene Expression`[c("CD3D","TCL1A","MS4A1"), 1:30]
```

The `.` values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0,  Seurat uses a sparse-matrix representation whenever possible. This results in significant memory and speed savings for Drop-seq/inDrop/10x data.



## Standard pre-processing workflow

The next steps encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 

### QC and selecting cells for further analysis

Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics [commonly used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community include

* The number of unique genes detected in each cell. 
    + Low-quality cells or empty droplets will often have very few genes
    + Cell doublets or multiplets may exhibit an aberrantly high gene count
* Total number of molecules detected within a cell (correlates strongly with unique genes)
* The percentage of reads that map to the mitochondrial genome
    + Low-quality / dying cells often exhibit extensive mitochondrial contamination
    + We calculate mitochondrial QC metrics with the `PercentageFeatureSet()` function, which calculates the percentage of counts originating from a set of features
    + We use the set of all genes starting with `MT-` as a set of mitochondrial genes


```{r mito}
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```


Where are QC metrics stored in Seurat?

* The number of unique genes and total molecules are automatically calculated during `CreateSeuratObject()`
    + You can find them stored in the object meta data
```{r qc, fig.height=7, fig.width=13}
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
```


In the example below, we visualize QC metrics, and use these to filter cells.

* We filter cells that have unique feature counts over 2,500 or less than 200
* We filter cells that have >5% mitochondrial counts
    
```{r qc2, fig.height=7, fig.width=13}
#Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt") #+geom_hline(yintercept=0.2)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") # +geom_hline(yintercept=1000)
plot1 + plot2


```
```{r qc3, fig.height=7, fig.width=13}
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
```


### Normalizing the data

After removing unwanted cells from the dataset, the next step is to normalize the data. 
By default, we employ a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in `pbmc[["RNA"]]@data`. These values are provided as default. 

```{r normalize.default, eval = FALSE}
pbmc <- NormalizeData(pbmc)
```


Another approach is using SCTransform. 


### Identification of highly variable features (feature selection)

We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). 
Focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets [see here](https://www.nature.com/articles/nmeth.2645).

The procedure in Seurat is described in detail [here](https://doi.org/10.1016/j.cell.2019.05.031), and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the `FindVariableFeatures()` function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.


```{r var_features, fig.height=5, fig.width=11}
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
```

### Scaling the data

Next, we apply a linear transformation ('scaling') that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The `ScaleData()` function:

* Shifts the expression of each gene, so that the mean expression across cells is 0
* Scales the expression of each gene, so that the variance across cells is 1
    + This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
* The results of this are stored in `pbmc[["RNA"]]@scale.data`

```{r regress, fig.height=7, fig.width=11, results='hide'}
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```



### Perform linear dimensional reduction

Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using `features` argument if you wish to choose a different subset.

```{r pca,results='hide'}
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

Seurat provides several useful ways of visualizing both cells and features that define the PCA, including `VizDimReduction()`, `DimPlot()`, and `DimHeatmap()`

```{r pca_viz, message=TRUE}
# Examine and visualize PCA results a few different ways
print(pbmc[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')
```

In particular `DimHeatmap()` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and features are ordered according to their PCA scores. Setting `cells` to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.

```{r single-heatmap}
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
```

```{r multi-heatmap, fig.height=15, fig.width=9}
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
```


## Clustering the cells

Seurat applies a graph-based clustering approach, building upon initial strategies in ([Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8)).
Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'. 

We first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). This step is performed using the `FindNeighbors()` function, and takes as input the previously defined dimensionality of the dataset (first 10 PCs).


To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [[SLM, Blondel *et al*., Journal of Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), to iteratively group cells together, with the goal of optimizing the standard modularity function. The `FindClusters()` function implements this procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters can be found using the `Idents()` function.

```{r cluster, fig.height=5, fig.width=7}
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
```


## Run non-linear dimensional reduction (UMAP/tSNE)

Seurat offers several non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets. The goal of these algorithms is to learn the underlying manifold of the data in order to place similar cells together in low-dimensional space. Cells within the graph-based clusters determined above should co-localize on these dimension reduction plots. As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis.

```{r umap, fig.height=5, fig.width=7}
pbmc <- RunUMAP(pbmc, dims = 1:10)
```

```{r umapplot, fig.height=5, fig.width=7}
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = 'umap')
```

```{r tsne, fig.height=5, fig.width=7}
pbmc <- RunTSNE(pbmc, dims = 1:10)
```

```{r tsneplot, fig.height=5, fig.width=7}
DimPlot(pbmc, reduction = 'tsne')
```
```{r markers - single} 
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
```
```{r markers - all } 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
```

### Visualizations
```{r markerplots, fig.height=10, fig.width=15}
# Violin plot of normalised data 
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# You can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = 'counts', log = TRUE)

# tSNE for each selected feature 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), reduction = "tsne")

# UMAP for each selected feature
FeaturePlot(pbmc, features = c("IL7R", "CCR7"), reduction = "tsne") 

```

```{r clusterHeatmap, fig.height=8, fig.width=15}
# Pick out top 10 markers per cluster 
pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10
# Plot heatmap 
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
```

```{r dotplot, fig.height=10, fig.width=8}
# Remove duplicate genes from top 10 markers per cluster list
goi = unique(top10$gene)
# PLot dotplot, flipping the x and y axes 
DotPlot(pbmc, features = goi ) + coord_flip()
```

### Cell annotation - manual 

```{r manual}
new.cluster.ids <- c("CD4 T", "Monocytes", "CD4 T", "CD8 T", "B", "NK", "Monocytes", "CD8 T", "CD8 T", "Platelets")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = 'tsne', label = TRUE, pt.size = 0.5) + NoLegend()
```

### Cell annotation - automatic
```{r save, eval = FALSE}
# Save object in format that Azimuth accepts. Here we are using .rds
saveRDS(pbmc, file="pbmcs.rds")
# Upload this file to the webtool
```

```{r cell annotation - azimuth online}
# Once completed, download the four files. Then import them into your object. 
imputed.assay <- readRDS('data/pbmcs10k/azimuth/azimuth_impADT.Rds')
pbmc <- pbmc[, Cells(imputed.assay)]
pbmc[['impADT']] <- imputed.assay

projected.umap <- readRDS('data/pbmcs10k/azimuth/azimuth_umap.Rds')
pbmc <- pbmc[, Cells(projected.umap)]
pbmc[['umap.proj']] <- projected.umap

predictions <- read.delim('data/pbmcs10k/azimuth/azimuth_pred.tsv', row.names = 1)
pbmc <- AddMetaData(object = pbmc, metadata = predictions)

```

```{r plot}
DimPlot(pbmc, reduction = 'tsne', label = TRUE, pt.size = 0.5, group.by = "predicted.celltype.l2")  
```

