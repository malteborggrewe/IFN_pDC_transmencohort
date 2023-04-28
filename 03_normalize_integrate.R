###############################################################################
#' Malte Borggrewe
#' Normalize + scale + regress + integrate to remove donor variation
###############################################################################

# Based on seurat integration tutorial
# https://satijalab.org/seurat/articles/integration_introduction.html

source("src/00_library.R")
seu <- readRDS(paste0(dirs["rds"],"seu.RDS"))

# split the dataset into a list of 3 seurat objects (1 for each donor)
seu.list <- SplitObject(seu, split.by = "sample")

# normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)
  x <- FindVariableFeatures(x,
                            selection.method = "vst",
                            nfeatures = seurat[["variable_features"]])
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seu.list)

# Perform integration
seu.anchors <- FindIntegrationAnchors(object.list = seu.list,
                                      anchor.features = features)

# this command creates an 'integrated' data assay
seu.combined <- IntegrateData(anchorset = seu.anchors)

# make integrated as default assay
DefaultAssay(seu.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seu.combined <- ScaleData(seu.combined,
                          vars.to.regress = c("percent.mt", "nFeature_RNA"),
                          verbose = TRUE)

# Save RDS
saveRDS(seu, paste0(dirs["rds"],"seu_corrected.RDS"))

