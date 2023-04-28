###############################################################################
#' Malte Borggrewe
#' Load all packages and define custom functions
###############################################################################

##### Load all packages -------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(clustree)
  library(MAST)
  library(VennDiagram)
  library(GSEABase)
  library(AUCell)
  library(pheatmap)
  library(EnhancedVolcano)
  library(enrichR)
  source("src/00_settings.R")
})


##### Small helpers -----------------------------------------------------------
# Helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}


##### Custom functions --------------------------------------------------------
#' Loads GEX and HTO data from H5 Matrix and demultiplexes based on HTO
#' Saves .csv with stats to qc folder and returns a ready-to-use Seurat object
#' Requires settings file (00_settings.R)
#' @param sample Name of the sample. A folder with the sample name needs to
#' exist in raw_dir folder.
#' @param raw_dir Folder with raw files containing h5 files (cellranger output)
#' @param output_dir Output dir for the QC files to be saved
load_demultiplex <- function(sample, raw_dir, output_dir) {
  # Read cell ranger output as H5 matrix (contains both GEX and HTO)
  gex_hto <- Read10X_h5(paste0(raw_dir, sample,
                               "_filtered_feature_bc_matrix.h5"))
  # Split the H5 into GEX and HTO
  gex <- gex_hto$`Gene Expression`
  hto <- gex_hto$`Multiplexing Capture`
  # Get colnames of all overlapping columns
  joint_bcs <- intersect(colnames(gex), colnames(hto))
  # Subset RNA and HTO counts by joint cell barcodes
  gex <- gex[, joint_bcs]
  hto <- hto[, joint_bcs]
  # Make Seurat object from GEX counts + add HTOs
  # Filter cells with low counts + features with low cell counts
  seu <- CreateSeuratObject(counts = gex)
  seu[["HTO"]] <- CreateAssayObject(counts = hto)
  # Normalize data (centered log transform for HTOs)
  seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR")
  # Demultiplex HTOs. Asign each cell as singlet, doublet, or negative
  # The default threshold is 0.99 positive quantile
  # The lower this setting, the "easier" an HTO is assigned to a cell, thus
  # you will get more doublets, less negatives and less singlets
  seu <- HTODemux(seu, assay = "HTO", positive.quantile = 0.99)
  # Table to record statistics of demultiplexing etc
  stats <- data.frame(
    n_bcs_gex = length(colnames(gex)),
    n_bcs_hto = length(colnames(hto)),
    n_bcs_overlap = length(joint_bcs),
    n_doublet = table(seu$HTO_classification.global)["Doublet"],
    n_negative = table(seu$HTO_classification.global)["Negative"],
    n_singlet = table(seu$HTO_classification.global)["Singlet"],
    n_HTO_1 = table(seu$hash.ID)["TotalSeq-C0251 anti-human Hashtag 1 Antibody"],
    n_HTO_2 = table(seu$hash.ID)["TotalSeq-C0251 anti-human Hashtag 2 Antibody"],
    n_HTO_3 = table(seu$hash.ID)["TotalSeq-C0251 anti-human Hashtag 3 Antibody"],
    n_HTO_4 = table(seu$hash.ID)["TotalSeq-C0251 anti-human Hashtag 4 Antibody"],
    row.names = sample)
  write.csv(stats,
            paste0(output_dir,"demultiplex_", sample,".csv"),
            row.names = TRUE)
  # Group cells based on the max HTO signal
  Idents(seu) <- "HTO_maxID"
  # Plot graphs showing distribution of HTO
  pdf(paste0(output_dir,"demultiplex_", sample,".pdf"), height = 18, width = 8)
  print(
    RidgePlot(seu, assay = "HTO", features = rownames(seu[["HTO"]]), ncol = 1))
  dev.off()
  # Quickly prep stats for ggplot and plot stats
  stats <- data.frame(t(stats)) %>% mutate(Statistics = names(stats))
  print(ggplot(stats, aes_string(x = "Statistics", y = sample)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          theme(axis.title = FALSE) + 
          theme_classic())
  ggsave(paste0(output_dir,"demultiplex_", sample,".png"))
  # Filter cells that have only 1 specific HTO (remove doublets + negatives)
  Idents(seu) <- "HTO_classification.global"
  seu <- subset(seu, idents = "Singlet")
  # Add sample information
  seu$sample <- sample
  # Return
  return(seu)
}


# Visualise cell filtering technique + returns suggested cell filter
# It calculates the log10 of counts per cell and uses median -/+ 3 mad
# MAD = Median absolute deviation
visualise_cell_filter <- function(seu) {
  # Calculate thresholds
  tot_cts <- log10(seu@meta.data$nFeature_RNA)
  tot_cts_median <- median(tot_cts)
  tot_cts_mad <- mad(tot_cts)
  print(paste("Median total features:", 10^tot_cts_median))
  print(paste("Mad total features:", 10^tot_cts_mad))
  # Set thresholds
  cell_thresholds <- list(low = tot_cts_median-3*tot_cts_mad,
                          high = tot_cts_median+3*tot_cts_mad)
  # Visualise using histogram
  hist(tot_cts,col="grey80",breaks = 100,
       xlab="log10(nFeature_RNA)", 
       main=paste("# of genes in ", seu@project.name),
       ylab="Number of cells"
  )
  abline(v = cell_thresholds$low, col = "red")
  abline(v = cell_thresholds$high, col = "red")
  # Return thresholds
  return(cell_thresholds)
}

# Create barplot with frequencies of cells per cluster by feature
plot_freq_plot <- function(clusters,
                           feature,
                           legend,
                           x_lab,
                           palette=NULL,
                           n_cells=FALSE,
                           table=FALSE) {
  # Create ggplot dataframe
  ggdata <- data.frame(table(clusters, feature))
  # Calculate percentage of cells per clusters for the provided featres
  for(i in 1:nrow(ggdata)) {
    sum <- sum(ggdata[ggdata[,"feature"]==ggdata[i,"feature"], "Freq"])
    ggdata[i,"Perc"] <- (ggdata[i,"Freq"]/sum)*100
  }
  # Plot graph with frequencies
  plot <- ggplot(ggdata, aes(x=feature, y=Perc, fill=clusters)) +
    geom_bar(stat="identity", position="stack") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(x_lab) + ylab("Frequency (%)") + labs(fill=legend)
  # Add custom color palette
  if(!is.null(palette)) plot <- plot + scale_fill_manual(values=palette)
  print(plot)
  # Plot second graph with ncells and not frequencies
  if(isTRUE(n_cells)) {
    print(ggplot(ggdata, aes(x=feature, y=Freq, fill=clusters)) +
            geom_bar(stat="identity", position="stack") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            scale_y_continuous(expand = c(0, 0)) +
            xlab(x_lab) + ylab("Number of cells") + labs(fill=legend))
  }
  # Return table
  if(isTRUE(table)) return(ggdata)
}


#' Differential pairwise testing using Seurats FindMarkers function
#' Requires settings file (00_settings.R)
#' @param seu Seurat object
#' @param column Which metadata column to run comparison on
#' @param condition1 First condition
#' @param condition2 Second condition
#' @param prefix Prefix used to save the csv file. NULL = CSV not saved.
#' @param ... Passed through to "FindMarkers" seurat function
differential_testing <- function(seu,
                                 column,
                                 condition1,
                                 condition2,
                                 prefix = "DEG_",
                                 ...) {
  # Save condition as idents
  Idents(seu) <- seu[[column]]
  # Do pairwise comparisons
  res <- FindMarkers(seu,
                     ident.1 = condition1,
                     ident.2 = condition2,
                     test.use = "MAST",
                     logfc.threshold = 0.25,
                     ...)
  res_sub <- res[res$p_val_adj<0.05,]
  # Write results
  if(!is.null(prefix))
    write.csv(res_sub, paste0(dirs[["deg"]], prefix, condition1, "_vs_",
                               condition2, ".csv"))
  return(res_sub)
}


#' Function for calculating AUC based on GSEA gene list and seurat object
#' Exports csv with AUCmatrix results
#' Requires Seurat, AUCell and GSEABase
#' nCores and output_dir are set based on settings lists
#' @param seu Seurat object with seurat_clusters and scaled data
#' @param gsea GSEA geneset data (getGmt from GSEABase)
#' @param title Title for plots and AUCMatrix CSV file
#' @param ncores Indicate how many cores should be used
#' @param output_dir Directory to store AUCmatrix as CSV output
run_AUC <- function(seuset,
                    gsea,
                    title,
                    ncores = settings[["ncores"]],
                    output_dir = dirs[["auc"]]) {
  if(length(gsea) < 2) stop("There must be more than 1 sets in a geneset.")
  # Make the expression matrix from the seurat object 
  exprMatrix <- as.matrix(seuset@assays$RNA@data)
  # Check how many genes of the geneset are in the expression matrix 
  gsea <- subsetGeneSets(gsea, rownames(exprMatrix)) 
  print(cbind(nGenes(gsea)))
  # Add the gene-set size into its name
  gsea <- setGeneSetNames(gsea,
                               newNames=paste(names(gsea),
                                              "_",
                                               nGenes(gsea), 
                                              "_genes",
                                              sep=""))
  # Rankings: Quantiles for the number of genes detected by cell: 
  cells_rankings <- AUCell_buildRankings(exprMatrix,
                                         nCores = ncores,
                                         plotStats = FALSE) 
  # Calculate enrichment for the gene signatures (AUC)
  cells_AUC <- AUCell_calcAUC(gsea,
                              cells_rankings,
                              nCores = ncores) 
  AUCmatrix <- getAUC(cells_AUC) %>% t(.)
  write.table(AUCmatrix, paste0(output_dir, title, ".csv"))
}


#' Function for plotting AUCmatrix output in UMAP + Heatmap
#' Requires Seurat, ggplot2, and pheatmap
#' @param seu Seurat object with seurat_clusters and scaled data
#' @param group_by groups to plot the heatmap by. Default: seurat_clusters
#' @param AUCmatrix AUCmatrix output from run_AUC function
#' @param scales Indicates different scaling to use for heatmap.
plot_AUC <- function(seuset,
                     group_by = "seurat_clusters",
                     AUCmatrix,
                     scales = c("none"),
                     cluster_cols = FALSE) {
  # Define colors
  umap_cols <- c("lightgrey", "lemonchiffon1", "tan1", "firebrick3")
  hm_cols <- colorRampPalette(c("steelblue", "slategray2", "white",
                                "tan1", "firebrick3"))(100)
  # Loop for each grouping to plot
  for(group in group_by) {
    # Plot the AUC in a UMAP using seurat clusters
    plot <- DimPlot(seuset,
                    group.by = group,
                    pt.size = 0.3,
                    label.size = 5, 
                    reduction = "umap")
    # Make ggplot compatible dataset with AUC values
    df <- as.data.frame(cbind(plot$data$UMAP_1,
                              plot$data$UMAP_2,
                              plot$data[group]),
                        row.names = rownames(plot$plot_env$data)) %>%
      setNames(c("UMAP1", "UMAP2", "GROUPING")) %>%
      merge(., AUCmatrix, by = 0, all = TRUE) %>%
      `rownames<-` (rownames(plot$plot_env$data)) %>%
      dplyr::select(., -1)
    # Plot UMAPs
    if(group == "seurat_clusters") {
      sapply(names(df[,-(1:3)]), FUN = function(set)
        print(ggplot(df, aes_string(x = "UMAP1", y= "UMAP2", colour = set)) +
                geom_point(size = 1) +
                scale_colour_gradientn("AUC", colors = umap_cols) +
                theme_classic() +
                xlab("UMAP1") +
                ylab("UMAP2") +
                ggtitle(paste("Comparison with", set)))
      )
    }
    # Plot the AUC in a heatmap with all groupings
    df2 <- aggregate(.~GROUPING, df[,-(1:2)], FUN=mean)
    rownames(df2) <- df2$GROUPING
    df2 <- as.matrix(df2[,-1])
    sapply(scales, FUN = function(scale)
      print(pheatmap(df2,
                     color=hm_cols,
                     scale = scale,
                     cluster_rows = F,
                     cluster_cols = cluster_cols,
                     main = paste("Scale:", scale),
                     border_color = NA))
    )
  }
}


#' Creates a list + plots of enriched GO terms and saves them to results
#' @param genes List of genes to run GO on. Minimum 5 genes. Symbol.
#' @param filename Output filename
#' @param folder Output folder. Default: current folder
#' @param pvalue pvalue threshold for enrichment analyis
#' @param qvalue qvalue threshold for enrichment analysis
#' @param orgDB Organism database, e.g. org.Hs.eg.db for Human
go_enrichment <- function(genes,
                          title,
                          folder="/",
                          pvalue=0.01,
                          qvalue=0.05,
                          orgDB = org.Hs.eg.db){
  # Only run script if the DEGlist is greater than 5 genes
  if (length(genes) < 5) {
    warning(paste("Not enough genes for", title,". List must be > 5."))
  } else {
    # Run GO analysis
    res <- enrichGO(gene = genes,
                    OrgDb = orgDB,
                    keyType = "SYMBOL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = pvalue,
                    qvalueCutoff = qvalue,
                    readable = FALSE)
    # Only plot graphs if GO terms are enriched
    if(is.null(res) || nrow(res) == 0)
      { 
        warning(paste("No GO terms enriched for", title))
    } else {
      # Write table
      write.csv(res, paste0(folder, "GO_", title, ".csv"))
      # Export GO graphs as pdf
      #pdf(paste0(folder, "GO_", title, ".pdf"), useDingbats = F)
      print(dotplot(res, showCategory = 20) + ggtitle(paste("Dotplot", title)))
      #print(barplot(res, font.size=6) + ggtitle(paste("Barplot", title)))
      #print(emapplot(res, font.size=6) + ggtitle(paste("Emapplot", title)))
      #print(cnetplot(res, font.size=6) + ggtitle(paste("Cnetplot", title)))
      #dev.off()
    }
  }
  return(res)
}


#' Enrichment analysis for genes in database using enrichr (human)
#' Plots a bargraph and exports csv
#' @param genes Array of genes.
#' @param title Title for plot and saving csv.
#' @param database Which database to use (https://maayanlab.cloud/Enrichr/)
enrich_enrichr <- function(genes,
                           title,
                           database,
                           bargraph = FALSE,
                           folder = NULL) {
  # Enrichment
  enriched <- enrichr(genes, database)[[1]]
  # write csv
  if(!is.null(folder)) write.csv(enriched, paste0(folder, title, ".csv"))
  # Plot bargraph
  if(isTRUE(bargraph)) print(plotEnrich(enriched) + ggtitle(title))
  # Return
  return(enriched)
}

#' Plot slingshot graphs using ggplot2
#' Slingshot uses baseR plotting functions to plot results. This is very
#' restrictive. This function converts slingshot output to ggcompatible format
#' and plots graph using ggplot.
#' @param sce_sling Provide sce_sling object containing getLineages + getCurves
#' @param colors Provide a color palette to color UMAPs by.
#' @param color_by Default would be clusters, e.g. seurat_clusters
#' @param dim Which dim reduction to use. Usually PCA or UMAP
#' @param pseudotime Indicate if each lineage should be plotted with pseudotime.
slingshot_to_ggplot <- function(sce_sling,
                                colors = NULL,
                                color_by,
                                dim = "PCA",
                                pseudotime = TRUE) {
  # Convert SCE PCA info to ggplot-compatible dataframe
  gg_dim <- as.data.frame(sce_sling@int_colData$reducedDims[[dim]][,1:2]) %>%
    cbind(sce_sling[[color_by]]) %>%
    setNames(c("dim1", "dim2", color_by))
  # Convert MST data from slingshot to ggplot-compatible dataframe
  # Rename columns from Dim.1 and Dim.2 to dim1 and dim2 to match gg_dim df
  gg_mst <- slingMST(sce_sling, as.df = TRUE) %>%
    rename("dim1" = names(.[1])) %>% rename("dim2" = names(.[2]))
  # Retrieve curves information from slingshot analysis
  # Rename columns from Dim.1 and Dim.2 to dim1 and dim2 to match gg_dim df
  gg_curves <- slingCurves(sce_sling, as.df = TRUE) %>%
    rename("dim1" = names(.[1])) %>% rename("dim2" = names(.[2]))
  # Plot graph + MST
  print(ggplot(gg_dim, aes(x = dim1, y = dim2)) +
          geom_point(aes_string(col = color_by), size = 0.4) +
          { if(!is.null(colors)) scale_color_manual(values = colors) } +
          theme_classic() +
          geom_point(data = gg_mst, size = 5) +
          geom_path(data = gg_mst %>% arrange(Order), aes(group = Lineage), size = 1))
  # Plot graph + Curves
  print(ggplot(gg_dim, aes(x = dim1, y = dim2)) +
          geom_point(aes_string(col = color_by), size = 0.4) +
          { if(!is.null(colors)) scale_color_manual(values = colors) } +
          theme_classic() +
          geom_path(data = gg_curves %>% arrange(Order),
                    aes(group = Lineage),
                    size = 1))
  # Plot a graph for each lineage with pseudotime
  if(pseudotime == TRUE) {
    for(i in unique(gg_curves$Lineage)) {
      # Add pseudotime information for lineage
      gg_dim$pseudotime <- sce_sling[[paste0("slingPseudotime_",i)]]
      # Plot graph
      print(ggplot(data = gg_dim, aes(x = dim1, y = dim2)) +
              geom_point(aes(color = pseudotime)) +
              scale_color_viridis() +
              theme_classic() +
              geom_path(data = gg_curves[gg_curves$Lineage==i,],
                        aes(group = Lineage),
                        size = 1) +
              ggtitle(paste0("Lineage ",
                             i, ": ",
                             toString(unlist(slingLineages(sce_sling)[i])))))
    }
  }
}


##### Save session info -------------------------------------------------------
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

