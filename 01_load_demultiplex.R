###############################################################################
#' Malte Borggrewe
#' Load and demultiplex with using HTO
###############################################################################


source("src/00_library.R")


# Load H5 matrix, convert to Seurat object, and demultiplex using HTO
seu_donor1 <- load_demultiplex("donor1", dirs[["raw"]], dirs[["demux"]])
seu_donor2 <- load_demultiplex("donor2", dirs[["raw"]], dirs[["demux"]])
seu_donor2 <- load_demultiplex("donor3", dirs[["raw"]], dirs[["demux"]])

# Merging will only use raw counts; normalization is removed
seu <- merge(seu_donor1,
             y = c(seu_donor2, seu_donor3),
             add.cell.ids = c("donor1", "donor2", "donor3"),
             project = "TLR7_pDC")

# Add some metadata
seu$condition <- seu$hash.ID %>%
  gsub("TotalSeq-C0251 anti-human Hashtag 1 Antibody", "T0 Unstim", .) %>%
  gsub("TotalSeq-C0251 anti-human Hashtag 2 Antibody", "T0 CL097", .) %>%
  gsub("TotalSeq-C0251 anti-human Hashtag 3 Antibody", "T3 Unstim", .) %>%
  gsub("TotalSeq-C0251 anti-human Hashtag 4 Antibody", "T3 CL097", .)
seu$time <- substr(seu$condition, 1, 2)
seu$treatment <- substr(seu$condition, 4, 10)


# Save RDS
saveRDS(seu, paste0(dirs["rds"],"seu_unfiltered.RDS"))

