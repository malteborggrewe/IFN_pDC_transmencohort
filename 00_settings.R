###############################################################################
#' Malte Borggrewe
#' Settings: defining global variables
###############################################################################

#' Initialize and define project directories
#' @param output_dir Change default output dir as desired.
#' Output_dir can be used upon function calling to indicate a subdir in output/.
project_dirs <- function(output_dir = "output/") {
  if(output_dir != "output/")
    output_dir = paste0("output/", output_dir, "/")
  dirs <- list(
    data = "data/",
    raw = "data/raw/",
    meta = "data/metadata/",
    rds = "data/rds/",
    output = output_dir,
    deg = paste0(output_dir, "DEG/"),
    auc = paste0(output_dir, "AUC/"),
    qc = "output/qc/",
    demux = "output/qc/demux/",
    reports = "reports/",
    src = "src/")
  for(d in dirs){
    dir.create(d, showWarnings = FALSE)
  }
  return(dirs)
}
dirs <- project_dirs()

# Seurat settings
seurat <- list(
  min_cells = 3, #Default: 3
  min_features = 870, #Default: 200
  max_features = 4670, #Default: 2500
  max_mt = 5, #Default: 5
  variable_features = 2000,
  ndims_use = 10,
  resolution_corrected = 0.8
)

# Colorpalette
colors <- list(
  seurat_clusters = c("1"="#DBD61F",
                      "2"="#08AEB2",
                      "3"="#FC7E87",
                      "4"="#C42380",
                      "5"="#D27FF4",
                      "6"="#5BF4D6",
                      "7"="#E0973F",
                      "8"="#55CE60",
                      "9"="#A129D3",
                      "10"="#E6194B",
                      "11"="#407EE2",
                      "12"="#ED2BA7",
                      "13"="#AAE012",
                      "14"="#42D4F4",
                      "15"="#72A379",
                      "16"="#F2612F",
                      "17"="#15377C"),
  condition = c("TP0 Unstim"="#C42380",
                "TP0 CL097"="#08AEB2",
                "TP3 Unstim"="#AAE012",
                "TP3 CL097"="#E0973F"),
  treatment = c("Unstim"="#C42380",
                "CL097"="#407EE2")
)

# Other settings
settings <- list(
  ncores = 20
)

