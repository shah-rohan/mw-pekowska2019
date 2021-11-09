# This script is executed using r_pekowska2019_rev_analysis conda env
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
#library(HGNChelper)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# By the time of the initial analysis (2018-10-05), HGNCHelper was not available (Published: 2019-10-24)
# Differences in output with and without symbol update are minor anyways.
outprefix <- "GOBP_Chen_Broad_no_symbol_update"

# Input data from:
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.3385/MediaObjects/41588_2015_BFng3385_MOESM25_ESM.xls
d <- read_xls("inp/xls/chen2015/ng.3385-S4.xls")
names(d) <- make.names(names(d))
# Annotations for representative samples from each tissue/cell line were manually gathered
anno <- read_xls("inp/xls/sample_annotations.xls")

m <- as.matrix(d[,anno$Chen2015_sample_name])
colnames(m) <- anno$Cell_origin

# Some Symbols are outdated since Chen2015 publication.
# Updating to current Symbol for more accurate clusterProfiler results
#d$genes <- checkGeneSymbols(d$genes)$Suggested.Symbol

# We only want to query genes associated with broad peaks,ie with peak around TSS > 4000bp
m <- m > 4000 

broad_genes <- apply(
  X = m,
  MARGIN =2,
  FUN = function(x) {
    d[x,]$genes
  }
)
res <- compareCluster(
  geneClusters = broad_genes,
  fun = "enrichGO",
  keyType = "SYMBOL",
  OrgDb = "org.Hs.eg.db",
  ont = "BP"
)
saveRDS(
  res,
  file = paste0(
    outprefix,
    "_res.rds"
  )
)
simplified_res <- simplify(
  res,
  cutoff = 0.1
)
saveRDS(
  simplified_res,
  file = paste0(
    outprefix,
    "_simplified_res.rds"
  )
)
#dotplot(simplified_res)

dc <- reshape2::dcast(
  simplified_res@compareClusterResult,
  Description ~ Cluster,
  value.var = "pvalue"
)

col_fun = colorRamp2(
  c(1      , 2    , 30),
  c("white", "red", "red4")
)
dm <- -log10(as.matrix(dc[,-1]))
rownames(dm) <- dc[,1]
dm[is.na(dm)] <- 0

# We want to keep ontology elements that have the lowest p-value in at least one of the samples
max_dm <- apply(
  dm, 
  2, 
  FUN = function(x) {
    x == max(x)
  }
)
dm <- dm[rowSums(max_dm) > 0,]
dm <- dm[
  do.call(
    order,
    as.data.frame(dm)
  ),
]
anno_cols_to_display <- c("Project_source", "Antibody")
pdf(
  file = paste0(
    outprefix,
    "_heatmap.pdf"
  ),
  width = 13.5,
  height = 6.5
)
set.seed(2) # seed is here to get reasonnably good annotation colors
Heatmap(
  dm,
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  heatmap_legend_param = list(
    col_fun = col_fun,
    title = "-log10(p-value)",
    at     = c(1   , 2  , 30),
    labels = c("<2", "2", ">30")
  ),
  bottom_annotation = HeatmapAnnotation(
    df = data.frame(anno[, anno_cols_to_display])
  )
)
dev.off()

# Alternative filtering where only the best GO term by sample is kept
dm <- apply(
  dm,
  2,
  FUN = function(x) {
    ifelse(
      x == max(x),
      x,
      0
    )
  }
)
dm <- dm[
  do.call(
    order,
    as.data.frame(dm)
  ),
]
pdf(
  file = paste0(
    outprefix,
    "_top1_heatmap.pdf"
  ),
  width = 13.5,
  height = 6.5
)
set.seed(2) # seed is here to get reasonnably good annotation colors
Heatmap(
  dm,
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  heatmap_legend_param = list(
    col_fun = col_fun,
    title = "-log10(p-value)",
    at     = c(1   , 2  , 30),
    labels = c("<2", "2", ">30")
  ),
  bottom_annotation = HeatmapAnnotation(
    df = data.frame(anno[, anno_cols_to_display])
  )
)
dev.off()

