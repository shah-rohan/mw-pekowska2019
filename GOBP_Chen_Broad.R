library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(HGNChelper)
library(reshape2)
library(ComplexHeatmap)
library(circlize)

# Input data from:
# https://static-content.springer.com/esm/art%3A10.1038%2Fng.3385/MediaObjects/41588_2015_BFng3385_MOESM25_ESM.xls
d <- read_xls("inp/xls/chen2015/ng.3385-S4.xls")
names(d) <- make.names(names(d))
# Annotations manually gathered by Salvatore Spicuglia
anno <- read_xls("inp/xls/sample_annotations.xls")

m <- as.matrix(d[,anno$Chen2015_sample_name])
colnames(m) <- anno$Cell_origin

# Some Symbols are outdated since Chen2015 publication.
# Updating to current Symbol for more accurate clusterProfiler results
d$genes <- checkGeneSymbols(d$genes)$Suggested.Symbol

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
simplified_res <- simplify(
  res,
  cutoff = 0.1
)

dotplot(simplified_res)

dc <- reshape2::dcast(
  res@compareClusterResult,
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
anno_cols_to_display <- c("Project_source", "Antibody")
pdf(
  file = "GOBP_Chen_Broad_heatmap.pdf",
  width = 13.5,
  height = 6.5

)
set.seed(2) # seed is here to get reasonnably good annotation colors
Heatmap(
  dm[rowSums(max_dm) > 0,],
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

