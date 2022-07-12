library(pheatmap)

combi_results <- function(folder, padj_thr, significant_minimum = 2){
  results = list()
  x = 0
  for(i in list.files(folder)){
    x = x+1
    data = read.table(paste0(folder,'/',i), sep = '\t', header  = 1, row.names = 1)
    selected = as.data.frame(data[, c('Stat..dist.dir.up.')])
    selected[!(data$p.adj..dist.dir.up. < padj_thr | data$p.adj..dist.dir.up. < padj_thr),] = 0
    rownames(selected) = rownames(data)
    colnames(selected) = c(gsub('.txt','',gsub('piano_','',i)))
    results[[x]] = selected
  }
  
  results = do.call(cbind, results)
  results = results[rowSums(results == 0) < dim(results)[2]-significant_minimum,]
}

folder = 'KEGG' ## Change this to the folder name where you store the functional analysis results.
padj_thr = 0.05 ## Adjust this threshold if you see too many significant terms in your results (aka unreadable)
significant_minimum = 2 ## This will show only terms that are significant in at least 2 comparisons. You can adjust this, 1 is best so we have the overall view, but sometimes it may results in too many terms
data = combi_results(folder, padj_thr, significant_minimum)

pheatmap(
  data,
  cluster_cols = FALSE,
  legend_labels = 'Gene-Level Statistics',
  breaks = seq(-max(abs(data)), max(abs(data)), length.out = 100)
  )
