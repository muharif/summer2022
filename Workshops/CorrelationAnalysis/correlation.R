install.packages("pheatmap")

# load package
library(pheatmap)

metadata_file = 'data/metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

tpm_file_nc = 'data/tpm_nc.txt' #adjust this based on your file location
tpm_nc = read.csv(tpm_file_nc,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

tpm_file_pc = 'data/tpm_pc.txt' #adjust this based on your file location
tpm_pc = read.csv(tpm_file_pc,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

ECS_genes = c('CNR1', 'CNR2', 'NAPEPLD',
              'GDE1', 'GDPD3',
              'PLA2G4E', 'PLAAT1', 'PLAAT2',
              'PLAAT3', 'PLAAT4', 'PLAAT5',
              'DAGLA', 'DAGLB', 'NAAA',
              'FAAH', 'MGLL', 'ABHD6', 'ABHD12')

## patient Correlation

#PC
high_exp = tpm_pc[rowMeans(tpm_pc) > 5,]

corr = cor(high_exp, method = 'spearman')

cutree_rows = 4
cutree_cols = 4

pheatmap(corr, show_rownames	= FALSE, 
         cutree_rows = cutree_rows, cutree_cols = cutree_cols,
         annotation_row = metadata,
         annotation_col = metadata,
)

#NC
high_exp = tpm_nc[rowMeans(tpm_nc) > 5,]

corr = cor(high_exp, method = 'spearman')

cutree_rows = 5
cutree_cols = 5

pheatmap(corr, show_rownames	= FALSE, 
         cutree_rows = cutree_rows, cutree_cols = cutree_cols,
         annotation_row = metadata,
         annotation_col = metadata,
)

## Heatmap Overview of ECS vs Non-Coding

corr = cor(t(tpm_nc), t(tpm_pc[ECS_genes,]))

#To remove the missing correlations
corr = corr[complete.cases(corr), ]

# Select top X most correlated genes in each
top_X = 200

selected_genes = c()
for(gene_ecs in ECS_genes){
  temp = names(corr[order(abs(corr[,gene_ecs]), decreasing = TRUE)[0:top_X], gene_ecs])
  selected_genes = append(selected_genes, temp)
}


selected_genes = unique(selected_genes)

cutree_rows = 8
cutree_cols = 5

p = pheatmap(corr[selected_genes,], show_rownames	= FALSE, cutree_rows = cutree_rows, cutree_cols = cutree_cols)

row_clust = as.factor(cutree(p$tree_row, cutree_rows))
row_clust = as.data.frame(row_clust)
col_clust = as.factor(cutree(p$tree_col, cutree_cols))
col_clust = as.data.frame(col_clust)

pheatmap(corr[selected_genes,], show_rownames	= FALSE, 
         cutree_rows = cutree_rows, cutree_cols = cutree_cols,
         annotation_row = row_clust,
         annotation_col = col_clust,
         breaks = seq(-1, 1, length.out = 100)
         )



## Single Gene Correlation

gene_pc = 'CNR1'
gene_nc = 'ENSG00000117242'

corr = cor.test(as.numeric(tpm_pc[gene_pc,]), as.numeric(tpm_nc[gene_nc,]), method = "spearman")
corr

## One-gene to all

gene_pc = 'CNR1'

results = c()

for(gene_nc in rownames(tpm_nc)){
  corr = cor.test(as.numeric(tpm_pc[gene_pc,]), as.numeric(tpm_nc[gene_nc,]), method = "spearman")
  results = append(results, list(c(gene_nc, corr$estimate, corr$p.value)))
  
}

results = t(as.data.frame(results))
colnames(results) = c('Gene', 'Correlation', 'P-Value')

FDR = p.adjust(results[,'P-Value'], method = 'fdr')

results = cbind(results, FDR)

write.table(results,file='Correlation_CNR1.txt',sep = '\t', na = '',row.names = F)

## All ECS genes to non-coding

results = c()
for(gene_pc in ECS_genes){
  for(gene_nc in rownames(tpm_nc)){
    corr = cor.test(as.numeric(tpm_pc[gene_pc,]), as.numeric(tpm_nc[gene_nc,]), method = "spearman")
    results = append(results, list(c(gene_pc, gene_nc, corr$estimate, corr$p.value)))
  }
}
results = t(as.data.frame(results))
colnames(results) = c('GenePC', 'GeneNC', 'Correlation', 'P-Value')

FDR = p.adjust(results[,'P-Value'], method = 'fdr')

results = cbind(results, FDR)

write.table(results,file='Correlation_ECS.txt',sep = '\t', na = '',row.names = F)

