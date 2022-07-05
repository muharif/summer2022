###Do this just once to install
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano', force = TRUE)
####

library(EnhancedVolcano)

#Loading the deseq file
deseq_file='deseq.txt' # Adjust this to your own file name if necessary
deseq = read.csv(deseq_file,sep='\t',stringsAsFactors = F,row.names = 1)

EnhancedVolcano(deseq,
                lab = rownames(deseq),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = rownames(deseq[order(deseq$padj), ])[0:10], ## To label top 10 lowest padj genes, you can change it if needed
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'P-Adjusted'),
                pCutoff = 0.05,
                pCutoffCol = 'padj',
                FCcutoff = 0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                # labFace = 'bold',
                # boxedLabels = TRUE,
                # colAlpha = 4/5,
                # legendPosition = 'right',
                # legendLabSize = 14,
                # legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black'
                )
## Save using the export button on top of the plot,
## my go to is usually size 8x9 potrait, but feel free to adjust