## Execute this just once
install.packages("factoextra")


## Draw PCA
library(factoextra)
library(FactoMineR) 

metadata_file = 'yourmetadata' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

tpm_file = 'yourtpm' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

#To run the PCA
res.pca <- PCA(t(tpm), scale = TRUE, graph =  FALSE)

# Visualizing the percentage of explained variance
fviz_eig(res.pca)

# Draw the PCA Plot
fviz_pca_ind(res.pca,
            habillage = as.factor(metadata$disease_state), 
            palette = c("green", "grey"), # Color, support RGB HEX too
             label = 'none',
            pointsize = 5, # Size
            invisible="quali",
            repel=TRUE,
            
)+
  labs(title ="", x = paste0("PC1 (", signif(res.pca$eig[1,2], 4) ,"%)"), y = paste0("PC2 (", signif(res.pca$eig[2,2], 4) ,"%)")) + 
  theme_minimal()



