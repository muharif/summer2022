### ForBen
metadata_file = 'yourmetadata' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

x = c()
for(i in rownames(metadata)){
  x = c(x,paste0('X',i))
}

rownames(metadata) = x

metadata = metadata[metadata$Strain != 'Unknown',]

tpm_file = 'yourtpm' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

library(ggplot2)
combined = cbind(t(tpm), metadata)

### BEN
#For Single Category
order =   c('Vehicle', 'TGFB1', 'Rimonabant', '1400W', 'MRI1867')

ggplot(combined, aes(x = as.factor(Treatment), y = NOS2,fill = Strain)) +  # Y is the selected gene
  scale_x_discrete(limits = c('Vehicle', 'TGFB1', 'Rimonabant', '1400W', 'MRI1867')) +
  geom_bar(stat="identity", position = "dodge") +
  labs(title ="", x = 'Conditions', y = 'TPM Values') +
  theme_minimal()


#For Multi Category
ggplot(combined, aes(x = as.factor(Strain), y = ABHD14B,fill = Treatment)) + # Y is the selected gene
  geom_bar(stat="identity", position = "dodge") +
  labs(title ="", x = 'Conditions', y = 'Expression Values') +
  theme_minimal()


### ANGELINA
metadata_file = 'yourmetadata' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

tpm_file = 'yourtpm' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

library(ggplot2)
combined = cbind(t(tpm), metadata)

#For Single Category

ggplot(combined, aes(x = as.factor(disease_state), y = ENSG00000256150)) +  # Y is the selected gene
  geom_bar(stat="identity", position = "dodge") +
  labs(title ="", x = 'Conditions', y = 'TPM Values') +
  theme_minimal()
