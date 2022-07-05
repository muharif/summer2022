## PCA

metadata_file = 'ncRNA1/metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

tpm_file = 'ncRNA1/tpm_lncRNA.txt' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

PCA=prcomp(t(tpm), scale=F)
plot(PCA$x,pch = 15,col=c('blue','blue','red','red','lightgreen','lightgreen','black','black'))

# SHAM_1D --> Blue
# MI_1D --> red
# SHAM_3D --> lightgreen
# MI_3D --> black

summary(PCA)

metadata_file = 'ncRNA1/metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

count_file = 'ncRNA1/count_lncRNA.txt' #adjust this based on your file location
count = read.csv(count_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

library('DESeq2')
conds=as.factor(metadata$disease_state)
coldata <- data.frame(row.names=rownames(metadata),conds)
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~conds)
dds <- DESeq(dds)

cond1 = 'IPF' #First Condition
cond2 = 'Control' #Reference Condition
res=results(dds,contrast=c('conds',cond1,cond2))
res=data.frame(res)
  
write.table(res,file='deseq_1D.txt',sep = '\t', na = '',row.names = T,col.names=NA)

cond1 = 'MI_3D' #First Condition
cond2 = 'SHAM_3D' #Reference Condition
res=results(dds,contrast=c('conds',cond1,cond2))
res=data.frame(res)

write.table(res,file='deseq_3D.txt',sep = '\t', na = '',row.names = T,col.names=NA)
``