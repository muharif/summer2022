## PCA

metadata_file = 'HPSPF/metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

x = c()
for(i in rownames(metadata)){
  x = c(x,paste0('X',i))
}

rownames(metadata) = x

tpm_file = 'HPSPF/tpm.txt' #adjust this based on your file location
tpm = read.csv(tpm_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

PCA=prcomp(t(tpm), scale=F)
plot(PCA$x,pch = 15,col=c('blue','blue','red','red','lightgreen','lightgreen','black','black'))

# SHAM_1D --> Blue
# MI_1D --> red
# SHAM_3D --> lightgreen
# MI_3D --> black

summary(PCA)

metadata_file = 'HPSPF/metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

x = c()
for(i in rownames(metadata)){
  x = c(x,paste0('X',i))
}

rownames(metadata) = x

count_file = 'HPSPF/count.txt' #adjust this based on your file location
count = read.csv(count_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

library('DESeq2')
conds=as.factor(metadata$conds)
coldata <- data.frame(row.names=rownames(metadata),conds)
dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~conds)
dds <- DESeq(dds)

cond1 = 'HPS_Vehicle' #First Condition
cond2 = 'Control_Vehicle' #Reference Condition
res=results(dds,contrast=c('conds',cond1,cond2))
res=data.frame(res)

write.table(res,file='deseq_1D.txt',sep = '\t', na = '',row.names = T,col.names=NA)
