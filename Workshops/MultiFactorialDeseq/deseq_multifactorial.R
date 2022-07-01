metadata_file = 'metadata.txt' #adjust this based on your file location
metadata = read.csv(metadata_file,sep='\t',stringsAsFactors = F,row.names = 1)

count_file = 'count_pc.txt' #adjust this based on your file location
count = read.csv(count_file,sep='\t',stringsAsFactors = F,row.names = 1)[,rownames(metadata)]

library('DESeq2')
coldata = metadata[,c('Batch', 'DiseaseStatus', 'tobacco', 'gender')]
coldata$Batch = as.factor(coldata$Batch)
coldata$DiseaseStatus = as.factor(coldata$DiseaseStatus)
coldata$tobacco = as.factor(coldata$tobacco)
coldata$gender = as.factor(coldata$gender)

dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(count)),colData=coldata,design=~Batch+gender+DiseaseStatus)
dds <- DESeq(dds)

col = 'DiseaseStatus'
cond1 = 'IPF'
cond2 = 'Normal'
res=results(dds,contrast=c(col,cond1,cond2))
res=data.frame(res)
