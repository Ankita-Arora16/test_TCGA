library(TCGAbiolinks)
library(SummarizedExperiment)


query= GDCquery(project='TCGA-BRCA',
                data.category='Transcriptome Profiling',
                data.type='Gene Expression Quantification',
                workflow.type = 'HTSeq - Counts')

# to select results based on cases 
samplesTable = getResults(query,cols=c('cases') )

#to view sample table 
#View(sampleTable) #1222 samples

dataTP= TCGAquery_SampleTypes(barcode= samplesTable, typesample = 'TP')
View(dataTP)
#1102 entries/ samples

dataNT= TCGAquery_SampleTypes(barcode = samplesTable, typesample='NT')
View(dataNP) #133 enteries

#to take out data for first 100 samples from both 

sample_dataTP= dataTP[1:100]
sample_dataNT = dataNT[1:100]

#preparing data table for the 
samplequery = GDCquery(project='TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       data.type='Gene Expression Quantification',
                       workflow.type='HTSeq - Counts',
                       barcode= c(sample_dataNT,sample_dataTP))

#to download the data
GDCdownload(query= samplequery)

#to prepare data 

prepdata = GDCprepare(query= samplequery, save= TRUE, save.filename = 'Geneexpdata.rda')


#Preparing data for preprocessing- removing outliners

data_preprocessed = TCGAanalyze_Preprocessing(object = prepdata, 
                                              cor.cut= 0.6,
                                              datatype = 'HTSeq - Counts')


#normalization of data

dataNorm = TCGAanalyze_Normalization(tabDF = data_preprocessed, 
                                     geneInfo = geneInfoHT,
                                     method = 'gcContent')

#to filter data out 

dataFilt = TCGAanalyze_Filtering(tabDF = dataNorm,
                                 method = 'quantile',
                                 qnt.cut = 0.25)

#applying filter to normal samples

samplesNT = TCGAquery_SampleTypes(colnames(dataFilt),
                                  typesample = c('NT'))
#applying filter to tumor samples

samplesTP = TCGAquery_SampleTypes(colnames(dataFilt),
                                  typesample = c('TP'))



#DEG Analysis using edgeR

dataDeg = TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                          mat2 = dataFilt[,samplesTP],
                          Cond1type = 'Normal',
                          Cond2type = 'Tumor',
                          fdr.cut = 0.01,
                          logFC.cut = 1,
                          method = 'glmLRT')

#DEG table with expression values in normal and tumor sample

dataDEG_table = TCGAanalyze_LevelTab(dataDeg, "Tumor",'Normal',
                                             dataFilt[,samplesTP],dataFilt[,samplesNT])

#writing a csv file for the table

write.csv(dataDEG_table, 'DEGtable.csv', quote = FALSE)

#filtering the genes with logFc greater than 6

DEG_filter = dataDeg[which(abs(dataDeg$logFC) >= 2), ]

#retrieving the table

write.csv(DEG_filter, 'Diff_expr_genes.csv', quote= FALSE)