#######################################################################
# The MIT License
#
# Copyright (c) 2017, Jérôme Audoux (jerome.audoux@inserm.fr); 2022, Ryan McMinds (mcmindsr@usf.edu)
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files 
# (the “Software”), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge,
# publish, distribute, sublicense, and/or sell copies of the Software,
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be 
# included in all copies or substantial portions of the Software.
#
# The Software is provided “as is”, without warranty of any kind, 
# express or implied, including but not limited to the warranties of 
# merchantability, fitness for a particular purpose and 
# noninfringement. In no event shall the authors or copyright holders
# be liable for any claim, damages or other liability, whether in an 
# action of contract, tort or otherwise, arising from, out of or in 
# connection with the software or the use or other dealings in the 
# Software. 
#######################################################################

library("data.table")
library("foreach")
library("doParallel")

args <- commandArgs(TRUE)

# Get parameters for the test
binary                    = args[1]#snakemake@input$binary
kmer_counts               = args[2]#snakemake@input$counts
sample_conditions         = args[3]#snakemake@input$sample_conditions
pvalue_threshold          = args[4]#snakemake@params$pvalue_threshold
log2fc_threshold          = args[5]#snakemake@params$log2fc_threshold
conditionA                = args[6]#snakemake@params$conditionA
conditionB                = args[7]#snakemake@params$conditionB
nb_core                   = args[8]#snakemake@threads
chunk_size                = args[9]#snakemake@params$chunk_size

# Get output files  
output_tmp                = args[10]#snakemake@output$tmp_dir
output_diff_counts        = args[11]#snakemake@output$diff_counts
output_pvalue_all         = args[12]#snakemake@output$pvalue_all
output_log                = args[13]#snakemake@log[[1]]


# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}

dir.create(output_tmp, showWarnings = FALSE)

logging(paste("Start binomial_diff_methods"))
logging(paste("pvalue_threshold", pvalue_threshold))
logging(paste("log2fc_threshold", log2fc_threshold))

## LOADING PRIOR KNOWN NORMALISATION FACTORS
logging(paste("loading colData:", date()))
colData = read.table(sample_conditions,header=T,row.names=1)

## LOAD KMER COUNTS
logging(paste("loading kmer counts:", date()))
kmer_count_data = read.table(kmer_counts,header=T,row.names=1)
kmer_count_data[kmer_count_data > 0] <- 1
logging(paste("finished loading kmer counts:", date()))

# Set the number of cores to use
registerDoParallel(cores=nb_core)

pvals <- numeric(nrow(kmer_count_data))
outdf <- data.frame(ID=rownames(kmer_count_data), 
                    pvalue=numeric(nrow(kmer_count_data)),
                    meanA=rowMeans(kmer_count_data[,rownames(colData)[colData$condition==conditionA]]), 
                    meanB=rowMeans(kmer_count_data[,rownames(colData)[colData$condition==conditionB]]),
                    log2FC=numeric(nrow(kmer_count_data)), ## actually logit; keeping name for compatibility
                    kmer_count_data)
## Binomial regression ANALYSIS ON EACH k-mer
logging(paste("starting parallel linear model fitting:", date()))
pv_log <- foreach(i=iter(kmer_count_data, by='row'), .combine=rbind) %dopar% {

  res <- glm(i ~ colData$condition + colData$normalization_factor, family=binomial(link='logit'))

  return(c(anova(res, test='LRT')$'Pr(>Chi)'[[2]], res$coefficients[[2]]))
  
}
logging(paste("finished parallel linear model fitting:", date()))

outdf$pvalue <- pv_log[,1]
outdf$log2FC <- pv_log[,2]

### print a table that can be used downstream
dir.create(dirname(output_pvalue_all))

outdf_filt <- outdf[,c('ID','pvalue')]
write.table(outdf_filt,
            file=gzfile(output_pvalue_all),
            sep="\t",
            quote=FALSE,
            col.names = FALSE,
            row.names = FALSE)

outdf_filt <- outdf[outdf$pvalue < pvalue_threshold, c('ID','pvalue','meanA','meanB','log2FC')]
colnames(outdf_filt)[1] <- 'tag'
write.table(outdf_filt,
            file=gzfile(output_diff_counts),
            sep="\t",
            quote=FALSE,
            col.names = TRUE,
            row.names = FALSE)

logging(paste("End binomial_diff_methods"))
