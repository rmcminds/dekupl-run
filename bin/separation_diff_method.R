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
library("MASS")
library("pscl")
library("doMPI")
## also requires RhpcBLASctl

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
chunk_size                = as.numeric(args[9])#snakemake@params$chunk_size
seed                      = args[14]#snakemake@params$seed

# Get output files  
output_tmp                = args[10]#snakemake@output$tmp_dir
output_diff_counts        = args[11]#snakemake@output$diff_counts
output_pvalue_all         = args[12]#snakemake@output$pvalue_all
output_log                = args[13]#snakemake@log[[1]]

# Temporary files
output_tmp_chunks         = paste(output_tmp,"/tmp_chunks/",sep="")
output_tmp_separ          = paste(output_tmp,"/tmp_separ/",sep="")
header_kmer_counts        = paste(output_tmp,"/header_kmer_counts.txt",sep="")
tmp_concat                = paste(output_tmp,"/tmp_concat.txt",sep="")
adj_pvalue                = paste(output_tmp,"/adj_pvalue.txt.gz",sep="")
dataseparAll              = paste(output_tmp,"/dataseparAll.txt.gz",sep="")
dataseparFiltered         = paste(output_tmp,"/dataseparFiltered.txt.gz",sep="")

# Function for logging to the output
logging <- function(str) {
  sink(file=paste(output_log), append=TRUE, split=TRUE)
  print(paste(Sys.time(),str))
  sink()
}
# Return the number of line in the last files of the directory
nbLineLastFile <- function(dir) {
  return(as.numeric(system(paste("zcat", paste(dir, "/$(ls ", dir, " |sort -n|grep subfile|tail -1)",sep=""), "| wc -l"), intern=TRUE)))
}
# Return the number of files in the directory
nbFiles <- function(dir) {
  return(as.numeric(system(paste("ls ", dir, "|grep subfile | wc -l"), intern=TRUE)))
}

# Set the number of cores to use
cl <- startMPIcluster() 
registerDoMPI(cl)

# Create directories
dir.create(output_tmp, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_chunks, showWarnings = FALSE, recursive = TRUE)
dir.create(output_tmp_separ, showWarnings = FALSE, recursive = TRUE)

logging(paste("Start separ_diff_methods"))
logging(paste("pvalue_threshold", pvalue_threshold))
logging(paste("log2fc_threshold", log2fc_threshold))

## LOADING PRIOR KNOWN NORMALISATION FACTORS
colData = read.table(sample_conditions,header=T,row.names=1)

# CLEAN THE TMP FOLDER
system(paste("rm -f ", output_tmp_chunks, "/*", sep=""))

# SAVE THE HEADER INTO A FILE
system(paste("zcat", kmer_counts, "| head -1 | cut -f2- >", header_kmer_counts))

# SHUFFLE AND SPLIT THE MAIN FILE INTO CHUNKS WITH AUTOINCREMENTED NAMES, ACCORDING TO SEED
logging(paste("Shuffling and splitting:", date()))
if(seed == 'fixed'){
    system(paste("zcat", kmer_counts, " >tmp_shuff; cat tmp_shuff| tail -n +2 | shuf --random-source=tmp_shuff | awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
             "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip --fast >\" output_tmp_chunks x}'"))
    system("rm tmp_shuff")
}else{
    system(paste("zcat", kmer_counts, "| tail -n +2 | shuf | awk -v", paste("chunk_size=", chunk_size,sep=""), "-v", paste("output_tmp_chunks=",output_tmp_chunks,sep=""),
                 "'NR%chunk_size==1{OFS=\"\\t\";x=++i\"_subfile.txt.gz\"}{OFS=\"\";print | \"gzip --fast >\" output_tmp_chunks x}'"))
}

logging("Shuffle and split done")

nb_line_last_file = nbLineLastFile(output_tmp_chunks)
nb_files = nbFiles(output_tmp_chunks)

# IF THE LAST FILE HAS LESS THAN HALF OF THE CHUNK SIZE
# CONCATENATE THE LAST 2 FILES AND THEN SPLIT
if(nb_files > 1 && nb_line_last_file < (chunk_size/2)) {

  ## CONCATENATE THE 2 FILES
  logging(paste("The last file has",nb_line_last_file,"line(s) it will be concatenated to the second last one"))

  before_last_file = paste(output_tmp_chunks, (nb_files - 1), "_subfile.txt.gz", sep="")
  last_file        = paste(output_tmp_chunks, (nb_files), "_subfile.txt.gz", sep="")

  #CONCATENATE THE LAST 2 FILES INTO A TMP FILE
  system(paste("cat", before_last_file, last_file, ">", tmp_concat, sep=" "))

  #NUMBER OF LINE OF THE TMP FILE
  nb_line_last_file = chunk_size + nb_line_last_file

  logging(paste("The last file has", nb_line_last_file, "line(s) it will be splitted in two"))

  ### DIVIDE IN TWO PARTS
  system(paste("zcat", tmp_concat, "| head -n", floor(as.integer(nb_line_last_file/2)), "| gzip --fast >", before_last_file))
  system(paste("zcat", tmp_concat, "| tail -n", paste("+", floor(as.integer(nb_line_last_file/2 + 1)), sep=""), "| gzip --fast >", last_file))
  system(paste("rm", tmp_concat))
}

## LOAD THE FILENAMES OF THE DIFFERENT CHUNKS
lst_files = system(paste("find",output_tmp_chunks,"-iname \"*_subfile.txt.gz\" | sort -n"), intern = TRUE)

logging("Split done")

## LOAD THE HEADER
header = as.character(unlist(read.table(file = header_kmer_counts, sep = "\t", header = FALSE)))

logging(paste("Foreach of the", length(lst_files),"files"))

## zero inflated negative binomial regression ANALYSIS ON EACH k-mer
invisible(foreach(i=1:length(lst_files)) %dopar% {
  
  RhpcBLASctl::omp_set_num_threads(1)
  RhpcBLASctl::blas_set_num_threads(1)
  
  bigTab = as.matrix(read.table(lst_files[[i]],header=F,stringsAsFactors=F,row.names=1))
  colnames(bigTab) <- header
  logicalTab <- bigTab > 0

  res <- apply(logicalTab, 1, function(presence) {
    
    isA <- colData$condition == conditionA
    isB <- colData$condition == conditionB
    alltreatment <- all(presence[isB])
    nocontrol <- !any(presence[isA])
    allcontrol <- all(presence[isA])
    notreatment <- !any(presence[isB])
  
    if(alltreatment & nocontrol) {
      return(c(1e-6,1))
    } else if (allcontrol & notreatment) {
      return(c(1e-6,-1))
    } else {
      return(c(1,0))
    }
    
  })
  
  outdf <- data.frame(ID=rownames(bigTab), 
                      pvalue=res[1,],
                      meanA=rowMeans(bigTab[,rownames(colData)[colData$condition==conditionA]]), 
                      meanB=rowMeans(bigTab[,rownames(colData)[colData$condition==conditionB]]),
                      log2FC=res[2,]) ## actually binary direction; keeping name for compatibility
                     
  outdf_filt <- outdf[,c('ID','pvalue')]
  write.table(outdf_filt,
              file=gzfile(paste(output_tmp_separ,i,"_pvalue_part_tmp.gz",sep="")),
              sep="\t",
              quote=FALSE,
              col.names = FALSE,
              row.names = FALSE)
  
  outdf_filt <- cbind(outdf[, c('ID','meanA','meanB','log2FC')], as.data.frame(bigTab))
  write.table(outdf_filt,
              file=gzfile(paste(output_tmp_separ,i,"_datasepar_part_tmp.gz", sep="")),
              sep="\t",
              quote=FALSE,
              col.names = FALSE,
              row.names = FALSE)
  
}) #END FOREACH

system(paste("rm -rf", output_tmp_chunks))

logging(paste("Foreach done:", date()))

### print a table that can be used downstream
dir.create(dirname(output_pvalue_all))

#MERGE ALL CHUNKS PVALUE INTO A FILE
system(paste("find", output_tmp_separ, "-name '*_pvalue_part_tmp.gz' | xargs cat >", output_pvalue_all))

logging(paste("Pvalues merged into",output_pvalue_all))

#MERGE ALL CHUNKS DESeq2 INTO A FILE
system(paste("find", output_tmp_separ, "-name '*_datasepar_part_tmp.gz' | xargs cat >", dataseparAll))

logging(paste("DESeq2 results merged into",dataseparAll))

# REMOVE DESeq2 CHUNKS RESULTS
system(paste("rm -rf", output_tmp_separ))
 
#LEFT JOIN INTO dataseparAll
#GET ALL THE INFORMATION (ID,MEAN_A,MEAN_B,LOG2FC,COUNTS) FOR DE KMERS
system(paste("echo \"LANG=en_EN join <(zcat ", output_pvalue_all," | LANG=en_EN sort -n -k1) <(zcat ", dataseparAll," | LANG=en_EN sort -n -k1) | awk 'function abs(x){return ((x < 0.0) ? -x : x)} {if (abs(\\$5) >=0 && \\$2 <= ", pvalue_threshold, ") print \\$0}' | tr ' ' '\t' | gzip > ", dataseparFiltered, "\" | bash", sep=""))

logging("Get counts for pvalues that passed the filter")

#CREATE THE FINAL HEADER USING ADJ_PVALUE AND dataseparAll ONES AND COMPRESS THE FILE
# CREATE THE HEADER FOR THE DESeq2 TABLE RESULT

#SAVE THE HEADER
system(paste("echo 'tag\tpvalue\tmeanA\tmeanB\tlog2FC' | paste - ", header_kmer_counts," | gzip > ", output_diff_counts))
system(paste("cat", dataseparFiltered, ">>", output_diff_counts))
system(paste("rm", dataseparFiltered))

logging("Analysis done")
