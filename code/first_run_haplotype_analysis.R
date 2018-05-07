#RAD-seq analysis
library(dbplyr)
library(tidyverse)

setwd("/Users/elizabethpermina/Documents/3_RADseq/code")
# firsst run haplotypes - illustrative data
hap_file<-"/Users/elizabethpermina/Documents/3_RADseq/scratch/haplotype/populations.haplotypes.tsv"
first_run_haplotypes<-read_tsv(hap_file)
(first_run_haplotypes$P1_A02.single)
