rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(here)

inputs <- read.table("input.txt", stringsAsFactors = FALSE)$V1
input_folder <- inputs[1]

count_gaps <- function(subseqs){
  str_count(subseqs, "-")
}

vp_good <- readRDS(here(input_folder, "vpgood"))

write.fasta(names = 1:nrow(vp_good), sequences = as.list(vp_good$AA), file.out = here(input_folder, "vpgood.fasta"))

#system("bash doMSA.sh") # if run locally
system(paste("/home/anggao/mafft/bin/mafft --thread 15 --auto", here(input_folder,"vpgood.fasta"), " > ", here(input_folder, "vp_aligned.fasta"), sep = " ")) # if run at eofe cluster
