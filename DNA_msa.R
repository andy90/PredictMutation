rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(here)
#library(doParallel)
#registerDoParallel(cores=16)
library(Rcpp)
sourceCpp("alignDNA.cpp")

inputs <- read.table("input.txt", stringsAsFactors = FALSE)$V1
input_folder <- inputs[1]

vp_good <- readRDS(here(input_folder, "vpgood"))
vp_msa <- read.fasta(here(input_folder, "vp_aligned.fasta"))

vp_dna <- str_split(vp_good$DNA, "")

a <- sapply(1:length(vp_msa), function(i){
  AAmsa <- vp_msa[[i]]
  DNAseq <- vp_dna[[i]]
  DNAmsa <- alignDNA(AAmsa, DNAseq)
})

a <- t(a)

saveRDS(a, here(input_folder, "DNAmsa"))