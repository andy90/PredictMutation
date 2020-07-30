rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(msa)
library(here)
#library(doParallel)
#registerDoParallel(cores=16)

inputs <- read.table("input.txt", stringsAsFactors = FALSE)$V1
input_folder <- inputs[1]

vp_good <- readRDS(here(input_folder, "vpgood"))

vp_dna_msa <- readRDS(here(input_folder, "DNAmsa"))
rownames(vp_dna_msa) <- NULL

vp_msa <-  do.call(rbind, read.fasta(here(input_folder, "vp_aligned.fasta")))
rownames(vp_msa) <-  NULL

gap_percent <- colSums(vp_msa == "-")/nrow(vp_msa)
sum(gap_percent > 0.125) # same criteria used by Ferguson
gap_position <- which(gap_percent > 0.125)

vp_msa_reduced <- vp_msa[, -gap_position]
dna_gap_position <- c( 3 * (gap_position - 1) + 1, 3 * (gap_position - 1) + 2, 3 * gap_position)
vp_dna_msa_reduced <- vp_dna_msa[, -dna_gap_position]

vp_dna_msa_collapsed <- 
  sapply(1:ncol(vp_msa_reduced), function(j){
    jbegin <- (3*(j-1) + 1)
    jmiddle <- 3*j - 1
    jend <- 3*j
    a <- paste(vp_dna_msa_reduced[, jbegin], vp_dna_msa_reduced[, jmiddle], vp_dna_msa_reduced[, jend], sep = "")
    a
    
  })

# also remove positions with no mutations
codon_nmuts <-
  apply(vp_dna_msa_collapsed, 2, function(seq){
    a_table <- sort(table(seq), decreasing = TRUE)
    length(a_table)
  })

ind_muts <- which(codon_nmuts != 1)

saveRDS(vp_dna_msa_collapsed[, ind_muts], here(input_folder, "vp_dna_msa_collapsed"))
saveRDS(vp_msa_reduced[, ind_muts], here(input_folder, "vp_msa_reduced"))