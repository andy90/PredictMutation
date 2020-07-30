# get the frequency of mutations at nucleotide level, ie, A -> G, A -> T
# include self to self transition probability
rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(DescTools)
library(here)

inputs <- read.table("input.txt", stringsAsFactors = FALSE)$V1
input_folder <- inputs[1]

vp_dna_msa_collapsed <- readRDS(here(input_folder, "vp_dna_msa_collapsed"))
dimnames(vp_dna_msa_collapsed) <- NULL

most_common_codon <- 
  apply(vp_dna_msa_collapsed, 2, function(seq){
    a <- names(sort(table(seq), decreasing = TRUE))[1] 
    a
  })

codon_frequencies <-
  apply(vp_dna_msa_collapsed, 2, function(seq){
    a_table <- sort(table(seq), decreasing = TRUE)
    aa_freq <- a_table/sum(a_table)
    data.frame(seq = names(aa_freq), Freq = as.numeric(aa_freq))
  })



get_mutations <- function(codon1, codon2){
  letter1 <- str_split(codon1, "")[[1]]
  letter2 <- str_split(codon2, "")[[1]]
  rbind(letter1, letter2)
}

all_mut_prob <-
  foreach(i = 1:length(codon_frequencies), .combine = "rbind") %do% {
    mut_prob <- matrix(data = 0, nrow = 5, ncol = 5)
    colnames(mut_prob) <- c("A", "T", "C", "G", "-")
    rownames(mut_prob) <- c("A", "T", "C", "G", "-")
    
    df_mut <- codon_frequencies[[i]]
    for (j in 1:nrow(df_mut)){
      muts <- get_mutations(df_mut[1, 1], df_mut[j,1])
      
      for (k in 1:ncol(muts)){
        mut_prob[muts[1,k], muts[2,k]] <- mut_prob[muts[1,k], muts[2,k]] + df_mut[j, 2]
      }
    }
    mut_prob
  }

reduced_mut_prob <- matrix(data = 0, nrow = 5, ncol = 5)
colnames(reduced_mut_prob) <- c("A", "T", "C", "G", "-")
rownames(reduced_mut_prob) <- c("A", "T", "C", "G", "-")

for( i in 1:(nrow(all_mut_prob)/5)){
  ind_begin <- (i-1)*5 + 1
  ind_end <- 5*i
  reduced_mut_prob <- reduced_mut_prob + all_mut_prob[ind_begin:ind_end, ]
}

reduced_mut_prob # this should be recovered in the new model
rowSums(reduced_mut_prob)
p_trans <- sweep(reduced_mut_prob[-5,-5], 1, rowSums(reduced_mut_prob[-5,-5]), FUN = "/")
saveRDS(p_trans, file = here(input_folder, "p_trans_new"))
