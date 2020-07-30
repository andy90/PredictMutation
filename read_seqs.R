rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(here)
substrRight <- function(x){
  substr(x, nchar(x)-2, nchar(x))
}

inputs <- read.table("input.txt", stringsAsFactors = FALSE)$V1
input_folder <- inputs[1]
input_file <- inputs[2]

vp_fasta_file <- here(input_folder, input_file)
vp_nt <- toupper(read.fasta(vp_fasta_file, seqtype = "DNA", as.string = TRUE, forceDNAtolower = FALSE))
vp_nt_seqs <- gsub("-", "", unlist(vp_nt, use.names = FALSE)) # remove gaps in each sequence

# if the seqs end with stop codon, then delete the codon
ind_stop <- substrRight(vp_nt_seqs) %in% c("TAA", "TAG", "TGA")
vp_nt_seqs[ind_stop] <-  gsub('.{3}$', '', vp_nt_seqs[ind_stop])

# translate nucleotide sequences to amino acid sequences, need fuzzy due to sequencing errors
vp_aa <- translate(DNAStringSet(vp_nt_seqs, use.names = FALSE), if.fuzzy.codon = "solve")

vp <- data.frame( "AA" = vp_aa, "DNA" = vp_nt_seqs,   stringsAsFactors = FALSE)

vp_good <-
  vp %>%
  filter(!grepl(pattern = "[^ATCG]", DNA)) %>% # remove those that have unwanted letters
  filter(!grepl("\\*", AA))

saveRDS(vp_good, here(input_folder, "vpgood"))
