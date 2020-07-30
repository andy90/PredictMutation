rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(DescTools)
library(expm)
data("BLOSUM62")
data("BLOSUM80")
data("BLOSUM45")
data("BLOSUM100")
data("PAM30")
data("PAM40")
data("PAM70")
data("PAM120")
data("PAM250")
AASimilarity <- BLOSUM100 #readRDS("sandberg_custom") # both sandberg and grantham do not work as well as Blosum

input_folder <- "HA"

p_trans <- readRDS(here(input_folder, "p_trans_new"))
diag(p_trans) <- diag(p_trans) -1


## read in the msa for the amino acid sequence, find the frequency of each
## amino acid at each position; find the most common AA at each position
vp_msa_reduced <- readRDS(here(input_folder, "vp_msa_reduced"))
dimnames(vp_msa_reduced) <- NULL

most_common_AA <- 
  apply(vp_msa_reduced, 2, function(seq){
    a <- names(sort(table(seq), decreasing = TRUE))[1] 
    a
  })

AA_frequecies <- 
  apply(vp_msa_reduced, 2, function(seq){
    a_table <- sort(table(seq), decreasing = TRUE)
    aa_freq <- a_table/sum(a_table)
    data.frame(aa_freq)
  })

## find the dna sequence of the msa aligned sequences
## find the most common codon
## find the frequency of each codon at each position
## find the translated amino acid from the most common codon
vp_dna_msa_collapsed <- readRDS(here(input_folder, "vp_dna_msa_collapsed"))
dimnames(vp_dna_msa_collapsed) <- NULL

most_common_codon <- 
  apply(vp_dna_msa_collapsed, 2, function(seq){
    a <- names(sort(table(seq), decreasing = TRUE))[1] 
    a
  })
most_common_AA_translated <- translate(DNAStringSet(most_common_codon))

codon_frequencies <-
  apply(vp_dna_msa_collapsed, 2, function(seq){
    a_table <- sort(table(seq), decreasing = TRUE)
    aa_freq <- a_table/sum(a_table)
    data.frame(aa_freq)
  })

## construct the codon table
## obtain the AA table corresponding to the codon table
nt <- c(1, 2, 3, 4)
nt_names <- c("A", "T", "C", "G")
names(nt) <- nt_names
codon_table <- c()
for (nt1 in nt_names ){
  for (nt2 in nt_names){
    for (nt3 in nt_names){
      codon_table <- c(codon_table, paste(c(nt1, nt2, nt3), collapse = ""))
    }
  }
}
AA_table <- translate(DNAStringSet(codon_table))

get_mutations <- function(codons, codon1){
  lapply(codons, function(codon){
    letter1 <- str_split(codon, "")[[1]]
    letter2 <- str_split(codon1, "")[[1]]
    unname(nt[unname(rbind(letter2, letter1))])
  })
  
}
df_mutcodon_pbdist <- 
  foreach(isite = 1:length(most_common_codon), .combine = "rbind") %do% {
    dfcodonp <- data.frame("codon" = codon_table, "AA" = AA_table)
    dfcodonp_reduced <-
      dfcodonp %>%
      filter(AA != most_common_AA_translated[isite]) %>%
      filter(AA != "*") %>%
      mutate(mutations = get_mutations(codon, most_common_codon[isite])) %>%
      mutate(bdist = sapply(AA, function(ai){
        dist <- AASimilarity[ai, as.character(most_common_AA_translated[isite])]
      })) %>%
      mutate(isite = isite)
    
    dfcodon_exp <- codon_frequencies[[isite]]
    dfcodon_exp_formated <- 
      dfcodon_exp %>%
      filter(seq != "---") %>%
      set_names("codon", "Freq") %>% 
      select(codon, Freq)
    
    dfcodonjoined <- left_join(dfcodonp_reduced, dfcodon_exp_formated, by = "codon")
    dfcodonjoined[is.na(dfcodonjoined)] <- 0
    dfcodonjoined <- 
      dfcodonjoined %>%
      mutate(Freqscaled = Freq/sum(Freq))
    
    dfcodonjoined
  }



get_corr_product_2 <- function(cab){ 
  
  a <- cab[1]
  b <- cab[2]
  p_trans <- expm(p_trans * b)
  df_new <- 
    df_mutcodon_pbdist %>%
    mutate(p = sapply(mutations, function(mutation){
      p_trans[mutation[1], mutation[2]]*p_trans[mutation[3], mutation[4]]*p_trans[mutation[5], mutation[6]]
    })) %>%
    group_by(isite) %>%
    mutate(pscaled = (p *exp(bdist*a) )/sum(p *exp(bdist*a)))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  # codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  # 
  df_new_AA <-
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_at(c("Freqscaled", "pscaled"), sum)

  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  # AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor
}

get_corr_product_2_othervalues <- function(cab){ 
  
  a <- cab[1]
  b <- cab[2]
  p_trans <- expm(p_trans * b)
  df_new <- 
    df_mutcodon_pbdist %>%
    mutate(p = sapply(mutations, function(mutation){
      p_trans[mutation[1], mutation[2]]*p_trans[mutation[3], mutation[4]]*p_trans[mutation[5], mutation[6]]
    })) %>%
    group_by(isite) %>%
    mutate(pscaled = (p *exp(bdist*a) )/sum(p *exp(bdist*a)))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  
  df_new_AA <-
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_at(c("Freqscaled", "pscaled"), sum)

  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  
  list(df_new, df_new_AA, AAcor)
}

res <- optim(c(0.2,0.2), get_corr_product_2)
print(res$value)
cab_max <- res$par
res_max <- get_corr_product_2_othervalues(cab_max)
