rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(DescTools)
library(here)
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


input_folder <- "p7"

sigmoid <- function(a, b, x){
  1/(1 + exp(-a*x + b))
}
## read in the nt mutation probability matrix
p_trans <- readRDS(here(input_folder, "p_trans_new"))
diag(p_trans) <- 1
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
nt_names <- c("A", "T", "C", "G")
codon_table <- c()
for (nt1 in nt_names ){
  for (nt2 in nt_names){
    for (nt3 in nt_names){
      codon_table <- c(codon_table, paste(c(nt1, nt2, nt3), collapse = ""))
    }
  }
}
AA_table <- translate(DNAStringSet(codon_table))

## obtain the predicted mutation probability for each codon at each position
p_predicted <-
  foreach(icodon = 1:length(most_common_codon), .combine = "cbind") %do% {
    codon <- most_common_codon[icodon]
    letter_codon <- str_split(codon, "")[[1]]
    p <- rep(1, length(nt_names)**3)
    ip <- 1
    for (nt1 in nt_names ){
      for (nt2 in nt_names){
        for (nt3 in nt_names){
          new_letter <- c(nt1, nt2, nt3)
          for (i in 1:3){
            p[ip] <- p[ip]*p_trans[letter_codon[i], new_letter[i]]
          }
          ip <- ip + 1
        }
      }
    }
    p
  }
dimnames(p_predicted) <- NULL

df_mutcodon_pbdist <- 
  foreach(isite = 1:length(most_common_codon), .combine = "rbind") %do% {
    p <- p_predicted[,isite]
    dfcodonp <- data.frame("codon" = codon_table, "AA" = AA_table, "p" = p)
    dfcodonp_reduced <-
      dfcodonp %>%
      filter(AA != most_common_AA_translated[isite]) %>%
      filter(AA != "*") %>%
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
    
    if (sum(dfcodonjoined$Freq) == 0){
      dfcodonjoined <- NULL
    }
    dfcodonjoined
  }

df_mutcodon_pbdist <-
  df_mutcodon_pbdist %>%
  mutate(bdist_scaled = (bdist - min(bdist))/(max(bdist) - min(bdist)))

get_corr <- function(alpha){ # maximum between 5 and 6
  
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = p*exp(alpha*bdist_scaled)/sum(p*exp(alpha*bdist_scaled))) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor
}

get_corr_bdistonly <- function(alpha){ # maximum between 5 and 6
  
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = exp(alpha*bdist_scaled)/sum(exp(alpha*bdist_scaled))) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor
}

get_corr_ponly <- function(alpha){ # maximum between 5 and 6
  
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = p**alpha/sum(p**alpha)) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor
}

get_corr_linear <- function(cab){
  a <- cab[1]
  b <- cab[2]
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = (p + a*bdist_scaled + b)/sum(p + a*bdist_scaled + b)) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor  
}

get_corr_product <- function(cab){
  a <- cab[1]
  b <- cab[2]
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = (p**b *bdist_scaled**a )/sum(p**b *bdist_scaled**a)) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")

  -AAcor  
}

get_corr_product_2 <- function(cab){
  a <- cab[1]
  b <- cab[2]
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = (p**b *exp(bdist_scaled*a) )/sum(p**b *exp(bdist_scaled*a))) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor  
}

get_corr_sigmoid <- function(cabcd){
  a <- cabcd[[1]]
  b <- cabcd[[2]]
  aa <- cabcd[[3]]
  bb <- cabcd[[4]]
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = sigmoid(a, b, p)*sigmoid(aa, bb, bdist_scaled)/sum(sigmoid(a, b, p)*sigmoid(aa, bb, bdist_scaled))) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  -AAcor  
}

optimize(get_corr, c(0,10)) 
                            
optimize(get_corr_bdistonly, c(0,20)) 

optimise(get_corr_ponly, c(0,10)) 
optim(c(1,1), get_corr_product_2) 
optim(c(1,1), get_corr_product) 
optim(c(0,0), get_corr_linear) 
optim(c(1,0,1,0), get_corr_sigmoid)

#####
# plot the best choices
#####
get_corr_product_2_for_plot <- function(cab){
  a <- cab[1]
  b <- cab[2]
  df_new <- 
    df_mutcodon_pbdist %>%
    group_by(isite) %>%
    mutate(pscaled = (p**b *exp(bdist_scaled*a) )/sum(p**b *exp(bdist_scaled*a))) %>%
    mutate(Freqscaled = Freq/sum(Freq))
  
  codoncorr <- cor(df_new$Freqscaled, df_new$pscaled)
  codoncorr_spear <- cor(df_new$Freqscaled, df_new$pscaled, method = "spearman")
  
  df_new_AA <- 
    df_new %>%
    ungroup() %>%
    group_by(isite, AA) %>%
    summarise_each(funs(sum), Freqscaled, pscaled)
  
  AAcor <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled)
  AAcor_spear <- cor(df_new_AA$Freqscaled, df_new_AA$pscaled, method = "spearman")
  
  list(df_new, df_new_AA)
  
}


df_AA <- get_corr_product_2_for_plot(c(3.0228598, 0.4581865))[[2]]
nexp <- 0.18
df_AA$Freqscaled <- df_AA$Freqscaled**nexp
df_AA$pscaled <- df_AA$pscaled**nexp
p <- ggplot(data = df_AA)
p <- p + geom_point(mapping = aes(x= pscaled, y = Freqscaled))
p <- p + coord_equal()
cust_breaks <- c(0, 10**(-5), 10**(-4), 10**(-3), 10**(-2), 10**(-1), 1)**nexp
cust_labels <- c("0", "10e-5", "10e-4", "10e-3", "10e-2", "10e-1", "1")
p <- p + scale_x_continuous(breaks = cust_breaks, 
                            labels = cust_labels,
                            name = "probability predicted",
                            limits = c(0,1))
p <- p + scale_y_continuous(breaks = cust_breaks, 
                            labels = cust_labels,
                            name = "probability observed",
                            limits = c(0,1))
p <- p + theme_classic()
p <- p + theme(axis.text.x= element_text(angle = 90))

#p <- p + ylim(0,1)
p
