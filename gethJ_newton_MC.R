rm(list = ls())
library(tidyverse)
library(seqinr)
library(foreach)
library(Biostrings)
library(msa)
library(DescTools)
library(Rcpp)
sourceCpp("MC/ising_rcpp.cpp")

p17_msa_reduced <- readRDS( "p17_msa_reduced")
dimnames(p17_msa_reduced) <- NULL

most_common_AA <- 
  apply(p17_msa_reduced, 2, function(seq){
    a <- names(sort(table(seq), decreasing = TRUE))[1] 
    a
  })

AA_mut <- sweep(p17_msa_reduced, 2, most_common_AA, FUN = "!=")
colSums(AA_mut)

nonsynony <- readRDS("nonsynony")
colSums(nonsynony)
which(colSums(AA_mut) != colSums(nonsynony)) # at this two positions, the most common AA does not corresponds to the most common codon

pi <- colSums(AA_mut)/nrow(AA_mut) 
pij <- t(AA_mut) %*% AA_mut / nrow(AA_mut) # this is our target

nsite <- length(pi)
h <- rep(0, nsite)
J <- matrix(0, nrow = nsite, ncol = nsite)

gammah <- 0.05 # if I use 0.5, the result will be really good for those not so small mutation probs
# for those small mutation probs I will just get 0 directly
gammaJ <- 0
for (i in 1:100){
  a <- get_distribution(J, h)
  pobi <- a$sav
  pobij <- a$scoupleav
  
  deltah <- -gammah * (pobi - pi)/(pi*(pi-1)  +  .Machine$double.xmin)
  deltaJ <- -gammaJ * (pobij - pij)/(pij*(pij-1) + 10**(-5))
  
  h <- h + deltah
  J <- J + deltaJ
  
  print(pobi)
}

cor(pobi, pi)
