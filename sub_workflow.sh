#!/bin/bash  
#SBATCH -N 1  
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --constraint=centos7
#SBATCH -p sched_mit_arupc_long

Rscript read_seqs.R
Rscript AA_msa.R
Rscript DNA_msa.R
Rscript remove_gaps.R
Rscript nt_trans_freq.R
#Rscript test_model.R
