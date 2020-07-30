#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
CharacterVector alignDNA(CharacterVector AAmsa, CharacterVector DNAseq) {
  int n = AAmsa.size();
  
  CharacterVector DNAmsa(3*n);
  int ind = 0;
  for (int j = 0; j < n; j++){
    int jbegin = 3*j;
    if (AAmsa[j] != "-"){
      int ind_begin = 3*ind;
      DNAmsa[jbegin] = DNAseq[ind_begin];
      DNAmsa[jbegin+1] = DNAseq[ind_begin+1];
      DNAmsa[jbegin+2] = DNAseq[ind_begin+2];
      ind = ind + 1;
    }else{
      DNAmsa[jbegin] = "-";
      DNAmsa[jbegin+1] = "-";
      DNAmsa[jbegin+2] = "-";
    }
  }
  return DNAmsa;
}


