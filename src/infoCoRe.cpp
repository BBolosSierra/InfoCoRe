#include "Headers.h"

//----------------
//NOTE:
// When editing this file, i.e. adding or removing functions, need to carry out the following tasks:
// 1) edit file NAMESPACE, adding/removing the function name in export(...)
// Next run following in R, inside the package:
// 2) $> cd rblock/
// 3) $> R
// 4)  > library(Rcpp)
// 5)  > Rcpp::compileAttributes()
// running 5) will create new RcppExports.R in R/ and
// new RcppExports.cpp in src/
// If need to run (to build Makevars file in src):
// 1) $> cd rblock/
// 2) $> autoconfig
// 3) $> ./configure
// BUILD PACKAGE:
// 1) $> R CMD build rblock/
// INSTALL PACKAGE (note 'INSTALL' should be in capital letters):
// 1) $> R CMD INSTALL rblock_1.0.tar.gz
//----------------


//Ref Papers:
// 1) P. Villegas et al. Laplacian paths in complex networks: Information core emerges from entropic transitions, PHYSICAL REVIEW RESEARCH 4, 033196 (2022).
// doi: DOI: 10.1103/PhysRevResearch.4.033196
// website: https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.033196


// [[Rcpp::export]]
arma::Mat<double> norm_laplacian( const arma::SpMat<double>& Adj, Rcpp::IntegerVector norm=1){
  //undirected network
  
  arma::uword N, option_norm;

  option_norm = norm[0];
  
  N = Adj.n_rows;

  // Degree Matrix
  arma::Mat<double> D(diagmat(sum(Adj,0)));
  
  
  // Laplacian
  arma::Mat<double> L(N,N); L=D-Adj;
  
  if( option_norm==1 ){
  
    D.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });

    // Normalised Laplacian
    L = eye(N,N) - D * Adj * D;
    
  }

  return L;
    
}

// [[Rcpp::export]]
arma::Mat<std::complex<double>> norm_laplacian_cx( const arma::SpMat<double>& Adj,
                                                   Rcpp::IntegerVector weighted=0,
                                                   Rcpp::IntegerVector norm=1){
  
  // Refs: 1) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 2) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 3) https://www.stats.ox.ac.uk/~cucuring/Hermitian_Clustering_AISTATS.pdf  
  
  arma::uword i,j,ii,jj,N,option_we,option_norm;

  option_we   = weighted[0];
  option_norm = norm[0];
  
  N = Adj.n_rows;  
  
  arma::Mat<std::complex<double>> H; H.zeros(N,N);

  if( option_we == 0 ){
    //unweighted    
  
    for(ii=0; ii<N; ii++){
      const arma::SpSubview_row<double> rindx = Adj.row(ii);
      const arma::uvec cindx = find(rindx);
      for(jj=0; jj<cindx.n_elem; jj++){

        // Define the Adj. matrix rows and column indices.
        i = ii;
        j = cindx(jj);

        // Find edge weights
        double we_ij=Adj.at(i,j);
        double we_ji=Adj.at(j,i);

        // Build Hermitain matrix
        if( we_ij > 0 && we_ji > 0 ){
          // edge in both directions
          std::complex<double> bidir(1,0);
          H.at(i,j) = bidir;
        } else {
          if( we_ij > 0 && we_ji == 0 ){
            // out-going edge
            std::complex<double> cWe_ij(0,1);
            H.at(i,j) = cWe_ij; H.at(j,i) = std::conj(cWe_ij);
          } else { 
            // in-coming edge
            std::complex<double> cWe_ji(0,1);        
            H.at(j,i) = cWe_ji; H.at(i,j) = std::conj(cWe_ji);
          }
        }
      }
    }

  } else {
    //weighted
    
     for(ii=0; ii<N; ii++){
      const arma::SpSubview_row<double> rindx = Adj.row(ii);
      const arma::uvec cindx = find(rindx);
      for(jj=0; jj<cindx.n_elem; jj++){

        // Define the Adj. matrix rows and column indices.
        i = ii;
        j = cindx(jj);

        // Find edge weights
        double we_ij=Adj.at(i,j);
        double we_ji=Adj.at(j,i);

        // Build Hermitain matrix
        if( we_ij > 0 && we_ji > 0 ){
          // edge in both directions
          std::complex<double> bidir((we_ij-we_ji),0);
          H.at(i,j) = bidir;
        } else {
          if( we_ij > 0 && we_ji == 0 ){
            // out-going edge
            std::complex<double> cWe_ij(0,we_ij);
            H.at(i,j) = cWe_ij; H.at(j,i) = std::conj(cWe_ij);
          } else { 
            // in-coming edge
            std::complex<double> cWe_ji(0,we_ji);        
            H.at(j,i) = cWe_ji; H.at(i,j) = std::conj(cWe_ji);
          }
        }
      }
    }
     
  }

  // Degree Matrix, using in- and out-degree
  arma::SpRow<double> d = sum(Adj,0) + sum(Adj,1).t();
  // Now cast from Sparse matrix 'd' to dense matrix D 
  arma::Mat<double> D(diagmat(d));

  // Laplacian
  arma::Mat<std::complex<double>> L(N,N); L=D-H;
  
  if( option_norm==1 ){
    
    // Normalised Laplacian
    D.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });
  
    L = eye(N,N) - D * H * D;

  }

  return L;
  

}


// [[Rcpp::export]]
void driver( const arma::SpMat<double>& Adj,
             Rcpp::IntegerVector weighted=0,
             Rcpp::IntegerVector directed=0,
             Rcpp::IntegerVector norm=1){

  arma::uword option_dir,option_we,option_norm;
  option_dir  = directed[0];
  option_we   = weighted[0];
  option_norm = norm[0];
  
  arma::Mat<double> L;
  arma::Mat<std::complex<double>> L_dir;
   
  if( option_dir == 0 ){
    cout << "> undirected Adj: " << endl;
    cout << "> calculate L... ";
    L = norm_laplacian(Adj,norm=option_norm);
    L.brief_print("L:");
  } else {
    cout << "> directed Adj: " << endl;
    cout << "> calculate L... ";
    L_dir = norm_laplacian_cx(Adj,weighted=option_we,norm=option_norm);
    L_dir.brief_print("L:");
  }
  
}
