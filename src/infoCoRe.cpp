#include "Headers.h"

//Ref Papers:
// 1) P. Villegas et al. Laplacian paths in complex networks: Information core emerges from entropic transitions, PHYSICAL REVIEW RESEARCH 4, 033196 (2022).
// doi: DOI: 10.1103/PhysRevResearch.4.033196
// website: https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.033196


// Numerical differentiation to calculate C(τ)
// [[Rcpp::export]]
void calculate_specific_heat(const arma::vec& tau_vec,
                             const arma::vec& ent_vector,
                             arma::vec& specific_heat,
                             arma::uvec& heat_peak_indices,
                             arma::vec& tau_peak_values) {
  // Check if the lengths of tau_vec and entropy_vector match
  if (tau_vec.n_elem != ent_vector.n_elem) {
    Rcpp::stop("Lengths of tau_vec and entropy_vector must be the same.");
  }
  // Calculate finite differences to approximate the derivative
  specific_heat.set_size(tau_vec.n_elem);
  for (arma::uword i = 1; i < tau_vec.n_elem; ++i) {
    double delta_tau = std::log(tau_vec[i]) - std::log(tau_vec[i - 1]);
    double delta_entropy = ent_vector[i] - ent_vector[i - 1];
    
    // Use finite difference formula to approximate the derivative
    specific_heat[i] = - (delta_entropy / delta_tau);
  }
  // The first element is left undefined, set it to zero or some default value
  specific_heat[0] = 0.0;
  
  // Detect peaks
  for (arma::uword i = 1; i < specific_heat.n_elem - 1; ++i) {
    if (specific_heat[i] > specific_heat[i - 1] && specific_heat[i] > specific_heat[i + 1]) {
      // i is a peak index
      heat_peak_indices.resize(heat_peak_indices.n_elem + 1);
      heat_peak_indices(heat_peak_indices.n_elem - 1) = i;
      
      // Store corresponding tau value
      tau_peak_values.resize(tau_peak_values.n_elem + 1);
      tau_peak_values(tau_peak_values.n_elem - 1) = tau_vec[i];
    }
  }
}

// [[Rcpp::export]]
Rcpp::List renormalize_graph(const arma::Mat<double>& x,
                            const arma::vec& tau_peak_values) {
  
  // If there is only one relevant tau, perform additional logic
  //if (tau_peak_values.size() == 1) {
  double tau = tau_peak_values(0);
  
  // Perform the logic from code
  arma::mat L1 = x; // Use the input matrix directly
  arma::mat K = exp(-tau * L1);
  double tr = arma::trace(K);
  arma::mat rho = K / tr;
  arma::mat adj2 = arma::zeros<arma::mat>(rho.n_rows, rho.n_cols);
  
  for (size_t i = 0; i < rho.n_rows; ++i) {
    for (size_t j = 0; j < rho.n_cols; ++j) {
      if (rho(i, j) >= rho(j, j) || rho(i, j) >= rho(i, i)) {
        adj2(i, j) = 1;
      }
    }
  }
  // Create a list
  Rcpp::List renorm;
  renorm["K"] = K;
  renorm["adj2"] = adj2;
  renorm["tr"] = tr;
  renorm["rho"] = rho;
  
  // Return the list
  return renorm;
}

// Function to calculate entropy
// Input: mu
// Output: ent - the calculated entropy value
// [[Rcpp::export]]
void entropy( const arma::vec& mu,
              double& ent ){
  
  // Declare variables
  arma::uword i;
  double N = mu.n_elem; // Number of elements in the vector mu
  double sum=0;
  arma::vec tmp; tmp.zeros(N); // Initialize a vector of zeros with the same size as mu

  // Calculate the sum of mu[i] * log(mu[i]) for each element in mu
  for( i=0; i<N; i++ ){
    sum += mu[i] * std::log(mu[i]);
  }
  
  // Calculate the entropy
  ent      = -sum/std::log(N);
  
}


// Variables:
//   lambda: Vector of input values
//   mu: Vector to store the output
//   tau: Parameter controlling the influence of the exponential operation
//        Default value is 1
// [[Rcpp::export]]
void net_operator( const arma::vec& lambda,
                   arma::vec& mu,
                   double tau=1 ){
  
  arma::uword i;
  double N    = lambda.n_elem;
  double norm = 0;
  double val;
  arma::vec tmp; tmp.zeros(N);
  

  for( i=0; i<N; i++ ){
    val    = std::exp(-1*lambda[i]*tau);
    tmp[i] = val; 
    norm  += val;
  }
    
  mu = tmp/norm;
}

// Function to calculate entropy for each column in a matrix of probability distributions
// Input: mu_matrix
// Output: ent_vector - a vector containing the calculated entropy values for each column
// [[Rcpp::export]]
void entropy_matrix(const arma::mat& mu_matrix, arma::vec& ent_vector) {
  
  arma::uword num_cols = mu_matrix.n_cols;
  ent_vector.set_size(num_cols);
  
  for (arma::uword j = 0; j < num_cols; ++j) {
    const arma::vec& mu = mu_matrix.col(j);
    
    // Declare variables
    double N = mu.n_elem; // Number of elements in the vector mu
    double sum = 0;
    arma::vec tmp;
    tmp.zeros(N); // Initialize a vector of zeros with the same size as mu
    
    // Calculate the sum of mu[i] * log(mu[i]) for each element in mu
    for (arma::uword i = 0; i < N; i++) {
      if (mu[i] > 0) {
        sum += mu[i] * std::log(mu[i]);
      } else {
        sum += 0;
      }
    }
    // Calculate the entropy
    //ent_vector(j) = -sum / std::log(N);
    // Calculate the entropy
    ent_vector(j) = (mu.n_elem > 0) ? -sum / std::log(N) : 0;
  }
}

// [[Rcpp::export]]
void net_operator_mat(const arma::vec& lambda,
                        arma::mat& mu_matrix,
                        const arma::vec& tau_vec){
  arma::uword i, j;
  double N = lambda.n_elem;
  double norm;
  double val;
  arma::vec tmp;
  mu_matrix.set_size(N, tau_vec.n_elem);
  
  for (j = 0; j < tau_vec.n_elem; ++j) {
    norm = 0;
    tmp.zeros(N);
    
    for (i = 0; i < N; i++) {
      val = std::exp(-1 * lambda[i] * tau_vec[j]);
      tmp[i] = val;
      norm += val;
    }
    mu_matrix.col(j) = tmp / norm;
  }
}

// Function to calculate the eigenvalues and/or eigenvectors of a dense matrix 'x'
// Input:
//   x: Input matrix
//   eigval: Vector to store the eigenvalues
//   eigvec: Matrix to store the eigenvectors
//   val_only: Optional parameter indicating whether to compute only eigenvalues (1) or both eigenvalues and eigenvectors (0)
//             Default is 0 (compute both)
//   order: Optional parameter indicating whether to order the eigenvalues in ascending order (0) or descending order (1)
//          Default is 1 (descending order)
// [[Rcpp::export]]
void get_eig( const arma::Mat<double>& x,
              arma::vec& eigval,
              arma::Mat<double>& eigvec,
              Rcpp::IntegerVector val_only=0,
              Rcpp::IntegerVector order=1 ){

  arma::uword i,option_valOnly,option_ord;

  option_valOnly = val_only[0];
  option_ord     = order[0];

  // Calculate the eigenvalues & vectors of dense matrix 'x'
  if( option_valOnly == 1 ){
    // Only eigenvalues
    arma::eig_sym(eigval, x);
  } else {
    // Both eigenvalues and eigenvectors
    arma::eig_sym(eigval, eigvec, x);
  }
  
  // Ordering of eigenvalues and coefficients
  if( option_ord == 1 ){
  
    // Swap the eigenvalues since they are ordered backwards (we need largest
    // to smallest).
    for( i=0; i<floor(eigval.n_elem / 2.0); ++i){
      eigval.swap_rows(i, (eigval.n_elem - 1) - i);
    }
    
    // Flip the coefficients to produce the same effect.
    if( option_valOnly == 0 ){
      eigvec = arma::fliplr(eigvec);
      
    }
  }
}



// [[Rcpp::export]]
void get_eig_cx( const arma::Mat<std::complex<double>>& x,
                 arma::vec& eigval,
                 arma::Mat<std::complex<double>>& eigvec,
                 Rcpp::IntegerVector val_only=0,
                 Rcpp::IntegerVector order=1 ){

  arma::uword i,option_valOnly,option_ord;

  option_valOnly = val_only[0];
  option_ord     = order[0];

  // Calculate the eigenvalues & vectors of dense complex matrix 'x'
  if( option_valOnly == 1 ){
    arma::eig_sym(eigval, x);
    //arma::eig_gen(eigval, x);
  } else {
    arma::eig_sym(eigval, eigvec, x); 
  }
   
  if( option_ord == 1 ){
  
    // Swap the eigenvalues since they are ordered backwards (we need largest
    // to smallest).
    for( i=0; i<floor(eigval.n_elem / 2.0); ++i){
      eigval.swap_rows(i, (eigval.n_elem - 1) - i);
    }
  
    // Flip the coefficients to produce the same effect.
    if( option_valOnly == 0 ){
      eigvec = arma::fliplr(eigvec);
    }
  }
  
}

// [[Rcpp::export]]
arma::Mat<double> laplacian( const arma::SpMat<double>& Adj, Rcpp::IntegerVector norm=1){
  //undirected network
  
  arma::uword N, option_norm;

  option_norm = norm[0];
  
  N = Adj.n_rows;

  // Degree Matrix
  arma::SpMat<double> D = diagmat(sum(Adj,0));
  
  // Laplacian
  arma::SpMat<double> Ltmp(N,N); Ltmp=D-Adj;
  
  if( option_norm==1 ){
  
    D.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });

    // Normalised Laplacian
    Ltmp = eye(N,N) - D * Adj * D;
    
  }

  // Cast from Sparse Matrix to Dense Matrix
  arma::Mat<double> L(Ltmp);
  
  return L;
    
}


// [[Rcpp::export]]
arma::Mat<std::complex<double>> laplacian_cx( const arma::SpMat<double>& Adj,
                                              Rcpp::IntegerVector weighted=0,
                                              Rcpp::IntegerVector norm=1){
  
  // Refs: 1) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 2) https://stackoverflow.com/questions/67189074/hermitian-adjacency-matrix-of-digraph
  // Refs: 3) https://www.stats.ox.ac.uk/~cucuring/Hermitian_Clustering_AISTATS.pdf  
  
  arma::uword i,j,ii,jj,N,option_we,option_norm;

  option_we   = weighted[0];
  option_norm = norm[0];
  
  N = Adj.n_rows;  
  
  arma::SpMat<std::complex<double>> H; H.zeros(N,N);

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
  arma::SpMat<double> D = diagmat(d);

  // Laplacian
  arma::SpMat<std::complex<double>> Ltmp(N,N); Ltmp=D-H;
  
  if( option_norm==1 ){
    
    // Normalised Laplacian
    D.transform(  [](double val) { return ( (val==0 ? val : pow(val,-0.5)) ); });
  
    Ltmp = eye(N,N) - D * H * D;

  }

  // Cast from Sparse Matrix to Dense Matrix   
  arma::Mat<std::complex<double>> L(Ltmp);
  
  return L;
  

}


// [[Rcpp::export]]
Rcpp::List driver( const arma::SpMat<double>& Adj,
             Rcpp::IntegerVector weighted=0,
             Rcpp::IntegerVector directed=0,
             Rcpp::IntegerVector norm=1,
             Rcpp::IntegerVector val_only=0,
             Rcpp::IntegerVector order=1,
             const arma::vec& custom_tau_vec = Rcpp::NumericVector::create(1.0, 10, 50, 100)){

  arma::uword option_dir,option_we,option_norm,option_valOnly,option_ord;
  option_dir     = directed[0];
  option_we      = weighted[0];
  option_norm    = norm[0];
  option_ord     = order[0];
  option_valOnly = val_only[0];
  
  arma::Mat<double> L;
  arma::Mat<std::complex<double>> L_dir;
  arma::vec tau_vec = custom_tau_vec;
  
  arma::Mat<double> eigvec;
  arma::vec         eigval;
  arma::vec         mu; // for mu = 1 (net_operator)
  arma::mat         mu_matrix; // Matrix for get_transition_tau
  arma::vec         ent_vector;
  arma::vec         specific_heat;
  arma::uvec        heat_peak_indices; 
  arma::vec         tau_peak_values;
  double            ent=-1;
  double            tau=1;
  
  arma::Mat<std::complex<double>> cx_eigvec;
  arma::vec         cx_eigval;
  
  Rcpp::List result;
  
  if( option_dir == 0 ){
    cout << "> undirected Adj: " << endl;
    
    cout << "> calculate L... ";
      
    L = laplacian(Adj, norm=option_norm);
    
    cout << "done!" << endl;
    
    L.brief_print("L:");
    
    result["L"] = Rcpp::wrap(L);
    
    cout << "> TEST: " << endl;
    
    get_eig(L, eigval, eigvec, val_only=option_valOnly, order=option_ord);

    eigval.brief_print("eigval:");

    eigvec.brief_print("eigvec:");

    result["eigval"] = Rcpp::wrap(eigval);
    
    result["eigvec"] = Rcpp::wrap(eigvec);
    
    cout << "> TEST2: " << endl;
    
    net_operator(eigval, mu, tau = 1);
    
    mu.brief_print("mu (net_operator):");
    
    result["mu"] = Rcpp::wrap(mu);
    
    cout << "> TEST3: " << endl;
    
    entropy(mu, ent);
    
    result["entropy"] = Rcpp::wrap(ent);
    
    cout << "> entropy: " << ent << endl;

    cout << "> TEST4: " << endl;
    
    net_operator_mat(eigval, mu_matrix, tau_vec);
    
    mu_matrix.brief_print("mu (get_transition_tau):");
    
    tau_vec.brief_print("taus:");
    
    result["mu_matrix"] = Rcpp::wrap(mu_matrix);
    
    result["tau_vec"] = Rcpp::wrap(tau_vec);
    
    cout << "> TEST5: " << endl;
    
    entropy_matrix(mu_matrix, ent_vector);
    
    ent_vector.brief_print("Entropy values per tau:");
    
    result["ent_vector"] = Rcpp::wrap(ent_vector);
    
    cout << "> TEST6: " << endl;
    
    calculate_specific_heat(tau_vec, ent_vector, specific_heat, heat_peak_indices, tau_peak_values);
    
    result["specific_heat"] = Rcpp::wrap(specific_heat);
    
    result["tau_peak_values"] = Rcpp::wrap(tau_peak_values);
    
    specific_heat.brief_print("Specific Heat values per tau:");
    
    tau_peak_values.brief_print("tau peak values:");
    
    cout << "> TEST7: " << endl;
    
    Rcpp::List renorm = renormalize_graph(L, tau_peak_values);
    
    result["K"] = Rcpp::wrap(renorm["K"]);
    result["adj2"] = Rcpp::wrap(renorm["adj2"]);
    result["tr"] = Rcpp::wrap(renorm["tr"]);
    result["rho"] = Rcpp::wrap(renorm["rho"]);
    
    // printing the matrix
    arma::mat K = renorm["K"];
    arma::mat adj2 = renorm["adj2"];
    double tr = renorm["tr"];
    arma::mat rho = renorm["rho"];
    
    K.brief_print("K: ");
    
    adj2.brief_print("adj2: ");
    
    Rcpp::Rcout << "tr: " << tr << std::endl;
    
    rho.brief_print("rho: ");
    
    cout << "> Done!! " << endl;
    
  } else {
    cout << "> directed Adj: " << endl;
    cout << "> calculate L... ";
    L_dir = laplacian_cx(Adj, weighted=option_we, norm=option_norm);
    cout << "done!" << endl;
    L_dir.brief_print("L:");

    cout << "> is Hermitian: " << L_dir.is_hermitian() << endl;
    
    cout << "> TEST: " << endl;
    get_eig_cx(L_dir, cx_eigval, cx_eigvec, val_only=option_valOnly, order=option_ord);
    
    cx_eigval.brief_print("cx_eigval:");

    cx_eigvec.brief_print("cx_eigvec:");

  }
  
  return Rcpp::wrap(result);
}

