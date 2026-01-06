

Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> Hess_logpxi(const Eigen::SparseMatrix<double> &Qv,
                                        const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv) {
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;
  
  return result;
}
')





Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> Hess_logpxi_NB(const Eigen::SparseMatrix<double> &Qv,
                                          const Eigen::SparseMatrix<double> & Wnb,
                                        const Eigen::SparseMatrix<double> &MVnb,
                                        const Eigen::SparseMatrix<double> &Xv) {
                                     
   // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;
 
  return result;
 
}
')




Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::MatrixXd solve_sparse_cholesky(const Eigen::SparseMatrix<double> &Qv,
                                        const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv,
                                     const Eigen::VectorXd &b) {
                                     
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;
  
  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);
 
}
')




Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::MatrixXd solve_sparse_cholesky_NB(const Eigen::SparseMatrix<double> &Qv,
                                          const Eigen::SparseMatrix<double> & Wnb,
                                        const Eigen::SparseMatrix<double> &MVnb,
                                        const Eigen::SparseMatrix<double> &Xv,
                                     const Eigen::VectorXd &b) {
                                     
   // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;
  
  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;
  
  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);
 
}
')


Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double determinant_Pv(const Eigen::SparseMatrix<double> &Pv) {
     // Perform LDLT decomposition
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(Pv);

  if (ldlt.info() != Eigen::Success) {
    stop("LDLT decomposition failed");
  }

  // log-determinant = sum(log(diagonal of D))
  Eigen::VectorXd D_diag = ldlt.vectorD();  // for SimplicialLDLT, this is the diagonal of D
  if ((D_diag.array() <= 0).any()) {
    stop("Non-positive values in D; log-determinant is not defined");
  }

  double logdet = D_diag.array().log().sum();
  return logdet;

}
')





Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> XWX_func (const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv) {
                                     
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;
  
  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }
  
  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;
  
  // Result = Xt_W_Xv 
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  
  return result;
}
')


Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;

  // [[Rcpp::depends(RcppEigen)]]

  // [[Rcpp::export]]
  
  
  Eigen::SparseMatrix<double> XWX_func_NB (const Eigen::SparseMatrix<double> &Xv,
                  const Eigen::SparseMatrix<double> &Wnb) {
                                     
    // Construct matrix
    Eigen::SparseMatrix<double> Pv = Xv.transpose() * Wnb*Xv;

  
  return Pv;
}
')



Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd Predict_se(const Eigen::SparseMatrix<double> &Prec,
                                      const Eigen::MatrixXd &Xpred,
                                      const Eigen::VectorXi &ind) {
                                      
  const int n = Xpred.rows();
  const int full_dim = Prec.rows();
  Eigen::VectorXd se(n);
  
  
  // Cholesky factorization
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    stop("Cholesky decomposition failed");
  }
  
  #pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    // Zero-padded vector of full_dim
    Eigen::VectorXd xfull = Eigen::VectorXd::Zero(full_dim);
    for (int j = 0; j < ind.size(); ++j) {
      xfull[ind[j]-1] = Xpred(i, j);  // put row i of Xpred into full-size vector
    }
    
    Eigen::VectorXd zi = chol.solve(xfull);
    se[i] = std::sqrt(xfull.dot(zi));
  }
  
  return se;
}

')




Rcpp::sourceCpp(code = '
  #include <RcppEigen.h>
  using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppEigen.h>
#include <omp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Eigen::MatrixXd Predict_se_all(
    const Eigen::SparseMatrix<double> &Prec,
    const Eigen::MatrixXd &Xpred,
    const Rcpp::IntegerVector &ind_cb,
    const Rcpp::List &ind_cb_spat_list) {

  const int n_areas = ind_cb_spat_list.size();
  const int n_preds = Xpred.rows();
  const int full_dim = Prec.rows();

  Eigen::MatrixXd se_matrix(n_areas, n_preds);

  // Cholesky decomposition (once!)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed");
  }

  for (int a = 0; a < n_areas; ++a) {
    Rcpp::IntegerVector ind_cb_spat = ind_cb_spat_list[a];

    // Combine indices
    std::vector<int> ind_combined(ind_cb.size() + ind_cb_spat.size());
    std::copy(ind_cb.begin(), ind_cb.end(), ind_combined.begin());
    std::copy(ind_cb_spat.begin(), ind_cb_spat.end(), ind_combined.begin() + ind_cb.size());

    for (int i = 0; i < n_preds; ++i) {
      Eigen::VectorXd xfull = Eigen::VectorXd::Zero(full_dim);
      for (int j = 0; j < ind_combined.size(); ++j) {
        xfull[ind_combined[j]-1] = Xpred(i, j);
      }

      Eigen::VectorXd zi = chol.solve(xfull);
      se_matrix(a, i) = std::sqrt(xfull.dot(zi));
    }
  }

  return se_matrix;
}

')



Rcpp::sourceCpp(code = '
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
double Hutchinson_pD(
     const Eigen::SparseMatrix<double> &Ipost,
    const Eigen::SparseMatrix<double> &Ilike,
    int R = 5000)  // number of random probes
{
  const int p = Ipost.rows();
  double trace_est = 0.0;
  
  // Sparse Cholesky factorization of Ipost
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Ipost);
  if (chol.info() != Eigen::Success) stop("Cholesky decomposition failed");
  
  Rcpp::RNGScope scope;
  for (int r = 0; r < R; r++) {
    // ±1 Rademacher random vector
    Eigen::VectorXd z(p);
    for (int i=0; i<p; i++) z[i] = (R::runif(0,1) < 0.5 ? -1.0 : 1.0);

    // b = I_like * z
    Eigen::VectorXd b = Ilike * z;

    // Solve I_post * x = b
    Eigen::VectorXd x = chol.solve(b);

    trace_est += z.dot(x);
  }

  trace_est /= R;
  return trace_est;
}

')







Rcpp::sourceCpp(code = '
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd simulate_from_precision(const SparseMatrix<double> &Prec,
                                        const int nsim = 500) {
  const int n = Prec.rows();
  Eigen::MatrixXd draws(n, nsim);

  // Sparse Cholesky factorization (LDLT)

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Prec);
  Eigen::VectorXd inv_sqrt_eigvals = es.eigenvalues().array().rsqrt();
  Eigen::MatrixXd Q_inv_sqrt = es.eigenvectors() * inv_sqrt_eigvals.asDiagonal() * es.eigenvectors().transpose();

  // For each simulation
  for (int k = 0; k < nsim; ++k) {
    // Step 1: sample z ~ N(0, I)
    Eigen::VectorXd z = VectorXd::NullaryExpr(n, [](){ return R::rnorm(0.0, 1.0); });

    // Step 2: x = Q^{-1/2} z
    Eigen::VectorXd x = Q_inv_sqrt * z;
    draws.col(k) = x;
  }

  return draws;
}
')



Rcpp::sourceCpp(code = '
#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::SparseMatrix;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SimplicialLLT;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
MatrixXd simulate_from_precision_correct(const SparseMatrix<double> &Prec,
                                         const int nsim = 500) {
  const int n = Prec.rows();
  MatrixXd draws(n, nsim);

  // 1) Sparse Cholesky decomposition: Prec = P^T L L^T P
  SimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(Prec);
  if (solver.info() != Eigen::Success)
    stop("Cholesky decomposition failed.");

  // 2) Generate standard normals
  MatrixXd Z = MatrixXd::NullaryExpr(n, nsim, [](){ return R::rnorm(0.0, 1.0); });

  // 3) Solve L^T Y = Z (upper-triangular solve)
  MatrixXd Y = solver.matrixU().triangularView<Eigen::Upper>().solve(Z);

  // 4) Apply permutation P^T to get correct order
  draws = solver.permutationP().transpose() * Y;

  return draws;
}
')



# Function to perform matrix-vector multiply for A = Q + XWX
matvec_A <- function(v, Q, XWX) {
  return(Q %*% v + XWX %*% v)  # Adjust based on structure, use sparse matmult if possible
}

# Lanczos tridiagonalization (simplified)
lanczos_tridiag <- function(Afun, n, k, v0) {
  alpha <- numeric(k)
  beta <- numeric(k)
  V <- matrix(0, n, k)
  v <- v0 / sqrt(sum(v0^2)) # Calculate unit vector from z
  V[,1] <- v # Starting vector w0 of unit 2-norm
  
  for (j in 1:k) {
    w <- Afun(v) # Calculate w1 = A*w0
    if (j > 1) w <- w - beta[j-1] * V[,j-1] # Calculate w2 = Av_1 - beta_1*v_1
    alpha[j] <- sum(v * w) # Calculate alpha1 = w1*v1
    w <- w - alpha[j] * v # Calculate w1 - alpha1 * v1
    if (j < k) {
      beta[j] <- sqrt(sum(w^2)) # Calculate norm of w1
      if (beta[j] == 0) break
      v <- w / beta[j] # Standardize
      V[,j+1] <- as.numeric(v) # This is v2
    }
  }
  return(list(alpha=alpha, beta=beta))
}

# Approximate logdet with m random vectors and k Lanczos steps
approx_logdet <- function(Q0, XWX, m=10, k=20) {
  
  set.seed(1)
  n <- nrow(Q0)
  Afun <- function(v) matvec_A(v, Q0, XWX) 
  estimates <- numeric(m)
  for (i in 1:m) {
    z <- sample(c(-1,1), n, replace=TRUE) # Sample Rademacher vector (vectors with +-1 entries of equal prob)
    lanczos_res <- lanczos_tridiag(Afun, n, k, z) # Apply m+1 steps of Lanczos to A with z as starting vector
    alpha <- lanczos_res$alpha
    beta <- lanczos_res$beta
    
    k <- length(alpha)
    T <- diag(alpha) # Construction of Lanczos matrix (alpha diagonal and beta off diagonal, starting from beta2)
    if (k > 1) {
      off_diag <- beta[1:(k-1)] 
      T[cbind(1:(k-1), 2:k)] <- off_diag
      T[cbind(2:k, 1:(k-1))] <- off_diag
    }
    eig <- eigen(T, symmetric = TRUE)
    vals <- eig$values
    vecs <- eig$vectors
    log_est <- sum((vecs[1, ])^2 * log(pmax(vals, 1e-10)))
    estimates[i] <- log_est * sum(z^2)
  }
  return(mean(estimates))
}




# DIC calculation


calculate_DIC_tot <- function(xi_mode, Prec, Qv_mode, Xv, y, offset, v_mode, model.family = "poisson",
                              per.area = F, ID_y = NULL, ind.comparison){
  # DIC of mean
  
  if(model.family == "poisson"){
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      return(sum((y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y))))
    }
    
    log_lik_ID <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      vec <- (y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y))
      return(vec) #tapply(vec, ID_y,sum)
    }


  } else{
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))
      
      return( sum((1/exp(v_mode[1]))*((y * gammanb) - bgammanb) +
                         lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1)))
    }
    
    log_lik_ID <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))
      
      vec <-(1/exp(v_mode[1]))*((y * gammanb) - bgammanb) + 
        lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1)
      
      return(vec)
    }

    
    
  }
  
  if(is.null(ind.comparison)){
  Dev_of_mean <- -2*log_lik(xi_mode, y, Xv, offset, v_mode)
  
  # pd
  pd_est <- dim(Prec)[1]-Hutchinson_pD(Prec, Qv_mode)
  
  DIC <- Dev_of_mean + 2*pd_est
  
  }else{
    betas = simulate_from_precision_correct(Prec, nsim = 500)

    pd_est<- DIC <- vector("list",length = length(ind.comparison))
    
    loglik_samples <- vector(mode = "list")
    for (s in 1:500){
      loglik_samples[[s]] <- log_lik_ID(xi_mode + betas[,s],y,
                                        Xv, offset, v_mode)
    }
    
    for (l in 1:length(ind.comparison)){
      Dev_of_mean <- -2*sum(log_lik_ID(xi_mode+rowMeans(betas), y[ind.comparison[[l]]], Xv[ind.comparison[[l]],],
                                          offset[ind.comparison[[l]]], v_mode))
      
      result <- lapply(loglik_samples, function(vec) {
        sum(vec[ind.comparison[[l]]])
      })    
      
      mean_of_Dev <- mean(-2*do.call(cbind, result))
      
      pd_est[[l]] = mean_of_Dev - Dev_of_mean
      DIC[[l]] = 2*mean_of_Dev - Dev_of_mean
      
    }
  }
  # 
  
  if (per.area){
    if(is.null(ind.comparison)){
      DIC_area = pd_area = NULL
      betas = simulate_from_precision_correct(Prec, nsim = 500)
      loglik_samples <- vector(mode = "list")
      for (s in 1:500){
          loglik_samples[[s]] <- log_lik_ID(xi_mode + betas[,s],y,
                                         Xv, offset, v_mode)
      }
  
      Dev_of_mean <- -2*tapply(log_lik_ID(xi_mode+rowMeans(betas), y, Xv,
                            offset, v_mode), ID_y,sum)
  
      
      result = lapply(loglik_samples, rowsum, group = ID_y)
      
      mean_of_Dev <-rowMeans(-2*do.call(cbind, result))
  
      pd_area = mean_of_Dev - Dev_of_mean
      DIC_area = 2*mean_of_Dev - Dev_of_mean
    }else{
      pd_area<- DIC_area <- vector("list",length = length(ind.comparison))
      
      for (l in 1:length(ind.comparison)){
        Dev_of_mean <- -2*tapply(log_lik_ID(xi_mode+rowMeans(betas), y[ind.comparison[[l]]], 
                                     Xv[ind.comparison[[l]],],
                                   offset[ind.comparison[[l]]], v_mode), ID_y[ind.comparison[[l]]], sum)
      
        result <- lapply(loglik_samples, function(vec) {
          rowsum(vec[ind.comparison[[l]]], ID_y[ind.comparison[[l]]])
        })        
        mean_of_Dev <-rowMeans(-2*do.call(cbind, result))
      
        pd_area[[l]] = mean_of_Dev - Dev_of_mean
        DIC_area[[l]] = 2*mean_of_Dev - Dev_of_mean    
      }
    }
  }else{
    DIC_area = NULL
    pd_area = NULL
  }
  

  
  return(list(DIC = DIC, pd = pd_est, DIC_area = DIC_area, pd_area = pd_area))
}


# WAIC
calculate_WAIC_tot <- function(xi_mode, Prec, Qv_mode, Xv, y, offset, v_mode, model.family = "poisson"){
  # WAIC of mean
  
  if(model.family == "poisson"){
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      return(y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y))
    }
    
  } else{
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))
      
      return((1/exp(v_mode[1]))*((y * gammanb) - bgammanb) +
                    lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1))
    }
    
    
    
  }
  
  
  betas = simulate_from_precision_correct(Prec, nsim = 500)
  loglik_samples <- matrix(0, nrow = nrow(Xv), ncol = 500)
  for (s in 1:500){
    loglik_samples[,s] <- log_lik(xi_mode + betas[,s],y, Xv, offset, v_mode)
  }

  
  
  lppd <- sum(log(rowMeans(exp(loglik_samples))))
  pWAIC2 <- sum(apply(loglik_samples, 1, var))
  WAIC <- -2 * (lppd - pWAIC2)
  
  return(list(WAIC = WAIC, pd = pWAIC2))
}


# CPO
calculate_CPO <- function(xi_mode, Prec, Qv_mode, Xv, y, offset, v_mode, model.family = "poisson"){
  # WAIC of mean
  
  if(model.family == "poisson"){
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      return(y * log(Cvxi_mean) - Cvxi_mean-lfactorial(y))
    }
    
  } else{
    log_lik <- function(xi_mode,y, Xv, offset, v_mode){
      Cvxi_mean <-  exp(as.numeric(Xv %*% xi_mode)+log(offset))
      varval <-  Cvxi_mean + (1 / exp(v_mode[1])) * (Cvxi_mean ^ 2)
      Vnb <- Matrix::Diagonal(x = Cvxi_mean * (1/varval - (Cvxi_mean/(varval^2)) *
                                                 (1 + 2*Cvxi_mean*(1/exp(v_mode[1])))))
      Wnb <- Matrix::Diagonal(x = ((Cvxi_mean) ^ 2) * (1 / varval))
      gammanb <- exp(v_mode[1]) * log(Cvxi_mean / (Cvxi_mean+ exp(v_mode[1])))
      bgammanb <- (-1) * (exp(v_mode[1])^2) * log(exp(v_mode[1])/(exp(v_mode[1]) + Cvxi_mean))
      
      return((1/exp(v_mode[1]))*((y * gammanb) - bgammanb) +
               lgamma(y + exp(v_mode[1])) - lgamma(exp(v_mode[1])) - lgamma(y+1))
    }
    
    
    
  }
  
  
  betas = simulate_from_precision_correct(Prec, nsim = 500)
  loglik_samples <- matrix(0, nrow = nrow(Xv), ncol = 500)
  for (s in 1:500){
    loglik_samples[,s] <- log_lik(xi_mode + betas[,s],y, Xv, offset, v_mode)
  }
  
  m <- apply(-loglik_samples, 1, max) # max per observation
  harmonic_mean_inv <- rowMeans(exp(-loglik_samples - m))
  logCPO <- -(log(harmonic_mean_inv) + m)
  
  CPO_result <- exp(logCPO)
  
  return(CPO_result)
}


