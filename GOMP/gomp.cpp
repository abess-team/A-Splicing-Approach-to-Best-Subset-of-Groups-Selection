#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppEigen)]]
#include <iostream>
#include <vector>
using namespace Rcpp;
using namespace std;

std::vector<Eigen::MatrixXd> Orthonormalize(Eigen::MatrixXd& X, Eigen::VectorXd& y, Eigen::VectorXi& index, Eigen::VectorXi& gsize, int n, int p, int N, Eigen::VectorXd& weights, Eigen::VectorXd& meanx, double& meany){
  int size;
  for (int i=0;i<n;i++) {
    X.row(i) = X.row(i)*sqrt(weights(i));
  }
  meanx = X.colwise().mean();
  X = X.rowwise()-meanx.transpose();
  std::vector<Eigen::MatrixXd> mat(N);
  for (int i=0;i<N;i++) {
    size = gsize(i);
    if (size == 1) {
      Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(1, 1);
      temp(0, 0) = X.col(index(i)).norm()/sqrt(n);
      mat[i] = temp;
      X.col(index(i)) = X.col(index(i))/temp(0, 0);
    }
    else
    {
      Eigen::MatrixXd X_ind = X.block(0, index(i), n, size);
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(X_ind.transpose()*X_ind/n, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::MatrixXd Sigma = Eigen::MatrixXd::Zero(size, size);
      Sigma.diagonal() =  svd.singularValues();
      Sigma = Sigma.ldlt().solve(Eigen::MatrixXd::Identity(size, size)).sqrt();
      Eigen::MatrixXd U = svd.matrixU();
      Eigen::MatrixXd temp = U*Sigma;
      X.block(0, index(i), n, size) = X_ind*temp;
      mat[i] = temp;
    }
  }
  y = y.cwiseProduct(weights);
  meany = y.mean();
  y = y.array() - meany;
  return mat;
}

Eigen::VectorXi find_ind(std::vector<int> L, Eigen::VectorXi& index, Eigen::VectorXi& gsize, int p, int N)
{
  unsigned int J = N;
  if (L.size() == J) {
    return Eigen::VectorXi::LinSpaced(p, 0, p-1);
  }
  else 
  {
    int mark = 0;
    Eigen::VectorXi ind = Eigen::VectorXi::Zero(p);
    for (unsigned int i=0;i<L.size();i++) {
      ind.segment(mark, gsize(L[i])) = Eigen::VectorXi::LinSpaced(gsize(L[i]), index(L[i]), index(L[i])+gsize(L[i])-1);
      mark = mark + gsize(L[i]);
    }
    return ind.head(mark);
  }
}

Eigen::MatrixXd X_seg(Eigen::MatrixXd& X, int n, Eigen::VectorXi& ind) 
{
  Eigen::MatrixXd X_new(n, ind.size());
  for (int k=0;k<ind.size();k++) {
    X_new.col(k) = X.col(ind[k]);
  }
  return X_new;
}

// [[Rcpp::export]]
List gompCpp(Eigen::MatrixXd x, Eigen::VectorXd y, Eigen::VectorXd weight, int N, int n, int p,
            Eigen::VectorXi index, Eigen::VectorXi gsize, double tol, int iter_num) {
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
  Eigen::VectorXd x_mean = Eigen::VectorXd::Zero(p);
  double y_mean = 0.0;
  std::vector<Eigen::MatrixXd> mat = Orthonormalize(x, y, index, gsize, n, p, N, weight, x_mean, y_mean);
  
  List A_out(iter_num);
  Eigen::MatrixXd beta_out(p, iter_num);
  vector<int> A;
  int ind;
  double maxcor;
  Eigen::VectorXd residual = y - x*beta;
  Eigen::VectorXd cor = Eigen::VectorXd::Zero(N);
  int iter = 0;
  for (iter;iter<iter_num;iter++) {
    for (int j=0;j < N;j++) {
      Eigen::MatrixXd X_temp = x.block(0, index(j), n, gsize(j));
      cor[j] = sqrt((X_temp.transpose()*residual).squaredNorm());
    }
    maxcor = cor.maxCoeff(&ind);
    // Rcout<<maxcor<<"\n";
    if (maxcor <= tol) break;
    A.push_back(ind);
    sort(A.begin(), A.end());
    A_out(iter) = A;
    cor = Eigen::VectorXd::Zero(N);
    Eigen::VectorXi ind_temp = find_ind(A, index, gsize, p, N);
    Eigen::MatrixXd X_temp = X_seg(x, n, ind_temp);
    Eigen::VectorXd beta_temp = X_temp.householderQr().solve(y);
    beta = Eigen::VectorXd::Zero(p);
    for (int m=0;m<ind_temp.size();m++) {
      beta(ind_temp(m)) = beta_temp(m);  
    }
    beta_out.col(iter) = beta;
    residual = y - X_temp*beta_temp;
  }
  Eigen::VectorXd coef_out = Eigen::VectorXd::Zero(iter_num);
  for (int i=0;i<N;i++) {
    beta_out.block(index(i), 0, gsize(i), iter_num) = mat[i] * beta_out.block(index(i), 0, gsize(i), iter_num);
  }
  coef_out = y_mean*Eigen::VectorXd::Ones(iter_num) - beta_out.transpose()*x_mean;
  return List::create(Named("beta")=beta_out, Named("coef")=coef_out,  Named("A_out")=A_out, Named("iter") = iter);
}


