
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
using namespace Numer;


// [[Rcpp::export]]
double C_erf(double x){
  double result = 2 * Rf_pnorm5(x * sqrt(2), 0.0, 1.0, 1, 0) - 1;
  return result;
}


// [[Rcpp::export]]
double Intermediate_LT_Halfnorm_det(double x, double y, double b, double m, double sigma){
  
  double m2p1 = pow(m,2) + 1;
  double bmxy = b + m * x - y;
  
  double val = ( -pow(2 * 3.1415926, 0.5) * pow(m2p1, 0.5) * bmxy * C_erf( bmxy / (pow(2, 0.5) * pow(m2p1, 0.5) * sigma) ) ) -
    ( 2 * m2p1 * sigma * exp( -pow(bmxy, 2) / (2 * m2p1 * pow(sigma, 2)) ) );
  
  return val;
}


// [[Rcpp::export]]
double LineTransect_HalfNorm_det(Rcpp::NumericVector xBounds, Rcpp::NumericVector yBounds,
                                 double b, double m, 
                                 double sigmaDet, double p0){
  
  double lineTransectIntgrl;
  
  if(m == 0){   
    
    lineTransectIntgrl = p0 * -pow((3.1415926/2), 0.5) * sigmaDet * (xBounds(1) - xBounds(0)) * 
      (C_erf((b - yBounds(1))/(pow(2,0.5) * sigmaDet)) - C_erf((b - yBounds(0))/(pow(2,0.5) * sigmaDet)));
    
  } else{     // might also want a case to handle undefined/infinite m / vertical lines
    
    //double sqrt_m2p1 = pow((pow(m,2) + 1), 0.5);
    
    lineTransectIntgrl = p0 * sigmaDet / (2 * m) * ( ( Intermediate_LT_Halfnorm_det(xBounds(1), yBounds(1), b, m, sigmaDet) -
      Intermediate_LT_Halfnorm_det(xBounds(1), yBounds(0), b, m, sigmaDet) ) - 
      ( Intermediate_LT_Halfnorm_det(xBounds(0), yBounds(1), b, m, sigmaDet) -
      Intermediate_LT_Halfnorm_det(xBounds(0), yBounds(0), b, m, sigmaDet) ) );
    
    // HAD ISSUES WITH THIS FORMULATION
    //lineTransectIntgrl = (p0 * pow(3.1415926/2, 0.5) * sqrt_m2p1 * sigmaDet / m) *
    
    //((((pow(2/3.1415926, 0.5) * sqrt_m2p1 * sigmaDet * -exp(-pow((b + m * xBounds(1) - yBounds(1)), 2) / (2 * sqrt_m2p1 * pow(sigmaDet,2)))) - 
    //((b + m * xBounds(1) - yBounds(1)) * C_erf((b + m * xBounds(1) - yBounds(1)) / (pow(2,0.5) * sqrt_m2p1 * sigmaDet)))) - 
    
    //((pow(2/3.1415926, 0.5) * sqrt_m2p1 * sigmaDet * -exp(-pow((b + m * xBounds(1) - yBounds(0)), 2) / (2 * sqrt_m2p1 * pow(sigmaDet,2)))) - 
    //((b + m * xBounds(1) - yBounds(0)) * C_erf((b + m * xBounds(1) - yBounds(0)) / (pow(2,0.5) * sqrt_m2p1 * sigmaDet))))) - 
    
    //(((pow(2/3.1415926, 0.5) * sqrt_m2p1 * sigmaDet * -exp(-pow((b + m * xBounds(0) - yBounds(1)), 2) / (2 * sqrt_m2p1 * pow(sigmaDet,2)))) - 
    //((b + m * xBounds(0) - yBounds(1)) * C_erf((b + m * xBounds(0) - yBounds(1)) / (pow(2,0.5) * sqrt_m2p1 * sigmaDet)))) - 
    
    //((pow(2/3.1415926, 0.5) * sqrt_m2p1 * sigmaDet * -exp(-pow((b + m * xBounds(0) - yBounds(0)), 2) / (2 * sqrt_m2p1 * pow(sigmaDet,2)))) - 
    //((b + m * xBounds(0) - yBounds(0)) * C_erf((b + m * xBounds(0) - yBounds(0)) / (pow(2,0.5) * sqrt_m2p1 * sigmaDet))))));
    
    //Rcpp::Rcout << "New 2 formulation = " << t2 << "\n";
    
  }
  
  return lineTransectIntgrl;
}


// [[Rcpp::export]]
double EndPoint_HalfNorm_det(Rcpp::NumericVector xBounds, Rcpp::NumericVector yBounds,
                             double epX, double epY, 
                             double sigmaDet, double p0){
  double endPointIntgrl = p0 * 0.5 * 3.1415926 * pow(sigmaDet, 2) * 
    (C_erf((epX - xBounds(1)) / (pow(2,0.5) * sigmaDet)) - C_erf((epX - xBounds(0)) / (pow(2,0.5) * sigmaDet))) * 
    (C_erf((epY - yBounds(1)) / (pow(2,0.5) * sigmaDet)) - C_erf((epY - yBounds(0)) / (pow(2,0.5) * sigmaDet)));
  
  return endPointIntgrl;
}


// functor for integration of hazard detection function
class HazardDet_TL: public MFunc{
private:
  const double sigmaDet;
  const double p0;
  const double beta;
  const double b;
  const double m;
  
  
  
public:
  HazardDet_TL(const double& sigmaDet_, const double& p0_, const double& beta_, 
               const double& b_, const double& m_) : 
  sigmaDet(sigmaDet_), p0(p0_), beta(beta_), b(b_), m(m_) {}
  
  double operator()(Constvec& x){
    double val1 = std::abs(b + m * x[0] - x[1]) / std::pow(m * m + 1, 0.5);
    double val2 = -1 * std::pow(val1/sigmaDet, -beta);
    double val3 = 1 - std::exp(val2);
    return p0 * val3; 
  }
};

// [[Rcpp::export]]
double HazardDetIntgrl_TL(double sigmaDet, double p0,
                          double beta, double b, double m,
                          Rcpp::NumericVector Lwr, Rcpp::NumericVector Upr){
  HazardDet_TL f(sigmaDet, p0, beta, b, m);
  Eigen::VectorXd lower(2);
  lower << Lwr(0), Lwr(1);
  Eigen::VectorXd upper(2);
  upper << Upr(0), Upr(1);
  double err_est;
  int err_code;
  const int maxeval = 5000;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval);
  return res;
}


// functor for integration of hazard detection function
class HazardDet_EP: public MFunc{
private:
  const double sigmaDet;
  const double p0;
  const double beta;
  const double xEP;
  const double yEP;
  
  
  
public:
  HazardDet_EP(const double& sigmaDet_, const double& p0_, const double& beta_, 
               const double& xEP_, const double& yEP_) : 
  sigmaDet(sigmaDet_), p0(p0_), beta(beta_), xEP(xEP_), yEP(yEP_) {}
  
  double operator()(Constvec& x){
    double val1 = std::pow( std::pow(x[0] - xEP, 2) + std::pow(x[1] - yEP, 2), 0.5);
    double val2 = -1 * std::pow(val1/sigmaDet, -beta);
    double val3 = 1 - std::exp(val2);
    return p0 * val3;   
  }
};

// [[Rcpp::export]]
double HazardDetIntgrl_EP(double sigmaDet, double p0,
                          double beta, double xEP, double yEP,
                          Rcpp::NumericVector Lwr, Rcpp::NumericVector Upr){
  HazardDet_EP f(sigmaDet, p0, beta, xEP, yEP);
  Eigen::VectorXd lower(2);
  lower << Lwr(0), Lwr(1);
  Eigen::VectorXd upper(2);
  upper << Upr(0), Upr(1);
  double err_est;
  int err_code;
  const int maxeval = 5000;
  const double res = integrate(f, lower, upper, err_est, err_code, maxeval);
  return res;
}



// [[Rcpp::export]]
arma::cube C_pMove_Normal_Full(int K, int nCells, Rcpp::NumericMatrix habitat, Rcpp::NumericMatrix distMat,
                                Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars){
  
  arma::cube pMove(nCells, nCells, K);
  arma::mat pMove_s_sums(nCells, K);
  arma::cube pMoveStd(nCells, nCells, K);
  double covarSum;
  //double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  
  if(nCovars == 0){
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //activity center cell
          pMove(i,j,k) = exp( -0.5 * pow(distMat(i,j) / betaMoves(0), 2) );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  } else{
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
        }
        for(int j = 0; j < nCells; j++){  //activity center cell
          pMove(i,j,k) = exp( -0.5 * pow(distMat(i,j) / betaMoves(0), 2) + covarSum );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  }
  
  for(int k = 0; k < K; k++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,k) = pMove(i,j,k) / pMove_s_sums(j,k);
      }
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
arma::cube C_pMove_SLTA_Full(int K, int nCells, Rcpp::NumericMatrix distMat, Rcpp::NumericMatrix angleMat,
                             double alpha, double beta, double mu, double gamma){
  
  arma::cube pMove(nCells, nCells, K);
  arma::mat pMove_s_sums(nCells, K);
  arma::cube pMoveStd(nCells, nCells, K);
  //double covarSum;
  pMove_s_sums.zeros();
  
  //if(alpha >= 1){
  //  for(int k = 0; k < K; k++){
  //    for(int i = 0; i < nCells; i++){    //movement outcome cell
  //      for(int j = 0; j < nCells; j++){  //activity center cell
  //        //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
  //        //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
  //        
  //        pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j)) *
  //          sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
  //        pMove_s_sums(j,k) += pMove(i,j,k);
  //      }
  //    }
  //  }
  //} else{
  //  for(int k = 0; k < K; k++){
  //    for(int i = 0; i < nCells; i++){    //movement outcome cell
  //      for(int j = 0; j < nCells; j++){  //activity center cell
  //        if(distMat(i,j) == 0){
  //          pMove(i,j,k) = 0;
  //        } else{
  //          //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
  //          //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
  //          pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j)) *
  //            sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
  //        }
  //        pMove_s_sums(j,k) += pMove(i,j,k);
  //      }
  //    }
  //  }
  //}
  
  if(gamma < 700){
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //activity center cell
          if(distMat(i,j) == 0){
            //pMove(i,j,k) = 0;
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(0.01) - beta * 0.01) *
              sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
          } else{
            //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
            //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j)) *
              sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
          }
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  } else{ // uniform angle -> wrapped cauchy approaches uniform at much smaller values of gamma too
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //activity center cell
          if(distMat(i,j) == 0){
            //pMove(i,j,k) = 0;
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(0.01) - beta * 0.01);
          } else{
            //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
            //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j));
          }
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  }

  for(int k = 0; k < K; k++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,k) = pMove(i,j,k) / pMove_s_sums(j,k);
      }
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
Rcpp::List C_pMove_SL_Full(int K, int nCells, Rcpp::NumericMatrix habitat, 
                           Rcpp::NumericMatrix distMat, Rcpp::NumericMatrix angleMat,
                           Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars,
                           double alpha, double beta, double pi){
  
  arma::cube pMove(nCells, nCells, K);
  arma::mat pMove_s_sums(nCells, K);
  arma::mat pMove_North_sums(nCells, K);
  arma::mat pMove_West_sums(nCells, K);
  arma::mat pMove_South_sums(nCells, K);
  arma::mat pMove_East_sums(nCells, K);
  arma::cube pMoveStd(nCells, nCells, K);
  arma::cube pMoveStd_North(nCells, nCells, K);
  arma::cube pMoveStd_West(nCells, nCells, K);
  arma::cube pMoveStd_South(nCells, nCells, K);
  arma::cube pMoveStd_East(nCells, nCells, K);
  double covarSum;
  pMove_s_sums.zeros();
  
  if(nCovars == 0){ // uniform angle -> wrapped cauchy approaches uniform at much smaller values of gamma too
    
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //activity center cell
          if(distMat(i,j) == 0){
            //pMove(i,j,k) = 0;
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(0.01) - beta * 0.01);
          } else{
            //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
            //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j));
          }
          pMove_s_sums(j,k) += pMove(i,j,k);
          
          if(angleMat(j,i) > pi/4 && angleMat(j,i) <= 3*pi/4){
            pMove_North_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) > 3*pi/4 && angleMat(j,i) <= 5*pi/4){
            pMove_West_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) > 5*pi/4 && angleMat(j,i) <= 7*pi/4){
            pMove_South_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) == 0){
            pMove_North_sums(j,k) += pMove(i,j,k)/4;
            pMove_West_sums(j,k) += pMove(i,j,k)/4;
            pMove_South_sums(j,k) += pMove(i,j,k)/4;
            pMove_East_sums(j,k) += pMove(i,j,k)/4;
          } else{
            pMove_East_sums(j,k) += pMove(i,j,k);
          }
          
        }
      }
    }
    
  } else{
    
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
        }
        
        for(int j = 0; j < nCells; j++){  //activity center cell
          if(distMat(i,j) == 0){
            //pMove(i,j,k) = 0;
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(0.01) - beta * 0.01 + covarSum);
          } else{
            //pMove(i,j,k) = pow(beta, alpha) * pow(distMat(i,j), alpha-1) * exp(-beta * distMat(i,j)) *
            //  sinh(gamma) / (cosh(gamma) - cos(angleMat(i,j) - mu));
            pMove(i,j,k) = exp( alpha * log(beta) + (alpha-1) * log(distMat(i,j)) - beta * distMat(i,j) + covarSum);
          }
          pMove_s_sums(j,k) += pMove(i,j,k);
          
          if(angleMat(j,i) > pi/4 && angleMat(j,i) <= 3*pi/4){
            pMove_North_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) > 3*pi/4 && angleMat(j,i) <= 5*pi/4){
            pMove_West_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) > 5*pi/4 && angleMat(j,i) <= 7*pi/4){
            pMove_South_sums(j,k) += pMove(i,j,k);
          } else if(angleMat(j,i) == 0){
            pMove_North_sums(j,k) += pMove(i,j,k)/4;
            pMove_West_sums(j,k) += pMove(i,j,k)/4;
            pMove_South_sums(j,k) += pMove(i,j,k)/4;
            pMove_East_sums(j,k) += pMove(i,j,k)/4;
          } else{
            pMove_East_sums(j,k) += pMove(i,j,k);
          }
          
        }
      }
    }
    
  }
  
  for(int k = 0; k < K; k++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,k) = pMove(i,j,k) / pMove_s_sums(j,k);
        
        if(angleMat(j,i) > pi/4 && angleMat(j,i) <= 3*pi/4){
          pMoveStd_North(i,j,k) = pMove(i,j,k) / pMove_North_sums(j,k);
        } else if(angleMat(j,i) > 3*pi/4 && angleMat(j,i) <= 5*pi/4){
          pMoveStd_West(i,j,k) = pMove(i,j,k) / pMove_West_sums(j,k);
        } else if(angleMat(j,i) > 5*pi/4 && angleMat(j,i) <= 7*pi/4){
          pMoveStd_South(i,j,k) = pMove(i,j,k) / pMove_South_sums(j,k);
        } else if(angleMat(j,i) == 0){
          pMoveStd_North(i,j,k) += pMove(i,j,k)/(4 * pMove_North_sums(j,k));
          pMoveStd_West(i,j,k) += pMove(i,j,k)/(4 * pMove_West_sums(j,k));
          pMoveStd_South(i,j,k) += pMove(i,j,k)/(4 * pMove_South_sums(j,k));
          pMoveStd_East(i,j,k) += pMove(i,j,k)/(4 * pMove_East_sums(j,k));
        } else{
          pMoveStd_East(i,j,k) = pMove(i,j,k) / pMove_East_sums(j,k);
        }
        
      }
    }
  }
  
  
  Rcpp::List output = Rcpp::List::create(pMoveStd, pMoveStd_North, pMoveStd_West, pMoveStd_South, pMoveStd_East);
  return output;
}




// [[Rcpp::export]]
arma::cube C_pMove_Laplace_Full(int K, int nCells, Rcpp::NumericMatrix habitat, Rcpp::NumericMatrix distMat,
                                Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars){
  
  arma::cube pMove(nCells, nCells, K);
  arma::mat pMove_s_sums(nCells, K);
  arma::cube pMoveStd(nCells, nCells, K);
  double covarSum;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  
  if(nCovars == 0){
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //activity center cell
          pMove(i,j,k) = exp( alphaMove * distMat(i,j) );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  } else{
    for(int k = 0; k < K; k++){
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
        }
        for(int j = 0; j < nCells; j++){  //activity center cell
          pMove(i,j,k) = exp( alphaMove * distMat(i,j) + covarSum );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
    }
  }
  
  for(int k = 0; k < K; k++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,k) = pMove(i,j,k) / pMove_s_sums(j,k);
      }
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
Rcpp::List C_pMove_Laplace_Open_move(int K, int nCells, int nCells_aug, arma::mat habitat, arma::mat distMat,
                                     arma::mat habitat_aug, arma::mat distMat_aug,
                                     arma::vec betaMoves, arma::mat covarCols, int nCovars){
  
  arma::cube pMove(nCells_aug, nCells, K);
  arma::mat pMove_s_sums(nCells, K);
  arma::mat pMove_s_sums_out(nCells, K);
  arma::mat pMove_s_sums_in(nCells, K);
  arma::cube pMoveStd(nCells, nCells, K);
  arma::mat exitMat(nCells, K);
  arma::mat enterDist(nCells, K);
  double covarSum;
  double covarSum_in;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  

  if(nCovars == 0){
    for(int k = 0; k < K; k++){
      
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //origin cell
          pMove(i,j,k) = exp( alphaMove * distMat(i,j) );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
      
      for(int i = 0; i < nCells_aug; i++){  //movement outcome cell
        for(int j = 0; j < nCells; j++){    //origin cell
          pMove_s_sums_out(j,k) += exp(alphaMove * distMat_aug(i,j) );
        }
      }
      
    }
    
    pMove_s_sums_in = pMove_s_sums_out;
    
  } else{
    for(int k = 0; k < K; k++){
      
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
        }
        for(int j = 0; j < nCells; j++){  //origin cell
          pMove(i,j,k) = exp( alphaMove * distMat(i,j) + covarSum );
          pMove_s_sums(j,k) += pMove(i,j,k);
        }
      }
      
      for(int i = 0; i < nCells_aug; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat_aug(i,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
        }
        for(int j = 0; j < nCells; j++){  //origin cell
          pMove(i,j,k) = exp( alphaMove * distMat_aug(i,j) + covarSum );
          pMove_s_sums_out(j,k) += pMove(i,j,k);
          
          covarSum_in = 0;
          for(int nc = 0; nc < nCovars; nc++){
            covarSum_in += betaMoves(nc+1) * habitat(j,covarCols(k,nc)-1); //passing the column index for R; subtract 1 for C++ column index
          }
          pMove_s_sums_in(j,k) += exp( alphaMove * distMat_aug(j,i) + covarSum_in );
        }
      }
      
    }
  }
  
  for(int k = 0; k < K; k++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,k) = pMove(i,j,k) / (pMove_s_sums(j,k) + pMove_s_sums_out(j,k));
      }
    }
  }
  
  for(int k = 0; k < K; k++){
    for(int j = 0; j < nCells; j++){
      exitMat(j,k) = pMove_s_sums_out(j,k) / (pMove_s_sums(j,k) + pMove_s_sums_out(j,k));
      enterDist(j,k) = pMove_s_sums_in(j,k) / accu(pMove_s_sums_in.col(k));
    }
  }
  
  Rcpp::List output = Rcpp::List::create(pMoveStd, exitMat, enterDist);
  return output;
}


// [[Rcpp::export]]
arma::cube C_pMove_RWAC_Laplace_Full(int nCells, int nACCells, Rcpp::NumericMatrix habitat,
                                Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars,
                                double rho){
  
  arma::cube pMove(nCells, nCells, nACCells);
  arma::mat pMove_s_sums(nCells, nACCells);
  arma::cube pMoveStd(nCells, nCells, nACCells);
  double covarSum;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();

  if(nCovars == 0){
    for(int ac = 0; ac < nACCells; ac++){   // activity center cell
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        for(int j = 0; j < nCells; j++){  //movement origin cell
          // Could switch the looping variables to j=outcome, i=origin and save the average location, x=rho * habitat(ac,0) + (1-rho) * habitat(i,0),
          // for each loop through outcome cells
          pMove(i,j,ac) = exp( alphaMove * pow( pow(rho * habitat(ac,0) + (1-rho) * habitat(j,0) - habitat(i,0), 2) + 
            pow(rho * habitat(ac,1) + (1-rho) * habitat(j,1) - habitat(i,1), 2), 0.5));
          //pMove(i,j,k) = exp( alphaMove * distMat(i,j) );
          pMove_s_sums(j,ac) += pMove(i,j,ac);
        }
      }
    }
  } else{
    for(int ac = 0; ac < nCells; ac++){   // activity center cell
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(0,nc)-1); //ASSUMES COVARS ARE CONSTANT - ALWAYS = FIRST OCCASION VALUES
        }
        for(int j = 0; j < nCells; j++){  //movement origin cell
          pMove(i,j,ac) = exp( alphaMove * pow( pow(rho * habitat(ac,0) + (1-rho) * habitat(j,0) - habitat(i,0), 2) + 
            pow(rho * habitat(ac,1) + (1-rho) * habitat(j,1) - habitat(i,1), 2), 0.5) + 
            covarSum );
          pMove_s_sums(j,ac) += pMove(i,j,ac);
        }
      }
    }
  }
  
  for(int ac = 0; ac < nCells; ac++){
    for(int i = 0; i < nCells; i++){
      for(int j = 0; j < nCells; j++){
        pMoveStd(i,j,ac) = pMove(i,j,ac) / pMove_s_sums(j,ac);
      }
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
arma::cube C_pMove_CorrelatedRW_Laplace_Full(int nCells, Rcpp::NumericMatrix habitat,
                                     Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars,
                                     arma::mat rotation){
  
  arma::cube pMove(nCells, nCells, nCells);
  arma::mat pMove_s_sums(nCells, nCells);
  arma::cube pMoveStd(nCells, nCells, nCells);
  double covarSum;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  arma::vec expected_location(2);
  
  if(nCovars == 0){
    for(int k_two = 0; k_two < nCells; k_two++){   // k-2 location
      for(int k_one = 0; k_one < nCells; k_one++){    // k-1 location 
        expected_location(0) = rotation(0,0) * (habitat(k_one,0) - habitat(k_two,0)) +
          rotation(0,1) * (habitat(k_one,1) - habitat(k_two,1)) +
          habitat(k_one,0);
        expected_location(1) = rotation(1,0) * (habitat(k_one,0) - habitat(k_two,0)) +
          rotation(1,1) * (habitat(k_one,1) - habitat(k_two,1)) +
          habitat(k_one,1);
        for(int i = 0; i < nCells; i++){  //movement outcome cell
          // Could switch the looping variables to j=outcome, i=origin and save the average location, x=rho * habitat(ac,0) + (1-rho) * habitat(i,0),
          // for each loop through outcome cells
          pMove(i,k_two,k_one) = exp( alphaMove * pow( pow(expected_location(0) - habitat(i,0), 2) + 
            pow(expected_location(1) - habitat(i,1), 2), 0.5));
          //pMove(i,j,k) = exp( alphaMove * distMat(i,j) );
          pMove_s_sums(k_two,k_one) += pMove(i,k_two,k_one);
        }
      }
    }
  } else{
    for(int k_two = 0; k_two < nCells; k_two++){   // activity center cell
      for(int k_one = 0; k_one < nCells; k_one++){    //movement outcome cell
        expected_location(0) = rotation(0,0) * (habitat(k_one,0) - habitat(k_two,0)) +
          rotation(0,1) * (habitat(k_one,1) - habitat(k_two,1)) +
          habitat(k_one,0);
        expected_location(1) = rotation(1,0) * (habitat(k_one,0) - habitat(k_two,0)) +
          rotation(1,1) * (habitat(k_one,1) - habitat(k_two,1)) +
          habitat(k_one,1);
        for(int i = 0; i < nCells; i++){  //movement origin cell
          covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
          for(int nc = 0; nc < nCovars; nc++){
            covarSum += betaMoves(nc+1) * habitat(i,covarCols(0,nc)-1); //ASSUMES COVARS ARE CONSTANT - ALWAYS = FIRST OCCASION VALUES
          }
          pMove(i,k_two,k_one) = exp( alphaMove * pow( pow(expected_location(0) - habitat(i,0), 2) + 
            pow(expected_location(1) - habitat(i,1), 2), 0.5) + 
            covarSum );
          pMove_s_sums(k_two,k_one) += pMove(i,k_two,k_one);
        }
      }
    }
  }
  
  for(int k_two = 0; k_two < nCells; k_two++){
    for(int k_one = 0; k_one < nCells; k_one++){
      for(int i = 0; i < nCells; i++){
        pMoveStd(i,k_two,k_one) = pMove(i,k_two,k_one) / pMove_s_sums(k_two,k_one);
      }
    }
  }
  
  return pMoveStd;
}



// [[Rcpp::export]]
arma::mat C_pMove_RWAC_Laplace_Initial(int nCells, int nACCells, Rcpp::NumericMatrix habitat, Rcpp::NumericMatrix distMat_AC,
                                     Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars,
                                     Rcpp::NumericVector p_init_state_ac){
  
  arma::mat pMove(nCells, nACCells);
  arma::vec pMove_s_sums(nACCells);
  arma::mat pMoveStd(nCells, nACCells);
  double covarSum;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  
  if(nCovars == 0){
    for(int ac = 0; ac < nACCells; ac++){   // activity center cell
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        pMove(i,ac) = exp( alphaMove * distMat_AC(i,ac));
        pMove_s_sums(ac) += pMove(i,ac);
      }
    }
  } else{
    for(int ac = 0; ac < nCells; ac++){   // activity center cell
      for(int i = 0; i < nCells; i++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(i,covarCols(0,nc)-1); //ASSUMES COVARS ARE CONSTANT - ALWAYS = FIRST OCCASION VALUES
        }
        pMove(i,ac) = exp(alphaMove * distMat_AC(i,ac) + covarSum);
        pMove_s_sums(ac) += pMove(i,ac);
      }
    }
  }
  
  for(int ac = 0; ac < nCells; ac++){
    for(int i = 0; i < nCells; i++){
        pMoveStd(i,ac) = p_init_state_ac(ac) * pMove(i,ac) / pMove_s_sums(ac);
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
arma::mat C_pMove_CorrelatedRW_Laplace_Initial(int nCells, Rcpp::NumericMatrix habitat, Rcpp::NumericMatrix distMat,
                                       Rcpp::NumericVector betaMoves, Rcpp::NumericMatrix covarCols, int nCovars,
                                       Rcpp::NumericVector p_init_state){
  
  arma::mat pMove(nCells, nCells);
  arma::vec pMove_s_sums(nCells);
  arma::mat pMoveStd(nCells, nCells);
  double covarSum;
  double alphaMove = -1 / betaMoves(0);
  pMove_s_sums.zeros();
  
  if(nCovars == 0){
    for(int k_two = 0; k_two < nCells; k_two++){   // k-2 location
      for(int k_one = 0; k_one < nCells; k_one++){    // k-1 location
        pMove(k_one,k_two) = exp( alphaMove * distMat(k_one,k_two));
        pMove_s_sums(k_two) += pMove(k_one,k_two);
      }
    }
  } else{
    for(int k_two = 0; k_two < nCells; k_two++){   // activity center cell
      for(int k_one = 0; k_one < nCells; k_one++){    //movement outcome cell
        covarSum = 0;                       // covariates * betaMoves for movement outcome cell i
        for(int nc = 0; nc < nCovars; nc++){
          covarSum += betaMoves(nc+1) * habitat(k_one,covarCols(0,nc)-1); //ASSUMES COVARS ARE CONSTANT - ALWAYS = FIRST OCCASION VALUES
        }
        pMove(k_one,k_two) = exp(alphaMove * distMat(k_one,k_two) + covarSum);
        pMove_s_sums(k_two) += pMove(k_one,k_two);
      }
    }
  }
  
  for(int k_two = 0; k_two < nCells; k_two++){
    for(int k_one = 0; k_one < nCells; k_one++){
      pMoveStd(k_one,k_two) = p_init_state(k_two) * pMove(k_one,k_two) / pMove_s_sums(k_two);
    }
  }
  
  return pMoveStd;
}


// [[Rcpp::export]]
arma::mat C_pDet_Full(int K_tot, int nCells, int nDetCells, arma::cube lines_arr,
                      arma::cube detGrid,
                      Rcpp::NumericVector sigmaDets, Rcpp::NumericVector g0s, 
                      Rcpp::NumericVector hSigmas, Rcpp::NumericVector hBetas,
                      double detArea, double truncDist,
                      Rcpp::NumericVector SDcols, int nSDcols,
                      Rcpp::NumericVector g0cols, int nG0cols,
                      Rcpp::NumericVector hscols, int nhscols,
                      Rcpp::NumericVector hbcols, int nhbcols){
  
  arma::mat pDet_arr(nCells, K_tot);
  Rcpp::NumericVector xBounds(2);
  Rcpp::NumericVector yBounds(2);
  Rcpp::NumericVector Lwr(2);
  Rcpp::NumericVector Upr(2);
  Rcpp::NumericVector pDet(nDetCells);
  Rcpp::NumericVector pDetSum(nCells);
  Rcpp::NumericVector pDetCount(nCells);
  Rcpp::NumericVector pDetAvg(nCells);
  double sigmaDet = exp(sigmaDets(0));
  double utSigmaDet;
  double g0 = 1 / (1 + exp(-g0s(0)));
  double utg0;
  double hSigma = exp(hSigmas(0));
  double utHS;
  double hBeta = exp(hBetas(0));
  double utHB;
  
  if(nSDcols == 0){                 
    sigmaDet = exp(sigmaDets(0));
  } 
  
  if(nG0cols == 0){                 
    g0 = 1 / (1 + exp(-g0s(0)));
  }
  
  if(nhscols == 0){
    hSigma = exp(hSigmas(0));
  }
  
  if(nhbcols == 0){
    hBeta = exp(hBetas(0));
  }
  
  for(int KT = 0; KT < K_tot; KT++){
    
    for(int i = 0; i < nCells; i++){
      pDetSum(i) = 0;
      pDetCount(i) = 0;
      pDetAvg(i) = 0;
    }
    
    for(int j = 0; j < nDetCells; j++){
      
      if(detGrid(j, 8, KT) >= truncDist){    //if grid cell is farther away than the truncation distance, detection probability is 0
        
        pDet(j) = 0;
        
      } else{
        
        // accomodate variable numbers of covariates in g0
        if(nG0cols > 0){
          utg0 = g0s(0);
          for(int k = 1; k < (nG0cols + 1); k++){
            utg0 += g0s(k) * detGrid(j, g0cols(k-1)-1, KT);
          }
          g0 = 1 / (1 + exp(-utg0));
        }
        
        if(detGrid(j, 10, KT) == 0){   //half-normal detection function
          
          xBounds(0) = detGrid(j, 3, KT);
          xBounds(1) = detGrid(j, 4, KT);
          yBounds(0) = detGrid(j, 5, KT);
          yBounds(1) = detGrid(j, 6, KT);
          
          // accomodate variable numbers of covariates in sigmaDet
          if(nSDcols > 0){
            utSigmaDet = sigmaDets(0);
            for(int k = 1; k < (nSDcols + 1); k++){
              utSigmaDet += sigmaDets(k) * detGrid(j, SDcols(k-1)-1, KT); 
            }
            sigmaDet = exp(utSigmaDet);
          }
          
          if(sigmaDet > 100000000){   //Deal with numerical overflow; Probably need to find seet of conditions to handle this for hazard model too
            
            pDet(j) = g0;
            
          } else{
            
            if(detGrid(j, 9, KT) == 1){  //line transect integral
              
              pDet(j) = LineTransect_HalfNorm_det(xBounds, yBounds, lines_arr(detGrid(j, 7, KT)-1, 5, KT), lines_arr(detGrid(j, 7, KT)-1, 4, KT), sigmaDet, g0) / detArea;
              
            } else if(detGrid(j, 9, KT) == 2){ //end point (1-start) integral 
              
              pDet(j) = EndPoint_HalfNorm_det(xBounds, yBounds, lines_arr(detGrid(j, 7, KT)-1, 0, KT), lines_arr(detGrid(j, 7, KT)-1, 2, KT), sigmaDet, g0) / detArea;
              
            } else if(detGrid(j,9, KT) == 3){ //end point (2-end) integral
              
              pDet(j) = EndPoint_HalfNorm_det(xBounds, yBounds, lines_arr(detGrid(j, 7, KT)-1, 1, KT), lines_arr(detGrid(j, 7, KT)-1, 3, KT), sigmaDet, g0) / detArea;
              
            }
            
          }
          
          
          
        } else if(detGrid(j, 10, KT) == 1){  //hazard detection function
          
          Lwr(0) = detGrid(j, 3, KT);
          Lwr(1) = detGrid(j, 5, KT);
          Upr(0) = detGrid(j, 4, KT);
          Upr(1) = detGrid(j, 6, KT);
          
          // accomodate variable numbers of covariates in hazard sigma
          if(nhscols > 0){
            utHS = hSigmas(0);
            for(int k = 1; k < (nhscols + 1); k++){
              utHS += hSigmas(k) * detGrid(j, hscols(k-1)-1, KT); 
            }
            hSigma = exp(utHS);
          }
          
          // accomodate variable numbers of covariates in hazard beta
          if(nhbcols > 0){
            utHB = hBetas(0);
            for(int k = 1; k < (nhbcols + 1); k++){
              utHB += hBetas(k) * detGrid(j, hbcols(k-1)-1, KT); 
            }
            hBeta = exp(utHB);
          }
          
          
          if(detGrid(j, 9, KT) == 1){  //line transect integral
            
            pDet(j) = HazardDetIntgrl_TL(hSigma, g0, hBeta, lines_arr(detGrid(j, 7, KT)-1, 5, KT), lines_arr(detGrid(j, 7, KT)-1, 4, KT), Lwr, Upr) / detArea;
            
          } else if(detGrid(j, 9, KT) == 2){   //end point (1) integral
            
            pDet(j) = HazardDetIntgrl_EP(hSigma, g0, hBeta, lines_arr(detGrid(j, 7, KT)-1, 0, KT), lines_arr(detGrid(j, 7, KT)-1, 2, KT), Lwr, Upr) / detArea;
            
          } else if(detGrid(j, 9, KT) == 3){   //end point (2) integral
            
            pDet(j) = HazardDetIntgrl_EP(hSigma, g0, hBeta, lines_arr(detGrid(j, 7, KT)-1, 1, KT), lines_arr(detGrid(j, 7, KT)-1, 3, KT), Lwr, Upr) / detArea;
            
          }
          
        }
        
        pDetSum(detGrid(j, 0, KT)-1) += pDet(j);   //only need to add if within truncDist, otherwise would just be adding 0
        
      }
      
      pDetCount(detGrid(j, 0, KT)-1) += 1;
      
    }
    
    
    for(int i = 0; i < nCells; i++){
      pDet_arr(i, KT) = pDetSum(i)/pDetCount(i);
    }
    
  }
  
  return pDet_arr;
  
}


// [[Rcpp::export]]
Rcpp::List C_pNoDet(int nCells, int K, Rcpp::NumericMatrix pDet, Rcpp::NumericMatrix pMove){
  
  Rcpp::NumericMatrix pNoDetOcc(nCells,K);
  Rcpp::NumericVector pNoDetGivenAC(nCells);
  
  for(int i = 0; i < nCells; i++){        //activity center
    pNoDetGivenAC(i) = 1;
    for(int k = 0; k < K; k++){           //occasion
      for(int j = 0; j < nCells; j++){    //movement outcome cell
        pNoDetOcc(i,k) += pMove(j,i) * (1 - pDet(j,k));
      }
      pNoDetGivenAC(i) *= pNoDetOcc(i,k);
    }
  }
  
  Rcpp::List output = Rcpp::List::create(pNoDetGivenAC, pNoDetOcc);
  return output;
}


// [[Rcpp::export]]
Rcpp::List C_pNoDet_RW(int K, Rcpp::NumericVector tr_b4_occ,
                         int nCells, arma::vec p_init_state,
                         arma::cube pMove, arma::mat p_no_det_occ){
  arma::mat p_no_det_tot(nCells, K);
  arma::vec temp_1(nCells);
  arma::vec temp_2(nCells);
  
  if(tr_b4_occ(0) > 0){
    p_no_det_tot.col(0) = pMove.slice(0) * p_init_state % p_no_det_occ.col(0);
  } else{
    p_no_det_tot.col(0) = p_init_state % p_no_det_occ.col(0);
  }
  
  for(int k = 1; k < K; k++){
    if(tr_b4_occ(k) > 1){
      temp_1 = p_no_det_tot.col(k-1);
      for(int j = 0; j < tr_b4_occ(k); j++){
        temp_2 = pMove.slice(k) * temp_1;     //SHOULD THIS BE .SLICE(k-1)? MATTERS IF TEMPORAL VARIATION IN TRANSITIONS
        temp_1 = temp_2;
      }
      p_no_det_tot.col(k) = temp_1 % p_no_det_occ.col(k);
    } else{
      p_no_det_tot.col(k) = pMove.slice(k) * p_no_det_tot.col(k-1) % p_no_det_occ.col(k);
    }
  }
  
  double pDot = 1 - accu(p_no_det_tot.col(K-1));
  
  Rcpp::List output = Rcpp::List::create(pDot, p_no_det_tot);
  return output;
}


// [[Rcpp::export]]
Rcpp::List C_pNoDet_RWAC(int K, Rcpp::NumericVector tr_b4_occ,
                         int nCells, int nACCells, arma::mat p_init_state_full,
                         arma::cube pMove, arma::mat p_no_det_occ){
  arma::cube p_no_det_tot(nCells, nACCells, K);
  arma::mat temp_1(nCells, nACCells);
  arma::mat temp_2(nCells, nACCells);

  if(tr_b4_occ(0) > 0){
    for(int i = 0; i < nACCells; i++){
      p_no_det_tot.slice(0).col(i) = pMove.slice(i) * p_init_state_full.col(i) % p_no_det_occ.col(0); 
    }
  } else{
    for(int i = 0; i < nACCells; i++){
      p_no_det_tot.slice(0).col(i) = p_init_state_full.col(i) % p_no_det_occ.col(0);
    }
  }
  
  for(int k = 1; k < K; k++){
    if(tr_b4_occ(k) > 1){
      temp_1 = p_no_det_tot.slice(k-1);
      for(int j = 0; j < tr_b4_occ(k); j++){
        for(int i = 0; i < nACCells; i++){
          temp_2.col(i) = pMove.slice(i) * temp_1.col(i);
        }
       temp_1 = temp_2;
      }
      for(int i = 0; i < nACCells; i++){
        p_no_det_tot.slice(k).col(i) = temp_1.col(i) % p_no_det_occ.col(k);
        //for(int j = 0; j < nCells; j++){
        //  p_no_det_tot(j,i,k) = temp_1(j,i) * p_no_det_occ(j,k);
        //}
      }
    } else{
      for(int i = 0; i < nACCells; i++){
        p_no_det_tot.slice(k).col(i) = pMove.slice(i) * p_no_det_tot.slice(k-1).col(i) % p_no_det_occ.col(k);
      }
    }
  }
  
  double pDot = 1 - accu(p_no_det_tot.slice(K-1));
  
  Rcpp::List output = Rcpp::List::create(pDot, p_no_det_tot);
  return output;
  
}


// [[Rcpp::export]]
Rcpp::List C_pNoDet_CorrelatedRW(int K, Rcpp::NumericVector tr_b4_occ,
                         int nCells, arma::mat p_init_state_full,
                         arma::cube pMove, arma::mat p_no_det_occ){
  arma::cube p_no_det_tot(nCells, nCells, K);
  arma::mat temp_1(nCells, nCells);
  arma::mat temp_2(nCells, nCells);
  
  if(tr_b4_occ(0) > 0){
    for(int i = 0; i < nCells; i++){
      p_no_det_tot.slice(0).col(i) = pMove.slice(i) * p_init_state_full.row(i).t() % p_no_det_occ.col(0); 
    }
  } else{
    for(int i = 0; i < nCells; i++){
      p_no_det_tot.slice(0).col(i) = p_init_state_full.col(i) % p_no_det_occ.col(0);
    }
  }
  
  for(int k = 1; k < K; k++){
    if(tr_b4_occ(k) > 1){
      temp_1 = p_no_det_tot.slice(k-1);
      for(int j = 0; j < tr_b4_occ(k); j++){
        for(int i = 0; i < nCells; i++){
          temp_2.col(i) = pMove.slice(k-1) * temp_1.row(i).t();
        }
        temp_1 = temp_2;
      }
      for(int i = 0; i < nCells; i++){
        p_no_det_tot.slice(k).col(i) = temp_1.col(i) % p_no_det_occ.col(k);
        //for(int j = 0; j < nCells; j++){
        //  p_no_det_tot(j,i,k) = temp_1(j,i) * p_no_det_occ(j,k);
        //}
      }
    } else{
      for(int i = 0; i < nCells; i++){
        p_no_det_tot.slice(k).col(i) = pMove.slice(i) * p_no_det_tot.slice(k-1).row(i).t() % p_no_det_occ.col(k);
      }
    }
  }
  
  double pDot = 1 - accu(p_no_det_tot.slice(K-1));
  
  Rcpp::List output = Rcpp::List::create(pDot, p_no_det_tot);
  return output;
  
}


// [[Rcpp::export]]
arma::vec C_indLklhd_RandomWalk(int N, int K, int Nstates, arma::vec tr_b4_occ,
                            arma::vec mu_init, 
                            arma::mat y, arma::mat y_pix, 
                            arma::mat detDists, arma::mat y_platform,
                            arma::cube y_ds_covar, arma::cube y_g0_covar,
                            arma::cube y_hs_covar, arma::cube y_hb_covar,
                            int n_ds_covars, int n_g0_covars,
                            int n_hs_covars, int n_hb_covars,
                            arma::vec g0s, arma::vec sigmaDets,
                            arma::vec Hsigmas, arma::vec Hbetas,
                            arma::mat p_no_det_occ, arma::cube t_mat){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        state_vec = t_mat.slice(occ) * state_vec_V2;
        state_vec_V2 = state_vec;
        occ += 1;
      }
    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(y_pix(i, 0) - 1) = state_vec(y_pix(i, 0) - 1) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(y_pix(i, 0) - 1) = state_vec(y_pix(i, 0) - 1) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        state_vec_V2 = t_mat.slice(occ) * state_vec;
        state_vec = state_vec_V2;
        occ += 1;
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(y_pix(i, kt) - 1) = state_vec(y_pix(i, kt) - 1) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(y_pix(i, kt) - 1) = state_vec(y_pix(i, kt) - 1) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
double C_indLklhd_AC(int N, int K, int Nstates, 
                  arma::vec mu_init, arma::mat y, arma::mat y_pix, arma::mat detDists, arma::mat y_platform,
                  arma::cube y_ds_covar, arma::cube y_g0_covar, arma::cube y_hs_covar, arma::cube y_hb_covar,
                  int n_ds_covars, int n_g0_covars, int n_hs_covars, int n_hb_covars,
                  arma::vec g0s, arma::vec sigmaDets, arma::vec Hsigmas, arma::vec Hbetas,
                  arma::mat p_no_det_occ, arma::mat p_move){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 

    state_vec = state_vec_V2;

    for(int kt = 0; kt < K; kt++){      // loop through remaining sampling occasions
      
      state_vec_V2.zeros();     // zero out while maintaining dimension,
      
      if(y(i,kt) == 1){        // DETECTED!
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          for(int st = 0; st < Nstates; st++){
            state_vec_V2(st) = g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))) * 
              p_move(y_pix(i,kt) - 1, st) * 
              state_vec(st);
          }

        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          for(int st = 0; st < Nstates; st++){
            state_vec_V2(st) = g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) * 
              p_move(y_pix(i,kt) - 1, st) * 
              state_vec(st);
          }
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
  }
  
  return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
arma::vec C_indLklhd_RWAC(int N, int K, int nCells, int nACCells, 
                          arma::vec tr_b4_occ,
                          arma::vec mu_init, arma::mat p_init_state_full,
                          arma::mat y, arma::mat y_pix, 
                          arma::mat detDists, arma::mat y_platform,
                          arma::cube y_ds_covar, arma::cube y_g0_covar,
                          arma::cube y_hs_covar, arma::cube y_hb_covar,
                          int n_ds_covars, int n_g0_covars,
                          int n_hs_covars, int n_hb_covars,
                          arma::vec g0s, arma::vec sigmaDets,
                          arma::vec Hsigmas, arma::vec Hbetas,
                          arma::mat p_no_det_occ, arma::cube t_mat){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  arma::mat state_mat(nCells,nACCells);
  double s_v_sum;
  arma::mat state_mat_V2(nCells, nACCells);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = p_init_state_full;
    state_mat_V2 = p_init_state_full;
    s_v_sum = accu(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    int occ = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        for(int ac = 0; ac < nACCells; ac++){
          state_mat.col(ac) = t_mat.slice(ac) * state_mat_V2.col(ac);
        }
        state_mat_V2 = state_mat;
        occ += 1;
      }
    } else{
      state_mat = state_mat_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_mat_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_mat_V2.row(y_pix(i, 0) - 1) = state_mat.row(y_pix(i, 0) - 1) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_mat_V2.row(y_pix(i, 0) - 1) = state_mat.row(y_pix(i, 0) - 1) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int ac = 0; ac < nACCells; ac++){
        state_mat_V2.col(ac) = state_mat.col(ac) % p_no_det_occ.col(0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = accu(state_mat_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_mat = state_mat_V2 / s_v_sum;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int ac = 0; ac < nACCells; ac++){
          state_mat_V2.col(ac) = t_mat.slice(ac) * state_mat.col(ac);
        }
        state_mat = state_mat_V2;
        occ += 1;
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_mat_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_mat_V2.row(y_pix(i, kt) - 1) = state_mat.row(y_pix(i, kt) - 1) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_mat_V2.row(y_pix(i, kt) - 1) = state_mat.row(y_pix(i, kt) - 1) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        for(int ac = 0; ac < nACCells; ac++){
          state_mat_V2.col(ac) = state_mat.col(ac) % p_no_det_occ.col(kt);
        }
        
      }
      
      s_v_sum = accu(state_mat_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_mat = state_mat_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_mat, state_mat_V2, s_v_sum);
  //return output;
  
}


// [[Rcpp::export]]
arma::vec C_indLklhd_CorrelatedRW(int N, int K, int nCells, 
                          arma::vec tr_b4_occ,
                          arma::vec mu_init, arma::mat p_init_state_full,
                          arma::mat y, arma::mat y_pix, 
                          arma::mat detDists, arma::mat y_platform,
                          arma::cube y_ds_covar, arma::cube y_g0_covar,
                          arma::cube y_hs_covar, arma::cube y_hb_covar,
                          int n_ds_covars, int n_g0_covars,
                          int n_hs_covars, int n_hb_covars,
                          arma::vec g0s, arma::vec sigmaDets,
                          arma::vec Hsigmas, arma::vec Hbetas,
                          arma::mat p_no_det_occ, arma::cube t_mat){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  arma::mat state_mat(nCells,nCells);
  double s_v_sum;
  arma::mat state_mat_V2(nCells, nCells);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = p_init_state_full;
    state_mat_V2 = p_init_state_full;
    s_v_sum = accu(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    int occ = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        for(int k = 0; k < nCells; k++){
          state_mat.col(k) = t_mat.slice(k) * state_mat_V2.row(k).t();
        }
        state_mat_V2 = state_mat;
        occ += 1;
      }
    } else{
      state_mat = state_mat_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_mat_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_mat_V2.row(y_pix(i, 0) - 1) = state_mat.row(y_pix(i, 0) - 1) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_mat_V2.row(y_pix(i, 0) - 1) = state_mat.row(y_pix(i, 0) - 1) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int k = 0; k < nCells; k++){
        state_mat_V2.col(k) = state_mat.col(k) % p_no_det_occ.col(0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = accu(state_mat_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_mat = state_mat_V2 / s_v_sum;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int k = 0; k < nCells; k++){
          state_mat_V2.col(k) = t_mat.slice(k) * state_mat.row(k).t();
        }
        state_mat = state_mat_V2;
        occ += 1;
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_mat_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_mat_V2.row(y_pix(i, kt) - 1) = state_mat.row(y_pix(i, kt) - 1) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_mat_V2.row(y_pix(i, kt) - 1) = state_mat.row(y_pix(i, kt) - 1) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        for(int k = 0; k < nCells; k++){
          state_mat_V2.col(k) = state_mat.col(k) % p_no_det_occ.col(kt);
        }
        
      }
      
      s_v_sum = accu(state_mat_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_mat = state_mat_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_mat, state_mat_V2, s_v_sum);
  //return output;
  
}



// [[Rcpp::export]]
arma::vec C_indLklhd_RandomWalk_RWGroup(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                        arma::vec mu_init, 
                                        arma::mat y, arma::mat y_pix, 
                                        arma::mat detDists, arma::mat y_platform,
                                        arma::cube y_ds_covar, arma::cube y_g0_covar,
                                        arma::cube y_hs_covar, arma::cube y_hb_covar,
                                        int n_ds_covars, int n_g0_covars,
                                        int n_hs_covars, int n_hb_covars,
                                        arma::vec g0s, arma::vec sigmaDets,
                                        arma::vec Hsigmas, arma::vec Hbetas,
                                        arma::mat p_no_det_occ, arma::cube t_mat,
                                        arma::vec g0s_wCalf, arma::vec y_group,
                                        arma::mat p_no_det_occ_wCalf, arma::cube t_mat_wCalf,
                                        arma::vec p_birth, arma::vec y_sex,
                                        arma::uvec other_states, arma::uvec preg_states,
                                        arma::uvec wCalf_states, arma::uvec male_states){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double g0_wCalf = 1 / (1 + exp(-g0s_wCalf(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  int Nstates_sub = Nstates/4;
  arma::vec state_vec_V2(Nstates);
  arma::vec fromPreg(Nstates_sub);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    
    if(y_sex(i) == 1){    // female
      
      state_vec(male_states).zeros();
      state_vec_V2(male_states).zeros();
      
      if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
        for(int tr = 0; tr < tr_b4_occ(0); tr++){
          state_vec(other_states) = t_mat.slice(occ) * state_vec_V2(other_states);
          fromPreg = t_mat.slice(occ) * state_vec_V2(preg_states);
          state_vec(preg_states) = fromPreg * (1 - p_birth(occ));
          state_vec(wCalf_states) = fromPreg * p_birth(occ) + t_mat_wCalf.slice(occ) * state_vec_V2(wCalf_states);
          state_vec_V2 = state_vec;
          occ += 1;
        }
      } else{
        state_vec(other_states) = state_vec_V2(other_states);
        state_vec(preg_states) = state_vec_V2(preg_states);
        state_vec(wCalf_states) = state_vec_V2(wCalf_states);
      }
      
      if(y(i,0) == 1){        // DETECTED!
        
        state_vec_V2.zeros();   // zero out while maintaining dimension
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,0) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(other_states(y_pix(i, 0) - 1)) = state_vec(other_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(preg_states(y_pix(i, 0) - 1)) = state_vec(preg_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, 0) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(other_states(y_pix(i, 0) - 1)) = state_vec(other_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          state_vec_V2(preg_states(y_pix(i, 0) - 1)) = state_vec(preg_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          
        }
        
      } else if(y(i,0) == 2){
        
        state_vec_V2.zeros();   // zero out while maintaining dimension
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s_wCalf(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s_wCalf(ng0) * y_g0_covar(i, 0, ng0-1);
          }
          g0_wCalf = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,0) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(wCalf_states(y_pix(i, 0) - 1)) = state_vec(wCalf_states(y_pix(i, 0) - 1)) * g0_wCalf * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, 0) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(wCalf_states(y_pix(i, 0) - 1)) = state_vec(wCalf_states(y_pix(i, 0) - 1)) * g0_wCalf * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          
        }
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates_sub; st++){
          state_vec_V2(other_states(st)) = state_vec(other_states(st)) * p_no_det_occ(st,0);
          state_vec_V2(preg_states(st)) = state_vec(preg_states(st)) * p_no_det_occ(st,0);
          state_vec_V2(wCalf_states(st)) = state_vec(wCalf_states(st)) * p_no_det_occ_wCalf(st,0);
        }
        
      }
      
      // add to log-likelihood
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
        
        for(int tr = 0; tr < tr_b4_occ(kt); tr++){
          state_vec_V2(other_states) = t_mat.slice(occ) * state_vec(other_states);
          fromPreg = t_mat.slice(occ) * state_vec(preg_states);
          state_vec_V2(preg_states) = fromPreg * (1 - p_birth(occ));
          state_vec_V2(wCalf_states) = fromPreg * p_birth(occ) + t_mat_wCalf.slice(occ) * state_vec(wCalf_states);
          state_vec = state_vec_V2;
          occ += 1;
        }
        
        if(y(i,kt) == 1){        // DETECTED!
          
          state_vec_V2.zeros();     // zero out while maintaining dimension,
          
          if(n_g0_covars > 0){    // construct detection g0 if there are covariates
            ut_g0 = g0s(0);
            for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
              ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
            }
            g0 = 1 / (1 + exp(-ut_g0));
          }
          
          if(y_platform(i,kt) == 0){   //Half-norm
            
            if(n_ds_covars > 0){    // construct detection sigma if there are covariates
              ut_sigmaDet = sigmaDets(0);
              for(int nds = 1; nds <= n_ds_covars; nds++){
                ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
              }
              sigmaDet = exp(ut_sigmaDet);
            }
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(other_states(y_pix(i, kt) - 1)) = state_vec(other_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            state_vec_V2(preg_states(y_pix(i, kt) - 1)) = state_vec(preg_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            
          } else if(y_platform(i, kt) == 1){   //Hazard function
            
            if(n_hs_covars > 0){    // construct detection sigma if there are covariates
              ut_Hsigma = Hsigmas(0);
              for(int nhs = 1; nhs <= n_hs_covars; nhs++){
                ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
              }
              Hsigma = exp(ut_Hsigma);
            }
            
            if(n_hb_covars > 0){    // construct detection sigma if there are covariates
              ut_Hbeta = Hbetas(0);
              for(int nhb = 1; nhb <= n_hb_covars; nhb++){
                ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
              }
              Hbeta = exp(ut_Hbeta);
            }
            
            //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
            //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
            //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(other_states(y_pix(i, kt) - 1)) = state_vec(other_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            state_vec_V2(preg_states(y_pix(i, kt) - 1)) = state_vec(preg_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            
          }  
          
        } else if(y(i,kt) == 2){
          
          state_vec_V2.zeros();     // zero out while maintaining dimension,
          
          if(n_g0_covars > 0){    // construct detection g0 if there are covariates
            ut_g0 = g0s_wCalf(0);
            for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
              ut_g0 += g0s_wCalf(ng0) * y_g0_covar(i, kt, ng0-1);
            }
            g0_wCalf = 1 / (1 + exp(-ut_g0));
          }
          
          if(y_platform(i,kt) == 0){   //Half-norm
            
            if(n_ds_covars > 0){    // construct detection sigma if there are covariates
              ut_sigmaDet = sigmaDets(0);
              for(int nds = 1; nds <= n_ds_covars; nds++){
                ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
              }
              sigmaDet = exp(ut_sigmaDet);
            }
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(wCalf_states(y_pix(i, kt) - 1)) = state_vec(wCalf_states(y_pix(i, kt) - 1)) * g0_wCalf * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            
          } else if(y_platform(i, kt) == 1){   //Hazard function
            
            if(n_hs_covars > 0){    // construct detection sigma if there are covariates
              ut_Hsigma = Hsigmas(0);
              for(int nhs = 1; nhs <= n_hs_covars; nhs++){
                ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
              }
              Hsigma = exp(ut_Hsigma);
            }
            
            if(n_hb_covars > 0){    // construct detection sigma if there are covariates
              ut_Hbeta = Hbetas(0);
              for(int nhb = 1; nhb <= n_hb_covars; nhb++){
                ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
              }
              Hbeta = exp(ut_Hbeta);
            }
            
            //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
            //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
            //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(wCalf_states(y_pix(i, kt) - 1)) = state_vec(wCalf_states(y_pix(i, kt) - 1)) * g0_wCalf * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            
          }
          
        }else{   // Not detected
          
          //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
          for(int st = 0; st < Nstates_sub; st++){
            state_vec_V2(other_states(st)) = state_vec(other_states(st)) * p_no_det_occ(st,kt);
            state_vec_V2(preg_states(st)) = state_vec(preg_states(st)) * p_no_det_occ(st,kt);
            state_vec_V2(wCalf_states(st)) = state_vec(wCalf_states(st)) * p_no_det_occ_wCalf(st,kt);
          }
          
        }
        
        s_v_sum = sum(state_vec_V2);
        Ind_logProb(i) += log(s_v_sum);
        state_vec = state_vec_V2 / s_v_sum;
        
        //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
        //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
        //}
        
        //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
        
      }
      
    } else if(y_sex(i) == 2){ // male
      
      state_vec(other_states).zeros();
      state_vec(preg_states).zeros();
      state_vec(wCalf_states).zeros();
      state_vec_V2(other_states).zeros();
      state_vec_V2(preg_states).zeros();
      state_vec_V2(wCalf_states).zeros();
      
      if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
        for(int tr = 0; tr < tr_b4_occ(0); tr++){
          state_vec(male_states) = t_mat.slice(occ) * state_vec_V2(male_states);
          state_vec_V2 = state_vec;
          occ += 1;
        }
      } else{
        state_vec(male_states) = state_vec_V2(male_states);
      }
      
      if(y(i,0) == 1){        // DETECTED!
        
        state_vec_V2.zeros();   // zero out while maintaining dimension
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,0) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(male_states(y_pix(i, 0) - 1)) = state_vec(male_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, 0) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(male_states(y_pix(i, 0) - 1)) = state_vec(male_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          
        }
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates_sub; st++){
          state_vec_V2(male_states(st)) = state_vec(male_states(st)) * p_no_det_occ(st,0);
        }
        
      }
      
      // add to log-likelihood
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
        
        for(int tr = 0; tr < tr_b4_occ(kt); tr++){
          state_vec_V2(male_states) = t_mat.slice(occ) * state_vec(male_states);
          state_vec = state_vec_V2;
          occ += 1;
        }
        
        if(y(i,kt) == 1){        // DETECTED!
          
          state_vec_V2.zeros();     // zero out while maintaining dimension,
          
          if(n_g0_covars > 0){    // construct detection g0 if there are covariates
            ut_g0 = g0s(0);
            for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
              ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
            }
            g0 = 1 / (1 + exp(-ut_g0));
          }
          
          if(y_platform(i,kt) == 0){   //Half-norm
            
            if(n_ds_covars > 0){    // construct detection sigma if there are covariates
              ut_sigmaDet = sigmaDets(0);
              for(int nds = 1; nds <= n_ds_covars; nds++){
                ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
              }
              sigmaDet = exp(ut_sigmaDet);
            }
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(male_states(y_pix(i, kt) - 1)) = state_vec(male_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            
          } else if(y_platform(i, kt) == 1){   //Hazard function
            
            if(n_hs_covars > 0){    // construct detection sigma if there are covariates
              ut_Hsigma = Hsigmas(0);
              for(int nhs = 1; nhs <= n_hs_covars; nhs++){
                ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
              }
              Hsigma = exp(ut_Hsigma);
            }
            
            if(n_hb_covars > 0){    // construct detection sigma if there are covariates
              ut_Hbeta = Hbetas(0);
              for(int nhb = 1; nhb <= n_hb_covars; nhb++){
                ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
              }
              Hbeta = exp(ut_Hbeta);
            }
            
            //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
            //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
            //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
            
            // only update new_state_vec where individual was detected -- probability of being in any other state is 0
            state_vec_V2(male_states(y_pix(i, kt) - 1)) = state_vec(male_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            
          }  
          
        } else{   // Not detected
          
          //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
          for(int st = 0; st < Nstates_sub; st++){
            state_vec_V2(male_states(st)) = state_vec(male_states(st)) * p_no_det_occ(st,kt);
          }
          
        }
        
        s_v_sum = sum(state_vec_V2);
        Ind_logProb(i) += log(s_v_sum);
        state_vec = state_vec_V2 / s_v_sum;
        
        //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
        //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
        //}
        
        //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
        
      }
      
    }
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
arma::vec C_indLklhd_3states(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                        arma::vec mu_init, 
                                        arma::mat y, arma::mat y_pix, 
                                        arma::mat detDists, arma::mat y_platform,
                                        arma::cube y_ds_covar, arma::cube y_g0_covar,
                                        arma::cube y_hs_covar, arma::cube y_hb_covar,
                                        int n_ds_covars, int n_g0_covars,
                                        int n_hs_covars, int n_hb_covars,
                                        arma::vec g0s, arma::vec sigmaDets,
                                        arma::vec Hsigmas, arma::vec Hbetas,
                                        arma::mat p_no_det_occ, arma::mat p_move_reside,
                                        arma::mat p_move_south, arma::mat p_move_north,
                                        arma::vec emigration, arma::vec first_loc_south,
                                        arma::vec gamma, arma::vec phi,
                                        arma::vec p_reside,
                                        arma::uvec south_states, arma::uvec reside_states,
                                        arma::uvec north_states, int out_states,
                                        int nCells){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  arma::vec temp_vec(nCells);
  arma::vec immigrating_temp(nCells);
  arma::vec residing_temp(nCells);
  arma::vec movingNorth_temp(nCells);
  
  arma::vec remainNorth(nCells);
  arma::vec remainSouth(nCells);
  for(int c = 0; c < nCells; c++){
    remainNorth(c) = 1 - emigration(c);
    remainSouth(c) = 1 - p_reside(c);
  }
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        for(int c = 0; c < nCells; c++){
          immigrating_temp(c) = state_vec_V2(out_states) * gamma(occ) * first_loc_south(c);
          residing_temp(c) = state_vec_V2(reside_states(c)) * phi(occ);
          movingNorth_temp(c) = state_vec_V2(reside_states(c)) * (1 - phi(occ));
        }
        
        state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
          accu(state_vec_V2(north_states) % emigration);
        state_vec(south_states) = p_move_south * (state_vec_V2(south_states) % remainSouth) + 
          immigrating_temp;
        state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
          p_move_reside * residing_temp;
        state_vec(north_states) = p_move_north * (state_vec_V2(north_states) % remainNorth) + 
          p_move_north * movingNorth_temp;
        
        //state_vec(south_states) = state_vec_V2(out_states) * gamma(occ) * first_loc_south + // may be an issue with differing lengths
        //  p_move_south * (state_vec_V2(south_states) % (1-p_reside));
        //state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
        //  p_move_reside * (state_vec_V2(reside_states) * phi(occ));
        //state_vec(north_states) = p_move_north * (state_vec_V2(reside_states) * (1 - phi(occ))) + 
        //  p_move_north * (state_vec_V2(north_states) % (1 - emigration));
        
        state_vec_V2 = state_vec;
        occ += 1;
      }
    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(south_states(y_pix(i, 0) - 1)) = state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(north_states(y_pix(i, 0) - 1)) = state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
        state_vec_V2(reside_states(y_pix(i, 0) - 1)) = state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(south_states(y_pix(i, 0) - 1)) = state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(north_states(y_pix(i, 0) - 1)) = state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(reside_states(y_pix(i, 0) - 1)) = state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < nCells; st++){
        state_vec_V2(south_states(st)) = state_vec(south_states(st)) * p_no_det_occ(st,0);
        state_vec_V2(north_states(st)) = state_vec(north_states(st)) * p_no_det_occ(st,0);
        state_vec_V2(reside_states(st)) = state_vec(reside_states(st)) * p_no_det_occ(st,0);
      }
      state_vec_V2(out_states) = state_vec(out_states);
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    state_vec_V2 = state_vec;
    
    //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
    
    //Rcpp::Rcout << "state_vec sum = " << accu(state_vec) << "\n";
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int c = 0; c < nCells; c++){
          immigrating_temp(c) = state_vec_V2(out_states) * gamma(occ) * first_loc_south(c);
          residing_temp(c) = state_vec_V2(reside_states(c)) * phi(occ);
          movingNorth_temp(c) = state_vec_V2(reside_states(c)) * (1 - phi(occ));
        }
        
        state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
          accu(state_vec_V2(north_states) % emigration);
        
        //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
        
        state_vec(south_states) = p_move_south * (state_vec_V2(south_states) % remainSouth) + 
          immigrating_temp;
        
        //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec(south_states)) << "\n";
        
        state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
          p_move_reside * residing_temp;
        
        //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec(reside_states)) << "\n";
        
        state_vec(north_states) = p_move_north * (state_vec_V2(north_states) % remainNorth) + 
          p_move_north * movingNorth_temp;
        
        //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec(north_states)) << "\n";
        
        //state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
        //  accu(state_vec_V2(north_states) % emigration);
        //state_vec(south_states) = state_vec_V2(out_states) * gamma(occ) * first_loc_south + // may be an issue with differing lengths
        //  p_move_south * (state_vec_V2(south_states) % (1-p_reside));
        //state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
        //  p_move_reside * (state_vec_V2(reside_states) * phi(occ));
        //state_vec(north_states) = p_move_north * (state_vec_V2(reside_states) * (1 - phi(occ))) + 
        //  p_move_north * (state_vec_V2(north_states) % (1 - emigration));
        
        state_vec_V2 = state_vec;
        occ += 1;
      }
      
      //Rcpp::Rcout << "state_vec " << kt+1 << " sum = " << accu(state_vec) << "\n";
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(south_states(y_pix(i, kt) - 1)) = state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(north_states(y_pix(i, kt) - 1)) = state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(reside_states(y_pix(i, kt) - 1)) = state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(south_states(y_pix(i, kt) - 1)) = state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(north_states(y_pix(i, kt) - 1)) = state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(reside_states(y_pix(i, kt) - 1)) = state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
        //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec(south_states)) << "\n";
        //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec(reside_states)) << "\n";
        //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec(north_states)) << "\n";
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < nCells; st++){
          state_vec_V2(south_states(st)) = state_vec(south_states(st)) * p_no_det_occ(st,kt);
          state_vec_V2(north_states(st)) = state_vec(north_states(st)) * p_no_det_occ(st,kt);
          state_vec_V2(reside_states(st)) = state_vec(reside_states(st)) * p_no_det_occ(st,kt);
        }
        state_vec_V2(out_states) = state_vec(out_states);
        
      }
      //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec_V2(out_states) << "\n";
      //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec_V2(south_states)) << "\n";
      //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec_V2(reside_states)) << "\n";
      //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec_V2(north_states)) << "\n";
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      state_vec_V2 = state_vec;
      
      //Rcpp::Rcout << "logProb " << kt+1 << " = " << Ind_logProb(i) << "\n";
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }

    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, s_v_sum, sigmaDet, g0);
  //return output;
  
}



// [[Rcpp::export]]
arma::vec C_indLklhd_3states_moveDirObs(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                        arma::vec mu_init, 
                                        arma::mat y, arma::mat y_pix, 
                                        arma::mat detDists, arma::mat y_platform,
                                        arma::cube y_ds_covar, arma::cube y_g0_covar,
                                        arma::cube y_hs_covar, arma::cube y_hb_covar,
                                        int n_ds_covars, int n_g0_covars,
                                        int n_hs_covars, int n_hb_covars,
                                        arma::vec g0s, arma::vec sigmaDets,
                                        arma::vec Hsigmas, arma::vec Hbetas,
                                        arma::mat p_no_det_occ, arma::mat p_move_reside,
                                        arma::mat p_move_south, arma::mat p_move_north,
                                        arma::vec emigration, arma::vec first_loc_south,
                                        arma::vec gamma, arma::vec phi,
                                        arma::vec p_reside,
                                        arma::uvec south_states, arma::uvec reside_states,
                                        arma::uvec north_states, int out_states,
                                        int nCells, double migrate_obs,
                                        double no_move_migrate_obs, double ew_migrate_obs, 
                                        double no_move_reside_obs, double news_reside_obs,
                                        arma::mat y_dir){
                             
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  arma::vec temp_vec(nCells);
  arma::vec immigrating_temp(nCells);
  arma::vec residing_temp(nCells);
  arma::vec movingNorth_temp(nCells);
  
  arma::vec remainNorth(nCells);
  arma::vec remainSouth(nCells);
  for(int c = 0; c < nCells; c++){
    remainNorth(c) = 1 - emigration(c);
    remainSouth(c) = 1 - p_reside(c);
  }
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        for(int c = 0; c < nCells; c++){
          immigrating_temp(c) = state_vec_V2(out_states) * gamma(occ) * first_loc_south(c);
          residing_temp(c) = state_vec_V2(reside_states(c)) * phi(occ);
          movingNorth_temp(c) = state_vec_V2(reside_states(c)) * (1 - phi(occ));
        }
        
        state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
          accu(state_vec_V2(north_states) % emigration);
        state_vec(south_states) = p_move_south * (state_vec_V2(south_states) % remainSouth) + 
          immigrating_temp;
        state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
          p_move_reside * residing_temp;
        state_vec(north_states) = p_move_north * (state_vec_V2(north_states) % remainNorth) + 
          p_move_north * movingNorth_temp;
        
        //state_vec(south_states) = state_vec_V2(out_states) * gamma(occ) * first_loc_south + // may be an issue with differing lengths
        //  p_move_south * (state_vec_V2(south_states) % (1-p_reside));
        //state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
        //  p_move_reside * (state_vec_V2(reside_states) * phi(occ));
        //state_vec(north_states) = p_move_north * (state_vec_V2(reside_states) * (1 - phi(occ))) + 
        //  p_move_north * (state_vec_V2(north_states) % (1 - emigration));
        
        state_vec_V2 = state_vec;
        occ += 1;
      }
    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        //state_vec_V2(south_states(y_pix(i, 0) - 1)) = state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        //state_vec_V2(north_states(y_pix(i, 0) - 1)) = state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
        //state_vec_V2(reside_states(y_pix(i, 0) - 1)) = state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
        
        if(y_dir(i,0) == 0){ 
          //observed direction = no direction
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = no_move_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = no_move_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = no_move_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          
        } else if(y_dir(i,0) == 1){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          
        } else if(y_dir(i,0) == 2){
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          
        } else if(y_dir(i,0) == 3){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          
        } else if(y_dir(i,0) == 4){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2)));
          
        }
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        //state_vec_V2(south_states(y_pix(i, 0) - 1)) = state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        //state_vec_V2(north_states(y_pix(i, 0) - 1)) = state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        //state_vec_V2(reside_states(y_pix(i, 0) - 1)) = state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
        if(y_dir(i,0) == 0){ 
          //observed direction = no direction
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = no_move_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = no_move_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = no_move_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          
        } else if(y_dir(i,0) == 1){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          
        } else if(y_dir(i,0) == 2){
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          
        } else if(y_dir(i,0) == 3){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
          state_vec_V2(north_states(y_pix(i, 0) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          
        } else if(y_dir(i,0) == 4){
          state_vec_V2(south_states(y_pix(i, 0) - 1)) = migrate_obs * state_vec(south_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          state_vec_V2(reside_states(y_pix(i, 0) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta)));
          
        }
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < nCells; st++){
        state_vec_V2(south_states(st)) = state_vec(south_states(st)) * p_no_det_occ(st,0);
        state_vec_V2(north_states(st)) = state_vec(north_states(st)) * p_no_det_occ(st,0);
        state_vec_V2(reside_states(st)) = state_vec(reside_states(st)) * p_no_det_occ(st,0);
      }
      state_vec_V2(out_states) = state_vec(out_states);
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    state_vec_V2 = state_vec;
    
    //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
    
    //Rcpp::Rcout << "state_vec sum = " << accu(state_vec) << "\n";
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int c = 0; c < nCells; c++){
          immigrating_temp(c) = state_vec_V2(out_states) * gamma(occ) * first_loc_south(c);
          residing_temp(c) = state_vec_V2(reside_states(c)) * phi(occ);
          movingNorth_temp(c) = state_vec_V2(reside_states(c)) * (1 - phi(occ));
        }
        
        state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
          accu(state_vec_V2(north_states) % emigration);
        
        //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
        
        state_vec(south_states) = p_move_south * (state_vec_V2(south_states) % remainSouth) + 
          immigrating_temp;
        
        //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec(south_states)) << "\n";
        
        state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
          p_move_reside * residing_temp;
        
        //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec(reside_states)) << "\n";
        
        state_vec(north_states) = p_move_north * (state_vec_V2(north_states) % remainNorth) + 
          p_move_north * movingNorth_temp;
        
        //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec(north_states)) << "\n";
        
        //state_vec(out_states) = state_vec_V2(out_states) * (1 - gamma(occ)) + 
        //  accu(state_vec_V2(north_states) % emigration);
        //state_vec(south_states) = state_vec_V2(out_states) * gamma(occ) * first_loc_south + // may be an issue with differing lengths
        //  p_move_south * (state_vec_V2(south_states) % (1-p_reside));
        //state_vec(reside_states) = p_move_reside * (state_vec_V2(south_states) % p_reside) + 
        //  p_move_reside * (state_vec_V2(reside_states) * phi(occ));
        //state_vec(north_states) = p_move_north * (state_vec_V2(reside_states) * (1 - phi(occ))) + 
        //  p_move_north * (state_vec_V2(north_states) % (1 - emigration));
        
        state_vec_V2 = state_vec;
        occ += 1;
      }
      
      //Rcpp::Rcout << "state_vec " << kt+1 << " sum = " << accu(state_vec) << "\n";
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          //state_vec_V2(south_states(y_pix(i, kt) - 1)) = state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          //state_vec_V2(north_states(y_pix(i, kt) - 1)) = state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          //state_vec_V2(reside_states(y_pix(i, kt) - 1)) = state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
          if(y_dir(i,kt) == 0){ 
            //observed direction = no direction
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = no_move_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = no_move_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = no_move_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            
          } else if(y_dir(i,kt) == 1){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            
          } else if(y_dir(i,kt) == 2){
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            
          } else if(y_dir(i,kt) == 3){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            
          } else if(y_dir(i,kt) == 4){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2)));
            
          }
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          //state_vec_V2(south_states(y_pix(i, kt) - 1)) = state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          //state_vec_V2(north_states(y_pix(i, kt) - 1)) = state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          //state_vec_V2(reside_states(y_pix(i, kt) - 1)) = state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
          if(y_dir(i,kt) == 0){ 
            //observed direction = no direction
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = no_move_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = no_move_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = no_move_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            
          } else if(y_dir(i,kt) == 1){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            
          } else if(y_dir(i,kt) == 2){
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            
          } else if(y_dir(i,kt) == 3){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
            state_vec_V2(north_states(y_pix(i, kt) - 1)) = ew_migrate_obs * state_vec(north_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            
          } else if(y_dir(i,kt) == 4){
            state_vec_V2(south_states(y_pix(i, kt) - 1)) = migrate_obs * state_vec(south_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            state_vec_V2(reside_states(y_pix(i, kt) - 1)) = news_reside_obs * state_vec(reside_states(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta)));
            
          }
          
        }  
        
      } else{   // Not detected
        
        //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec(out_states) << "\n";
        //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec(south_states)) << "\n";
        //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec(reside_states)) << "\n";
        //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec(north_states)) << "\n";
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < nCells; st++){
          state_vec_V2(south_states(st)) = state_vec(south_states(st)) * p_no_det_occ(st,kt);
          state_vec_V2(north_states(st)) = state_vec(north_states(st)) * p_no_det_occ(st,kt);
          state_vec_V2(reside_states(st)) = state_vec(reside_states(st)) * p_no_det_occ(st,kt);
        }
        state_vec_V2(out_states) = state_vec(out_states);
        
      }
      //Rcpp::Rcout << "state_vec(out_states) sum = " << state_vec_V2(out_states) << "\n";
      //Rcpp::Rcout << "state_vec(south_states) sum = " << accu(state_vec_V2(south_states)) << "\n";
      //Rcpp::Rcout << "state_vec(reside_states) sum = " << accu(state_vec_V2(reside_states)) << "\n";
      //Rcpp::Rcout << "state_vec(north_states) sum = " << accu(state_vec_V2(north_states)) << "\n";
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      state_vec_V2 = state_vec;
      
      //Rcpp::Rcout << "logProb " << kt+1 << " = " << Ind_logProb(i) << "\n";
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
arma::vec C_indLklhd_RandomWalk_Direction(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                          arma::vec mu_init, 
                                          arma::mat y, arma::mat y_pix, 
                                          arma::mat detDists, arma::mat y_platform,
                                          arma::cube y_ds_covar, arma::cube y_g0_covar,
                                          arma::cube y_hs_covar, arma::cube y_hb_covar,
                                          int n_ds_covars, int n_g0_covars,
                                          int n_hs_covars, int n_hb_covars,
                                          arma::vec g0s, arma::vec sigmaDets,
                                          arma::vec Hsigmas, arma::vec Hbetas,
                                          arma::mat p_no_det_occ, 
                                          arma::cube pMove_North, arma::cube pMove_West,
                                          arma::cube pMove_South, arma::cube pMove_East,
                                          arma::uvec north_cells, arma::uvec west_cells,
                                          arma::uvec south_cells, arma::uvec east_cells,
                                          int out1_cell, int out2_cell, arma::uvec cov_open_transitions,
                                          double p_nsns, double p_nsew, double p_nssn,
                                          double p_ewew, double p_ewns, double p_ewwe,
                                          arma::vec phi, arma::vec gamma, arma::vec D_bars,
                                          arma::mat D_xs, arma::vec open_transitions){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  int n_open_transitions = open_transitions.n_elem;
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    int next_open_transition = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        if(occ == open_transitions(next_open_transition)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = state_vec_V2(out1_cell) * (1 - gamma(occ));
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition += 1;
          
        } else{
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe);
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew);
          
          state_vec(out1_cell) = state_vec_V2(out1_cell);
          
          state_vec(out2_cell) = state_vec_V2(out2_cell);
          
        }
        
        state_vec_V2 = state_vec;
        occ += 1;
      }


    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(north_cells(y_pix(i, 0) - 1)) = state_vec(north_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(west_cells(y_pix(i, 0) - 1)) = state_vec(west_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(south_cells(y_pix(i, 0) - 1)) = state_vec(south_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(east_cells(y_pix(i, 0) - 1)) = state_vec(east_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(north_cells(y_pix(i, 0) - 1)) = state_vec(north_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(west_cells(y_pix(i, 0) - 1)) = state_vec(west_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(south_cells(y_pix(i, 0) - 1)) = state_vec(south_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(east_cells(y_pix(i, 0) - 1)) = state_vec(east_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    state_vec_V2 = state_vec;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        if(occ == open_transitions(next_open_transition) && next_open_transition < (n_open_transitions - 1)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = state_vec_V2(out1_cell) * (1 - gamma(occ));
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition += 1;


        } else if(occ == open_transitions(next_open_transition) && next_open_transition == (n_open_transitions - 1)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = 0;
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition = 0;
          
          
        } else{
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe);
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew);
          
          state_vec(out1_cell) = state_vec_V2(out1_cell);
          
          state_vec(out2_cell) = state_vec_V2(out2_cell);
          
        }
        
        state_vec_V2 = state_vec;
        occ += 1;
        //Rcpp::Rcout << "k = " << kt << "; occ = " << occ << "; next_open_transition = " << next_open_transition << "\n";
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(north_cells(y_pix(i, kt) - 1)) = state_vec(north_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(west_cells(y_pix(i, kt) - 1)) = state_vec(west_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(south_cells(y_pix(i, kt) - 1)) = state_vec(south_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(east_cells(y_pix(i, kt) - 1)) = state_vec(east_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(north_cells(y_pix(i, kt) - 1)) = state_vec(north_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(west_cells(y_pix(i, kt) - 1)) = state_vec(west_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(south_cells(y_pix(i, kt) - 1)) = state_vec(south_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(east_cells(y_pix(i, kt) - 1)) = state_vec(east_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      state_vec_V2 = state_vec;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
      
    }
    
    //Rcpp::Rcout << "next_open_transition: " << next_open_transition << "\n";
    
  }
  
  return Ind_logProb;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
double C_indLklhd_OpenAC(int N, int K, int Nstates, 
                         arma::vec mu_init, arma::mat y, arma::mat y_pix, arma::mat detDists, arma::mat y_platform,
                         arma::cube y_ds_covar, arma::cube y_g0_covar, arma::cube y_hs_covar, arma::cube y_hb_covar,
                         int n_ds_covars, int n_g0_covars, int n_hs_covars, int n_hb_covars,
                         arma::vec g0s, arma::vec sigmaDets, arma::vec Hsigmas, arma::vec Hbetas,
                         arma::mat p_no_det_occ, arma::cube t_arr, arma::uvec K_sec, arma::cube p_move, arma::uvec K_covar,
                         int nCells){
                     
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    
    state_vec = state_vec_V2;
    
    state_vec_V2.zeros();     // zero out while maintaining dimension,
    
    int kt = 0;
    
    if(y(i,kt) == 1){        // DETECTED!
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,kt) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        for(int st = 0; st < nCells; st++){
          state_vec_V2(st) = g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))) * 
            p_move(y_pix(i,kt) - 1, st, 0) * 
            state_vec(st);
        }
        
      } else if(y_platform(i, kt) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
        //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
        //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
        
        for(int st = 0; st < nCells; st++){
          state_vec_V2(st) = g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) * 
            p_move(y_pix(i,kt) - 1, st, 0) * 
            state_vec(st);
        }
        
      }  
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
      }
      
    }
    
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      if(K_sec(kt) != K_sec(kt-1)){
        state_vec_V2 = state_vec;
        state_vec = t_arr.slice(K_sec(kt-1)) * state_vec_V2;
      }
      
      state_vec_V2.zeros();     // zero out while maintaining dimension,
      
      if(y(i,kt) == 1){        // DETECTED!
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          for(int st = 0; st < nCells; st++){ 
            state_vec_V2(st) = g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))) * 
              p_move(y_pix(i,kt) - 1, st, K_covar(kt)) * 
              state_vec(st);
          }
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          for(int st = 0; st < nCells; st++){
            state_vec_V2(st) = g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) * 
              p_move(y_pix(i,kt) - 1, st, K_covar(kt)) * 
              state_vec(st);
          }
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
  }
  
  return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}



// [[Rcpp::export]]
arma::mat C_RandomWalk_ReturnState(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                    arma::vec mu_init, 
                                    arma::mat y, arma::mat y_pix, 
                                    arma::mat detDists, arma::mat y_platform,
                                    arma::cube y_ds_covar, arma::cube y_g0_covar,
                                    arma::cube y_hs_covar, arma::cube y_hb_covar,
                                    int n_ds_covars, int n_g0_covars,
                                    int n_hs_covars, int n_hb_covars,
                                    arma::vec g0s, arma::vec sigmaDets,
                                    arma::vec Hsigmas, arma::vec Hbetas,
                                    arma::mat p_no_det_occ, arma::cube t_mat){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  arma::mat State_Mat(Nstates, N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        state_vec = t_mat.slice(occ) * state_vec_V2;
        state_vec_V2 = state_vec;
        occ += 1;
      }
    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(y_pix(i, 0) - 1) = state_vec(y_pix(i, 0) - 1) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(y_pix(i, 0) - 1) = state_vec(y_pix(i, 0) - 1) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        state_vec_V2 = t_mat.slice(occ) * state_vec;
        state_vec = state_vec_V2;
        occ += 1;
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(y_pix(i, kt) - 1) = state_vec(y_pix(i, kt) - 1) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(y_pix(i, kt) - 1) = state_vec(y_pix(i, kt) - 1) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
    for(int state = 0; state < Nstates; state++){
      State_Mat(state, i) = state_vec(state);
    }
    
  }
  
  return State_Mat;
  
}



// [[Rcpp::export]]
arma::mat C_OpenAC_ReturnState(int N, int K, int Nstates, 
                               arma::vec mu_init, arma::mat y, arma::mat y_pix, arma::mat detDists, arma::mat y_platform,
                               arma::cube y_ds_covar, arma::cube y_g0_covar, arma::cube y_hs_covar, arma::cube y_hb_covar,
                               int n_ds_covars, int n_g0_covars, int n_hs_covars, int n_hb_covars,
                               arma::vec g0s, arma::vec sigmaDets, arma::vec Hsigmas, arma::vec Hbetas,
                               arma::mat p_no_det_occ, arma::cube t_arr, arma::uvec K_sec, arma::cube p_move, arma::uvec K_covar,
                               int nCells){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  arma::mat State_Mat(Nstates, N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    
    state_vec = state_vec_V2;
    
    state_vec_V2.zeros();     // zero out while maintaining dimension,
    
    int kt = 0;
    
    if(y(i,kt) == 1){        // DETECTED!
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,kt) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        for(int st = 0; st < nCells; st++){
          state_vec_V2(st) = g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))) * 
            p_move(y_pix(i,kt) - 1, st, 0) * 
            state_vec(st);
        }
        
      } else if(y_platform(i, kt) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
        //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
        //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
        
        for(int st = 0; st < nCells; st++){
          state_vec_V2(st) = g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) * 
            p_move(y_pix(i,kt) - 1, st, 0) * 
            state_vec(st);
        }
        
      }  
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
      }
      
    }
    
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      if(K_sec(kt) != K_sec(kt-1)){
        state_vec_V2 = state_vec;
        state_vec = t_arr.slice(K_sec(kt-1)) * state_vec_V2;
      }
      
      state_vec_V2.zeros();     // zero out while maintaining dimension,
      
      if(y(i,kt) == 1){        // DETECTED!
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          for(int st = 0; st < nCells; st++){ 
            state_vec_V2(st) = g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))) * 
              p_move(y_pix(i,kt) - 1, st, K_covar(kt)) * 
              state_vec(st);
          }
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          for(int st = 0; st < nCells; st++){
            state_vec_V2(st) = g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) * 
              p_move(y_pix(i,kt) - 1, st, K_covar(kt)) * 
              state_vec(st);
          }
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
    for(int st = 0; st < Nstates; st++){
      State_Mat(st,i) = state_vec(st);
    }
  }
  
  return State_Mat;
  
}


// [[Rcpp::export]]
arma::mat C_RandomWalk_Direction_ReturnState(int N, int K, int Nstates, arma::vec tr_b4_occ,
                                             arma::vec mu_init, 
                                             arma::mat y, arma::mat y_pix, 
                                             arma::mat detDists, arma::mat y_platform,
                                             arma::cube y_ds_covar, arma::cube y_g0_covar,
                                             arma::cube y_hs_covar, arma::cube y_hb_covar,
                                             int n_ds_covars, int n_g0_covars,
                                             int n_hs_covars, int n_hb_covars,
                                             arma::vec g0s, arma::vec sigmaDets,
                                             arma::vec Hsigmas, arma::vec Hbetas,
                                             arma::mat p_no_det_occ, 
                                             arma::cube pMove_North, arma::cube pMove_West,
                                             arma::cube pMove_South, arma::cube pMove_East,
                                             arma::uvec north_cells, arma::uvec west_cells,
                                             arma::uvec south_cells, arma::uvec east_cells,
                                             int out1_cell, int out2_cell, arma::uvec cov_open_transitions,
                                             double p_nsns, double p_nsew, double p_nssn,
                                             double p_ewew, double p_ewns, double p_ewwe,
                                             arma::vec phi, arma::vec gamma, arma::vec D_bars,
                                             arma::mat D_xs, arma::vec open_transitions){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec = mu_init;
  double s_v_sum;
  arma::vec state_vec_V2(Nstates);
  int n_open_transitions = open_transitions.n_elem;
  arma::mat State_Mat(Nstates, N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    Ind_logProb(i) = log(s_v_sum);    
    state_vec_V2 = state_vec / s_v_sum;  //Standardize state_vec to sum to one 
    int occ = 0;
    int next_open_transition = 0;
    
    if(tr_b4_occ(0) > 0){                       // Transitions before the first sampling occasion
      for(int tr = 0; tr < tr_b4_occ(0); tr++){
        if(occ == open_transitions(next_open_transition)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = state_vec_V2(out1_cell) * (1 - gamma(occ));
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition += 1;
          
        } else{
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe);
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew);
          
          state_vec(out1_cell) = state_vec_V2(out1_cell);
          
          state_vec(out2_cell) = state_vec_V2(out2_cell);
          
        }
        
        state_vec_V2 = state_vec;
        occ += 1;
      }
      
      
    } else{
      state_vec = state_vec_V2;
    }
    
    if(y(i,0) == 1){        // DETECTED!
      
      state_vec_V2.zeros();   // zero out while maintaining dimension
      
      if(n_g0_covars > 0){    // construct detection g0 if there are covariates
        ut_g0 = g0s(0);
        for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
          ut_g0 += g0s(ng0) * y_g0_covar(i, 0, ng0-1);
        }
        g0 = 1 / (1 + exp(-ut_g0));
      }
      
      if(y_platform(i,0) == 0){   //Half-norm
        
        if(n_ds_covars > 0){    // construct detection sigma if there are covariates
          ut_sigmaDet = sigmaDets(0);
          for(int nds = 1; nds <= n_ds_covars; nds++){
            ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, 0, nds-1);
          }
          sigmaDet = exp(ut_sigmaDet);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(north_cells(y_pix(i, 0) - 1)) = state_vec(north_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(west_cells(y_pix(i, 0) - 1)) = state_vec(west_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(south_cells(y_pix(i, 0) - 1)) = state_vec(south_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        state_vec_V2(east_cells(y_pix(i, 0) - 1)) = state_vec(east_cells(y_pix(i, 0) - 1)) * g0 * exp(-pow(detDists(i, 0), 2) / (2*pow(sigmaDet, 2))); 
        
      } else if(y_platform(i, 0) == 1){   //Hazard function
        
        if(n_hs_covars > 0){    // construct detection sigma if there are covariates
          ut_Hsigma = Hsigmas(0);
          for(int nhs = 1; nhs <= n_hs_covars; nhs++){
            ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, 0, nhs-1);
          }
          Hsigma = exp(ut_Hsigma);
        }
        
        if(n_hb_covars > 0){    // construct detection sigma if there are covariates
          ut_Hbeta = Hbetas(0);
          for(int nhb = 1; nhb <= n_hb_covars; nhb++){
            ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, 0, nhb-1);
          }
          Hbeta = exp(ut_Hbeta);
        }
        
        // only update new_state_vec where individual was detected -- probability of being in any other state is 0
        state_vec_V2(north_cells(y_pix(i, 0) - 1)) = state_vec(north_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(west_cells(y_pix(i, 0) - 1)) = state_vec(west_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(south_cells(y_pix(i, 0) - 1)) = state_vec(south_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        state_vec_V2(east_cells(y_pix(i, 0) - 1)) = state_vec(east_cells(y_pix(i, 0) - 1)) * g0 * (1 - exp(-pow(detDists(i, 0) / Hsigma, -Hbeta))); 
        
      }
      
    } else{   // Not detected
      
      //state_vec_V2 = state_vec * p_no_det_occ(allStates,occ);
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,0);
      }
      
    }
    
    // add to log-likelihood
    s_v_sum = sum(state_vec_V2);
    Ind_logProb(i) += log(s_v_sum);
    state_vec = state_vec_V2 / s_v_sum;
    state_vec_V2 = state_vec;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        if(occ == open_transitions(next_open_transition) && next_open_transition < (n_open_transitions - 1)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell)) * gamma(occ) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = state_vec_V2(out1_cell) * (1 - gamma(occ));
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition += 1;
          
          
        } else if(occ == open_transitions(next_open_transition) && next_open_transition == (n_open_transitions - 1)){
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns) * phi(occ) +
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew) * phi(occ) + 
            D_xs.col(cov_open_transitions(next_open_transition)) * 
            ((state_vec_V2(out1_cell) + state_vec_V2(out2_cell) * gamma(occ)) / (4*D_bars(cov_open_transitions(next_open_transition))));
          
          state_vec(out1_cell) = 0;
          
          state_vec(out2_cell) = state_vec_V2(out2_cell) * (1 - gamma(occ)) + 
            (1 - phi(occ)) * (accu(state_vec_V2(north_cells)) + accu(state_vec_V2(west_cells)) + 
            accu(state_vec_V2(south_cells)) + accu(state_vec_V2(east_cells)));
          
          next_open_transition = 0;
          
          
        } else{
          
          state_vec(north_cells) = pMove_North.slice(occ) * (state_vec_V2(north_cells) * p_nsns + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nssn + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(west_cells) = pMove_West.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewew +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewwe);
          
          state_vec(south_cells) = pMove_South.slice(occ) * (state_vec_V2(north_cells) * p_nssn + state_vec_V2(west_cells) * p_ewns +
            state_vec_V2(south_cells) * p_nsns + state_vec_V2(east_cells) * p_ewns);
          
          state_vec(east_cells) = pMove_East.slice(occ) * (state_vec_V2(north_cells) * p_nsew + state_vec_V2(west_cells) * p_ewwe +
            state_vec_V2(south_cells) * p_nsew + state_vec_V2(east_cells) * p_ewew);
          
          state_vec(out1_cell) = state_vec_V2(out1_cell);
          
          state_vec(out2_cell) = state_vec_V2(out2_cell);
          
        }
        
        state_vec_V2 = state_vec;
        occ += 1;
        //Rcpp::Rcout << "k = " << kt << "; occ = " << occ << "; next_open_transition = " << next_open_transition << "\n";
      }
      
      if(y(i,kt) == 1){        // DETECTED!
        
        state_vec_V2.zeros();     // zero out while maintaining dimension,
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kt, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kt) == 0){   //Half-norm
          
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kt, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(north_cells(y_pix(i, kt) - 1)) = state_vec(north_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(west_cells(y_pix(i, kt) - 1)) = state_vec(west_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(south_cells(y_pix(i, kt) - 1)) = state_vec(south_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          state_vec_V2(east_cells(y_pix(i, kt) - 1)) = state_vec(east_cells(y_pix(i, kt) - 1)) * g0 * exp(-pow(detDists(i, kt), 2) / (2*pow(sigmaDet, 2))); 
          
        } else if(y_platform(i, kt) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kt, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
          
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kt, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Rcpp::Rcout << "Hbeta = " << Hbeta << "; Hsigma = " << Hsigma << "\n";
          //Rcpp::Rcout << "Det mod: " << g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))) << "\n";
          //Rcpp::Rcout << "Movement mod: " << state_vec(y_pix(i, kt) - 1) << "\n";
          
          // only update new_state_vec where individual was detected -- probability of being in any other state is 0
          state_vec_V2(north_cells(y_pix(i, kt) - 1)) = state_vec(north_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(west_cells(y_pix(i, kt) - 1)) = state_vec(west_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(south_cells(y_pix(i, kt) - 1)) = state_vec(south_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          state_vec_V2(east_cells(y_pix(i, kt) - 1)) = state_vec(east_cells(y_pix(i, kt) - 1)) * g0 * (1 - exp(-pow(detDists(i, kt) / Hsigma, -Hbeta))); 
          
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        for(int st = 0; st < Nstates; st++){
          state_vec_V2(st) = state_vec(st) * p_no_det_occ(st,kt);
        }
        
      }
      
      s_v_sum = sum(state_vec_V2);
      Ind_logProb(i) += log(s_v_sum);
      state_vec = state_vec_V2 / s_v_sum;
      state_vec_V2 = state_vec;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
      
    }
    
    for(int st = 0; st < Nstates; st++){
      State_Mat(st,i) = state_vec(st);
    }
    
    //Rcpp::Rcout << "next_open_transition: " << next_open_transition << "\n";
    
  }
  
  return State_Mat;
  
}
