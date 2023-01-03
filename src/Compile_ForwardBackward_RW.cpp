
#include <RcppArmadillo.h>
#include <RcppNumerical.h>
using namespace Numer;


// [[Rcpp::export]]
Rcpp::List C_ForwardBackward_AC(int N, int K, int Nstates, 
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
  arma::cube Forward(Nstates,K,N);
  arma::cube Backward(Nstates,K,N);
  arma::mat sv_sum_f(K,N);
  arma::mat sv_sum_b(K,N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    s_v_sum = sum(state_vec);
    //Ind_logProb(i) = log(s_v_sum);    
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
              p_move(y_pix(i,kt), st) * 
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
      
      Forward.slice(i).col(kt) = state_vec;
      sv_sum_f(kt,i) = s_v_sum;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
  }
  
  for(int i = 0; i < N; i++){
    
    int kt = K - 1;
    
    for(int st = 0; st < Nstates; st++){
      state_vec(st) = 1;
      Backward(st,kt,i) = 1;
    }
    
    for(int kt = K-1; kt > 0; kt--){
      
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
              p_move(y_pix(i,kt), st) * 
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
      state_vec = state_vec_V2 / s_v_sum;
      
      Backward.slice(i).col(kt-1) = state_vec;
      sv_sum_b(kt-1,i) = s_v_sum;
      
    }
    
    
  }
  
  //return sum(Ind_logProb);
  Rcpp::List output = Rcpp::List::create(Ind_logProb, Forward, Backward, sv_sum_f, sv_sum_b);
  return output;
  
}


// [[Rcpp::export]]
Rcpp::List C_ForwardBackward_RW(int N, int K, int Nstates, arma::vec tr_b4_occ,
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
  arma::cube Forward(Nstates,K,N);
  arma::cube Backward(Nstates,K,N);
  arma::cube Back_temp_cube(Nstates, Nstates, K);
  arma::cube Back_temp_cube2(Nstates, Nstates, K);
  arma::mat Back_temp_mat(Nstates, Nstates);
  arma::mat ColOnes(Nstates,1);
  
  ColOnes.ones();
  
  
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
    
    Forward.slice(i).col(0) = state_vec;
    //sv_sum_f(0,i) = s_v_sum;
    
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
      
      Forward.slice(i).col(kt) = state_vec;
      //sv_sum_f(kt,i) = s_v_sum;
      
    }
    
    for(int kb = 1; kb < K; kb++){
      arma::mat Det_mat(Nstates, Nstates);
      Det_mat.zeros();
      if(y(i,kb) == 1){        // DETECTED!
        
        if(n_g0_covars > 0){    // construct detection g0 if there are covariates
          ut_g0 = g0s(0);
          for(int ng0 = 1; ng0 <= n_g0_covars; ng0++){
            ut_g0 += g0s(ng0) * y_g0_covar(i, kb, ng0-1);
          }
          g0 = 1 / (1 + exp(-ut_g0));
        }
        
        if(y_platform(i,kb) == 0){   //Half-norm
        
          if(n_ds_covars > 0){    // construct detection sigma if there are covariates
            ut_sigmaDet = sigmaDets(0);
            for(int nds = 1; nds <= n_ds_covars; nds++){
              ut_sigmaDet += sigmaDets(nds) * y_ds_covar(i, kb, nds-1);
            }
            sigmaDet = exp(ut_sigmaDet);
          }
          
          //Det_mat(y_pix(i, kb) - 1, y_pix(i, kb) - 1) = g0 * exp(-pow(detDists(i, kb), 2) / (2*pow(sigmaDet, 2))); 
          
          Back_temp_cube.slice(kb-1).zeros();
          for(int st = 0; st < Nstates; st++){
            Back_temp_cube.slice(kb-1).col(y_pix(i,kb)-1).row(st) = t_mat.slice(kb-1).row(y_pix(i,kb)-1).col(st) * g0 * exp(-pow(detDists(i, kb), 2) / (2*pow(sigmaDet, 2)));
          }
          //Back_temp_cube.slice(kb-1).col(y_pix(i,kb)-1) = (t_mat.slice(kb-1).row(y_pix(i,kb)-1).t() * g0 * exp(-pow(detDists(i, kb), 2) / (2*pow(sigmaDet, 2)))).t(); //???
            //.row instead of .t().col(), because t_mat is column stochastic, but this formulation of backward algorithm requires row stochastic transition matrix
            
        } else if(y_platform(i, kb) == 1){   //Hazard function
          
          if(n_hs_covars > 0){    // construct detection sigma if there are covariates
            ut_Hsigma = Hsigmas(0);
            for(int nhs = 1; nhs <= n_hs_covars; nhs++){
              ut_Hsigma += Hsigmas(nhs) * y_hs_covar(i, kb, nhs-1);
            }
            Hsigma = exp(ut_Hsigma);
          }
        
          if(n_hb_covars > 0){    // construct detection sigma if there are covariates
            ut_Hbeta = Hbetas(0);
            for(int nhb = 1; nhb <= n_hb_covars; nhb++){
              ut_Hbeta += Hbetas(nhb) * y_hb_covar(i, kb, nhb-1);
            }
            Hbeta = exp(ut_Hbeta);
          }
          
          //Det_mat(y_pix(i, kb) - 1, y_pix(i, kb) - 1) = g0 * (1 - exp(-pow(detDists(i, kb) / Hsigma, -Hbeta))); 
          
          Back_temp_cube.slice(kb-1).zeros();
          for(int st = 0; st < Nstates; st++){
            Back_temp_cube.slice(kb-1).col(y_pix(i,kb)-1).row(st) = t_mat.slice(kb-1).row(y_pix(i,kb)-1).col(st) * g0 * (1 - exp(-pow(detDists(i, kb) / Hsigma, -Hbeta)));
          }
          //Back_temp_cube.slice(kb-1).col(y_pix(i,kb)-1) = (t_mat.slice(kb-1).row(y_pix(i,kb)-1) * g0 * (1 - exp(-pow(detDists(i, kb) / Hsigma, -Hbeta)))).t(); //???
            //.row instead of .t().col(), because t_mat is column stochastic, but this formulation of backward algorithm requires row stochastic transition matrix
        }  
        
      } else{   // Not detected
        
        //state_vec_V2 = state_vec % p_no_det_occ(allStates,occ);
        //for(int st = 0; st < Nstates; st++){
        //  Det_mat(st, st) = p_no_det_occ(st,kb);
        //}
        
        for(int st = 0; st < Nstates; st++){
          for(int st2 = 0; st2 < Nstates; st2++){
            Back_temp_cube.slice(kb-1).row(st).col(st2) = t_mat.slice(kb-1).col(st).row(st2) * p_no_det_occ.col(kb).row(st);
          }
          //Back_temp_cube.slice(kb-1).row(st) = (t_mat.slice(kb-1).col(st) % p_no_det_occ.col(kb)).t(); //???
          //.col(st) instead of .t().row(st), because t_mat is column stochastic, but this formulation of backward algorithm requires row stochastic transition matrix
        }
        
        
      }
      
      //Back_temp_cube.slice(kb-1) = t_mat.slice(kb-1) * Det_mat;
      
      //if(kb > 1){
      //  Back_temp_mat = Back_temp_cube.slice(kb-2) * t_mat.slice(kb-1) * Det_mat;
      //  Back_temp_cube.slice(kb-1) = Back_temp_mat / accu(Back_temp_mat);
      //} else{
      //  Back_temp_mat = t_mat.slice(kb-1) * Det_mat;
      //  Back_temp_cube.slice(kb-1) = Back_temp_mat / accu(Back_temp_mat); 
      //}
      
    }
    
    Backward.slice(i).col(K-1) = ColOnes;
    for(int kb = K-2; kb > -1; kb--){
      Backward.slice(i).col(kb) = Back_temp_cube.slice(kb) * Backward.slice(i).col(kb+1);
    }
    
  }
  
  //return sum(Ind_logProb);
  Rcpp::List output = Rcpp::List::create(Ind_logProb, Forward, Backward);//, Back_temp_cube, Back_temp_cube2, ColOnes);
  return output;
  
}


// [[Rcpp::export]]
arma::cube C_ForwardBackward_RWAC(int N, int K, int nCells, int nACCells, 
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
  
  arma::cube Forward_AC(nACCells,K,N);
  arma::cube Forward_RW(nCells,K,N);
  arma::cube Backward_AC(nACCells,K,N);
  arma::cube Backward_RW(nCells,K,N);
  //arma::cube Back_temp_cube(Nstates, Nstates, K);
  //arma::cube Back_temp_cube2(Nstates, Nstates, K);
  //arma::mat Back_temp_mat(Nstates, Nstates);
  //arma::mat ColOnes(Nstates,1);
  
  //ColOnes.ones();
  
  arma::cube Forward(nCells, nACCells, N);
  
  
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
      
      if(kt == K-1){
        Forward.slice(i) = state_mat;
      }
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return Forward;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_mat, state_mat_V2, s_v_sum);
  //return output;
  
}



// [[Rcpp::export]]
arma::cube C_ForwardBackward_CRW(int N, int K, int nCells, 
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
  
  arma::cube Forward(nCells, nCells, N);
  
  
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
      
      if(kt == K-1){
        Forward.slice(i) = state_mat;
      }
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return Forward;
  
}




// [[Rcpp::export]]
arma::mat C_ForecastObs_AC(int N, int K, int Nstates, 
                           arma::vec mu_init, arma::mat y, arma::mat y_pix, arma::mat detDists, arma::mat y_platform,
                           arma::cube y_ds_covar, arma::cube y_g0_covar, arma::cube y_hs_covar, arma::cube y_hb_covar,
                           int n_ds_covars, int n_g0_covars, int n_hs_covars, int n_hb_covars,
                           arma::vec g0s, arma::vec sigmaDets, arma::vec Hsigmas, arma::vec Hbetas,
                           arma::mat p_no_det_occ, arma::mat p_move,
                           arma::mat FBProbs){
  
  arma::vec Ind_logProb(N);
  double sigmaDet = exp(sigmaDets(0));
  double g0 = 1 / (1 + exp(-g0s(0)));
  double Hsigma = exp(Hsigmas(0));
  double Hbeta = exp(Hbetas(0));
  double ut_sigmaDet;
  double ut_g0;
  double ut_Hsigma;
  double ut_Hbeta;
  arma::vec state_vec(Nstates);
  arma::mat s_v_sum(N, K);
  arma::vec state_vec_V2(Nstates);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = FBProbs.col(i);    // Don't update state_vec; calculate forecast conditional on FB probs, do not continue forward algorithm
    state_vec_V2 = FBProbs.col(i);   
    
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
      
      s_v_sum(i,kt) = sum(state_vec_V2);
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
  }
  
  return s_v_sum;
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_vec, state_vec_V2, s_v_sum, sigmaDet, g0);
  //return output;
  
}


// [[Rcpp::export]]
arma::mat C_ForecastObs_RW(int N, int K, int Nstates, arma::vec tr_b4_occ,
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
                           arma::mat FBProbs){
  
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
  arma::mat s_v_sum(N,K);
  arma::vec state_vec_V2(Nstates);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = FBProbs.col(i);
    state_vec_V2 = FBProbs.col(i);
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
    s_v_sum(i,0) = sum(state_vec_V2);
    state_vec_V2 = state_vec;
    
    
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
      
      s_v_sum(i,kt) = sum(state_vec_V2);
      state_vec_V2 = state_vec;
      
    }
  }
  
  return s_v_sum;
  
}



// [[Rcpp::export]]
arma::mat C_ForecastObs_RWAC(int N, int K, int nCells, int nACCells, 
                                  arma::vec tr_b4_occ,
                                  arma::vec mu_init, arma::cube FBProbs,
                                  arma::mat y, arma::mat y_pix, 
                                  arma::mat detDists, arma::mat y_platform,
                                  arma::cube y_ds_covar, arma::cube y_g0_covar,
                                  arma::cube y_hs_covar, arma::cube y_hb_covar,
                                  int n_ds_covars, int n_g0_covars,
                                  int n_hs_covars, int n_hb_covars,
                                  arma::vec g0s, arma::vec sigmaDets,
                                  arma::vec Hsigmas, arma::vec Hbetas,
                                  arma::mat p_no_det_occ, arma::cube t_mat){
  
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
  arma::mat s_v_sum(N,K);
  arma::mat state_mat_V2(nCells, nACCells);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = FBProbs.slice(i);
    state_mat_V2 = FBProbs.slice(i);
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
    s_v_sum(i,0) = accu(state_mat_V2);
    state_mat_V2 = state_mat;
    
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
      
      s_v_sum(i,kt) = accu(state_mat_V2);
      state_mat_V2 = state_mat;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return s_v_sum;
  //return sum(Ind_logProb);
  //Rcpp::List output = Rcpp::List::create(Ind_logProb, state_mat, state_mat_V2, s_v_sum);
  //return output;
  
}



// [[Rcpp::export]]
arma::mat C_ForecastObs_CRW(int N, int K, int nCells, 
                                  arma::vec tr_b4_occ,
                                  arma::vec mu_init, arma::cube FBProbs,
                                  arma::mat y, arma::mat y_pix, 
                                  arma::mat detDists, arma::mat y_platform,
                                  arma::cube y_ds_covar, arma::cube y_g0_covar,
                                  arma::cube y_hs_covar, arma::cube y_hb_covar,
                                  int n_ds_covars, int n_g0_covars,
                                  int n_hs_covars, int n_hb_covars,
                                  arma::vec g0s, arma::vec sigmaDets,
                                  arma::vec Hsigmas, arma::vec Hbetas,
                                  arma::mat p_no_det_occ, arma::cube t_mat){
  
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
  arma::mat s_v_sum(N,K);
  arma::mat state_mat_V2(nCells, nCells);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = FBProbs.slice(i);
    state_mat_V2 = FBProbs.slice(i);
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
    s_v_sum(i,0) = accu(state_mat_V2);
    state_mat_V2 = state_mat;
    
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
      
      s_v_sum(i,kt) = accu(state_mat_V2);
      state_mat_V2 = state_mat;
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
    
  }
  
  return s_v_sum;
  
}




// [[Rcpp::export]]
Rcpp::List C_ForecastState_AC(int N, int K, int Nstates, 
                             arma::mat y_pix, arma::mat p_move,
                             arma::mat FBProbs){
  
  arma::vec state_vec(Nstates);
  arma::mat s_v_sum(N, K);
  arma::vec state_vec_V2(Nstates);
  arma::mat Forecast_state(Nstates,N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = FBProbs.col(i);    // Don't update state_vec; calculate forecast conditional on FB probs, do not continue forward algorithm
    state_vec_V2 = FBProbs.col(i);   
    
    for(int st = 0; st < Nstates; st++){        //activity center cell
      for(int st2 = 0; st2 < Nstates; st2++){   //actual location cell
        Forecast_state(st2,i) += state_vec(st) * p_move(st2, st);
      }
    }
     
    
    for(int kt = 0; kt < K; kt++){      // loop through remaining sampling occasions
      
      state_vec_V2.zeros();     // zero out while maintaining dimension,
      
      for(int st = 0; st < Nstates; st++){
        state_vec_V2(st) = p_move(y_pix(i,kt) - 1, st) * state_vec(st);
      }
      
      s_v_sum(i,kt) = sum(state_vec_V2);
      
      //if(Rcpp::traits::is_nan<REALSXP>(s_v_sum)){
      //  Rcpp::Rcout << "s_v_sum is NaN at " << kt+1 << "\n";
      //}
      
      //Rcpp::Rcout << "kt = " << kt + 1 << " Ind_logProb = " << Ind_logProb(i) << "\n";
      
    }
  }
  
  //return s_v_sum;
  Rcpp::List output = Rcpp::List::create(s_v_sum, Forecast_state);
  return output;
  
}



// [[Rcpp::export]]
Rcpp::List C_ForecastState_RW(int N, int K, int Nstates, 
                             arma::vec tr_b4_occ,
                             arma::mat y_pix, arma::cube t_mat,
                             arma::mat FBProbs){
  
  arma::vec state_vec(Nstates);
  arma::mat s_v_sum(N, K);
  arma::vec state_vec_V2(Nstates);
  arma::cube Forecast_state(Nstates,K,N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = FBProbs.col(i);
    state_vec_V2 = FBProbs.col(i);
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
    
    Forecast_state.slice(i).col(0) = state_vec;
    
    s_v_sum(i, 0) = state_vec(y_pix(i,0) - 1);
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        state_vec_V2 = t_mat.slice(occ) * state_vec;
        state_vec = state_vec_V2;
        occ += 1;
      }
      
      Forecast_state.slice(i).col(kt) = state_vec;
      
      s_v_sum(i,kt) = state_vec(y_pix(i,kt) - 1);
      
    }
  }
  
  //return s_v_sum;
  Rcpp::List output = Rcpp::List::create(s_v_sum, Forecast_state);
  return output;
  
}


// [[Rcpp::export]]
Rcpp::List C_ForecastState_RWAC(int N, int K, int nCells, int nACCells, 
                             arma::vec tr_b4_occ,
                             arma::vec mu_init, arma::cube FBProbs,
                             arma::mat y_pix, arma::cube t_mat){
  
  arma::vec state_vec = mu_init;
  arma::mat state_mat(nCells,nACCells);
  arma::mat s_v_sum(N,K);
  arma::mat state_mat_V2(nCells, nACCells);
  arma::cube Forecast_state_RW(nCells,K,N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = FBProbs.slice(i);
    state_mat_V2 = FBProbs.slice(i);
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
    
    //Forecast of actual location - marginalized over activity centers
    Forecast_state_RW.slice(i).col(0) = sum(state_mat_V2, 1);
    
    // add to log-likelihood
    s_v_sum(i,0) = accu(state_mat_V2.row(y_pix(i,0) - 1));
    //Ind_logProb(i) += log(s_v_sum);
    state_mat_V2 = state_mat;
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int ac = 0; ac < nACCells; ac++){
          state_mat_V2.col(ac) = t_mat.slice(ac) * state_mat.col(ac);
        }
        state_mat = state_mat_V2;
        occ += 1;
      }
      
      //Forecast of actual location - marginalized over activity centers
      Forecast_state_RW.slice(i).col(kt) = sum(state_mat_V2, 1);
      
      s_v_sum(i,kt) = accu(state_mat_V2.row(y_pix(i,kt) - 1));;
      //Ind_logProb(i) += log(s_v_sum);
      state_mat_V2 = state_mat;
      
    }
    
  }
  
  Rcpp::List output = Rcpp::List::create(s_v_sum, Forecast_state_RW);
  return output;
  
}


// [[Rcpp::export]]
Rcpp::List C_ForecastState_CRW(int N, int K, int nCells, 
                            arma::vec tr_b4_occ,
                            arma::vec mu_init, arma::cube FBProbs,
                            arma::mat y_pix, arma::cube t_mat){
  

  arma::vec state_vec = mu_init;
  arma::mat state_mat(nCells,nCells);
  arma::mat s_v_sum(N,K);
  arma::mat state_mat_V2(nCells, nCells);
  arma::cube Forecast_state_RW(nCells,K,N);
  
  
  for(int i = 0; i < N; i++){
    
    state_vec = mu_init;
    state_mat = FBProbs.slice(i);
    state_mat_V2 = FBProbs.slice(i);
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

    //Forecast of actual location - marginalized over activity centers
    Forecast_state_RW.slice(i).col(0) = sum(state_mat_V2, 1);
    
    // add to log-likelihood
    s_v_sum(i,0) = accu(state_mat_V2.row(y_pix(i,0) - 1));
    //Ind_logProb(i) += log(s_v_sum);
    state_mat_V2 = state_mat;
    
    
    for(int kt = 1; kt < K; kt++){      // loop through remaining sampling occasions
      
      for(int tr = 0; tr < tr_b4_occ(kt); tr++){
        for(int k = 0; k < nCells; k++){
          state_mat_V2.col(k) = t_mat.slice(k) * state_mat.row(k).t();
        }
        state_mat = state_mat_V2;
        occ += 1;
      }
      
      //Forecast of actual location - marginalized over activity centers
      Forecast_state_RW.slice(i).col(kt) = sum(state_mat_V2, 1);
      
      // add to log-likelihood
      s_v_sum(i,kt) = accu(state_mat_V2.row(y_pix(i,0) - 1));
      //Ind_logProb(i) += log(s_v_sum);
      state_mat_V2 = state_mat;
      
    }
    
  }
  
  Rcpp::List output = Rcpp::List::create(s_v_sum, Forecast_state_RW);
  return output;
  
}



// [[Rcpp::export]]
Rcpp::List C_pDots_RWAC(int K, Rcpp::NumericVector tr_b4_occ,
                         int nCells, int nACCells, arma::mat p_init_state_full,
                         arma::cube pMove, arma::mat p_no_det_occ){
  arma::cube p_no_det_tot(nCells, nACCells, K);
  arma::cube state_distribution(nCells, nACCells, K);
  arma::mat temp_1(nCells, nACCells);
  arma::mat temp_2(nCells, nACCells);
  arma::vec pDot(K);
  
  if(tr_b4_occ(0) > 0){
    for(int i = 0; i < nACCells; i++){
      temp_1.col(i) = pMove.slice(i) * p_init_state_full.col(i);
      //p_no_det_tot.slice(0).col(i) = pMove.slice(i) * p_init_state_full.col(i) % p_no_det_occ.col(0); 
    }
  } else{
    temp_1 = p_init_state_full;
    //for(int i = 0; i < nACCells; i++){
    //  p_no_det_tot.slice(0).col(i) = p_init_state_full.col(i) % p_no_det_occ.col(0);
    //}
  }
  
  for(int i = 0; i < nACCells; i++){
    p_no_det_tot.slice(0).col(i) = temp_1.col(i) % p_no_det_occ.col(0);
  }
  
  state_distribution.slice(0) = temp_1;
  pDot(0) = 1 - accu(p_no_det_tot.slice(0));
  
  for(int k = 1; k < K; k++){
    if(tr_b4_occ(k) > 1){
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
        temp_2.col(i) = pMove.slice(i) * temp_1.col(i);
        //p_no_det_tot.slice(k).col(i) = pMove.slice(i) * p_no_det_tot.slice(k-1).col(i) % p_no_det_occ.col(k);
      }
      temp_1 = temp_2;
      for(int i = 0; i < nACCells; i++){
        p_no_det_tot.slice(k).col(i) = temp_1.col(i) % p_no_det_occ.col(k);
      }
    }
    
    pDot(k) = 1 - accu(p_no_det_tot.slice(k));
    state_distribution.slice(k) = temp_1;
  }
  
  Rcpp::List output = Rcpp::List::create(pDot, p_no_det_tot, state_distribution);
  return output;
  
}



// [[Rcpp::export]]
Rcpp::List C_pDots_CRW(int K, Rcpp::NumericVector tr_b4_occ,
                                 int nCells, arma::mat p_init_state_full,
                                 arma::cube pMove, arma::mat p_no_det_occ){
  arma::cube p_no_det_tot(nCells, nCells, K);
  arma::mat temp_1(nCells, nCells);
  arma::mat temp_2(nCells, nCells);
  arma::vec pDot(K);
  
  if(tr_b4_occ(0) > 0){
    for(int i = 0; i < nCells; i++){
      temp_1.col(i) = pMove.slice(i) * p_init_state_full.row(i).t(); 
    }
  } else{
    for(int i = 0; i < nCells; i++){
      temp_1.col(i) = p_init_state_full.col(i);
    }
  }
  
  for(int i = 0; i < nCells; i++){
    p_no_det_tot.slice(0).col(i) = temp_1.col(i) % p_no_det_occ.col(0);
  }
  
  pDot(0) = 1 - accu(p_no_det_tot.slice(0));
  
    
  for(int k = 1; k < K; k++){
    if(tr_b4_occ(k) > 1){
      for(int j = 0; j < tr_b4_occ(k); j++){
        for(int i = 0; i < nCells; i++){
          temp_2.col(i) = pMove.slice(k-1) * temp_1.row(i).t();
        }
        temp_1 = temp_2;
      }
      for(int i = 0; i < nCells; i++){
        p_no_det_tot.slice(k).col(i) = temp_1.col(i) % p_no_det_occ.col(k);
      }
    } else{
      for(int i = 0; i < nCells; i++){
        temp_2.col(i) = pMove.slice(i) * temp_1.row(i).t();
      }
      temp_1 = temp_2;
      for(int i = 0; i < nCells; i++){
        p_no_det_tot.slice(k).col(i) = temp_1.col(i) % p_no_det_occ.col(k);
      }
    }
    pDot(k) = 1 - accu(p_no_det_tot.slice(k));
  }
  

  Rcpp::List output = Rcpp::List::create(pDot, p_no_det_tot);
  return output;
  
}

