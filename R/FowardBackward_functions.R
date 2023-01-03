#####################################################################################
### Forward and Backward Probabilities and Forecasting for Random Walk SCR Models ###
#####################################################################################

################################################################################
################################################################################
### Calculate foward-backward probabilities for closed activity center model ###
################################################################################
################################################################################
#' Calculate the forward-backward probabilities for the location of each individual on each occasion under a spatial capture-recapture model with an activity center movement model
#' 
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "AC"
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "AC"
#'
#' @export
ForwardBackward_Closed_AC = function(data, par){
  y = data$y
  detDists = data$detDists 
  habitat = data$habitat
  detGrid = data$detGrid
  lines_arr = data$lines.arr
  distMat = data$distMat
  y.pix = data$y.pix
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  nCells = data$nCells
  nDetCells = data$nDetCells
  N = data$N
  K = data$K
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  
  y.g0.covar = data$y.g0.covar
  y.ds.covar = data$y.ds.covar
  y.hs.covar = data$y.hs.covar
  y.hb.covar = data$y.hb.covar
  y.platform = data$y.platform
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  #gammaMod = data$gammaMod
  #phiMod = data$phiMod
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  
  if(n.d.covars == 0){
    p.move = C_pMove_Laplace_Full(K = 1, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars)
  }else{
    p.move = C_pMove_Laplace_Full(K = 1, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars)
  }
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  p.no.det.list = C_pNoDet(nCells = nCells, K = K, pDet = p.det.arr, pMove = p.move[,,1])
  p.no.det.given.ac = p.no.det.list[[1]]
  p.no.det.occ = p.no.det.list[[2]]
  
  #probability of being detected at least once
  p.dot = 1 - sum(p.init.state * p.no.det.given.ac)
  
  if(is.na(p.dot) | p.dot < 0){
    #print(paste0("p.dot = ", p.dot, " < 0; par = "))
    #print(par)
    p.dot = 0
  }
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  # expected number of unobserved individuals
  EN.unobs = sum(D.x * p.no.det.given.ac)
  
  
  ### Calculate sum(log-likelihood of individuals encounter histories) ###
  FB_Probs = C_ForwardBackward_AC(N = N+1, K = K, Nstates = N.states,  
                                  mu_init = mu.init, 
                                  y = rbind(y, rep(0,K)),
                                  y_pix = rbind(y.pix, rep(0,K)),
                                  detDists = rbind(detDists, rep(0,K)),
                                  y_platform = rbind(y.platform, rep(NA,K)),
                                  y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                  n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                  g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                  p_no_det_occ = p.no.det.occ, p_move = p.move[,,1])
  
  FB.arr = FB_Probs[[2]]*FB_Probs[[3]]
  for(i in 1:(N+1)){
    for(k in 1:K){
      FB.arr[,k,i] = FB.arr[,k,i]/sum(FB.arr[,k,i])
    }
  }
  
  FB.mat = FB.arr[,K,]
  
  return(list(ForwardBackward_Full = FB.arr, FB.T = FB.mat))
}


############################################################################
############################################################################
### Calculate foward-backward probabilities for closed random walk model ###
############################################################################
############################################################################
#' Calculate the forward-backward probabilities for the location of each individual on each occasion under a spatial capture-recapture model with a simple random walk movement model
#' 
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW"
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW"
#'
#' @export
ForwardBackward_Closed_SimpleRW  = function(data, par){  
  
  y = data$y
  detDists = data$detDists 
  habitat = data$habitat
  detGrid = data$detGrid
  lines_arr = data$lines.arr
  distMat = data$distMat
  y.pix = data$y.pix
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  nCells = data$nCells
  nDetCells = data$nDetCells
  N = data$N
  K = data$K
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K - 1 #data$N.transitions # number of transitions, number of primary periods - 1
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ
  
  y.g0.covar = data$y.g0.covar
  y.ds.covar = data$y.ds.covar
  y.hs.covar = data$y.hs.covar
  y.hb.covar = data$y.hb.covar
  y.platform = data$y.platform
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  #gammaMod = data$gammaMod
  #phiMod = data$phiMod
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  
  if(n.d.covars == 0){
    p.move = C_pMove_Laplace_Full(K = K, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars)
  }else{
    p.move = C_pMove_Laplace_Full(K = K, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  # transition matrices between periods #
  t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  for(i in 2:K){  
    # use 2:N.transitions
    # habitat at t = 1 determines initial distribution
    # then habitat at t = 2:N.transitions determines movement to those locations in those sampling periods
    
    t.mat[,,i-1] = p.move[,,i]
    
  }
  
  
  ### Could make p.no.det.tot a vector and update it each iteration if that is faster? ###
  p.no.det.tot = array(dim = c(N.states, K, 2))
  
  # if the study begins with a first transition, the dimension of t.mat will need to be changed; it isn't that flexible yet #
  if(tr_b4_occ[1] > 0){
    p.no.det.tot[,1,1] = t.mat[,,1] %*% p.init.state * p.no.det.occ[,1]
    p.no.det.tot[,1,2] = t.mat[,,1] %*% D.x * p.no.det.occ[,1]
  } else{
    p.no.det.tot[,1,1] = p.init.state * p.no.det.occ[,1]
    p.no.det.tot[,1,2] = D.x * p.no.det.occ[,1]
  }
  
  for(i in 2:K){
    if(tr_b4_occ[i] > 1){
      temp = p.no.det.tot[,i-1,1]
      temp.d = p.no.det.tot[,i-1,2]
      for(j in 1:(tr_b4_occ[i]-1)){
        temp = t.mat[,,i-1] %*% temp   # this assumes same transition for each day without observation
        temp.d = t.mat[,,i-1] %*% temp.d
        ### MAY WANT TO BUILD IN FLEXIBILITY INTO t.mat TO ACCOMODATE DIFFERENT TRANSITIONS HERE ###
      }
      p.no.det.tot[,i,1] = temp * p.no.det.occ[,i]
      p.no.det.tot[,i,2] = temp.d * p.no.det.occ[,i]
      
    } else{
      p.no.det.tot[,i,1] = t.mat[,,i-1] %*% p.no.det.tot[,i-1,1] * p.no.det.occ[,i]
      p.no.det.tot[,i,2] = t.mat[,,i-1] %*% p.no.det.tot[,i-1,2] * p.no.det.occ[,i]
    }
  }
  
  #probability of being detected at least once
  p.dot = 1 - sum(p.no.det.tot[,K,1])
  
  EN.unobs = sum(p.no.det.tot[,K,2])
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  ### Set this up following book; but t.mat is row stochastic in book, whereas mine is column stochastic, so input the transpose  ###
  ### the transpose is row stochastic, but should still double check the indexing is correct t(p.init.state) %*% t(t.mat) =? t.mat %*% p.init.state
  ### transposition works but t.mat is array, so transpose slices in c++ function
    ### ONLY FOR BACKWARDS PROBS - can use existing code to get forward probs easily ###
  t1 = Sys.time()
  FB_Probs = C_ForwardBackward_RW(N = N + 1, K = K, Nstates = N.states, tr_b4_occ = tr_b4_occ,
                                  mu_init = mu.init, 
                                  y = rbind(y, rep(0,K)), y_pix = rbind(y.pix, rep(0,K)), 
                                  detDists = rbind(detDists, rep(0,K)), y_platform = rbind(y.platform, rep(0,K)),
                                  y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar,
                                  y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                  n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars,
                                  n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                  g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                  p_no_det_occ = p.no.det.occ, t_mat = t.mat)
  t2 = Sys.time()
  t2 - t1
  
  
  FB.arr = FB_Probs[[2]]*FB_Probs[[3]]
  for(i in 1:(N+1)){
    for(k in 1:K){
      FB.arr[,k,i] = FB.arr[,k,i]/sum(FB.arr[,k,i])
    }
  }
  
  
  return(list(FBProbs = FB.arr, EN.unobs = EN.unobs))
}




##################################################################################################
##################################################################################################
### Calculate foward-backward probabilities for closed random walk with activity center model ###
##################################################################################################
##################################################################################################
#' Calculate the forward-backward probabilities for the location of each individual on each occasion under a spatial capture-recapture model with a biased random walk movement model
#' 
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_BiasedRW 
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_BiasedRW 
#'
#' @export
ForwardBackward_Closed_BiasedRW = function(data, par){  
  
  y = data$y
  detDists = data$detDists
  detGrid = data$detGrid
  lines_arr = data$lines.arr
  y.pix = data$y.pix
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  
  nCells = data$nCells
  nDetCells = data$nDetCells
  nACCells = data$nACCells
  
  habitat = data$habitat
  habitat.AC = data$habitat.AC
  
  distMat = data$distMat
  distMat.AC = data$distMat.AC
  
  N = data$N
  K = data$K
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K - 1 #data$N.transitions # number of transitions, number of primary periods - 1
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ
  
  y.g0.covar = data$y.g0.covar
  y.ds.covar = data$y.ds.covar
  y.hs.covar = data$y.hs.covar
  y.hb.covar = data$y.hb.covar
  y.platform = data$y.platform
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  rho.init.indices = data$rho.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  rho = 1 / (1 + exp(-par[rho.init.indices]))
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  # have to integrate across all possible activity centers
  # probability of moving from cell i to cell j given the activity center is in cell k
  # u_t ~ N(u_t-1 * p + s * (1 - p), sigma^2I), where 0 <= p <= 1
  # need to loop over these calcs for each N.states s -> do this in C++
  if(n.d.covars == 0){
    p.move = C_pMove_RWAC_Laplace_Full(nCells = nCells, nACCells = nACCells, habitat = habitat,
                                       betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars, 
                                       rho = rho)
  }else{
    p.move = C_pMove_RWAC_Laplace_Full(nCells = nCells, nACCells = nACCells, habitat = habitat, 
                                       betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                       rho = rho)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat.AC[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat.AC[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  
  # transition matrices between periods #
  # t.mat should be 4d, cell to, cell from, activity center cell, transitions
  # likelihood and calculation of thinning rate (p.dot) need to be integrated over all possible activity centers
  #t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  
  ### THE ABOVE (mu.init) IS THE POINT PROCESS FOR ACTIVITY CENTERS ###
  ### ALSO NEED TO DESCRIBE INITIAL LOCATION ABOUT ACTIVITY CENTERS ###
  ### USE A MATRIX INSTEAD OF A VECTOR TO DESCRIBE ACTIVITY CENTER AND LOCATION STATE ###
  ### LAPLACE INITIAL v2 -> changed from marginalizing over ac's to marginilizing over all previous locations and ac's ###
  ### v2 ALSO DID NOT WORK ###
  ### i.e., treat the situation like the individual's location was uniformly distributed during the time step before the first occasion ###
  ### having issues estimating rho, think it is because there is too much influence from first occasion on activity center location ###
  ### could also use a model where the first location is distributed uniformly across the state space, regardless of the activity center location ###
  if(n.d.covars == 0){
    p.init.state.full = C_pMove_RWAC_Laplace_Initial(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
                                                     betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                     p_init_state_ac = p.init.state)
    #p.init.state.full = C_pMove_RWAC_Laplace_Initial_v2(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
    #                                                 betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
    #                                                 p_init_state_ac = p.init.state, rho = rho)
    
    # try uniformly distributed initial location and activity center
    #p.init.state.full = matrix(rep(p.init.state/nCells, nCells), nrow = nCells, ncol = nACCells, byrow = T)
    
  } else{
    p.init.state.full = C_pMove_RWAC_Laplace_Initial(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
                                                     betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                     p_init_state_ac = p.init.state)
    #p.init.state.full = C_pMove_RWAC_Laplace_Initial_v2(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
    #                                                 betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
    #                                                 p_init_state_ac = p.init.state, rho = rho)
    
    # try uniformly distributed initial location
    #p.init.state.full = matrix(rep(p.init.state/nCells, nCells), nrow = nCells, ncol = nACCells, byrow = T)
  }
  
  # p.init.state.full: activity center is column index; movement cell is row index 
  
  
  
  # probability of having not yet been detected on each occasion given initial activity center and locations  
  pNoDet_RWAC_out = C_pNoDet_RWAC(K = K, tr_b4_occ = tr_b4_occ,
                                  nCells = nCells, nACCells = nACCells, p_init_state_full = p.init.state.full,
                                  pMove = p.move, p_no_det_occ = p.no.det.occ)
  # probability of being detected at least once 
  p.dot = pNoDet_RWAC_out[[1]]
  p.no.det.tot = pNoDet_RWAC_out[[2]]
  
  if(is.na(p.dot) | p.dot < 0){
    #print(paste0("p.dot = ", p.dot, " < 0; par = "))
    #print(par)
    p.dot = 0
  }
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  EN.unobs = (1 - p.dot) * D.bar
  
  
  ### Calculate sum(log-likelihood of individuals encounter histories) ###
  Forward.K = C_ForwardBackward_RWAC(N = N + 1, K = K, nCells = nCells, nACCells = nACCells,
                                   tr_b4_occ = tr_b4_occ, 
                                   mu_init = mu.init, p_init_state_full = p.init.state.full, 
                                   y = rbind(y, rep(0,K)), y_pix = rbind(y.pix, rep(0,K)), 
                                   detDists = rbind(detDists, rep(0,K)), y_platform = rbind(y.platform, rep(0,K)),
                                   y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                   n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                   g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                   p_no_det_occ = p.no.det.occ, t_mat = p.move)
  
  
  return(list(FBProbs = Forward.K, EN.unobs = EN.unobs))
}



#######################################################################################
#######################################################################################
### Calculate foward-backward probabilities for closed correlated random walk model ###
#######################################################################################
#######################################################################################
#' Calculate the forward-backward probabilities for the location of each individual on each occasion under a spatial capture-recapture model with a correlated random walk movement model
#' 
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_CorrelatedRW 
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_CorrelatedRW 
#'
#' @export
ForwardBackward_Closed_CorrelatedRW = function(data, par){  
  
  y = data$y
  detDists = data$detDists
  detGrid = data$detGrid
  lines_arr = data$lines.arr
  y.pix = data$y.pix
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  
  nCells = data$nCells
  nDetCells = data$nDetCells
  nACCells = data$nACCells
  
  habitat = data$habitat
  habitat.AC = data$habitat.AC
  
  distMat = data$distMat
  distMat.AC = data$distMat.AC
  
  N = data$N
  K = data$K
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K - 1 #data$N.transitions # number of transitions, number of primary periods - 1
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ
  
  y.g0.covar = data$y.g0.covar
  y.ds.covar = data$y.ds.covar
  y.hs.covar = data$y.hs.covar
  y.hb.covar = data$y.hb.covar
  y.platform = data$y.platform
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  gamma.init.indices = data$gamma.init.indices
  beta.init.indices = data$beta.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  gamma = 1 / (1 + exp(-par[gamma.init.indices]))
  beta = (2*pi / (1 + exp(-par[beta.init.indices]))) - pi
  
  rotation = gamma * matrix(c(cos(beta), sin(beta), -sin(beta), cos(beta)), nrow = 2)
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  # have to integrate across all possible activity centers
  # probability of moving from cell i to cell j given the activity center is in cell k
  # u_t ~ N(u_t-1 * p + s * (1 - p), sigma^2I), where 0 <= p <= 1
  # need to loop over these calcs for each N.states s -> do this in C++
  if(n.d.covars == 0){
    p.move = C_pMove_CorrelatedRW_Laplace_Full(nCells = nCells, habitat = habitat,
                                               betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars, 
                                               rotation = rotation)
  }else{
    p.move = C_pMove_CorrelatedRW_Laplace_Full(nCells = nCells, habitat = habitat, 
                                               betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                               rotation = rotation)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat.AC[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat.AC[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  
  # transition matrices between periods #
  # t.mat should be 4d, cell to, cell from, activity center cell, transitions
  # likelihood and calculation of thinning rate (p.dot) need to be integrated over all possible activity centers
  #t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  
  ### THE ABOVE (mu.init) IS THE POINT PROCESS FOR ACTIVITY CENTERS ###
  ### ALSO NEED TO DESCRIBE INITIAL LOCATION ABOUT ACTIVITY CENTERS ###
  ### USE A MATRIX INSTEAD OF A VECTOR TO DESCRIBE ACTIVITY CENTER AND LOCATION STATE ###
  if(n.d.covars == 0){
    p.init.state.full = C_pMove_CorrelatedRW_Laplace_Initial(nCells = nCells, habitat = habitat, distMat = distMat,
                                                             betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                             p_init_state = p.init.state)
    
    
  } else{
    p.init.state.full = C_pMove_CorrelatedRW_Laplace_Initial(nCells = nCells, habitat = habitat, distMat = distMat,
                                                             betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                             p_init_state = p.init.state)
    
  }
  
  # p.init.state.full: activity center is column index; movement cell is row index 
  
  
  
  # probability of having not yet been detected on each occasion given initial activity center and locations  
  pNoDet_CorrelatedRW_out = C_pNoDet_CorrelatedRW(K = K, tr_b4_occ = tr_b4_occ,
                                                  nCells = nCells, p_init_state_full = p.init.state.full,
                                                  pMove = p.move, p_no_det_occ = p.no.det.occ)
  # probability of being detected at least once 
  p.dot = pNoDet_CorrelatedRW_out[[1]]
  p.no.det.tot = pNoDet_CorrelatedRW_out[[2]]

  
  if(is.na(p.dot) | p.dot < 0 | p.dot > 1){
    #print(paste0("p.dot = ", p.dot, " < 0; par = "))
    #print(par)
    p.dot = 0
  }
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  EN.unobs = (1 - p.dot) * D.bar
  
  
  ### Calculate sum(log-likelihood of individuals encounter histories) ###
  Forward.K = C_ForwardBackward_CRW(N = N + 1, K = K, nCells = nCells, 
                                           tr_b4_occ = tr_b4_occ, 
                                           mu_init = mu.init, p_init_state_full = p.init.state.full, 
                                           y = rbind(y, rep(0,K)), 
                                    y_pix = rbind(y.pix, rep(0,K)), 
                                    detDists = rbind(detDists, rep(0,K)), 
                                    y_platform = rbind(y.platform, rep(0,K)),
                                           y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                           n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                           g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                           p_no_det_occ = p.no.det.occ, t_mat = p.move)
  
  
  
  
  return(list(FBProbs = Forward.K, EN.unobs = EN.unobs))
}



#############################################################################################################
#############################################################################################################
### Calculate likelihood of future observations and states (locations) given forecast for closed AC model ###
#############################################################################################################
#############################################################################################################
#' Calculate the likelihood of future encounter histories and movement trajectories under the forecast for a spatial capture-recapture model with an activity center movement model
#'
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "AC"
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "AC"
#' @param FBProbs Array of forward-backward probabilities; passed from ForwardBackward_Closed_AC
#'
#' @export
Forecast_Closed_AC = function(data, par, FBProbs){
  y = data$y.forecast
  detDists = data$detDists.forecast 
  habitat = data$habitat
  detGrid = data$detGrid.forecast
  lines_arr = data$lines.arr.forecast
  distMat = data$distMat
  y.pix = data$y.pix.forecast
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  nCells = data$nCells
  nDetCells = data$nDetCells
  N = data$N.forecast
  K = data$K.forecast
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  
  y.g0.covar = data$y.g0.covar.forecast
  y.ds.covar = data$y.ds.covar.forecast
  y.hs.covar = data$y.hs.covar.forecast
  y.hb.covar = data$y.hb.covar.forecast
  y.platform = data$y.platform.forecast
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  #gammaMod = data$gammaMod
  #phiMod = data$phiMod
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  
  if(n.d.covars == 0){
    p.move = C_pMove_Laplace_Full(K = 1, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars)
  }else{
    p.move = C_pMove_Laplace_Full(K = 1, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars)
  }
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  p.no.det.list = C_pNoDet(nCells = nCells, K = K, pDet = p.det.arr, pMove = p.move[,,1])
  p.no.det.given.ac = p.no.det.list[[1]]
  p.no.det.occ = p.no.det.list[[2]]
  
  #probability of being detected at least once
  p.dot = 1 - sum(p.init.state * p.no.det.given.ac)
  
  if(is.na(p.dot) | p.dot < 0){
    #print(paste0("p.dot = ", p.dot, " < 0; par = "))
    #print(par)
    p.dot = 0
  }
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  # pad data with an unobserved encounter history - representative of all individuals still unobserved through forecast sampling
  Forecast_Distribution_prevObs = C_ForecastObs_AC(N = data$N, K = K, Nstates = N.states,  
                                                   mu_init = mu.init, 
                                                   y = y[1:data$N,],
                                                   y_pix = y.pix[1:data$N,],
                                                   detDists = detDists[1:data$N,],
                                                   y_platform = y.platform[1:data$N,],
                                                   y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                   n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                   g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                   p_no_det_occ = p.no.det.occ, p_move = p.move[,,1], FBProbs = FBProbs[,1:data$N])
  
  
  # expected number of unobserved individuals
  EN.unobs = sum(D.x * p.no.det.given.ac)
  D.x.unobs = FBProbs[,data$N+1] * EN.unobs
  mu.x.unobs = D.x.unobs * (1-p.no.det.occ)
  
  if(N - data$N == 0){
    Occ_dist_unobs = rep(-mu.bar, K)
  } else{
    Forecast_Distribution_Unobs_v1 = C_ForecastObs_AC(N = N - data$N, K = K, Nstates = N.states,  
                                                      mu_init = mu.init, 
                                                      y = matrix(y[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_pix = matrix(y.pix[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      detDists = matrix(detDists[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_platform = matrix(y.platform[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                      n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                      g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                      p_no_det_occ = p.no.det.occ, p_move = p.move[,,1], 
                                                      FBProbs = matrix(rep(D.x.unobs, N - data$N), nrow = nrow(FBProbs), ncol = N - data$N))
    # above FBProbs is D.x.unobs following the Poisson point process likelihood listed in Glennie et al.
    
    Occ_dist_unobs = numeric(K)
    for(k in 1:K){
      Occ_dist_unobs[k] = sum(log(Forecast_Distribution_Unobs_v1[y[(data$N+1):N,k]==1,k])) - sum(mu.x.unobs[,k]) - ifelse(sum(y[(data$N+1):N,k]) == 0, 0, sum(log(1:sum(y[(data$N+1):N,k]))))
    }
    
  }
  
  Occasion_distribution = colSums(log(Forecast_Distribution_prevObs)) + Occ_dist_unobs
  
  Forecast_State_prevObs_list = C_ForecastState_AC(N = data$N, K = K, Nstates = N.states, 
                                              y_pix = data$y.pix.all.obs[,(data$K+1):(data$K+data$K.forecast)], p_move = p.move[,,1],
                                              FBProbs[,1:data$N]) # treat detection as perfect no matter the location
  
  Forecast_State_prevObs = Forecast_State_prevObs_list[[1]]
  Forecast_prevObs = Forecast_State_prevObs_list[[2]]
  
  # Brier score #
  Brier.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Quantiles.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Size.mat.prevObs = matrix(nrow = data$N, ncol = K)
  for(i in 1:data$N){
    for(k in 1:K){
      Brier.mat.prevObs[i,k] = (1 - Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], i])^2 + 
        sum(Forecast_prevObs[-data$y.pix.all.obs[i,k+data$K], i]^2)
      #Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k], i]) / N.states
      Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[,i][Forecast_prevObs[,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], i]])
      Size.mat.prevObs[i,k] = sum(Forecast_prevObs[,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], i])
    }
  }
  
  
  Forecast_State_Unobs_list = C_ForecastState_AC(N = nrow(data$y.pix.all.unobs), K = K, Nstates = N.states,
                                            y_pix = data$y.pix.all.unobs[,(data$K+1):(data$K+data$K.forecast)], p_move = p.move[,,1],
                                            FBProbs = matrix(rep(D.x.unobs, nrow(data$y.pix.all.unobs)), 
                                                             nrow = nrow(FBProbs), ncol = nrow(data$y.pix.all.unobs)))
  
  Forecast_State_Unobs = Forecast_State_Unobs_list[[1]]
  Forecast_unobs = Forecast_State_Unobs_list[[2]] / EN.unobs
  
  # Brier score #
  Brier.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Quantiles.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Size.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  for(i in 1:nrow(data$y.pix.all.unobs)){
    for(k in 1:K){
      Brier.mat.unobs[i,k] = (1 - Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], i])^2 + 
        sum(Forecast_unobs[-data$y.pix.all.unobs[i,k+data$K], i]^2)
      #Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[,i] >= Forecast_unobs[data$y.pix.all.unobs[i,k], i]) / N.states
      Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[,i][Forecast_unobs[,i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], i]])
      Size.mat.unobs[i,k] = sum(Forecast_unobs[,i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], i])
    }
  }
  
  Occ_state_dist_unobs = numeric(K)
  for(k in 1:K){
    Occ_state_dist_unobs[k] = sum(log(Forecast_State_Unobs[,k])) - sum(D.x.unobs) - sum(log(1:nrow(data$y.pix.all.unobs)))
  }
  
  Occasion_state_distribution = colSums(log(Forecast_State_prevObs)) + Occ_state_dist_unobs
  
  return(list(Observations = Occasion_distribution, States = Occasion_state_distribution,
              Observations_prevObs = colSums(log(Forecast_Distribution_prevObs)), States_prevObs = colSums(log(Forecast_State_prevObs)),
              Observations_unobs = Occ_dist_unobs, States_unobs = Occ_state_dist_unobs,
              Forecast_prevObs = Forecast_prevObs, Forecast_unobs = Forecast_unobs,
              Brier.prevObs = colMeans(Brier.mat.prevObs), Brier.unobs = colMeans(Brier.mat.unobs),
              Quantiles.prevObs = Quantiles.mat.prevObs, Quantiles.unobs = Quantiles.mat.unobs,
              Size.mat.prevObs = Size.mat.prevObs, Size.mat.unobs = Size.mat.unobs))
  
}




#############################################################################################################
#############################################################################################################
### Calculate likelihood of future observations and states (locations) given forecast for closed RW model ###
#############################################################################################################
#############################################################################################################
#' Calculate the likelihood of future encounter histories and movement trajectories under the forecast for a spatial capture-recapture model with a simple random walk movement model
#'
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW"
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW"
#' @param FBProbs Array of forward-backward probabilities; passed from ForwardBackward_Closed_SimpleRW
#'
#' @export
Forecast_Closed_SimpleRW = function(data, par, FBProbs, EN.unobs){

  y = data$y.forecast
  detDists = data$detDists.forecast 
  habitat = data$habitat
  detGrid = data$detGrid.forecast
  lines_arr = data$lines.arr.forecast
  distMat = data$distMat
  y.pix = data$y.pix.forecast
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  nCells = data$nCells
  nDetCells = data$nDetCells
  N = data$N.forecast
  K = data$K.forecast
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K # one transition between last occasion of y.training and first occasion of y.forecast
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ.forecast
  tr_b4_occ[1] = ifelse(tr_b4_occ[1] < 1, 1, tr_b4_occ[1])
  
  y.g0.covar = data$y.g0.covar.forecast
  y.ds.covar = data$y.ds.covar.forecast
  y.hs.covar = data$y.hs.covar.forecast
  y.hb.covar = data$y.hb.covar.forecast
  y.platform = data$y.platform.forecast
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  #gammaMod = data$gammaMod
  #phiMod = data$phiMod
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  
  if(n.d.covars == 0){
    p.move = C_pMove_Laplace_Full(K = K, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars)
  }else{
    p.move = C_pMove_Laplace_Full(K = K, nCells = nCells, habitat = habitat, distMat = distMat, 
                                  betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  # transition matrices between periods #
  t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  for(i in 1:K){  
    
    t.mat[,,i] = p.move[,,i]
    
  }
  
  ### Probability of not being detected given not detected during initial sampling period (y.training) ###
  ### Could make p.no.det.tot a vector and update it each iteration if that is faster? ###
  p.no.det.tot = array(dim = c(N.states, K, 2))
  
  # if the study begins with a first transition, the dimension of t.mat will need to be changed; it isn't that flexible yet #
  if(tr_b4_occ[1] > 0){
    p.no.det.tot[,1,1] = t.mat[,,1] %*% FBProbs[,data$N + 1]    # holds state distribution
    p.no.det.tot[,1,2] = p.no.det.tot[,1,1] * p.no.det.occ[,1]  # holds observation distribution
  } else{
    p.no.det.tot[,1,1] = FBProbs[,data$N + 1] 
    p.no.det.tot[,1,2] = p.no.det.tot[,1,1] * p.no.det.occ[,1]
  }
  
  for(i in 2:K){
    if(tr_b4_occ[i] > 1){
      temp = p.no.det.tot[,i-1,1]
      
      for(j in 1:(tr_b4_occ[i]-1)){
        temp = t.mat[,,i] %*% temp   # this assumes same transition for each day without observation
        
        ### MAY WANT TO BUILD IN FLEXIBILITY INTO t.mat TO ACCOMODATE DIFFERENT TRANSITIONS HERE ###
      }
      p.no.det.tot[,i,1] = temp
      p.no.det.tot[,i,2] = temp * p.no.det.occ[,i]
      
    } else{
      p.no.det.tot[,i,1] = t.mat[,,i] %*% p.no.det.tot[,i-1,1] 
      p.no.det.tot[,i,2] = p.no.det.tot[,i,1] * p.no.det.occ[,i]
    }
  }
  
  # probability of individual that is detected on an occasion given it was not detected during initial sampling period
  p.dot = 1 - sapply(1:K, function(X) sum(p.no.det.tot[,X,2]))
  
  # thinned Poisson point process rate parameter for each forecasted occasion # expected number of individuals detected that were not detected during initial sampling 
  mu.bar = p.dot * EN.unobs
  
  # pad data with an unobserved encounter history - representative of all individuals still unobserved through forecast sampling
  Forecast_Distribution_prevObs = C_ForecastObs_RW(N = data$N, K = K, Nstates = N.states, 
                                                   tr_b4_occ = tr_b4_occ,
                                                   mu_init = mu.init, 
                                                   y = y[1:data$N,],
                                                   y_pix = y.pix[1:data$N,],
                                                   detDists = detDists[1:data$N,],
                                                   y_platform = y.platform[1:data$N,],
                                                   y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                   n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                   g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                   p_no_det_occ = p.no.det.occ, t_mat = t.mat, FBProbs = FBProbs[,1:data$N])
  
  
  # expected number of unobserved individuals
  D.x.unobs = FBProbs[,data$N+1] * EN.unobs
  
  if(N - data$N == 0){
    Occ_dist_unobs = -mu.bar
  } else{
    Forecast_Distribution_Unobs_v1 = C_ForecastObs_RW(N = N - data$N, K = K, Nstates = N.states, 
                                                      tr_b4_occ = tr_b4_occ,
                                                      mu_init = mu.init, 
                                                      y = matrix(y[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_pix = matrix(y.pix[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      detDists = matrix(detDists[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_platform = matrix(y.platform[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                      y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                      n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                      g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                      p_no_det_occ = p.no.det.occ, t_mat = t.mat, 
                                                      FBProbs = matrix(rep(D.x.unobs, N - data$N), nrow = nrow(FBProbs), ncol = N - data$N))
    # above FBProbs is D.x.unobs following the Poisson point process likelihood listed in Glennie et al.
    
    Occ_dist_unobs = numeric(K)
    for(k in 1:K){
      Occ_dist_unobs[k] = sum(log(Forecast_Distribution_Unobs_v1[y[(data$N+1):N,k]==1,k])) - mu.bar[k] - ifelse(sum(y[(data$N+1):N,k]) == 0, 0, sum(log(1:sum(y[(data$N+1):N,k]))))
    }
  }
  
  Occasion_distribution = colSums(log(Forecast_Distribution_prevObs)) + Occ_dist_unobs
  
  Forecast_State_prevObs_list = C_ForecastState_RW(N = data$N, K = K, Nstates = N.states, tr_b4_occ = tr_b4_occ,
                                              y_pix = data$y.pix.all.obs[,(data$K+1):(data$K+data$K.forecast)], t_mat = t.mat,
                                              FBProbs[,1:data$N]) # treat detection as perfect no matter the location
  
  Forecast_State_prevObs = Forecast_State_prevObs_list[[1]]
  Forecast_prevObs = Forecast_State_prevObs_list[[2]]
  
  # Brier score #
  Brier.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Quantiles.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Size.mat.prevObs = matrix(nrow = data$N, ncol = K)
  for(i in 1:data$N){
    for(k in 1:K){
      Brier.mat.prevObs[i,k] = (1 - Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])^2 + 
        sum(Forecast_prevObs[-data$y.pix.all.obs[i,k+data$K], k, i]^2)
      #Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[, k, i] >= Forecast_prevObs[data$y.pix.all.obs[i,k], k, i]) / N.states
      Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i][Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i]])
      Size.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])
    }
  }
  
  Forecast_State_Unobs_list = C_ForecastState_RW(N = nrow(data$y.pix.all.unobs), K = K, Nstates = N.states, tr_b4_occ = tr_b4_occ,
                                            y_pix = data$y.pix.all.unobs[,(data$K+1):(data$K+data$K.forecast)], t_mat = t.mat,
                                            FBProbs = matrix(rep(D.x.unobs, nrow(data$y.pix.all.unobs)), 
                                                             nrow = nrow(FBProbs), ncol = nrow(data$y.pix.all.unobs)))
  Forecast_State_Unobs = Forecast_State_Unobs_list[[1]]
  Forecast_unobs = Forecast_State_Unobs_list[[2]] / EN.unobs
  
  # Brier score #
  Brier.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Quantiles.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Size.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  for(i in 1:nrow(data$y.pix.all.unobs)){
    for(k in 1:K){
      Brier.mat.unobs[i,k] = (1 - Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])^2 + 
        sum(Forecast_unobs[-data$y.pix.all.unobs[i,k+data$K], k, i]^2)
      #Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k], k, i]) / N.states
      Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[,k,i][Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i]])
      Size.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])
    }
  }
  
  Occ_state_dist_unobs = numeric(K)
  for(k in 1:K){
    Occ_state_dist_unobs[k] = sum(log(Forecast_State_Unobs[,k])) - sum(D.x.unobs) - sum(log(1:nrow(data$y.pix.all.unobs)))
  }
  
  Occasion_state_distribution = colSums(log(Forecast_State_prevObs)) + Occ_state_dist_unobs
  
  return(list(Observations = Occasion_distribution, States = Occasion_state_distribution,
              Observations_prevObs = colSums(log(Forecast_Distribution_prevObs)), States_prevObs = colSums(log(Forecast_State_prevObs)),
              Observations_unobs = Occ_dist_unobs, States_unobs = Occ_state_dist_unobs,
              Forecast_prevObs = Forecast_prevObs, Forecast_unobs = Forecast_unobs,
              Brier.prevObs = colMeans(Brier.mat.prevObs), Brier.unobs = colMeans(Brier.mat.unobs),
              Quantiles.prevObs = Quantiles.mat.prevObs, Quantiles.unobs = Quantiles.mat.unobs,
              Size.mat.prevObs = Size.mat.prevObs, Size.mat.unobs = Size.mat.unobs))
  
}



###############################################################################################################
###############################################################################################################
### Calculate likelihood of future observations and states (locations) given forecast for closed RWAC model ###
###############################################################################################################
###############################################################################################################
#' Calculate the likelihood of future encounter histories and movement trajectories under the forecast for a spatial capture-recapture model with a biased random walk movement model
#'
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_BiasedRW 
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_BiasedRW 
#' @param FBProbs Array of forward-backward probabilities; passed from ForwardBackward_Closed_BiasedRW
#'
#' @export
Forecast_Closed_BiasedRW = function(data, par, FBProbs, EN.unobs){
  
  y = data$y.forecast
  detDists = data$detDists.forecast
  detGrid = data$detGrid.forecast
  lines_arr = data$lines.arr.forecast
  y.pix = data$y.pix.forecast
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  
  nCells = data$nCells
  nDetCells = data$nDetCells
  nACCells = data$nACCells
  
  habitat = data$habitat
  habitat.AC = data$habitat.AC
  
  distMat = data$distMat
  distMat.AC = data$distMat.AC
  
  N = data$N.forecast
  K = data$K.forecast
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K - 1 #data$N.transitions # number of transitions, number of primary periods - 1
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ.forecast
  tr_b4_occ[1] = ifelse(tr_b4_occ[1] < 1, 1, tr_b4_occ[1])
  
  
  y.g0.covar = data$y.g0.covar.forecast
  y.ds.covar = data$y.ds.covar.forecast
  y.hs.covar = data$y.hs.covar.forecast
  y.hb.covar = data$y.hb.covar.forecast
  y.platform = data$y.platform.forecast
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  rho.init.indices = data$rho.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  rho = 1 / (1 + exp(-par[rho.init.indices]))
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  # have to integrate across all possible activity centers
  # probability of moving from cell i to cell j given the activity center is in cell k
  # u_t ~ N(u_t-1 * p + s * (1 - p), sigma^2I), where 0 <= p <= 1
  # need to loop over these calcs for each N.states s -> do this in C++
  if(n.d.covars == 0){
    p.move = C_pMove_RWAC_Laplace_Full(nCells = nCells, nACCells = nACCells, habitat = habitat,
                                       betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars, 
                                       rho = rho)
  }else{
    p.move = C_pMove_RWAC_Laplace_Full(nCells = nCells, nACCells = nACCells, habitat = habitat, 
                                       betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                       rho = rho)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat.AC[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat.AC[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  
  # transition matrices between periods #
  # t.mat should be 4d, cell to, cell from, activity center cell, transitions
  # likelihood and calculation of thinning rate (p.dot) need to be integrated over all possible activity centers
  #t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  
  ### THE ABOVE (mu.init) IS THE POINT PROCESS FOR ACTIVITY CENTERS ###
  ### ALSO NEED TO DESCRIBE INITIAL LOCATION ABOUT ACTIVITY CENTERS ###
  ### USE A MATRIX INSTEAD OF A VECTOR TO DESCRIBE ACTIVITY CENTER AND LOCATION STATE ###
  ### LAPLACE INITIAL v2 -> changed from marginalizing over ac's to marginilizing over all previous locations and ac's ###
  ### v2 ALSO DID NOT WORK ###
  ### i.e., treat the situation like the individual's location was uniformly distributed during the time step before the first occasion ###
  ### having issues estimating rho, think it is because there is too much influence from first occasion on activity center location ###
  ### could also use a model where the first location is distributed uniformly across the state space, regardless of the activity center location ###
  if(n.d.covars == 0){
    p.init.state.full = C_pMove_RWAC_Laplace_Initial(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
                                                     betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                     p_init_state_ac = p.init.state)
    #p.init.state.full = C_pMove_RWAC_Laplace_Initial_v2(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
    #                                                 betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
    #                                                 p_init_state_ac = p.init.state, rho = rho)
    
    # try uniformly distributed initial location and activity center
    #p.init.state.full = matrix(rep(p.init.state/nCells, nCells), nrow = nCells, ncol = nACCells, byrow = T)
    
  } else{
    p.init.state.full = C_pMove_RWAC_Laplace_Initial(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
                                                     betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                     p_init_state_ac = p.init.state)
    #p.init.state.full = C_pMove_RWAC_Laplace_Initial_v2(nCells = nCells, nACCells = nACCells, habitat = habitat, distMat_AC = distMat.AC,
    #                                                 betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
    #                                                 p_init_state_ac = p.init.state, rho = rho)
    
    # try uniformly distributed initial location
    #p.init.state.full = matrix(rep(p.init.state/nCells, nCells), nrow = nCells, ncol = nACCells, byrow = T)
  }
  
  # p.init.state.full: activity center is column index; movement cell is row index 
  
  # probability of being detected on an occasion given the individual was not detected during the initial sampling (training set)
  p.dots = C_pDots_RWAC(K = K, tr_b4_occ = tr_b4_occ,
                        nCells = nCells, nACCells = nACCells, p_init_state_full = FBProbs[,,data$N+1],
                        pMove = p.move, p_no_det_occ = p.no.det.occ)


  
  # thinned Poisson point process rate parameter for each forecasted occasion # expected number of individuals detected that were not detected during initial sampling
  mu.bar = p.dots[[1]] * EN.unobs
  
  
  # pad data with an unobserved encounter history - representative of all individuals still unobserved through forecast sampling
  Forecast_Distribution_prevObs = C_ForecastObs_RWAC(N = data$N, K = K, nCells = nCells, nACCells = nACCells,
                                                   tr_b4_occ = tr_b4_occ,
                                                   mu_init = mu.init, FBProbs = FBProbs[,,1:data$N],
                                                   y = y[1:data$N,],
                                                   y_pix = y.pix[1:data$N,],
                                                   detDists = detDists[1:data$N,],
                                                   y_platform = y.platform[1:data$N,],
                                                   y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                   n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                   g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                   p_no_det_occ = p.no.det.occ, t_mat = p.move)
  
  
  # expected number of unobserved individuals
  D.x.unobs = FBProbs[,,data$N+1] * EN.unobs
  
  FBProbs.unobs = array(dim = c(dim(FBProbs)[1:2],100 - data$N))
  for(i in 1:(100-data$N)){
    FBProbs.unobs[,,i] = D.x.unobs
  }
  
  if(N - data$N == 0){
    Occ_dist_unobs = -mu.bar
  } else{
    Forecast_Distribution_Unobs_v1 = C_ForecastObs_RWAC(N = N - data$N, K = K, nCells = nCells, nACCells = nACCells, 
                                                        tr_b4_occ = tr_b4_occ,
                                                        mu_init = mu.init, FBProbs = FBProbs.unobs,
                                                        y = matrix(y[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_pix = matrix(y.pix[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        detDists = matrix(detDists[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_platform = matrix(y.platform[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                        n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                        g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                        p_no_det_occ = p.no.det.occ, t_mat = p.move)
    # above FBProbs is D.x.unobs following the Poisson point process likelihood listed in Glennie et al.
    
    Occ_dist_unobs = numeric(K)
    for(k in 1:K){
      Occ_dist_unobs[k] = sum(log(Forecast_Distribution_Unobs_v1[y[(data$N+1):N,k]==1,k])) - mu.bar[k] - ifelse(sum(y[(data$N+1):N,k]) == 0, 0, sum(log(1:sum(y[(data$N+1):N,k]))))
    }
  }
  
  Occasion_distribution = colSums(log(Forecast_Distribution_prevObs)) + Occ_dist_unobs
  
  Forecast_State_prevObs_list = C_ForecastState_RWAC(N = data$N, K = K, nCells = nCells, nACCells = nACCells, 
                                                   tr_b4_occ = tr_b4_occ, mu_init = mu.init,
                                                   FBProbs = FBProbs[,,1:data$N],
                                                   y_pix = data$y.pix.all.obs[,(data$K+1):(data$K+data$K.forecast)], 
                                                   t_mat = p.move) # treat detection as perfect no matter the location
  

  
  Forecast_State_prevObs = Forecast_State_prevObs_list[[1]]
  Forecast_prevObs = Forecast_State_prevObs_list[[2]]
  
  # Brier score #
  Brier.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Quantiles.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Size.mat.prevObs = matrix(nrow = data$N, ncol = K)
  for(i in 1:data$N){
    for(k in 1:K){
      Brier.mat.prevObs[i,k] = (1 - Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])^2 + 
        sum(Forecast_prevObs[-data$y.pix.all.obs[i,k+data$K], k, i]^2)
      #Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[, k, i] >= Forecast_prevObs[data$y.pix.all.obs[i,k], k, i]) / N.states
      Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i][Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i]])
      Size.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])
    }
  }
  
  Forecast_State_Unobs_list = C_ForecastState_RWAC(N = nrow(data$y.pix.all.unobs), K = K, nCells = nCells, nACCells = nACCells,
                                                   tr_b4_occ = tr_b4_occ, mu_init = mu.init,
                                                   FBProbs = FBProbs.unobs,
                                                   y_pix = data$y.pix.all.unobs[,(data$K+1):(data$K+data$K.forecast)], 
                                                   t_mat = p.move)
  
  Forecast_State_Unobs = Forecast_State_Unobs_list[[1]]
  Forecast_unobs = Forecast_State_Unobs_list[[2]] / EN.unobs
  
  # Brier score #
  Brier.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Quantiles.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Size.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  for(i in 1:nrow(data$y.pix.all.unobs)){
    for(k in 1:K){
      Brier.mat.unobs[i,k] = (1 - Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])^2 + 
        sum(Forecast_unobs[-data$y.pix.all.unobs[i,k+data$K], k, i]^2)
      #Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k], k, i]) / N.states
      Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[,k,i][Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i]])
      Size.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])
    }
  }
  
  Occ_state_dist_unobs = numeric(K)
  for(k in 1:K){
    Occ_state_dist_unobs[k] = sum(log(Forecast_State_Unobs[,k])) - sum(D.x.unobs) - sum(log(1:nrow(data$y.pix.all.unobs)))
  }
  
  Occasion_state_distribution = colSums(log(Forecast_State_prevObs)) + Occ_state_dist_unobs
  
  return(list(Observations = Occasion_distribution, States = Occasion_state_distribution,
              Observations_prevObs = colSums(log(Forecast_Distribution_prevObs)), States_prevObs = colSums(log(Forecast_State_prevObs)),
              Observations_unobs = Occ_dist_unobs, States_unobs = Occ_state_dist_unobs,
              Forecast_prevObs = Forecast_prevObs, Forecast_unobs = Forecast_unobs,
              Brier.prevObs = colMeans(Brier.mat.prevObs), Brier.unobs = colMeans(Brier.mat.unobs),
              Quantiles.prevObs = Quantiles.mat.prevObs, Quantiles.unobs = Quantiles.mat.unobs,
              Size.mat.prevObs = Size.mat.prevObs, Size.mat.unobs = Size.mat.unobs))
  
}



##############################################################################################################
##############################################################################################################
### Calculate likelihood of future observations and states (locations) given forecast for closed CRW model ###
##############################################################################################################
##############################################################################################################
#' Calculate the likelihood of future encounter histories and movement trajectories under the forecast for a spatial capture-recapture model with a correlated random walk movement model
#'
#' @param data Named list of data; passed from the output of a call to scr_lt_Closed_CorrelatedRW 
#' @param par Vector of maximum likelihood parameter estimates; passed from the output of a call to scr_lt_Closed_CorrelatedRW 
#' @param FBProbs Array of forward-backward probabilities; passed from ForwardBackward_Closed_CorrelatedRW
#'
#' @export
Forecast_Closed_CorrelatedRW = function(data, par, FBProbs, EN.unobs){  
  
  y = data$y.forecast
  detDists = data$detDists.forecast
  detGrid = data$detGrid.forecast
  lines_arr = data$lines.arr.forecast
  y.pix = data$y.pix.forecast
  Area = data$Area          # state space pixel area
  detArea = data$detArea    # detection grid pixel area; will need to be included in detGrid if cell size varies
  truncDist = data$truncDist  # how far from survey effort should detection probability be fixed to 0
  
  nCells = data$nCells
  nDetCells = data$nDetCells
  nACCells = data$nACCells
  
  habitat = data$habitat
  habitat.AC = data$habitat.AC
  
  distMat = data$distMat
  distMat.AC = data$distMat.AC
  
  N = data$N.forecast
  K = data$K.forecast
  log.fact.N = data$log.fact.N
  
  N.states = data$N.states  # number of spatial states
  N.transitions = K - 1 #data$N.transitions # number of transitions, number of primary periods - 1
  SecOcc = data$SecOcc
  tr_b4_occ = data$tr_b4_occ.forecast
  tr_b4_occ[1] = ifelse(tr_b4_occ[1] < 1, 1, tr_b4_occ[1])
  
  y.g0.covar = data$y.g0.covar.forecast
  y.ds.covar = data$y.ds.covar.forecast
  y.hs.covar = data$y.hs.covar.forecast
  y.hb.covar = data$y.hb.covar.forecast
  y.platform = data$y.platform.forecast
  
  d.init.indices = data$d.init.indices
  m.init.indices = data$m.init.indices
  g0.init.indices = data$g0.init.indices
  ds.init.indices = data$ds.init.indices
  hs.init.indices = data$hs.init.indices
  hb.init.indices = data$hb.init.indices
  gamma.init.indices = data$gamma.init.indices
  beta.init.indices = data$beta.init.indices
  
  n.d.covars = data$n.d.covars
  n.ds.covars = data$n.ds.covars
  n.g0.covars = data$n.g0.covars
  n.hs.covars = data$n.hs.covars
  n.hb.covars = data$n.hb.covars
  
  d.covar.cols = data$d.covar.cols
  ds.covar.cols = data$ds.covar.cols
  g0.covar.cols = data$g0.covar.cols
  hs.covar.cols = data$hs.covar.cols
  hb.covar.cols = data$hb.covar.cols
  
  ### parameters ###
  ## need parameters: FirstOccIn, gamma vector, phi vector ##
  
  betas = par[d.init.indices]   # density
  
  betaMove = exp(par[m.init.indices])
  
  g0s = par[g0.init.indices]  # g0
  
  if(ds.init.indices[1] == -1){
    sigmaDets = 1 # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    sigmaDets = par[ds.init.indices]  # detection scale
  }
  
  if(hs.init.indices[1] == -1){
    hSigma = 1    # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hSigma = par[hs.init.indices]
  }
  
  if(hb.init.indices[1] == -1){
    hBeta = 1   # not in the model, fix for C_pDet as it requires it as an argument
  } else{
    hBeta = par[hb.init.indices]
  }
  
  gamma = 1 / (1 + exp(-par[gamma.init.indices]))
  beta = (2*pi / (1 + exp(-par[beta.init.indices]))) - pi
  
  rotation = gamma * matrix(c(cos(beta), sin(beta), -sin(beta), cos(beta)), nrow = 2)
  
  
  p.det.arr = C_pDet_Full(K_tot = K, nCells = nCells, nDetCells = nDetCells, lines_arr = lines_arr, detGrid = detGrid, 
                          sigmaDets = c(sigmaDets), g0s = c(g0s), hSigmas = c(hSigma), hBetas = c(hBeta), 
                          detArea = detArea, truncDist = truncDist,
                          SDcols = ds.covar.cols, nSDcols = n.ds.covars,
                          g0cols = g0.covar.cols, nG0cols = n.g0.covars,
                          hscols = hs.covar.cols, nhscols = n.hs.covars,
                          hbcols = hb.covar.cols, nhbcols = n.hb.covars)
  
  
  
  
  # have to integrate across all possible activity centers
  # probability of moving from cell i to cell j given the activity center is in cell k
  # u_t ~ N(u_t-1 * p + s * (1 - p), sigma^2I), where 0 <= p <= 1
  # need to loop over these calcs for each N.states s -> do this in C++
  if(n.d.covars == 0){
    p.move = C_pMove_CorrelatedRW_Laplace_Full(nCells = nCells, habitat = habitat,
                                               betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars, 
                                               rotation = rotation)
  }else{
    p.move = C_pMove_CorrelatedRW_Laplace_Full(nCells = nCells, habitat = habitat, 
                                               betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                               rotation = rotation)
  }
  
  
  
  
  # probability of not being detected on each primary and secondary occasion for each activity center
  #  p.no.det.given.ac = matrix(nrow = N.states, ncol = K.prim)
  #  p.no.det.occ = matrix(nrow = N.states, ncol = K.tot) # issue if variable number of secondry occasion
  ### CHANGING SO THAT THERE ARE ONLY SPATIAL STATES ### 
  #p.no.det.given.ac = matrix(nrow = nCells, ncol = K.prim)
  #p.no.det.occ = matrix(nrow = nCells, ncol = K.tot) # issue if variable number of secondry occasion
  ### THESE ARE THE SAME WHEN THERE IS ONLY ONE OCCASION PER PRIMARY PERIOD
  p.no.det.occ = 1 - p.det.arr
  
  
  # expected density in each cell at start of each season
  if(n.d.covars == 0){
    p.D <- exp(betas[1])
  } else if(n.d.covars > 0){
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat.AC[,d.covar.cols[1,x]]))
    p.D <- exp(ut.D)
  }
  
  # expected activity center abundance in each cell
  D.x = p.D * habitat.AC[,"Area"]
  
  # expected total abundance
  D.bar = sum(D.x)
  
  # initial spatial distribution for each individual
  p.init.state = D.x / D.bar
  
  # intensity of point process
  mu.init = D.x
  
  
  # transition matrices between periods #
  # t.mat should be 4d, cell to, cell from, activity center cell, transitions
  # likelihood and calculation of thinning rate (p.dot) need to be integrated over all possible activity centers
  #t.mat = array(dim = c(N.states, N.states, N.transitions))   
  spatialStates = 1:N.states
  
  
  ### THE ABOVE (mu.init) IS THE POINT PROCESS FOR ACTIVITY CENTERS ###
  ### ALSO NEED TO DESCRIBE INITIAL LOCATION ABOUT ACTIVITY CENTERS ###
  ### USE A MATRIX INSTEAD OF A VECTOR TO DESCRIBE ACTIVITY CENTER AND LOCATION STATE ###
  if(n.d.covars == 0){
    p.init.state.full = C_pMove_CorrelatedRW_Laplace_Initial(nCells = nCells, habitat = habitat, distMat = distMat,
                                                             betaMoves = c(betaMove), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                             p_init_state = p.init.state)
    
    
  } else{
    p.init.state.full = C_pMove_CorrelatedRW_Laplace_Initial(nCells = nCells, habitat = habitat, distMat = distMat,
                                                             betaMoves = c(betaMove, betas[2:(n.d.covars+1)]), covarCols = d.covar.cols, nCovars = n.d.covars,
                                                             p_init_state = p.init.state)
    
  }
  
  # p.init.state.full: activity center is column index; movement cell is row index 
  
  
  
  # probability of having not yet been detected on each occasion given initial activity center and locations  
  pNoDet_CorrelatedRW_out = C_pNoDet_CorrelatedRW(K = K, tr_b4_occ = tr_b4_occ,
                                                  nCells = nCells, p_init_state_full = p.init.state.full,
                                                  pMove = p.move, p_no_det_occ = p.no.det.occ)
  # probability of being detected at least once 
  p.dot = pNoDet_CorrelatedRW_out[[1]]
  p.no.det.tot = pNoDet_CorrelatedRW_out[[2]]
  
  
  ### LEFT OFF HERE - REMOVE ABOVE p.dot AND CALCULATE PROBABILITY OF BEING DETECTED EACH OCCASION GIVEN NOT DETECTED DURING TRAINING SET
  # probability of being detected on an occasion given the individual was not detected during the initial sampling (training set)
  p.dots = C_pDots_CRW(K = K, tr_b4_occ = tr_b4_occ,
                        nCells = nCells, p_init_state_full = FBProbs[,,data$N+1],
                        pMove = p.move, p_no_det_occ = p.no.det.occ)
  
  
  
  # thinned Poisson point process rate parameter for each forecasted occasion # expected number of individuals detected that were not detected during initial sampling
  mu.bar = p.dots[[1]] * EN.unobs
  
  
  
  # pad data with an unobserved encounter history - representative of all individuals still unobserved through forecast sampling
  Forecast_Distribution_prevObs = C_ForecastObs_CRW(N = data$N, K = K, nCells = nCells, 
                                                     tr_b4_occ = tr_b4_occ,
                                                     mu_init = mu.init, FBProbs = FBProbs[,,1:data$N],
                                                     y = y[1:data$N,],
                                                     y_pix = y.pix[1:data$N,],
                                                     detDists = detDists[1:data$N,],
                                                     y_platform = y.platform[1:data$N,],
                                                     y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                     n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                     g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                     p_no_det_occ = p.no.det.occ, t_mat = p.move)
  
  
  # expected number of unobserved individuals
  D.x.unobs = FBProbs[,,data$N+1] * EN.unobs
  
  FBProbs.unobs = array(dim = c(dim(FBProbs)[1:2],100 - data$N))
  for(i in 1:(100-data$N)){
    FBProbs.unobs[,,i] = D.x.unobs
  }
  
  if(N - data$N == 0){
    Occ_dist_unobs = -mu.bar
  } else{
    Forecast_Distribution_Unobs_v1 = C_ForecastObs_CRW(N = N - data$N, K = K, nCells = nCells,  
                                                        tr_b4_occ = tr_b4_occ,
                                                        mu_init = mu.init, FBProbs = FBProbs.unobs,
                                                        y = matrix(y[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_pix = matrix(y.pix[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        detDists = matrix(detDists[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_platform = matrix(y.platform[(data$N+1):N,], ncol = K, nrow = N - data$N),
                                                        y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                                        n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                                        g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                                        p_no_det_occ = p.no.det.occ, t_mat = p.move)
    # above FBProbs is D.x.unobs following the Poisson point process likelihood listed in Glennie et al.
    
    Occ_dist_unobs = numeric(K)
    for(k in 1:K){
      Occ_dist_unobs[k] = sum(log(Forecast_Distribution_Unobs_v1[y[(data$N+1):N,k]==1,k])) - mu.bar[k] - ifelse(sum(y[(data$N+1):N,k]) == 0, 0, sum(log(1:sum(y[(data$N+1):N,k]))))
    }
  }
  
  Occasion_distribution = colSums(log(Forecast_Distribution_prevObs)) + Occ_dist_unobs
  
  Forecast_State_prevObs_list = C_ForecastState_CRW(N = data$N, K = K, nCells = nCells,  
                                                     tr_b4_occ = tr_b4_occ, mu_init = mu.init,
                                                     FBProbs = FBProbs[,,1:data$N],
                                                     y_pix = data$y.pix.all.obs[,(data$K+1):(data$K+data$K.forecast)], 
                                                     t_mat = p.move) # treat detection as perfect no matter the location
  
  
  
  Forecast_State_prevObs = Forecast_State_prevObs_list[[1]]
  Forecast_prevObs = Forecast_State_prevObs_list[[2]]
  
  # Brier score #
  Brier.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Quantiles.mat.prevObs = matrix(nrow = data$N, ncol = K)
  Size.mat.prevObs = matrix(nrow = data$N, ncol = K)
  for(i in 1:data$N){
    for(k in 1:K){
      Brier.mat.prevObs[i,k] = (1 - Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])^2 + 
        sum(Forecast_prevObs[-data$y.pix.all.obs[i,k+data$K], k, i]^2)
      #Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[, k, i] >= Forecast_prevObs[data$y.pix.all.obs[i,k], k, i]) / N.states
      Quantiles.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i][Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i]])
      Size.mat.prevObs[i,k] = sum(Forecast_prevObs[,k,i] >= Forecast_prevObs[data$y.pix.all.obs[i,k+data$K], k, i])
    }
  }
  
  Forecast_State_Unobs_list = C_ForecastState_CRW(N = nrow(data$y.pix.all.unobs), K = K, nCells = nCells, 
                                                   tr_b4_occ = tr_b4_occ, mu_init = mu.init,
                                                   FBProbs = FBProbs.unobs,
                                                   y_pix = data$y.pix.all.unobs[,(data$K+1):(data$K+data$K.forecast)], 
                                                   t_mat = p.move)
  
  Forecast_State_Unobs = Forecast_State_Unobs_list[[1]]
  Forecast_unobs = Forecast_State_Unobs_list[[2]] / EN.unobs
  
  # Brier score #
  Brier.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Quantiles.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  Size.mat.unobs = matrix(nrow = nrow(data$y.pix.all.unobs), ncol = K)
  for(i in 1:nrow(data$y.pix.all.unobs)){
    for(k in 1:K){
      Brier.mat.unobs[i,k] = (1 - Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])^2 + 
        sum(Forecast_unobs[-data$y.pix.all.unobs[i,k+data$K], k, i]^2)
      #Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k], k, i]) / N.states
      Quantiles.mat.unobs[i,k] = sum(Forecast_unobs[,k,i][Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i]])
      Size.mat.unobs[i,k] = sum(Forecast_unobs[, k, i] >= Forecast_unobs[data$y.pix.all.unobs[i,k+data$K], k, i])
    }
  }
  
  Occ_state_dist_unobs = numeric(K)
  for(k in 1:K){
    Occ_state_dist_unobs[k] = sum(log(Forecast_State_Unobs[,k])) - sum(D.x.unobs) - sum(log(1:nrow(data$y.pix.all.unobs)))
  }
  
  Occasion_state_distribution = colSums(log(Forecast_State_prevObs)) + Occ_state_dist_unobs
  
  return(list(Observations = Occasion_distribution, States = Occasion_state_distribution,
              Observations_prevObs = colSums(log(Forecast_Distribution_prevObs)), States_prevObs = colSums(log(Forecast_State_prevObs)),
              Observations_unobs = Occ_dist_unobs, States_unobs = Occ_state_dist_unobs,
              Forecast_prevObs = Forecast_prevObs, Forecast_unobs = Forecast_unobs,
              Brier.prevObs = colMeans(Brier.mat.prevObs), Brier.unobs = colMeans(Brier.mat.unobs),
              Quantiles.prevObs = Quantiles.mat.prevObs, Quantiles.unobs = Quantiles.mat.unobs,
              Size.mat.prevObs = Size.mat.prevObs, Size.mat.unobs = Size.mat.unobs))
}

