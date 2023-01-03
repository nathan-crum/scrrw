#############################################
#############################################
#### Function to optimize the likelihood ####
#############################################
#############################################
#' Fit a closed population spatial capture-recapture model with a simple random walk and line transect sampling
#'  
#' @param data Named list of input data
#' @param density.movement Formulation for the initial population density and movement models
#' @param g0 Formulation for the detection model - probability of detection for individuals on the trackline
#' @param det.scale Formulation for the detection model - scale parameter of half-normal detection function; -1 if half-normal detection function not used
#' @param hSigma Formulation for the detection model - parameter of hazard detection function; -1 if hazard detection function not used
#' @param hBeta Formulation for the detection model - parameter of hazard detection function; -1 if hazard detection function not used
#' @param MOVEMENT.MODEL Character string, either "RW" for simple random walk or "AC" for the traditional activty center model
#'
#' @export
scr_lt_Closed_SimpleRW = function(data, density.movement = ~1, g0 = ~1, det.scale = ~1, hSigma = ~-1, hBeta = ~-1, MOVEMENT.MODEL = "RW"){
  
  # interpret formulas - only handling additive covariate effects, supply interactive and polynomial covariates if part of the model
  
  ### Density/Movement ###
  dt = terms(density.movement)
  d.covars = rownames(attr(dt, "factors"))
  n.d.covars = length(d.covars)
  if(n.d.covars > 0){
    d.covar.cols = matrix(nrow = data$K.prim, ncol = n.d.covars)
    d.matches = match(d.covars, colnames(data$habitat))
    d.unmatched = d.covars[is.na(d.matches)]
    if(any(!is.na(d.matches))){
      for(i in 1:sum(!is.na(d.matches))){
        d.covar.cols[,i] = d.matches[!is.na(d.matches)][i]
      }
      nextcol = sum(!is.na(d.matches)) + 1
    } else{
      nextcol = 1
    }
    
    if(length(d.unmatched) > 0){
      for(i in 1:length(d.unmatched)){
        for(k in 1:data$K.prim){
          covar_occ = paste0(d.unmatched[i], "_", k)
          d.covar.cols[k,nextcol] = match(covar_occ, colnames(data$habitat))
        }
        nextcol = nextcol + 1
      }
    }
    
    if(any(is.na(d.covar.cols))){
      print("Error: Density covariates incorrectly specified")
      return(NULL)
    }
    
  } else{
    d.covar.cols = -1
  }
  
  
  ### g0 ###
  gt = terms(g0)                            ### g0 ###
  g0.covars = rownames(attr(gt, "factors"))
  g0.covar.cols = match(g0.covars, colnames(data$detGrid))
  n.g0.covars = length(g0.covars)
  
  
  ### det.scale ###
  dst = terms(det.scale)                    ### det.scale ###
  if(attr(dst, "intercept") == 0){  # removed from model
    ds.covars = NULL
    ds.covar.cols = -1
    n.ds.covars = 0
  } else{
    ds.covars = rownames(attr(dst, "factors"))
    ds.covar.cols = match(ds.covars, colnames(data$detGrid))
    n.ds.covars = length(ds.covars)
  }
  
  hst = terms(hSigma)
  if(attr(hst, "intercept") == 0){
    hs.covars = NULL
    hs.covar.cols = -1
    n.hs.covars = 0
  } else{
    hs.covars = rownames(attr(hst, "factors"))
    hs.covar.cols = match(hs.covars, colnames(data$detGrid))
    n.hs.covars = length(hs.covars)
  }
  
  hbt = terms(hBeta)
  if(attr(hbt, "intercept") == 0){
    hb.covars = NULL
    hb.covar.cols = -1
    n.hb.covars = 0
  } else{
    hb.covars = rownames(attr(hbt, "factors"))
    hb.covar.cols = match(hb.covars, colnames(data$detGrid))
    n.hb.covars = length(hb.covars)
  }

  
  
  
  # add in the number of covariates for each submodel #
  data$n.d.covars = n.d.covars
  data$n.g0.covars = n.g0.covars
  data$n.ds.covars = n.ds.covars
  data$n.hs.covars = n.hs.covars
  data$n.hb.covars = n.hb.covars
  
  # add in the columns corresponding to each covariate
  data$d.covar.cols = d.covar.cols
  data$g0.covar.cols = g0.covar.cols
  data$ds.covar.cols = ds.covar.cols
  data$hs.covar.cols = hs.covar.cols
  data$hb.covar.cols = hb.covar.cols
  
  ### get initial values for parameters ###
  par.out = SCR_lt_Closed_SimpleRW_inits(data)
  
  # update covariate columns if there aren't any covariates
  if(n.d.covars == 0){
    data$d.covar.cols = matrix(0, nrow = 1, ncol = 1)
  }
  
  if(n.ds.covars == 0){
    data$ds.covar.cols = 0
  }
  
  if(n.g0.covars == 0){
    data$g0.covar.cols = 0
  }
  
  if(n.hs.covars == 0){
    data$hs.covar.cols = 0
  }
  
  if(n.hb.covars == 0){
    data$hb.covar.cols = 0
  }
  
  
  ### add the indices for each set of parameters to data ###
  data$ds.init.indices = par.out$ds.init.indices
  data$g0.init.indices = par.out$g0.init.indices
  data$m.init.indices = par.out$m.init.indices
  data$d.init.indices = par.out$d.init.indices
  data$hs.init.indices = par.out$hs.init.indices
  data$hb.init.indices = par.out$hb.init.indices
  
  if(MOVEMENT.MODEL == "RW"){
    results.BFGS <- optim(par = par.out$par, fn = SCR_lt_Closed_SimpleRW_loglklhd, data = data, hessian = T,
                          method = "BFGS", control = list(maxit = 10000))
  } else if(MOVEMENT.MODEL == "AC"){
    results.BFGS <- optim(par = par.out$par, fn = SCR_lt_Closed_AC_loglklhd, data = data, hessian = T,
                          method = "BFGS", control = list(maxit = 10000))
  } else if(MOVEMENT.MODEL == "RWAC"){
    results.BFGS <- optim(par = par.out$par, fn = SCR_lt_Closed_RWAC_loglklhd, data = data, hessian = T,
                          method = "BFGS", control = list(maxit = 10000))
  }

  
  
  return(list(results = results.BFGS, data = data))
  
}



#########################################################
#########################################################
#### Generate initial values for likelihood function ####
#########################################################
#########################################################
#' Generate initial values for the spatial capture-recapture simple random walk model
#' 
#' @param data Named list of input data and parameter indices; passed from scr_lt_Closed_SimpleRW
#'
#' @export
SCR_lt_Closed_SimpleRW_inits = function(data){
  
  # initial values for betas
  par = numeric()
  par[1] = log(data$N/sum(data$habitat[,"Area"]))
  if(data$n.d.covars > 0){
    par[2:(data$n.d.covars+1)] = 0
  }
  d.init.indices = 1:length(par)
  
  # movement scale parameter - 1 per group
  m.init.indices = length(par) + 1
  par[m.init.indices] = log(data$mean.ac.move)
  
  # g0 initial values
  if(data$n.g0.covars == 0){
    g0.init.indices = length(par) + 1
    par[g0.init.indices] = 0  # qlogis(0.5)
  } else if(data$n.g0.covars > 0){
    g0.init.indices = (length(par)+1):(length(par)+data$n.g0.covars+1)
    par[g0.init.indices] = 0
  }
  
  # half normal scale initial values
  if(data$n.ds.covars == 0){
    if(length(data$ds.covar.cols) == 0){
      ds.init.indices = length(par) + 1
      par[ds.init.indices] = log(quantile(data$detDists, 0.68, na.rm=T))
    } else if(data$ds.covar.cols == -1){   # removed from the model
      ds.init.indices = -1
    } 
  } else if(data$n.ds.covars > 0){
    ds.init.indices = (length(par) + 1):(length(par) + data$n.ds.covars + 1)
    par[ds.init.indices[1]] = log(quantile(data$detDists, 0.68, na.rm=T))
    par[ds.init.indices[2:length(ds.init.indices)]] = 0
  }
  
  # hazard scale initial values
  if(data$n.hs.covars == 0){
    if(length(data$hs.covar.cols) == 0){
      hs.init.indices = length(par) + 1
      par[hs.init.indices] = log(quantile(data$detDists, 0.68, na.rm=T))
    } else if(data$hs.covar.cols == -1){
      hs.init.indices = -1
    }
  } else if(data$n.hs.covars > 0){
    hs.init.indices = (length(par) + 1):(length(par) + data$n.hs.covars + 1)
    par[hs.init.inidices[1]] = log(quantile(data$detDists, 0.68, na.rm=T))
    par[hs.init.indices[2:length(hs.init.indices)]] = 0
  }
  
  # hazard beta initial values
  if(data$n.hb.covars == 0){
    if(length(data$hb.covar.cols) == 0){
      hb.init.indices = length(par) + 1
      par[hb.init.indices] = log(1)
    } else if(data$hb.covar.cols == -1){
      hb.init.indices = -1
    } 
  } else if(data$n.hb.covars > 0){
    hb.init.indices = (length(par) + 1):(length(par) + data$n.hb.covars + 1)
    par[hb.init.inidices[1]] = log(1)
    par[hb.init.indices[2:length(hb.init.indices)]] = 0
  }

  
  
  return(list(par = par, d.init.indices = d.init.indices, m.init.indices = m.init.indices, 
              g0.init.indices = g0.init.indices, ds.init.indices = ds.init.indices,
              hs.init.indices = hs.init.indices, hb.init.indices = hb.init.indices))
  
}



##################################
##################################
### Activity Center Likelihood ###
##################################
##################################
#' Likelihood function for the closed population, activity center, spatial capture-recapture model
#' 
#' @param data Named list of input data; passed from scr_lt_Closed_SimpleRW
#' @param par Numeric vector of parameters; passed from scr_lt_Closed_SimpleRW
#'
#' @export
SCR_lt_Closed_AC_loglklhd = function(data, par){  
  
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
  
  
  ### Calculate sum(log-likelihood of individuals encounter histories) ###
  ind.likelihood = C_indLklhd_AC(N = N, K = K, Nstates = N.states,  
                                 mu_init = mu.init, y = y, y_pix = y.pix, detDists = detDists, y_platform = y.platform,
                                 y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                 n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                 g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                 p_no_det_occ = p.no.det.occ, p_move = p.move[,,1])
  
  
  ### DIFFERENCE FROM SCR BOOK FORMULATION ###
  #log.Likelihood = sum(logProb.ind) + N * log(D * area.state.space) + (-D * area.state.space * p.star)
  # Glennie likelihood of thinned Poisson point process (thinned by the detection function):
  # exp(-mu.bar)/factorial(number of individuals detected) * product of encounter history likelihoods
  # this is essentially a Poisson pmf: exp(-lambda) * lambda^x/factorial(x)
  # lambda is the rate parameter of the Poisson distribution, equivalent to mu.bar, the rate of the thinned Poisson point process
  # x is the number of detected individuals
  # lambda^x is the product of x encounter histories
  #log.Likelihood = sum(logProb.ind) - mu.bar - log.fact.N
  log.Likelihood = ind.likelihood - mu.bar - log.fact.N
  
  if(is.na(log.Likelihood)){
    print(par)
  }
  
  return(-log.Likelihood)
}




##############################
##############################
### Random Walk Likelihood ###
##############################
##############################
#' Likelihood function for the closed population, simple random walk, spatial capture-recapture model
#' 
#' @param data Named list of input data; passed from scr_lt_Closed_SimpleRW
#' @param par Numeric vector of parameters; passed from scr_lt_Closed_SimpleRW
#'
#' @export
SCR_lt_Closed_SimpleRW_loglklhd = function(data, par){  
  
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
  p.no.det.tot = matrix(nrow = N.states, ncol = K)

  # if the study begins with a first transition, the dimension of t.mat will need to be changed; it isn't that flexible yet #
  if(tr_b4_occ[1] > 0){
    p.no.det.tot[,1] = t.mat[,,1] %*% p.init.state * p.no.det.occ[,1]
  } else{
    p.no.det.tot[,1] = p.init.state * p.no.det.occ[,1]
  }
  
  for(i in 2:K){
    if(tr_b4_occ[i] > 1){
      temp = p.no.det.tot[,i-1]
      for(j in 1:(tr_b4_occ[i]-1)){
        temp = t.mat[,,i-1] %*% temp   # this assumes same transition for each day without observation
        ### MAY WANT TO BUILD IN FLEXIBILITY INTO t.mat TO ACCOMODATE DIFFERENT TRANSITIONS HERE ###
      }
      p.no.det.tot[,i] = temp * p.no.det.occ[,i]
      
    } else{
      p.no.det.tot[,i] = t.mat[,,i-1] %*% p.no.det.tot[,i-1] * p.no.det.occ[,i]
    }
  }
  
  #probability of being detected at least once
  p.dot = 1 - sum(p.no.det.tot[,K])
  
  if(is.na(p.dot) | p.dot < 0){
    #print(paste0("p.dot = ", p.dot, " < 0; par = "))
    #print(par)
    p.dot = 0
  }
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  
  ### Calculate sum(log-likelihood of individuals encounter histories) ###
  ind.likelihood = C_indLklhd_RandomWalk(N = N, K = K, Nstates = N.states, tr_b4_occ = tr_b4_occ, 
                                        mu_init = mu.init, y = y, y_pix = y.pix, detDists = detDists, y_platform = y.platform,
                                        y_ds_covar = y.ds.covar, y_g0_covar = y.g0.covar, y_hs_covar = y.hs.covar, y_hb_covar = y.hb.covar,
                                        n_ds_covars = n.ds.covars, n_g0_covars = n.g0.covars, n_hs_covars = n.hs.covars, n_hb_covars = n.hb.covars,
                                        g0s = c(g0s), sigmaDets = c(sigmaDets), Hsigmas = c(hSigma), Hbetas = c(hBeta),
                                        p_no_det_occ = p.no.det.occ, t_mat = t.mat)
  
  #ind.likelihood.r = numeric(N)
  #for(i in 1:N){
  #  ind.likelihood.r[i] = log(sum(mu.init))
  #  state.vec = p.init.state
  #  if(y[i,1] == 1){
  #    state.vec[-y.pix[i,1]] = 0
  #    det.prob = plogis(g0s) * exp(-detDists[i,1]^2 / (2*exp(sigmaDets)^2))
  #    ind.likelihood.r[i] = ind.likelihood.r[i] + log(sum(state.vec * det.prob))
  #  } else{
  #    ind.likelihood.r[i] = ind.likelihood.r[i] + log(sum(state.vec * p.no.det.occ[,1]))
  #  }
  #  state.vec = state.vec / sum(state.vec)
  #  
  #  for(j in 1:N.transitions){
  #    state.vec = t.mat[,,j] %*% state.vec
  #    
  #    if(y[i,j+1] == 1){
  #      state.vec[-y.pix[i,j+1]] = 0
  #      det.prob = plogis(g0s) * exp(-detDists[i,j+1]^2 / (2*exp(sigmaDets)^2))
  #      ind.likelihood.r[i] = ind.likelihood.r[i] + log(sum(state.vec * det.prob))
  #    } else{
  #      ind.likelihood.r[i] = ind.likelihood.r[i] + log(sum(state.vec * p.no.det.occ[,j+1]))
  #    }
  #    state.vec = state.vec / sum(state.vec)
  #  }
  #}

  
  ### DIFFERENCE FROM SCR BOOK FORMULATION ###
  #log.Likelihood = sum(logProb.ind) + N * log(D * area.state.space) + (-D * area.state.space * p.star)
  # Glennie likelihood of thinned Poisson point process (thinned by the detection function):
  # exp(-mu.bar)/factorial(number of individuals detected) * product of encounter history likelihoods
  # this is essentially a Poisson pmf: exp(-lambda) * lambda^x/factorial(x)
  # lambda is the rate parameter of the Poisson distribution, equivalent to mu.bar, the rate of the thinned Poisson point process
  # x is the number of detected individuals
  # lambda^x is the product of x encounter histories
  #log.Likelihood = sum(logProb.ind) - mu.bar - log.fact.N
  log.Likelihood = sum(ind.likelihood) - mu.bar - log.fact.N
  
  if(is.na(log.Likelihood)){
    print(par)
  }
  
  return(-log.Likelihood)
}







