#################################
#################################
#### make the detection grid ####
#################################
#################################
#' Make a detection grid; an array that holds information about the distance to the closest transect line and the detection model
#' 
#' @param habitat A matrix with columns named x and y holding the x and y coordinates for each point in the movement model grid
#' @param lines.arr An array containing information on the transect lines surveyed on each occasion
#' @param detGrid.reso A number describing the resolution of the detection grid
#' @param surveyCovars A character vector of the names of survey covariates; needs to correspond to column names in lines.arr
#'
#' @export
makeDetGrid = function(habitat, lines.arr, detGrid.reso, surveyCovars = NULL){
  
  # habitat is state space data frame that has columns for the center x and y coordinates of each pixel, Area of pixel
  # lines.arr is an array of survey line data, will need covariates if there are such line specific covariates
  # detGrid.reso is resolution of detection grid -> resolution of habitat/detGrid.reso should be an integer
  # surveyCovars is a vector of the names of survey covariates corresponding to column names in lines.arr
  
  ###############################################################
  ### NEED TO INCORPORATE POSSIBILITY OF DETECTION COVARIATES ###
  ###############################################################
  
  dat = as.data.frame(habitat)
  
  reso = sqrt(habitat[1,"Area"])
  
  detGrid = array(NA, dim = c(nrow(dat) * (reso/detGrid.reso)^2, 12 + length(surveyCovars), dim(lines.arr)[3]))
  dimnames(detGrid)[[2]] = c("datCell", "x", "y", "xlim1", "xlim2", "ylim1", "ylim2", "closestLine", "dist", "lt_eps_epe", "hn_haz", "allSame", surveyCovars)
  
  surveyCovarCols = match(surveyCovars, colnames(lines.arr[,,1]))
  
  reso.scale = reso
  det.reso = detGrid.reso
  
  range = 1:((reso/detGrid.reso)^2)
  
  numRange = length(range)
  
  for(i in 1:nrow(dat)){  # loop over larger state space cells
    
    # center of detection grid cells within dat cell i
    xs = seq(dat[i,]$x - reso.scale/2 + det.reso/2, dat[i,]$x + reso.scale/2 - det.reso/2, by = det.reso)
    ys = seq(dat[i,]$y - reso.scale/2 + det.reso/2, dat[i,]$y + reso.scale/2 - det.reso/2, by = det.reso)
    xs.s = rep(xs, reso.scale/det.reso)
    ys.s = rep(ys, reso.scale/det.reso)
    ys.s = ys.s[order(ys.s)]
    
    for(j in 1:dim(lines.arr)[3]){    # loop over occasions
      
      detGrid[range,"datCell",j] = i
      detGrid[range,"x",j] = xs.s
      detGrid[range,"y",j] = ys.s
      detGrid[range,"xlim1",j] = xs.s - det.reso/2
      detGrid[range,"xlim2",j] = xs.s + det.reso/2
      detGrid[range,"ylim1",j] = ys.s - det.reso/2
      detGrid[range,"ylim2",j] = ys.s + det.reso/2
      
      # distances to lines
      lines = !is.na(lines.arr[,1,j])
      for(k in range){
        
        XonLine = (detGrid[k,"x",j] + lines.arr[lines,"slope",j] * detGrid[k,"y",j] - lines.arr[lines,"slope",j] * lines.arr[lines,"intercept",j]) / (lines.arr[lines,"slope",j]^2 + 1)
        YonLine = lines.arr[lines,"slope",j] * XonLine + lines.arr[lines,"intercept",j]
        
        LT_EP = (XonLine >= lines.arr[lines,"x.start",j] & XonLine <= lines.arr[lines,"x.end",j]) | (XonLine <= lines.arr[lines,"x.start",j] & XonLine >= lines.arr[lines,"x.end",j])
        
        dist.lt = sqrt((XonLine - detGrid[k,"x",j])^2 + (YonLine - detGrid[k,"y",j])^2)
        dist.ep.strt = sqrt((lines.arr[lines,"x.start",j] - detGrid[k,"x",j])^2 + (lines.arr[lines,"y.start",j] - detGrid[k,"y",j])^2)
        dist.ep.end = sqrt((lines.arr[lines,"x.end",j] - detGrid[k,"x",j])^2 + (lines.arr[lines,"y.end",j] - detGrid[k,"y",j])^2)
        
        EP = sapply(1:sum(lines), function(x) which(c(dist.ep.strt[x], dist.ep.end[x]) == min(c(dist.ep.strt[x], dist.ep.end[x]))))
        dist.ep = ifelse(EP == 1, dist.ep.strt, dist.ep.end)
        
        dist = ifelse(LT_EP, dist.lt, dist.ep)
        
        detGrid[k,"dist",j] = min(dist)
        detGrid[k,"closestLine",j] = which(dist == min(dist))
        detGrid[k,"lt_eps_epe",j] = ifelse(LT_EP[detGrid[k,"closestLine",j]], 1, ifelse(EP[detGrid[k,"closestLine",j]] == 1, 2, 3))  # 1 = lt, 2 = ep start, 3 = ep end
        ### not filling in hn_haz because everything will be treated as hazard without other info about platform in lines.arr
        
        if(length(surveyCovars) > 0){
          for(z in surveyCovars){
            detGrid[k,z,j] = lines.arr[detGrid[k,"closestLine",j],z,j]
          }
        }
        
        detGrid[k,"hn_haz",j] = lines.arr[detGrid[k,"closestLine",j],"detMod",j]
        
      }
      
      if(all(detGrid[range,"closestLine",j] == detGrid[range[1],"closestLine",j]) & all(detGrid[range,"lt_eps_epe",j] == detGrid[range[1],"lt_eps_epe",j])){
        detGrid[range,"allSame",j] = 1
      } else{
        detGrid[range,"allSame",j] = 0
      }
      
    }
    
    range = range + numRange
    
  }
  
  detGrid[,"hn_haz",][is.na(detGrid[,"hn_haz",])] = 0   # half-normal if anything is left without specification
  
  return(detGrid)
  
}


###############################################
###############################################
#### calculate mean maximum distance moved ####
###############################################
###############################################
#' Calculate the mean of the maximum pariwise distance between each individuals' observations
#' 
#' @param y.pix A matrix; NA when individual -row- not observed on -occasion-; cell ID where observed otherwise
#' @param distMat A matrix; holds pairwise distances between all grid cells in the habitat matrix 
#'
#' @export
meanMaxMove = function(y.pix, distMat){
  maxMoves = numeric(nrow(y.pix))
  for(i in 1:nrow(y.pix)){
    pixs = y.pix[i,][!is.na(y.pix[i,])]
    maxMoves[i] = max(distMat[pixs,pixs])
  }
  return(mean(maxMoves[maxMoves > 0]))
}


###########################################
###########################################
### Activity Center Abundance Estimates ###
###########################################
###########################################
#' Calculate the expected and realized abundance estimate from the closed population, activity center, spatial capture-recapture model
#' 
#' @param data Named list; output from `scr_lt_Closed_SimpleRW()` with MOVEMENT.MODEL = "AC" 
#' @param results Named list; output from `scr_lt_Closed_SimpleRW()` with MOVEMENT.MODEL = "AC" 
#'
#' @export
SCR_lt_Closed_AC_Abundance = function(data, results){  
  
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
  
  par = results$par
  
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
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  # expected activity center abundance in each cell
  mu = D.x
  
  # expected number of unobserved individuals
  EN.unobs = sum(D.x * p.no.det.given.ac)
  
  VCV = solve(results$hessian)
  
  if(n.d.covars == 0){
    
    # Expected abundance and 95% CI
    E.Abund = exp(par[d.init.indices]) * sum(habitat[,"Area"])
    var.E.Abund = sum(exp(par[data$d.init.indices]) * habitat[,"Area"])^2 * VCV[d.init.indices,d.init.indices]
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  } else if(n.d.covars > 0){
    
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[x]]))
    E.Abund.per.site <- exp(ut.D) * habitat[,"Area"]
    E.Abund <- sum(E.Abund.per.site)
    
    gradient = numeric(n.d.covars + 1)
    gradient[1] = E.Abund
    gradient[2:length(gradient)] = sapply(1:n.d.covars, function(x) sum(E.Abund.per.site * habitat[,d.covar.cols[x]]))
    
    var.E.Abund = t(gradient) %*% VCV[d.init.indices, d.init.indices] %*% gradient
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  }
  
  var.EN.unobs = var.E.Abund - E.Abund
  C = exp(1.96 * sqrt(log(1 + (var.EN.unobs / (EN.unobs^2)))))
  R.Abund = N + EN.unobs
  R.Abund.95.CI = c(N + EN.unobs/C, N + EN.unobs*C)
  
  output = data.frame(Estimate = c(E.Abund, R.Abund), lcb = c(E.Abund.95.CI[1], R.Abund.95.CI[1]), ucb = c(E.Abund.95.CI[2], R.Abund.95.CI[2]))
  rownames(output) = c("Expected Abundance", "Realized Abundance")
  
  return(output)
}


#######################################
#######################################
### Random Walk Abundance Estimates ###
#######################################
#######################################
#' Calculate the expected and realized abundance estimate from the closed population, simple random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW" 
#' @param results Named list; output from scr_lt_Closed_SimpleRW with MOVEMENT.MODEL = "RW" 
#'
#' @export
SCR_lt_Closed_SimpleRW_Abundance = function(data, results){  
  
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
  
  par = results$par
  
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
    p.no.det.tot[,1,1] = t.mat[,,1] %*% p.init.state * p.no.det.given.ac[,1]
    p.no.det.tot[,1,2] = t.mat[,,1] %*% D.x * p.no.det.given.ac[,1]
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
  
  # expected activity center abundance in each cell
  mu = D.x
  
  VCV = solve(results$hessian)
  
  if(n.d.covars == 0){
    
    # Expected abundance and 95% CI
    E.Abund = exp(par[d.init.indices]) * sum(habitat[,"Area"])
    var.E.Abund = sum(exp(par[data$d.init.indices]) * habitat[,"Area"])^2 * VCV[d.init.indices,d.init.indices]
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  } else if(n.d.covars > 0){
    
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[x]]))
    E.Abund.per.site <- exp(ut.D) * habitat[,"Area"]
    E.Abund <- sum(E.Abund.per.site)
    
    gradient = numeric(n.d.covars + 1)
    gradient[1] = E.Abund
    gradient[2:length(gradient)] = sapply(1:n.d.covars, function(x) sum(E.Abund.per.site * habitat[,d.covar.cols[x]]))
    
    var.E.Abund = t(gradient) %*% VCV[d.init.indices, d.init.indices] %*% gradient
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  }
  
  var.EN.unobs = var.E.Abund - E.Abund
  C = exp(1.96 * sqrt(log(1 + (var.EN.unobs / (EN.unobs^2)))))
  R.Abund = N + EN.unobs
  R.Abund.95.CI = c(N + EN.unobs/C, N + EN.unobs*C)
  
  output = data.frame(Estimate = c(E.Abund, R.Abund), lcb = c(E.Abund.95.CI[1], R.Abund.95.CI[1]), ucb = c(E.Abund.95.CI[2], R.Abund.95.CI[2]))
  rownames(output) = c("Expected Abundance", "Realized Abundance")
  
  return(output)
  
}


########################################
########################################
### Movement and Detection Estimates ###
########################################
########################################
#' Calculate the parameter estimates from a closed population, activity center or simple random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_SimpleRW
#' @param results Named list; output from scr_lt_Closed_SimpleRW
#'
#' @export
Parameter_Estimates_SimpleRW = function(data, results){
  
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
  
  par = results$par
  
  VCV = solve(results$hessian)
  
  ut.move = par[m.init.indices]
  ut.move.var = VCV[m.init.indices, m.init.indices]
  #delta.move.var = exp(ut.move)^2*ut.move.var
  #c(-1.96, 1.96)*sqrt(delta.move.var) + exp(ut.move)
  
  move.ci = exp(c(-1.96,1.96)*sqrt(ut.move.var) + ut.move)
  
  ut.det.sd = par[ds.init.indices]
  ut.det.sd.var = VCV[ds.init.indices, ds.init.indices]
  det.sd.ci = exp(c(-1.96,1.96)*sqrt(ut.det.sd.var) + ut.det.sd)
  
  ut.g0 = par[g0.init.indices]
  ut.g0.var = VCV[g0.init.indices, g0.init.indices]
  g0.ci = plogis(c(-1.96,1.96)*sqrt(ut.g0.var) + ut.g0)
  
  
  output = data.frame(Estimate = c(exp(ut.move), plogis(ut.g0), exp(ut.det.sd)),
                      lcb = c(move.ci[1], g0.ci[1], det.sd.ci[1]),
                      ucb = c(move.ci[2], g0.ci[2], det.sd.ci[2]))
  rownames(output) = c("Movement Scale", "g0", "Detection Scale")
  
  return(output)
  
}


#######################################
#######################################
### Random Walk Abundance Estimates ###
#######################################
#######################################
#' Calculate the expected and realized abundance estimate from the closed population, biased random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_BiasedRW
#' @param results Named list; output from scr_lt_Closed_BiasedRW
#'
#' @export
SCR_lt_Closed_BiasedRW_Abundance = function(data, results){  
  
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
  
  par = results$par
  
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
  
  
  EN.unobs = D.bar * (1 - p.dot)
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  # expected activity center abundance in each cell
  mu = D.x
  
  VCV = solve(results$hessian)
  
  if(n.d.covars == 0){
    
    # Expected abundance and 95% CI
    E.Abund = exp(par[d.init.indices]) * sum(habitat[,"Area"])
    var.E.Abund = sum(exp(par[data$d.init.indices]) * habitat[,"Area"])^2 * VCV[d.init.indices,d.init.indices]
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  } else if(n.d.covars > 0){
    
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[x]]))
    E.Abund.per.site <- exp(ut.D) * habitat[,"Area"]
    E.Abund <- sum(E.Abund.per.site)
    
    gradient = numeric(n.d.covars + 1)
    gradient[1] = E.Abund
    gradient[2:length(gradient)] = sapply(1:n.d.covars, function(x) sum(E.Abund.per.site * habitat[,d.covar.cols[x]]))
    
    var.E.Abund = t(gradient) %*% VCV[d.init.indices, d.init.indices] %*% gradient
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  }
  
  var.EN.unobs = var.E.Abund - E.Abund
  C = exp(1.96 * sqrt(log(1 + (var.EN.unobs / (EN.unobs^2)))))
  R.Abund = N + EN.unobs
  R.Abund.95.CI = c(N + EN.unobs/C, N + EN.unobs*C)
  
  output = data.frame(Estimate = c(E.Abund, R.Abund), lcb = c(E.Abund.95.CI[1], R.Abund.95.CI[1]), ucb = c(E.Abund.95.CI[2], R.Abund.95.CI[2]))
  rownames(output) = c("Expected Abundance", "Realized Abundance")
  
  return(output)
  
}


########################################
########################################
### Movement and Detection Estimates ###
########################################
########################################
#' Calculate the parameter estimates from a closed population, biased random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_BiasedRW
#' @param results Named list; output from scr_lt_Closed_BiasedRW
#'
#' @export
Parameter_Estimates_BiasedRW = function(data, results){
  
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
  
  par = results$par
  
  VCV = solve(results$hessian)
  
  ut.move = par[m.init.indices]
  ut.move.var = VCV[m.init.indices, m.init.indices]
  #delta.move.var = exp(ut.move)^2*ut.move.var
  #c(-1.96, 1.96)*sqrt(delta.move.var) + exp(ut.move)
  
  move.ci = exp(c(-1.96,1.96)*sqrt(ut.move.var) + ut.move)
  
  ut.det.sd = par[ds.init.indices]
  ut.det.sd.var = VCV[ds.init.indices, ds.init.indices]
  det.sd.ci = exp(c(-1.96,1.96)*sqrt(ut.det.sd.var) + ut.det.sd)
  
  ut.g0 = par[g0.init.indices]
  ut.g0.var = VCV[g0.init.indices, g0.init.indices]
  g0.ci = plogis(c(-1.96,1.96)*sqrt(ut.g0.var) + ut.g0)
  
  ut.rho = par[rho.init.indices]
  ut.rho.var = VCV[rho.init.indices, rho.init.indices]
  rho.ci = plogis(c(-1.96,1.96)*sqrt(ut.rho.var) + ut.rho)
  
  
  output = data.frame(Estimate = c(exp(ut.move), plogis(ut.g0), exp(ut.det.sd), plogis(ut.rho)),
                      lcb = c(move.ci[1], g0.ci[1], det.sd.ci[1], rho.ci[1]),
                      ucb = c(move.ci[2], g0.ci[2], det.sd.ci[2], rho.ci[2]))
  rownames(output) = c("Movement Scale", "g0", "Detection Scale", "rho")
  
  return(output)
  
}



##################################################
##################################################
### Correlated Random Walk Abundance Estimates ###
##################################################
##################################################
#' Calculate the expected and realized abundance estimate from the closed population, correlated random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_CorrelatedRW
#' @param results Named list; output from scr_lt_Closed_CorrelatedRW
#'
#' @export
SCR_lt_Closed_CorrelatedRW_Abundance = function(data, results){  
  
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
  
  par = results$par
  
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
  
  
  # expected number of unobserved individuals
  EN.unobs = D.bar * (1 - p.dot)
  
  # thinned Poisson point process rate parameter
  mu.bar = p.dot * D.bar
  
  # expected activity center abundance in each cell
  mu = D.x
  
  VCV = solve(results$hessian)
  
  if(n.d.covars == 0){
    
    # Expected abundance and 95% CI
    E.Abund = exp(par[d.init.indices]) * sum(habitat[,"Area"])
    var.E.Abund = sum(exp(par[data$d.init.indices]) * habitat[,"Area"])^2 * VCV[d.init.indices,d.init.indices]
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  } else if(n.d.covars > 0){
    
    ut.D <- betas[1] + rowSums(sapply(1:n.d.covars, function(x) betas[x+1] * habitat[,d.covar.cols[x]]))
    E.Abund.per.site <- exp(ut.D) * habitat[,"Area"]
    E.Abund <- sum(E.Abund.per.site)
    
    gradient = numeric(n.d.covars + 1)
    gradient[1] = E.Abund
    gradient[2:length(gradient)] = sapply(1:n.d.covars, function(x) sum(E.Abund.per.site * habitat[,d.covar.cols[x]]))
    
    var.E.Abund = t(gradient) %*% VCV[d.init.indices, d.init.indices] %*% gradient
    C.prime = exp(1.96 * sqrt(log(1 + ( var.E.Abund / E.Abund^2 ))))
    E.Abund.95.CI = c(E.Abund/C.prime, E.Abund*C.prime)
    
  }
  
  var.EN.unobs = var.E.Abund - E.Abund
  C = exp(1.96 * sqrt(log(1 + (var.EN.unobs / (EN.unobs^2)))))
  R.Abund = N + EN.unobs
  R.Abund.95.CI = c(N + EN.unobs/C, N + EN.unobs*C)
  
  output = data.frame(Estimate = c(E.Abund, R.Abund), lcb = c(E.Abund.95.CI[1], R.Abund.95.CI[1]), ucb = c(E.Abund.95.CI[2], R.Abund.95.CI[2]))
  rownames(output) = c("Expected Abundance", "Realized Abundance")
  
  return(output)
}






########################################
########################################
### Movement and Detection Estimates ###
########################################
########################################
#' Calculate the parameter estimates from a closed population, correlated random walk, spatial capture-recapture model
#' 
#' @param data Named list; output from scr_lt_Closed_CorrelatedRW
#' @param results Named list; output from scr_lt_Closed_CorrelatedRW
#'
#' @export
Parameter_Estimates_CorrelatedRW = function(data, results){
  
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
  
  par = results$par
  
  VCV = solve(results$hessian)
  
  ut.move = par[m.init.indices]
  ut.move.var = VCV[m.init.indices, m.init.indices]
  #delta.move.var = exp(ut.move)^2*ut.move.var
  #c(-1.96, 1.96)*sqrt(delta.move.var) + exp(ut.move)
  
  move.ci = exp(c(-1.96,1.96)*sqrt(ut.move.var) + ut.move)
  
  ut.det.sd = par[ds.init.indices]
  ut.det.sd.var = VCV[ds.init.indices, ds.init.indices]
  det.sd.ci = exp(c(-1.96,1.96)*sqrt(ut.det.sd.var) + ut.det.sd)
  
  ut.g0 = par[g0.init.indices]
  ut.g0.var = VCV[g0.init.indices, g0.init.indices]
  g0.ci = plogis(c(-1.96,1.96)*sqrt(ut.g0.var) + ut.g0)
  
  ut.gamma = par[gamma.init.indices]
  ut.gamma.var = VCV[gamma.init.indices, gamma.init.indices]
  gamma.ci = plogis(c(-1.96,1.96)*sqrt(ut.gamma.var) + ut.gamma)
  
  ut.beta = par[beta.init.indices]
  ut.beta.var = VCV[beta.init.indices, beta.init.indices]
  beta.ci = 2 * pi * plogis(c(-1.96,1.96)*sqrt(ut.beta.var) + ut.beta) - pi
  
  
  output = data.frame(Estimate = c(exp(ut.move), plogis(ut.g0), exp(ut.det.sd), plogis(ut.gamma), 2*pi*plogis(ut.beta) - pi),
                      lcb = c(move.ci[1], g0.ci[1], det.sd.ci[1], gamma.ci[1], beta.ci[1]),
                      ucb = c(move.ci[2], g0.ci[2], det.sd.ci[2], gamma.ci[2], beta.ci[2]))
  rownames(output) = c("Movement Scale", "g0", "Detection Scale", "gamma", "beta")
  
  return(output)
  
}