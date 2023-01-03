############################################################################
## Data simulation function - Location-t ~ Laplace(Location-t-1, move.sd) ##
############################################################################
#' Simulate a spatial capture-recapture data set where individuals move according to a simple random walk
#'  
#' @param K.train A number > 0; The number of sampling occasions in the encounter history that the model will be fit too
#' @param K.forecast A number >= 0; The number of sampling occasions in the encounter history that will be compared to model forecasts
#' @param det.sd A number > 0; The scale parameter of the half-normal detection function
#' @param p0 A number > 0 and <= 1; Detection probability for an individual on a transect line
#' @param move.sd A number > 0; The scale parameter for the movement model
#' @param moveDist A character string, either "Laplace" or "Normal"; The distribution used to simulate movement
#' @param A.scale A number > 0; The scale of the grid the model uses to compute the movement of individuals
#' @param detgrid.scale A number > 0; The scale of the grid the model uses to compute detection probability
#' @param makeDetGrid Logical
#'
#' @export
simData_SRW <- function(K.train, K.forecast, det.sd, p0, move.sd, moveDist = "Laplace", A.scale = 2, detgrid.scale = 0.5, makeDetGrid = T){
  
  ########################
  ### STATE SPACE DATA ###
  ########################
  
  K = K.train + K.forecast
  
  # survey limits
  xlim_surveys = c(10,30)
  ylim_surveys = c(11,29)
  
  # state space limits for SCR
  xlim = c(0,40)
  ylim = c(0,40)
  
  Area = (xlim[2]-xlim[1])*(ylim[2]-ylim[1])
  
  # northings and east/west endpoints for tracklines
  # 10 tracklines, 20 units long
  northings = seq(ylim_surveys[1],ylim_surveys[2],by=2)
  ew_bounds = c(xlim_surveys[1],xlim_surveys[2])
  
  # habitat is a grid with grid cell dimension of CellDim X CellDim units
  CellDim = 1
  ycoords = seq(ylim[1]+CellDim/2, ylim[2]-CellDim/2, by=CellDim)
  xcoords = ycoords
  habitat = matrix(ncol = 3, nrow = length(ycoords)^2)
  colnames(habitat) = c("x", "y", "Area")
  
  habitat[,1] = rep(xcoords,length(ycoords))
  habitat[,2] = rep(ycoords, length(xcoords))[order(rep(ycoords, length(xcoords)))]
  habitat[,3] = CellDim^2
  
  distMat = matrix(nrow = nrow(habitat), ncol = nrow(habitat))
  
  for(i in 1:nrow(habitat)){
    distMat[i,] = sqrt((habitat[i,1] - habitat[,1])^2 + (habitat[i,2] - habitat[,2])^2)
  }

  
  #########################
  ### s AND u LOCATIONS ###
  #########################
  
  # population size to simulate
  N = 100   # 0.0625 individuals per unit area

  #actual location on each occasion
  Ux <- matrix(NA, N, K)
  Uy <- matrix(NA, N, K)
  
  if(moveDist == "Laplace"){
    for(i in 1:N){
      UPix = which(rmultinom(1, 1, rep(1, nrow(habitat))) == 1)
      Ux[i,1] = runif(1, habitat[UPix, "x"] - CellDim/2, habitat[UPix, "x"] + CellDim/2)
      Uy[i,1] = runif(1, habitat[UPix, "y"] - CellDim/2, habitat[UPix, "y"] + CellDim/2)
      
      for(k in 2:K){
        cellProbs.pre = (1/(2*move.sd)) * exp(-distMat[UPix,]/move.sd)
        cellProbs.std = cellProbs.pre / sum(cellProbs.pre)
        UPix = which(rmultinom(1, 1, cellProbs.std) == 1)
        Ux[i,k] = runif(1, habitat[UPix, "x"] - CellDim/2, habitat[UPix, "x"] + CellDim/2)
        Uy[i,k] = runif(1, habitat[UPix, "y"] - CellDim/2, habitat[UPix, "y"] + CellDim/2)
      }
    }
  } else if(moveDist == "Normal"){
    for(i in 1:N){
      UPix = which(rmultinom(1, 1, rep(1, nrow(habitat))) == 1)
      Ux[i,1] = runif(1, habitat[UPix, "x"] - CellDim/2, habitat[UPix, "x"] + CellDim/2)
      Uy[i,1] = runif(1, habitat[UPix, "y"] - CellDim/2, habitat[UPix, "y"] + CellDim/2)
      
      for(k in 2:K){
        cellProbs.pre = (1/(2*move.sd*sqrt(2*pi)) * exp(-0.5 * (distMat[UPix,]/move.sd)^2))
        cellProbs.std = cellProbs.pre / sum(cellProbs.pre)
        UPix = which(rmultinom(1, 1, cellProbs.std) == 1)
        Ux[i,k] = runif(1, habitat[UPix, "x"] - CellDim/2, habitat[UPix, "x"] + CellDim/2)
        Uy[i,k] = runif(1, habitat[UPix, "y"] - CellDim/2, habitat[UPix, "y"] + CellDim/2)
      }
    }
  }
  
  #plot(Uy[1,] ~ Ux[1,], ylim = ylim, xlim = xlim)
  #text(Uy[1,] ~ Ux[1,], labels = 1:10)
  #for(i in 2:5){
  #  points(Ux[i,], Uy[i,], col = i)
  #  text(Uy[i,] ~ Ux[i,], col = i)
  #}
  
  #plot(Ux,Uy)

  #############################
  ### DETECTION PROBABILITY ###
  #############################
  
  # treat surveys as contnuous horizontal lines between eastern and western most points
  min.d = matrix(nrow=N, ncol=K)
  ytr = seq(ylim_surveys[1],ylim_surveys[2],by=2)
  northings = ytr
  ew_bounds = c(xlim_surveys[1],xlim_surveys[2])
  for(i in 1:N){
    for(k in 1:K){
      if(Ux[i,k] >= xlim_surveys[1] & Ux[i,k] <= xlim_surveys[2]){
        ### if U is within east west boundaries, 
        # min distance is min of distances between survey line northings and U
        min.d[i,k] = min(abs(northings - Uy[i,k]))
      } else {
        ### if U is outside east west boundaries,
        # min distance is distance to from U to end point of line with closest northing
        closest_northing = northings[which(abs(northings - Uy[i,k]) == min(abs(northings - Uy[i,k])))]
        min.d[i,k] = min(sapply(ew_bounds, function(z) sqrt((z-Ux[i,k])^2 + (closest_northing-Uy[i,k])^2)))
      }										
    }
  }
  
  
  # detection probability
  p <- p0 * exp(-0.5 * (min.d/det.sd)^2) #to store detection probs
  
  
  
  #########################
  ### ENCOUNTER HISTORY ###
  #########################
  
  #observations
  y <- matrix(NA, N, K) #to store detections
  y.x <- matrix(NA, N, K) #to store observed locations
  y.y <- matrix(NA, N, K) #to store observed locations
  for(i in 1:N){ # Loop over individuals
    for(k in 1:K){ # Loop over occasions
      y[i,k] <- rbinom(1, 1, p[i,k])
      if (y[i,k]==1) { #if detected, save locations
        y.x[i,k] <- Ux[i,k]
        y.y[i,k] <- Uy[i,k]
      } #if
    } #k
  } #i
  
  ##############################
  ### FORMAT DATA FOR MODELS ###
  ##############################
  
  #remove individuals never detected
  captured.test <- which(rowSums(y[,1:K.train])>0)
  y.test <- y[captured.test,1:K.train]
  y.x.test <- y.x[captured.test,1:K.train]
  y.y.test <- y.y[captured.test,1:K.train]
  detDists.test <- min.d[captured.test,1:K.train]
  
  captured <- which(rowSums(y) > 0)
  newly.captured <- captured[!(captured %in% captured.test)]
  y.forecast <- rbind(y[captured.test, (K.train+1):(K.train+K.forecast)],
                      y[newly.captured, (K.train+1):(K.train+K.forecast)])
  y.x.forecast <- rbind(y.x[captured.test,(K.train+1):(K.train+K.forecast)],
                        y.x[newly.captured, (K.train+1):(K.train+K.forecast)])
  y.y.forecast <- rbind(y.y[captured.test,(K.train+1):(K.train+K.forecast)],
                        y.y[newly.captured, (K.train+1):(K.train+K.forecast)])
  detDists.forecast <- rbind(min.d[captured.test,(K.train+1):(K.train+K.forecast)],
                             min.d[newly.captured, (K.train+1):(K.train+K.forecast)])
  
  
  # analysis movement grid #
  A.ylim = ylim    #c(0,40)
  A.xlim = xlim    #c(0,40)
  A.ycoords = seq(A.ylim[1]+A.scale/2, A.ylim[2]-A.scale/2, A.scale)
  A.xcoords = A.ycoords
  
  A.grid = matrix(nrow = length(A.ycoords)^2, ncol = 3)
  colnames(A.grid) = c("x","y","Area")
  A.grid[,1] = rep(A.xcoords, length(A.ycoords))
  A.grid[,2] = rep(A.ycoords, length(A.xcoords))[order(rep(A.ycoords, length(A.xcoords)))]
  A.grid[,3] = A.scale^2
  
  A.distMat = matrix(nrow = nrow(A.grid), ncol = nrow(A.grid))
  
  for(i in 1:nrow(A.grid)){
    A.distMat[i,] = sqrt((A.grid[i,1] - A.grid[,1])^2 + (A.grid[i,2] - A.grid[,2])^2)
  }
  
  A.angleMat = matrix(nrow = nrow(A.grid), ncol = nrow(A.grid))
  
  for(i in 1:nrow(A.grid)){
    A.angleMat[i,] = atan((A.grid[,2]-A.grid[i,2]) / (A.grid[,1]-A.grid[i,1]))
  }
  
  #A.distCube = array(dim = c(nrow(A.grid), nrow(A.grid), nrow(A.grid)))
  #for(i in 1:nrow(A.grid)){
  #  for(j in 1:nrow(A.grid)){
  #    A.distCube[i,,j] = sqrt(( (A.grid[i,1] + A.grid[j,1])/2 - A.grid[,1] )^2 + ( (A.grid[i,2] + A.grid[j,2])/2 - A.grid[,2] )^2)
  #  }
  #}
  
  y.pix.test = matrix(0, nrow = nrow(y.test), ncol = ncol(y.test))
  for(i in 1:nrow(y.pix.test)){
    for(j in 1:ncol(y.pix.test)){
      if(y.test[i,j] == 1){
        dists = sqrt((y.x.test[i,j] - A.grid[,"x"])^2 + (y.y.test[i,j] - A.grid[,"y"])^2)
        y.pix.test[i,j] = which(dists == min(dists))
      }
    }
  }
  
  y.pix.forecast = matrix(0, nrow = nrow(y.forecast), ncol = ncol(y.forecast))
  for(i in 1:nrow(y.pix.forecast)){
    for(j in 1:ncol(y.pix.forecast)){
      if(y.forecast[i,j] == 1){
        dists = sqrt((y.x.forecast[i,j] - A.grid[,"x"])^2 + (y.y.forecast[i,j] - A.grid[,"y"])^2)
        y.pix.forecast[i,j] = which(dists == min(dists))
      }
    }
  }
  
  y.pix.all = matrix(0, nrow = nrow(y), ncol = ncol(y))
  for(i in 1:nrow(y)){
    for(j in 1:ncol(y)){
        dists = sqrt((Ux[i,j] - A.grid[,"x"])^2 + (Uy[i,j] - A.grid[,"y"])^2)
        y.pix.all[i,j] = which(dists == min(dists))
    }
  }
  
  y.pix.all.obs = y.pix.all[captured.test,]
  y.pix.all.unobs = y.pix.all[-captured.test,]
  
  # survey line data
  tlines <- matrix(NA, length(ytr), 7) 
  colnames(tlines) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept", "detMod")
  for(n in 1:length(ytr)){ # Loop over transect lines
    tlines[n,1:4] <- c(xlim_surveys[1], xlim_surveys[2], ytr[n], ytr[n]) #x.start, x.end, y.start, y.end
    tlines[n,5] <- (tlines[n,4] - tlines[n,3]) / (tlines[n,2] - tlines[n,1]) #slope = rise/run
    tlines[n,6] <- tlines[n,3] - (tlines[n,1]*tlines[n,5]) #intercept = y.start -x.start*slope
    tlines[n,7] <- 0
  }
  lines.arr.test <- array(tlines, dim=c(dim(tlines), K.train))
  colnames(lines.arr.test) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept", "detMod")
  lines.arr.forecast <- array(tlines, dim=c(dim(tlines), K.forecast))
  colnames(lines.arr.forecast) = c("x.start", "x.end", "y.start", "y.end", "slope", "intercept", "detMod")
  
  if(makeDetGrid){
    detGrid.test = makeDetGrid(habitat = A.grid, lines.arr = lines.arr.test, detGrid.reso = detgrid.scale)
    detGrid.forecast = makeDetGrid(habitat = A.grid, lines.arr = lines.arr.forecast, detGrid.reso = detgrid.scale)
  }
  
  mean.ac.move = meanMaxMove(y.pix = y.pix.test, distMat = A.distMat)
  
  y.g0.covar.test = array(y.test, dim = c(nrow(y.test), ncol(y.test), 1))
  y.g0.covar.test[y.g0.covar.test == 0] = NA
  y.g0.covar.test[y.g0.covar.test == 1] = 0
  y.ds.covar.test = y.g0.covar.test
  y.hs.covar.test = y.g0.covar.test
  y.hb.covar.test = y.g0.covar.test
  y.platform.test = y.g0.covar.test[,,1]
  
  y.g0.covar.forecast = array(y.forecast, dim = c(nrow(y.forecast), ncol(y.forecast), 1))
  y.g0.covar.forecast[y.g0.covar.forecast == 0] = NA
  y.g0.covar.forecast[y.g0.covar.forecast == 1] = 0
  y.ds.covar.forecast = y.g0.covar.forecast
  y.hs.covar.forecast = y.g0.covar.forecast
  y.hb.covar.forecast = y.g0.covar.forecast
  y.platform.forecast = y.g0.covar.forecast[,,1]
  
  if(makeDetGrid){
    out = list(y = y.test, 
               y.pix = y.pix.test,
               y.g0.covar = y.g0.covar.test,
               y.ds.covar = y.ds.covar.test,
               y.hs.covar = y.hs.covar.test,
               y.hb.covar = y.hb.covar.test,
               y.platform = y.platform.test,
               detDists = detDists.test, 
               habitat = A.grid,
               detGrid = detGrid.test,
               distMat = A.distMat,
               angleMat = A.angleMat,
               lines.arr = lines.arr.test,
               Area = A.scale^2,
               detArea = detgrid.scale^2,
               truncDist = 5,
               nCells = nrow(A.grid),
               nDetCells = nrow(detGrid.test),
               N = nrow(y.test),
               K = K.train,
               N.states = nrow(A.grid),
               N.transitions = K.train - 1,
               tr_b4_occ = c(0,rep(1,K.train-1)),
               log.fact.N = sum(log(1:nrow(y.test))),
               mean.ac.move = mean.ac.move,
               y.forecast = y.forecast, 
               y.pix.forecast = y.pix.forecast,
               y.g0.covar.forecast = y.g0.covar.forecast,
               y.ds.covar.forecast = y.ds.covar.forecast,
               y.hs.covar.forecast = y.hs.covar.forecast,
               y.hb.covar.forecast = y.hb.covar.forecast,
               y.platform.forecast = y.platform.forecast,
               detDists.forecast = detDists.forecast,
               detGrid.forecast = detGrid.forecast,
               lines.arr.forecast = lines.arr.forecast,
               N.forecast = nrow(y.forecast),
               K.forecast = K.forecast,
               tr_b4_occ.forecast = c(0,rep(1,K.forecast-1)),
               y.pix.all.obs = y.pix.all.obs,
               y.pix.all.unobs = y.pix.all.unobs)
  } else{
    out = list(y = y.test, 
               y.pix = y.pix.test,
               y.g0.covar = y.g0.covar.test,
               y.ds.covar = y.ds.covar.test,
               y.hs.covar = y.hs.covar.test,
               y.hb.covar = y.hb.covar.test,
               y.platform = y.platform.test,
               detDists = detDists.test, 
               habitat = A.grid,
               #detGrid = detGrid.test,
               distMat = A.distMat,
               angleMat = A.angleMat,
               lines.arr = lines.arr.test,
               Area = A.scale^2,
               detArea = detgrid.scale^2,
               truncDist = 5,
               nCells = nrow(A.grid),
               #nDetCells = nrow(detGrid),
               N = nrow(y.test),
               K = K.train,
               N.states = nrow(A.grid),
               N.transitions = K.train - 1,
               tr_b4_occ = c(0,rep(1,K.train-1)),
               log.fact.N = sum(log(1:nrow(y.test))),
               mean.ac.move = mean.ac.move,
               y.forecast = y.forecast, 
               y.pix.forecast = y.pix.forecast,
               y.g0.covar.forecast = y.g0.covar.forecast,
               y.ds.covar.forecast = y.ds.covar.forecast,
               y.hs.covar.forecast = y.hs.covar.forecast,
               y.hb.covar.forecast = y.hb.covar.forecast,
               y.platform.forecast = y.platform.forecast,
               detDists.forecast = detDists.forecast,
               #detGrid.forecast = detGrid.forecast,
               lines.arr.forecast = lines.arr.forecast,
               N.forecast = nrow(y.forecast),
               K.forecast = K.forecast,
               tr_b4_occ.forecast = c(0,rep(1,K.forecast-1)),
               y.pix.all.obs = y.pix.all.obs,
               y.pix.all.unobs = y.pix.all.unobs)
  }

  

  
  return(out)
}
