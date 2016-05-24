## Authors 
## Martin Schlather, schlather@math.uni-mannheim.de
##
##
## Copyright (C) 2015 Martin Schlather
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  



RFempiricalvariogram <- function(
 x, y = NULL, z = NULL, T = NULL, data, grid, bin = NULL, 
 phi=NULL,  ## phi, number of anglular segments per PI
 theta = NULL, ## aehnlich
 deltaT = NULL, ##  deltaT[1] max abstand, deltaT[2] : gitterabstand
 distances, vdim, ...
 ) {

  ## repetition is last dimension

  ## bin centers will be a vector of scalar distances  (in cylinder coord, e.g.)
  ## for the angles: start always with the first on negative angle, continue
  ##                 counter clockwise [0, 2pi]
  ## in 3 d the third angle is zero if vector in the (x, y) plane, positive
  ##                 angle starts with points above plane
  
  ## make sure that exactly one negative value appears, and that zero is
  ## added if bin starts with a positive value
#Print("OK");
  stopifnot(length(theta) <= 1, length(phi) <= 1)
  
  if ((is(data, "RFsp") || isSpObj(data)) && !missing(x))
    stop("x, y, z, T may not be given if 'data' is of class 'RFsp' or an 'sp' object")
#Print("OK1"); 

  ## to do: distances
  if (!missing(distances) && length(distances)>0)
    stop("option distances not programmed yet.")
  
  RFoptOld <- internal.rfoptions(...)
  on.exit(RFoptions(LIST=RFoptOld[[1]]))
  RFopt <- RFoptOld[[2]]
  varunits <- RFopt$coords$varunits
  call <- match.call()

  Z <- StandardizeData(x=x, y=y, z=z, T=T, distances=distances, grid=grid,
                       RFopt = RFopt,
                       data=data,  allowFirstCols=FALSE,
                       vdim=if (missing(vdim)) NULL else vdim)
  grid <- sapply(Z$coord, function(z) z$grid)
  fft <- RFopt$empvario$fft && grid[1] && all(grid == grid[1])
  time <- Z$Zeit

#  Print(Z, deltaT)
  if (time && !fft) {
    stop("currently time components are not possible")
    if (grid[1])
      stop(if (all(grid == grid[1])) "Failed for 'fft=FALSE'. Try 'fft=TRUE'."
           else "currently time components are not possible")
    if (Z$spatialdim == 3) stop("currently time components are not possible for spatial dimensions greater than 2")
    warning("time component treated as spatial component")
    time <- Z$Zeit <- FALSE
    Z$xdimOZ <- Z$xdimOZ + 1
    Z$spatialdim <- Z$spatialdim + 1
    rangeT <- NULL
    for (i in 1:length(Z$coord)) {
      T <- Z$coord[[i]]$T
      Z$coord[[i]]$T <- NULL
      T <- seq(T[1], step=T[2], length=T[3])
      rangeT <- range(rangeT, T)
      #Print(T, rep(T, each=nrow(Z$coord[[i]]$x)))
      Z$coord[[i]]$x <- cbind(matrix(ncol=ncol(Z$coord[[i]]$x),
                                     rep(Z$coord[[i]]$x, length(T))),
                              rep(T, each=nrow(Z$coord[[i]]$x)))
      if (length(Z$coord[[i]]$y) > 0)
        Z$coord[[i]]$y <- cbind(matrix(ncol=ncol(Z$coord[[i]]$y),
                                       rep(Z$coord[[i]]$y, length(T))),
                                rep(T, each=nrow(Z$coord[[i]]$y)))
      Z$coord[[i]]$spatialdim <- Z$coord[[i]]$spatialdim + 1
      Z$coord[[i]]$Zeit <- FALSE
    }
    Z$rangex <- cbind(Z$rangex, rangeT)
    ##  ## to do
    ## to do multivariate;

    if (Z$spatialdim == 1 || length(phi) == 0) phi <- 10
    if (Z$spatialdim == 2 && length(theta)==0) theta <- 10
    deltaT <- NULL
    #Print(Z)
    
  }

  if (Z$dist.given) stop("option distances not programmed yet.")
  
  if (missing(vdim) || length(vdim) == 0) {
    vdim <- if (!is.na(Z$vdim)) Z$vdim else 1
  } else {
    if (!is.na(Z$vdim) && vdim!=Z$vdim)
      warning("given multivariate dimension 'vdim' does not match multivariate dimension of the data")
  }

  data <- RFboxcox(Z$data)
   
  
  restotal <- sapply(Z$coord, function(z) z$restotal)
  spatialdim <- Z$spatialdim
 
  len.data <- sapply(data, length)
  repetitions <- as.integer(len.data / (restotal * vdim))  
  if (any(repetitions)==0) stop("no data given")
  if (any(len.data != restotal * vdim * repetitions))
    stop("number of data does not match coordinates")
  sets <- length(Z$data)
#Print("OK2"); 

  for (i in 1:sets) {
    dim.data <- c(restotal[i], vdim, repetitions[i])
    dim(data[[i]]) <- dim.data

    if (vdim > 1 && repetitions[i] > 1) {
      dataX <- aperm(data[[i]], c(1, 3, 2)) ## now: coord, repet, vdim
      dim(dataX) <- c(dim.data[1] * dim.data[3], dim.data[2])
      variance <- cov(dataX)
      rm(dataX)
    } else {
      dim(data[[i]]) <- if (vdim == 1) prod(dim.data) else dim.data[1:2]    
      variance <- var(data[[i]])
      dim(data[[i]]) <- dim.data
    }
  }
 
  if(is.null(bin) || length(bin)==0) bin <- 20

  if (length(bin) == 1) {
    ## automatic bin depending on coords
    xx <- Z$coord[[1]]$x
    if(grid[1])
      bin <- seq(0, max(xx[2, ] * xx[3, ]) / 2, len = bin) 
    else {
      bin <- seq(0, sqrt(sum((apply(xx, 2, max)-apply(xx, 2, min))^2))/2,
                 len = bin)
    }
    if (RFopt$general$printlevel >= PL_SUBIMPORTANT)
      message("Bins in RFempiricalvariogram are chosen automatically:\n", 
              paste(signif(bin, 2), collapse=" "))
  }

  pseudo <- RFopt$empvario$pseudovariogram
  phi0 <- RFopt$empvario$phi0 # 0 if automatic
  theta0 <- RFopt$empvario$theta0 # 0 if automatic

  thetagiven <- length(theta)>0 && spatialdim > 2 && theta > 1
  phigiven <-  length(phi)>0 && spatialdim > 1 && phi > 1 
  deltaTgiven <- length(deltaT)>0 && all(deltaT > 0)
  basic <- !(time || phigiven || thetagiven)
     
  if(time && pseudo)
    stop("Time component is not compatible with Pseudo variogram") # to do
  if(!fft && pseudo) ## to do
    stop("Pseudo variogram is only implemented for grids with the FFT method")
  ## IS THE FFT FLAG SET

  #fft <- fft && repetitions == 1 # to do ! fft should allow for repetitions  
  
  bin <- prepareBin(bin)
  stopifnot(length(bin)>=2, all(is.finite(bin)))
  if (any(diff(bin)<=0)) stop("bin must be a strictly increasing sequence")
   ##  is.null(bin) in fft : see version 3.0.12 or earlier ! to do ?! 
  
  centers <- pmax(0, (bin[-1] + bin[-length(bin)])/2)
  n.bins <- length(bin) - 1

#Print(phi0, phigiven)
#  Print(deltaT)

  
  phi <- if (!phigiven) c(0, 0) else c(phi0, phi)        
  theta <- if (!thetagiven) c(0, 0) else c(theta0, theta)
  if (!deltaTgiven) deltaT <- c(0, 0)
  stopifnot(0 <= phi[1], 2 * pi > phi[1], 
            0 <= theta[1], 2 * pi > theta[1], 
            phi[2] >= 0,  phi[2] == as.integer(phi[2]), 
            theta[2] >= 0, theta[2] == as.integer(theta[2]),
            all(is.finite(deltaT)), all(deltaT >= 0))

  if (time) {
    T.start  <- sapply(Z$coord, function(x) x$T[1])
    T.step <- sapply(Z$coord, function(x) x$T[2])
    T.len  <- sapply(Z$coord, function(x) x$T[3])
    if (sets > 1) {
      if (any(abs(diff(diff(T.step))) > 1e-15))
        stop("only data sets with the same time step allowed") # generalise todo
    }
    T <-  c(0, T.step[1], max(T.len))
  } else {
    T <-  c(1, 1, 1)
  }

  if (length(deltaT) == 1) deltaT <- c(deltaT, 1)
  realdelta <- deltaT[2] * T[2]
  
  NotimeComponent <- !deltaTgiven || T[3]==1
  stepT <-  deltaT[2] / T[2]

  #Print(deltaT, T, stepT, T[2] * (T[3]-1), max(T[2], realdelta) )
  
  if (stepT != as.integer(stepT)) {
    #Print(T, stepT, deltaT)
    stop("deltaT not multiple of distance of temporal grid")
  }
  stepT <- max(1, stepT)
  nstepT <- as.integer(min(deltaT[1], T[2] * (T[3]-1)) / max(T[2], realdelta))

  #Print(nstepT)
  
  n.theta <- max(1, theta[2])
  n.delta <- 1 + nstepT
  n.phi <- max(1, phi[2])
  if (!fft && !basic) {
    n.phibin <- 2 * n.phi
  } else {
    n.phibin <-
      if (!pseudo && NotimeComponent) max(1, n.phi)
      else if (phi[2]==0) 1 else 2 * n.phi
  }

  totalbinsOhnevdim <- as.integer(n.bins * n.phibin * n.theta * n.delta)
  totalbins <- totalbinsOhnevdim * vdim^2

  phibins <- thetabins <- Tbins <- NULL

 # Print(nstepT, realdelta)
  
  if (!NotimeComponent) Tbins <- (0:nstepT) * realdelta
  if (phi[2] > 0) phibins <- phi[1] + 0 : (n.phibin - 1) * pi / n.phi

  if (n.theta > 1)
    thetabins <- theta[1] + (0 : (n.theta-1) + 0.5) * pi / n.theta
 
  dims <- c(bins=n.bins, phi=n.phibin, theta=n.theta, delta=n.delta,
            vdim=rep(vdim, 2))

  emp.vario.sd <- NULL
#Print("OK4", fft, vdim, basic, thetagiven, phigiven); 

  if (fft) {
    ## to do: das liest sich alles irgendwie komisch
    maxspatialdim <- 3
   
    if (Z$spatialdim > maxspatialdim)
      stop("fft does not work yet for spatial dimensions greater than ",
           maxspatialdim)
 
    emp.vario <- n.bin <- 0
    for (i in 1:sets) {
      xx <- Z$coord[[i]]$x
      if (ncol(xx)<maxspatialdim) { # not matrix(0, ...) here!
        ##                              since x is a triple
        xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=maxspatialdim-ncol(xx)))
      }
      T3 <- if (length(Z$coord[[i]]$T) == 0) 1 else Z$coord[[i]]$T[3]
      neudim <- c(xx[3, ], if (time) T3)
      
      ## last: always repetitions
      ## last but: always vdim
      ## previous ones: coordinate dimensions

     ## Print(data, xx, T3, neudim, c(neudim, vdim, length(data[[i]]) / vdim / prod(neudim)))
      
      dim(data[[i]]) <- c(neudim, vdim, length(data[[i]]) / vdim / prod(neudim))

      
    
      
      ## to achieve a reflection in x and z instead of y we transpose the
      ## array
      crossvar <- doVario(X=data[[i]], asVector=TRUE, pseudo=pseudo, time=time)
      sumvals <- crossvar[[1]]
      nbvals <- crossvar[[2]]
    
      back <- .Call("fftVario3D", as.double(xx), 
                    as.double(sumvals), as.double(nbvals), 
                    as.double(bin), as.integer(n.bins), 
                    as.integer(T3), 
                    as.integer(stepT), as.integer(nstepT),       
                    as.double(phi), 
                    as.double(theta), 
                    as.integer(repetitions[i]),
                    as.integer(vdim),
                    totalbinsOhnevdim,
                    as.logical(pseudo), 
                    PACKAGE="RandomFields")       
   
    ## the results are now reformatted into arrays
    ## the angles are given in clear text
 
      emp.vario <- emp.vario + back[, 1]
      n.bin <- n.bin + back[, 2]
    }
    emp.vario <- emp.vario / n.bin ## might cause 0/0, but OK
    n.bin <- as.integer(round(n.bin))   
  } else {
    
    ## #####################################################################
    ##
    ## MARTINS CODE WENN FFT == FALSE
    ##
    ## #####################################################################

    if (vdim > 1) stop("multivariat only progrmmed for fft up to now")
    
    if (basic) {
      n.bin <- emp.vario.sd <- emp.vario <- 0

      
      for (i in 1:sets) {

 #       Print(i, Z$coord[[i]]$x,  Z$coord[[i]]$l, data[[i]], repetitions[i],
  #            grid[i])
        
        back <- .C("empiricalvariogram", 
                   as.double(Z$coord[[i]]$x), ## Z definition
                   as.integer(spatialdim),
                   as.integer(Z$coord[[i]]$l), 
                   as.double(data[[i]]),
                   as.integer(repetitions[i]), as.integer(grid[i]), 
                   as.double(bin), as.integer(n.bins), as.integer(0), 
                   emp.vario = double(totalbins),
                   emp.vario.sd=double(totalbins), 
                   n.bin= integer(totalbins), PACKAGE="RandomFields")        
        n.bin <- n.bin + back$n.bin
        dummy <- back$emp.vario.sd
        dummy[back$n.bin == 0] <- 0
        emp.vario.sd <- emp.vario.sd + dummy^2 * back$n.bin
 #       print(dummy)
        dummy <- back$emp.vario      
        dummy[(is.na(dummy) & (centers==0)) | back$n.bin == 0] <- 0
        emp.vario <- emp.vario + dummy * back$n.bin
  #      print(dummy)
     }
      emp.vario <- emp.vario / n.bin
      emp.vario.sd <- sqrt(emp.vario.sd / n.bin)      
      rm("back")
    } else { ## anisotropic space-time
      ## always transform to full 3 dimensional space-time coordinates
      ## with all angles given. Otherwise there would be too many special
      ## cases to treat in the c program. However, there is some lost
      ## of speed in the calculations...

      for (i in 1:sets) {        
        coord <-  Z$coord[[i]]
        xx <- coord$x
        stopifnot(is.matrix(xx))
        if (ncol(xx)<3)  # not matrix(0, ...) here! since x could be a triple
          xx <- cbind(xx, matrix(1, nrow=nrow(xx), ncol=3-ncol(xx)))

        ## x fuer grid und nicht-grid: spalte x, y, bzw z
        n.bin <- emp.vario.sd <- emp.vario <- 0

        back <-
          .C("empvarioXT", ## to do : write as .Call
             as.double(xx),
             as.double(if (length(coord$T)>0) coord$T else rep(1,3)),
             as.integer(Z$coord[[i]]$l), 
             as.double(data[[i]]),
             as.integer(repetitions[i]),
             as.integer(grid[i]), 
             as.double(bin), as.integer(n.bins), 
             as.double(c(phi[1], phi[2])), 
             as.double(c(theta[1], theta[2])), 
             as.integer(c(stepT, nstepT)), 
             ## input : deltaT[1] max abstand, deltaT[2] : echter gitterabstand,
             ##   c   : delta[1] : index gitterabstand, deltaT[2] : # of bins -1
             ##                   (zero is the additional distance)
             emp.vario = double(totalbins),
             emp.vario.sd = double(totalbins),  
             n.bin = integer(totalbins),
             PACKAGE="RandomFields")
        n.bin <- n.bin + back$n.bin
        emp.vario.sd <- emp.vario.sd + back$emp.vario.sd
        emp.vario <- emp.vario + back$emp.vario
        rm("back")
      }
      
      if (!time && vdim == 1) {
        ## vario is symmetric in phi;
        ## so the number of phi's can be halfened in this case
        dim(emp.vario) <- dims
        dim(n.bin) <- dims
        dim(emp.vario.sd) <- dims
         
        if (dims[2] > 1) {
          dims[2] <- as.integer(dims[2] / 2)
          half <- 1 : dims[2]
          n.bin <- n.bin[, half,,,,, drop=FALSE] +n.bin[, -half,,,,, drop=FALSE]
          emp.vario <- emp.vario[, half, , , , , drop=FALSE] +
            emp.vario[, -half, , , , , drop=FALSE]
          emp.vario.sd <- emp.vario.sd[, half, , , , , drop=FALSE] +
            emp.vario.sd[, -half, , , , , drop=FALSE]
          phibins <- phibins[half]
        }
      }

      emp.vario <- emp.vario / n.bin ## might cause 0/0, but OK
    
      idx <- n.bin > 1 & emp.vario != 0    
      evsd <- emp.vario.sd[idx] / (n.bin[idx] - 1) -
        n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2
      if (any(evsd < -1e-14)) {
        Print(idx, n.bin[idx] - 1, emp.vario.sd[idx], #
              emp.vario.sd[idx] / (n.bin[idx] - 1), #
              emp.vario.sd[idx] / (n.bin[idx] - 1) -
              n.bin[idx] / (n.bin[idx] -1) * emp.vario[idx]^2,
              emp.vario)
        warning(paste(evsd))
      }
      evsd[evsd < 0] <- 0
      emp.vario.sd[idx] <- sqrt(evsd)   
      emp.vario.sd[!idx] <- NaN
    }

    ## ################################################################
    ##
    ## END OF MARPINS CODE WENN FFT == FALSE
    ##
    ## ################################################################


  } # !fft
      
  
  dim(emp.vario) <- dims
  dim(n.bin) <- dims
  if (!is.null(emp.vario.sd)) dim(emp.vario.sd) <- dims
 
  name <- list()
  namedim <- names(dims)
  for (i in 1:length(dims)) {
    name[[i]] <-
      if (namedim[i] %in% c("vdim1", "vdim2")) {
          if (length(Z$varnames) == 0) NULL
          else rep(Z$varnames, length.out=dims[i])
      } else if (namedim[i] != "bins") paste(namedim[i], 1:dims[i], sep="")  
  }
  dimnames(emp.vario) <- name
#  } else names(emp.vario) <- Z$varnames[1]

  if (RFopt$general$spConform)
    l <- new("RFempVariog",
            centers=centers,
            emp.vario=emp.vario,
            var=variance,
            sd= emp.vario.sd,
            n.bin=n.bin,
            phi.centers=phibins,
            theta.centers=thetabins,
            T=Tbins,
            vdim = vdim,
            coordunits = Z$coordunits,
            varunits = varunits,
            call=call)
  else {
    l <- list(centers=centers,
              emp.vario=emp.vario,
              var=variance,
              sd= emp.vario.sd,
              n.bin=n.bin,
              phi.centers=phibins,
              theta.centers=thetabins,
              T=Tbins,
              vdim = vdim,
              coordunits =  Z$coordunits,
              varunits = varunits
           )
    class(l) <- "RF_empVariog"
  }

#  Print(l)
 # print(emp.vario)
  return(l)
   
} # function RFempiricalvariogram


## ############################################
## END OF MAIN FUNCTION 
## ############################################



doVario <- function(X, asVector=FALSE, pseudo=FALSE, time=FALSE) {
  dimX <- dim(X)
  idx.repet <- length(dimX) 
  idx.vdim <- length(dimX) - 1
  
  d <- length(dimX) - 2## last two dimensions are repet & vdim
  twoD <- dimX[3] == 1
  n <- d + pseudo 
  len<- 2^(n-1)

  
  numbers <- cubes <- array(dim=c(dimX[1:d], len, dimX[idx.repet],
                              rep(dimX[idx.vdim], 2)))
  X_list <- as.list(rep(NA, len))
  X_list[[1]] <- X
  
  ##reflect the data, carefully with time reflection
  refl.order <- if(time && !pseudo) c(1,3,4) else c(1,3,2)

  j <- 2
  for (i in 1:(n-1)) {
    for (k in 1:(2^(i-1))) {
      X_list[[j]] <- reflection(X_list[[k]], refl.order[i])
      j <- j + 1
    }      
  }

  ## to do the crossvariogram
  
  ## decide which blocks are needed
  blockidx <- rep(FALSE, 8)
  if(!time && !pseudo){
    if(twoD)  ## 2-dim case
      blockidx[1:2] <- TRUE
    else                ## 3-dim case
      blockidx[1:4] <- TRUE
  } else if(time && pseudo) {
    stop("Time component is not compatible with Pseudo variogram")
  } else { # ((time && !pseudo) || (!time && pseudo))
    if(twoD)  ## 2-dim case
      blockidx[c(1:2, 5:6)] <- TRUE
    else                ## 3-dim case
      blockidx[1:8] <- TRUE
  }

 
  
  for (i in c(1:len)){
    crossvar <- crossvario(X_list[[i]], pseudo=pseudo, dummy=!blockidx[i])
    if (time) {
      cubes[,,,,i ,,,] <- crossvar[[1]]
      numbers[,,,,i ,,,] <- crossvar[[2]]
    } else {
      cubes[,,,i ,,,] <- crossvar[[1]]
      numbers[,,,i ,,,] <- crossvar[[2]]
    }
  }

  if(asVector) return(list(as.vector(cubes), as.vector(numbers)))
 
  ##revert the reflection ## currently not used as asVector
  cubes <- crossvar[[1]]
  numbers <- crossvar[[2]]
  i<- n - 1
  for (i in (n-1):1) {
    parts<- len / (2^i)      
    positions <- 2^(i - 1)       
    for (j in 1:parts) {
      for (k in 1:positions) {
        idx <- 2* positions * j- positions + k
        if (time) {
          cubes[,,,,idx ,,,] <- reflection(cubes[,,,,idx ,,,], i)
          numbers[,,,,idx ,,,] <- reflection(numbers[,,,,idx ,,,], i)
        } else {
          cubes[,,,idx ,,,] <- reflection(cubes[,,,idx ,,,], i)
          numbers[,,,idx ,,,] <- reflection(numbers[,,,idx ,,,], i)
        }
      }
    }
  }
  return(list(cubes, numbers))
} 

crossvario<-function(f, pseudo = FALSE, dummy = FALSE) {
  d <- dim(f)
  idx.repet <- length(d) 
  idx.vdim <- length(d) - 1
  repetvdim <- c(idx.vdim, idx.repet)
  vdim <- d[idx.vdim]
  repet <- d[idx.repet]
  CVd <- c(d[-repetvdim], repet, vdim, vdim)
  if(dummy) return(list(array(1, dim=CVd), array(1, dim=CVd)))

  
  idx <- rep(TRUE, length(d) - 2)
  idx.data <- paste("[", paste(1, ":", d, collapse=", "), "]")
  idx.vario <- paste("[", paste(rep(",", length(d)-2), collapse=""), "r, i, j]")
  idx.w <- paste("[", paste(1, ":", d[-repetvdim], collapse=", "), "]")
 
  dim.coord <- 2 * d[-repetvdim]-1
  F <- If <- array(0, dim=c(dim.coord, d[repetvdim]))
  eval(parse(text=paste("If", idx.data, "<- !is.na(f)")))
  f[is.na(f)] <- 0
  eval(parse(text=paste("F", idx.data,  "<- f")))
  LIf <- list(If)
  LF <- list(F)

  nbvals <- Crossvario <- array(0, CVd)
   
  for (i in 1:vdim) {
    for (j in 1:vdim) {
      for (r in 1:repet) {
        #
        
        
        If <- do.call("[", c(LIf, idx, i, r))
        dim(If) <- dim.coord
        Ig <- do.call("[", c(LIf, idx, j, r))        
        dim(Ig) <- dim.coord
        F <- do.call("[", c(LF, idx, i, r))
        dim(F) <- dim.coord
        G <- do.call("[", c(LF, idx, j, r))
        dim(G) <- dim.coord
        if (!pseudo) {    
          fftIfIg <- fft(If * Ig)
          fftFG <- fft(F * G)
          fftIfG <- fft(G * If)
          fftIgF <- fft(F * Ig)   
          z <- fft(Conj(fftFG) * fftIfIg
                   + Conj(fftIfIg) * fftFG
                   - Conj(fftIgF) * fftIfG
                   - Conj(fftIfG) * fftIgF, inverse=TRUE)
          N <- fft( Conj(fftIfIg) * fftIfIg, inverse=TRUE )
        } else {
          F2 <- F^2
          G2 <- G^2
          fftIf <- fft(If)
          fftIg <- fft(Ig)
          z <- fft( Conj(fft(F2))* fftIg
                   + Conj(fftIf) * fft(G2)
                   - 2* Conj(fft(F)) * fft(G), inverse=TRUE)
          ## N <- 2* fft(Conj(fftIf)*fftIg, inverse=TRUE)
          N <- fft(Conj(fftIf)*fftIg, inverse=TRUE)
        }

        w <- Re(z) / (2 * prod(dim(N))) # sumvals
  
        
        eval(parse(text=paste("Crossvario", idx.vario, "<- w", idx.w)))
        eval(parse(text=paste("nbvals", idx.vario,
                     "<- Re(N", idx.w, ") / prod(dim(N))")))
      }
    }
  }  
  return(list(Crossvario, as.array(round(nbvals))))
}


prepareBin <- function(bin)
{
  if(missing(bin)) return(NULL)
  if (bin[1] > 0) {
      if (RFoptions()$general$printlevel>1)
	message("empirical variogram: left bin border 0 added\n")
      bin <- c(0, bin)
    }
  if (bin[1]==0) bin <- c(-1, bin)
  if (bin[1] < 0) bin <- c(bin[1], bin[bin>=0])
  
  bin
}

