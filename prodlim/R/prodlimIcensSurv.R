prodlimIcensSurv <- function(response,
                             grid,
                             tol=7,
                             maxiter,
                             ml=FALSE,
                             exact=TRUE){

  # {{{ data
  ntol <- 10^{-tol}
  L <- response[,"L"]
  N <- length(L)
  R <- response[,"R"]
  status <- response[,"status"]
  # }}}
  # {{{ one-step idea

  if (ml==FALSE) {
      # right censored observations
      # are defined by status
      R[status==0] <- L[status==0]
      if (missing(grid))
          grid <- sort(unique(c(L,R)))
      else
          if (exact)
              grid <- sort(unique(c(min(L,R),grid)))
          else
              grid <- sort(unique(grid))
      ## need an extra grid point before the smallest
      ## `L' to catch right censored and exact
      ## event times that match this smallest `L'
    
    stopifnot(all(grid >= 0))
    if (grid[1]==0)
      grid <- c(-1,grid)
    else
      grid <- c(0,grid)
    
    indexR <- sindex(jump.times=grid,eval.times=R)
    indexL <- sindex(jump.times=grid,eval.times=L)

      ## indexR <- match(R,grid)
      ## indexL <- match(L,grid)
      NS <- length(grid)
      Ind <- iindex(L,R,grid)
       ## fit <- list("icens_prodlim",
       ## as.double(L),
       ## as.double(R),
       ## as.double(grid),
       ## as.integer(indexL),
       ## as.integer(indexR),
       ## as.integer(Ind$iindex),
       ## as.integer(c(Ind$imax,0)),
       ## as.integer(status),
       ## as.double(N),
       ## as.double(NS),
       ## nrisk=double(NS),
       ## nevent=double(NS),
       ## ncens=double(NS),
       ## hazard=double(NS),
       ## varhazard=double(NS),
       ## surv=double(NS),
       ## oldsurv=double(NS),
       ## as.double(ntol),
       ## as.integer(maxiter),
       ## n.iter=integer(1),
       ## package="prodim")
      fit <- .C("icens_prodlim",
                as.double(L),
                as.double(R),
                as.double(grid),
                as.integer(indexL),
                as.integer(indexR),
                as.integer(Ind$iindex),
                as.integer(c(Ind$imax,0)),
                as.integer(status),
                as.double(N),
                as.double(NS),
                nrisk=double(NS),
                nevent=double(NS),
                ncens=double(NS),
                hazard=double(NS),
                varhazard=double(NS),
                surv=double(NS),
                oldsurv=double(NS),
                as.double(ntol),
                as.integer(maxiter),
                n.iter=integer(1),
                package="prodim")
      ## rename the extra grid point before the smallest `L'
      ## if it is negative
      if (grid[1]<0) grid[1] <- 0
      res <- list("time"=rbind(c(0,grid[-length(grid)]),c(grid)),
                  "n.risk"=round(pmax(0,fit$nrisk),tol),
                  "n.event"=round(pmax(0,fit$nevent),tol),
                  "n.lost"=round(fit$ncens,tol),
                  "hazard"=round(fit$hazard,tol),
                  "surv"=round(pmax(0,fit$surv),tol),
                  "maxtime"=max(grid),
                  "n.iter"=fit$n.iter,
                  "tol"=ntol,
                  "model"="survival")
      #    res <- list("time"=rbind(c(0,0,grid[-length(grid)]),c(0,grid)),"n.risk"=c(N,round(pmax(0,fit$nrisk),tol)),"n.event"=c(0,round(pmax(0,fit$nevent),tol)),"n.lost"=c(0,round(fit$ncens,tol)),"hazard"=c(0,round(fit$hazard,tol)),"surv"=c(1,round(pmax(0,fit$surv),tol)),"maxtime"=max(grid),"n.iter"=fit$n.iter,"tol"=ntol,"model"="survival")
  }
  else{
      # }}}
      # {{{ npmle 

    
    ## artificial closure of right censored intervals 
    ## R[Rna] <- max(c(L,R)) + 1
    R[status==0] <- max(c(L,R[status!=0])) + 1
    ##     R[status==0] <- max(c(L,R)) + 1
    ##     print(R[status==0])
    peto.intervals  <-  PetoInt(L,R,status)
    indices <- IntIndex(x=peto.intervals,L=L,R=R)
    Mindex <- indices$Mindex
    Mstrata <- indices$Mstrata
    Iindex <- indices$Iindex
    Istrata <- indices$Istrata
    M <- length(Mstrata)
    N <- length(Istrata)
    ## Zsurv <- predictSurv(prodlimIcensSurv(response=response,grid=grid,tol=tol,maxiter=1,ml=FALSE))
    Z <- rep(1/M,M)
    fit  <- .C('GMLE',as.integer(c(0,Mstrata)),as.integer(c(0,Istrata)),as.integer(Mindex),as.integer(Iindex),as.integer(N),as.integer(M),Z=as.double(Z),double(length(Z)),as.double(ntol),as.integer(maxiter),steps=integer(1),package="prodlim")
    n.event <- c(0,fit$Z*M)
    surv <- c(1,1-cumsum(fit$Z))
    hazard <- c(0,fit$Z)/surv
    res <- list("time"=cbind(c(0,0),peto.intervals),"n.risk"=N-n.event,"n.event"=n.event,"n.lost"= c(0,rep(0,M)),"hazard"=round(hazard,tol),"surv"=round(surv,tol),"maxtime"=max(c(peto.intervals)),"n.iter"=fit$steps,"tol"=ntol,"model"="survival")
  }
  # }}}

  class(res) <- "prodlim"
  res  
}
