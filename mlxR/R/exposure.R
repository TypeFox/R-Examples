#' Computation of AUC, Cmax and Cmin
#'
#' Compute the area under the curve, the maximum and minimum values of a function of time over a given interval or at steady state
#' 
#' Input arguments are the input arguments of Simulx (http://simulx.webpopix.org)
#'  
#' Specific input arguments can be also used for computing the exposure at steady state,
#' i.e. after the administration of an "infinite" number of doses.              
#' See http://simulx.webpopix.org/exposure/ for more details. 
#' 
#' @param output a list  with fields: 
#' \itemize{
#'   \item \code{name}: a vector of output names
#'   \item \code{time}: = 'steady.state' 
#'   \item \code{ntp}: number of time points used for computing the exposure (default=100) 
#'   \item \code{tol}: tolerance number, between 0 and 1, for approximating steaty-state  (default=0.01) 
#'   \item \code{ngc}: number of doses used for estimating the convergence rate to steaty-state  (default=5) 
#' }
#' @param treatment a list with fields
#' \itemize{
#'   \item \code{tfd} : time of first dose (default=0),
#'   \item \code{ii} : inter dose interval (mandatory),
#'   \item \code{amount} : the amount for each dose,
#'   \item \code{rate} : the infusion rate  (default=\code{Inf}),
#'   \item \code{tinf} : the time of infusion  (default=0),
#'   \item \code{type} : the type of input (default=1),
#'   \item \code{target} : the target compartment (default=NULL). 
#' }
#' 
#' @param model a \code{Mlxtran} or \code{PharmML} model used for the simulation
#' @param group a list, or a list of lists, with fields: 
#' \itemize{
#'   \item \code{size} : size of the group (default=1),
#'   \item \code{level} : level(s) of randomization,
#'   \item \code{parameter} : if different parameters per group are defined,
#'   \item \code{output} : if different outputs per group are defined,
#'   \item \code{treatment} : if different treatements per group are defined,
#'   \item \code{regressor} : if different regression variables per group are defined.
#' }
#' @param parameter a vector of parameters with their names and values
#' @param data a list
#' @param project the name of a Monolix project
#' @param settings a list of optional settings
#' \itemize{
#'   \item \code{record.file} : name of the datafile where the simulated data is written (string),
#'   \item \code{seed} : initialization of the random number generator (integer),
#'   \item \code{load.design} : TRUE/FALSE (if load.design is not defined, a test is automatically performed to check if a new design has been defined),
#'   \item \code{data.in} : TRUE/FALSE (default=FALSE)
#'   \item \code{id.out}  : add columns id (when N=1) and group (when #group=1), TRUE/FALSE (default=FALSE)
#'   \item \code{Nmax} : maximum group size used in a single call of mlxCompute (default=100)
#' }   
#' @param regressor a list, or a list of lists, with fields
#' \itemize{
#'   \item \code{name} : a vector of regressor names,
#'   \item \code{time} : a vector of times,
#'   \item \code{value} : a vector of values.
#' }
#' @param varlevel a list, or a list of lists, with fields
#' \itemize{
#'   \item \code{name} : a vector of names of variability levels,
#'   \item \code{time} : a vector of times that define the occasions.
#' }
#' @return A list of data frames. One data frame per output is created with columns \code{id} (if number of subject >1),
#' \code{group} (if number of groups >1), \code{t1} (beginning of time interval), \code{t2} (end of time interval),
#' \code{step} (time step), \code{auc} (area under the curve), \code{tmax} (time of maximum value), \code{cmax} (maximum value),
#' \code{tmin} (time of minimum value), \code{cmin} (minimum value).
#' 
#' @importFrom utils tail
#' @export
exposure <- function(model,output, group=NULL,treatment=NULL,parameter=NULL,
                     data=NULL, project=NULL, settings=NULL, regressor=NULL, varlevel=NULL)
{  
  if (!is.null(group)){
    if (!is.null(group$output))
      stop("\n'output' cannot be defined in 'group' with exposure\n")
  }
  if (identical(output$time,"steady.state")){
    
    if (is.null(output$ntp)) {ntp<-100}  else {ntp<-output$ntp}
    if (is.null(output$ngc)) {ngc<-5}    else {ngc<-output$ngc}
    if (is.null(output$tol)) {tol<-0.01} else {tol<-output$tol}
        
    # groups and treatment
    if (is.null(group))
      group <- list(NULL)
    
    if (!is.null(names(group)))  
      group <- list(group) 
    G <- length(group)
    
    if (!is.null(treatment)){
      for (g in (1:G))
        group[[g]]$treatment <- treatment
      treatment <- NULL
    }
    for (g in (1:G)){
      if (!is.null(names(group[[g]]$treatment)))  
        group[[g]]$treatment <- list(group[[g]]$treatment) 
    } 
    
    # ppcm
    pmtau=1
    tfd <- Inf
    
    for (g in (1:G)){
      for (k in seq(1,length(group[[g]]$treatment))){
        trtk <- group[[g]]$treatment[[k]]
        if (is.null(trtk$ii))
          stop("inter dose interval (ii) should be an element of treatment when exposure at steady state is computed")
        pmtau <- ppcm(pmtau,trtk$ii) 
        if(!any("tfd" %in% names(trtk)))
          group[[g]]$treatment[[k]]$tfd <- 0
        tfd <- min(tfd, group[[g]]$treatment[[k]]$tfd)
      }
    }
    #--------------------------------
    mgc <- 2
    T <- ngc*pmtau
    gr1 <- group
    for (g in (1:G)){
      trtg <- group[[g]]$treatment
      trt1 <- trtg
      for (k in seq(1,length(trtg))){
        trtk <- trtg[[k]]
        trtk$time <- seq(trtk$tfd,T+trtk$tfd,by=trtk$ii)   
        trtk$ii <- NULL
        trtk$tfd <- NULL
        trt1[[k]] <- trtk
      }
      gr1[[g]]$treatment <- trt1
    }
    out1 <- list(name=output$name, time=seq(tfd,T+tfd-pmtau/mgc,by=pmtau/mgc))
    if (length(gr1)==1){
      trt1 <- gr1[[1]]$treatment
      gr1[[1]]$treatment <- NULL
      if (length(gr1[[1]])==0)
        gr1 <- NULL
    }else{
      trt1 <- NULL
    }   
    r1 <- simulx(model=model,group=gr1,treatment=trt1, parameter=parameter,
                 output=out1,settings=settings,varlevel=varlevel)
    r1$treatment <- NULL
    r.names <- names(r1)
    alpha <- Inf
    for (k in seq(1:length(r1))){
      namek <- r.names[k]
      rk <- r1[[k]]
      if (any("time" %in% names(rk))){
        i.id <- any("id" %in% names(rk))
        if (!i.id) {rk$id <- factor(1)}
        N <- nlevels(rk$id)
        i2 <- 0
        for (i in (1:N)){
          i1 <- i2 + 1 
          i2 <- i2 + ngc*mgc
          rki <- rk[i1:i2,]
          fki <- matrix(rki[namek][[1]],ncol=ngc)
          dki <- diff(apply(fki,2,mean))
          alphai <- diff(log(abs(dki)))
          if (max(alphai) <0){
            alpha <- min(alpha,-tail(alphai,n=1))
          }else if (min(alphai)>0){
            stop("\n\nSorry... exposure was not able to estimate the rate of convergence to steady state...
Try increasing ngc, or fix the number of doses")
          }
        }
      }
    }
    M <- ceiling(-log(tol)/alpha)+1
    #----------------------------------
    
    T <- M*pmtau + tfd
    for (g in (1:G)){
      trtg <- group[[g]]$treatment
      for (k in seq(1,length(trtg))){
        trtk <- trtg[[k]]
        trtk$time <- seq(trtk$tfd,T+trtk$tfd,by=trtk$ii)   
        trtk$ii <- NULL
        trtk$tfd <- NULL
        trtg[[k]] <- trtk
      }
      group[[g]]$treatment <- trtg
    }
    output <- list(name=output$name, time=seq(T-pmtau,T,length=ntp))
    if (length(group)==1){
      treatment <- group[[1]]$treatment
      group[[1]]$treatment <- NULL
      if (length(group[[1]])==0)
        group <- NULL
    }
  }
  t <- output$time
  t.min <- min(t)
  t.max <- max(t)
  t.n   <- length(t)
  t.step <- (t.max - t.min)/(t.n-1)
  t.time <- seq(t.min,t.max, by=t.step)
  output$time <- t.time
  
  r.simul <- simulx(model=model,group=group,treatment=treatment,parameter=parameter,
                    output=output,data=data,project=project,settings=settings,
                    regressor=regressor,varlevel=varlevel)
  r.simul$treatment <- NULL
  r.names <- names(r.simul)
  res <- list()
  kk <- 0
  for (k in seq(1:length(r.simul))){
    namek <- r.names[k]
    rk <- r.simul[[k]]
    if (any("time" %in% names(rk))){
      kk <- kk+1
      i.id <- any("id" %in% names(rk))
      i.group <- any("group" %in% names(rk))
      if (!i.id){
        rk$id <- 1
      }
      N <- nlevels(rk$id)
      cmin <- vector(length=N)
      tmin <- vector(length=N)
      cmax <- vector(length=N)
      tmax <- vector(length=N)
      auc  <- vector(length=N)
      gr   <- vector(length=N)
      i2 <- 0
      for (i in (1:N)){
        i1 <- i2 + 1 
        i2 <- i2 + t.n
        rki <- rk[i1:i2,]
        tki <- rki["time"][[1]]
        fki <- rki[namek][[1]]
        imin <- which.min(fki)
        tmin[i] <- tki[imin]
        cmin[i] <- fki[imin]
        imax <- which.max(fki)
        tmax[i] <- tki[imax]
        cmax[i] <- fki[imax]
        auc[i] <- ((fki[1]+fki[t.n])/2+sum(fki[2:(t.n-1)]))*t.step
        if (i.group){ gr[i] <- rki$group[1]  }
      }
      if (!i.id){
        res[[kk]] <- data.frame(t1=t.min, t2=t.max, step=t.step, auc, 
                                tmax, cmax, tmin, cmin)
      }else if (!i.group){
        res[[kk]] <- data.frame(id=factor(1:N), t1=t.min, t2=t.max, step=t.step, auc, 
                                tmax, cmax, tmin, cmin)
      }else{
        res[[kk]] <- data.frame(id=factor(1:N), group=as.factor(gr), t1=t.min, t2=t.max, step=t.step, auc, 
                                tmax, cmax, tmin, cmin)        
      }
      names(res)[kk] <- namek
    }
  }
  res$output <- r.simul
  return(res)
}


pgcd <- function(a,b){ 
  if (abs(b-0)<1e-7) { 
    pgs <- a 
  } 
  else { 
    r <- a%%b 
    pgs <- pgcd(b,r) 
  } 
  return(pgs) 
} 

ppcm<-function(a,b) a*b/(pgcd(a,b))

