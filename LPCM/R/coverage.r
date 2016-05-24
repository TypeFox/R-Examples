
coverage.raw <-function(X, vec, tau, weights=1, plot.type="p", print=FALSE, label=NULL,...){
     X<- as.matrix(X)
     p <- dim(vec)[1]
     n <- dim(X)[1]
     d <- dim(X)[2]

     if (n %% length(weights) !=0){
         stop("The length of the vector of weights is not a multiple of the sample size.")
     } else {
      weights <-rep(weights, n %/%   length(weights))
     }   

     min.distance  <- rep(0,n)
     ins.distance <- rep(1,n)

     for (i in 1: n){
        #if (i %%10==0){ print(i)}
        if (is.null(label)){
             min.distance[i] <- mindist(vec, X[i,])$mindist
        } else {
            if (!is.matrix(vec)){vec<-matrix(vec, nrow=1)}
            if (as.character(label[i]) %in%  dimnames(vec)[[1]]){
                  min.distance[i] <- mindist(matrix(vec[as.character(label[i]),],nrow=1), X[i,])$mindist   # 11/10/10 experimental, for clustering
            } else { 
                min.distance[i]<- tau+1   #01/11/10 if data point not allocatable to any centre.
            }
        }    
        ins.distance[i] <- (min.distance[i] <= tau) #indicator for being inside/outside the tube
           }
     ci<- weighted.mean(min.distance <= tau, w=weights) 
     if (plot.type %in% c("p","l")) {
        plot(X, col=ins.distance+1,...)
        if (plot.type=="p"){points(vec, col=3,lwd=2 )} else if (plot.type=="l"){lines(vec, col=3,lwd=2 )}
        }
     if (print){print(c(tau,ci))}
     return(list(tau=tau, coverage= ci, min=min.distance, inside= ins.distance))
 }



coverage<-function(X, vec,  taumin=0.02, taumax,  gridsize=25, weights=1, plot.type="o", print=FALSE,...){
  if (missing(taumax)){
  m <-colMeans(X)
  Xm <- sweep(X, 1, m, "-")
  Xs<- apply(X,1, sd)
  taumax<-max(Xs)
 }
 all.taus      <- seq(taumin, taumax, length=gridsize)
 all.coverages <- rep(0, gridsize)

 for (j in 1:gridsize){
    all.coverages[j]<- coverage.raw(X, vec, all.taus[j],  weights, plot.type=0, print)$coverage
 }

 if (plot.type!=0){
   plot(all.taus, all.coverages, type=plot.type, xlab=expression(tau), ylab=expression(C(tau)),...)
 }
 return(list(tau=all.taus, coverage=all.coverages))
}


lpc.coverage<-function(object, taumin=0.02, taumax, gridsize=25,  quick=TRUE, plot.type="o", print=FALSE, ...){
 
 if (class(object)=="lpc"){
   X <-object$data
   scaled<-object$scaled
   weights <- object$weights
   if (quick){
       lpc.vec<- object$LPC
   } else {
       lpc.vec<-lpc.spline(object, project=TRUE)$closest.coords
   }
 } else if (class(object)=="lpc.spline"){
    X <-object$lpcobject$data
    scaled <- object$lpcobject$scaled
    weights <- object$lpcobject$weights
    if (object$closest.coords[1] != "none") {
       lpc.vec<- object$closest.coords
     } else if (quick){
       lpc.vec<- object$lpcobject$LPC
     } else {
       lpc.vec<- lpc.spline(object$lpcobject,project=TRUE)$closest.coords
    }
  } else {
    stop("Invalid lpc object.")
  }
    
   if (missing(taumax)){   
      if (scaled){taumax<-0.5} else {
         m <-colMeans(X)
         Xm <- sweep(X, 2, m, "-")
         Xs<- apply(X,2, sd)
        taumax<-max(Xs)
      } 
    }
 
 
  result <- coverage(X, lpc.vec, taumin, taumax, gridsize,  weights, plot.type=0, print)
  if (plot.type!=0){
   plot(result$tau,result$coverage, type=plot.type, xlab=expression(tau), ylab=expression(C(tau)),...)
 }
  return(result)
}  
  
  

lpc.self.coverage <-
function (X, taumin = 0.02, taumax = 0.5, gridsize = 25, x0=1, 
    way = "two", scaled = TRUE, weights = 1, pen = 2, 
    depth = 1, control = lpc.control(boundary = 0, cross = FALSE), 
    quick = TRUE, plot.type = "o", print = FALSE, ...) 
    {
    if (class(X) %in% c("lpc", "lpc.spline")) {
        stop("Invalid data matrix.")
    }
   
    
    Xi <- as.matrix(X)
    N <- dim(Xi)[1]
    d <- dim(Xi)[2]
    
    s1       <- apply(Xi, 2, function(dat){ diff(range(dat))})  # range
    
    mult <- control$mult
   
    if (length(x0)==1){
       ms.sub<-control$ms.sub
       if (!is.null(control$ms.h)){
          ms.h<-control$ms.h
       } else {
          if (!scaled){ms.h   <- s1/10} else {ms.h<- 0.1}
       }   
    }    

   
    if (length(x0)==1 && x0==1){   
      n  <- sample(N,1);
      X  <-  if (scaled){ sweep(Xi, 2, s1, "/")} else {Xi}     
      x0 <- matrix(ms.rep(X, X[n,],ms.h, plotms=0)$final, nrow=1)
      x0 <- if (scaled) x0*s1 else x0   # unscales again
      rm(X)
     }  else if (length(x0)==1 && x0==0 ){
          if (N<= ms.sub) {
              sub<- 1:ms.sub
          }  else {    
              Nsub <- min(max(ms.sub, floor(ms.sub*N/100)), 10*ms.sub)
              #print(Nsub)
              sub <- sample(1:N, Nsub)
          }
          x0<- suppressWarnings(unscale(ms(Xi, ms.h, subset=sub, plotms=0, scaled=scaled)))$cluster.center
     }  else {
        if (is.null(x0)){
             if (is.null(mult)){stop("One needs to allow for at least one starting point.")}
             x0<- matrix(0, nrow=0, ncol=d)
        } else {    
             x0 <- matrix(x0, ncol=d, byrow=TRUE)
        }   
    } 
    if(!is.null(control$mult)){
     #  stop("It is not permitted to modify mult for this operation.")
      if (dim(x0)[1] < mult) {n <- runif(mult-dim(x0)[1],1,N+1)%/%1; x0 <- rbind(x0,Xi[n,])} #
      if (dim(x0)[1] > mult) {x0<-x0[1:mult]} 
    }
    x0       <- as.matrix(x0)    # putting in matrix format; just in case....

    #print(x0)
    if ((!scaled) && taumax < 1) {
        warning("Please adjust the range (taumin, taumax) of tube widths by hand, as the data are not scaled.")
    }
    Pm <- NULL
    h0 <- taumin
    h1 <- taumax
    h <- seq(h0, h1, length = gridsize)
    #n <- gridsize
    cover <- matrix(0, gridsize, 2)
    for (i in 1:gridsize) {
        new.h0 <- h[i]
        fit <- lpc(Xi, h = new.h0, t0 = new.h0, x0 = x0,
            way = way, scaled = scaled, weights = weights, pen = pen, 
            depth = depth, control)
        if (!quick) {
            fit.spline <- lpc.spline(fit, project = TRUE)
        }
        if (quick) {
            Pm[[i]] <- fit$LPC
        }
        else {
            Pm[[i]] <- fit.spline$closest.coords
        }
        cover[i, ] <- as.numeric(coverage.raw(fit$data, Pm[[i]], 
            new.h0, weights, plot.type = 0, print = print)[1:2])
    }

   
    select <- select.self.coverage(self = cover,  
        smin = 2/3, plot.type = plot.type)
    result <- list(self.coverage.curve = cover, select = select, x0.unscaled=suppressWarnings(unscale(fit))$start,
        type = "lpc")
    class(result) <- "self"
    result
}


select.self.coverage <-
function (self,  smin, plot.type = "o", plot.segments=NULL) 
{
    if (class(self) == "self") {
        cover <- self$self.coverage.curve
    }
    else {
        cover <- self
    }
    if (missing(smin)) {
        if (class(self) == "self") {
            smin <- switch(self$type, lpc = 2/3, ms = 1/3)
        }
        else stop("Please specify `smin' argument.")
    }
    n <- dim(cover)[1]
    diff1 <- diff2 <- rep(0, n)
    diff1[2:n] <- cover[2:n, 2] - cover[1:(n - 1), 2]
    diff2[2:(n - 1)] <- diff1[3:n] - diff1[2:(n - 1)]
    select <- select.coverage <- select.2diff <- NULL
       
    if (plot.type != 0) {
        plot(cover, type = plot.type, xlab = "h", ylab = "S(h)", 
            ylim = c(0, 1))
    }
    for (i in (3:(n - 1))) {
        if (diff2[i] < 0 && cover[i, 2] > max(smin, cover[1:(i - 
            1), 2])) {
            select <- c(select, cover[i, 1])
            select.coverage <- c(select.coverage, cover[i,2])
            select.2diff <- c(select.2diff, diff2[i])
            #if (plot.type != 0) {
            #    segments(cover[i, 1], 0, cover[i, 1], cover[i, 
            #      2], col = scol[i], lty = slty[i], lwd=slwd[i])
            #}
        }
    }
             
    selected <-select[order(select.2diff)]
    covered<- select.coverage[order(select.2diff)]          
    
    if (plot.type != 0) {
      d<-length(selected)
      slty <- slwd <-scol<-rep(0,d)
      if (is.null(plot.segments)){
         slty[1:3]<- c(1,2,3)
         slwd[1:3] <-c(2,1,1)
         scol[1:3] <- c(3,3,3)
      } else {
          r<-max(length(plot.segments$lty), length(plot.segments$lwd), length(plot.segments$col))
            slty[1:r] <- slwd[1:r] <-scol[1:r]<-1
            if (length(plot.segments$lty)>0){ slty[1:length(plot.segments$lty)]<- plot.segments$lty}
             if (length(plot.segments$lwd)>0){slwd[1:length(plot.segments$lwd)]<- plot.segments$lwd}
            if (length(plot.segments$col)>0){ scol[1:length(plot.segments$col)]<- plot.segments$col}
       }          
      for (j in 1:d){
           segments(selected[j], 0, selected[j], covered[j], col = scol[j], lty = slty[j], lwd=slwd[j])
         }
      }         
             
    return(list("select"=selected, "select.coverage"=covered, "select.2diff"=select.2diff[order(select.2diff)]))
}

