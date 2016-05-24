CodomMarker <-
function(y, ng=5,
                      mutype=0, sdtype="sd.const", ptype="p.free",
                      clus=TRUE, mu.start=NA, sd.start=rep(0.075,ng), sd.fixed=0.075, p=NA,
                      maxiter=500, maxn.bin=200, nbin=200,
                      plothist=TRUE, nbreaks=40,
                      maintitle=NULL, subtitle=NULL, xlabel=NULL, xaxis="s") {
tryCatch( {
  if (ng < 1) simpleError("ng<1")
  if (length(y) < 10*ng) simpleError("y should have at least 10*ng elements")
  if (length(y[is.na(y)]) > 0) simpleError("y may not contain NA")
  if (length(y[y<0 | y>1]) > 0) simpleError("all y values should be between 0 and 1 (inclusive)")
  
  if (is.na(p[1]) || sum(p)==0) p <- rep(1/ng,ng)
  else p <- p/sum(p)

  if (length(mu.start) != ng) mu.start<-NA
  if (length(sd.start) != ng) sd.start<-NA
  
  # transform data to arcsine(square root):
  yw <- asin(sqrt(y))
  if (!is.na(mu.start[1])) mu.start<-asin(sqrt(mu.start))

  minyw <- min(yw, na.rm=T)
  maxyw <- max(yw, na.rm=T)

  #determine the number of free parameters of the model:
  if (mutype %in% 0:10) mupar <- getMuModelNpar(mutype,ng)
  else  simpleError("invalid mutype")
  sigmapar = switch (sdtype,
    sd.free = ng,
    sd.const = 1,
    sd.fixed = 0,
    simpleError("invalid sdtype")
  )
  ppar = switch (ptype,
    p.free = ng-1,
    p.HW = 1,
    p.fixed = 0,
    simpleError("invalid ptype")
  )
  npar <- mupar+sigmapar+ppar

  result<-list(message="Error in CodomMarker") #default, to be replaced

  #check for ng==1 (1 peak):
  if (ng==1) {
    result <- fitOneDist(yw, sdtype, sd.fixed, npar)
    }
  else {
    if (is.na(mu.start[1])) {  # no starting values for mu were give, so calculate these first

      # initialize parameters through clustering if asked for
      if (clus) {
        # original: using hclust; very slow for large sample numbers
        #isna <- is.na(yw); yh <- yw[!isna]; n <- length(yh)
        # priorclass <- cutree(hclust(dist(yh), method="average"), ng)
        # mu.start <- aggregate(data.frame(yh,priorclass), list(priorclass), mean)[,2]
        # o <- order(mu.start); mu.start <- mu.start[o]
        # sd.start <- aggregate(data.frame(yh,priorclass), list(priorclass), sd)[,2]
        # sd.start <- sd.start[o]
      
        # alternative (3-5-2012): use kmeans clustering, much faster 
        yh <- yw[!is.na(yw)]
        # .Random.seed (affects and is affected by kmeans) is saved, set and restored in fitTetra

        clus.init <- ClusterInit(yh, ng) # returns mu's and sd on transformed scale

        mu.start <- clus.init$clus.mu
        #sd.start <- rep(clus.init$clus.sd, ng)  # start with equal standard deviations
        sd.start <- clus.init$clus.sd #no need for rep, already ng values
        sd.start[sd.start<0.01] <- 0.01  # don't allow sds smaller than 0.01
        
        if (sdtype=="sd.const") sd.start <- rep(mean(sd.start),ng)
        else if (sdtype=="sd.fixed") sd.start <- rep(sd.fixed,ng)
      }
      else { # mu.start==NA and clus==FALSE; spread evenly over range
        ylo <- asin(sqrt(0.02))
        yhi <- asin(sqrt(0.98))
        mu.start <- seq(ylo, yhi, by=(yhi-ylo)/(ng-1))
      } 
    } else {
      # mu.start given, check if sd.start is available:
      if (is.na(sd.start[1])) {
        sd.start=rep(0.075,ng) # as the default parameter value
      }  
    }

    result <- EMGaussMix(yw, ng=ng, ptype=ptype, mutype=mutype, sdtype=sdtype, sd.fixed=sd.fixed,
                         p.start=p, mu.start=mu.start, sd.start=sd.start,
                         npar=npar, maxiter=maxiter, maxn.bin=maxn.bin, nbin=nbin)
  } # if ng==1 else
  if ("message" %in% names(result)) {
    if (result$message=="") {
      # transform back to original scale
      mu.back <- sin(result$psi$mu)^2 + 0.5*(2*(cos(result$psi$mu)^2 - sin(result$psi$mu)^2))*(result$psi$sigma^2)   # Use E(f(y)) ~ f(Ey) + 0.5*f''(Ey) var(y)
      sigma.back <- 2*sin(result$psi$mu)*cos(result$psi$mu)*result$psi$sigma
      back <- list(mu.back=mu.back, sigma.back=sigma.back)
      result$back <- back

      # prepare for histogram visualization
      if (plothist==TRUE)
        PlotHistDensity(y, result, trafo="asr",
                        maintitle=maintitle, subtitle=subtitle, xlabel=xlabel,
                        nbreaks=nbreaks)

      result
    }
    else { #result$message!=""
      simpleError(paste("Error1 in CodomMarker: ",result$message,sep=""))
    }
  } #"message" in names(result)
  else simpleError("Error3 in CodomMarker")
}, error = function(ex) { simpleError(paste("Error2 in CodomMarker: ",ex,sep="")) } )
}
