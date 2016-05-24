#SMACOF for individual differences (list of dissimilarity matrices)

smacofIndDiff <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
                          constraint = c("indscal", "idioscal", "identity"),
                          weightmat = NULL, init = "torgerson", ties = "primary", 
                          verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6,
                          spline.degree = 2, spline.intKnots = 2)
  
# delta ... list of input objects: either of class dist() or a symmetric matrix
# contstraint ... either NULL, "identity", "diagonal", "idioscal"
# ties ... primary, secondary, tertiary
{
  
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  constraint <- match.arg(constraint, c("indscal", "idioscal", "identity"), several.ok = FALSE)
  
  diss <- delta
  p <- ndim
  if (constraint == "indscal") constraint <- "diagonal"
  constr <- constraint
  if (!is.list(diss)) diss <- list(diss)
  if ((is.matrix(diss[[1]])) || (is.data.frame(diss[[1]]))) diss <- lapply(diss, strucprep)
  checkdiss(diss)           ## sanity check
  
  ## --- weight matrix
  if (is.null(weightmat)) wgths <- initWeights(diss) else wgths <- weightmat
  if (!is.list(wgths)) {
    wgths <- list(wgths)
    if (length(wgths) != length(diss)) wgths <- sapply(diss, function(wwr) return(wgths))
  }
  if ((is.matrix(wgths[[1]])) || (is.data.frame(wgths[[1]]))) wgths <- lapply(wgths, strucprep)  
  
  ## --- Prepare for optimal scaling
  trans <- type
  if (trans=="ratio"){
    trans <- "none"
  } else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
  } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
  } else if(trans=="spline"){
    trans <- "mspline"
  }
  disobj <- list()
  for (i in 1:length(diss)) {
      disobj[[i]] <- transPrep(diss[[i]],trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
  }
  ## --- end optimal scaling prep
  
  
  n <- attr(diss[[1]],"Size")
  if (p > (n - 1)) stop("Maximum number of dimensions is n-1!")
  
  nn <- n*(n-1)/2
  m <- length(diss)
  itel <- 1
  
  if (is.null(attr(diss[[1]], "Labels"))) {
     for (i in 1:m) attr(diss[[i]], "Labels") <- paste(1:n)
  }
        
  dr <- list()
  wr <- list()
  vr <- list()
  dh <- list()
      
  for (j in 1:m) {                              #initialize weights, V, norm d as lists
	  wr <- appendList(wr,vmat(wgths[[j]]))              
	  vr <- appendList(vr,myGenInv(wr[[j]]))
	  dh <- appendList(dh,normDissN(diss[[j]],wgths[[j]],1))
  }
  xr <-list()                                  #configurations as list
  sold <- sf1 <- sf2 <- 0                      #stress init

  ## --- starting values
  aconf <- initConf(init, diss, n, p, inddiff = TRUE)
    
  
  bconf <- repList(diag(p),m)                  #1-matrix
  for (j in 1:m) {                             
    xr[[j]] <- aconf%*%bconf[[j]]  #same starting values for all ways
    dr[[j]] <- dist(xr[[j]])                          #configuration distances
    sf1 <- sf1 + sum(wgths[[j]]*dr[[j]]*dh[[j]])
    sf2 <- sf2 + sum(wgths[[j]]*dr[[j]]^2)
  }
  
  lb <- sf1/sf2                           #normalization constant
  aconf <- lb*aconf   
  for (j in 1:m) {                             #normalize X, D, compute stress      
    #aconf <- lb*aconf                     
  	xr[[j]] <- lb*xr[[j]]
    dr[[j]] <- lb*dr[[j]]
	  sold <- sold + sum(wgths[[j]]*(dh[[j]]-dr[[j]])^2)
  }

  #--------------- begin majorization ------------------
  repeat
  {
    br <- list()
    yr <- list()
    er <- list()
    sunc <- 0
    for (j in 1:m) {                          #compute B, Y, 
    	br <- appendList(br,bmat(dh[[j]],wgths[[j]],dr[[j]]))
	    yr <- appendList(yr, vr[[j]] %*% (br[[j]] %*% xr[[j]]))
	    er <- appendList(er,dist(yr[[j]]))
	    sunc <- sunc + sum(wgths[[j]]*(dh[[j]]-er[[j]])^2)
    }
    scon<-sunc

    #--------- impose constraints ---------
    if (!is.null(constr)) {
	     scon <- 0
       er <- list()

       #-- same configurations across ways, configuration weights I
       if (constr=="identity") {
		    z <- matrix(0,n,p)
        u <- matrix(0,n,n)
		    for (j in 1:m) {
			   z<-z+wr[[j]]%*%yr[[j]]
			   u<-u+wr[[j]]
		    }
		    aconf<-myGenInv(u)%*%z
        yr<-repList(aconf,m)
	     }

       #-- configuration weights diagonal INDSCAL
       if (constr=="diagonal") {
		    aux0<-matrix(0,n,p)
		    for (j in 1:m) {
			   aux1<-diag(crossprod(aconf,wr[[j]]%*%yr[[j]]))
			   aux2<-diag(crossprod(aconf,wr[[j]]%*%aconf))
			   bconf[[j]]<-diag(aux1/aux2)
			   aux0<-aux0+(wr[[j]]%*%yr[[j]]%*%bconf[[j]])		    			    
		    } 
		    for (s in 1:p) {
			   aux1<-matrix(0,n,n)
			   for (j in 1:m) aux1<-aux1+(bconf[[j]][s,s]^2)*wr[[j]]
				 aconf[,s]<-myGenInv(aux1)%*%aux0[,s]
		  	}
		    for (j in 1:m) yr[[j]]<-aconf%*%bconf[[j]]
			 }

       #-- no constraints, idioscal ---- 
       if (constr=="idioscal") {
		    aux0<-matrix(0,n,p); auxk<-matrix(0,(n*p),(n*p))
		    for (j in 1:m) {
			    aux1<-crossprod(aconf,wr[[j]]%*%yr[[j]])
			    aux2<-crossprod(aconf,wr[[j]]%*%aconf)
			    auxb<-solve(aux2,aux1)
          bconf[[j]]<-auxb
			    auxc<-crossprod(t(auxb))
			    aux0<-aux0+(wr[[j]]%*%yr[[j]]%*%t(auxb))	
			    auxk<-auxk+kronecker(auxc,wr[[j]])    			    
		    }
		    auxv<-kronecker(diag(p),matrix((1/n),n,n))
		    aconf<-matrix(solve(auxk+auxv,as.vector(aux0)),n,p)
		    for (j in 1:m)  yr[[j]]<-aconf%*%bconf[[j]]
	     }

      for (j in 1:m) {
		    er <- appendList(er,dist(yr[[j]]))
		    scon <- scon+sum(wgths[[j]]*(dh[[j]]-er[[j]])^2) #constraint stress computation
      }
    }
    #-------- end constraints -----------
    
    snon <- scon   
    
    #dhat2 <- transform(e, disobj, w = wgths, normq = nn)  ## dhat update
    #dh <- dhat2$res
    snon<-0
    dh <- list()
    for (j in 1:m) {
      do <- transform(er[[j]], disobj[[j]], w = wgths[[j]], normq = nn)$res
      dh <- appendList(dh, do)
      snon <- snon+sum(wgths[[j]]*(dh[[j]]-er[[j]])^2)
    }
    #-------- nonmetric MDS ----------
#     if (type == "ordinal") {
# 	    if ((itel%%modulus) == 0) {
#         snon<-0
#         dh<-list()
# 	       for(j in 1:m) {
# 		            ds <- diss[[j]]
#                 es <- er[[j]]
#                 ws <- wgths[[j]]
# 	            	if (ties=="primary") do <- monregP(ds,es,ws)
# 	            	if (ties=="secondary") do <- monregS(ds,es,ws)
# 	            	if (ties=="tertiary") do <- monregT(ds,es,ws)
# 		            dh <- appendList(dh,normDissN(do,ws,1))
# 		            snon <- snon+sum(ws*(dh[[j]]-es)^2)
# 	       }
#       }
#     }
#     
#     if (type == "interval") {
#       snon<-0
#       dh<-list()
#       for(j in 1:m) {
#         ds <- diss[[j]]
#         es <- er[[j]]
#         ws <- wgths[[j]]
#         Amat <- cbind(1, as.vector(ds), as.vector(ds)^2) 
#         do <- nnlsPred(Amat, as.vector(es), as.vector(ws))$pred
#         dh <- appendList(dh,normDissN(do,ws,1))
#         snon <- snon+sum(ws*(dh[[j]]-es)^2)
#       }
#     }
#     
    if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress (not normalized): ", formatC(c(snon),digits=8,width=12,format="f"),"\n")

    if (((sold-snon)<eps) || (itel == itmax)) break()         #convergence
    
    xr <- yr
    dr <- er
    sold <- snon
    itel <- itel+1
  }
  #------------ end majorization -----------------
  
  
  names(dh) <- names(er) <- names(yr) <- names(delta)
  cnames <- paste("D", 1:p, sep = "")
  for (i in 1:length(yr)) {
    colnames(yr[[i]]) <- cnames
    rownames(yr[[i]]) <- labels(diss[[i]])
    rownames(bconf[[i]]) <- colnames(bconf[[i]]) <- cnames
    dh[[i]] <- structure(dh[[i]], Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
    attr(dh[[i]], "Labels") <- attr(er[[i]], "Labels") <- labels(diss[[i]])
  }
  
  
  colnames(aconf) <- cnames
  rnames <- rownames(as.matrix(delta[[1]]))
  rownames(aconf) <- rnames
  names(bconf) <- names(dh)
  
  snon <- (snon/m)/nn                   #stress normalization  nn <- n*(n-1)/2, m number lists
  stress <- sqrt(snon)
  
  
  confdiss <- rep(list(NULL), m)
  for (j in 1:m) {                              #initialize weights, V, norm d as lists
	 confdiss[[j]] <- normDissN(er[[j]], wgths[[j]], 1)
  }

    
  ## stress-per-point 
  spoint <- list()
  sps <- matrix(0, m, n)
  rss <- NULL
  for (j in 1:m) {
     spointj <- spp(dh[[j]], confdiss[[j]], wgths[[j]])    
     sps[j,] <- spointj$spp                                   ## SPP per subject
     rss[j] <- sum(spointj$resmat[lower.tri(spointj$resmat)]) ## RSS per subject    
  }   
  colnames(sps) <- rnames
  rownames(sps) <- names(delta)
  spp <- colMeans(sps)                                      ## SPP
  rss <- sum(rss)                                           ## total RSS
  
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
  
  #return configurations, configuration distances, normalized observed distances 
  result <- list(delta = diss, dhat = dh, confdiss = confdiss, conf = yr, gspace = aconf, cweights = bconf,
                 stress = stress, spp = spp, weightmat = wgths, resmat = spoint$resmat, rss = rss, sps = sps, ndim = p, 
                 model = "Three-way SMACOF", niter = itel, nobj = n, type = type, call = match.call()) 
  class(result) <- "smacofID"
  result 
}



