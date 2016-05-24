#' Individual Differences Scaling of provenance data
#'
#' Performs 3-way Multidimensional Scaling analysis using Carroll and
#' Chang (1970)'s INdividual Differences SCALing method as implemented
#' using De Leeuw and Mair (2011)'s stress majorization algorithm.
#' @param ... a sequence of datasets of class \code{distributional} or
#' \code{compositional}
#' @param type is either "ratio" or "ordinal"
#' @return an object of class \code{INDSCAL}, i.e. a list containing
#' the following items:
#' 
#' delta: Observed dissimilarities
#'
#' obsdiss: List of observed dissimilarities, normalized
#'
#' confdiss: List of configuration dissimilarities
#'
#' conf: List of matrices of final configurations
#' 
#' gspace: Joint configurations aka group stimulus space
#' 
#' cweights: Configuration weights
#'
#' stress: Stress-1 value
#'
#' spp: Stress per point
#' 
#' sps: Stress per subject (matrix)
#' 
#' ndim: Number of dimensions
#'
#' model: Type of smacof model
#' 
#' niter: Number of iterations
#'
#' nobj: Number of objects
#' @author
#' Jan de Leeuw and Patrick Mair
#' @references
#' de Leeuw, J., & Mair, P. (2009). Multidimensional scaling using
#' majorization: The R package smacof. Journal of Statistical
#' Software, 31(3), 1-30, < http://www.jstatsoft.org/v31/i03/>
#' @examples
#' data(Namib)
#' plot(indscal(Namib$DZ,Namib$HM))
#' @export
indscal <- function(...,type='ordinal'){
    slist <- list(...)
    dnames <- get.data.names(slist)
    names(slist) <- dnames
    disslist <- getdisslist(slist)
    out <- smacofIndDiff(disslist, constraint = "indscal",type=type)
    class(out) <- "INDSCAL"
    return(out)
}

smacofIndDiff <- function (delta, ndim = 2, type = c("ratio", "ordinal"), 
    constraint = c("indscal", "idioscal", "identity"), weightmat = NULL, 
    init = NULL, ties = "primary", verbose = FALSE, modulus = 1, 
    itmax = 1000, eps = 1e-06) {
    type <- match.arg(type, c("ratio", "ordinal"), 
        several.ok = FALSE)
    constraint <- match.arg(constraint, c("indscal", "idioscal", 
        "identity"), several.ok = FALSE)
    diss <- delta
    p <- ndim
    if (constraint == "indscal") 
        constraint <- "diagonal"
    constr <- constraint
    if (!is.list(diss)) 
        diss <- list(diss)
    if ((is.matrix(diss[[1]])) || (is.data.frame(diss[[1]]))) 
        diss <- lapply(diss, strucprep)
    checkdiss(diss)
    if (is.null(weightmat)) 
        wgths <- initWeights(diss)
    else wgths <- weightmat
    if (!is.list(wgths)) {
        wgths <- list(wgths)
        if (length(wgths) != length(diss)) 
            wgths <- sapply(diss, function(wwr) return(wgths))
    }
    if ((is.matrix(wgths[[1]])) || (is.data.frame(wgths[[1]]))) 
        wgths <- lapply(wgths, strucprep)
    n <- attr(diss[[1]], "Size")
    if (p > (n - 1)) 
        stop("Maximum number of dimensions is n-1!")
    nn <- n * (n - 1)/2
    m <- length(diss)
    itel <- 1
    if (is.null(attr(diss[[1]], "Labels"))) {
        for (i in 1:m) attr(diss[[i]], "Labels") <- paste(1:n)
    }
    dr <- list()
    wr <- list()
    vr <- list()
    dh <- list()
    for (j in 1:m) {
        wr <- appendList(wr, vmat(wgths[[j]]))
        vr <- appendList(vr, myGenInv(wr[[j]]))
        dh <- appendList(dh, normDissN(diss[[j]], wgths[[j]], 
            1))
    }
    xr <- list()
    sold <- sf1 <- sf2 <- 0
    aconf <- initConf(init, diss, n, p, inddiff = TRUE)
    bconf <- repList(diag(p), m)
    for (j in 1:m) {
        xr[[j]] <- aconf %*% bconf[[j]]
        dr[[j]] <- stats::dist(xr[[j]])
        sf1 <- sf1 + sum(wgths[[j]] * dr[[j]] * dh[[j]])
        sf2 <- sf2 + sum(wgths[[j]] * dr[[j]]^2)
    }
    lb <- sf1/sf2
    aconf <- lb * aconf
    for (j in 1:m) {
        xr[[j]] <- lb * xr[[j]]
        dr[[j]] <- lb * dr[[j]]
        sold <- sold + sum(wgths[[j]] * (dh[[j]] - dr[[j]])^2)
    }
    repeat {
        br <- list()
        yr <- list()
        er <- list()
        sunc <- 0
        for (j in 1:m) {
            br <- appendList(br, bmat(dh[[j]], wgths[[j]], dr[[j]]))
            yr <- appendList(yr, vr[[j]] %*% br[[j]] %*% xr[[j]])
            er <- appendList(er, stats::dist(yr[[j]]))
            sunc <- sunc + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2)
        }
        scon <- sunc
        if (!is.null(constr)) {
            scon <- 0
            er <- list()
            if (constr == "identity") {
                z <- matrix(0, n, p)
                u <- matrix(0, n, n)
                for (j in 1:m) {
                  z <- z + wr[[j]] %*% yr[[j]]
                  u <- u + wr[[j]]
                }
                aconf <- myGenInv(u) %*% z
                yr <- repList(aconf, m)
            }
            if (constr == "diagonal") {
                aux0 <- matrix(0, n, p)
                for (j in 1:m) {
                  aux1 <- diag(crossprod(aconf, wr[[j]] %*% yr[[j]]))
                  aux2 <- diag(crossprod(aconf, wr[[j]] %*% aconf))
                  bconf[[j]] <- diag(aux1/aux2)
                  aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% bconf[[j]])
                }
                for (s in 1:p) {
                  aux1 <- matrix(0, n, n)
                  for (j in 1:m) aux1 <- aux1 + (bconf[[j]][s, 
                    s]^2) * wr[[j]]
                  aconf[, s] <- myGenInv(aux1) %*% aux0[, s]
                }
                for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
            }
            if (constr == "idioscal") {
                aux0 <- matrix(0, n, p)
                auxk <- matrix(0, (n * p), (n * p))
                for (j in 1:m) {
                  aux1 <- crossprod(aconf, wr[[j]] %*% yr[[j]])
                  aux2 <- crossprod(aconf, wr[[j]] %*% aconf)
                  auxb <- solve(aux2, aux1)
                  bconf[[j]] <- auxb
                  auxc <- crossprod(t(auxb))
                  aux0 <- aux0 + (wr[[j]] %*% yr[[j]] %*% t(auxb))
                  auxk <- auxk + kronecker(auxc, wr[[j]])
                }
                auxv <- kronecker(diag(p), matrix((1/n), n, n))
                aconf <- matrix(solve(auxk + auxv, as.vector(aux0)), 
                  n, p)
                for (j in 1:m) yr[[j]] <- aconf %*% bconf[[j]]
            }
            for (j in 1:m) {
                er <- appendList(er, stats::dist(yr[[j]]))
                scon <- scon + sum(wgths[[j]] * (dh[[j]] - er[[j]])^2)
            }
        }
        snon <- scon
        if (type == "ordinal") {
            if ((itel%%modulus) == 0) {
                snon <- 0
                dh <- list()
                for (j in 1:m) {
                  ds <- diss[[j]]
                  es <- er[[j]]
                  ws <- wgths[[j]]
                  if (ties == "primary") 
                    do <- monregP(ds, es, ws)
                  if (ties == "secondary") 
                    do <- monregS(ds, es, ws)
                  if (ties == "tertiary") 
                    do <- monregT(ds, es, ws)
                  dh <- appendList(dh, normDissN(do, ws, 1))
                  snon <- snon + sum(ws * (dh[[j]] - es)^2)
                }
            }
        }
        #if (type == "interval") {
        #    snon <- 0
        #    dh <- list()
        #    for (j in 1:m) {
        #        ds <- diss[[j]]
        #        es <- er[[j]]
        #        ws <- wgths[[j]]
        #        Amat <- cbind(1, as.vector(ds), as.vector(ds)^2)
        #        do <- nnlsPred(Amat, as.vector(es), as.vector(ws))$pred
        #        dh <- appendList(dh, normDissN(do, ws, 1))
        #        snon <- snon + sum(ws * (dh[[j]] - es)^2)
        #    }
        #}
        if (verbose) 
            cat("Iteration: ", formatC(itel, width = 3, format = "d"), 
                " Stress (not normalized): ", formatC(c(snon), 
                  digits = 8, width = 12, format = "f"), "\n")
        if (((sold - snon) < eps) || (itel == itmax)) 
            (break)()
        xr <- yr
        dr <- er
        sold <- snon
        itel <- itel + 1
    }
    names(dh) <- names(er) <- names(yr) <- names(delta)
    cnames <- paste("D", 1:p, sep = "")
    for (i in 1:length(yr)) {
        colnames(yr[[i]]) <- cnames
        rownames(yr[[i]]) <- labels(diss[[i]])
        rownames(bconf[[i]]) <- colnames(bconf[[i]]) <- cnames
        dh[[i]] <- structure(dh[[i]], Size = n, call = quote(as.dist.default(m = b)), 
            class = "dist", Diag = FALSE, Upper = FALSE)
        attr(dh[[i]], "Labels") <- attr(er[[i]], "Labels") <- labels(diss[[i]])
    }
    colnames(aconf) <- cnames
    rownames(aconf) <- labels(diss[[1]])
    names(bconf) <- names(dh)
    snon <- (snon/m)/nn
    stress <- sqrt(snon)
    confdiss <- rep(list(NULL), m)
    for (j in 1:m) {
        confdiss[[j]] <- normDissN(er[[j]], wgths[[j]], 1)
    }
    spoint <- spp(dh, confdiss, wgths)
    reslist <- mapply(function(ldh, lcd) {
        (as.matrix(ldh) - as.matrix(lcd))^2
    }, dh, confdiss, SIMPLIFY = FALSE)
    sps <- sapply(reslist, mean)
    if (itel == itmax) 
        warning("Iteration limit reached! Increase itmax argument!")
    result <- list(delta = diss, dhat = dh, confdiss = confdiss, 
        conf = yr, gspace = aconf, cweights = bconf, stress = stress, 
        spp = spoint$spp, weightmat = wgths, resmat = spoint$resmat, 
        sps = sps, ndim = p, model = "Three-way SMACOF", niter = itel, 
        nobj = n, type = type, call = match.call())
    class(result) <- "smacofID"
    result
}

## sanity checks for dissimilatiries (>=0)
checkdiss <- function(diss) {
  if (any(sapply(diss, function(d0) any(d0 < 0, na.rm = TRUE))))
      stop("Dissimilarities should be non-negative!")
}

initWeights <- function(diss) {
  if (!is.list(diss)) {
	  n <- attr(diss,"Size")
    ww <- matrix(1, n, n)
    ww[is.na(as.matrix(diss))] <- 0 ## blank out missings
	  return(stats::as.dist(ww))
  } else {
  n <- attr(diss[[1]],"Size")
  m <- length(diss)
  ww <- repList(matrix(1,n,n),m)
  for (i in 1:m) {
    wwi <- ww[[i]]
    wwi[is.na(as.matrix(diss[[i]]))] <- 0
    ww[[i]] <- stats::as.dist(wwi)
  }
  return(ww)
  }
}

#used in initWeights()
repList<-function(x,n) {
  z <- list()
  for (i in 1:n)
    z<-c(z,list(x))
  return(z)
}

appendList<-function(x,a) {
return(c(x,list(a)))
}

`vmat` <- function(wgths) {
  v <- as.matrix(wgths)
  r <- rowSums(v)  #row margins of weight matrix
  return(diag(r)-v)
}

`myGenInv` <- function(x) {
  n <- dim(x)[1]
  nn <-1/n
  return(solve(x+nn)-nn)
}

normDissN <- function(diss,wghts,m) {
  N <- length(diss)*m
  dissnorm <- diss/sqrt(sum(wghts*diss^2, na.rm = TRUE))*sqrt(N)
  return(dissnorm)
}

initConf <- function(init, diss, n, p, inddiff = FALSE) { 
  if (inddiff) diss <- stats::as.dist(apply(simplify2array(lapply(diss, as.matrix)),
                                     c(1,2), sum, na.rm = TRUE))
    
  if (is.null(init)) {
    meandiss <- mean(diss, na.rm = TRUE)   ## mean dissimilarity for imputation
    diss1 <- as.matrix(diss)
    diss1[is.na(diss1)] <- meandiss
    x <- torgerson(diss1, p = p) 
    init <- "dummy"
  }
  if (init[1] == "random") {
    x <- matrix(stats::runif(n*p, min = -1), ncol = p)
  } 
  if (is.matrix(init)) {
    x <- as.matrix(init)
    if (any(dim(x) != c(n,p)))
        stop(paste0("Dimension of the starting configuration matrix needs to be ",
                    n, " times ", p, "!"))
  }
  return(x)
}

`bmat` <- function(diss, wgths, d, eps = 1E-12) {
  z <- ifelse(d < eps, 1, 0)
  b <- as.matrix((wgths*diss*(1-z))/(d+z))
  r <- rowSums(b) 
  return(diag(r)-b)
}

# primary approach to ties
`monregP` <- function(x, y, w = rep(1,length(x)), block = stats::weighted.mean) {
#x ... observed distances (diss)
#y ... Guttman transformed distances
#w ... weights

  o <- order(x,y)
  r <- order(o)
  return(pavasmacof(y[o],w[o])[r])
}

`pavasmacof` <- function(x, w = rep(1,length(x)), block=stats::weighted.mean) {

#---------------------- subroutines ---------------------------
  merge.block.up <- function(blocklist, blockvalues, x, w, i, block) {
    n <- length(blockvalues)
    nn <- 1:n
    ii <- which(i+1!=nn)
    blocklist[i,] <- c(blocklist[i,1],blocklist[i+1,2])
    indi <- blocklist[i,1]:blocklist[i+1,2]
    blockvalues[i] <- block(x[indi],w[indi])
    blocklist <- blocklist[ii,]
    if (length(ii) == 1) dim(blocklist)<-c(1,2)
    blockvalues <- blockvalues[ii]
    list(v = blockvalues, l = blocklist)
  }

  put.back <- function(n, blocklist, blockvalues) {
    x <- rep(0,n)
    nb <- length(blockvalues)
    for (i in 1:nb) x[blocklist[i,1]:blocklist[i,2]] <- blockvalues[i]
    return(x)
  }

  is.up.satisfied <- function(x,i) (i == length(x))||(x[i]<=x[i+1])
  is.down.satisfied <- function(x,i) (i == 1)||(x[i-1]<=x[i])

  weighted.median <- function(x,w=rep(1,length(x))) {
    ox <- order(x)
    x <- x[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0,w))
    up<-sum(w)-low
    df<-low-up
    repeat {
      if (df[k] < 0) k<-k+1
	else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
	       else return(x[k-1])
    }
  }

  weighted.pth.fractile <- function(x,w=rep(1,length(x)),a=1,b=1) {
    ox <- order(x)
    x <- x[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0,w))
    up <- sum(w)-low
    df <- a*low-b*up
    repeat {
	if (df[k] < 0) k<-k+1
	  else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
		 else return(x[k-1])
    }
  }

#-------------------- end subroutines -------------------------

nblock <- n <-length(x) 
blocklist<-array(1:n,c(n,2)) 
blockvalues<-x
active<-1
repeat{	if (!is.up.satisfied(blockvalues,active)) {
		blockmerge<-merge.block.up(blocklist,blockvalues,x,w,active,block)
		blockvalues<-blockmerge$v; blocklist<-blockmerge$l
		nblock<-nblock-1
		while (!is.down.satisfied(blockvalues,active)) {
			blockmerge<-merge.block.up(blocklist,blockvalues,x,w,active-1,block)
			blockvalues<-blockmerge$v; blocklist<-blockmerge$l; 
			nblock<-nblock-1; active<-active-1;
			}
		}
	else if (active == nblock) break() else active<-active+1
	}	
put.back(n,blocklist,blockvalues)

}

## compute stress per point
spp <- function(dhat, confdiss, wgths) {
  if (!is.list(dhat)) {                       ## all except smacofIndDiff
    resmat <- as.matrix(wgths)*as.matrix(dhat - confdiss)^2    #point stress
    diag(resmat) <- NA
    spp <- colMeans(resmat, na.rm = TRUE)
    spp <- spp/sum(spp)*100
    names(spp) <- colnames(resmat) <- rownames(resmat) <- attr(dhat, "Labels") 
  } else {                                   ## smacofIndDiff
    resmat <- as.matrix(sumList(dhat) - sumList(confdiss))^2 ##point stress
    diag(resmat) <- NA                    
    spp <- colMeans(resmat, na.rm = TRUE)
    spp <- spp/sum(spp)*100
    names(spp) <- colnames(resmat) <- rownames(resmat) <- attr(dhat[[1]], "Labels") 
  }
  return(list(spp = spp, resmat = resmat))
}

#sum of distances over 
sumList <- function(x) {
  m <- length(x)
  z <- x[[1]]
  if (m == 1) return(z)
  for (j in 2:m) z <- z+x[[j]]
  return(z)
}

torgerson <- function(delta, p = 2) {  
#diss ... dissimilarity matrix
#p ... number of dimensions
#------------------- begin subroutines -------------------
#Torgerson's double centering
  doubleCenter <- function(x) {
    n <- dim(x)[1]
    m <- dim(x)[2]
    s <- sum(x)/(n*m)
    xr <- rowSums(x)/m
    xc <- colSums(x)/n
    return((x-outer(xr,xc,"+"))+s)
  }
#-------------------- end subroutines --------------------
  diss <- as.matrix(delta)
  
  z <- eigen(-doubleCenter(as.matrix(diss)^2)/2,symmetric=TRUE)
  v <- pmax(z$values,0)
  if (p == 1) normdiag <- cbind(sqrt(v[1])) else normdiag <- diag(sqrt(v[1:p]))
  conf <- z$vectors[,1:p]%*%normdiag
  rownames(conf) <- rownames(diss)
  return(conf)
}

`strucprep` <- function(x) {
  distvec <- as.vector(x[lower.tri(x)])
  n <- dim(x)[1]
  dissim <- structure(distvec, Size = n, call = quote(as.dist.default(m=b)),
                      class = "dist", Diag = FALSE, Upper = FALSE)
}

# secondary approach to ties
`monregS` <- function(x, y, w = rep(1,length(x)), block = stats::weighted.mean) {
#x ... observed distance matrix
#y ... Guttman transformed distances
#w ... weights

  wag <- tapply(w,x,sum)
  yag <- tapply(y,x,mean)
  xag <- tapply(x,x,mean)
  o <- order(xag)
  r <- order(o)
  e <- pavasmacof(yag[o],wag[o])[r] #call pava
  return(ifelse(outer(x,xag,"=="),1,0)%*%e)
}

#tertiary approach to ties
`monregT` <- function(x,y,w=rep(1,length(x)),block = stats::weighted.mean) {
#x ... observed distance matrix
#y ... Guttman transformed distances
#w ... weights
  wag <- tapply(w,x,sum)
  yag <- tapply(y,x,mean)
  xag <- tapply(x,x,mean)
  o <- order(xag)
  r <- order(o)
  e <- pavasmacof(yag[o],wag[o])[r]  #call pava
  return(y+ifelse(outer(x,xag,"=="),1,0)%*%(e-yag[o]))
}
