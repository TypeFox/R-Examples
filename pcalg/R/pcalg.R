#### Main functions pc(), fci(), fciPlus(), ...
####                ----  -----
### Auxiliaries       -->  ./Aaux.R
### Classes & Methods -->  ./AllClasses.R

trueCov <- function(dag, back.compatible = FALSE)
{
  ##  as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
  ## weightMatrix <- t(as(dag,"matrix"))
  wm <- if(back.compatible) wgtMatrix.0(dag) else wgtMatrix(dag)
  p <- length(dag@nodes)
  ## SS'  where S = (I - W)^{-1} :
  tcrossprod(solve(diag(p) - wm))
}

## buggy randomDAG:
## randomDAG <- function (n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
## {
##     stopifnot(n >= 2, is.numeric(prob), length(prob) == 1,
## 	      0 <= prob, prob <= 1,
## 	      is.numeric(lB), is.numeric(uB), lB <= uB)
##     edL <- vector("list", n)
##     nmbEdges <- 0L
##     for (i in seq_len(n - 2)) {
##         listSize <- rbinom(1, n - i, prob)
##         nmbEdges <- nmbEdges + listSize
##         edgeList <- sample(seq(i + 1, n), size = listSize)
##         weightList <- runif(length(edgeList), min = lB, max = uB)
##         edL[[i]] <- list(edges = edgeList, weights = weightList)
##     }
##     if (nmbEdges > 0) {
## 	edL[[n-1]] <-
## 	    if (rbinom(1, 1, prob) == 1)
## 		list(edges = n,
## 		     weights = runif(1, min = lB, max = uB))
## 	    else
## 		list(edges = integer(0), weights = numeric(0))
## 	edL[[n]] <- list(edges = integer(0), weights = numeric(0))
## 	names(edL) <- V
## 	new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
##     }
##     else
## 	new("graphNEL", nodes = V, edgemode = "directed")
## }

randomDAG <- function (n, prob, lB = 0.1, uB = 1, V = as.character(1:n))
{
    stopifnot(n >= 2, is.numeric(prob), length(prob) == 1,
	      0 <= prob, prob <= 1,
	      is.numeric(lB), is.numeric(uB), lB <= uB)
    edL <- vector("list", n)
    nmbEdges <- 0L
    for (i in seq_len(n - 2)) {
        listSize <- rbinom(1, n - i, prob)
        nmbEdges <- nmbEdges + listSize
        edgeList <- sample(seq(i + 1, n), size = listSize)
        weightList <- runif(length(edgeList), min = lB, max = uB)
        edL[[i]] <- list(edges = edgeList, weights = weightList)
    }
    ## i=n-1 separately
    ## (because of sample(7,1) is actually sample(1:7,1) and not 7)
    listSize <- rbinom(1, 1, prob)
    if (listSize > 0) {
        nmbEdges <- nmbEdges + 1
        edgeList <- n
        weightList <- runif(1, min = lB, max = uB)
    } else {
          edgeList <- integer(0)
          weightList <- numeric(0)
      }
    edL[[n-1]] <- list(edges = edgeList, weights = weightList)
    if (nmbEdges > 0) {
	edL[[n]] <- list(edges = integer(0), weights = numeric(0))
	names(edL) <- V
	new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
    }
    else
	new("graphNEL", nodes = V, edgemode = "directed")
}


## A version of this is also in	/u/maechler/R/MM/Pkg-ex/graph/weightmatrix.R
## another on in  Matrix/R/sparseMatrix.R  function graph.wgtMatrix() :
## No longer in use __apart__ for  rmvDAG(..., back.compatible=TRUE)
wgtMatrix.0 <- function(g, transpose = TRUE)
{
  ## Purpose: work around "graph" package's  as(g, "matrix") bug
  ## ----------------------------------------------------------------------
  ## ACHTUNG: mat_[i,j]==1 iff j->i,
  ## whereas with as(g,"matrix") mat_[i,j]==1 iff i->j
  ## ----------------------------------------------------------------------
  ## Arguments: g: an object inheriting from (S4) class "graph"
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, based on Seth Falcon's code;  Date: 12 May 2006

  ## MM: another buglet for the case of  "no edges":
  if(numEdges(g) == 0) {
    p <- length(nd <- nodes(g))
    return( matrix(0, p,p, dimnames = list(nd, nd)) )
  }
  ## Usual case, when there are edges:
  if(!("weight" %in% names(edgeDataDefaults(g))))
    edgeDataDefaults(g, "weight") <- 1L
  w <- unlist(edgeData(g, attr = "weight"))
  ## we need the *transposed* matrix typically:
  tm <- if(transpose) t(as(g, "matrix")) else as(g, "matrix")
  ## now is a 0/1 - matrix (instead of 0/wgts) with the 'graph' bug
  if(any(w != 1)) ## fix it
    tm[tm != 0] <- w
  ## tm_[i,j]==1 iff i->j
  tm
}

wgtMatrix <- function(g, transpose = TRUE) {
  res <- as(g, "matrix") # from 'graph' package, now reliable (we hope)
  if (transpose) ## default!
      t(res) else res
}

rmvDAG <-
  function(n, dag,
	   errDist = c("normal", "cauchy", "t4", "mix", "mixt3", "mixN100"),
	   mix = 0.1, errMat = NULL, back.compatible = FALSE,
           use.node.names = !back.compatible)
{
  ## Purpose: Generate data according to a given DAG (with weights) and
  ## given node distribution (rows: number of samples; cols: node values in
  ## topological order)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n     : Number of samples
  ## - dag   : Graph object containing the DAG and weights
  ## - errDist: "normal" or "mix" for pure standard normal node distribution
  ##           or mixing with standard cauchy
  ## - mix   : Percentage of points sampled from standard cauchy
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006;  Martin Maechler

  ## check input &  initialize variables
  stopifnot(is(dag, "graph"),
            (p <- length(nodes(dag))) >= 2)

  ##  as(.,"matrix") now {for some versions of 'graph' pkg} is 0/1
  ## weightMatrix <- t(as(dag,"matrix"))
  weightMatrix <- if(back.compatible) wgtMatrix.0(dag) else wgtMatrix(dag)

  ## check if top. sorted
  nonZeros <- which(weightMatrix != 0, arr.ind = TRUE)
  if (nrow(nonZeros) > 0) {
    if (any(nonZeros[,1] - nonZeros[,2] < 0) || any(diag(weightMatrix) != 0))
      stop("Input DAG must be topologically ordered!")
  }

  errDist <- match.arg(errDist)
  if(grepl("^mix", errDist))
    eMat <- function(outs) { # (n,p)
      X <- c(rnorm(n*p - length(outs)), outs)
      matrix(sample(X), nrow = n)
    }
  if(is.null(errMat)) {
    ## generate errors e_i
    errMat <-
      switch(errDist,
             "normal" = matrix(rnorm  (n*p),  nrow = n),
             "cauchy" = matrix(rcauchy(n*p),  nrow = n),
             "t4" =     matrix(rt(n*p, df = 4), nrow = n),
             "mix"    = eMat(rcauchy(round(mix*n*p))),
             "mixt3"  = eMat(     rt(round(mix*n*p), df = 3)),
             "mixN100"= eMat(  rnorm(round(mix*n*p), sd = 10)))
  }
  else { ## check & use 'errMat' argument:
    stopifnot(!is.null(dim.eM <- dim(errMat)),
              dim.eM == c(n,p), is.numeric(errMat))
  }
  if(use.node.names)
    colnames(errMat) <- nodes(dag) # == colnames(weightMatrix)

  ## compute X matrix X_i
  if (sum(weightMatrix) > 0) {
    X <- errMat
    for (j in 2:p) { ## uses X[*, 1:(j-1)] -- "recursively" !
      ij <- 1:(j-1)
      X[,j] <- X[,j] + X[, ij, drop = FALSE] %*% weightMatrix[j, ij]
    }
    X
  }
  else
    errMat
}


pcSelect <- function(y, dm, alpha, corMethod = "standard",
                     verbose = FALSE, directed = FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## ----------------------------------------------------------------------
  ## Value: List
  ## - G: boolean vector with connected nodes
  ## - zMin: Minimal z values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 27.4.07

  stopifnot((n <- nrow(dm)) >= 1,
            (p <- ncol(dm)) >= 1)
  vNms <- colnames(dm)
  ## cl <- match.call()

  zMin <- c(0,rep.int(Inf,p))
  C <- mcor(cbind(y,dm), method = corMethod)
  cutoff <- qnorm(1 - alpha/2)
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- c(FALSE,rep.int(TRUE,p))
  seq_p <- seq_len(p+1L) # = 1:(p+1)

  done <- FALSE
  ord <- 0
  while (!done && any(G)) {
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G)
    remEdges <- length(ind)
    if(verbose >= 1)
      cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
    for (i in 1:remEdges) {
      if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=",i,"|iMax=",remEdges,"\n")
      y <- 1
      x <- ind[i]

      if (G[x]) {
        nbrsBool <- G
        nbrsBool[x] <- FALSE
        nbrs <- seq_p[nbrsBool]
        ## neighbors of y without itself and x
        length_nbrs <- length(nbrs)

        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)

          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat {
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            z <- zStat(x,y, nbrs[S], C,n)
            if(abs(z) < zMin[x]) zMin[x] <- abs(z)
            if (verbose >= 2)
              cat(paste("x:",vNms[x-1],"y:",(ytmp <- round((p+1)/2)),"S:"),
                  c(ytmp,vNms)[nbrs[S]],paste("z:",z,"\n"))
            if (abs(z) <= cutoff) {
              G[x] <- FALSE
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          } ## {repeat}
        }
      } ## end if( G )
    } ## end for(i ..)
    ord <- ord+1
  } ## end while

  ## return
  list(G = setNames(G[-1L], vNms),
       zMin = zMin[-1])
}## pcSelect

zStat <- function(x,y, S, C, n)
{
  ## Purpose: Fisher's z-transform statistic of partial corr.(x,y | S)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y cond. indep. given S?
  ## - C: Correlation matrix among nodes
  ## - n: Samples used to estimate correlation matrix
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, 22 May 2006; Markus Kalisch

  ## for C use:
  ## dyn.load("/u/kalisch/cCode/pcAlgo/parcorC.so")
  ##  res <- 0

  ##   if (length(S) < 4) {
  r <- pcorOrder(x,y, S, C)
  ##  } else {
  ##    k <- solve(C[c(x,y,S),c(x,y,S)])
  ##    r <- -k[1,2]/sqrt(k[1,1]*k[2,2])
  ##      r <- .C("parcorC",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(dim(C)[1]),as.double(as.vector(C)))[[1]]
  ##  }

  res <- sqrt(n- length(S) - 3) * 0.5*log.q1pm(r)
  if (is.na(res)) 0 else res
}

condIndFisherZ <- function(x,y,S,C,n, cutoff,
                           verbose = isTRUE(getOption("verbose.pcalg.condIFz")))
{
  ## Purpose: Return boolean result on conditional independence using
  ## Fisher's z-transform
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y cond. indep. given S?
  ## - C: Correlation matrix among nodes
  ## - n: Samples used to estimate correlation matrix
  ## - cutoff: Cutoff for significance level for individual
  ##           partial correlation tests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:32

  ## R Variante
  r <- pcorOrder(x,y, S, C)

  ## C Variante
  ## C res <- 0
  ## C p <- dim(C)[1]
  ## C r <- .C("parcor",as.double(res),as.integer(x-1),as.integer(y-1),as.integer(S-1),as.integer(length(S)),as.integer(p),as.double(C))[[1]]

  T <- sqrt(n-length(S)-3)* 0.5*log.q1pm(r)

  ## cat(" (",x,",",y,") | ",S," : T = ",T,"\n", sep='')

  ## MM: T is only NA when 'r' already is (r = +- 1  <==>  T = +- Inf) -- better check there (FIXME?)
  ## is.na(T) <==>  T <- 0 # if problem, delete edge: be conservative
  is.na(T) || abs(T) <= cutoff
}

pcorOrder <- function(i,j, k, C, cut.at = 0.9999999) {
  ## Purpose: Compute partial correlation
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - i,j,k: Partial correlation of i and j given k
  ## - C: Correlation matrix among nodes
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  if (length(k) == 0) {
    r <- C[i,j]
  } else if (length(k) == 1) {
    r <- (C[i, j] - C[i, k] * C[j, k])/sqrt((1 - C[j, k]^2) * (1 - C[i, k]^2))
  } else { ## length(k) >= 2
    PM <- pseudoinverse(C[c(i,j,k), c(i,j,k)])
    r <- -PM[1, 2]/sqrt(PM[1, 1] * PM[2, 2])
  }
  if(is.na(r)) 0 else min(cut.at, max(-cut.at, r))
}


compareGraphs <- function(gl,gt) {
  ## Purpose: Return TPR, FPR and TDR of comparison of two undirected graphs
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gl: Estimated graph (may be directed, but the direction will
  ##       be dropped internally)
  ## - gt: True graph (may be directed, but the direction will
  ##       be dropped internally)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:35

  ## When 'graph' returns the 'weight matrix' again:
  ##   ml <- as(ugraph(gl), "matrix")
  ##   mt <- as(ugraph(gt), "matrix")
  ##   p <- dim(ml)[1]
  ml <- wgtMatrix(ugraph(gl))
  mt <- wgtMatrix(ugraph(gt))
  p <- dim(ml)[2]

  mt[mt != 0] <- rep(1,sum(mt != 0))
  ml[ml != 0] <- rep(1,sum(ml != 0)) ## inserted to fix bug

  ## FPR :=  #{misplaced edges} / #{true gaps}
  diffm <- ml-mt
  nmbTrueGaps <- (sum(mt == 0)-p)/2
  ##  print(list(p=p,sum=sum(mt==0),mt=mt,ml=ml,nmbTrueGaps=nmbTrueGaps,diffm=diffm))
  fpr <- if (nmbTrueGaps == 0) 1 else (sum(diffm > 0)/2)/nmbTrueGaps

  ## TPR := #{correctly found edges} / #{true edges}
  diffm2 <- mt-ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- if (nmbTrueEdges == 0) 0 else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges

  ## TDR := #{correctly found edges} / #{found edges}
  trueEstEdges <- (nmbTrueEdges-sum(diffm2 > 0)/2) ## #{true edges} - #{not detected}
  tdr <-
    if (sum(ml == 1) == 0) { ## no edges detected
      if (trueEstEdges == 0) 1 ## no edges in true graph
      else 0
    } else trueEstEdges/(sum(ml == 1)/2)

  ## return named vector:
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}

getNextSet <- function(n,k,set) {
  ## Purpose: Generate the next set in a list of all possible sets of size
  ##          k out of 1:n;
  ##  Also returns a boolean whether this set was the last in the list.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - n,k: Choose a set of size k out of numbers 1:n
  ## - set: previous set in list
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006, 17:37

  ## chInd := changing Index
  chInd <- k - (zeros <- sum((seq(n-k+1,n)-set) == 0))
  wasLast <- (chInd == 0)
  if (!wasLast) {
    set[chInd] <- s.ch <- set[chInd] + 1
    if (chInd < k)
      set[(chInd+1):k] <- seq(s.ch +1L, s.ch +zeros)
  }
  list(nextSet = set, wasLast = wasLast)
}

mcor <- function(dm, method = c("standard", "Qn", "QnStable",
				"ogkScaleTau2", "ogkQn", "shrink"))
{
  ## Purpose: Compute correlation matrix (perhaps elementwise)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix; rows: samples, cols: variables
  ## - method: "Qn" or "standard" (default) envokes robust (based on Qn
  ##           scale estimator) or standard correlation estimator, respectively.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 27 Jan 2006

  p <- ncol(dm)
  method <- match.arg(method)
  switch(method,
	 "Qn" = {
	   res <- diag(nrow = p)
	   Qn. <- apply(dm, 2L, Qn)
	   for (i in seq_len(p-1)) {
	     for (j in i:p) {
	       qnSum  <- Qn(dm[,i] + dm[,j])
	       qnDiff <- Qn(dm[,i] - dm[,j])
               res[j,i] <- res[i,j] <-
                    max(-1,
                        min(1, (qnSum^2 - qnDiff^2) / (4*Qn.[i]*Qn.[j])))
	     }
	   }
	   res
	 },
         "QnStable" = {
	   res <- diag(nrow = p)
	   Qn. <- apply(dm, 2L, Qn) # Qn(.) for all columns
           ## xQn := each column divided by its Qn(.)
           xQn <- dm / rep(Qn., each=nrow(dm))
	   for (i in seq_len(p-1)) {
             xQn.i <- xQn[,i]
	     for (j in i:p) {
	       qnSum  <- Qn(xQn.i + xQn[,j])
	       qnDiff <- Qn(xQn.i - xQn[,j])
               res[j,i] <- res[i,j] <-
                    max(-1,
                        min(1, (qnSum^2 - qnDiff^2) / (qnSum^2 + qnDiff^2)))
	     }
	   }
	   res
	 },
	 "ogkScaleTau2" = {
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = scaleTau2,
			  weight.fn = hard.rejection)$cov)
	 },
	 "ogkQn" = {
	   cov2cor(covOGK(dm, n.iter = 2, sigmamu = s_Qn,
			  weight.fn = hard.rejection)$cov)
	 },
	 "standard" = cor(dm),
         "shrink" = {
           CM <- cor(dm)
           n <- nrow(dm)
           p <- ncol(dm)
           S1 <- sum(CM[1,-1]^2)
           S2 <- sum((1-CM[1,-1]^2)^2)

           g <- S1 / (S1 + S2/n) # is in [0,1]
           scor3 <- CM
           for(i in 2:p) {
             scor3[1,i] <- g*CM[1,i]
             scor3[i,1] <- g*CM[i,1]
           }
           scor3
         }
         )# {switch}
}

pcSelect.presel <- function(y, dm, alpha, alphapre, corMethod = "standard",
                            verbose = 0, directed = FALSE)
{
  ## Purpose: Find columns in dm, that have nonzero parcor with y given
  ## any other set of columns in dm with use of some kind of preselection
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - y: Response Vector (length(y)=nrow(dm))
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - alpha   : Significance level of individual partial correlation tests
  ## - alphapre: Significance level in preselective use of pcSelect
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## ----------------------------------------------------------------------
  ## Author: Philipp Ruetimann, Date: 5 Mar 2008

  tmp <- pcSelect(y, dm, alpha=alphapre,
                  corMethod=corMethod, verbose=verbose, directed=directed)
  pcs <- tmp$G
  zmi <- tmp$zMin
  ppcs <- which(pcs == 1L)
  Xnew <- dm[, ppcs, drop=FALSE]
  lang <- length(pcs)
  tmp2 <- pcSelect(y, dm=Xnew, alpha=alpha,
                   corMethod=corMethod, verbose=verbose, directed=directed)
  pcSnew <- tmp2$G
  zminew <- tmp2$zMin
  zmi[ppcs] <- zminew
  ## MM FIXME -- do without for() :
  k <- 1
  for (i in 1:lang) {
    if(pcs[i] == 1) {
      pcs[i] <- pcSnew[k]
      k <- k+1
    }
  }
  list(pcs = pcs, Xnew = Xnew, zMin = zmi)
}


corGraph <- function(dm, alpha = 0.05, Cmethod = "pearson")
{
  ## Purpose: Computes a correlation graph. Significant correlations are
  ## shown using the given correlation method.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix with rows as samples and cols as variables
  ## - alpha: Significance level for correlation test
  ## - Cmethod: a character string indicating which correlation coefficient
  ##          is to be  used for the test.  One of '"pearson"',
  ##          '"kendall"', or '"spearman"', can be abbreviated.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 30 Jan 2006; Martin Maechler (preserve cns)

  stopifnot(is.numeric(p <- ncol(dm)), p >= 2)
  if(is.null(cns <- colnames(dm))) cns <- as.character(1:p)
  mat <- matrix(0, p,p)
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      mat[i,j] <- cor.test(dm[,i], dm[,j], alternative = "two.sided",
                           method = Cmethod)$p.value < alpha
    }
  }
  mat <- mat + t(mat)
  dimnames(mat) <- list(cns,cns)
  as(mat, "graphNEL")
}

##################################################
## CPDAG
##################################################

##################################################
## dag2cpdag
##################################################

dag2cpdag <- function(g)
{
  ## Purpose: Compute the (unique) completed partially directed graph (CPDAG)
  ## that corresponds to the input DAG; result is a graph object
  ## In the current implementation, this function is just a wrapper
  ## function for 'dag2essgraph'
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dag: input DAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Alain Hauser, Date: 13 Mar 2015
  nn <- nodes(g) ## *: to keep node labels
  dag <- as(g, "GaussParDAG")
  res <- as(dag2essgraph(dag), "graphNEL")
  nodes(res) <- nn ## *
  res
}


#dag2cpdag <- function(g)
#{
#  ## Purpose: Compute the (unique) completed partially directed graph (CPDAG)
#  ## that corresponds to the input DAG; result is a graph object
#  ## ----------------------------------------------------------------------
#  ## Arguments:
#  ## - dag: input DAG (graph object)
#  ## ----------------------------------------------------------------------
#  ## Author: Diego Colombo, Date: 10 Jun 2013, 11:06
#
#    amat <- as(g, "matrix")
#    amat[amat != 0] <- 1
#    skel.amat <- amat + t(amat)
#    skel.amat[skel.amat == 2] <- 1
#    cpdag <- skel.amat
#
#    ## search the v-structures in the DAG
#    ind <- which((amat == 1 & t(amat) == 0), arr.ind = TRUE)
#    tripleMatrix <- matrix(,0,3)
#    ## Go through all edges
#    for (i in seq_len(nrow(ind))) { ## MM(FIXME): growth of tripleMatrix
#        x <- ind[i,1]
#        y <- ind[i,2]
#        indY <- setdiff(which((amat[,y] == 1 & amat[y,] == 0), arr.ind = TRUE),x) ## x-> y <- z
#	if(length(newZ <- indY[amat[x,indY] == 0])) ## deparse.l.=0: no colnames
#	  tripleMatrix <- rbind(tripleMatrix, cbind(x, y, newZ, deparse.level=0),
#				deparse.level=0)
#    }
#    if ((m <- nrow(tripleMatrix)) > 0) {
#        deleteDupl <- logical(m)# all FALSE
#        for (i in seq_len(m))
#            if (tripleMatrix[i,1] > tripleMatrix[i,3])
#                deleteDupl[i] <- TRUE
#        if(any(deleteDupl))
#          tripleMatrix <- tripleMatrix[!deleteDupl,, drop=FALSE]
#
#        ## orient the v-structures in the CPDAG
#        for (i in seq_len(nrow(tripleMatrix))) {
#                x <- tripleMatrix[i,1]
#                y <- tripleMatrix[i,2]
#                z <- tripleMatrix[i,3]
#                cpdag[x,y] <- cpdag[z,y] <- 1
#                cpdag[y,x] <- cpdag[y,z] <- 0
#        }
#    }
#
#    ## orient the edges with the 3 orientation rules
#    repeat {
#        old_cpdag <- cpdag
#        ## Rule 1
#        ind <- which((cpdag == 1 & t(cpdag) == 0), arr.ind = TRUE)
#        for (i in seq_len(nrow(ind))) {
#            a <- ind[i, 1]
#            b <- ind[i, 2]
#            isC <- ((cpdag[b, ] == 1 & cpdag[, b] == 1) &
#                    (cpdag[a, ] == 0 & cpdag[, a] == 0))
#            if (any(isC)) {
#                indC <- which(isC)
#                cpdag[b, indC] <- 1
#                cpdag[indC, b] <- 0
#            }
#        }
#        ## Rule 2
#        ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
#        for (i in seq_len(nrow(ind))) {
#            a <- ind[i, 1]
#            b <- ind[i, 2]
#            isC <- ((cpdag[a, ] == 1 & cpdag[, a] == 0) &
#                    (cpdag[, b] == 1 & cpdag[b, ] == 0))
#            if (any(isC)) {
#                cpdag[a, b] <- 1
#                cpdag[b, a] <- 0
#            }
#        }
#        ## Rule 3
#        ind <- which((cpdag == 1 & t(cpdag) == 1), arr.ind = TRUE)
#        for (i in seq_len(nrow(ind))) {
#            a <- ind[i, 1]
#            b <- ind[i, 2]
#            indC <- which((cpdag[a, ] == 1 & cpdag[, a] == 1) &
#                          (cpdag[, b] == 1 & cpdag[b, ] == 0))
#            if (length(indC) >= 2) {
#                cmb.C <- combn(indC, 2)
#                cC1 <- cmb.C[1, ]
#                cC2 <- cmb.C[2, ]
#                for (j in seq_along(cC1)) {
#                    c1 <- cC1[j]
#                    c2 <- cC2[j]
#                    if (c1 != c2 && cpdag[c1, c2] == 0 && cpdag[c2,c1] == 0) {
#                        cpdag[a, b] <- 1
#                        cpdag[b, a] <- 0
#                        break
#                    }
#                }
#            }
#        }
#        if (all(cpdag == old_cpdag))
#            break
#    }
#    as(cpdag,"graphNEL")
#}


## dag2cpdag <- function(dag) {
  ## Purpose: Compute the (unique) completed partially directed graph (CPDAG)
  ## that corresponds to the input DAG; result is a graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dag: input DAG (graph object)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:30

##  p <- numNodes(dag)
##  ## transform DAG to adjacency matrix if any edges are present
##  if (numEdges(dag)==0) {
##    cpdag.res <- dag
##  } else {
##    dag <- as(dag,"matrix")
##    dag[dag!=0] <- 1

##    ## dag is adjacency matrix
##    e.df <- labelEdges(dag)
##    cpdag <- matrix(0, p,p)
##    for (i in seq_len(nrow(e.df))) {
##      if (e.df$label[i]) {
##        cpdag[e.df$tail[i],e.df$head[i]] <- 1
##      } else {
##        cpdag[e.df$tail[i],e.df$head[i]] <- cpdag[e.df$head[i],e.df$tail[i]] <- 1
##      }
##    }
##    rownames(cpdag) <- colnames(cpdag) <- as.character(seq(1,p))
##    cpdag.res <- as(cpdag,"graphNEL")
##  }
##  cpdag.res
##}

## make.edge.df <- function(amat) {
  ## Purpose: Generate a data frame describing some properties of a DAG
  ## (for extending to a CPDAG)
  ## The output contains xmin,xmax,head,tail,order (NA or number),
  ## type (1="d",0="u") in lexikographic order
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG [x_ij=1 means i->j]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:43

  ## INPUT: Adjacency matrix
##  stopifnot(sum(amat)>0)
##  e <- which(amat==1,arr.ind=TRUE)
##  e.dup <- duplicated(t(apply(e,1,sort)))
##  nmb.edges <- sum(!e.dup)
##  res <- data.frame(xmin=rep(NA,nmb.edges),xmax=rep(NA,nmb.edges),
##                    tail=rep(NA,nmb.edges),head=rep(NA,nmb.edges),
##                    order=rep(NA,nmb.edges),type=rep(1,nmb.edges))
##  pure.edges <- e[!e.dup,]
##  if(length(pure.edges)==2) dim(pure.edges) <- c(1,2)
##  for (i in seq_len(nrow(pure.edges))) {
##    if (all(amat[pure.edges[i,1],pure.edges[i,2]]==
##            amat[pure.edges[i,2],pure.edges[i,1]])) {
##      res$type[i] <- 0
##      res$head[i] <- NA
##      res$tail[i] <- NA
##    } else {
##      res$head[i] <- pure.edges[i,2]
##      res$tail[i] <- pure.edges[i,1]
##    }
##  }
##  s.pure.edges <- t(apply(pure.edges,1,sort))
##  ii <- order(s.pure.edges[,1],s.pure.edges[,2])
##  res <- res[ii,]
##  res$xmin <- s.pure.edges[ii,1]
##  res$xmax <- s.pure.edges[ii,2]
##  res
##}

## orderEdges <- function(amat) {
  ## Purpose: Order the edges of a DAG according to Chickering
  ## (for extension to CPDAG)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:42

##  stopifnot(isAcyclic(amat))
##  ordered.nodes <- topOrder(amat) ## parents before children
##  edge.df <- make.edge.df(amat)

##  eOrder <- 0
##  while(any(unOrdered <- is.na(edge.df$order))) {
##    counter <- 0
    ## find y
##    y <- NA
##    found <- FALSE
##    while(!found) {
##      counter <- counter+1
##      node <- ordered.nodes[counter]
      ## which edges are incident to node?
##      nbr.nodes <- which(amat[,node]==1)
##      if(length(nbr.nodes)>0) {
##        unlabeled <- rep.int(FALSE, length(nbr.nodes))
##        for(i in seq_along(nbr.nodes)) {
##          x <- nbr.nodes[i]
          ## is edge edge x-y unlabeled?
##	  unlabeled[i] <- length(intersect(which(edge.df$xmin==min(node,x) &
##						 edge.df$xmax==max(node,x)),
##					   which(unOrdered))) > 0
##        }
        ## choose unlabeled edge with highest order node
##        if(any(unlabeled)) {
##          nbr.unlab <- nbr.nodes[unlabeled] # nbrnodes w. unlabeled edges
##          tmp <- ordered.nodes[ordered.nodes %in% nbr.unlab]
##          y <- tmp[length(tmp)]
          ## y <- last(ordered.nodes[which(ordered.nodes %in% nbr.unlab)])
##          edge.df$order[edge.df$xmin==min(node,y) &
##                        edge.df$xmax==max(node,y)] <- eOrder
##          eOrder <- eOrder+1
##          found <- TRUE
##        }
##      }

##    } ## while !found

##  } ## while any(unOrdered)
##  edge.df
##}


## labelEdges <- function(amat) {
  ## Purpose: Label the edges in a DAG with "compelled" and "reversible"
  ## (for extension to a CPDAG)
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix of DAG
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;  prettified: MMaechler

  ## label=TRUE -> compelled
  ## label=FALSE -> reversible
##  edge.df <- orderEdges(amat)
##  lab <- rep(NA,dim(edge.df)[1])
##  edge.df <- edge.df[order(edge.df$order),]
##  Head <- edge.df$head
##  Tail <- edge.df$tail

##  while(any(ina <- is.na(lab))) {
##    x.y <- which(ina)[1]
##    x <- Tail[x.y]
##    y <- Head[x.y]
##    y.is.head <- Head == y
##    e1 <- which(Head == x & lab)
##    for(ee in e1) {
##      w <- Tail[ee]
##      if (any(wt.yh <- w == Tail & y.is.head))
##        lab[wt.yh] <- TRUE
##      else {
##        lab[y.is.head] <- TRUE
##        break
##      }
##    }
    ## edges going to y not starting from x
##    cand <- which(y.is.head  &  Tail != x)
##    if (length(cand) > 0) {
##      valid.cand <- rep(FALSE,length(cand))
##      for (iz in seq_along(cand)) {
##        z <- Tail[cand[iz]]
##        if (!any(Tail==z & Head==x)) ## NOT.parent.of.x :
##          valid.cand[iz] <- TRUE
##      }
##      cand <- cand[valid.cand]
##    }
##    lab[which(y.is.head & is.na(lab))] <- (length(cand) > 0)
##  }
##  edge.df$label <- lab
##  edge.df
##}

##################################################
## pdag2dag
##################################################
find.sink <- function(gm) {
  ## Purpose: Find sink of an adj matrix; return numeric(0) if there is none;
  ## a sink may have incident undirected edges, but no directed ones
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix (gm_i_j is edge from j to i)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;  speedup: Martin Maechler, Dec.2013

  ## treat undirected edges
  gm[gm == t(gm) & gm == 1] <- 0
  ## treat directed edges
  which(colSums(gm) == 0)
}

adj.check <- function(gm,x) {
  ## Purpose:  Return "TRUE", if:
  ## For every vertex y, adj to x, with (x,y) undirected, y is adjacent to
  ## all the other vertices which are adjacent to x
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: adjacency matrix of graph
  ## - x: node number (number)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006;
  ## several smart speedups: Martin Maechler, Dec.2013

  gm.1 <- (gm == 1)
  xr <- gm.1[x,]
  xc <- gm.1[,x]
  nx <- which(xr | xc)
  ## undirected neighbors of x
  un <- which(xr & xc)
  for(y in un) {
      adj.x <- setdiff(nx, y)
      adj.y <- setdiff(which(gm.1[y,] | gm.1[,y]), x)
      if(!all(adj.x %in% adj.y))
          return(FALSE)
  }
  TRUE
}


amat2dag <- function(amat) {
  ## Purpose: Transform the adjacency matrix of an PDAG to the adjacency
  ## matrix of a SOME DAG in the equiv. class
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: adjacency matrix; x -> y if amat[x,y]=1,amat[y,x]=0
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 31 Oct 2006, 15:23

  p <- dim(amat)[1]
  ## amat to skel
  skel <- amat+t(amat)
  skel[which(skel > 1)] <- 1

  ## permute skel
  ord <- sample.int(p)
  skel <- skel[ord,ord]

  ## skel to dag
  for (i in 2:p) {
    for (j in 1:(i-1)) {
      if(skel[i,j] == 1) skel[i,j] <- 0
    }
  }
  ## inverse permutation
  i.ord <- order(ord)
  skel[i.ord,i.ord]
}

##################################################
## udag2pdag
##################################################
udag2pdag <- function(gInput, verbose = FALSE) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03

  res <- gInput
  if (numEdges(gInput@graph) > 0) {
    g <- as(gInput@graph,"matrix") ## g_ij if i->j
    p <- as.numeric(dim(g)[1])
    pdag <- g
    ind <- which(g == 1,arr.ind = TRUE)

    ## Create minimal pattern
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,] == 1),x) ## x-y-z
      for (z in allZ) {
        if (g[x,z] == 0  &&
            !(y %in% gInput@sepset[[x]][[z]] ||
              y %in% gInput@sepset[[z]][[x]])) {
          if (verbose) {
            cat("\n",x,"->",y,"<-",z,"\n")
            cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
          }
          pdag[x,y] <- pdag[z,y] <- 1
          pdag[y,x] <- pdag[y,z] <- 0
        }
      }
    }

    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))

    if (res2$success) {
      ## Convert to complete pattern: use rules by Pearl
      old_pdag <- matrix(0, p,p)
      while (!all(old_pdag == pdag)) {
        old_pdag <- pdag
        ## rule 1
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[b,] == 1 & pdag[,b] == 1) & (pdag[a,] == 0 & pdag[,a] == 0))
          if (length(indC) > 0) {
            pdag[b,indC] <- 1
            pdag[indC,b] <- 0
            if (verbose)
              cat("\nRule 1:",a,"->",b," and ",b,"-",indC,
                  " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule1")

        ## rule 2
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[a,] == 1 & pdag[,a] == 0) & (pdag[,b] == 1 & pdag[b,] == 0))
          if (length(indC) > 0) {
            pdag[a,b] <- 1
            pdag[b,a] <- 0
            if (verbose) cat("\nRule 2: Kette ",a,"->",indC,"->",
                  b,":",a,"->",b,"\n")
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule2")

        ## rule 3
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE) ## a - b
        for (i in seq_len(nrow(ind))) {
          a <- ind[i,1]
          b <- ind[i,2]
          indC <- which( (pdag[a,] == 1 & pdag[,a] == 1) & (pdag[,b] == 1 & pdag[b,] == 0))
          if (length(indC) >= 2) {
            ## cat("R3: indC = ",indC,"\n")
            g2 <- pdag[indC,indC]
            ## print(g2)
            if (length(g2) <= 1) {
              g2 <- 0
            } else {
              diag(g2) <- rep(1,length(indC)) ## no self reference
            }
            if (any(g2 == 0)) { ## if two nodes in g2 are not connected
              pdag[a,b] <- 1
              pdag[b,a] <- 0
              if (verbose) cat("\nRule 3:",a,"->",b,"\n")
            }
          }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule3")

        ## rule 4
        ##-         ind <- which((pdag==1 & t(pdag)==1), arr.ind=TRUE) ## a - b
        ##-         if (length(ind)>0) {
        ##-           for (i in seq_len(nrow(ind))) {
        ##-             a <- ind[i,1]
        ##-             b <- ind[i,2]
        ##-             indC <- which( (pdag[a,]==1 & pdag[,a]==1) & (pdag[,b]==0 & pdag[b,]==0))
        ##-             l.indC <- length(indC)
        ##-             if (l.indC>0) {
        ##-               found <- FALSE
        ##-               ic <- 0
        ##-               while(!found & (ic < l.indC)) {
        ##-                 ic <- ic + 1
        ##-                 c <- indC[ic]
        ##-                 indD <- which( (pdag[c,]==1 & pdag[,c]==0) & (pdag[,b]==1 & pdag[b,]==0))
        ##-                 if (length(indD)>0) {
        ##-                   found <- TRUE
        ##-                   pdag[b,a] = 0
        ##-                   if (verbose) cat("Rule 4 applied \n")
        ##-                 }
        ##-               }
        ##-             }
        ##-           }
        ##-         }

      }
      res@graph <- as(pdag,"graphNEL")
    } else {
      ## was not extendable; random DAG chosen
      res@graph <- res2$graph
      ## convert to CPDAG
      res@graph <- dag2cpdag(res@graph)
    }
  }
  return(res)
} ## udag2pdag

shd <- function(g1,g2)
{
  ## Purpose: Compute Structural Hamming Distance between graphs g1 and g2
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g1, g2: Input graphs
  ## (graph objects;connectivity matrix where m[x,y]=1 iff x->1
  ## and m[x,y]=m[y,x]=1 iff x-y; pcAlgo-objects)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  1 Dec 2006, 17:21

                                        ## Idea: Transform g1 into g2
                                        ## Transform g1 and g2 into adjacency matrices
  if (is(g1, "pcAlgo")) g1 <- g1@graph
  if (is(g2, "pcAlgo")) g2 <- g2@graph

  if (is(g1, "graphNEL")) {
    m1 <- wgtMatrix(g1, transpose = FALSE)
    m1[m1 != 0] <- 1
  }
  if (is(g2, "graphNEL")) {
    m2 <- wgtMatrix(g2, transpose = FALSE)
    m2[m2 != 0] <- 1
  }

  shd <- 0
                                        ## Remove superfluous edges from g1
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
                                        ## Add missing edges to g1
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
                                        ## Compare Orientation
  d <- abs(m1-m2)
  ## return
  shd + sum((d + t(d)) > 0)/2
}

################################################################################
## New in V8 ; uses  vcd  package
################################################################################
ci.test <- function(x,y, S = NULL, dm.df) {
  stopifnot(is.data.frame(dm.df), ncol(dm.df) > 1)
  tab <- table(dm.df[,c(x,y,S)])
  if (any(dim(tab) < 2))
    1
  else if (length(S) == 0)
    fisher.test(tab, simulate.p.value = TRUE)$p.value
  else
    vcd::coindep_test(tab,3:(length(S)+2))$p.value
}

pcAlgo <- function(dm = NA, C = NA, n = NA, alpha, corMethod = "standard",
                   verbose = FALSE, directed = FALSE,
                   G = NULL, datatype = 'continuous', NAdelete = TRUE,
                   m.max = Inf, u2pd = "rand", psepset = FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdagu
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo

  .Deprecated(msg = "pcAlgo() is deprecated and only kept for backward compatibility.
 Please use skeleton, pc, or fci instead\n")
  cl <- match.call()

  if (any(is.na(dm))) {
    stopifnot(all(!is.na(C)),!is.na(n), (p <- ncol(C)) > 0)
  } else {
    n <- nrow(dm)
    p <- ncol(dm)
  }
  n <- as.integer(n)

  if (is.null(G)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else if (!(identical(dim(G),c(p,p))))
      stop("Dimensions of the dataset and G do not agree.")

  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  zMin <- matrix(Inf, p,p)
  n.edgetests <- numeric(1)# final length = max { ord}
  done <- FALSE
  ord <- 0

  if (datatype == 'continuous') {
    diag(zMin) <- 0
    if (any(is.na(C))) C <- mcor(dm, method = corMethod)
    cutoff <- qnorm(1 - alpha/2)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord+1] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[,1]), ]
      remEdges <- nrow(ind)
      if(verbose)
        cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
      for (i in 1:remEdges) {
        if(verbose && i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n")
        x <- ind[i,1]
        y <- ind[i,2]
        if (G[y,x]) {
          nbrsBool <- G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE
            S <- seq(length = ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord+1] <- n.edgetests[ord+1]+1
              z <- zStat(x,y, nbrs[S], C,n)
              if (verbose) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("z:",z,"\n"))
              if(abs(z) < zMin[x,y]) zMin[x,y] <- abs(z)
              if (abs(z) <= cutoff) {
                G[x,y] <- G[y,x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if(nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            }
          }
        } ## end if(!done)

      } ## end for(i ..)
      ord <- ord+1
      ##    n.edgetests[ord] <- remEdges
    } ## while

    for (i in 1:(p-1)) {
      for (j in 2:p) {
        zMin[i,j] <- zMin[j,i] <- min(zMin[i,j],zMin[j,i])
      }
    }
  }
  else {
    ##
    ##
    ## DISCRETE DATA ######################################################
    ##
    if (datatype == 'discrete') {
      dm.df <- as.data.frame(dm)
      while (!done && any(G) && ord <= m.max) {
        n.edgetests[ord+1] <- 0
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ## For comparison with C++ sort according to first row
        ind <- ind[order(ind[,1]), ]
        remEdges <- nrow(ind)
        if(verbose)
          cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
        for (i in 1:remEdges) {
          if(verbose) { if(i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n") }
          x <- ind[i,1]
          y <- ind[i,2]
          if (G[y,x]) {
            nbrsBool <- G[,x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord) done <- FALSE
              S <- seq(length = ord)
              repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
                n.edgetests[ord+1] <- n.edgetests[ord+1]+1
                prob <- ci.test(x,y, nbrs[S], dm.df)
                if (verbose) cat("x=",x," y=",y," S=",nbrs[S],":",prob,"\n")
                if (is.na(prob)) prob <- if(NAdelete) 1 else 0
                if(prob >= alpha) { # independent
                  G[x,y] <- G[y,x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                } else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if(nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          } ## end if(!done)

        } ## end for(i ..)
        ord <- ord+1
        ##    n.edgetests[ord] <- remEdges
      } ## while
    } else
      stop("Datatype must be 'continuous' or 'discrete'.")
  }

  if (psepset) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,] == 1),x) ## x-y-z

      for (z in allZ) {
        if (amat[x,z] == 0 &&
            !((y %in% sepset[[x]][[z]]) ||
              (y %in% sepset[[z]][[x]]))) {
          if (verbose >= 2) {
            cat("\n",x,"*->",y,"<-*",z,"\n")
            cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
          }

          ## x o-> y <-o z
          amat[x,y] <- amat[z,y] <- 2

        } ## for
      } ## if
    } ## for

    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,] != 0)) {
        tf1 <- setdiff(reach(x,-1,-1,amat), x)
        for (y in seq_p[amat[x,] != 0]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf) > 0) {
            az <- abs(zStat(x,y,tf,C,n))
            if (az < zMin[x,y]) zMin[x,y] <- az
            if (az <= cutoff) {
              ## delete x-y
              amat[x, y] <- amat[y, x] <- 0
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
            }
            if (verbose >= 2)
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - |z| = ",az,"\n")
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
  } ## end if(psepset)

  if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

  ## transform matrix to graph object (if not deprecated anyway: FIX to use correct node names!)
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = as.character(seq_p))
  } else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    as(G,"graphNEL")
  }

  res <- new("pcAlgo", graph = Gobject,
             call = cl, n = n, max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = zMin)
  if (directed)
    switch (u2pd,
            "rand"    = udag2pdag       (res),
            "retry"   = udag2pdagSpecial(res)$pcObj,
            "relaxed" = udag2pdagRelaxed(res))
  else
    res
} ## {pcAlgo} __ deprecated __

flipEdges <- function(amat,ind) {
  res <- amat
  for (i in seq_len(nrow(ind))) {
    x <- ind[i,]
    res[x[1],x[2]] <- amat[x[2],x[1]]
    res[x[2],x[1]] <- amat[x[1],x[2]]
  }
  res
}

pdag2dag <- function(g, keepVstruct = TRUE) {
  ## Purpose: Generate a consistent extension of a PDAG to a DAG; if this
  ## is not possible, a random extension of the skeleton is returned and
  ## a warning is issued.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g: PDAG (graph object)
  ## - keepVstruct: TRUE - vStructures are kept
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006; tweaks: Martin M

  not.yet <- FALSE
  if (numEdges(g) == 0) {
    graph <- g
  } else {
    gm <- wgtMatrix(g) ## gm_i_j is edge from j to i
    storage.mode(gm) <- "integer"
    gm[which(gm > 0 & gm != 1)] <- 1L

    a <- gm2 <- gm
    cn2 <- colnames(gm2)
    go.on <- TRUE
    while(go.on && length(a) > 1 && sum(a) > 0) {
      not.yet <- TRUE
      sinks <- find.sink(a)
      if (length(sinks) > 0) {
        counter <- 1L
        while(not.yet && counter <= length(sinks)) {
          x <- sinks[counter]
          if (!keepVstruct || adj.check(a,x)) {
            not.yet <- FALSE
            ## orient edges
            inc.to.x <- which(a[,x] == 1L & a[x,] == 1L) ## undirected
            if (length(inc.to.x) > 0) {
              ## map var.names to col pos in orig adj matrix
              ## bug: real.inc.to.x <- as.numeric(rownames(a)[inc.to.x])
              real.inc.to.x <- which(cn2 %in% rownames(a)[inc.to.x])
              real.x        <- which(cn2 %in% rownames(a)[x])
              gm2[real.x, real.inc.to.x] <- 1L
              gm2[real.inc.to.x, real.x] <- 0L
            }
            ## remove x and all edges connected to it
            a <- a[-x,-x]
          }
          counter <- counter + 1L
        }
      }
      go.on <- !not.yet
    }## { while }

    graph <- if (not.yet) {
               ## warning("PDAG not extendible: Random DAG on skeleton drawn")
               as(amat2dag(gm), "graphNEL")
             } else ## success :
               as(t(gm2), "graphNEL")
  }
  list(graph = graph, success = !not.yet)
}

udag2pdagSpecial <- function(gInput, verbose = FALSE, n.max = 100) {
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PDAG using
  ## the rules of Pearl. The output is again a pcAlgo-object. Ambiguous
  ## v-structures are reoriented until extendable or max number of tries
  ## is reached. If still not extendable, a DAG is produced starting from the
  ## current PDAG even if introducing new v-structures.
  ##
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gInput: pcAlgo object
  ## - verbose: 0 - no output, 1 - detailed output
  ## - n.max: Maximal number of tries to reorient v-strucutres
  ## ----------------------------------------------------------------------
  ## Values:
  ## - pcObj: Oriented pc-Object
  ## - evisit: Matrix counting the number of orientation attemps per edge
  ## - xtbl.orig: Is original graph with v-structure extendable
  ## - xtbl: Is final graph with v-structure extendable
  ## - amat0: Adj.matrix of original graph with v-structures
  ## - amat1: Adj.matrix of graph with v-structures after reorienting
  ##          edges from double edge visits
  ## - status:
  ##   0: original try is extendable
  ##   1: reorienting double edge visits helps
  ##   2: orig. try is not extendable; reorienting double visits don't help;
  ##      result is acyclic, has orig. v-structures, but perhaps
  ##      additional v-structures
  ## - counter: Number of reorientation tries until success or max.tries
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Sep 2006, 15:03
  counter <- 0
  res <- gInput
  status <- 0
  p <- length(nodes(res@graph))
  evisit <- amat0 <- amat1 <- matrix(0,p,p)
  xtbl <- xtbl.orig <- TRUE
  if (numEdges(gInput@graph) > 0) {
    g <- as(gInput@graph,"matrix") ## g_ij if i->j
    p <- dim(g)[1]
    pdag <- g
    ind <- which(g == 1,arr.ind = TRUE)
    ## ind <- unique(t(apply(ind,1,sort)))

    ## Create minimal pattern
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(g[y,] == 1),x) ## x-y-z
      for(z in allZ) {
        if ((g[x,z] == 0) &&
            !(y %in% gInput@sepset[[x]][[z]] ||
              y %in% gInput@sepset[[z]][[x]])) {
          if (verbose) {
            cat("\n",x,"->",y,"<-",z,"\n")
            cat("Sxz=",gInput@sepset[[z]][[x]],"Szx=",gInput@sepset[[x]][[z]])
          }
          ## check if already in other direction directed
          if (pdag[x,y] == 0 && pdag[y,x] == 1) {
            evisit[x,y] <- evisit[x,y] + 1
            evisit[y,x] <- evisit[y,x] + 1
          }
          if (pdag[z,y] == 0 && pdag[y,z] == 1) {
            evisit[z,y] <- evisit[z,y] + 1
            evisit[y,z] <- evisit[y,z] + 1
          }
          pdag[x,y] <- pdag[z,y] <- 1
          pdag[y,x] <- pdag[y,z] <- 0
        } ## if
      } ## for
    } ## for ( i )

    amat0 <- pdag
    ## Test whether this pdag allows a consistent extension
    res2 <- pdag2dag(as(pdag,"graphNEL"))
    xtbl <- res2$success
    xtbl.orig <- xtbl

    if (!xtbl && (max(evisit) > 0)) {
      tmp.ind2 <- unique(which(evisit > 0,arr.ind = TRUE))
      ind2 <- unique(t(apply(tmp.ind2,1,sort)))
      ## print(ind2)
      n <- nrow(ind2)
      n.max <- min(2^n-1,n.max)
      counter <- 0
      ## xtbl is FALSE because of if condition
      while((counter < n.max) & !xtbl) {
        ## if (counter%%100 == 0) cat("\n counter=",counter,"\n")
        counter <- counter + 1
        dgBase <- digitsBase(counter)
        dgBase <- dgBase[length(dgBase):1]
        ## print(dgBase)
        indBase <- matrix(0,1,n)
        indBase[1,seq_along(dgBase)] <- dgBase
        ## indTmp <- ind2[ss[[counter]],,drop=FALSE]
        indTmp <- ind2[(indBase == 1),,drop = FALSE]
        ## print(indTmp)
        pdagTmp <- flipEdges(pdag,indTmp)
        resTmp <- pdag2dag(as(pdagTmp,"graphNEL"))
        xtbl <- resTmp$success
      }
      pdag <- pdagTmp
      status <- 1
    }
    amat1 <- pdag

    if (xtbl) {
      ## Convert to complete pattern: use rules by Pearl
      old_pdag <- matrix(0, p,p)
      while (any(old_pdag != pdag)) {
        old_pdag <- pdag
        ## rule 1
        ind <- which((pdag == 1 & t(pdag) == 0), arr.ind = TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[b,] == 1 & pdag[,b] == 1) & (pdag[a,] == 0 & pdag[,a] == 0))
            if (length(indC) > 0) {
              pdag[b,indC] <- 1
              pdag[indC,b] <- 0
              if (verbose)
                cat("\nRule 1:",a,"->",b," and ",b,"-",indC,
                    " where ",a," and ",indC," not connected: ",b,"->",indC,"\n")
            }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule1")

        ## rule 2
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,] == 1 & pdag[,a] == 0) & (pdag[,b] == 1 & pdag[b,] == 0))
            if (length(indC) > 0) {
              pdag[a,b] <- 1
              pdag[b,a] <- 0
              if (verbose) cat("\nRule 2: Kette ",a,"->",indC,"->",
                    b,":",a,"->",b,"\n")
            }
        }
        ## x11()
        ## plot(as(pdag,"graphNEL"), main="After Rule2")

        ## rule 3
        ind <- which((pdag == 1 & t(pdag) == 1), arr.ind = TRUE) ## a - b
        for (i in seq_len(nrow(ind))) {
            a <- ind[i,1]
            b <- ind[i,2]
            indC <- which( (pdag[a,] == 1 & pdag[,a] == 1) & (pdag[,b] == 1 & pdag[b,] == 0))
            if (length(indC) >= 2) {
              ## cat("R3: indC = ",indC,"\n")
              g2 <- pdag[indC,indC]
              ## print(g2)
              if (length(g2) <= 1) {
                g2 <- 0
              } else {
                diag(g2) <- rep(1,length(indC)) ## no self reference
              }
              if (any(g2 == 0)) { ## if two nodes in g2 are not connected
                pdag[a,b] <- 1
                pdag[b,a] <- 0
                if (verbose) cat("\nRule 3:",a,"->",b,"\n")
              }
          }
        }
      }
      res@graph <- as(pdag,"graphNEL")
    } else {
      res@graph <- dag2cpdag(pdag2dag(as(pdag,"graphNEL"),keepVstruct = FALSE)$graph)
      status <- 2
      ## res@graph <- res2$graph
    }
  }
  list(pcObj = res, evisit = evisit, xtbl = xtbl, xtbl.orig = xtbl.orig,
       amat0 = amat0, amat1 = amat1, status = status, counter = counter)
}

udag2pdagRelaxed <- function(gInput, verbose = FALSE, unfVect = NULL, solve.confl = FALSE, orientCollider = TRUE, rules = rep(TRUE, 3))
{

##################################################
  ## Internal functions
##################################################

  ## replace 'else if' branch in 'if( !solve.confl )' statement
  orientConflictCollider <- function(pdag, x, y, z) { ## x - y - z
    ## pdag: amat, pdag[x,y] = 1 and pdag[y,x] = 0 means x -> y
    ## x,y,z: colnumber of nodes in pdag
    ## only used if conflicts should be solved

    ## orient x - y
    if (pdag[x,y] == 1) {
      ## x --- y, x --> y => x --> y
      pdag[y,x] <- 0
    } else {
      ## x <-- y, x <-> y => x <-> y
      pdag[x,y] <- pdag[y,x] <- 2
    }

    ## orient z - y
    if (pdag[z,y] == 1) {
      ## z --- y, z --> y => z --> y
      pdag[y,z] <- 0
    } else {
      ## z <-- y, z <-> y => z <-> y
      pdag[z,y] <- pdag[y,z] <- 2
    }

    pdag
  }

  ## TODO: include correct VERBOSE statements
  rule1 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ## Rule 1: a -> b - c goes to a -> b -> c
    ## Interpretation: No new collider is introduced
    ## Out: Updated pdag
    search.pdag <- pdag
    ind <- which(pdag == 1 & t(pdag) == 0, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      ## find all undirected neighbours of b not adjacent to a
      isC <- which(search.pdag[b, ] == 1 & search.pdag[, b] == 1 &
                   search.pdag[a, ] == 0 & search.pdag[, a] == 0)
      if (length(isC) > 0) {
        for (ii in seq_along(isC)) {
          c <- isC[ii]
          ## if the edge between b and c has not been oriented previously,
          ## orient it using normal R1
          if (!solve.confl | (pdag[b,c] == 1 & pdag[c,b] == 1) ) { ## no conflict
            ## !! before, we checked search.pdag, not pdag !!
            if (!is.null(unfVect)) { ## deal with unfaithful triples
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                ## if unfaithful triple, don't orient
                pdag[b, c] <- 1
                pdag[c, b] <- 0
              }
            } else {
              ## don't care about unfaithful triples -> just orient
              pdag[b, c] <- 1
              pdag[c, b] <- 0
              ## cat("Rule 1\n")
            }
            if (verbose)
              cat("\nRule 1':", a, "->", b, " and ",
                  b, "-", c, " where ", a, " and ", c,
                  " not connected and ", a, b, c, " faithful triple: ",
                  b, "->", c, "\n")
          } else if (pdag[b,c] == 0 & pdag[c,b] == 1) {
            ## conflict that must be solved
            ## solve conflict: if the edge is b <- c because of a previous
            ## orientation within for loop then output <->
            if (!is.null(unfVect)) { ## deal with unfaithful triples
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pdag[b, c] <- 2
                pdag[c, b] <- 2
                if (verbose)
                  cat("\nRule 1':", a, "->", b, "<-",
                      c, " but ", b, "->", c, "also possible and",
                      a, b, c, " faithful triple: ", a,"->", b, "<->", c,"\n")
              }
            } else {
              ## don't care about unfaithful triples -> just orient
              pdag[b, c] <- 2
              pdag[c, b] <- 2
              if (verbose)
                cat("\nRule 1':", a, "->", b, "<-",
                    c, " but ", b, "->", c, "also possible and",
                    a, b, c, " faithful triple: ", a,"->", b, "<->", c,"\n")
            } ## unfVect: if else
          } ## conflict: if else
        } ## for c
      } ## if length(isC)
      if (!solve.confl) search.pdag <- pdag
    } ## for ind
    pdag
  }

  rule2 <- function(pdag, solve.confl = FALSE) {
    ## Rule 2: a -> c -> b with a - b: a -> b
    ## Interpretation: Avoid cycle
    ## normal version = conservative version
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      isC <- which(search.pdag[a, ] == 1 & search.pdag[, a] == 0 &
                   search.pdag[, b] == 1 & search.pdag[b, ] == 0)
      for (ii in seq_along(isC)) {
        c <- isC[ii]
        ## if the edge has not been oriented yet, orient it with R2
        ## always do this if you don't care about conflicts
        if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1) ) {
          pdag[a, b] <- 1
          pdag[b, a] <- 0
          if (verbose)
            cat("\nRule 2: Chain ", a, "->", c,
                "->", b, ":", a, "->", b, "\n")
        }
        ## else if the edge has been oriented as a <- b by a previous R2
        else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
          pdag[a, b] <- 2
          pdag[b, a] <- 2
          if (verbose)
            cat("\nRule 2: Chain ", a, "->", c,
                "->", b, ":", a, "<->", b, "\n")
        }
      }
      if (!solve.confl) search.pdag <- pdag
    }
    pdag
  }

  rule3 <- function(pdag, solve.confl = FALSE, unfVect = NULL) {
    ## Rule 3: a-b, a-c1, a-c2, c1->b, c2->b but c1 and c2 not connected;
    ## then a-b => a -> b
    search.pdag <- pdag
    ind <- which(search.pdag == 1 & t(search.pdag) == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      a <- ind[i, 1]
      b <- ind[i, 2]
      c <- which(search.pdag[a, ] == 1 & search.pdag[, a] == 1 &
                 search.pdag[, b] == 1 & search.pdag[b, ] == 0)
      if (length(c) >= 2) {
        cmb.C <- combn(c, 2)
        cC1 <- cmb.C[1, ]
        cC2 <- cmb.C[2, ]
        for (j in seq_along(cC1)) {
          c1 <- cC1[j]
          c2 <- cC2[j]
          if (search.pdag[c1, c2] == 0 && search.pdag[c2,c1] == 0) {
            if (!is.null(unfVect)) {
              if (!any(unfVect == triple2numb(p, c1, a, c2), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c2, a, c1), na.rm = TRUE)) {
                ## if the edge has not been oriented yet, orient it with R3
                if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1) ) {
                  pdag[a, b] <- 1
                  pdag[b, a] <- 0
                  if (!solve.confl) search.pdag <- pdag
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                        a, "->", b, "\n")
                  break
                }
                ## else if: we care about conflicts and  the edge has been oriented as a <- b by a previous R3
                else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
                  pdag[a, b] <- pdag[b, a] <- 2
                  if (verbose)
                    cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                        a, "<->", b, "\n")
                  break
                } ## if solve conflict
              } ## if unf. triple found
            } else { ## if care about unf. triples; else don't care
              if (!solve.confl | (pdag[a, b] == 1 & pdag[b, a] == 1) ) {
                pdag[a, b] <- 1
                pdag[b, a] <- 0
                if (!solve.confl) search.pdag <- pdag
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                      a, "->", b, "\n")
                break
              }
              ## else if: we care about conflicts and  the edge has been oriented as a <- b by a previous R3
              else if (pdag[a, b] == 0 & pdag[b, a] == 1) {
                pdag[a, b] <- pdag[b, a] <- 2
                if (verbose)
                  cat("\nRule 3':", a, c1, c2, "faithful triple: ",
                      a, "<->", b, "\n")
                break
              } ## if solve conflict
            } ## if care about unf. triples
          } ## if c1 and c2 are not adjecent
        } ## for all pairs of c's
      } ## if at least two c's are found
    } ## for all undirected edges
    pdag
  }
##################################################
  ## Main
##################################################

  ## prepare adjacency matrix of skeleton
  if (numEdges(gInput@graph) == 0)
    return(gInput)
  g <- as(gInput@graph, "matrix")
  p <- nrow(g)
  pdag <- g

  ## orient collider
  if (orientCollider) {
    ind <- which(g == 1, arr.ind = TRUE)
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(g[y, ] == 1), x) ## x - y - z
      for (z in allZ) {
        ## check collider condition
        if (g[x, z] == 0 &&
            !((y %in% gInput@sepset[[x]][[z]]) ||
              (y %in% gInput@sepset[[z]][[x]]))) {
          if (length(unfVect) == 0) { ## no unfaithful triples
            if (!solve.confl) { ## don't solve conflicts
              pdag[x, y] <- pdag[z, y] <- 1
              pdag[y, x] <- pdag[y, z] <- 0
            } else { ## solve conflicts
              pdag <- orientConflictCollider(pdag, x, y, z)
            }
          } else { ## unfaithful triples are present
            if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
              if (!solve.confl) { ## don't solve conflicts
                pdag[x, y] <- pdag[z, y] <- 1
                pdag[y, x] <- pdag[y, z] <- 0
              } else { ## solve conflicts
                pdag <- orientConflictCollider(pdag, x, y, z)
              }
            }
          }
        }
      } ## for z
    } ## for i
  } ## end: Orient collider

  ## Rules 1 - 3
  repeat {
    old_pdag <- pdag
    if (rules[1]) {
      pdag <- rule1(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (rules[2]) {
      pdag <- rule2(pdag, solve.confl = solve.confl)
    }
    if (rules[3]) {
      pdag <- rule3(pdag, solve.confl = solve.confl, unfVect = unfVect)
    }
    if (all(pdag == old_pdag))
      break
  } ## repeat

  gInput@graph <- as(pdag, "graphNEL")
  gInput

}

## DEPRECATED! -- use  ida() --
beta.special <- function(dat = NA, x.pos, y.pos, verbose = 0, a = 0.01,
                         myDAG = NA, myplot = FALSE, perfect = FALSE,
                         method = "local", collTest = TRUE, pcObj = NA, all.dags = NA, u2pd = "rand")
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj and all DAGs
  ## can be precomputed
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dat: data
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - verbose: 0=no comments, 1=progress in BB, 2=detail on estimates
  ## - a: significance level of tests for finding CPDAG
  ## - myDAG: needed if bootstrp==FALSE
  ## - myplot: plot estimated graph
  ## - perfect: True cor matrix is calculated from myDAG
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - collTest: True - Exclude orientations of undirected edges that
  ##   introduce a new collider
  ## - pcObj: Fit of PC Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  ## Covariance matrix: Perfect case / standard case
  if (perfect) {
    if(!is(myDAG, "graphNEL")) stop("For perfect-option the true DAG is needed!")
    mcov <- trueCov(myDAG)
    mcor <- cov2cor(mcov)
  } else {
    mcov <- cov(dat)
  }

  ## estimate skeleton and CPDAG of given data
  res <-
    if (is(pcObj, "pcAlgo"))
      pcObj
    else if(perfect)
      pcAlgo.Perfect(mcor, corMethod = "standard",directed = TRUE,u2pd = u2pd)
    else
      pcAlgo(dat, alpha = a, corMethod = "standard",directed = TRUE,u2pd = u2pd)

  ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
  amat <- ad.res <- wgtMatrix(res@graph)
  amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1

  if (method == "local") {
##############################
    ## local method
    ## Main Input: mcov
##############################
    ## find unique parents of x
    wgt.est <- ad.res
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x.pos,] != 0)
    if (y.pos %in% pa1) {
      ## x is parent of y -> zero effect
      beta.hat <- 0
    } else { ## y.pos not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x.pos,] != 0)
      if (verbose == 2) {
        cat("\n\nx=",x.pos,"y=",y.pos,"\n")
        cat("pa1=",pa1,"\n")
        cat("pa2=",pa2,"\n")
      }

      ## estimate beta
      if (length(pa2) == 0) {
        beta.hat <- lm.cov(mcov, y.pos, c(x.pos,pa1))
        if (verbose == 2)
          cat("Fit - y:",y.pos, "x:",c(x.pos,pa1), "|b.hat=", beta.hat)
      } else {
        beta.hat <- NA
        ii <- 1
        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        ## check for new collider
        if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
          beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
          if (verbose == 2)
            cat("\ny:",y.pos,"x:",c(x.pos,pa1),"|b.hat=", beta.hat[ii])
        }## else {
        ##   cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
        ## }
        ## exactly one member of pa2
        for (i2 in seq_along(pa2)) {
          ## check for new collider
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
            ii <-  ii+1
            if (y.pos %in% pa2.t) {
              ## cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
              beta.hat[ii] <- 0
            } else {
              beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
            }
            if (verbose == 2) { cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa2[i2]),
                  "|b.hat=",beta.hat[ii])
}
          } else {
            ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
          }
        }
        ## higher order subsets
        if (length(pa2) > 1) {
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2, i, simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))) {
              pa2.t <- pa.tmp[,j]
              pa2.f <- setdiff(pa2,pa2.t)
              ## teste auf neuen collider
              if (!collTest || !has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
                ii <- ii+1
                if (y.pos %in% pa2.t) {
                  cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                  beta.hat[ii] <- 0
                } else {
                  beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1,pa2.t))
                }
                if (verbose == 2) { cat("\ny:",y.pos,"x:",c(x.pos,pa1,pa2.t),
                      "|b.hat=",beta.hat[ii])
}
              } else {
                ## cat("\nx:",x.pos," pa1:",pa1," pa2.t:",pa2.t," pa2.f:",pa2.f)
              }
            }
          }
        }
      } ## if pa2
    } ## if y in pa1
  } else {
##############################
    ## global method
    ## Main Input: mcov
##############################
    p <- numNodes(res@graph)
    am.pdag <- ad.res
    am.pdag[am.pdag != 0] <- 1
    ## find all DAGs if not provided externally
    ad <- if (is.na(all.dags)) allDags(am.pdag,am.pdag,NULL) else all.dags
    n.dags <- nrow(ad)
    beta.hat <- rep.int(NA,n.dags)
    if (n.dags > 0) {
      if (myplot) {
        ## x11()
        par(mfrow = c(ceiling(sqrt(n.dags)), round(sqrt(n.dags)) ))
      }
      for (i in 1:n.dags) {
        ## compute effect for every DAG
        gDag <- as(matrix(ad[i,],p,p),"graphNEL")
        if (myplot) Rgraphviz::plot(gDag)
        ## path from y to x
        rev.pth <- RBGL::sp.between(gDag,as.character(y.pos),
                                    as.character(x.pos))[[1]]$path
        if (length(rev.pth) > 1) {
          ## if reverse path exists, beta=0
          beta.hat[i] <- 0
        } else {
          ## path from x to y
          pth <- RBGL::sp.between(gDag,as.character(x.pos),
                                  as.character(y.pos))[[1]]$path
          if (length(pth) < 2) {
            ## sic! There is NO path from x to y
            beta.hat[i] <- 0
          } else {
            ## There is a path from x to y
            wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
            pa1 <- which(wgt.unique[x.pos,] != 0)
            if (y.pos %in% pa1) {
              cat("Y in Parents: ",y.pos," in ",pa1,"\n")
              beta.hat[i] <- 0
            } else {
              beta.hat[i] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
            }
            if (verbose == 2)
              cat("Fit - y:",y.pos,"x:",c(x.pos,pa1), "|b.hat=",beta.hat,"\n")
          } ## if length(pth)
        } ## if rev.pth
      } ## for n.dags
    } ## if n.dags
  } ## if method
  beta.hat
} ## {beta.special}



## DEPRECATED! -- use  ida() / idafast() --
beta.special.pcObj <- function(x.pos,y.pos,pcObj,mcov = NA,amat = NA,amatSkel = NA,
                               t.amat = NA)
{
  ## Purpose: Estimate the causal effect of x on y; the pcObj has to be
  ## precomputed. This method is intended to be a fast version of
  ##
  ## beta.special(dat=NA,x.pos,y.pos,verbose=0,a=NA,myDAG=NA,myplot=FALSE,
  ## perfect=FALSE,method="local",collTest=TRUE,pcObj=pcObj,all.dags=NA,u2pd="relaxed")
  ##
  ## Thus, this is a faster version for the local method given a
  ## precomputed PC-Algo Object (relaxed udag2pdag, so CPDAG might not
  ## be a real CPDAG; this does not matter, since we try not to extend).
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - pcObj: Fit of pc Algorithm (semidirected); if this is available, no
  ##   new fit is done
  ## - mcov: covariance matrix of pcObj fit
  ## - amat,amatSkel,g2,t.amat are variants of the adjacency matrix that
  ##   are used internally but can be precomputed; the relevant code
  ##   is commented out
  ## ----------------------------------------------------------------------
  ## Value: List with two elements
  ## - beta.res: beta.causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 21 Nov 2007, 11:18

  cat("This function is deprecated and is only kept for backward compatibility.
Please use ida or idaFast instead\n")

  if (is.na(amat) | is.na(amatSkel) | is.na(t.amat)) {
    ## Code for computing precomputable variables
    ## prepare adjMatrix and skeleton {MM FIXME : can be improved}
    amat <- wgtMatrix(pcObj@graph)
    amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
    t.amat <- t(amat)
    amatSkel <- amat + t.amat
    amatSkel[amatSkel != 0] <- 1
  }

  ## find unique parents of x
  tmp <- amat-t.amat
  tmp[which(tmp < 0)] <- 0
  wgt.unique <- tmp
  pa1 <- which(wgt.unique[x.pos,] != 0)
  if (y.pos %in% pa1) {
    cat("Y in Parents: ",y.pos," in ",pa1,"\n")
    beta.hat <- 0
  } else { ## y.pos not in pa1
    ## find ambiguous parents of x
    wgt.ambig <- amat-wgt.unique
    pa2 <- which(wgt.ambig[x.pos,] != 0)
    pa2 <- setdiff(pa2,y.pos)
    ## estimate beta
    if (length(pa2) == 0) {
      beta.hat <- lm.cov(mcov,y.pos,c(x.pos,pa1))
    } else {
      beta.hat <- NA
      ii <- 1
      ## no member of pa2
      ## check for new collider
      pa2.f <- pa2
      pa2.t <- NA
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        beta.hat[ii] <- lm.cov(mcov,y.pos,c(x.pos,pa1))
      }
      ## exactly one member of pa2
      for (i2 in seq_along(pa2)) {
        ## check for new collider
        pa2.f <- pa2[-i2]
        pa2.t <- pa2[i2]
        if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
          ii <-  ii+1
          beta.hat[ii] <-
            if (y.pos %in% pa2.t) {
              cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
              0
            } else lm.cov(mcov,y.pos,c(x.pos,pa1,pa2[i2]))
        }
      }
      ## higher order subsets
      if (length(pa2) > 1) {
        for (i in 2:length(pa2)) {
          pa.tmp <- combn(pa2, i, simplify = TRUE)
          for (j in seq_len(ncol(pa.tmp))) {
            ## teste auf neuen collider
            pa2.t <- pa.tmp[,j]
            pa2.f <- setdiff(pa2,pa2.t)
            if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
              ii <- ii+1
              beta.hat[ii] <-
                if (y.pos %in% pa2.t) {
                  cat("Y in Parents: ",y.pos," in ",pa2.t,"\n")
                  0
                } else lm.cov(mcov,y.pos,c(x.pos,pa1,pa2.t))
            }
          }
        }
      } ## if pa2
    } ## length(pa2)
  } ## y.pos %in% pa2
  beta.hat
} ## {beta.special.pcObj}

##' @title
##' @param C covariance matrix
##' @param y column of response
##' @param x columns of expl. vars
##' @return
lm.cov <- function (C, y, x) {
  solve(C[x, x], C[x, y, drop = FALSE])[1, ]
}

causalEffect <- function(g,y,x) {
  ## Compute true causal effect of x on y in g
  wmat <- wgtMatrix(g)
  p <- ncol(wmat)
  vec <- matrix(0,p,1)
  vec[x] <- 1
  ## compute and return  beta_{true} :
  if(y-x > 1) {
    for (i in (x+1):y) vec[i] <- wmat[i,]%*%vec
    vec[y]
  } else {
    wmat[y,x]
  }
}

has.new.coll <- function(amat,amatSkel, x, pa1, pa2.t, pa2.f) {
  ## Check if undirected edges that are pointed to x create a new v-structure
  ## Additionally check, if edges that are pointed away from x create
  ## new v-structure; i.e. x -> pa <- papa would be problematic
  ## pa1 are definit parents of x
  ## pa2 are undirected "parents" of x
  ## pa2.t are the nodes in pa2 that are directed towards pa2
  ## pa2.f are the nodes in pa2 that are directed away from pa2
  ## Value is TRUE, if new collider is introduced
  res <- FALSE
  if (length(pa2.t) > 0 && !all(is.na(pa2.t))) {
    ## check whether all pa1 and all pa2.t are connected;
    ## if no, there is a new collider
    if (length(pa1) > 0 && !all(is.na(pa1))) {
      res <- min(amatSkel[pa1, pa2.t]) == 0 ## TRUE if new collider
    }
    ## in addition, all pa2.t have to be connected
    if (!res && length(pa2.t) > 1) {
      A2 <- amatSkel[pa2.t,pa2.t]
      diag(A2) <- 1
      res <- min(A2) == 0 ## TRUE if new collider
    }
  }
  if (!res && length(pa2.f) > 0 && !all(is.na(pa2.f))) {
    ## consider here only the DIRECTED Parents of pa2.f
    ## remove undirected edges
    A <- amat-t(amat)
    A[A < 0] <- 0
    ## find parents of pa2.f
    cA <- colSums(A[pa2.f,,drop = FALSE])
    papa <- setdiff(which(cA != 0), x)
    ## if any node in papa is not directly connected to x, there is a new
    ## collider
    if (length(papa) > 0)
      res <- min(amatSkel[x,papa]) == 0 ## TRUE if new collider
  }
  res
}

allDags <- function(gm,a,tmp, verbose = FALSE)
{
  ## Purpose: Find all DAGs for a given PDAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - gm: Adjacency matrix of initial PDAG; only 0-1 entries
  ##   i -> j iff gm(j,i)=1
  ## - a: copy of gm
  ## - tmp: "current set of DAGs", initially NULL
  ## ----------------------------------------------------------------------
  ## Value:
  ## - one 0/1 adj.matrix per row
  ## Reversion to graph: as(matrix(res[i,],p,p),"graphNEL")
  ## Reversion to wgtMatrix (i->j iff a[j,i]=1): t(matrix(res[i,],p,p))
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  7 Apr 2008, 14:08
  if (sum(a) == 0) {
    if (verbose) {
      cat("Last Call - Final Graph: \n")
      print(gm)
      cat("#################### \n")
    }
    tmp2 <- rbind(tmp,c(t(gm)))
    if (all(!duplicated(tmp2))) tmp <- tmp2
  } else {
    sinks <- find.sink(a)
    if (verbose) {
      cat("Main Call: ################## \n")
        print(gm)
      print(a)
      cat("Sinks: ",sinks,"\n")
    }
    for(x in sinks) {
      if (verbose) cat("Try removing", x," in a.\n")
      gm2 <- gm
      a2 <- a
      if (adj.check(a,x)) {
        inc.to.x <- a[, x] == 1 & a[x, ] == 1
        if (any(inc.to.x)) {
          real.inc.to.x <- as.numeric(rownames(a)[inc.to.x])
          real.x <- as.numeric(rownames(a)[x])
          gm2[real.x, real.inc.to.x] <- 1
          gm2[real.inc.to.x, real.x] <- 0
        }
        a2 <- a[-x,-x]
        if (verbose) {
          cat("Removed sink",as.numeric(rownames(a)[x]),
              "in g (", x,"in a).\n")
          cat("New graphs: \n")
          print(gm2)
          print(a)
        }
        tmp <- allDags(gm2, a2, tmp, verbose)
        ##     ------- *recursively*
      }
    }
  }
  tmp
}

pcAlgo.Perfect <- function(C, cutoff = 1e-8, corMethod = "standard", verbose = 0,
                           directed = FALSE, u2pd = "rand", psepset = FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - C: True Correlation matrix
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - verbose: 0-no output, 1-small output, 2-details
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modification: Diego Colombo, Sept 2009
  ## backward compatibility
  stopifnot((p <- nrow(C)) >= 2)
  cl <- match.call()
  seq_p <- seq_len(p)# 1:p
  pcMin <- matrix(Inf, p,p)
  diag(pcMin) <- 0
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  n.edgetests <- numeric(1)# final length = max { ord}
  ## G := complete graph :
  G <- matrix(TRUE, p,p)
  diag(G) <- FALSE

  done <- FALSE
  ord <- 0L
  while (!done && any(G)) {
    n.edgetests[ord+1] <- 0
    done <- TRUE
    ind <- which(G, arr.ind = TRUE)
    ## For comparison with C++ sort according to first row
    ind <- ind[order(ind[,1]), ]
    remEdges <- nrow(ind)
    if(verbose >= 1)
      cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
    for (i in 1:remEdges) {
      if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=",i,"|iMax=",remEdges,"\n")
      x <- ind[i,1]
      y <- ind[i,2]
      ##      done <- !G[y,x] # i.e. (x,y) was not already deleted in its (y,x) "version"
      ##      if(!done) {
      if (G[y,x]) {
        nbrsBool <- G[,x]
        nbrsBool[y] <- FALSE
        nbrs <- seq_p[nbrsBool]
        length_nbrs <- length(nbrs)
        ## if(verbose)

        ##        done <- length_nbrs < ord
        ##        if (!done) {
        if (length_nbrs >= ord) {
          if (length_nbrs > ord) done <- FALSE
          S <- seq_len(ord)

          ## now includes special cases (ord == 0) or (length_nbrs == 1):
          repeat { ## condition w.r.to all  nbrs[S] of size 'ord' :
            ##  if (condIndFisherZ(x,y, nbrs[S], C,n, cutoff,verbose)) {
            ## MM: want to use this: --- but it changes the result in some cases!!
            ##            cat("X=",x,"|Y=",y,"|ord=",ord,"|nbrs=",nbrs[S],"|iMax=",remEdges,"\n")
            n.edgetests[ord+1] <- n.edgetests[ord+1]+1
            pc.val <- pcorOrder(x,y,nbrs[S],C)
            if (abs(pc.val) < pcMin[x,y]) pcMin[x,y] <- abs(pc.val)
            if (verbose >= 2) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("pc:",pc.val,"\n"))
            if (abs(pc.val) <= cutoff) {
              ##              ##  pnorm(abs(z), lower.tail = FALSE) is the p-value
              G[x,y] <- G[y,x] <- FALSE
              sepset[[x]][[y]] <- nbrs[S]
              break
            }
            else {
              nextSet <- getNextSet(length_nbrs, ord, S)
              if(nextSet$wasLast)
                break
              S <- nextSet$nextSet
            }
          }

        } } ## end if(!done)

    } ## end for(i ..)
    ord <- ord+1L
    ##    n.edgetests[ord] <- remEdges
  }

  if (psepset) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,] == 1L),x) ## x-y-z

      if (length(allZ) > 0) {
        for (j in seq_along(allZ)) {
          z <- allZ[j]
          if ((amat[x,z] == 0L) && !((y %in% sepset[[x]][[z]]) | (y %in% sepset[[z]][[x]]))) {
            if (verbose == 2) {
              cat("\n",x,"*->",y,"<-*",z,"\n")
              cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
            }

            ## x o-> y <-o z
            amat[x,y] <- amat[z,y] <- 2L

          } ## if
        } ## for
      } ## if
    } ## for

    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,] != 0L)) {
        tf1 <- setdiff(reach(x,-1,-1,amat),x)
        for (y in seq_p[amat[x,] != 0L]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf) > 0) {
            pc.val <- pcorOrder(x,y,tf,C)
            if (abs(pc.val) < pcMin[x,y]) {
              pcMin[x,y] <- abs(pc.val)
            }
            if (abs(pc.val) <= cutoff) {
              ## delete x-y
              amat[x,y] <- amat[y,x] <- 0L
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
              if (verbose == 2) {
                cat("Delete edge",x,"-",y,"\n")
              }
            }
            if (verbose == 2) {
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - pc = ",pc.val,"\n")
            }
          }
        }
      }
    }

    G[amat == 0L] <- FALSE
    G[amat == 1L] <- TRUE
  } ## end if(psepset)

  if(verbose >= 1) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }

  ## transform matrix to graph object :
  if (sum(G) == 0) {
    Gobject <- new("graphNEL", nodes = as.character(seq_p))
  }
  else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    Gobject <- as(G,"graphNEL")
  }

  for (i in 1:(p-1)) {
    for (j in 2:p) {
      pcMin[i,j] <- pcMin[j,i] <- min(pcMin[i,j],pcMin[j,i])
    }
  }

  res <- new("pcAlgo",
             graph = Gobject,
             call = cl, n = as.integer(1), max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = pcMin)

  if (directed)
    switch(u2pd,
           "rand" = udag2pdag(res),
           "retry" = udag2pdagSpecial(res)$pcObj,
           "relaxed" = udag2pdagRelaxed(res))
  else
    res
}## pcAlgo.Perfect

### reach(): currently only called from  pcAlgo() and pcAlgo.Perfect()
### -------  and only in "possibledsep" version
## Function that computes the Possible d-sepset, done by Spirtes
reach <- function(a,b,c,adjacency)
{
  ## reachable      set of vertices;
  ## edgeslist      array[1..maxvertex] of list of edges
  ## numvertex      integer
  ## labeled        array (by depth) of list of edges that have been labeled

  makeedge <- function(x,y) list(list(x,y))

  legal.pdsep <- function(r,s) {
    ## Modifying global 'edgeslist'
    if ((adjacency[r[[1]],r[[2]]] == 2 &&
         adjacency[s,     r[[2]]] == 2 && r[[1]] != s) ||
        (adjacency[r[[1]],s] != 0 && r[[1]] != s)) {
      edgeslist[[r[[2]]]] <<- setdiff(edgeslist[[r[[2]]]],s)
      makeedge(r[[2]],s)
    }
  }

  initialize.pdsep <- function(x,y) mapply(makeedge, x = x, y = y)

  labeled <- list()
  numvertex <- dim(adjacency)[1]
  edgeslist <- list()
  for (i in 1:numvertex)
    edgeslist <- c(edgeslist,list(which(adjacency[,i] != 0)))
  labeled[[1]] <- initialize.pdsep(a, edgeslist[[a]])
  edgeslist[[a]] <- list()
  depth <- 2
  repeat {
    labeled[[depth]] <- list()
    for (i in seq_along(labeled[[depth-1]])) {
      lab.i <- labeled[[depth-1]][[i]]
      edgestemp <- edgeslist[[lab.i[[2]]]]
      if (length(edgestemp) == 0) break
      for (j in seq_along(edgestemp))
        labeled[[depth]] <- union(legal.pdsep(lab.i, edgestemp[[j]]),
                                  labeled[[depth]])
    }
    if (length(labeled[[depth]]) == 0)
      break
    ## else :
    depth <- depth  + 1
  }
  unique(unlist(labeled))
}


plotAG <- function(amat)
{
  ## Purpose: Plot ancestral graph
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: Adjacency matrix
  ##   amat[i,j]=3 & amat[j,i]=1 iff i 1-3 j
  ##   "0": no edge; "1": circle; "2": arrow; "3": tail
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 16 Feb 2009, 18:01
  check.Rgraphviz()

  g <- as(amat,"graphNEL")
  nn <- nodes(g)
  p <- numNodes(g)
  n.edges <- numEdges(g)
  ah.list <- at.list <- rep("none",n.edges)
  counter <- 0
  list.names <- NULL
  amat[amat == 1] <- "odot"
  amat[amat == 2] <- "normal"
  amat[amat == 3] <- "none"
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      x <- nn[i]
      y <- nn[j]
      if (amat[x,y] != 0) {
        counter <- counter + 1
        ah.list[[counter]] <- amat[x,y]
        at.list[[counter]] <- amat[y,x]
        list.names <- c(list.names,paste(x,"~",y,sep = ""))
      }
    }
  }
  names(ah.list) <- names(at.list) <- list.names

  edgeRenderInfo(g) <- list(arrowhead = ah.list, arrowtail = at.list)
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}

skeleton <- function(suffStat, indepTest, alpha, labels, p,
		     method = c("stable", "original", "stable.fast"), m.max = Inf,
		     fixedGaps = NULL, fixedEdges = NULL,
		     NAdelete = TRUE, verbose = FALSE)
{
  ## Purpose: Perform undirected part of PC-Algorithm, i.e.,
  ## estimate skeleton of DAG given data
  ## Order-independent version! NEU
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat: List containing all necessary elements for the conditional
  ##             independence decisions in the function "indepTest".
  ## - indepTest: predefined function for testing conditional independence
  ## - alpha: Significance level of individual partial correlation tests
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G, sepset, pMax, ord, n.edgetests
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 09.12.2009
  ## Modification: Diego Colombo; Martin Maechler

  ## x,y,S konstruieren
  ##-   tst <- try(indepTest(x,y,S, obj))
  ##-   if(inherits(tst, "try-error"))
  ##-     stop("the 'indepTest' function does not work correctly with 'obj'")

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p))
      p <- length(labels)
    else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    ## Don't want message, in case this is called e.g. from fciPlus():
    ## else
    ##   message("No need to specify 'p', when 'labels' is given")
  }
  seq_p <- seq_len(p)
  method <- match.arg(method)
  ## C++ version still has problems under Windows; will have to check why
                                        ##  if (method == "stable.fast" && .Platform$OS.type == "windows") {
                                        ##    method <- "stable"
                                        ##    warning("Method 'stable.fast' is not available under Windows; using 'stable' instead.")
                                        ##  }

  ## G := !fixedGaps, i.e. G[i,j] is true  iff  i--j  will be investigated
  if (is.null(fixedGaps)) {
    G <- matrix(TRUE, nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedGaps), c(p, p)))
    stop("Dimensions of the dataset and fixedGaps do not agree.")
  else if (!identical(fixedGaps, t(fixedGaps)) )
    stop("fixedGaps must be symmetric")
  else
    G <- !fixedGaps

  diag(G) <- FALSE

  if (any(is.null(fixedEdges))) { ## MM: could be sparse
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  }
  else if (!identical(dim(fixedEdges), c(p, p)))
    stop("Dimensions of the dataset and fixedEdges do not agree.")
  else if (!identical(fixedEdges, t(fixedEdges)) )
    stop("fixedEdges must be symmetric")

  if (method == "stable.fast") {
    ## Do calculation in C++...
    if (identical(indepTest, gaussCItest))
      indepTestName <- "gauss"
    else
      indepTestName <- "rfun"
    options <- list(
      verbose = as.integer(verbose),
      m.max = as.integer(ifelse(is.infinite(m.max), p, m.max)),
      NAdelete = NAdelete)
    res <- .Call("estimateSkeleton", G, suffStat, indepTestName, indepTest, alpha, fixedEdges, options);
    G <- res$amat
    ## sepset <- res$sepset
    sepset <- lapply(seq_p, function(i) c(
      lapply(res$sepset[[i]], function(v) if(identical(v, as.integer(-1))) NULL else v),
      vector("list", p - length(res$sepset[[i]])))) # TODO change convention: make sepset triangular
    pMax <- res$pMax
    n.edgetests <- res$n.edgetests
    ord <- length(n.edgetests) - 1L
  }
  else {
    ## Original R version

    pval <- NULL
    sepset <- lapply(seq_p, function(.) vector("list",p))# a list of lists [p x p]
    ## save maximal p value
    pMax <- matrix(-Inf, nrow = p, ncol = p)
    diag(pMax) <- 1
    done <- FALSE
    ord <- 0L
    n.edgetests <- numeric(1)# final length = max { ord}
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord1 <- ord+1L] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[, 1]), ]
      remEdges <- nrow(ind)
      if (verbose)
        cat("Order=", ord, "; remaining edges:", remEdges,"\n",sep = "")
      if(method == "stable") {
        ## Order-independent version: Compute the adjacency sets for any vertex
        ## Then don't update when edges are deleted
        G.l <- split(G, gl(p,p))
      }
      for (i in 1:remEdges) {
        if(verbose && (verbose >= 2 || i%%100 == 0)) cat("|i=", i, "|iMax=", remEdges, "\n")
        x <- ind[i, 1]
        y <- ind[i, 2]
        if (G[y, x] && !fixedEdges[y, x]) {
          nbrsBool <- if(method == "stable") G.l[[x]] else G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord)
              done <- FALSE
            S <- seq_len(ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord1] <- n.edgetests[ord1] + 1
              pval <- indepTest(x, y, nbrs[S], suffStat)
              if (verbose)
                cat("x=", x, " y=", y, " S=", nbrs[S], ": pval =", pval, "\n")
              if(is.na(pval))
                pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
              if (pMax[x, y] < pval)
                pMax[x, y] <- pval
              if(pval >= alpha) { # independent
                G[x, y] <- G[y, x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              }
              else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if (nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            } ## repeat
          }
        }
      }# for( i )
      ord <- ord + 1L
    } ## while()
    for (i in 1:(p - 1)) {
      for (j in 2:p)
        pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j,i])
    }
  }

  ## transform matrix to graph object :
  Gobject <-
    if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
    } else {
      colnames(G) <- rownames(G) <- labels
      as(G,"graphNEL")
    }

  ## final object
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = as.integer(ord - 1), n.edgetests = n.edgetests,
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}## end{ skeleton }




pc <- function(suffStat, indepTest, alpha, labels, p,
               fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE, m.max = Inf,
               u2pd = c("relaxed", "rand", "retry"),
               skel.method = c("stable", "original", "stable.fast"),
               conservative = FALSE, maj.rule = FALSE,
               solve.confl = FALSE, verbose = FALSE)
{
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdag
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - gTrue: Graph suffStatect of true DAG
  ## - conservative: If TRUE, conservative PC is done
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo, Markus Kalisch

  ## Initial Checks
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  u2pd <- match.arg(u2pd)
  skel.method <- match.arg(skel.method)
  if(u2pd != "relaxed") {
    if (conservative || maj.rule)
      stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")

    if (solve.confl)
      stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n can only be run with 'u2pd = relaxed'")
  }

  if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")

  ## Skeleton
  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  skel@call <- cl # so that makes it into result

  ## Orient edges
  if (!conservative && !maj.rule) {
    switch (u2pd,
            "rand" = udag2pdag(skel),
            "retry" = udag2pdagSpecial(skel)$pcObj,
            "relaxed" = udag2pdagRelaxed(skel, verbose = verbose, solve.confl = solve.confl))
  }
  else { ## u2pd "relaxed" : conservative _or_ maj.rule

    ## version.unf defined per default
    ## Tetrad CPC works with version.unf=c(2,1)
    ## see comment on pc.cons.intern for description of version.unf
    pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                          version.unf = c(2,1), maj.rule = maj.rule, verbose = verbose)
    udag2pdagRelaxed(pc.$sk, verbose = verbose,
                     unfVect = pc.$unfTripl, solve.confl = solve.confl)
  }
} ## {pc}


gSquareBin <- function(x, y, S, dm, adaptDF = FALSE, n.min = 10*df,
                       verbose = FALSE)
{
  ## Purpose: G^2 statistic to test for (conditional) independence
  ##          of *binary* variables   X and Y given S  --> ../man/binCItest.Rd
  ## -------------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y conditionally independent given S (S can be NULL)?
  ## - dm: data matrix (rows: samples, columns: variables) with binary entries
  ## - verbose: if TRUE, some additional info is outputted during the
  ##            computations
  ## - adaptDF: lower the degrees of freedom by one for each zero count.
  ##            The value for the DF cannot go below 1.
  ## -------------------------------------------------------------------------

  stopifnot((n <- nrow(dm)) >= 1) # nr of samples
  if(!all(as.logical(dm) == dm))
    stop("'dm' must be binary, i.e. with values in {0,1}")

  if(verbose) cat('Edge ',x,'--',y, ' with subset S =', S,'\n')

  lenS <- length(S)
  ## degrees of freedom assuming no structural zeros
  df <- 2^lenS

  if (n < n.min) { ## not enough samples to perform the test:
    warning(gettextf(
      "n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)",
                     n, n.min), domain = NA)
    return( 1 )
  }
  ## else --  enough data to perform the test
  d.x1 <- dm[,x] + 1L
  d.y1 <- dm[,y] + 1L
  if(lenS <= 5) { # bei gSquareDis lenS <= 4
    n12 <- 1:2
    switch(lenS+1L, {
      ## lenS == 0 ----------------------------------
      nijk <- array(0L, c(2,2)) # really 'nij', but 'nijk' is "global" name
      for (i in n12) {
        d.x.i <- d.x1 == i
	for (j in n12)
	  nijk[i,j] <- sum(d.x.i & d.y1 == j)
      }
      ## marginal counts
      ## compute G^2
      t.log <- n*(nijk / tcrossprod(rowSums(nijk), colSums(nijk)))
    } ,

    { ## lenS == 1 ----------------------------------

      ## a.pos <- sort(c(x,y,S))
      dmS.1 <- dm[,S] + 1L
      nijk <- array(0L, c(2,2,2))
      for(i in n12) {
        d.x.i <- d.x1 == i
        for(j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for(k in n12)
            nijk[i,j,k] <- sum(d.x.i.y.j & dmS.1 == k)
        }
      }
      alt <- c(x,y,S)
      c <- which(alt == S)
      nik <- apply(nijk,c,rowSums)
      njk <- apply(nijk,c,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,2))
      if(c == 3) {
        for (k in n12)
          t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k], njk[,k]))
      } else if(c == 1) {
        for (k in n12)
          t.log[k,,] <- nijk[k,,]*(nk[k] / tcrossprod(nik[,k], njk[,k]))
      } else { ## c == 2 (?)
        for (k in n12)
          t.log[,k,] <- nijk[,k,]*(nk[k] / tcrossprod(nik[,k], njk[,k]))
      }
    } ,

    { ## lenS == 2 ----------------------------------

      ## a.pos <- sort(c(x,y,S))
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      nijk2 <- array(0L, c(2,2,2,2))
      for(i in n12) {
        d.x.i <- d.x1 == i
	for(j in n12) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for(k in n12) {
            d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
            for(l in n12)
              nijk2[i,j,k,l] <- sum(d.x.y.S1 & dmS2.1 == l)
          }
        }
      }

      nijk <- array(0L, c(2,2,4))
      for(i in n12) {
	for(j in n12)
          nijk[,,2*(i-1)+j] <- nijk2[,,i,j]
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,4))
      for (k in 1:4)
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k], njk[,k]))
    } ,

    { ## lenS == 3 ----------------------------------
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      dmS3.1 <- dm[,S[3]] + 1L
      nijk <- array(0L, c(2,2,8))
      for(i1 in n12) {
        d.x.i <- d.x1 == i1
	for(i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for(i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for(i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for(i5 in n12)
                nijk[i1,i2,4*(i3-1)+2*(i4-1)+i5] <-
			sum(d.xy.S1S2 & dmS3.1 == i5)
            }
          }
        }
      }

      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,8))
      for (k in 1:8)
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k], njk[,k]))
    } ,

    { ## lenS == 4 ----------------------------------
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      dmS3.1 <- dm[,S[3]] + 1L
      dmS4.1 <- dm[,S[4]] + 1L
      nijk <- array(0L, c(2,2,16))
      for(i1 in n12) {
        d.x.i <- d.x1 == i1
	for(i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for(i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for(i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for(i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == i5
                for(i6 in n12)
                  nijk[i1,i2,8*(i3-1)+4*(i4-1)+2*(i5-1)+i6] <-
			sum(d.xy.S1S2S3 & dmS4.1 == i6)
              }
            }
          }
        }
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0, c(2,2,16))
      for (k in 1:16)
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k], njk[,k]))

    } ,

    { ## lenS == 5 ----------------------------------
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      dmS3.1 <- dm[,S[3]] + 1L
      dmS4.1 <- dm[,S[4]] + 1L
      dmS5.1 <- dm[,S[5]] + 1L
      nijk <- array(0L, c(2,2,32))
      for(i1 in n12) {
        d.x.i <- d.x1 == i1
	for(i2 in n12) {
          d.x.y.i12 <- d.x.i & d.y1 == i2
          for(i3 in n12) {
            d.x.y.S1 <- d.x.y.i12 & dmS1.1 == i3
            for(i4 in n12) {
              d.xy.S1S2 <- d.x.y.S1 & dmS2.1 == i4
              for(i5 in n12) {
                d.xy.S1S2S3 <- d.xy.S1S2 & dmS3.1 == i5
                for(i6 in n12) {
                  d.xy.S1S2S3S4 <- d.xy.S1S2S3 & dmS4.1 == i6
                  for(i7 in n12)
                    nijk[i1,i2,16*(i3-1)+8*(i4-1)+4*(i5-1)+2*(i6-1)+i7] <-
 			sum(d.xy.S1S2S3S4 & dmS5.1 == i7)
                }
              }
            }
          }
        }
      }

      nik <- apply(nijk,3,rowSums)
      njk <- apply(nijk,3,colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0,c(2,2,32))
      for (k in 1:32)
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k], njk[,k]))

    })# end {lenS = 5}, end{switch}
  }
  else { # --- |S| >= 6 -------------------------------------------------
    nijk <- array(0L, c(2,2,1))
    ## first sample 'by hand' to avoid if/else in the for-loop
    i <- d.x1[1]
    j <- d.y1[1]
    ## create directly a list of all k's  -- MM{FIXME}: there must be a better way
    k <- NULL
    lapply(as.list(S), function(x) { k <<- cbind(k,dm[,x]+1); NULL })
    ## first set of subset values
    parents.count <- 1L ## counter variable corresponding to the number
    ## of value combinations for the subset varibales
    ## observed in the data
    parents.val <- t(k[1,])
    nijk[i,j,parents.count] <- 1L        # cell counts

    ## Do the same for all other samples. If there is already a table
    ## for the subset values of the sample, increase the corresponding
    ## cell count. If not, create a new table and set the corresponding
    ## cell count to 1.
    for (it.sample in 2:n) {
      new.p <- TRUE
      i <- d.x1[it.sample]
      j <- d.y1[it.sample]
      ## comparing the current values of the subset variables to all
      ## already existing combinations of subset variables values
      t.comp <- t(parents.val[1:parents.count,]) == k[it.sample,]
      ## Have to be careful here. When giving dimension to a list,
      ## R fills column after column, and NOT row after row.
      dim(t.comp) <- c(lenS,parents.count)
      for (it.parents in 1:parents.count) {
        ## check if the present combination of value alreay exists
        if(all(t.comp[,it.parents])) {
          ## if yes, increase the corresponding cell count
          nijk[i,j,it.parents] <- nijk[i,j,it.parents] + 1L
          new.p <- FALSE
          break
        }
      }
      ## if the combination of subset values is new...
      if (new.p) {
        ## ...increase the number of subset 'types'
        parents.count <- parents.count + 1L
        if (verbose >= 2)
          cat(sprintf(' adding new parents (count = %d) at sample %d\n',
                      parents.count, it.sample))
        ## ...add the new subset to the others
        parents.val <- rbind(parents.val, k[it.sample,])
        ## ...update the cell counts (add new array)
        nijk <- abind(nijk, array(0,c(2,2,1)))
        nijk[i,j,parents.count] <- 1L
      }## if(new.p)
    }## end for(it.sample ..)

    nik <- apply(nijk,3,rowSums)
    njk <- apply(nijk,3,colSums)
    nk <- colSums(njk)
    ## compute G^2
    t.log <- array(0, c(2,2,parents.count))
    for (k in 1:parents.count)
      t.log[,,k] <- nijk[,,k] * (nk[k] / tcrossprod(nik[,k],njk[,k]))

  } ## |S| >= 6

  G2 <- sum(2 * nijk * log(t.log), na.rm = TRUE)

  if (adaptDF && lenS > 0) {
    ## lower the degrees of freedom according to the amount of zero
    ## counts; add zero counts corresponding to the number of parents
    ## combinations that are missing
    zero.counts <- sum(nijk == 0L) + 4*(2^lenS-dim(nijk)[3])
    ndf <- max(1, df-zero.counts)
    if(verbose) cat("adaptDF: (df=",df,", zero.counts=",zero.counts,
                    ") ==> new df = ", ndf, "\n", sep = "")
    df <- ndf
  }

  pchisq(G2, df, lower.tail = FALSE)# i.e. == 1 - P(..)

}## gSquareBin()

gSquareDis <- function(x, y, S, dm, nlev, adaptDF = FALSE, n.min = 10*df,
                       verbose = FALSE) {

  ## Purpose: G^2 statistic to test for (conditional) independence of
  ##          __discrete__ X and Y given S  --> ../man/disCItest.Rd
  ## -------------------------------------------------------------------
  ## Arguments:
  ## - x,y,S: Are x,y conditionally independent given S (S can be NULL)?
  ## - dm: data matrix (rows: samples, columns: variables) with
  ##       discrete entries
  ## - nlev: vector with numbers of levels for each variable
  ## - verbose: if TRUE, some additional info is outputted during the
  ##            computations
  ## - adaptDF: lower the degrees of freedom by one for each zero count.
  ##            The value for the DF cannot go below 1.
  ## -------------------------------------------------------------------

  stopifnot((n <- nrow(dm)) >= 1, # nr of samples
            (p <- ncol(dm)) >= 2) # nr of variables or nodes
  if(!all(1 <= c(x,y,S) & c(x,y,S) <= p))
    stop("x, y, and S must all be in {1,..,p}, p=",p)
  if(any(as.integer(dm) != dm))
    stop("'dm' must be discrete, with values in {0,1,..}")
  if(!any(dm == 0))
    stop("'dm' must have values in {0,1,..} with at least one '0' value")

  if(verbose) cat('Edge ', x,'--',y, ' with subset S =', S,'\n')

  lenS <- length(S)
  if(missing(nlev) || is.null(nlev))
    nlev <- vapply(seq_len(p),
		   function(j) length(levels(factor(dm[,j]))), 1L)
  else
    stopifnot(is.numeric(nlev), length(nlev) == p, !is.na(nlev))
  if(!all(nlev >= 2))
    stop("Each variable, i.e., column of 'dm', must have at least two different values")
  nl.x <- nlev[x]
  nl.y <- nlev[y]
  nl.S <- nlev[S]
  ## degrees of freedom assuming no structural zeros
  df <- (nl.x-1)*(nl.y-1)*prod(nl.S)
  if (n < n.min) { ## not enough samples to perform the test:
    warning(gettextf(
      "n=%d is too small (n < n.min = %d ) for G^2 test (=> treated as independence)",
                     n, n.min), domain = NA)
    return( 1 )   ## gerster prob=0
  }
  ## else --  enough data to perform the test
  i.ny <- seq_len(nl.y) # = 1:nl.y
  lenS <- length(S)
  d.x1 <- dm[,x] + 1L
  d.y1 <- dm[,y] + 1L
  if(lenS <= 4) { # bei gSquareBin lenS <= 5
    switch(lenS+1L, { ## lenS == 0 -------------------
      nij <- array(0L, c(nl.x,nl.y))
      for (i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for (j in i.ny)
          nij[i,j] <- sum(d.x.i & d.y1 == j)
      }
      ## marginal counts
      t.X <- rowSums(nij)
      t.Y <- colSums(nij)
      ## compute G^2
      t.log <- n*(nij/ tcrossprod(t.X, t.Y))
      t.G2 <- 2*nij * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    } ,
    { ## lenS == 1 ----------------------------------
      in.S <- seq_len(nl.S)
      dmS.1 <- dm[,S] + 1L
      nijk <- array(0L, c(nl.x,nl.y,nl.S))
      for(i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for(j in i.ny) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for(k in in.S)
            nijk[i,j,k] <- sum(d.x.i.y.j & dmS.1 == k)
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0, c(nl.x,nl.y,prod(nl.S)))
      for (k in 1:prod(nl.S))
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k],njk[,k]))

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    } ,
    { ## lenS == 2 ----------------------------------
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      nijk <- array(0L, c(nl.x,nl.y,nl.S1*nl.S2))
      for(i in 1:nl.x) {
        d.x.i <- d.x1 == i
        for(j in i.ny) {
          d.x.i.y.j <- d.x.i & d.y1 == j
          for(k in in.S1) {
            d.x.y.S1 <- d.x.i.y.j & dmS1.1 == k
            for(l in in.S2)
              nijk[i,j,nl.S2*(k-1)+l] <- sum(d.x.y.S1 & dmS2.1 == l)
          }
        }
      }

      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0, c(nl.x,nl.y,prod(nl.S)))
      for (k in 1:prod(nl.S))
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k],njk[,k]))

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    } ,
    { ## lenS == 3 ----------------------------------
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      in.S3 <- seq_len(nl.S3 <- nl.S[3])
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      dmS3.1 <- dm[,S[3]] + 1L
      nijk <- array(0L, c(nl.x,nl.y,prod(nl.S)))
      for(i1 in 1:nl.x) {
        d.x.i1 <- d.x1 == i1
        for(i2 in i.ny) {
          d.x.i1.y.i2 <- d.x.i1 & d.y1 == i2
          for(i3 in in.S1) {
            d.x.y.S1 <- d.x.i1.y.i2 & dmS1.1 == i3
            for(i4 in in.S2) {
              d.x.y.S1.2 <- d.x.y.S1 & dmS2.1 == i4
              for(i5 in in.S3)
                nijk[i1,i2, nl.S3*nl.S2*(i3-1)+nl.S3*(i4-1)+i5] <-
                   sum(d.x.y.S1.2 & dmS3.1 == i5)
            }
          }
        }
      }

      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0, c(nl.x,nl.y,prod(nl.S)))
      for (k in 1:prod(nl.S))
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k],njk[,k]))

      t.G2 <- 2 * nijk * log(t.log)
      t.G2[which(is.nan(t.G2),arr.ind = TRUE)] <- 0
      G2 <- sum(t.G2)
    } ,
    { ## lenS == 4 ----------------------------------
      in.S1 <- seq_len(nl.S1 <- nl.S[1])
      in.S2 <- seq_len(nl.S2 <- nl.S[2])
      in.S3 <- seq_len(nl.S3 <- nl.S[3])
      in.S4 <- seq_len(nl.S4 <- nl.S[4])
      dmS1.1 <- dm[,S[1]] + 1L
      dmS2.1 <- dm[,S[2]] + 1L
      dmS3.1 <- dm[,S[3]] + 1L
      dmS4.1 <- dm[,S[4]] + 1L
      nijk <- array(0L, c(nl.x, nl.y, prod(nl.S)))
      for(i1 in 1:nl.x) {
        d.x.i1 <- d.x1 == i1
        for(i2 in i.ny) {
          d.x.i1.y.i2 <- d.x.i1 & d.y1 == i2
          for(i3 in in.S1) {
            d.x.y.S1 <- d.x.i1.y.i2 & dmS1.1 == i3
            for(i4 in in.S2) {
              d.x.y.S1.2 <- d.x.y.S1 & dmS2.1 == i4
              for(i5 in in.S3) {
                d.x.y.S1.2.3 <- d.x.y.S1.2 & dmS3.1 == i5
                for(i6 in in.S4)
                  nijk[i1,i2, nl.S4*nl.S3*nl.S2*(i3-1) +
                              nl.S4*nl.S3*(i4-1)+
                              nl.S4*(i5-1)+i6] <- sum(d.x.y.S1.2.3 & dmS4.1 == i6)
              }
            }
          }
        }
      }
      nik <- apply(nijk, 3, rowSums)
      njk <- apply(nijk, 3, colSums)
      nk <- colSums(njk)

      ## compute G^2
      t.log <- array(0, c(nl.x,nl.y,prod(nl.S)))
      for (k in 1:prod(nl.S))
        t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k],njk[,k]))
      t.G2 <- 2 * nijk * log(t.log)
      t.G2[is.nan(t.G2)] <- 0
      G2 <- sum(t.G2)
    }) # end lens == 4, end{switch}
  }
  else { #  |S| = lenS >= 5  (gSquareBin: lenS >= 6)
    nijk <- array(0L, c(nl.x, nl.y, 1L))
    ## first sample 'by hand' to avoid if/else in the for-loop
    i <- d.x1[1]
    j <- d.y1[1]
    ## create directly a list of all k's  -- MM{FIXME}: there must be a better way
    k <- NULL
    lapply(as.list(S), function(x) { k <<- cbind(k,d.x1); NULL })
    ## first set of subset values
    parents.count <- 1L ## counter variable corresponding to the number
    ## of value combinations for the subset varibales
    ## observed in the data
    parents.val <- t(k[1,])
    nijk[i,j,parents.count] <- 1L # cell counts

    ## Do the same for all other samples. If there is already a table
    ## for the subset values of the sample, increase the corresponding
    ## cell count. If not, create a new table and set the corresponding
    ## cell count to 1.
    for (it.sample in 2:n) {
      flag <- 0
      i <- d.x1[it.sample]
      j <- d.y1[it.sample]
      ## comparing the current values of the subset variables to all
      ## already existing combinations of subset variables values
      t.comp <- t(parents.val[1:parents.count,]) == k[it.sample,]
      ## Have to be careful here. When giving dimension to a list,
      ## R fills column after column, and NOT row after row.
      dim(t.comp) <- c(lenS,parents.count)
      for (it.parents in 1:parents.count) {
        ## check if the present combination of value alreay exists
        if(all(t.comp[,it.parents])) {
          ## if yes, increase the corresponding cell count
          nijk[i,j,it.parents] <- nijk[i,j,it.parents] + 1L
          flag <- 1
          break
        }
      }# end for(it.parents...)
      ## if the combination of subset values is new...
      if (flag == 0) {
        ## ...increase the number of subset 'types'
        parents.count <- parents.count + 1L
        if (verbose >= 2)
          cat(sprintf(' adding new parents (count = %d) at sample %d\n',
                      parents.count, it.sample))
        ## ...add the new subset to the others
        parents.val <- rbind(parents.val,k[it.sample,])
        ## ...update the cell counts (add new array)
        nijk <- abind(nijk, array(0L, c(nl.x,nl.y,1)))
        nijk[i,j,parents.count] <- 1L
      } # end if(flag==0)
    } ## for(it in 2:n)
    if (verbose && verbose < 2)
      cat(sprintf(" added a total of %d new parents\n", parents.count))

    nik <- apply(nijk,3,rowSums)
    njk <- apply(nijk,3,colSums)
    nk <- colSums(njk)
    ## compute G^2
    t.log <- array(0,c(nl.x,nl.y,parents.count))
    for (k in 1:parents.count)
      t.log[,,k] <- nijk[,,k]*(nk[k] / tcrossprod(nik[,k],njk[,k]))
    t.G2 <- 2 * nijk * log(t.log)
    t.G2[which(is.nan(t.G2),arr.ind = TRUE)] <- 0
    G2 <- sum(t.G2)
  }

  if (adaptDF && lenS > 0) {
    ## lower the degrees of freedom according to the amount of
    ## zero counts
    zero.counts <-
      if(lenS == 0)
        length(which(nij == 0))
      else
        length(which(nijk == 0)) + 4*(2^lenS - dim(nijk)[3])
    ## add zero counts corresponding to the number of parents
    ## combinations that are missing
    df <- max((df-zero.counts),1)
  } # end adaptDF

  pchisq(G2, df, lower.tail = FALSE)# i.e. == 1 - P(..)

}## gSquareDis()


gaussCItest <- function(x,y,S,suffStat) {
  ## suffStat$C: correlation matrix
  ## suffStat$n: sample size
  z <- zStat(x,y,S, C = suffStat$C, n = suffStat$n)
  2*pnorm(abs(z), lower.tail = FALSE)
}


dsep <- function(a,b, S = NULL, g, john.pairs = NULL)
{
  ## Purpose: Are the set a and the set b d-separeted given the set S?
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - a,b,S: vectors of node names
  ## - g: graphNEL object
  ## - john.pairs: matrix from johnson.all.pairs.sp
  ## ----------------------------------------------------------------------
  ## Value:
  ## Boolean decision
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch

  ## Check that g is a DAG
  amatTmp <- wgtMatrix(g)
  amatTmp[amatTmp != 0] <- 1
  if (max(amatTmp+t(amatTmp)) > 1) stop("dsep: Undirected edge in input graph!")
  p <- numNodes(g)
  ## build node union of a,b,S
  if(is.null(john.pairs)) john.pairs <- johnson.all.pairs.sp(g)
  nodeUnion <- if(length(S) > 0) c(a,b,S) else c(a,b)
  my.nodes <- nodes(g)

  ## find ancestor graph of nodeUnion
  anc.set <- NULL
  for (i in seq_len(p)) {
    desc.nodes <- my.nodes[which(john.pairs[i,] < Inf)]
    if (any(desc.nodes %in% nodeUnion)) anc.set <- c(anc.set,my.nodes[i])
  } ## for (i in 1:p)
  gS <- subGraph(anc.set,g)

  ## Moralize in amatM
  amat <- wgtMatrix(gS, transpose = FALSE)
  if(all(a0 <- amat == 0))
    ## if no edge in graph, nodes are d-separated
    return( TRUE )
  ## else :
  amat[!a0] <- 1
  amatM <- amat
  ind <- which(amat == 1,arr.ind = TRUE)

  for (i in seq_len(nrow(ind))) {
    ## input is guaranteed to be directed
    x <- ind[i,1]
    y <- ind[i,2] ## x -> y
    allZ <- setdiff(which(amat[y,] == 0 & amat[,y] == 1), x) ## x -> y <- z
    for (z in allZ)
      if (amat[x,z] == 0 && amat[z,x] == 0)
        amatM[x,z] <- 1 ## moralize
  } ## for (i in seq_len(nrow(ind)))

  ## make undirected graph
  ## (up to now, there is NO undirected edge -> just add t(amat))
  gSM <- as(amatM+t(amatM),"graphNEL")

  if (length(S) > 0) { ## check separation
    separates(a,b,S,gSM)
  } else {
    b %nin% (if(is.list(bfs. <- bfs(gSM,a))) bfs.[[1]] else bfs.)
  }
} ## {dsep}


## Orakel:
dsepTest <- function(x,y, S = NULL, suffStat) {
  ## suffStat$ g: True graph (graphNEL suffStatect)
  ## suffStat$jp: johnson all pairs
  ## Returns "p-value" P =
  ##	0: keep edge / d-connected
  ##    1: drop edge / d-separated

  if( x == y || x %in% S || y %in% S) {
    0
  } else {
    stopifnot(is(g <- suffStat$g, "graph"))
    jp <- suffStat$jp
    V <- nodes(g)
    ## return  0 / 1
    as.numeric(dsep(a = V[x], b = V[y], S = V[S], g = g, john.pairs = jp))
  }
} ## {dsepTest}

disCItest <- function(x,y,S,suffStat) {
  if(is.data.frame(dm <- suffStat$dm)) dm <- data.matrix(dm)
  else stopifnot(is.matrix(dm))
  nlev <- suffStat$nlev
  adaptDF <- suffStat$adaptDF
  ## p-value:
  gSquareDis(x = x, y = y, S = S, dm = dm, nlev = nlev, adaptDF = adaptDF,
             verbose = FALSE)
}

binCItest <- function(x,y,S,suffStat) {
  if(is.data.frame(dm <- suffStat$dm)) dm <- data.matrix(dm)
  else stopifnot(is.matrix(dm))
  adaptDF <- suffStat$adaptDF
  ## p-value:
  gSquareBin(x = x, y = y, S = S, dm = dm, adaptDF = adaptDF, verbose = FALSE)
}



ida <- function(x.pos, y.pos, mcov, graphEst, method = c("local","global"),
                y.notparent = FALSE, verbose = FALSE, all.dags = NA)
{
  ## Purpose: Estimate the causal effect of x on y; the graphEst and correlation
  ## matrix have to be precomputed; all DAGs can be precomputed;
  ## Orient undirected edges at x in a way so that no new collider
  ## is introduced
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## - method: "local" - local (all combinations of parents in regr.)
  ##           "global" - all DAGs
  ## - y.notparent: if TRUE, the effect of x <- y is ignored;
  ##                (remove y from all parents set pa1 or pa2)
  ##                if FALSE, the effect of x <- y is set to zero
  ## - verbose: if TRUE, details on regressions that were used
  ## - all.dags: All DAGs in the format of function allDags; if this is
  ##   available, no new function call allDags is done
  ## ----------------------------------------------------------------------
  ## Value: causal values
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010; tweaks: Martin Maechler

  stopifnot(x.pos == (x <- as.integer(x.pos)),
            y.pos == (y <- as.integer(y.pos)),
            length(x) == 1, length(y) == 1)
  method <- match.arg(method)

  ## prepare adjMatrix and skeleton
  amat <- ad.g <- wgtMatrix(graphEst)
  amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1

  if (method == "local") {
##############################
    ## local method
    ## Main Input: mcov, graphEst
##############################
    ## find unique parents of x
    wgt.est <- (ad.g != 0)
    if (y.notparent) {
      ## Direct edge btw. x and y towards y
      wgt.est[x, y] <- FALSE
    }
    tmp <- wgt.est-t(wgt.est)
    tmp[which(tmp < 0)] <- 0
    wgt.unique <- tmp
    pa1 <- which(wgt.unique[x,] != 0)
    if (y %in% pa1) {
      ## x is parent of y -> zero effect
      beta.hat <- 0
    } else { ## y not in pa1
      ## find ambiguous parents of x
      wgt.ambig <- wgt.est-wgt.unique
      pa2 <- which(wgt.ambig[x,] != 0)
      if (verbose)
        cat("\n\nx=", x, "y=",y, "\npa1=",pa1, "\npa2=",pa2,"\n")

      ## estimate beta
      if (length(pa2) == 0) {
        beta.hat <- lm.cov(mcov,y,c(x,pa1))
        if (verbose) cat("Fit - y:",y,"x:",c(x,pa1), "|b.hat=",beta.hat,"\n")
      } else {
        ## at least one undirected parent
        beta.hat <- NA
        ii <- 1

        ## no member of pa2
        pa2.f <- pa2
        pa2.t <- NA
        if (!has.new.coll(amat,amatSkel,x,pa1,pa2.t,pa2.f)) {
          beta.hat[ii] <- lm.cov(mcov,y,c(x,pa1))
          if (verbose) cat("Fit - y:",y,"x:",c(x,pa1),
                           "|b.hat=",beta.hat[ii],"\n")
        }
        ## exactly one member of pa2
        for (i2 in seq_along(pa2)) {
          pa2.f <- pa2[-i2]
          pa2.t <- pa2[i2]
          if (!has.new.coll(amat,amatSkel,x,pa1,pa2.t,pa2.f)) {
            ii <-  ii+1
            if (y %in% pa2.t) {
              beta.hat[ii] <- 0
            } else {
              beta.hat[ii] <- lm.cov(mcov,y,c(x,pa1,pa2[i2]))
              if (verbose) cat("Fit - y:",y,"x:",c(x,pa1,pa2[i2]),
                               "|b.hat=",beta.hat[ii],"\n")
            }
          } ## if (!check..)
        } ## for (i2 in seq_along(pa2))

        ## higher order subsets of pa2
        if (length(pa2) > 1)
          for (i in 2:length(pa2)) {
            pa.tmp <- combn(pa2,i,simplify = TRUE)
            for (j in seq_len(ncol(pa.tmp))) {
              pa2.t <- pa.tmp[,j]
              pa2.f <- setdiff(pa2,pa2.t)
              if (!has.new.coll(amat,amatSkel,x,pa1,pa2.t,pa2.f)) {
                ii <- ii+1
                if (y %in% pa2.t) {
                  beta.hat[ii] <- 0
                } else {
                  beta.hat[ii] <- lm.cov(mcov,y,c(x,pa1,pa2.t))
                  if (verbose)
                    cat("Fit - y:",y,"x:",c(x,pa1,pa2.t),
                        "|b.hat=",beta.hat[ii],"\n")
                }
              } ## if (!check..)
            } ## for (j ...)
          } ## for (i ...)
      } ## if (length(pa2) ..)
    } ## if (y %in% pa1)

  } else {
##############################
    ## global method
    ## Main Input: mcov, graphEst
##############################
    p <- numNodes(graphEst)
    am.pdag <- ad.g
    am.pdag[am.pdag != 0] <- 1
    if (y.notparent) {
      ## Direct edge btw. x and y towards y
      am.pdag[x, y] <- 0
    }

    ## find all DAGs if not provided externally
    ad <- if(is.na(all.dags)) allDags(am.pdag,am.pdag,NULL) else all.dags
    n.dags <- nrow(ad)
    beta.hat <- rep(NA,n.dags)
    for (i in 1:n.dags) {
      ## compute effect for every DAG
      ## gDag <- as(matrix(ad[i,],p,p),"graphNEL")
      ## path from y to x
      ## rev.pth <- RBGL::sp.between(gDag,as.character(y),
      ##                    as.character(x))[[1]]$path
      ## if (length(rev.pth)>1) {
      ## if reverse path exists, beta=0
      ##  beta.hat[i] <- 0
      ## } else {
      ## path from x to y
      ##       pth <- RBGL::sp.between(gDag,as.character(x),
      ##                       as.character(y))[[1]]$path
      ##   if (length(pth)<2) {
      ## sic! There is NO path from x to y
      ##   beta.hat[i] <- 0
      ## } else {
      ## There is a path from x to y
      wgt.unique <- t(matrix(ad[i,],p,p)) ## wgt.est is wgtMatrix of DAG
      pa1 <- which(wgt.unique[x,] != 0)
      if (y %in% pa1) {
        beta.hat[i] <- 0
      } else {
        beta.hat[i] <- lm.cov(mcov, y, c(x,pa1))
        if (verbose) cat("Fit - y:",y,"x:",c(x,pa1),
                         "|b.hat=",beta.hat[i],"\n")
      }
    } ## for ( i  n.dags)
  } ## else : method = "global"
  beta.hat
} ## {ida}

idaFast <- function(x.pos, y.pos.set, mcov, graphEst)
{
  ## Purpose: Estimate the causal effect of x on each element in the
  ## set y using the local method; graphEst and correlation matrix
  ## have to be precomputed; orient
  ## undirected edges at x in a way so that no new collider is
  ## introduced; if there is an undirected edge between x and y, both directions are considered;
  ## i.e., y might be partent of x in which case the effect is 0.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x.pos, y.pos: Column of x and y in d.mat
  ## - mcov: Covariance matrix that was used to estimate graphEst
  ## - graphEst: Fit of PC Algorithm (semidirected)
  ## ----------------------------------------------------------------------
  ## Value: list of causal values; one list element for each element of
  ## y.pos.set
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 7 Jan 2010, 11:18

  ## prepare adjMatrix and skeleton
  amat <- ad.g <- wgtMatrix(graphEst)
  amat[which(amat != 0)] <- 1 ## i->j if amat[j,i]==1
  amatSkel <- amat + t(amat)
  amatSkel[amatSkel != 0] <- 1

  ## find unique and ambiguous parents of x
  wgt.est <- (ad.g != 0)
  tmp <- wgt.est-t(wgt.est)
  tmp[which(tmp < 0)] <- 0
  wgt.unique <- tmp
  wgt.ambig <- wgt.est-wgt.unique
  pa1 <- which(wgt.unique[x.pos,] != 0)
  pa2 <- which(wgt.ambig[x.pos,] != 0)

  ## estimate beta
  if (length(pa2) == 0) {
    beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
    beta.tmp[y.pos.set %in% pa1] <- 0
    beta.hat <- cbind(beta.tmp)
  } else {    ## at least one undirected parent
    ## no member of pa2
    pa2.f <- pa2
    pa2.t <- NA
    beta.hat <-
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
	beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1)) ####
	beta.tmp[y.pos.set %in% pa1] <- 0
	cbind(beta.tmp)
      } # else NULL

    ## exactly one member of pa2
    for (i2 in seq_along(pa2)) {
      pa2.f <- pa2[-i2]
      pa2.t <- pa2[i2]
      if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
        beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
        beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
        beta.hat <- cbind(beta.hat, beta.tmp)
      }
    } ## for (i2 in seq_along(pa2))

    ## higher order subsets of pa2
    if (length(pa2) > 1)
      for (i in 2:length(pa2)) {
        pa.tmp <- combn(pa2,i,simplify = TRUE)
        for (j in seq_len(ncol(pa.tmp))) {
          pa2.t <- pa.tmp[,j]
          pa2.f <- setdiff(pa2, pa2.t)
          if (!has.new.coll(amat,amatSkel,x.pos,pa1,pa2.t,pa2.f)) {
            beta.tmp <- lm.cov(mcov,y.pos.set,c(x.pos,pa1,pa2.t)) ####
            beta.tmp[y.pos.set %in% c(pa1,pa2.t)] <- 0
            beta.hat <- cbind(beta.hat, beta.tmp)
          }
        } ## for (j )
      } ## for (i )

  } ## if .. else length(pa2) > 0)

  ## MM: for now, maybe in the future get sensible column names:
  colnames(beta.hat) <- NULL
  if (nrow(beta.hat) > 0) rownames(beta.hat) <- as.character(y.pos.set)
  beta.hat
}

## -> ../man/legal.path.Rd
## only called in  qreach()  with only 'c' varying
legal.path <- function(a,b,c, amat)
{
  ## Purpose: Is path a-b-c legal (either collider in b or a,b,c is triangle)
  ## !! a-b-c must be in a path !! this is not checked !!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - a, b, c: nodes
  ## - amat: adj matrix (coding 0,1,2 for no edge, circle, arrowhead)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 29 Oct 2009; Martin Maechler
  if(a == c || (a.b <- amat[a,b]) == 0 || amat[b,c] == 0)
    return(FALSE)
  ## else  a != c  and  amat[a,b] != 0  and   amat[b,c] != 0
  ## return  TRUE iff
  (amat[a,c] != 0 || ## triangle
   ## need not check [c,a], since there must be SOME edgemark !=0 at [a,c], if
   ## edge is present
   (a.b == 2 && amat[c,b] == 2)) ## a collider
}



##' Plot a subgraph for a specified starting node and a given graph
##' @param graphObj Graph object
##' @param y        Starting node
##' @param dist     Distance of nodes included in subgraph from starting node y
##' @param amat     Adjacency matrix of skeleton graph (optional)
##' @param directed Boolean; should the plotted subgraph be directed?
##' @param main
##' --------------- see ../man/plotSG.Rd
plotSG <- function(graphObj, y, dist, amat = NA, directed = TRUE,
                   main = paste("Subgraph of ", deparse(substitute(graphObj)),
                     "\nfrom ", y, " with distance ", dist ))
{
  ## Author: Daniel Stekhoven, Date: 26 Jan 2010, 08:56

  check.Rgraphviz()
  stopifnot(dist >= 1)

  ## Extract adjacency matrix (if necessary)
  if (any( is.na(amat)))
    amat <- wgtMatrix(graphObj)

  ## Diagonalise (has no effect if already diagonal)
  amat[amat != 0] <- 1
  amat <- amat + t(amat)
  diag(amat) <- 0 # can that happen anyway??
  amat[amat == 2] <- 1

  ## Find connected nodes hierarchically
  nh <- which( amat[y,] == 1 )
  rel.nodes <- c( y, nh )
  for (i in seq_len(dist-1L)) {
    nh <-
      if ( length(nh) == 1 )
        which( amat[nh,] == 1 )
      else if (length(nh) != 0)
        which(amat[nh,] == 1, arr.ind = TRUE)[,"col"]
    ## else NULL
    rel.nodes <- unique( c( rel.nodes, nh ) )
  }

  ## Name nodes
  if(is(graphObj, "graphNEL"))
    names(rel.nodes) <- graphObj@nodes[rel.nodes]

  ## subgraph - distinguish between directed edges or not
  sg <- if (directed)
    subGraph(as.character(rel.nodes), graphObj)
  else
    as(amat[rel.nodes, rel.nodes], "graphNEL")

  ## Plot subgraph
  Rgraphviz::plot( sg )
  if(!is.null(main))
    title(main = main)
  invisible(sg)
}


##########################################################################
##
## Functions for the conservative versions (PC, FCI (all), and RFCI)
##
#########################################################################

## Functions used by all algorithms
##########################################################
pc.cons.intern <- function(sk, suffStat, indepTest, alpha,
                           version.unf = c(NA,NA), maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose:  For any unshielded triple A-B-C, consider all subsets D of
  ## the neighbors of A and of the neighbors of C, and record the sets
  ## D for which A and C are conditionally independent given D. If B
  ## is in none of these sets, do nothing (it is a
  ## v-structure) and also delete B from sepset(A,C) if present (so we are
  ## sure that a v-structure will be created). If B is in all sets, do nothing
  ## (it is not a v-structure) and also add B to sepset(A,C) if not present
  ## (so we are sure that a v-structure will not be created). If maj.rule=FALSE
  ## the normal conservative version is applied, hence if B is in
  ## some but not all sets, mark the triple as "ambiguous". If maj.rule=TRUE
  ## we mark the triple as "ambiguous" if B is in exactly 50% of the cases,
  ## if less than 50% define it as a v-structure, and if in more than 50%
  ## no v-structure.
  ## ----------------------------------------------------------------------
  ## Arguments: - sk: output returned by function "skeleton"
  ##            - suffStat: Sufficient statistics for independent tests
  ##            - indepTest: Function for independence test
  ##            - alpha: Significance level of test
  ##            - version.unf[1]: 1 it checks if b is in some sepsets,
  ##                              2 it also checks if there exists a sepset
  ##                              which is a subset of the neighbours.
  ##            - version.unf[2]: 1 same as in Tetrad (do not consider
  ##                              the initial sepset), 2 it also considers
  ##                              the initial sepset
  ##            - maj.rule: FALSE/TRUE if the majority rule idea is applied
  ## ----------------------------------------------------------------------
  ## Value: - unfTripl: Triple that were marked as unfaithful
  ##        - vers: vector containing the version (1 or 2) of the
  ##                corresponding triple saved in unfTripl (1=normal
  ##                unfaithful triple that is B is in some sepsets;
  ##                2=triple coming from version.unf[1]==2
  ##                that is a and c are indep given the initial sepset
  ##                but there doesn't exist a subset of the neighbours
  ##                that d-separates them)
  ##        - sk: updated skelet object, sepsets might have been updated
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 10:43
  ## Modifications: Diego Colombo

  g <- as(sk@graph,"matrix")
  stopifnot(all(g == t(g))) ## g is guaranteed to be symmetric
  p <- as.numeric(dim(g)[1])
  unfTripl <- vers <- rep(NA,min(p*p,100000))
  counter <- 0
  if (sum(g) > 0) {
    ind <- which(g == 1, arr.ind = TRUE)
    tripleMatrix <- NULL
    ## Go through all edges
    for (i in seq_len(nrow(ind))) {
      a <- ind[i,1]
      b <- ind[i,2]
      allC <- setdiff(which(g[b,] == 1),a) ## a-b-c
      newC <- allC[g[a,allC] == 0]
      tmpMatrix <- cbind(rep(a,length(newC)),rep(b,length(newC)),newC)
      tripleMatrix <- rbind(tripleMatrix,tmpMatrix)
      colnames(tripleMatrix) <- c("","","")
    }
    if ((m <- nrow(tripleMatrix)) > 0) {
      deleteDupl <- logical(m)# all FALSE
      for (i in seq_len(m))
        if (tripleMatrix[i,1] > tripleMatrix[i,3])
          deleteDupl[i] <- TRUE
      if(any(deleteDupl))
        tripleMatrix <- tripleMatrix[!deleteDupl,, drop = FALSE]

      for (i in seq_len(nrow(tripleMatrix))) {
        ## pay attention to the size of counter
        if (counter+1L == length(unfTripl)) {
          n.xtra <- min(p*p, 100000)
          new.len <- counter+1L + n.xtra
          length(unfTripl) <- new.len
          length(vers)     <- new.len
        }
        a <- tripleMatrix[i,1]
        b <- tripleMatrix[i,2]
        c <- tripleMatrix[i,3]
        nbrsA <- which(g[,a] != 0) ## G symm; c no nbr of a
        nbrsC <- which(g[,c] != 0)
        if (verbose) {
          cat("\nTriple:", a,b,c,"and sepset by skelet:",
              unique(sk@sepset[[a]][[c]],sk@sepset[[c]][[a]]),"\n")
        }
        r.abc <- checkTriple(a, b, c, nbrsA, nbrsC,
                             sk@sepset[[a]][[c]], sk@sepset[[c]][[a]],
                             suffStat = suffStat, indepTest = indepTest, alpha = alpha,
                             version.unf = version.unf, maj.rule = maj.rule, verbose = verbose)
        ## 1: in NO set; 2: in ALL sets; 3: in SOME but not all
        ## Take action only if case "3"
        if (r.abc$decision == 3) {
          ## record ambiguous triple
          counter <- counter + 1
          unfTripl[counter] <- triple2numb(p,a,b,c)
          vers[counter] <- r.abc$version
        }
        ## can happen the case in Tetrad, so we must save the triple
        ## as ambiguous:
        ## a and c independent given S but not given subsets of the
        ## adj(a) or adj(c)
        if ((version.unf[1] == 2) && (r.abc$version == 2) && (r.abc$decision != 3)) {
          counter <- counter + 1
          unfTripl[counter] <- triple2numb(p,a,b,c)
          vers[counter] <- r.abc$version
        }
        sk@sepset[[a]][[c]] <- r.abc$SepsetA
        sk@sepset[[c]][[a]] <- r.abc$SepsetC
      }
    }
  }
  length(unfTripl) <- length(vers) <- counter
  list(unfTripl = unfTripl, vers = vers, sk = sk)
}


## Called both from pc.cons.intern() and rfci.vStruc() :
checkTriple <- function(a, b, c, nbrsA, nbrsC, sepsetA, sepsetC,
                        suffStat, indepTest, alpha, version.unf = c(NA,NA),
                        maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose: For each subset of nbrsA and nbrsC where a and c are cond.
  ## independent, it is checked if b is in the conditioning set.
  ## ----------------------------------------------------------------------
  ## Arguments: - a,b,c: Nodes (positions in adjacency matrix)
  ##            - nbrsA: Neighbors of a
  ##            - nbrsC: Neighbors of c
  ##            - sepsetA: sepset(a,c)
  ##            - sepsetC: sepset(c,a)
  ##            - suffStat: Sufficient statistics for independent tests
  ##            - indepTest: Function for independence test
  ##            - alpha: Significance level of test
  ##            - version.unf[1]: 1 it checks if b is in some sepsets,
  ##                              2 it also checks if there exists a sepset
  ##                              which is a subset of the neighbours.
  ##            - version.unf[2]: 1 same as Tetrad (do not consider the initial
  ##                              sepset), 2 consider if b is in sepsetA
  ##                              or sepsetC
  ##            - maj.rule: FALSE/TRUE if the majority rule idea is applied
  ## ----------------------------------------------------------------------
  ## Value: - decision: res
  ##          res = 1: b is in NO sepset (-> v-structure)
  ##          res = 2: b is in ALL sepsets (-> no v-structure)
  ##          res = 3: b is in SOME but not all sepsets (-> ambiguous triple)
  ##        - version: version (1 or 2) of the ambiguous triple
  ##                (1=normal ambiguous triple that is b is in some sepsets;
  ##                2=triple coming from version.unf[1]==2 that is a and c are
  ##                indep given the initial sepset but there doesn't exist a
  ##                subset of the neighbours that d-separates them)
  ##        - sepsetA and sepsetC: updated separation sets
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 12 Feb 2010, 12:13
  ## Modifications: Diego Colombo, Martin Maechler


  ## loop through all subsets of parents

  nr.indep <- 0
  stopifnot(length(version.unf) == 2, version.unf %in% 1:2)
  ## Tetrad
  tmp <- if (version.unf[2] == 2) ## our version
    (b %in% sepsetA || b %in% sepsetC) ## else NULL = Tetrad version
  version <- 0
  ## start with the neighbours of a
  if ((nn <- length(nbrsA)) > 0) {
    allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
    ## loop through all subsets of neighbours
    for (i in 1:nrow(allComb)) { ## == 1:(2^nn)
      S <- nbrsA[which(allComb[i,] != 0)]
      pval <- indepTest(a, c, S, suffStat)
      ## save the pval and the set that produced this pval
      if (verbose) cat("a: S =",S," - pval =",pval,"\n")
      if (pval >= alpha) {
        nr.indep <- nr.indep + 1
        ## is b in set?
        tmp <- c(tmp, b %in% S)
        version <- 1
      }
    }
  }
  ## now with the neighbours of c
  if ((nn <- length(nbrsC)) > 0) {
    allComb <- expand.grid(lapply(integer(nn), function(.) 0:1))
    ## loop through all subsets of neighbours
    for (i in 1:nrow(allComb)) { ## == 1:(2^nn)
      S <- nbrsC[which(allComb[i,] != 0)]
      pval <- indepTest(a, c, S, suffStat)
      ## save the pval and the set that produced this pval
      if (verbose) cat("c: S =",S," - pval =",pval,"\n")
      if (pval >= alpha) {
        nr.indep <- nr.indep + 1
        ## is b in set?
        tmp <- c(tmp, b %in% S)
        version <- 1
      }
    }
  }
  if (version.unf[1] == 2  && nr.indep == 0) {
    version <- 2
  }
  if (is.null(tmp)) tmp <- FALSE

  if (all(tmp)) {
    res <- 2 ## in ALL sets
    ## therefore a - b - c is not a v-structure, hence add b to sepset(a,c)
    ## and sepset(c,a)
    ## for example it can happen that b is not in sepset(a,c) or sepset(c,a)
    ## but now it is in each set that separates a and c given the neighbours
    if (b %nin% sepsetA) sepsetA <- c(sepsetA, b)
    if (b %nin% sepsetC) sepsetC <- c(sepsetC, b)
  } else {
    if (all(!tmp)) {
      res <- 1 ## in NO set
      ## therefore a - b - c is a v-structure, hence delete b from
      ## sepset(a,c) and sepset(c,a)
      ## for example it can happen that b is in sepset(a,c) or sepset(c,a)
      ## but now it is in no set that separates a and c given the neighbours
      sepsetA <- setdiff(sepsetA,b)
      sepsetC <- setdiff(sepsetC,b)
    } else {
      ## normal conservative PC, b is in some sets hence the triple
      ## is unfaithful
      if (!maj.rule) {
        res <- 3 ## in SOME sets
      } else {
        ## use the majority rule to test if the triple is faithful
        ## or not and then decide if it is a v-structure or not accordigly
        ## NEW check the percentage of b in the conditioning sets that
        ## make a and c independent
        if (sum(tmp)/length(tmp) < 0.5) {
          ## we accept that b is in NO set
          res <- 1 ## in NO set
          ## therefore a - b - c is a v-structure, hence delete b
          ## from sepset(a,c) and sepset(c,a)
          ## for example it can happen that b is in sepset(a,c) or
          ## sepset(c,a) but now it is in no set that
          ## separates a and c given the neighbours
          sepsetA <- setdiff(sepsetA,b)
          sepsetC <- setdiff(sepsetC,b)
        } else if (sum(tmp)/length(tmp) > 0.5) {
          ## we accept that b is in ALL set
          res <- 2 ## in ALL sets
          ## therefore a - b - c is not a v-structure, hence add b
          ## to sepset(a,c) and sepset(c,a)
          ## for example it can happen that b is not in sepset(a,c)
          ## or sepset(c,a) but now it is in each set that
          ## separates a and c given the neighbours
          if (b %nin% sepsetA) sepsetA <- c(sepsetA,b)
          if (b %nin% sepsetC) sepsetC <- c(sepsetC,b)
        } else if (sum(tmp)/length(tmp) == 0.5) {
          ## define the triple as unfaithful, because half of the
          ## times b is in the set and half of them in not in
          res <- 3 ## in SOME sets
        }
      }
    }
  }
  if (verbose && res == 3) cat("Triple ambiguous\n")

  ## if you save a variable <- NULL into a list it will delete this element!
  ## The following also transforms NULL sepset* to integer(0):
  lapply(list(decision = res, version = version, SepsetA = sepsetA, SepsetC = sepsetC),
         as.integer)
} ## {checkTriple}

## For R5 and R9-R10
######################################################
faith.check <- function(cp, unfVect, p, boolean = TRUE)
{
  ## Purpose: check if every triple on the circle path  cp is unambiguous
  ## ----------------------------------------------------------------------
  ## Arguments: cp: circle path to check for unambiguity
  ##  boolean:  if TRUE,  return TRUE iff there is no ambiguity, i.e. "faithful"
  ##            if FALSE, return the number of ambiguities.
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 May 2010; 'boolean' etc by Martin Maechler
  if(!boolean) res <- 0L
  n <- length(cp)
  ## MM: FIXME speedup by precomputing (l%%n)+1, ((l+1)%%n)+1, and  ((l+2)%%n)+1 as *integers*
  ## ok, in steps: first "slowly" but surely correct
  ii <- 0:(n-1L)
  i1 <- (ii      %% n)+1L
  i2 <- ((ii+1L) %% n)+1L
  i3 <- ((ii+2L) %% n)+1L
  for (l in ii) {
    if (any(unfVect == triple2numb(p,cp[i1[l]],cp[i2[l]],cp[i3[l]]), na.rm = TRUE) ||
        any(unfVect == triple2numb(p,cp[i3[l]],cp[i2[l]],cp[i1[l]]), na.rm = TRUE)) {
      if(boolean) return(FALSE) # the first time we found one: not faithful
      ## else count
      res <- res + 1L
    }
  }
  if(boolean) TRUE else res
}


###############################################################################
##
## FCI, RFCI, and fast FCI-oracle
##
###############################################################################


## Internal function used by FCI and RFCI
##############################################################################

triple2numb <- function(p,i,j,k)
{
  ## Purpose: transform a triple i-j-k into a unique number
  ## ----------------------------------------------------------------------
  ## Arguments:-p:number of nodes;-triple: i-j-k
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 11 May 2010;  Martin Maechler
  ## as.numeric(.): not integer arithmetic which easily overflows
  k + (p. <- as.numeric(p)) * (j + p.*i)
}

## Functions for R4, R5, and R9-R10 for FCI(all) and RFCI
############################################################################

updateList <- function(path, set, old.list)
{
  ## Purpose: update the list of all paths in the iterative functions
  ## minDiscrPath, minUncovCircPath and minUncovPdPath
  ## ----------------------------------------------------------------------
  ## Arguments: - path: the path under investigation
  ##            - set: (integer) index set of variables to be added to path
  ##            - old.list: the list to update
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2011; Without for() by Martin Maechler
  c(old.list, lapply(set, function(s) c(path,s)))
}

## R9-R10
minUncovPdPath <- function(p, pag, a,b,c, unfVect, verbose = FALSE)
{
  ## Purpose: find a minimal uncovered pd path for a,b,c saved in path.
  ## Check also for the conservative case that it is unambiguous
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - p: number of nodes in the graph
  ##            - pag: adjacency matrix
  ##            - a,b,c : nodes under interest
  ##            - unfVect: vector containing the ambiguous triples
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 19 Oct 2011; small changes: Martin Maechler

  visited <- rep(FALSE, p)
  visited[c(a,b,c)] <- TRUE
  min.upd.path <- NA
  ## find all neighbours of b not visited yet
  indD <- which((pag[b,] == 1 | pag[b,] == 2) &
                (pag[,b] == 1 | pag[,b] == 3) &
                (pag[,a] == 0) & !visited)
  if (length(indD) > 0) {
    path.list <- updateList(b, indD, NULL)
    done <- FALSE
    while ((length(path.list) > 0) && (!done)) {
      ## next element in the queue
      mpath <- path.list[[1]]
      m <- length(mpath)
      d <- mpath[m]
      path.list[[1]] <- NULL
      visited[d] <- TRUE
      if (any(pag[d,c] == 1:2) && any(pag[c,d] == c(1,3))) {
        ## pd path found
        mpath <- c(a, mpath, c)
        n <- length(mpath)
        ## check the path to be uncovered
        uncov <- TRUE
	for (l in seq_len(n - 2)) {
	  if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
		pag[mpath[l + 2], mpath[l]] == 0)) {

	    uncov <- FALSE
	    break ## speed up!
	  }
	}
        ## if it is uncovered
        if (uncov)
          if (length(unfVect) == 0 || ## <<- normal version: just save
              ## conservative version, check the path to be faithful:
              faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.upd.path <- mpath
            done <- TRUE
          }
      }
      else {
        ## d and c are either not connected or connected with a "wrong" edge -----> search iteratively
        ## find all neighbours of d not visited yet
        indR <- which((pag[d,] == 1 | pag[d,] == 2) &
                      (pag[,d] == 1 | pag[,d] == 3) & !visited)
        if (length(indR) > 0) {
          ## update the queues
          path.list <- updateList(mpath, indR, path.list)
        }
      }
    } ## {while}
  }
  min.upd.path
} ## {minUncovPdPath}

## R5
minUncovCircPath <- function(p, pag, path, unfVect, verbose = FALSE)
{
  ## Purpose: find a minimal uncovered circle path for a,c,d,b saved in path.
  ## Check also for the conservative case that it is unambiguous
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - p: number of nodes in the graph
  ##            - pag: adjacency matrix
  ##            - path: a,c,d,b under interest
  ##            - unfVect: vector containing the unfaithful triples
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 19 Oct 2011, 13:11

  visited <- rep(FALSE, p)
  visited[path] <- TRUE # (a,b,c,d) all 'visited'
  a <- path[1]
  c <- path[2]
  d <- path[3]
  b <- path[4]
  min.ucp.path <- NA
  ## find all neighbours of c not visited yet
  indX <- which(pag[c,] == 1 & pag[,c] == 1 & !visited) ## c o-o x
  if (length(indX) > 0) {
    path.list <- updateList(c, indX, NULL)
    done <- FALSE
    while (!done && length(path.list) > 0) {
      ## next element in the queue
      mpath <- path.list[[1]]
      x <- mpath[length(mpath)]
      path.list[[1]] <- NULL
      visited[x] <- TRUE
      if (pag[x,d] == 1 && pag[d,x] == 1) {
        ## circle path found
        mpath <- c(a, mpath, d, b)
        n <- length(mpath)
        ## check the path to be uncovered
        uncov <- TRUE
	for (l in seq_len(n - 2)) {
	  if (!(pag[mpath[l], mpath[l + 2]] == 0 &&
		pag[mpath[l + 2], mpath[l]] == 0)) {

	    uncov <- FALSE
	    break ## speed up!
	  }
	}
        ## if it is uncovered
        if (uncov)
          if (length(unfVect) == 0 || ## <<- normal version: just save
              ## conservative version, check the path to be faithful:
              faith.check(mpath, unfVect, p)) {
            ## save the path to be returned
            min.ucp.path <- mpath
            done <- TRUE
          }
      }
      else {
        ## x and d are either not connected or connected with an edge which is not o--o -----> search iteratively
        ## find all neighbours of x not visited yet
        indR <- which(pag[x,] == 1 & pag[,x] == 1 & !visited) ## x o--o r
        if (length(indR) > 0) {
          ## update the queues
          path.list <- updateList(mpath, indR, path.list)
        }
      }
    }## {while}
  }
  return(min.ucp.path)
}

## R4
minDiscrPath <- function(pag, a,b,c, verbose = FALSE)
{
  ## Purpose: find a minimal discriminating path for a,b,c.
  ## If a path exists this is the output, otherwise NA
  ## ----------------------------------------------------------------------
  ## Arguments: - pag: adjacency matrix
  ##            - a,b,c: node positions under interest
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 Jan 2011; speedup: Martin Maechler

  p <- as.numeric(dim(pag)[1])
  visited <- rep(FALSE, p)
  visited[c(a,b,c)] <- TRUE # {a,b,c} "visited"
  ## find all neighbours of a  not visited yet
  indD <- which(pag[a,] != 0 & pag[,a] == 2 & !visited) ## d *-> a
  if (length(indD) > 0) {
    path.list <- updateList(a, indD, NULL)
    while (length(path.list) > 0) {
      ## next element in the queue
      mpath <- path.list[[1]]
      m <- length(mpath)
      d <- mpath[m]
      if (pag[c,d] == 0 & pag[d,c] == 0)
	## minimal discriminating path found :
	return( c(rev(mpath), b,c) )

      ## else :
      pred <- mpath[m-1]
      path.list[[1]] <- NULL
      visited[d] <- TRUE

      ## d is connected to c -----> search iteratively
      if (pag[d,c] == 2 && pag[c,d] == 3 && pag[pred,d] == 2) {
        ## find all neighbours of d not visited yet
        indR <- which(pag[d,] != 0 & pag[,d] == 2 & !visited) ## r *-> d
        if (length(indR) > 0)
          ## update the queues
          path.list <- updateList(mpath, indR, path.list)
      }
    } ## {while}
  }
  ## nothing found:  return
  NA
} ## {minDiscrPath}


###############################################################################
## FCI
###############################################################################

fci <- function(suffStat, indepTest, alpha, labels, p,
                skel.method = c("stable", "original", "stable.fast"),
                type = c("normal", "anytime", "adaptive"),
                fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                m.max = Inf, pdsep.max = Inf, rules = rep(TRUE, 10),
                doPdsep = TRUE, biCC = FALSE, conservative = FALSE,
                maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose: Perform FCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximum size of conditioning set
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - rules: array of length 10 wich contains TRUE or FALSE corresponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - doPdsep: compute possible dsep
  ## - biCC: TRUE or FALSE variable containing if biconnected components are
  ##         used to compute pdsep
  ## - conservative: TRUE or FALSE defining if
  ##          the v-structures after the pdsep
  ##          have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes
  ## - type: it specifies the version of the FCI that has to be used.
  ##         Per default it is normal, the normal FCI algorithm. It can also be
  ##         anytime for the Anytime FCI and in this cas m.max must be specified;
  ##         or it can be adaptive for Adaptive Anytime FCI and in this case
  ##         m.max must not be specified.

  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: Dec 2009; update: Diego Colombo, 2012; Martin Maechler, 2013

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  ## Check that the type is a valid one
  type <- match.arg(type)
  if (type == "anytime" && m.max == Inf)
    stop("To use the Anytime FCI you must specify a finite 'm.max'.")
  if (type == "adaptive" && m.max != Inf)
    stop("To use the Adaptive Anytime FCI you must not specify 'm.max'.")

  if (conservative && maj.rule)
    stop("Choose either conservative FCI or majority rule FCI")

  cl <- match.call()
  if (verbose) cat("Compute Skeleton\n================\n")

  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  skel@call <- cl # so that makes it into result
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  tripleList <- NULL

  if (doPdsep) {
    if (verbose) cat("\nCompute PDSEP\n=============\n")
    pc.ci <- pc.cons.intern(skel, suffStat, indepTest,
                            alpha = alpha, version.unf = c(1,1),
                            maj.rule = FALSE, verbose = verbose)
    ## Recompute (sepsets, G, ...):
    pdsepRes <- pdsep(skel@graph, suffStat, indepTest = indepTest, p = p,
                      sepset = pc.ci$sk@sepset, alpha = alpha, pMax = pMax,
                      m.max = if (type == "adaptive") max.ordSKEL else m.max,
                      pdsep.max = pdsep.max, NAdelete = NAdelete,
                      unfVect = pc.ci$unfTripl, # "tripleList.pdsep"
                      biCC = biCC, verbose = verbose)

    ## update the graph & sepset :
    G <- pdsepRes$G
    sepset <- pdsepRes$sepset
    pMax <- pdsepRes$pMax
    allPdsep <- pdsepRes$allPdsep
    n.edgetestsPD <- pdsepRes$n.edgetests
    max.ordPD <- pdsepRes$max.ord
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      tmp.pdsep <- new("pcAlgo", graph = as(G, "graphNEL"), call = cl,
                       n = integer(0), max.ord = as.integer(max.ordSKEL),
                       n.edgetests = n.edgetestsSKEL, sepset = sepset,
                       pMax = pMax, zMin = matrix(NA, 1, 1))
      sk. <- pc.cons.intern(tmp.pdsep, suffStat, indepTest, alpha,
                            verbose = verbose, version.unf = c(1, 1),
                            maj.rule = maj.rule)
      tripleList <- sk.$unfTripl
      ## update the sepsets
      sepset <- sk.$sk@sepset
    }
  }
  else {## !doPdsep : "do not Pdsep"
    n.edgetestsPD <- 0
    max.ordPD <- 0
    allPdsep <- vector("list", p)
    if (conservative || maj.rule) {
      if (verbose)
        cat("\nCheck v-structures conservatively\n=================================\n")
      nopdsep <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                                verbose = verbose, version.unf = c(2, 1),
                                maj.rule = maj.rule)
      tripleList <- nopdsep$unfTripl
      ## update the sepsets
      sepset <- nopdsep$sk@sepset
    }
  }
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  res <- udag2pag(pag = G, sepset, rules = rules, unfVect = tripleList,
                  verbose = verbose)
  colnames(res) <- rownames(res) <- labels
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),
      max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} ## {fci}

## only called in pdsep()
qreach <- function(x,amat,verbose = FALSE)
{
  ## Purpose: Compute possible-d-sep(x) ("psep")
  ## !! The non-zero entries in amat must be symmetric !!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x: node of which psep berechnet werden soll
  ## - amat: adjacency matrix
  ##         amat[i,j] = 0 iff no edge btw i,j
  ##         amat[i,j] = 1 iff i *-o j
  ##         amat[i,j] = 2 iff i *-> j
  ## - verbose: Show checked node sequence
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 29 Oct 2009, 11:54
  ## Stopping:
  ## =========
  ## At every iteration, Q get's reduced by one. It is only increased by
  ## at least one, if there are edges in A. and then, at least one
  ## edge in A. is removed. Edges are never inserted into A..
  ## Thus, either Q or A. becomes empty and the loop stops.
  ## Runtime:
  ## ========
  ## At least O(|V|), since look up in adjacency matrix is made. Assume O(|E|)>O
  ## Every edge can be visited at most twice. At each visit, there are no
  ## more than max(deg(V_i)) neighboring edges to look at. Thus, the runtime is
  ## O(2*|E| * max(deg(V_i))) = O(|E|^2) [worst case]; O(|E|) if sparse in the
  ## sense that max(deg(V_i)) is constant.
  ## Correctness:
  ## ============
  ## (A) All nodes in PSEP have a path of legal triples from x.
  ## (B) If there is a legal path from x to y in amat, at least one of them
  ## is found and y is recorded in PSEP:
  ## Suppose there is a node y != x that has a legal path from x, but is not in
  ## PSEP. y cannot be in nbrs(x), because they are added and nothing is
  ## deleted from PSEP. Hence, there must be a node z which
  ## is in PSEP but has a neighbor w in amat that has a legal path from
  ## x but is not in PSEP.
  ## Assuming that the function legal is correct, and noting that at (*) all
  ## neighbors of z in A. (tmp!) are analyzed, it follows that w is not
  ## in adj(z) in A. (but in amat). Thus, w must have been removed
  ## before from adj(z). Because of (+), it could only have been removed if
  ## u-z-w was found legal at some point. But then, w would have been added
  ## to PSEP. This is a contradiction to the assumption that w is not in PSEP.

  ## check: quadratic; x in V; edgemarks ok; non-zeroes symmetric
  stopifnot((ncol(amat) == nrow(amat)),x <= ncol(amat),all(amat %in% c(0,1,2)),
            all((amat != 0) == (t(amat != 0))))

  A. <- (amat != 0) ## A.[i,j] is true  <===>  edge i--j
  PSEP <- Q <- nb <- which(A.[x,])
  P <- rep.int(x, length(Q))
  A.[x,nb] <- FALSE ## delete edge to nbrs

  while(length(Q) > 0) {
    ## Invariants:
    ## ===========
    ## (A1) length(Q) == length(P) > 0
    ## (A2) non-zero in A. -> non-zero in amat [no non-zero added]
    ## (A3) Q[i] and P[i] are adjacent in amat [using (A2)]
    if (verbose) {
      cat("\n-------------","\n")
      cat("Queue Q:",Q,"\n")
      cat("Queue P:",P,"\n")
    }
    a <- Q[1]
    Q <- Q[-1]
    pred <- P[1] ## not empty because of (A1)
    P <- P[-1]
    if (verbose) cat("Select",pred,"towards",a,"\n")
    nb <- which(A.[a,]) ## (*)
    if (verbose) cat("Check nbrs",nb,"\nLegal:")

    for (b in nb) {
      ## Guaranteed: pred-a-b are a path because of (A3)
      if (lres <- legal.path(pred,a,b,amat)) {
        A.[a,b] <- FALSE ## remove b out of adj(a) in A. (+)
        Q <- c(Q,b)
        P <- c(P,a)
        PSEP <- c(PSEP,b)
      }
      if (verbose) cat(if(lres)"T" else ".", "")
    }
  }
  sort(setdiff(PSEP,x))
} ## {qreach}

## only called in fci() [by default:  doPdsep=TRUE]
pdsep <- function (skel, suffStat, indepTest, p, sepset, alpha, pMax, m.max = Inf,
                   pdsep.max = Inf, NAdelete = TRUE, unfVect = NULL,
                   biCC = FALSE, verbose = FALSE) ## FIXME: verbose : 2 --> qreach(verbose)
{
  ## Purpose: Compute Possible-D-SEP for each node, perform the condittional
  ##          independent tests and adapt graph accordingly
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - skel: Graph object returned by function skeleton
  ## - suffStat, indepTest: infofor the independence tests
  ## - p: number of nodes in the graph
  ## - sepset: Sepset that was used for finding the skeleton
  ## - alpha: niveau for the tests
  ## - pMax: Maximal p-values during estimation of skeleton
  ## - m.max: maximal size of the conditioning sets
  ## - pdsep.max: maximaum size of conditioning set for Possible-D-SEP
  ## - unfVect: vector containing the unfaithful triples, used for the
  ##   conservative orientation of the v-structures
  ## - biCC: if the biconnected components have to be used
  ## ----------------------------------------------------------------------
  ## Value:
  ## - G: Updated boolean adjacency matrix
  ## - sepset: Updated sepsets
  ## - pMax: Updated pMax
  ## - allPdsep: Possible d-sep for each node [list]
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  9 Dec 2009
  ## Modification: Diego Colombo; Martin Maechler

  G <- (as(skel, "matrix") != 0)
  n.edgetests <- rep(0, 1000)
  ord <- 0L
  allPdsep.tmp <- vector("list", p)
  if(biCC)
    conn.comp <- lapply(biConnComp(skel), as.numeric)
  if (any(G)) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    if (verbose) cat("\nCompute collider:\n")
    for (i in seq_len(nrow(ind))) {
      x <- ind[i, 1]
      y <- ind[i, 2]
      allZ <- setdiff(which(amat[y, ] != 0), x)
      for (z in allZ) {
        if (amat[x, z] == 0 &&
            !(y %in% sepset[[x]][[z]] ||
              y %in% sepset[[z]][[x]])) {

          if (length(unfVect) == 0) { ## normal version -------------------
            amat[x, y] <- amat[z, y] <- 2
            if (verbose) cat("\n",x,"*->", y, "<-*", z, "\n")
          }
          else { ## conservative version : check if x-y-z is faithful
            if (!any(unfVect == triple2numb(p,x,y,z), na.rm = TRUE) &&
                !any(unfVect == triple2numb(p,z,y,x), na.rm = TRUE)) {
              amat[x, y] <- amat[z, y] <- 2
              if (verbose)
                cat("\n",x,"*->", y, "<-*", z, "\n")
            }
          }
        }
      } ## for( z )
    } ## for( i  )
    allPdsep <- lapply(1:p, qreach, amat = amat)# verbose = (verbose >= 2)
    allPdsep.tmp <- vector("list", p)
    for(x in 1:p) {
      if(verbose) cat("\nPossible D-Sep of", x, "is:", allPdsep[[x]], "\n")
      if (any(an0 <- amat[x, ] != 0)) {
        tf1 <- setdiff(allPdsep[[x]], x)
        adj.x <- which(an0)
        for (y in adj.x) {
          if(verbose) cat(sprintf("\ny = %3d\n.........\n", y))
          tf <- setdiff(tf1, y)
          diff.set <- setdiff(tf, adj.x)
          ## bi-connected components
          if (biCC) {
            for(cci in conn.comp) {
              if (x %in% cci && y %in% cci)
                break ## found it
            }
            bi.conn.comp <- setdiff(cci, c(x,y))
            tf <- intersect(tf, bi.conn.comp)
            if (verbose) {
              cat("There is an edge between",x,"and",y,"\n")
              cat("Possible D-Sep of", x,
                  "intersected with the biconnected component of",x,"and",y,
                  "is:", tf, "\n")
            }
          } ## if(biCC)
          allPdsep.tmp[[x]] <- c(tf,y) ## you must add y to the set
          ## for the large scale simulations, we need to stop the algorithm if
          ## it takes to much time, i.e. sepset>25
          if (length(tf) > pdsep.max) {
            if(verbose)
              cat("Size of Possible-D-SEP bigger than",pdsep.max,
                  ". Break the search for the edge between", x,"and",y,"\n")
          } else if (length(diff.set) > 0) {
            done <- FALSE
            ord <- 0L
            while (!done && ord < min(length(tf), m.max)) {
              ord <- ord + 1L
              if(verbose) cat("ord = ", ord, "\n")
              if (ord == 1) {
                for (S in diff.set) {
                  pval <- indepTest(x, y, S, suffStat)
                  n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                  if (is.na(pval))
                    pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                  if (pval > pMax[x, y])
                    pMax[x, y] <- pval
                  if (pval >= alpha) {
                    amat[x, y] <- amat[y, x] <- 0
                    sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                    done <- TRUE
                    if (verbose)
                      cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                    break
                  }
                }
              }
              else { ## ord > 1
                tmp.combn <- combn(tf, ord) ## has  choose( |tf|, ord ) columns
                if (ord <= length(adj.x)) {
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    if (!all(S %in% adj.x)) {
                      n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                      pval <- indepTest(x, y, S, suffStat)
                      if (is.na(pval))
                        pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                      if(pMax[x, y] < pval)
                        pMax[x, y] <- pval
                      if (pval >= alpha) {
                        amat[x, y] <- amat[y, x] <- 0
                        sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                        done <- TRUE
                        if (verbose)
                          cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                        break
                      }
                    }
                  } ## for(k ..)
                }
                else { ## ord > |adj.x| :
                  ## check all combinations; no combination has been tested before
                  for (k in seq_len(ncol(tmp.combn))) {
                    S <- tmp.combn[, k]
                    n.edgetests[ord + 1] <- n.edgetests[ord + 1] + 1
                    pval <- indepTest(x, y, S, suffStat)
                    if (is.na(pval))
                      pval <- as.numeric(NAdelete) ## = if(NAdelete) 1 else 0
                    if(pMax[x, y] < pval)
                      pMax[x, y] <- pval
                    if (pval >= alpha) {
                      amat[x, y] <- amat[y, x] <- 0
                      sepset[[x]][[y]] <- sepset[[y]][[x]] <- S
                      done <- TRUE
                      if (verbose)
                        cat("x=", x, " y=", y, " S=", S, ": pval =", pval, "\n")
                      break
                    }
                  } ## for(k ..)
                } ## else: { ord > |adj.x| }
              } ## else

            } ## while(!done ..)
          }

        } ## for(y ..)

      } ## if(any( . ))

    } ## for(x ..)
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
    G[amat == 2] <- TRUE

  } ## if(any(G))

  list(G = G, sepset = sepset, pMax = pMax, allPdsep = allPdsep.tmp,
       max.ord = ord, n.edgetests = n.edgetests[1:(ord + 1)])
} ## {pdsep}


udag2pag <- function(pag, sepset, rules = rep(TRUE,10), unfVect = NULL, verbose = FALSE, orientCollider = TRUE)
{
  ## Purpose: Transform the Skeleton of a pcAlgo-object to a PAG using
  ## the rules of Zhang. The output is an adjacency matrix.
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - pag: adjacency matrix of size pxp
  ## - sepset: list of all separation sets
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - unfVect: Vector with ambiguous triples (coded as number using triple2numb)
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 6 Mar 2009; cleanup: Martin Maechler, 2010
  ## update: Diego Colombo, 2012

  ## Notation:
  ## ----------------------------------------------------------------------
  ## 0: no edge
  ## 1: -o
  ## 2: -> (arrowhead)
  ## 3: - (tail)
  ## a=alpha
  ## b=beta
  ## c=gamma
  ## d=theta

  stopifnot(is.logical(rules), length(rules) == 10)

  if (any(pag != 0)) {
    p <- as.numeric(dim(pag)[1])

    ## orient collider
    if (orientCollider) {
      ind <- which(pag == 1, arr.ind = TRUE)
      for (i in seq_len(nrow(ind))) {
        x <- ind[i, 1]
        y <- ind[i, 2]
        allZ <- setdiff(which(pag[y, ] != 0), x)
        for (z in allZ) {
          if (pag[x, z] == 0 && !((y %in% sepset[[x]][[z]]) ||
                   (y %in% sepset[[z]][[x]]))) {
            if (length(unfVect) == 0) {
              if (verbose) {
                cat("\n", x, "*->", y, "<-*", z, "\n")
                cat("Sxz=", sepset[[z]][[x]], "and",
                    "Szx=", sepset[[x]][[z]], "\n")
              }
              pag[x, y] <- pag[z, y] <- 2
            }
            else {
              if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
                if (verbose) {
                  cat("\n", x, "*->", y, "<-*", z, "\n")
                  cat("Sxz=", sepset[[z]][[x]], "and",
                      "Szx=", sepset[[x]][[z]], "\n")
                }
                pag[x, y] <- pag[z, y] <- 2
              }
            }
          }
        }
      }
    } ## end: Orient collider

    old_pag1 <- matrix(0, p, p)
    while (any(old_pag1 != pag)) {
      old_pag1 <- pag
      ##-- R1 ------------------------------------------------------------------
      if (rules[1]) {
        ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          b <- ind[i, 2]
          indC <- which((pag[b, ] != 0 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          indC <- setdiff(indC, a)
          if (length(indC) > 0) {
            if (length(unfVect) == 0) {
              pag[b, indC] <- 2
              pag[indC, b] <- 3
              if (verbose)
                cat("\nRule 1",
                    "\nOrient:", a, "*->", b, "o-*", indC,
                    "as:", b, "->", indC, "\n")
            }
            else {
              for (c in indC) {
                if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                    !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                  if (verbose)
                    cat("\nRule 1",
                        "\nConservatively orient:", a, "*->", b, "o-*",
                        c, "as:", b, "->", c, "\n")
                }
              } ## for( c )
            }
          }
        } ## for( i )
      }
      ##-- R2 ------------------------------------------------------------------
      if (rules[2]) {
        ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          indB <- which((pag[a, ] == 2 & pag[, a] == 3 & pag[c, ] != 0 & pag[, c] == 2) | (pag[a, ] == 2 & pag[, a] != 0 & pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) > 0) {
            pag[a, c] <- 2
            if (verbose) {
              cat("\nRule 2","\n")
              cat("Orient:", a, "->", indB, "*->", c, "or", a, "*->", indB, "->", c, "with", a, "*-o", c, "as:", a, "*->", c, "\n")
            }
          }
        }
      }
      ##-- R3 ------------------------------------------------------------------
      if (rules[3]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          d <- ind[i, 2]
          indAC <- which((pag[b, ] != 0 & pag[, b] == 2) & (pag[, d] == 1 & pag[d, ] != 0))
          if (length(indAC) >= 2) {
            if (length(unfVect) == 0) {
              counter <- 0
              while ((counter < (length(indAC) - 1)) && (pag[d, b] != 2)) {
                counter <- counter + 1
                ii <- counter
                while (ii < length(indAC) && pag[d, b] != 2) {
                  ii <- ii + 1
                  if (pag[indAC[counter], indAC[ii]] == 0 && pag[indAC[ii], indAC[counter]] == 0) {
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Orient:", d, "*->", b, "\n")
                    }
                    pag[d, b] <- 2
                  }
                }
              }
            }
            else {
              comb.indAC <- combn(indAC, 2)
              for (j in 1:dim(comb.indAC)[2]) {
                a <- comb.indAC[1, j]
                c <- comb.indAC[2, j]
                if (pag[a, c] == 0 && pag[c, a] == 0 && c != a) {
                  if (!any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                    pag[d, b] <- 2
                    if (verbose) {
                      cat("\nRule 3","\n")
                      cat("Conservatively orient:", d, "*->", b, "\n")
                    }
                  }
                }
              }
            }
          }
        }
      }
      ##-- R4 ------------------------------------------------------------------
      if (rules[4]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all a s.t. a -> c and a <-* b
          indA <- which((pag[b, ] == 2 & pag[, b] != 0) &
                        (pag[c, ] == 3 & pag[, c] == 2))
          ## chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while (length(indA) > 0 && pag[c,b] == 1) {
            a <- indA[1]
            indA <- indA[-1]
            ## path is the initial triangle
            ## abc <- c(a, b, c)
            ## Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
### MM: FIXME?? Isn't  Done  set to TRUE in *any* case inside the following
### while(.), the very first time already ??????????
            while (!Done && pag[a,b] != 0 && pag[a,c] != 0 && pag[b,c] != 0) {
              ## find a minimal discriminating path for a,b,c
              md.path <- minDiscrPath(pag, a,b,c, verbose = verbose)
              ## if the path doesn't exists, we are done with this triangle
              if ((N.md <- length(md.path)) == 1) {
                Done <- TRUE
              }
              else {
                ## a path exists
                ## if b is in sepset
                if ((b %in% sepset[[md.path[1]]][[md.path[N.md]]]) ||
                    (b %in% sepset[[md.path[N.md]]][[md.path[1]]])) {
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is in Sepset of",
                        c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                  pag[b, c] <- 2
                  pag[c, b] <- 3
                }
                else {
                  ## if b is not in sepset
                  if (verbose)
                    cat("\nRule 4",
                        "\nThere is a discriminating path between",
                        md.path[1], "and", c, "for", b, ",and", b, "is not in Sepset of",
                        c, "and", md.path[1], ". Orient:", a, "<->", b, "<->",
                        c, "\n")
                  pag[a, b] <- pag[b, c] <- pag[c, b] <- 2
                }
                Done <- TRUE
              }
            }
          }
        }
      }
      ##-- R5 ------------------------------------------------------------------
      if (rules[5]) {
        ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE) ## a o-o b
        while (length(ind) > 0) {
          a <- ind[1, 1]
          b <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all c s.t. a o-o c and c is not connected to b
          indC <- which((pag[a, ] == 1 & pag[, a] == 1) & (pag[b, ] == 0 & pag[, b] == 0))
          ## delete b since it is surely in indC
          indC <- setdiff(indC, b)
          ## find all d s.t. b o-o d and d is not connected to a
          indD <- which((pag[b, ] == 1 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
          ## delete a since it is surely in indD
          indD <- setdiff(indD, a)
          if (length(indC) > 0 && length(indD) > 0) {
            counterC <- 0
            while ((counterC < length(indC)) && pag[a, b] == 1) {
              counterC <- counterC + 1
              c <- indC[counterC]
              counterD <- 0
              while ((counterD < length(indD)) && pag[a, b] == 1) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if (pag[c, d] == 1 && pag[d, c] == 1) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[a, b] <- pag[b, a] <- 3
                    pag[a, c] <- pag[c, a] <- 3
                    pag[c, d] <- pag[d, c] <- 3
                    pag[d, b] <- pag[b, d] <- 3
                    if (verbose)
                      cat("\nRule 5",
                          "\nThere exists an uncovered circle path between", a, "and", b,
                          ". Orient:", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                  }
                  else { ## conservative: check that every triple on the circle is faithful
                    path2check <- c(a,c,d,b)
                    if (faith.check(path2check, unfVect, p)) {
                      pag[a, b] <- pag[b, a] <- 3
                      pag[a, c] <- pag[c, a] <- 3
                      pag[c, d] <- pag[d, c] <- 3
                      pag[d, b] <- pag[b, d] <- 3
                      if (verbose)
                        cat("\nRule 5",
                            "\nThere exists a faithful uncovered circle path between",
                            a, "and", b, ". Conservatively orient:",
                            a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                    }
                  }
                }
                ## search with a breitensuche a minimal uncovered circle path
                else {
                  ## Find a minimal uncovered circle path for these a,b,c, and d.
                  ## This path has already been checked to be uncovered and
                  ## to be faithful for the conservative case
		  ucp <- minUncovCircPath(p, pag = pag, path = c(a,c,d,b),
					  unfVect = unfVect, verbose = verbose)
                  ## there is a path ---> orient
                  if (length(ucp) > 1) {
                    ## orient every edge on the path as --
                    n <- length(ucp)
                    pag[ucp[1], ucp[n]] <- pag[ucp[n], ucp[1]] <- 3 ## a--b
                    for (j in 1:(length(ucp)-1)) ## each edge on the path --
                      pag[ucp[j], ucp[j + 1]] <- pag[ucp[j + 1], ucp[j]] <- 3
                  }
                }
              }
            }
          }
        }
      }
      ##-- R6 ------------------------------------------------------------------
      if (rules[6]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)## b o-* c
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          if (any(pag[b, ] == 3 & pag[, b] == 3)) {
            pag[c, b] <- 3
            if (verbose)
              cat("\nRule 6",
                  "\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
          }
        }
      }
      ##-- R7 ------------------------------------------------------------------
      if (rules[7]) {
        ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          b <- ind[i, 1]
          c <- ind[i, 2]
          indA <- which((pag[b, ] == 3 & pag[, b] == 1) & (pag[c, ] == 0 & pag[, c] == 0))
          indA <- setdiff(indA, c)
          if (length(indA) > 0) {
            if (length(unfVect) == 0) {
              pag[c, b] <- 3
              if (verbose)
                cat("\nRule 7",
                    "\nOrient:", indA, "-o", b, "o-*",
                    c, "as", b, "-*", c, "\n")
            }
            else for (a in indA)
              if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) &&
                  !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                pag[c, b] <- 3
                if (verbose)
                  cat("\nRule 7",
                      "\nConservatively orient:", a, "-o", b, "o-*",
                      c, "as", b, "-*", c, "\n")
              }
          }
        }
      }
      ##-- R8 ------------------------------------------------------------------
      if (rules[8]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
          a <- ind[i, 1]
          c <- ind[i, 2]
          indB <- which(pag[, a] == 3 & (pag[a, ] == 2 | pag[a, ] == 1) &
                        pag[c, ] == 3 & pag[, c] == 2)
          if (length(indB) > 0) {
            pag[c, a] <- 3
            if (verbose)
              cat("\nRule 8",
                  "\nOrient:", a, "->", indB, "->", c,
                  "or", a, "-o", indB, "->", c, "with", a,
                  "o->", c, "as", a, "->", c, "\n")
          }
        }
      }
      ##-- R9 ------------------------------------------------------------------
      if (rules[9]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. a (o-)--(o>) b and b and c are not connected
          indB <- which((pag[a, ] == 2 | pag[a, ] == 1) &
                        (pag[, a] == 1 | pag[, a] == 3) &
                        (pag[c, ] == 0 & pag[, c] == 0))
          ## delete c from indB since it is surely inside
          indB <- setdiff(indB, c)
          ## chose one b s.t. the initial structure exists and the edge hasn't been oriented yet
          while ((length(indB) > 0) && (pag[c,a] == 1)) {
            b <- indB[1]
            indB <- indB[-1]
	    ## find a minimal uncovered pd path from initial (a,b,c) :
	    upd <- minUncovPdPath(p, pag, a,b,c,
				  unfVect = unfVect, verbose = verbose)
	    ## there is a path ---> orient it
	    if (length(upd) > 1) {
	      pag[c, a] <- 3
	      if (verbose)
		cat("\nRule 9",
		    "\nThere exists an uncovered potentially directed path between", a, "and", c,
		    ". Orient:", a, " ->",c, "\n")
	    }
	  }
	}
      }
      ##-- R10 ------------------------------------------------------------------
      if (rules[10]) {
        ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. b --> c
          indB <- which((pag[c, ] == 3 & pag[, c] == 2))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (pag[c, a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (pag[c, a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if ((pag[a, b] == 1 || pag[a, b] == 2) &&
                    (pag[b, a] == 1 || pag[b, a] == 3) &&
                    (pag[a, d] == 1 || pag[a, d] == 2) &&
                    (pag[d, a] == 1 || pag[d, a] == 3) && pag[d, b] == 0 && pag[b, d] == 0) {
                  if (length(unfVect) == 0) { ## normal version
                    pag[c, a] <- 3
                    if (verbose)
                      cat("\nRule 10 [easy]",
                          "\nOrient:", a, "->", c, "\n")
                  }
                  else ## conservative version: check faithfulness of b-a-d
                    if (!any(unfVect == triple2numb(p,b,a,d), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,d,a,b), na.rm = TRUE)) {
                      pag[c, a] <- 3
                      if (verbose)
                        cat("\nRule 10 [easy]",
                            "\nConservatively orient:", a, "->", c, "\n")
                    }
                }
                ## search with a breitensuche two minimal uncovered circle paths
                else {
                  ## find all x s.t. a (o-)--(o>) x
                  indX <- which((pag[a, ] == 1 | pag[a, ] == 2) &
                                (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    counterX1 <- 0
                    while (counterX1 < length(indX) && pag[c, a] == 1) {
                      counterX1 <- counterX1 + 1
                      first.pos <- indA[counterX1]
                      indX2 <- setdiff(indX, first.pos)
                      counterX2 <- 0
                      while (counterX2 < length(indX2) && pag[c, a] == 1) {
                        counterX2 <- counterX2 + 1
                        sec.pos <- indX2[counterX2]
                        t1 <- minUncovPdPath(p, pag, a, first.pos, b,
                                             unfVect = unfVect, verbose = verbose)
                        if (length(t1) > 1) { # otherwise, can skip next minUnc..()
                          t2 <- minUncovPdPath(p, pag, a, sec.pos, d,
                                               unfVect = unfVect, verbose = verbose)
                          if (length(t2) > 1 &&
                              first.pos != sec.pos && pag[first.pos, sec.pos] == 0) {
                            ## we found 2 uncovered pd paths
                            if (length(unfVect) == 0) { ## normal version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10", "\nOrient:", a, "->", c, "\n")
                            }
                            else if(!any(unfVect == triple2numb(p,first.pos, a, sec.pos), na.rm = TRUE) &&
                                    !any(unfVect == triple2numb(p,sec.pos, a, first.pos), na.rm = TRUE)) {
                              ## conservative version
                              pag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10",
                                    "\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                      } #  # while ( counterX2 .. )
                    }
                  }
                } # else
              } # while ( counterD .. )
            } # while ( counterB .. )
          } # if (length(indB) .)
        }
      } ## if (rules[10] ..)
    }
  }
  pag
} ## udag2pag()



################################################################################
## RFCI
################################################################################

rfci <- function(suffStat, indepTest, alpha, labels, p,
                 skel.method = c("stable", "original", "stable.fast"),
                 fixedGaps = NULL, fixedEdges = NULL,
                 NAdelete = TRUE, m.max = Inf, rules = rep(TRUE, 10),
                 conservative = FALSE, maj.rule = FALSE,
                 verbose = FALSE)
{
  ## Purpose: Perform RFCI-Algorithm, i.e., estimate PAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - p: number of nodes in the graph
  ## - alpha: Significance level of individual partial correlation tests
  ## - verbose: 0 - no output, 1 - detailed output
  ## - fixedGaps: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical); gaps fixed here are not changed
  ## - fixedEdges: Edges marked here are not changed (logical)
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - conservative: TRUE or FALSE defining if the v-structures after
  ##                 the skeleton have to be oriented conservatively or nor
  ## - maj.rule: TRUE or FALSE variable containing if the majority rule is
  ##             used instead of the normal conservative
  ## - labels: names of the variables or nodes

  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, 2011; modifications: Martin Maechler

  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  if (conservative && maj.rule)
      stop("Can only choose one of conservative or majority rule RFCI")
  if (verbose) cat("Compute Skeleton\n================\n")

  skel <- skeleton(suffStat, indepTest, alpha, labels = labels, method = skel.method,
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, verbose = verbose)
  sk.A <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  ## the list of all ordered unshielded triples (the graph g does not change it is just a search!)
  u.t <- find.unsh.triple(sk.A, check = FALSE)
  ## check and orient v-structures recursively
  r.v. <- rfci.vStruc(suffStat, indepTest, alpha, sepset, sk.A,
                      unshTripl = u.t$unshTripl, unshVect = u.t$unshVect,
                      conservative = (conservative || maj.rule),
                      version.unf = c(1,1), maj.rule = maj.rule, verbose = verbose)
  A <- r.v.$amat
  sepset <- r.v.$sepset

  ## orient as many edge marks as possible
  if (verbose)
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules), "\n")

  res <- udag2apag(A, suffStat, indepTest, alpha, sepset,
		   rules = rules, unfVect = r.v.$unfTripl, verbose = verbose)
  Amat <- res$graph
  colnames(Amat) <- rownames(Amat) <- labels
  new("fciAlgo", amat = Amat, call = cl, n = integer(0),
        max.ord = as.integer(skel@max.ord), max.ordPDSEP = 0L,
        n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
        sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
} ## {rfci}

find.unsh.triple <- function(g, check = TRUE)
{
  ## Purpose: find the ordered (<x,y,z> with x<z) list of all the unshielded
  ##          triples in the graph
  ## ----------------------------------------------------------------------
  ## Arguments: g: adjacency matrix
  ## ----------------------------------------------------------------------
  ## Values: - unshTripl: matrix with 3 rows containing in each column
  ##                      an unshielded triple
  ##         - unshVect: containing the unique number for each column
  ##                     in unshTripl
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010; clean- and speed-up: Martin Maechler

  if(check) stopifnot( all(g == t(g)) )
  m <- 0L
  unshTripl <- matrix(integer(), 3, m)
  if (any(g != 0)) { # if( at least one edge) -- find all unshielded triples in g
    p <- nrow(g)
    indS <- which(g == 1, arr.ind = TRUE) ## x-y
    for (i in seq_len(nrow(indS))) {
      xy <- indS[i,]
      x <- xy[1]
      y <- xy[2]
      allZ <- setdiff(which(g[y, ] == 1), x) ## x-y-z
      for (z in allZ) {
        if (g[x,z] == 0 && g[z,x] == 0) {
          ## save the matrix
          unshTripl <- cbind(unshTripl, c(xy, z))
        }
      }
    }
    ## delete duplicates in the matrix
    if ((m <- ncol(unshTripl)) > 0) {
      deleteDupl <- logical(m)# all FALSE
      for (i in seq_len(m)) ## FIXME -- make faster!
        if (unshTripl[1,i] > unshTripl[3,i])
          deleteDupl[i] <- TRUE
      if(any(deleteDupl))
        m <- ncol(unshTripl <- unshTripl[,!deleteDupl, drop = FALSE])
    }
  }
  ## define the vector with the unique number for each triple
  unshVect <- vapply(seq_len(m), function(k)
                     triple2numb(p, unshTripl[1,k], unshTripl[2,k], unshTripl[3,k]),
                     numeric(1))
  list(unshTripl = unshTripl, unshVect = unshVect)
} ## {find.unsh.triple}

##' called only from  rfci()  and  checkEdges()
rfci.vStruc <- function(suffStat, indepTest, alpha, sepset, g.amat,
                        unshTripl, unshVect, conservative = FALSE,
                        version.unf = c(2,1), maj.rule = FALSE, verbose = FALSE)
{
  ## Purpose: check the unshielded triples in unshTripl as v-structures
  ##          recursively check if new unshielded triples have been found
  ##          then save and orient the final ones in finalList and finalVect
  ## ----------------------------------------------------------------------
  ## Arguments: - suffStat, indepTest, p,alpha: arguments from the algorithm
  ##            - sepset, g.amat: output of the function skeleton
  ##            - unshTripl, unshVect: list/numbers of the unshielded triples
  ##                                   in graph
  ##            - conservative: TRUE or FALSE
  ##            - version.unf=c(2,1): conservative version
  ##            - maj.rule: TRUE or FALSE variable containing if the majority
  ##             rule is used instead of the normal conservative
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset, graph (g.amat) and list of unfaithful v-structures
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010, 14:15

  ## check every column (unshielded triple) of unshTripl
  ## NOTE: what's done is done!!!! I don't look back in the matrix

  stopifnot(is.matrix(unshTripl), nrow(unshTripl) == 3)
  nT <- ncol(unshTripl)
  A <- g.amat
  ## list of all unfaithful v-structures
  unfTripl <- NULL
  if (nT) {
    if (verbose)
      cat("\nCheck unshielded triples for conditional dependence
================================================\n")
    ## check all the columns in unshTripl as v-structure or not
    checktmp <- dep.triple(suffStat, indepTest, alpha, sepset = sepset, apag = g.amat,
                           unshTripl = unshTripl, unshVect = unshVect,
                           ## per default every triple is defined as a v-structure :
                           trueVstruct = rep(TRUE, nT), verbose = verbose)
    ## save the updated objects
    A          <- checktmp$apag
    sepset     <- checktmp$sepset
    ## note that no column is deleted from the matrix, we can add new triples instead
    unshTripl  <- checktmp$triple
    unshVect   <- checktmp$vect
    trueVstruct <- checktmp$trueVstruct

    if (verbose) cat(
      if (conservative)
      "\nOrient the v-structures conservatively
=====================================\n"
      else "\nOrient the v-structures\n=======================\n")

    ## if there is at least one triple with the desired properties
    if (any(trueVstruct)) {
      for (i in 1:dim(unshTripl)[2]) {
        if (trueVstruct[i]) {
          x <- unshTripl[1, i]
          y <- unshTripl[2, i]
          z <- unshTripl[3, i]
          if (!conservative) {
            if (A[x, z] == 0 && A[x,y] != 0 && A[z,y] != 0 &&
                !((y %in% sepset[[x]][[z]]) || (y %in% sepset[[z]][[x]]))) {
              ## this is to avoid the problem of:
              ## assume that <a,b,c> was saved in finalList because a "true" unshielded triple
              ## but after in the matrix appears <b,c,d> and the edge b-c is deleted
              ## of course the triple <a,b,c> stays in the finalList but we cannot orient the edge
              ## b-c because it doesn'e exist anymore

              if (verbose)
                cat("\n", x, "*->", y, "<-*", z,
                  "\nSxz=", sepset[[z]][[x]], "and",
                    "Szx=", sepset[[x]][[z]], "\n")
              A[c(x,z), y] <- 2
            }
          } else { ## conservative version
            ## check if x-y-z is faithful
            ## find neighbours of x and z
            nbrsX <- which(A[,x] != 0) ## G symm; x no nbr of z
            nbrsZ <- which(A[,z] != 0)
            if (verbose)
              cat("\nTriple:",x,y,z,"and sepset by skelet:",
                  unique(sepset[[x]][[z]],sepset[[z]][[x]]),"\n")
            r.abc <- checkTriple(x, y, z, nbrsX, nbrsZ,
                                 sepset[[x]][[z]], sepset[[z]][[x]],
                                 suffStat = suffStat, indepTest = indepTest, alpha = alpha,
                                 version.unf = version.unf, maj.rule = maj.rule,
                                 verbose = verbose)
            p <- nrow(A)
            ## 1: in NO set; 2: in ALL sets; 3: in SOME but not all
            if (r.abc$decision == 3 || ## <- take action if case "3", or
                ## can happen the case in Tetrad, so we must save the triple as unfaithful
                ## a and c independent given S but not given subsets of the adj(x) or adj(z)
                (version.unf[1] == 2 && r.abc$version == 2)) {

              ## record unfaithful triple
              unfTripl <- c(unfTripl, triple2numb(p, x, y, z))
              if (verbose >= 2) cat("new unfTriple:", x, y, z, "\n")
            }
            sepset[[x]][[z]] <- r.abc$SepsetA
            sepset[[z]][[x]] <- r.abc$SepsetC
            if (A[x,z] == 0 && A[x,y] != 0 && A[z,y] != 0 &&
                !((y %in% sepset[[x]][[z]]) ||
                  (y %in% sepset[[z]][[x]])) &&
                !any(unfTripl == triple2numb(p, x, y, z), na.rm = TRUE) &&
                !any(unfTripl == triple2numb(p, z, y, x), na.rm = TRUE)) {
              if (verbose)
                cat("\nOrient:", x, "*->", y, "<-*", z,
                  "\nSxz=", sepset[[z]][[x]], "and",
                    "Szx=", sepset[[x]][[z]], "\n")
              A[c(x,z), y] <- 2
            }
          } ## else : conservative
        }
      } ## for (i ...)
    }
  }
  list(sepset = sepset, amat = A, unfTripl = unfTripl)
} ## {rfci.vStruc}


## only called from rfci.vStruc() - with trueVstruct a vector of TRUE
dep.triple <- function(suffStat, indepTest, alpha, sepset, apag,
                       unshTripl, unshVect, trueVstruct, verbose = FALSE)
{
  ## Purpose: test the two edges of any unshielded triple in unshTripl (column)
  ##          for dependence given sepset. If independent find the minimal
  ##          sepset, delete the edge, define the triple as FALSE in
  ##          trueVstruct and search in the graph for new triple or destroyed
  ##          ones. Otherwise do nothing.
  ## ----------------------------------------------------------------------
  ## Arguments: - suffStat, indepTest, alpha, sepset, apag: skeleton parameters / result
  ##            - unshTripl: matrix containing the unshielded triples (columns)
  ##            - unshVect: triple2numb of unshTripl
  ##            - trueVstruct: vector containing T/F for the v-structures
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset, apag, unshTripl, unshVect, trueVstruct
  ## Author: Diego Colombo, Date: 16 Aug 2010, 15:07

  p <- nrow(apag) ## == length(sepset)
  for(k in seq_len(ncol(unshTripl))) {
    ## Note that trueVstruct[.] is changed *inside* !
    if (trueVstruct[k]) { ## triple under inspection x-y-z
      x <- unshTripl[1, k]
      y <- unshTripl[2, k]
      z <- unshTripl[3, k]
      SepSet <- setdiff(unique(c(sepset[[x]][[z]],
				 sepset[[z]][[x]])), y)
      nSep <- length(SepSet)
      if (verbose)
        cat("\nTriple:", x, y, z, "and sepSet (of size", nSep,") ", SepSet,"\n")
      if (nSep != 0) {
        x. <- min(x,y)
        y. <- max(x,y)
        y_ <- min(y,z)
        z_ <- max(y,z)

        ## Check x-y ---------------------------------------------------------------

        del1 <- FALSE
        ## check first x and y given the whole sepset(x,z)
        if (indepTest(x, y, SepSet, suffStat) >= alpha) {
          ## afterwards we have to define this triple as FALSE in trueVstruct
          del1 <- TRUE
          ## x and y are independent, then find the minimal sepset
          done <- FALSE
          ord <- 0L
          while (!done && ord < nSep) {
            ord <- ord + 1L
            ## all combinations of SepSet of size ord
            S.j <- if (ord == 1 && nSep == 1) matrix(SepSet,1,1) else combn(SepSet, ord)
            for (i in seq_len(ncol(S.j))) {
              pval <- indepTest(x, y, S.j[,i], suffStat)
              if (verbose) cat("x=", x, " y=", y, "S=", S.j[,i], ": pval =", pval, "\n")
              if (pval >= alpha) {
                ## delete edge and save set in sepset
                apag[x, y] <- apag[y, x] <- 0
                sepset[[x]][[y]] <- sepset[[y]][[x]] <- S.j[,i]
                done <- TRUE
                break
              }
            }
          }

          ## case 1: before we had a triangle and now it is an unshielded triple x-m-y
          indM <- which((apag[x, ] == 1 & apag[, x] == 1) & (apag[y, ] == 1 & apag[, y] == 1))
          indM <- setdiff(indM, c(x,y,z)) ## just to be sure
          for (m in indM) {
            ## in the matrix the first column is always smaller than the third
            unshTripl <- cbind(unshTripl,        c(x., m, y.))
            unshVect <- c(unshVect, triple2numb(p, x., m, y.))
            ## per default this new triple is set as TRUE in trueVstruct
            trueVstruct <- c(trueVstruct, TRUE)
          }
          ## case 2: an existent unshielded triple has been destroyed
          ## case 2.a): we had q-x-y or y-x-q in the matrix but not anymore
          indQ <- which((apag[x, ] == 1 & apag[, x] == 1) & (apag[y, ] == 0 & apag[, y] == 0))
          indQ <- setdiff(indQ, c(x,y,z)) ## just to be sure
          for (q in indQ) {
            ## define the triple as FALSE in trueVstruct
            delTripl <- unshVect == (
              if (q < y) triple2numb(p, q, x, y) else triple2numb(p, y, x, q))
            ## if we haven't checked the triple yet, define it as FALSE
            if (any(delTripl))
              trueVstruct[which.max(delTripl)] <- FALSE
          }

          ## case 2.b): we had r-y-x or x-y-r in the matrix but not anymore
          indR <- which((apag[x, ] == 0 & apag[, x] == 0) & (apag[y, ] == 1 & apag[, y] == 1))
          indR <- setdiff(indR, c(x,y,z)) ## just to be sure
          ##                ^^^ (had typo here: 'indQ' !!)
          for (r in indR) {
            ## define the triple as FALSE in trueVstruct
            delTripl <- unshVect == (
              if (r < x) triple2numb(p, r, y, x) else triple2numb(p, x, y, r))
            if (any(delTripl))
              trueVstruct[which.max(delTripl)] <- FALSE
          }

        } ## if (pv1 .. indepTest(..) )

        ## Check z-y ---------------------------------------------------------------

        del2 <- FALSE
        ## check first z and y given the whole sepset(x,z)
        if (indepTest(z, y, SepSet, suffStat) >= alpha) {
          del2 <- TRUE ## will have to set this triple as FALSE in trueVstruct
          ## z and y are independent, then find the minimal sepset
          Done <- FALSE
          Ord <- 0L
          while (!Done && Ord < nSep) {
            Ord <- Ord + 1L
            ## all combinations of SepSet of size Ord
            S.j <- if (Ord == 1 && nSep == 1) matrix(SepSet,1,1) else combn(SepSet, Ord)
            for (i in seq_len(ncol(S.j))) {
              pval <- indepTest(z, y, S.j[,i], suffStat)
              if (verbose) cat("x=", z, " y=", y, " S=", S.j[,i], ": pval =", pval, "\n")
              if (pval >= alpha) {
                ## delete edge and save set in sepset
                apag[z, y] <- apag[y, z] <- 0
                sepset[[z]][[y]] <- sepset[[y]][[z]] <- S.j[,i]
                Done <- TRUE
                break
              }
            }
          } ## while( . )

          ## case 1: before we had a triangle and now it is an unshielded triple z-m-y
          indM <- which((apag[z, ] == 1 & apag[, z] == 1) & (apag[y, ] == 1 & apag[, y] == 1))
          indM <- setdiff(indM,c(x,y,z)) ## just to be sure
          for (m in indM) {
            ## in the matrix the first column is always smaller than the third
            unshTripl <- cbind(unshTripl,        c(y_, m, z_))
            unshVect <- c(unshVect, triple2numb(p, y_, m, z_))
            ## per default the triple is defined as TRUE
            trueVstruct <- c(trueVstruct,TRUE)
          }

          ## case 2: an existent unshielded triple has been destroyed
          ## case 2.a): we had q-z-y or y-z-q in the matrix but not anymore
          indQ <- which((apag[z, ] == 1 & apag[, z] == 1) & (apag[y, ] == 0 & apag[, y] == 0))
          indQ <- setdiff(indQ, c(x,y,z)) ## just to be sure
          for (q in indQ) {
            ## define the triple as FALSE in trueVstruct
            delTripl <- unshVect == (if (q < y) triple2numb(p, q, z, y) else triple2numb(p, y, z, q))
            ## if we haven't checked the triple yet, save it as FALSE
            if (any(delTripl))
              trueVstruct[which.max(delTripl)] <- FALSE
          }
          ## case 2.b): we had r-y-z or z-y-r in the matrix but not anymore
          indR <- which((apag[z, ] == 0 & apag[, z] == 0) & (apag[y, ] == 1 & apag[, y] == 1))
          indR <- setdiff(indR, c(x,y,z)) ## just to be sure
          ##                ^^^ (had typo here: 'indQ' !!)
          for (r in indR) {
            ## define the triple as FALSE in trueVstruct
            delTripl <- unshVect == (
              if (r < z) triple2numb(p, r, y, z) else triple2numb(p, z, y, r))
            if (any(delTripl))
              trueVstruct[which.max(delTripl)] <- FALSE
          }

        }
        ## if at least one edge has been deleted this is not a future v-structure
        if (any(del1, del2)) trueVstruct[k] <- FALSE

      } ## if (nSep != 0)
      ##
      ## ELSE: the sepset is empty
      ##   so surely <x,y,z> is an unshielded triple because they cannot be indep given the empty set
      ##   nothing changed in sepset and in graph
      ## else{
      ##   sepset <- sepset
      ##   apag <- apag
      ## }

      ## recursion on the next column of unshTripl
      ## rec.res <- dep.triple(suffStat, indepTest, p, alpha, sepset, apag, unshTripl, unshVect, trueVstruct, k, verbose=verbose)
      ## save the modified objects
      ## unshTripl <- rec.res$triple
      ## unshVect <- rec.res$vect
      ## sepset <- rec.res$sepset
      ## apag <- rec.res$graph
      ## trueVstruct <- rec.res$trueVstruct
    }
  }## for ( k )

  list(triple = unshTripl, vect = unshVect, sepset = sepset,
       apag = apag, trueVstruct = trueVstruct)
} ## {dep.triple}

##' Check the dependence of the edges for R4
##' only called from  udag2apag()
checkEdges <- function(suffStat, indepTest, alpha, apag, sepset, path,
                       unfVect = NULL, verbose = FALSE)
{
  ## Purpose: check if every edge on the path should exist in R4
  ## ----------------------------------------------------------------------
  ## Values: - updated sepset and apag
  ##         - deleted==FALSE no edge has been deleted on the path
  ##                  ==TRUE the discriminating path doesn't exist anymore
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 17 Aug 2010, 16:21

  stopifnot((n.path <- length(path)) >= 2)
  ## did we delete an edge?
  found <- FALSE
  ## conservative v-structures or not?
  conservative <- (length(unfVect) > 0)

  ## set that at the beginning there are no new unfaithful v-structures
  unfTripl <- NULL
  ## define the sepset
  SepSet.tot <- unique(c(sepset[[path[1]]][[path[n.path]]],
                         sepset[[path[n.path]]][[path[1]]]))
  if (length(SepSet.tot) != 0) {
    if (verbose)
      cat("\nCheck discriminating path:", path,
          "for dependence of any edge given sepset", SepSet.tot,"\n")
    p <- nrow(apag)
    ## check every edge on the path for independence given every possible subset of SepSet
    for (i in seq_len(n.path-1)) {
      x <- path[i]
      y <- path[i+1]
      SepSet <- setdiff(SepSet.tot, c(x,y))
      x. <- min(x,y)
      y. <- max(x,y)
      if (verbose >= 2)
        cat("Edge: ",x,"*-*",y, "; Sepset=", SepSet, "; |S|=", length(SepSet),"\n")
      if (length(SepSet) != 0) {
        j <- 0
        while (!found && j < length(SepSet)) {
          j <- j + 1
          ## all combinations of SepSet of size j
          S.j <- if (j == 1 && length(SepSet) == 1) matrix(SepSet,1,1) else combn(SepSet, j)
          ii <- 0
          while (!found && ii < ncol(S.j)) {
            ii <- ii + 1
            pval <- indepTest(x, y, S.j[,ii], suffStat)
            if (verbose)
              cat("x=", x, " y=", y, " S=", S.j[,ii], ": pval =", pval,"\n")
            if (pval >= alpha) {
              if (verbose) cat("Independence found: delete edge between",x,"and",y,"\n")
              found <- TRUE
              ## delete edge and save set in sepset
              apag[x, y] <- apag[y, x] <- 0
              sepset[[x]][[y]] <- sepset[[y]][[x]] <- S.j[,ii]
              ## before we had a triangle and now it is an unshielded triple x-m-y
              indM <- setdiff(which(apag[x, ] != 0 & apag[, x] != 0 &
                                    apag[y, ] != 0 & apag[, y] != 0),
                              c(x,y))## just to be sure
              ## create the list with all the new unshielded triples to be tested
              if ((nI <- length(indM)) > 0) {
                triplM <- matrix(integer(), 3, nI)
                newVect <- numeric(nI)
                for (jj in seq_len(nI)) {
                  m <- indM[jj]
                  triplM[,jj] <- c(x.,m,y.)
                  newVect[jj] <- triple2numb(p, x.,m,y.)
                }
                ## new unshielded triple to be tested
                r.v <- rfci.vStruc(suffStat, indepTest, alpha, sepset, apag,
                                   unshTripl = triplM, unshVect = newVect,
                                   conservative = conservative, verbose = verbose)
                ## save the modified graph g in apag
                apag <- r.v$amat
                ## save the new sepset
                sepset <- r.v$sepset
                ## save the new unfTripl, since we tested new v-structures and some can be unfaithful
                ## for the conservative version
                unfTripl <- r.v$unfTripl
              }
            }
          } ## while(!found && i < *)
        } ## while(!found && j < *)
      }
    } ## for(i ....)
  }
  ## if SepSet is the empty set do nothing because surely the vertices are dependent
  list(deleted = found, apag = apag, sepset = sepset, unfTripl = unfTripl)
}

##' called only from  rfci()
udag2apag <- function (apag, suffStat, indepTest, alpha, sepset,
                       rules = rep(TRUE, 10), unfVect = NULL,
                       verbose = FALSE)
### FIXME part of this is cut-n-paste from udag2pag()
{
  ## Purpose: use the 10 orientation rules to orient the skeleton. R4 about
  ##          discriminating paths also checks every edge for dependence
  ##          given all subsets of the separating set.
  ##          Note that the preliminaries v-structures are already oriented
  ##          in apag
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - apag: adjacency matrix outputs of the function rfci.vStruc
  ## - suffStat, indepTest, alpha: arguments required for the tests in R4
  ## - sepset: list of all separation sets
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - unfVect: Vector with ambiguous triples (coded as number using triple2numb)
  ## - verbose: 0 - no output, 1 - detailed output
  ## ----------------------------------------------------------------------
  ## Values: updated apag (oriented) and sepset
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 21 Oct 2010, 15:13

  if (any(apag != 0)) {
    p <- ncol(apag)
    old_apag1 <- matrix(0, nrow = p, ncol = p)
    while (any(old_apag1 != apag)) {
      old_apag1 <- apag
      ##-- R1 ------------------------------------------------------------------
      if (rules[1]) {
        ind <- which((apag == 2 & t(apag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            indC <- which(apag[b, ] != 0 & apag[, b] == 1 &
                          apag[a, ] == 0 & apag[, a] == 0)
            indC <- setdiff(indC, a)
            if (length(indC) > 0) {
              ## normal version
              if (length(unfVect) == 0) {
                apag[b, indC] <- 2
                apag[indC, b] <- 3
                if (verbose) {
                  cat("\nRule 1","\n")
                  cat("Orient:", a, "*->", b, "o-*", indC,
                      "as:", b, "->", indC, "\n")
                }
              }
              ## conservative
              else
                for (j in seq_along(indC)) {
                  c <- indC[j]
                  ## check that a-b-c faithful
                  if (!any(unfVect == triple2numb(p,a,b,c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p,c,b,a), na.rm = TRUE)) {
                    apag[b, c] <- 2
                    apag[c, b] <- 3
                    if (verbose) {
                      cat("\nRule 1","\n")
                      cat("Conservatively orient:", a, "*->", b, "o-*",
                          c, "as:", b, "->", c, "\n")
                    }
                  }
                }
            }
          }
      }
      ##-- R2 ------------------------------------------------------------------
      if (rules[2]) {
        ind <- which((apag == 1 & t(apag) != 0), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which(((apag[a, ] == 2 & apag[, a] == 3) &
                           (apag[c, ] != 0 & apag[, c] == 2)) |
                          ((apag[a, ] == 2 & apag[, a] != 0) &
                           (apag[c, ] == 3 & apag[, c] == 2)))
            if (length(indB) > 0) {
              apag[a, c] <- 2
              if (verbose) {
                cat("\nRule 2","\n")
                cat("Orient:", a, "->", indB, "*->",
                    c, "or", a, "*->", indB, "->", c, "with",
                    a, "*-o", c, "as:", a, "*->", c, "\n")
              }
            }
        }
      }
      ##-- R3 ------------------------------------------------------------------
      if (rules[3]) {
        ind <- which(apag != 0 & t(apag) == 1, arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            b <- ind[i, 1]
            d <- ind[i, 2]
            indAC <- which(apag[b, ] != 0 & apag[, b] == 2 &
                           apag[, d] == 1 & apag[d, ] != 0)
            if (length(indAC) >= 2) {
              if (length(unfVect) == 0) { ## normal version

                counter <- 0
                while ((counter < (length(indAC) - 1)) &&
                       (apag[d, b] != 2)) {
                  counter <- counter + 1
                  ii <- counter
                  while ((ii < length(indAC)) && (apag[d, b] != 2)) {
                    ii <- ii + 1
                    if (apag[indAC[counter], indAC[ii]] == 0 &&
                        apag[indAC[ii], indAC[counter]] == 0) {
                      apag[d, b] <- 2
                      if (verbose)
                        cat("\nRule 3","\nOrient:", d, "*->", b, "\n")
                    }
                  }
                }
              }
              else { ## conservative version
                comb.indAC <- combn(indAC,2)
                for (j in seq_len(ncol(comb.indAC))) {
                  a <- comb.indAC[1,j]
                  c <- comb.indAC[2,j]
                  if (apag[a,c] == 0 && apag[c,a] == 0 && c != a) {
                    ## check faithfulness a-d-c
                    if (!any(unfVect == triple2numb(p,a,d,c), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,c,d,a), na.rm = TRUE)) {
                      apag[d, b] <- 2
                      if (verbose)
                        cat("\nRule 3","\nConservatively orient:",  d, "*->", b, "\n")
                    }
                  }
                }
              }
            }
        }
      }
      ##-- R4 ------------------------------------------------------------------
      if (rules[4]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)## b o-* c
        while (length(ind) > 0) {
          b <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all a s.t. a -> c and a <-* b
          indA <- which((apag[b, ] == 2 & apag[, b] != 0) & (apag[c, ] == 3 & apag[, c] == 2))
          ## chose one a s.t. the initial triangle structure exists and the edge hasn't been oriented yet
          while (length(indA) > 0 && apag[c,b] == 1) {
            a <- indA[1]
            indA <- indA[-1]
            ## Done is TRUE if either we found a minimal path or no path exists for this triangle
            Done <- FALSE
            while (!Done && apag[a,b] != 0 && apag[a,c] != 0 && apag[b,c] != 0) {
              ## find a minimal disciminating path for a,b,c
              md.path <- minDiscrPath(apag, a,b,c, verbose = verbose)
              N.md <- length(md.path)
              ## if the path doesn't exists, we are done with this triangle
              if (N.md == 1) {
                Done <- TRUE
              }
              else {
                ## a path exists and needs to be checked and maybe oriented
                ## first check every single edge for independence
		chkE <- checkEdges(suffStat, indepTest, alpha = alpha,
				   apag = apag, sepset = sepset, path = md.path,
				   unfVect = unfVect, verbose = verbose)
                ## save updated graph, sepset,and UnfVect
                sepset <- chkE$sepset
                apag <- chkE$apag
                unfVect <- c(unfVect, chkE$unfTripl)
                ## no edge deleted ----> orient the edges
                if (!chkE$deleted) {
                  ## if b is in sepset
                  if (b %in% sepset[[md.path[1]]][[md.path[N.md]]] ||
                      b %in% sepset[[md.path[N.md]]][[md.path[1]]]) {
                    if (verbose) {
                      cat("\nRule 4","\n")
                      cat("There is a discriminating path between",
                          md.path[1], "and", c, "for", b, ",and", b, "is in Sepset of",
                          c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                    }
                    apag[b, c] <- 2
                    apag[c, b] <- 3
                  }
                  else {
                    ## if b is not in sepset
                    if (verbose) {
                      cat("\nRule 4","\n")
                      cat("There is a discriminating path between:",
                          md.path[1], "and", c, "for", b, ",and", b, "is not in Sepset of",
                          c, "and", md.path[1], ". Orient", a, "<->", b, "<->",
                          c, "\n")
                    }
                    apag[a, b] <- apag[b, c] <- apag[c, b] <- 2
                  }
                  Done <- TRUE
                }
              }
            } ## while(!Done && ...)
          }
        }
      }
      ##-- R5 ------------------------------------------------------------------
      if (rules[5]) {
        ind <- which((apag == 1 & t(apag) == 1), arr.ind = TRUE) ## a o-o b
        while (length(ind) > 0) {
          a <- ind[1, 1]
          b <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all c s.t. a o-o c and c is not connected to b
          indC <- which((apag[a, ] == 1 & apag[, a] == 1) & (apag[b, ] == 0 & apag[, b] == 0))
          ## delete b since it is surely in indC
          indC <- setdiff(indC, b)
          ## find all d s.t. b o-o d and d is not connected to a
          indD <- which((apag[b, ] == 1 & apag[, b] == 1) & (apag[a, ] == 0 & apag[, a] == 0))
          ## delete a since it is surely in indD
          indD <- setdiff(indD, a)
          if (length(indC) > 0 && length(indD) > 0) {
            counterC <- 0
            while ((counterC < length(indC)) && apag[a, b] == 1) {
              counterC <- counterC + 1
              c <- indC[counterC]
              counterD <- 0
              while ((counterD < length(indD)) && apag[a, b] == 1) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if (apag[c, d] == 1 && apag[d, c] == 1) {
                  ## normal version
                  if (length(unfVect) == 0) {
                    apag[a, b] <- apag[b, a] <- 3
                    apag[a, c] <- apag[c, a] <- 3
                    apag[c, d] <- apag[d, c] <- 3
                    apag[d, b] <- apag[b, d] <- 3
                    if (verbose) {
                      cat("\nRule 5","\n")
                      cat("There exists an uncovered circle path between",
                          a, "and", b, ". Orient", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                    }
                  }
                  ## conservative
                  else {
                    ## check that every triple on the circle is faithful
                    path2check <- c(a,c,d,b)
                    if (faith.check(path2check, unfVect, p)) {
                      apag[a, b] <- apag[b, a] <- 3
                      apag[a, c] <- apag[c, a] <- 3
                      apag[c, d] <- apag[d, c] <- 3
                      apag[d, b] <- apag[b, d] <- 3
                      if (verbose) {
                        cat("\nRule 5","\n")
                        cat("There exists a faithful uncovered circle path between", a, "and", b,
                            ". Conservatively orient:", a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                      }
                    }
                  }
                }
                ## search with a breitensuche a minimal uncovered circle path
                else {
                  ## Find a minimal uncovered circle path for these a,b,c, and d.
                  ## This path has already been checked to be uncovered and
                  ## to be faithful for the conservative case
		  ucp <- minUncovCircPath(p, pag = apag, path = c(a,c,d,b),
					  unfVect = unfVect, verbose = verbose)
                  ## there is a path ---> orient
                  if (length(ucp) > 1) {
                    ## orient every edge on the path as --
                    n <- length(ucp)
                    apag[ucp[1], ucp[n]] <- apag[ucp[n], ucp[1]] <- 3 ## a--b
                    for (j in seq_len(length(ucp)-1)) {
                      apag[ucp[j], ucp[j + 1]] <- apag[ucp[j + 1], ucp[j]] <- 3 ## each edge on the path --
                    }
                  }
                }
              }
            }
          }
        }
      }
      ##-- R6 ------------------------------------------------------------------
      if (rules[6]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)## b o-* c
        for (i in seq_len(nrow(ind))) {
            b <- ind[i, 1]
            c <- ind[i, 2]
            indA <- which(apag[b, ] == 3 & apag[, b] == 3)
            if (length(indA) > 0) {
              apag[c, b] <- 3
              if (verbose)
                cat("\nRule 6","\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
            }
        }
      }
      ##-- R7 ------------------------------------------------------------------
      if (rules[7]) {
        ind <- which((apag != 0 & t(apag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            b <- ind[i, 1]
            c <- ind[i, 2]
            indA <- which((apag[b, ] == 3 & apag[, b] == 1) &
                          (apag[c, ] == 0 & apag[, c] == 0))
            indA <- setdiff(indA, c)
            if (length(indA) > 0) {
              if (length(unfVect) == 0) { ## normal version

                apag[c, b] <- 3
                if (verbose) {
                  cat("\nRule 7","\n")
                  cat("Orient", indA, "-o", b, "o-*", c, "as", b, "-*", c, "\n")
                }
              }
              else { ## conservative
                for (j in seq_along(indA)) {
                  a <- indA[j]
                  ## check faithfulness of a-b-c
                  if (!any(unfVect == triple2numb(p,a,b,c), na.rm = TRUE) &&
                      !any(unfVect == triple2numb(p,c,b,a), na.rm = TRUE)) {
                    apag[c, b] <- 3
                    if (verbose) {
                      cat("\nRule 7","\n")
                      cat("Conservatively orient:", a, "-o", b, "o-*",
                          c, "as", b, "-*", c, "\n")
                    }
                  }
                }
              }
            }
          }
      }
      ##-- R8 ------------------------------------------------------------------
      if (rules[8]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE)
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            c <- ind[i, 2]
            indB <- which(((apag[a, ] == 2 & apag[, a] == 3) |
                           (apag[a, ] == 1 & apag[, a] == 3)) &
                          (apag[c, ] == 3 & apag[, c] == 2))
            if (length(indB) > 0) {
              apag[c, a] <- 3
              if (verbose) {
                cat("\nRule 8","\n")
                cat("Orient:", a, "->", indB, "->", c,
                    "or", a, "-o", indB, "->", c, "with",
                    a, "o->", c, "as", a, "->", c, "\n")
              }
            }
        }
      }
      ##-- R9 ------------------------------------------------------------------
      if (rules[9]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. a (o-)--(o>) b and b and c are not connected
          indB <- which((apag[a, ] == 2 | apag[a, ] == 1) & (apag[, a] == 1 | apag[, a] == 3) & (apag[c, ] == 0 & apag[, c] == 0))
          ## delete c from indB since it is surely inside
          indB <- setdiff(indB, c)
          ## chose one b s.t. the initial structure exists and the edge hasn't been oriented yet
          while ((length(indB) > 0) && (apag[c,a] == 1)) {
            b <- indB[1]
            indB <- indB[-1]
            ## find a minimal uncovered pd path from initial (a, b, c) :
	    upd <- minUncovPdPath(p, apag, a, b, c,
				  unfVect = unfVect, verbose = verbose)
            ## there is a path ---> orient it
            if (length(upd) > 1) {
              apag[c, a] <- 3
              if (verbose) {
                cat("\nRule 9", "\n")
                cat("There exists an uncovered potentially directed path between", a, "and", c, ". Orient:", a, " ->",c, "\n")
              }
            }
          }## while
        }
      }
      ##-- R10 ------------------------------------------------------------------
      if (rules[10]) {
        ind <- which((apag == 2 & t(apag) == 1), arr.ind = TRUE) ## a o-> c
        while (length(ind) > 0) {
          a <- ind[1, 1]
          c <- ind[1, 2]
          ind <- ind[-1,, drop = FALSE]
          ## find all b s.t. b --> c
          indB <- which((apag[c, ] == 3 & apag[, c] == 2))
          if (length(indB) >= 2) {
            counterB <- 0
            while (counterB < length(indB) && (apag[c, a] == 1)) {
              counterB <- counterB + 1
              b <- indB[counterB]
              indD <- setdiff(indB, b)
              counterD <- 0
              while ((counterD < length(indD)) && (apag[c, a] == 1)) {
                counterD <- counterD + 1
                d <- indD[counterD]
                ## this is the easiest one
                if ((apag[a, b] == 1 || apag[a, b] == 2) &&
                    (apag[b, a] == 1 || apag[b, a] == 3) &&
                    (apag[a, d] == 1 || apag[a, d] == 2) &&
                    (apag[d, a] == 1 || apag[d, a] == 3) &&
                    (apag[d, b] == 0 && apag[b, d] == 0)) {
                  if (length(unfVect) == 0) { ## normal version
                    apag[c, a] <- 3
                    if (verbose)
                      cat("\nRule 10","\nOrient:", a, "->", c, "\n")
                  }
                  else { ## conservative version
                    ## check faithfulness of b-a-d
                    if (!any(unfVect == triple2numb(p,b,a,d), na.rm = TRUE) &&
                        !any(unfVect == triple2numb(p,d,a,b), na.rm = TRUE)) {
                      apag[c, a] <- 3
                      if (verbose)
                        cat("\nRule 10","\nConservatively orient:", a, "->", c, "\n")
                    }
                  }
                }
                ## search with a breitensuche two minimal uncovered circle paths
                else {
                  ## find all x s.t. a (o-)--(o>) x
                  indX <- which((apag[a, ] == 1 | apag[a, ] == 2) &
                                (apag[, a] == 1 | apag[, a] == 3), arr.ind = TRUE)
                  indX <- setdiff(indX, c)
                  if (length(indX >= 2)) {
                    i1 <- 0
                    while (i1 < length(indX) && apag[c, a] == 1) {
                      i1 <- i1 + 1
                      pos.1 <- indA[i1]
                      indX2 <- setdiff(indX, pos.1)
                      i2 <- 0
                      while (i2 < length(indX2) && apag[c, a] == 1) {
                        i2 <- i2 + 1
                        pos.2 <- indX2[i2]
			tmp1 <- minUncovPdPath(p, apag, a, pos.1, b,
					       unfVect = unfVect, verbose = verbose)
			tmp2 <- minUncovPdPath(p, apag, a, pos.2, d,
					       unfVect = unfVect, verbose = verbose)
                        ## we found 2 uncovered pd paths
                        if (length(tmp1) > 1 && length(tmp2) > 1 &&
                            pos.1 != pos.2 && apag[pos.1, pos.2] == 0) {
                          if (length(unfVect) == 0) { ## normal version
                            apag[c, a] <- 3
                            if (verbose)
                              cat("\nRule 10","\nOrient:", a, "->", c, "\n")
                          }
                          else { ## conservative version
                            if (!any(unfVect == triple2numb(p,pos.1, a, pos.2), na.rm = TRUE) &&
                                !any(unfVect == triple2numb(p,pos.2, a, pos.1), na.rm = TRUE)) {
                              apag[c, a] <- 3
                              if (verbose)
                                cat("\nRule 10","\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      } ## rules[10]
    }
  }
  list(graph = apag, sepset = sepset)

} ## {udag2apag}

######################################################################
## Oracle FCI (fast)
#######################################################################

ancTS <- function(g) {

  ## Purpose: compute the ancestors of each node in the graph
  ## ----------------------------------------------------------------------
  ## Arguments: - g: graph in topological order (e.g. produced with randomDAG)
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date:  6 Aug 2012, 10:57

  stopifnot((p <- numNodes(g)) >= 2)
  m <- as(g, "matrix")
  an <- pa <- vector("list", p)
  for (i in 2:p) {
    pa[[i]] <- pa.i <- which(m[1:(i-1),i] != 0)
    if (length(pa.i) > 0) {
      tmp <- c(pa.i)
      ## FIXME: unlist(lapply(pa.i, function(.) an[[.]])) :
      for (p.ij in pa.i)
        tmp <- c(tmp, an[[p.ij]])
      an[[i]] <- sort(unique(tmp))
    }
  }
  an
}

dag2pag <- function(suffStat, indepTest, graph, L, alpha, rules = rep(TRUE,10), verbose = FALSE) {
  ## Purpose: Perform fat oracle FCI-Algorithm, i.e., estimate PAG from
  ##          the true DAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - suffStat, indepTest: info for the tests
  ## - graph: the true DAG
  ## - L: array containing the latent variables (columns in the adjacency matrix)
  ## - alpha: Significance level of individual partial correlation tests
  ## - rules: array of length 10 wich contains TRUE or FALSE corrsponding
  ##          to each rule. TRUE means the rule will be applied.
  ##          If rules==FALSE only R0 (minimal pattern) will be used
  ## - verbose: 0 - no output, 1 - detailed output

  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, 2011

  p <- numNodes(graph)
  cl <- match.call()
  ## find the ancestor sets
  ancList <- ancTS(graph)
  if (verbose) {
    cat("Compute Skeleton\n================\n")
  }
  ## find the skeleton
  skel <- skeleton.dag2pag(suffStat, indepTest, graph,
                           ancList, L, alpha, verbose = verbose)
  G <- as(skel@graph, "matrix")
  sepset <- skel@sepset
  pMax <- skel@pMax
  n.edgetestsSKEL <- skel@n.edgetests
  max.ordSKEL <- skel@max.ord
  allPdsep <- NA
  ## tripleList <- NULL
  n.edgetestsPD <- 0
  max.ordPD <- 0
  allPdsep <- vector("list", p-length(L))

  ## if (sum(G) == 0) {
  ##   Gobject <- new("graphNEL", nodes = as.character(1:(p-length(L))))
  ## }
  ## else {
  ##   colnames(G) <- rownames(G) <- as.character(1:(p-length(L)))
  ##   Gobject <- as(G, "graphNEL")
  ## }
  if (sum(G) != 0)
    colnames(G) <- rownames(G) <- as.character(seq_len(p-length(L)))
  if (verbose) {
    cat("\nDirect egdes:\n=============\nUsing rules:", which(rules),
        "\nCompute collider:\n")
  }
  res <- if (numEdges(skel@graph) > 0)
    udag2pag(pag = G, sepset, rules = rules, verbose = verbose)
  else G
  new("fciAlgo", amat = res, call = cl, n = integer(0),
      max.ord = as.integer(max.ordSKEL),  max.ordPDSEP = as.integer(max.ordPD),
      n.edgetests = n.edgetestsSKEL, n.edgetestsPDSEP = n.edgetestsPD,
      sepset = sepset, pMax = pMax, allPdsep = allPdsep)
} # {dag2pag}

##' Perform undirected part of oracle FCI-Algorithm, i.e.,
##' find skeleton of PAG given the true DAG
##'
##' @title
##' @param suffStat  info for the tests
##' @param indepTest info for the tests
##' @param graph the true DAG
##' @param ancList list containing the ancestors of each node in the graph
##' @param L array containing the latent variables (columns in the adjacency matrix)
##' @param alpha Significance level of individual partial correlation tests
##' @param NAdelete delete edge if pval=NA (for discrete data)
##' @param verbose 0 - no output, 1 - detailed output
##' @return "pcAlgo" object : (G, sepset, pMax, ord, n.edgetests)
##' @author  Diego Colombo, 2011; speedup: Martin Maechler, Nov.2013
skeleton.dag2pag <- function(suffStat, indepTest, graph, ancList, L, alpha,
                             NAdelete = TRUE, verbose = FALSE)
{
  ## Currently exported but undocumented, used only in  dag2pag()

  new.ord <- function(start, old, L) {
    seq_start <- 1:start
    tmp <- setdiff(seq_start,L) ## <- independent of old -- FIXME speed
    new <- rep(0, length(old))
    for (i in seq_along(old)) {
      new[i] <- which(tmp == old[i])
    }
    return(new)
  }

  cl <- match.call()
  pval <- NULL
  p <- numNodes(graph)
  stopifnot(2 <= (new.p <- p-length(L)))
  sepset <- vector("list", new.p)
  G <- matrix(TRUE, nrow = new.p, ncol = new.p)
  diag(G) <- FALSE
  for (iList in 1:new.p) sepset[[iList]] <- vector("list", new.p)
  pMax <- matrix(-Inf, nrow = new.p, ncol = new.p)
  diag(pMax) <- 1
  ## test independence given the ancestors
  for (i in setdiff(1:(p-1), L)) { ## i.e. i is not in L
    iL <- new.ord(p,i,L)
    for (j in setdiff((i+1):p, L)) { ## i.e. j is not in L
        cond.set <- c(ancList[[i]],ancList[[j]])
        cond.set <- setdiff(cond.set,c(i,j,L))
        pval <- indepTest(i, j, cond.set, suffStat)
        jL <- new.ord(p,j,L)
        if (verbose) {
          cat("x=", iL, " y=", jL,
              " S=", new.ord(p,cond.set,L),": pval =", pval, "\n")
        }
        if(is.na(pval)) pval <- as.numeric(NAdelete) ## == if(NAdelete) 1 else 0
        if (pval > pMax[iL, jL]) {
          pMax[iL, jL] <- pval
        }
        if (pval >= alpha) {
          G[iL, jL] <- G[jL, iL] <- FALSE
          sepset[[iL]][[jL]] <- new.ord(p,cond.set,L)
        }
    }
  }
  for (i in 1:(new.p - 1)) {
    for (j in 2:new.p) {
      pMax[i, j] <- pMax[j, i] <- max(pMax[i, j], pMax[j, i])
    }
  }
  Gobject <- if (sum(G) == 0)
      new("graphNEL", nodes = as.character(1:new.p))
  else {
      colnames(G) <- rownames(G) <- as.character(1:new.p)
      as(G, "graphNEL")
  }
  new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
      max.ord = integer(0), n.edgetests = integer(0),
      sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}

iplotPC <- function(pc.fit, labels = NULL) {
  ## igraph alternative to Rgraphviz (only for pcAlgo objects)
  adjm <- wgtMatrix(getGraph(pc.fit), transpose = FALSE)
  if (!is.null(labels)) dimnames(adjm) <- list(labels, labels)
  g1 <- graph.adjacency( adjm )
  plot.igraph(g1)
}

showEdgeList <- function(object, labels = NULL)
{
    cat("\nEdge List: \n")
    g <- getGraph(object)
    if (is.null(labels)) labels <- nodes(g)
    isDir <- isDirected(g)
    wm <- wgtMatrix(g)

    if(isDir && !is(object,"pcAlgo"))
      stop("implementation only for 'pcAlgo', currently ...") ##-> Markus

    ## MM: in general, maybe rather look at each pair (i,j) in the lower triangular part of [W, t(W)]
    wmU <- wm + t(wm)   # weights for undirected e.
    wmD <- t(wm - t(wm))# weights for   directed e.
    u <- which(wmU == 2 & upper.tri(wmU), arr.ind = TRUE)
    cat("\n", if(isDir) "Undirected Edges" else "Undirected graph with Edges",
	":\n", sep = "")
    for (i in seq_len(nrow(u)))
      cat(" ", paste(labels[u[i,1]], " --- ", labels[u[i,2]], "\n"))
    if(isDir) {
      d <- which(wmD == 1, arr.ind = TRUE)
      d <- d[order(d[,1]),]
      cat("\nDirected Edges:\n")
      for (i in seq_len(nrow(d)))
	cat(" ", paste(labels[d[i,1]], " --> ", labels[d[i,2]], "\n"))
    } else {
      d <- matrix(0,0,0)
    }
    invisible(list(undir = u, direct = d))
}

showAmat <- function(object) {
  g <- getGraph(object)
  cat("\nAdjacency Matrix G:",
      "G[i,j] = 1/2 if edge mark of edge i-j at j is head/tail.",
      "", sep = "\n")
  wm <- wgtMatrix(g)
  mTmp <- t(wm - 2*t(wm))
  mTmp[ mTmp < 0 ] <- 2
  mTmp
}


#######################################################################
## Generalized backdoor criterion
#######################################################################

visibleEdge <- function(amat, x, z)
{
  ## Purpose: check if the directed edge from x to z in a MAG or in a PAG
  ##          is visible or not
  ## ----------------------------------------------------------------------
  ## Arguments: amat, x, z
  ## ----------------------------------------------------------------------
  ## Value: T/F
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 25 Apr 2012;  simplification: Martin Maechler

  ## 1. scenario: there exists a vertex not adjacent to z with *---> x
  ## c *-> x --> z and c and z not connected :
  hasC1 <- any(amat[x,] != 0 & amat[,x] == 2 & amat[z,] == 0)
  ## if there is at least one vertex
  if (any(hasC1)) ## the edge is visible
      return(TRUE)

  ## Otherwise:
  ## 2. scenario: there exists a collider path that is into x and
  ##		  every vertex on the path is a parent of z

  ## c <--> x --> z and c is a parent of z :
  indC2 <- which(amat[x,] == 2 & amat[,x] == 2 &
		 amat[z,] == 3 & amat[,z] == 2)
  for(c in indC2) {
    ## find a minimal discriminating path for c,x,z
    if (length(minDiscrPath(amat, c,x,z)) > 1)
      ## a path exists
      return(TRUE)
  }
  ## nothing found:  return
  FALSE
}

possibleDe <- function(amat,x)
{
    ## Purpose: in a DAG, CPDAG, MAG, or PAG determine which nodes are
    ##          possible descendants of x on definite status paths
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## - amat: matrix corresponding to the DAG, CPDAG, MAG, or PAG
    ## - x: node of interest
    ## ----------------------------------------------------------------------
    ## Value:
    ## - de.list: array containing the possible descendants of x
    ## ----------------------------------------------------------------------
    ## Author: Diego Colombo, Date: 26 Apr 2012, 16:58

    stopifnot(is.matrix(amat))
    p <- nrow(amat)
    is.de <- rep.int(FALSE, p) ##
    ## 1. case: x is a possible child of itself
    is.de[x] <- TRUE
    ## 2. case: find all the possible children of x
    indD <- which(amat[x,] != 0  & amat[,x] != 2 & !is.de) ## x (o,-)-* d
    i.pr <- rep(x,length(indD))
    while (length(indD) > 0) {
        ## next element in the queue
        d <- indD[1]
        indD <- indD[-1]
        pred <- i.pr[1]
        i.pr <- i.pr[-1]
        is.de[d] <- TRUE
        a.d <- amat[,d]
        a.d.p <- a.d[pred]
        ## find all possible children of d not visited yet
        indR <- which(amat[d,] != 0 & a.d != 2 & !is.de) ## d (o,-)-* r
        for(j in seq_along(indR)) {
          ## check that the triple <pred,d,r> is of a definite status
          ## 1. d is a collider on this subpath; this is impossible
          ##    because the edge between d and r cannot be into d
          ## 2. d is a definite non-collider
          r <- indR[j]
          if (a.d.p == 3 || a.d[r] == 3 ||
              (a.d.p == 1 && a.d[r] == 1 && amat[pred,r] == 0)) {
            ## update the queues
            indD <- c(indD, r)
            i.pr <- c(i.pr, d)
          }
        }
    }
    ## return 'de.list' :
    which(is.de)

} ## {possibleDe}


dreach <- function(x, y, amat, verbose = FALSE)
{
  ## Purpose: Compute d-sep(x,y) ("dsep")
  ## !!!!! WE DO NOT ASSUME SELECTION VARIABLES HENCE NO --- EDGES IN A MAG!!!!
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - x, y: nodes of which dsep must be computed
  ## - amat: adjacency matrix of a MAG
  ##         amat[i,j] = 0 iff no edge btw i,j
  ##         amat[i,j] = 3 iff i *-- j
  ##         amat[i,j] = 2 iff i *-> j
  ## - verbose: Show checked node sequence
  ## ----------------------------------------------------------------------
  ## Value:
  ## - DSEP: set containing the nodes in D-SEP(x,y)
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 31 Jan 2013, 11:04

  ## check: quadratic; x in V; edgemarks ok; non-zeroes symmetric
  stopifnot(is.matrix(amat), (p <- nrow(amat)) == ncol(amat),
            x <= p, amat %in% c(0,2,3),
            ## The non-zero entries in amat must be symmetric:
            (aN0 <- amat != 0) == t(aN0))

  ## compute the ancestors of x and y
  amat.an <- amat
  amat.t <- amat + t(amat)
  ## because we have only ---> and <--> edges, in amat.t we have either 5s (--->) or 4s (<-->)
  ## delete all the bi-directed edges
  amat.an[which(amat.t == 4)] <- 0
  ## transform the directed edges from 3 and 2 in 1 and 0
  amat.an[which(amat.an == 3)] <- 0
  amat.an[which(amat.an == 2)] <- 1
  dir.graph <- as(amat.an, "graphNEL")
  john.pairs <- johnson.all.pairs.sp(dir.graph)
  anc.y <- anc.x <- rep(FALSE, p)
  anc.y[which(john.pairs[,y] < Inf)] <- TRUE
  anc.x[which(john.pairs[,x] < Inf)] <- TRUE

  ## the set containing the ancestors of either x or y
  anc <- anc.x | anc.y
  anc[x] <- FALSE

  ## prepare the setting for the search
  amat.tmp <- amat
  nb <- which(amat[x,] != 0 & anc)
  Q <- nb
  P <- rep.int(x,length(Q))
  DSEP <- nb

  amat.tmp[x,nb] <- 0 ## delete edge to nbrs

  while(length(Q) > 0) {
    ## Invariants:
    ## ===========
    ## (A1) length(Q) == length(P) > 0
    ## (A2) non-zero in amat.tmp -> non-zero in amat [no non-zero added]
    ## (A3) Q[i] and P[i] are adjacent in amat [using (A2)]
    if (verbose) {
      cat("\n-------------","\n")
      cat("Queue Q:",Q,"\n")
      cat("Queue P:",P,"\n")
    }
    a <- Q[1]
    Q <- Q[-1]
    pred <- P[1] ## not empty because of (A1)
    P <- P[-1]
    if (verbose) cat("Select",pred,"towards",a,"\n")
    nb <- which(amat.tmp[a,] != 0 & anc) ## (*)

    for (i in seq_along(nb)) {
        b <- nb[i]
        ## Guaranteed: pred-a-b are a path because of (A3)
        ## if a is a collider, b is in D-SEP of x
        if (amat[pred,a] == 2 && amat[b,a] == 2) {
          amat.tmp[a,b] <- 0 ## remove b out of adj(a) in amat.tmp (+)
          Q <- c(Q,b)
          P <- c(P,a)
          DSEP <- c(DSEP,b)
        }
    }

  }
  sort(unique(DSEP))

} ## {dreach}

pag2magAM <- function(amat.pag, x, max.chordal = 10, verbose = FALSE)
{
  ## Purpose: transform a PAG/CPDAG into a valid MAG/DAG using the arrowhead
  ##          augmentation from Zhang, without orienting additional edges
  ##          into x and then orient any chordal component into a DAG
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat.pag: adjacency matrix of the PAG/CPDAG
  ## - x: node of interest
  ## ----------------------------------------------------------------------
  ## Value:
  ## - valid.DAG.mat: adjacency matrix corresponding to the valid MAG
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 12 Apr 2013, 14:24;
  ## Tweaks by Martin Maechler

  ## 1. step: arrowhead augmentation
  ## find all o-> edges in the PAG
  indA <- which((amat.pag == 2 & t(amat.pag) == 1), arr.ind = TRUE)
  for (i in seq_len(nrow(indA))) {
      a <- indA[i, 1]
      b <- indA[i, 2]
      amat.pag[b,a] <- 3L
  }
  ## because we don't allow selection variables in the generalized backdoor
  ## criterion these kind of edges should not be present BUT YOU NEVER KNOW...
  ## find all o-- edges in the PAG
  indB <- which(amat.pag == 3 & t(amat.pag) == 1, arr.ind = TRUE)
  for (i in seq_len(nrow(indB))) {
      a <- indB[i, 1]
      b <- indB[i, 2]
      amat.pag[b,a] <- 2L
  }
  ## 2.step: find all the chordal components in the updated graph
  ## find the adjacency matrix that contains
  amat.undir <- amat.dir <- amat.pag
  amat.t <- amat.pag + t(amat.pag)
  ## we are interested only in the o-o edges
  ## delete all the other edges
  amat.undir[amat.t != 2] <- 0L
  amat.dir  [amat.t == 2] <- 0L
  g.undir <- as(amat.undir,"graphNEL")

  ## Get all connected components [maybe should rather keep labels ??]
  conn.comp <- lapply(connectedComp(g.undir), as.numeric)

  if(verbose) {
      lc <- vapply(conn.comp, length, 1L)
      t2 <- if(any(lc == 1)) gettextf(" and %d singletons", sum(lc == 1)) else ""
      cat(gettextf("The undirected graph has %d connected components%s",
                   sum(lc > 1), t2),"\n")
  }
  valid.DAG.mat <- amat.undir

  ## get all semilocal extensions of the undirected cpdag
  for (cci in conn.comp) if (length(cci) > 1) {
      if(length(cci) > max.chordal) return(NULL)
      ## check whether the component is chordal
      am.u.ii <- amat.undir[cci,cci]
      if(!is.chordal(graph.adjacency(am.u.ii, "undirected"),
		     fillin = TRUE)$chordal)
        return(NULL)
      ## else: not extendable or too large to handle
      am. <- my.SpecialDag(amat.undir, a=am.u.ii, X=x, verbose=verbose)
      valid.DAG.mat <- valid.DAG.mat * am.
    } # if() ..end for

  ## add directed edges and return "valid.DAG.mat":
  valid.DAG.mat + amat.dir

} ## {pag2magAM}

##' Auxiliary for  pag2magAM()
my.SpecialDag <- function (gm, a, X, verbose = FALSE)
{
    while (sum(a) != 0 ) {
        sinks <- find.sink(a)
        if (verbose) {
            cat("Main Call: ################## \n")
            print(gm)
            print(a)
            cat("Sinks: ", sinks, "\n")
        }
        for (x in sinks) {
            if (verbose)
                cat("Try removing", x, " in a.\n")
            if (adj.check(a, x) & as.numeric(rownames(a)[x]) != X) {
                inc.to.x <- a[, x] == 1 & a[x, ] == 1
                if (any(inc.to.x)) {
                    real.inc.to.x <- as.numeric(rownames(a)[inc.to.x])
                    real.x        <- as.numeric(rownames(a)[x])
                    gm[real.x, real.inc.to.x] <- 3
                    gm[real.inc.to.x, real.x] <- 2
                }
                if (verbose) {
                    cat("Removed sink", as.numeric(rownames(a)[x]),
                        "in g (", x, "in a).\n")
                    cat("New graphs: \n")
                    print(gm)
                    print(a)
                }
                a <- a[-x, -x]
                break
            }
        }
    }
    gm
} ## {my.SpecialDag}


backdoor <- function(amat, x, y, type = "pag", max.chordal = 10, verbose = FALSE)
{
  ## Purpose: for a given pair of nodes (x,y) and a given graph (DAG, CPDAG,
  ##          MAG, or PAG) estimate if the causal effect of x on y using the
  ##          generalized backdoor criterion is identifiable or not. In a
  ##          first step we check if the effect is identifiable or not. If
  ##          it is identifiable we compute the set W that satisfies the
  ##          generalized backdoor criterion. If it is not identifiable
  ##          the output is NA
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - amat: adjacency matrix of the corresponding graph
  ## - x, y: nodes for which we want to estimate the causal effect of x on y
  ## - type: specifies the type of graph of the adjacency matrix amat. It
  ##         can be a DAG (type="dag"), a CPDAG (type="cpdag"), a MAG
  ##         (type="mag"), or a PAG (type="pag")
  ## ----------------------------------------------------------------------
  ## Value:
  ## - cEffect: causal effect of x on y (either NA not identifiable or a number)
  ## ----------------------------------------------------------------------
  ## Author: Diego Colombo, Date: 12 Apr 2013, 16:06

    stopifnot (dim(amat)[1] > 0, x != y, length(type) == 1)
    set.w <- NA
    if (type == "dag" || type == "cpdag") {
        ## transfor each directed edge from 0-1 to 2-3 in the adjacency matrix
        ind <- which((amat == 1 & t(amat) == 0), arr.ind = TRUE) ## a -> b
        for (i in seq_len(nrow(ind))) {
            a <- ind[i, 1]
            b <- ind[i, 2]
            amat[a,b] <- 2
            amat[b,a] <- 3
        }
    }

    ## check that x and y belong to the same connected component
    ## compute the connected components of the whole graph
    conn.comp <- lapply(connectedComp(as(amat, "graphNEL")),
                        as.numeric)
    found.conncomp <- FALSE
    for(ic in seq_along(conn.comp))
        if (x %in% conn.comp[[ic]] & y %in% conn.comp[[ic]]) {
            found.conncomp <- TRUE
            break ## we found one {FIXME: should use info conn.comp[[ic]] ! }
        }

    ## if x and y are not in the same connected component ---> empty set
    if (!found.conncomp)
        return( NULL )

    ## else :

    ## 1. generate a MAG/DAG belonging to the equivalence class of the
    ## PAG or a CPDAG without orienting additional edges into X
    amat.r <- if (type == "cpdag" || type == "pag")
        pag2magAM(amat, x, max.chordal = max.chordal, verbose = verbose)
    else amat

    ## 2. compute the truncated corresponding graph
    ## 2.a) find all the visible edges in the original graph that are out of X
    ## 2.b) all the definitely visible edges out of X have to be deleted
    ## in the graph constructed previously
    indD <- which(amat[x,] == 2 & amat[,x] == 3) ## x--> d

    ## CASE A: in a CPDAG or DAG every directed edge out of X is visible
    ##        then delete them all
    if (type == "cpdag" || type == "dag") {
	for (z in indD)
	    ## delete visible edges out of x edges
	    amat.r[z,x] <- amat.r[x,z] <- 0
    } else {
        ## CASE B: in a MAG or PAG the directed edges out of X need
        ##        to be checked for visibility
        if (length(indD) > 0) {
	    del.edges <- vapply(indD, function(z) {
		## check each directed edge out of X for visibility
		visibleEdge(amat, x, z)
	    }, NA)
            ## delete visible edges out of x edges
            amat.r[indD[del.edges],x] <- amat.r[x,indD[del.edges]] <- 0
            ## transform invisible edges out of x by bi-directed
            ## amat.r[indD[!del.edges],x] <- amat.r[x,indD[!del.edges]] <- 2
        }
    }

    ## 3. compute possible descendants of x along definite status paths
    list.de <- possibleDe(amat, x)

    ## 4. compute D-SEP(x,y)_path in the truncated graph
    dsep.set <- dreach(x, y, amat.r)

    ## 5. check that no possible descendants of x along definite status
    ## paths are in D-SEP(x,y,amat.r) or y is in adj(x, amat.r)
    if (amat.r[x,y] == 0 && !any(list.de %in% dsep.set))
        ## D-SEP(x,y,amat.r) closes all the paths
        set.w <- dsep.set

    set.w
} ## {backdoor}

##' Utility that transforms a CPDAG matrix into a FCI matrix
trafoCPDmat <- function(M)
{
  n <- ncol(M)
  if(n >= 2) for(i in 2:n)
    for(j in seq_len(i-1L)) { ## 1 <= j < i <= n
      if (M[i,j] == 1L) {
        ## if the edge is i -> j
        if (M[j,i] == 0L) {
          M[i,j] <- 2L  # arrow head
          M[j,i] <- 1L  # arrow circle
        }
        ## if the edge is i - j  --- nothing to do: M[i,j] = M[j,i] = 1 already
        ## else if (M[j,i] == 1)
        ## {
        ##   M[i,j] <- 1  # arrow circle
        ##   M[j,i] <- 1  # arrow circle
        ## }
      }
      ## if the edge is i <- j
      else if (M[i,j] == 0L && M[j,i] == 1L)
      {
        M[i,j] <- 1L  # arrow circle
        M[j,i] <- 2L  # arrow head
      }
    }
  M
}

## this function finds ancestors of a node x in the subgraph of G+
## when some confounders are removed,
## for i=1..p (p - number of nodes in the full graph)
## vector[i] = position of the node i in the subgraph (if the node is not in
## the subgraph vector[i]=0)
## m is the adjacency matrix of the full graph
###___ CURRENTLY UNUSED ___,  instead  NAA_ancestors() below
NAAancestors.1 <- function(x, m, vector)
{
  ## q denotes unvisited nodes/ nodes in queue
  ## v denotes visited nodes
  p <- nrow(m)
  q <- v <- integer(p) # all 0
  q[1L] <- x
  i <- k <- 1L
  while(k <= i && q[k] != 0)
  {
    t <- q[k]
    ## mark t as visited
    v[k] <- t
    k <- k+1L
    ## in this for cycle adds all nodes that have an undirected
    ## edge with node t and all parents of node t to queue
    for(j in 1:p)
      if (m[j,t] == 2 && m[t,j] == 3)
        ## only add nodes that haven't been added
        if (vector[j] != 0 && j %nin% q)
        {
          i <- i+1L
          q[i] <- j
        }
  }
  ## Not Against Arrow Ancestors

  ## If node a is in NAAA that means that there is a not against arrowhead path
  ## from a to x
  setdiff(v, c(0L,x))
}

##' Finds the "not against arrowhead ancestors" of the node x in G+
##' using the adjacency matrix m of G+, and breadth first search.
##' "Not Against Arrowhead Ancestors" (NAAA) are all the nodes u in G+ that
##' have a path u.. -> x in G+  which doesn't go against arrowhead
##'
##' Currently only used as auxiliary for  PosDsepLinks()
##'
##' @title Find "Not Against Arrowhead" Ancestors (NAAA)
##' @param x scalar integer in 1..p, denoting a node in graph G+
##' @param m p x p adjacency matrix of graph G+
##' @return a (possibly empty) vector of integers, the NAAA nodes of x in G+
NAA_ancestors <- function(x, m)
{
  ## q denotes unvisited nodes/ nodes in queue
  ## v denotes visited nodes
  p <- ncol(m); i.p <- 1:p
  q <- v <- integer(p) # all 0
  q[1L] <- x
  i <- k <- 1L
  while(k <= i && q[k] != 0) {
    t <- q[k]
    ## mark t as visited
    v[k] <- t
    k <- k+1L
    ## in this for cycle adds all nodes that have an undirected
    ## edge with node t and all parents of node t to queue
    for(j in i.p)
      if ((t == x && m[j,t] == 2 && m[t,j] == 1) ||
          (t != x && m[t,j] != 2 && m[t,j] != 0))
        ## only add nodes that haven't been added
        if (j %nin% q) {
          i <- i+1L
          q[i] <- j
        }
  }
  ## return  Not Against Arrow Ancestors (NAAA):
  ## A node a is in NAAA iff there is a "not against arrowhead" path from a to x
  setdiff(v, c(0L,x))
}

##' Find possible Dsep links, given the adjacency matrix,
##' as the edges that satisfy the pattern described in Lemma 4.
##'
##' The edge X <-> Y is a PosDsep link in G+ if there exist nodes, U,V such that
##'  U <-> X <-> Y <-> V in G+, and U and V are not adjacent and paths V.. -> X
##' U.. -> Y exist in G+ and they do not go against arrowhead
##' @title Find POSsible D-Separation LINKS
##' @param m adjacency matrix
##' @return a (possibly empty) data frame with columns "x" and "y"
##'	 where  \code{x[i] *-* y[i]} is a possible Dsep link for all i.
PosDsepLinks <- function(m, verbose=TRUE)
{
  stopifnot(1 <= (p <- ncol(m)))
  i.p <- 1:p
  NAAA <- vector("list", p)
  x <- y <- integer()
  ## for all pairs  1 <= j < i <= p
  for(i in i.p[-1L]) {
    ## nodes that have a  "not against arrowhead" path to i :
    NAAA[[i]] <- NAAA_i <- NAA_ancestors(i, m)
    for (j in 1:i)
      ## first find a bidirected i <-> j edge in G+
      if (m[i,j] == 2 && m[j,i] == 2) {
        u <- v <- integer()
        ## Then find bidirected edges u <-> i and j <-> v
        for (k in i.p) if (k != i && k != j)
        {
          if (m[i,k] == 2 && m[k,i] == 2)
            u <- union(u,k)
          if (m[j,k] == 2 && m[k,j] == 2)
            v <- union(v,k)
        }

        ## find nodes u (v) that have both a bidirected edge with i (j) and
        ## a not against arrowhead path to j (i)
        u <- intersect(u, NAAA[[j]]) ## as j <= i, this is already computed
        v <- intersect(v, NAAA_i)

        ## check if there is a pair of nodes in (u,v) that is not adjacent
        found.i.j <- FALSE
        for(k in seq_along(u)) {
          u.k <- u[k]
          for(r in seq_along(v))
            if (u.k != v[r] && m[u.k,v[r]] == 0) {
              found.i.j <- TRUE
              break
            }
          if(found.i.j)
            break
        }

        ## found.i.j is true iff we found a pair of nodes fulfilling all 3 conditions.
        ## In that case, the node pair (i,j)  is added to the  PosDsepLinks set:
        if (found.i.j) {
          if(verbose) cat(sprintf("Added PosDsepLink  %2d *-* %2d\n", i,j))
          x <- c(x, i)
          y <- c(y, j)
        }
      }
  }

  ## returns data frame with possible Dsep links x[i] *-* y[i] for all i
  data.frame(x = x, y = y)
}

## this function finds a minimal Dsep set, given a Dsep set by
## going through all the subsets of the Dsep set, in the same way it is
## done as in the skeleton function

## x and y are the nodes separated by the sepset sep,
## while suffStat is a sufficient statistic, and
## indepTest is an independence test
MinimalDsep <- function(x,y, sep, suffStat,indepTest, alpha = 0.01)
{
  stopifnot((n <- length(sep)) >= 1)
  for(i in 1:n)
  {
    S <- seq_len(i)
    wasLast <- FALSE
    while (!wasLast) {
      if (indepTest(x,y,sep[S],suffStat) > alpha) ## had 'alpha= 0.01' hardcoded
        return(sep[S])

      z <- getNextSet(n,i,S)
      S <- z$nextSet
      wasLast <- z$wasLast
    }
  }
}

##' Find the hierarchy HIE(X,I) given a vector of nodes x and a list of
##' sepsets, according to the definition of a hierarchy HIE(X,I).
##'
##' This is a recursive function, so even though vector x represents the
##' original vector given to the function, it also represents the hierarchy
##' we are in the process of constructing
##' @title Find Hierarchy HIE(X,I)
##' @param x a vector of nodes.
##' @param sepset a list of \dQuote{sepsets}, the "independence set".
##' @return
HIE <- function(x, sepset)
{

  ## flag1 is used to detect wheteher the hierarchy
  ## is complete or if more nodes should be added
  flag1 <- FALSE

  ## first check if all the nodes have been added
  ## if there are nodes that still need to be added
  ## set flag1 to TRUE and exit the loop
  i.x <- seq_along(x)
  for(i in i.x)
  {
    if (flag1) break
    for(j in seq_along(x)) {
      ss.ij <- sepset[[x[i]]][[x[j]]]
      if (!is.null(ss.ij) && length(ss.ij) != 0 && length(setdiff(ss.ij, x)) != 0) {
        flag1 <- TRUE
        break
      }
    }
  }

  ## if flag1 is still FALSE all nodes have been added
  ## and we can exit the function
  if (!flag1)
    return(x)

  ## else (flag1 is TRUE) there are still nodes to be added

  ## flag2 is an indicator that prevents the recursion
  ## from doing some calculations more than once
  flag2 <- FALSE
  for(i in i.x) {
    if (flag2) break
    for(j in i.x)
    {
      if (flag2)
        break

      ## else :
      ## if there are still nodes in the separating set of two nodes
      ## from the hierarchy that are not in the hierarchy
      ## add those nodes to the hierarchy
      ss.ij <- sepset[[x[i]]][[x[j]]]
      if (!is.null(ss.ij) && length(ss.ij) != 0 && length(setdiff(ss.ij, x)) != 0) {
        ## recurse
        return(HIE(union(x,ss.ij), sepset))
      }
    }
  }
} ## {HIE}

AugmentGraph <- function(M, suffStat, sepsets, indepTest, alpha = 0.01)
{
  ## first transform matrix so that it differentiates between edge marks
  stopifnot((p <- ncol(M)) >= 1)
  if(hasNms <- !is.null(dns <- dimnames(M))) dimnames(M) <- NULL # speedup
  i.p <- 1:p
  for(i in i.p) {
    seps.i <- sepsets[[i]]
    for(j in i.p) ## go through all the sepsets
      if (!is.null(sep <- seps.i[[j]])) {
        ## only check neighbors of i, j and sepset
        adjacent <- union(which(M[i,] != 0),
                          which(M[j,] != 0))
        for (r in seq_along(sep))
          adjacent <- union(adjacent, which(M[sep[r],] != 0))

        ## delete i,j, sepset from the set of adjacent nodes
        del.nodes <- union(c(i,j), sep)
        adjacent <- setdiff(adjacent, del.nodes)

        ## if adding a node to the sepset breaks the cond. independence
        ## orient according to lemma 2. (i)
        for (adj in adjacent) {
          if (indepTest(i, j, union(sep,adj), suffStat) < alpha) {
            for(node in del.nodes)
              if ((M.k <- M[node,adj]) != 0 && M.k != 2) {
                M[node,adj] <- 2
                ## cat("oriented this node:",node," towards this node:",adj,"\n")
              }
          }
        }
      }
  }
  ## return newly oriented matrix
  if(hasNms) structure(M, dimnames = dns) else M
}

## The final fciplus function, it uses all other functions and returns
## an adjacency matrix; however, it does not go through the 10 orientation rules.
##
## the input for the function is a pc fitted graph - pc.fit
## as well as a sufficient statistic and an independence test
fciplus.intern <- function(pc.fit, alpha = 0.01, suffStat, indepTest, verbose=TRUE)
{
  sepsets <- pc.fit@sepset
  cpdag <- pc.fit@graph
  m <- trafoCPDmat(as(cpdag, "matrix"))

  ## first run the augment graph function to orient invariant
  ## arrowheads according to Lemma 2. (1)
  mat <- AugmentGraph(m, suffStat, sepsets, indepTest)
  ## which(mat[,27]==2)
  ## Check if there are any possible Dsep links
  link <- PosDsepLinks(mat)

  ## if there are possible Dsep links it should be checked
  ## which of these are actually edges that should be removed
  ## an edge X - Y is a Dsep link if X and Y are Dseparated by the hierarchy
  ## of the parents of the nodes X and Y (excluding X and Y)

  ## a  counter to help us go through the data frame of all possible Dsep lins
  count <- 1L
  while (count <= length(link$x)) {

    x <- link$x[count]
    y <- link$y[count]

    ## basex and basey are vectors of the nodes neighboring nodes x,y
    basex <- setdiff(which(mat[x,] != 0), y)
    basey <- setdiff(which(mat[y,] != 0), x)
    all <- union(basex,basey)

    ## to determine whether a posDsep link is an actual Dsep link
    ## it is necessarry to go through all the subsets of the sets basex and
    ## basey (since we don't know the parents) and build a hierarchy
    ## for each combination of those subsets

    ## since it is complicated to go through all the subets of both basex and basey,
    ## I decided to combine them into one set and go through subsets of that set
    ## but only building a hierarchy using a subset that contatins neighbors of
    ## both x and y

    ## bound is the number of neighbors of x and y
    bound <- length(all)

    ## flag is an indicator that will be given the value TRUE if
    ## an actual Dsep link is discovered so we can use it to break
    ## from the for cycle and return to the while cycle
    flag <- FALSE
    for(i in seq_len(bound)) {
      if (flag) break
      S <- seq_len(i)
      exit <- FALSE
      while (!exit && !flag)
      {
        a.S <- all[S]
        ## check if there are more than 1 and less than k neighbors of x and y in the subset
        if (1 <= (l.S.bx <- length(intersect(a.S, basex))) && l.S.bx <= i &&
            1 <= (l.S.by <- length(intersect(a.S, basey))) && l.S.by <= i)
        {
          ## build a hierarchy
          potential <- HIE(union(c(x,y), a.S), sepsets)
          ## remove x and y from the hierarchy
          potential <- setdiff(potential, c(x,y))

          ## if this set is actually a separating set it should be added to
          ## the list of sepsets and we should begin the process again (AugmentGraph, PosDsepLink)
          if (indepTest(x,y, potential, suffStat) > alpha)
          {
            ## find the minimal sepset
            mindsep <- MinimalDsep(x,y,potential,suffStat,indepTest)

            ## update sepsets
            sepsets[[x]][[y]] <- mindsep

            ## update matrix by removing the Dsep link
            mat[x,y] <- mat[y,x] <- 0

            cat("Found Dsep Link for x=",x,"y=",y,"sepset=",mindsep,"\n")

            ## orient invariant arrowheads
            mat <- AugmentGraph(mat,suffStat,sepsets,indepTest)

            ## search for new PosDsepLinks
            link <- PosDsepLinks(mat)

            ## reset counter and set flag to TRUE
            count <- 1L
            flag <- TRUE
          }
        }

        z <- getNextSet(bound, i, S)
        S <- z$nextSet
        exit <- z$wasLast
      } ## while(!exit ..)
    } ## for(i .. )
    ## if we've exited the for cycle and flag is still false it means x-y is not
    ## a Dseplink so we should move on to the next posDseplink
    if (!flag)
      count <- count+1L
  } ## while(count <= ..)
  ## return the adjusted adjacency matrix and sepset
  list(mat = mat, sepset = sepsets)
} ## {fciplus.intern}

fciPlus <- function(suffStat, indepTest, alpha, labels, p, verbose=TRUE)
{
  ## Author: Markus Kalisch, Date:  7 Jul 2014, 12:08
  cl <- match.call()
  if(!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  if(missing(labels)) {
    ## FIXME ---- make function for this and use in skeleton(), pc(), fci(), ....
    ## if(is.matrix(C <- suffStat[["C"]]) && (d <- dim(C))[1] == d[2]) {
    ##   ## we can derive 'p'  *and* 'labels' -- in 99% of cases !
    ##   p <- d[1]
    ##   if(is.null(labels <- colnames(C)))
    ##      labels <- as.character(seq_len(p))
    ## } else {
    if(missing(p)) stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
    ## }
  } else { ## use labels ==> p  from it
    stopifnot(is.character(labels))
    if(missing(p)) {
      p <- length(labels)
    } else if(p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else
      message("No need to specify 'p', when 'labels' is given")
  }

  skel <- skeleton(suffStat = suffStat, indepTest = indepTest, alpha = alpha,
                   labels = labels, p = p)
  fit1 <- udag2pdagRelaxed(gInput = skel, orientCollider = FALSE)
  fcip <- fciplus.intern(pc.fit = fit1, alpha=alpha, suffStat=suffStat,
                         indepTest=indepTest, verbose=verbose)
  fciplus.amat <- udag2pag(pag = fcip$mat, sepset = fcip$sepset,
                           orientCollider = FALSE)
  colnames(fciplus.amat) <- rownames(fciplus.amat) <- labels
  new("fciAlgo", amat = fciplus.amat, call = cl, n = integer(0),
      max.ord = integer(0),
      max.ordPDSEP = integer(0),
      n.edgetests = integer(0), n.edgetestsPDSEP = integer(0),
      sepset = list(), pMax = matrix(0,1,1), allPdsep = list())
}

###  MM: (ess-set-style 'DEFAULT) : we have much nesting ==> only indent by 2
## Local Variables:
## eval: (ess-set-style 'DEFAULT 'quiet)
## delete-old-versions: never
## End:

