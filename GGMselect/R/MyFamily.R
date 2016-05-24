selectMyFam <- function(X,MyFamily,K=2.5, min.ev=10**(-8)){
  # ---------------------------------------------------------------
  # FUNCTION
  #   Select a graph within a given family of graphs.
  #   End-user function.
  # INPUT
  # X : nxp data matrix
  # MyFamily : list of pxp adjacency matrix  
  # K : scalar or vector larger than 1
  # (K=0 works but is not advised)
  # min.ev : small positive real number
  # OUPUT
  #   A list with components:
  #       crit.min, Neighb, G: see selectFast
  # SEEALSO
  #   selectQE selectFast
  # ---------------------------------------------------------------
  # Liberer la memoire on exit
  # car celle-ci n'est pas liberee en cas d'interruption
#  VOIR on.exit(gc(verbose=FALSE), add=TRUE)
  

# Input checks
  if (any(K < 0))
    stop("K must be greater or equal to 0")
# (K=0 works but is not advised)
  if (!is.matrix(X))
    stop("X must be a matrix")
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (p<2)
    stop("X must have more than one column")
  if (n<4)
    stop("X must have more than three rows")
  
  if (!is.list(MyFamily))
    stop("MyFamily must be a list")
  if (length(MyFamily) == 0) 
	stop("MyFamily is an empty list")
        
  if ((min.ev < 0) || (min.ev > 0.01))
    stop("bad value of min.ev; must be in [0, 0.01] ")

  # Calculate Dmax: maximum degree of the nodes of the graph
  Dmax <- checkinit(p, n,MyFamily)

  # Calculate the penalty
  pen <- penalty(p, n,  Dmax, K)

  # Centrer les donnees
  X <- scale(X,center=TRUE,scale=FALSE)

  # Main calculations
  res <- mincrit(X,MyFamily,n,p,pen,min.ev, max(Dmax))
  
  # Suppress the last dimension of G and Neighb when equal to 1
  if (ncol(pen) ==1) {
    res$G <- as.matrix(res$G[,,1])
    res$Neighb <- as.matrix(res$Neighb[,,1])
  }
    
  return(res)
} # fin selectMyFam


	
checkinit <- function(p, n, MyFamily){
  # ---------------------------------------------------------------
  # FUNCTION
  #   Verify the arguments of the main function
  # INPUT
  #   p, n : integers
  #   MyFamily : list of pxp adjacency matrices
  # OUTPUT
  #   Dmax: integer vector of dimension p
  # CALLED BY
  #  selectMyFam
  # ---------------------------------------------------------------
	dmax <- 0
	for (g in 1:length(MyFamily)) {
          G <- MyFamily[[g]]

		# qq verifications concernant G
               if (!is.matrix(G))
                 stop( "the ", g, "st family is not a matrix" )
		if (dim(G)[1]!=p || dim(G)[2]!=p) 
			stop(
             "the ", g, "st adjacency matrix is not a pxp matrix. Please provide square pxp matrices, with p=ncol(X)=", p)

		if (!all(G %in% c(0,1))) 
			stop(
             "the ", g, "st adjacency matrix has values different from 0 and 1. Please provide matrices with values in 0 and 1")

		if (!isSymmetric(G)) 
			stop(
             "the ", g, "st adjacency matrix is not symmetric. Please provide symmetric matrices")



		# calcul du degre et degre max
		deg <- max(apply(G,2,sum))
		dmax <- max(dmax,deg)
		} # fin G
        
	if (dmax>=n-2) 
		stop(
     "some graphs in MyFamily have a degree equal or larger than n-2. Please provide graphs with degree smaller than n-2, with n=nrow(X)=",n)

# retourner dmax repete p fois
	return(array(dmax,p))
} # fin checkinit

# ----------------------------------
# Fonction appelee par mincrit via apply
# -------------------------------------
fung2 <- function(G) {
  return (sum(G==1))
} # fin fung2

# ----------------------------------
# Fonction appelee par mincrit via apply
# -------------------------------------
fung1<- function(G, maxnvois) {
  l = grep(1,G)
  
  return ( c(l, rep(0, maxnvois-length(l))))
  } # fin fung1



# ----------------------------------
# Fonction appelee par mincrit via apply
# -------------------------------------
funGGMSCRa<- function(a,
            n, p, X, min.ev, NVois, sumX2,
            GG, scr, iwork, work, svdMd, r1,
            W1, M, W2, W3, W4, vu, svdMv, xvals, Pr) {

# 10/09/2014: remove the option DUP=FALSE (deprecated)  
          res <- .C("GGMSCRa",
               as.integer(a), as.integer(n), as.integer(p),
               as.double(X), as.double(min.ev),
               as.integer(NVois), as.double(sumX2),
               as.integer(GG), scr=as.double(scr),
               as.integer(iwork),as.double(work),
               as.double(svdMd), as.double(r1),
               as.double(W1), as.double(M),
               as.double(W2), as.double(W3), as.double(W4),
		    as.double(vu), as.double(svdMv),
		    as.double(xvals),
		    as.double(Pr))$scr

          return(res[a])
        } # fin funGGMSCRa

  # ----------------------------------


mincrit <- function(X,MyFamily,n,p,pen,min.ev, maxNVois){
  # ---------------------------------------------------------------
  # FUNCTION
  #   Calculate the minimum criterion over all the families
  #   and the corresponding G
  # INPUT
  #   See selectMyFam
  # OUTPUT
  #    crit.min, G: see selectMyFam
  #    ind.min: index of the corresponding family 
  # CALLED BY
  #  selectMyFam
  # ---------------------------------------------------------------
  lK <- ncol(pen)
  # Allocation et initialisation des sorties
  critmin <- rep(Inf, lK)
  indmin <- array(Inf, lK)
  Dmaxmax <- nrow(pen)-1
  Gmin <- array(0,c(p,p, lK))
  Neighbmin <- array(0, c(p,   Dmaxmax, lK))
# Pour l'appel au C
  NVois <- array(0, p)
  work <- array(0, n*maxNVois)
 
  
#CHANGE 5/12/2011 taille iwork <- array(0,p)
#CHANGE 17/04/2012  iwork <- array(0,max(n,p))
  iwork <- array(0, n*p)

  
          svdMd<- array(0,p*p)
          r1<- array(0, n*p)
          W1<- array(0, n*p)
          W2<- array(0, p*p)
          W3<- array(0, p*p)
          M<- array(0, p*p)
          W4<- array(0,n)
          vu<- array(0, p*p)
          svdMv<- array(0, p*p)
          xvals<- array(0, p*p)
          Pr<- array(0,n)
          
          sumX2 <- apply(X**2,2, sum)
          scr <- array(0,p)
          Neighb <- array(0, c(p,   Dmaxmax, lK))
          crit <- rep(Inf, lK)

  # Boucle sur les familles
  for (g in 1:length(MyFamily)) {
          G <- MyFamily[[g]]
	if (!all(diag(G) ==0)) {
			warning(
"the diagonal of the ", g, "st adjacency matrix is supposed to be null")
                        diag(G) <- 0
                      }


           NVois <- apply(G, 2, fung2)
     #     GG <-  matrix(0, ncol=maxNVois, nrow=p)
           GG= matrix(apply(G, 1, fung1, maxNVois), ncol=maxNVois, byrow=TRUE)

            # GG[a,] = les indices des voisins de a

          vectp <- 1:p
           scr = sapply(vectp, funGGMSCRa,
            n, p, X, min.ev, NVois, sumX2,
            GG, scr, iwork, work, svdMd, r1,
            W1, M, W2, W3, W4, vu, svdMv, xvals, Pr)
          
          err <-0
          Neighb[,,] <- 0
          crit[] <- Inf
# 10/09/2014: remove the option DUP=FALSE (deprecated)  
             res <- .C("GGMGrMin",
                    as.integer(n), as.integer(p),
                    as.integer(lK), as.integer(ncol(GG)),
	       as.integer(Dmaxmax), as.double(scr), as.double(pen),
                    as.integer(GG),
	       as.integer(NVois), 
	      crit=as.double(crit),
               Neighb= as.integer(Neighb),
                    err=as.integer(err),  NAOK=TRUE)
          

          if (res$err != 0)
            warning(
 "some values of penalty are greater than 1e+08: the criterion have been set to Inf at least once during the process (family ", g,")")
          res$Neighb <- array(res$Neighb, c(p,   Dmaxmax, lK))
                              
		for (iK in 1:lK) {
			if (res$crit[iK] < critmin[iK]) {
			Gmin[,,iK] <- G
			critmin[iK] <- res$crit[iK]
                        indmin[iK] <- g
                        
                        Neighbmin[,,iK] <- res$Neighb[,,iK]
                      } # fin (res$critmin[iK] <= critmin[iK])
			} # fin iK
		} # fin boucle sur les familles
  
          # Labeller les sorties
          dimnames(Neighbmin) <- list(as.character(1:p),
                           NULL,  dimnames(pen)[[2]])
          dimnames(Gmin) <- list(dimnames(Neighbmin)[[1]],
                                    dimnames(Neighbmin)[[1]],
                                    dimnames(pen)[[2]])
          names(critmin) <- dimnames(pen)[[2]]
           # ---------------------

	return(list(crit.min=critmin, G=Gmin, Neighb=Neighbmin,
                    ind.min=indmin))
	}
	
	
