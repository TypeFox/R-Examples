verifyArg <- function( X, dmax, K, min.ev, max.iter) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Verify the arguments common to the main functions
  # selectFast and selectQE
  # INPUT
  #   See selectFast
  # OUTPUT
  #   A list with components:
  # n: integer scalar. Number of rows of X
  # p:  integer scalar. Number of columns of X
  # Dmax: integer vector of length p.
  #     (Dmax[j] = maximum degree of node j)
  #     si dmax est un scalaire: Dmax=rep(dmax,p)
  #      sinon Dmax=dmax
  # CALLED BY
  #  selectQE, selectFast
  # ---------------------------------------------------------------

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

  if (max.iter<=0)
    stop("bad value of max.iter: should be positive")
  
  # dmax verification
  ldmax <-length(dmax) 
  if ( any(dmax<1) ||  any(dmax > (n-3)) || any(dmax > (p-1)))
#//  if ( any(dmax<1) ||  any(dmax > (n-3)) || any(dmax > p))
    stop(paste(
         "values in dmax must be greater than 0, less or equal than the number of rows of X - 3 (",
               (n-3),
               "), less  or equal than the number of columns of X (",
               (p-1), ")"))
    if ( any(round(dmax) != dmax)) {
      stop("dmax must be integer")
    }
    if ( (ldmax != 1) && (ldmax != p ))
      stop("dmax must be of length 1 or p")

  if (ldmax ==1) Dmax <- rep(dmax,p)
  else
    Dmax <- dmax
  
  # min.ev verification
  if ((min.ev < 0) || (min.ev > 0.01))
    stop("bad value of min.ev: should be in [0, 0.01]")

return(list(n=n, p=p, Dmax=Dmax))
} # fin verifyArg


selectFast <- function(X,
                  dmax=min(floor(nrow(X)/3),nrow(X)-3,ncol(X)-1),
                  K=2.5, family="EW",
                  min.ev=10**(-8),
                  max.iter=200, eps=0.01,
                  beta=nrow(X)*nrow(X)/2,
                  tau=1/sqrt(nrow(X)*(ncol(X)-1)),
                  h=0.001, 
                  T0=10,
                  verbose=FALSE) {
  
  # ---------------------------------------------------------------
  # FUNCTION
  #   Main function for the LA, EW, C01 families
  #    of the package GGMselect
  #   End-user function.
  #   Estimate the graph of partial correlation from data
  # INPUT
  # Same input as selectQE:
  #   X : n x p matrix (data set)
  #   dmax : scalar or p dimensional vector
  #         (maximum degree of the nodes of the graph)
  #         positive integers <=  n-3, p-1
  #         Default value: min(c(3,n-3,p-1))
  #   K : scalar or vector (tuning parameter)
  #   family : scalar or vector taking value in "LA", "EW", "C01"
  #   (methods for building the collection of candidate graphs)
  #   min.ev : positive real number (for inversion)
  # Entrees specifiques a selectFast:
  #   max.iter : positive integer number (maximal number of
  #              iterations)
  #   eps  : positive real number (precision)
  #   beta, tau, h, T0: parameters of ModLasso
  # Other:
  #   verbose: flag for printing intermediary results in real-time
  # OUTPUT
  #   A list with components: EW, LA,C01, C01.LA,  C01.LA.EW
  #     (according to the family argument)
  #    Each one is the result obtained with a family.
  #    It is a list with components:
  #       Neighb:  array p x max(dmax) x length(K)
  #         Neighb[j, , iK ]: nodes connected to j for K[iK]
  #       crit.min: vector of dimension length(K).
  #         It gives the minimal values of the criterion
  #         for each value of K
  #       G: adjacency matrix with dimension pxpx length(K)
  #        Each slice iK (for iK=1 to lK) is the adjacency
  #         matrix of the  "directed" graph for K=K[iK]
  #     (0 = no edge, 1 = an edge)
  # Careful: Neighb and G are matrices if length(K)=1
  # SEEALSO
  #   selectQE selectMyFam
  # ---------------------------------------------------------------
  # Liberer la memoire on exit
  # car celle-ci n'est pas liberee en cas d'interruption
  #  par l'utilisateur
  on.exit(gc(verbose=FALSE), add=TRUE)
  
  # Verifier les arguments qui sont communs a selectQE
  res <- verifyArg( X, dmax, K, min.ev, max.iter)
  # Deplier la liste retournee.
  n <- res$n; p <- res$p;  Dmax <- res$Dmax

 # Verifier les autres arguments
  validmethods <- c("LA", "EW", "C01")
  # Chaque methode ne doit apparaitre qu'une fois
  if (any((unique(family) != family)==TRUE) ||
      any(!is.element(family, validmethods)))
    stop ("bad value of family")
  if ((eps<=0) || (eps>1))
    stop("bad value of eps")
  
  # Centrer les donnees
  X <- scale(X,center=TRUE,scale=FALSE)
  XNorm <- scale(X,center=FALSE)/sqrt(n-1)
  
  # Calculate the penalty
  pen <- penalty(p, n,  Dmax, K)
  
  # Main calculations
  return( calcLarsNEW(X, XNorm, Dmax,  pen, family,
                      min.ev, max.iter, eps,
                      beta, tau,h, T0,
                      verbose))
}




calcLarsNEW <- function(X, XNorm, Dmax, pen, family,
                        min.ev,
                        max.iter, eps,
                        beta, tau, h,  T0,
                        verbose) {

  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute the estimated graph for each family
  # INPUT
  #   X :  n x p data  matrix
  #   XNorm: n x p data matrix after rescaling (columns of norm 1)
  #   Dmax : p dimensional vector
  #     (Dmax[j] = maximum degree of node j)
  #   pen: max(Dmax) x lK array (penalty)
  #   family: vector of values in "EW", "LA", "C01"
  #   min.ev: (small) positive real number
  #   max.iter : positive integer number (maximal number of
  #              iterations)
  #   eps  : positive real number (precision)
  #   other: see selectFast
  # OUTPUT
  #   same as selectFast
  # CALLED BY selectFast
  # ---------------------------------------------------------------

  p <- dim(X)[2]
  n <- dim(X)[1]
  lK <- dim(pen)[2]

  # Initialisations
  output <- NULL
  scr.init <- apply(X**2,2,sum) # residual sum of square (no edge)
  crit.min.init <- rep(sum(scr.init),lK) # selection criterion (no edge)
  names(crit.min.init) <- dimnames(pen)[[2]]
  Dmaxmax <- max(Dmax) # Degree of the graph
# liste des voisins de chaque sommet, pour chaque valeur de K:

  Neighb <- array(0,c(p, Dmaxmax, lK)) 
  dimnames(Neighb) <- list(as.character(1:p),
                           NULL,  dimnames(pen)[[2]])
  # liste des voisins de chaque sommet, pour le graphe courant:
  NVoisGraph <- array(0,p)
  Graph<- array(0,c(p,Dmaxmax))  # current graph

  # Structures de travail, pour l'appel au C:
 #CHANGE 5/12/2011  iwork <- array(0,p)
#CHANGE 5/09/2012  iwork <- array(0,max(n,p))
    iwork <- array(0,n*p)
  work <- array(0, n*Dmaxmax)
  
  svdMd<- array(0,p)
   r1<- array(0, n*p)
   W1<- array(0, n*p)
  W2<- array(0, p*p)
  W3<- array(0, p*p)
  M<- array(0, p*p)
   W4<- array(0,Dmaxmax )
    vu<- array(0, p*p)
    svdMv<- array(0, p*p)
    xvals<- array(0, p*p)
   Pr<- array(0,n)

  
  if ( (any(family == "LA")) || (any(family == "EW")))  {
    max.st <- min(n,p-1)
    Gr <- array(0,c(p,max.st))  # current memory of positive actions
    NVoisGr <- array(0,p)
  }
  

  
  for (afamily in family)
    {
      if (verbose==TRUE)
        cat("*** Run family", afamily,"\n")
      # Reinitialisation
      Neighb[, ,] <- 0
        NVoisGraph[] <- 0
        Graph[,] <- 0

      scr <- scr.init
      crit.min <- crit.min.init
      if ( (afamily == "LA") || (afamily == "EW"))  {
        Gr[,]  <- 0
        NVoisGr[] <- 0
        GrGlob  <- calcModLasso(XNorm, Dmax, afamily, max.st,
                                max.iter, eps, beta, tau,h, T0, verbose)

      } 
      switch(afamily,
             LA =
             {
           # Collection "LA"
res <- .C("GGMloopAND",
   as.integer(n), as.integer(p), as.integer(lK),
   as.integer(nrow(GrGlob)), as.integer(ncol(GrGlob)),
as.integer(GrGlob), as.integer(Dmax),
   as.double(min.ev), as.double(X),  as.double(scr.init),
   as.double(pen), as.integer(ncol(Gr)), as.integer(ncol(Graph)),
as.integer(NVoisGraph), as.integer(NVoisGr),
	       as.integer(Graph), as.integer(Gr), as.integer(Dmaxmax),
	       as.double(scr),as.integer(iwork), as.double(work),
	       as.double(svdMd), as.double(r1),
		    as.double(W1), as.double(M),
		    as.double(W2), as.double(W3), as.double(W4),
		    as.double(vu), as.double(svdMv),
		    as.double(xvals),
	       as.double(Pr),
	       crit.min=as.double(crit.min),
             Neighb=as.integer(Neighb),
          NAOK=TRUE)
names(res$crit.min) <-   names(crit.min)      
res$Neighb <- array(res$Neighb,c(p, Dmaxmax, lK))
dimnames(res$Neighb) <- dimnames(Neighb)
               output$LA$Neighb <-res$Neighb
               output$LA$crit.min <- res$crit.min
               output$LA$G <- convGraph(res$Neighb)
             }, #end family=LA
             
             EW =
             {

res <- .C("GGMloopEWOR",
   as.integer(n), as.integer(p), as.integer(lK),
   as.integer(nrow(GrGlob)), as.integer(ncol(GrGlob)),
as.integer(GrGlob), as.integer(Dmax),
   as.double(min.ev), as.double(X),  as.double(scr.init),
   as.double(pen), as.integer(ncol(Gr)), as.integer(ncol(Graph)),
as.integer(NVoisGraph), as.integer(NVoisGr),
	       as.integer(Graph), as.integer(Gr), as.integer(Dmaxmax),
	       as.double(scr),as.integer(iwork), as.double(work),
	       as.double(svdMd), as.double(r1),
		    as.double(W1), as.double(M),
		    as.double(W2), as.double(W3), as.double(W4),
		    as.double(vu), as.double(svdMv),
		    as.double(xvals),
	       as.double(Pr),
	       crit.min=as.double(crit.min),
             Neighb=as.integer(Neighb),
          NAOK=TRUE)

names(res$crit.min) <-   names(crit.min)      
res$Neighb <- array(res$Neighb,c(p, Dmaxmax, lK))
dimnames(res$Neighb) <- dimnames(Neighb)
               output$EW$Neighb <- res$Neighb
               output$EW$crit.min <- res$crit.min
               output$EW$G <- convGraph(res$Neighb)
             }, # end family EW
             C01 =
             {
               GrGlobC01  <- calcModC01(XNorm)
               res <- .C("GGMloopC01",
   as.integer(n), as.integer(p), as.integer(lK),
   as.integer(nrow(GrGlobC01)), as.integer(ncol(GrGlobC01)),
as.integer(GrGlobC01), as.integer(Dmax),
   as.double(min.ev), as.double(X),  as.double(scr.init),
   as.double(pen),  as.integer(ncol(Graph)),
as.integer(NVoisGraph),  as.integer(Graph),
           as.integer(Dmaxmax),
	       as.double(scr),as.integer(iwork), as.double(work),
	       as.double(svdMd), as.double(r1),
		    as.double(W1), as.double(M),
		    as.double(W2), as.double(W3), as.double(W4),
		    as.double(vu), as.double(svdMv),
		    as.double(xvals),
	       as.double(Pr),
	       crit.min=as.double(crit.min),
             Neighb=as.integer(Neighb),
          NAOK=TRUE)
               
               names(res$crit.min) <-   names(crit.min)      
               res$Neighb <- array(res$Neighb,c(p, Dmaxmax, lK))
               dimnames(res$Neighb) <- dimnames(Neighb)
               output$C01$Neighb <- res$Neighb
               output$C01$crit.min <- res$crit.min
               output$C01$G <- convGraph(res$Neighb)
             }, # end family C01
             stop("Internal error : bad value of 'family'")
             ) # end of switch

      
      if (verbose==TRUE)
      cat("*** End family", afamily,"\n")
    } # end of loop afamily
  
  # MIXTE: LA+C01
  if (is.element("C01",family) && is.element("LA",family)) {
    crit <- cbind(output$LA$crit.min,output$C01$crit.min)
    # Pour chaque valeur de K, on prend le plus petit
    # critere parmi celui calcule par LA et celui calcule par C01
    ind.min <- apply(crit,1,which.min)
    output$C01.LA$crit.min <- 0*output$LA$crit.min
    output$C01.LA$Neighb <- 0*output$C01$Neighb
    for (iK in 1:lK) {
      output$C01.LA$crit.min[iK] <- crit[iK,ind.min[iK]]

      switch(ind.min[iK],
             # =1 (le 1ier critere est celui calcule par LA)
             output$C01.LA$Neighb[,,iK] <- output$LA$Neighb[,,iK],
             # =2 (le 2ieme critere est celui calcule par C01)
             output$C01.LA$Neighb[,,iK] <- output$C01$Neighb[,,iK],
             stop("Internal error: bad value ind.min")
             ) # fin switch
    }

    # Calculer la matrice d'adjacence a partir du tableau
    # des voisins
    output$C01.LA$G <- convGraph(output$C01.LA$Neighb)
  }
  #MIXT LA+C01+EW
  if
  (is.element("C01",family) && is.element("LA",family) && is.element("EW",family)) {
#    output$C01.LA.EW <- list(NULL)
    crit <- cbind(output$LA$crit.min,output$C01$crit.min,output$EW$crit.min)
    ind.min <- apply(crit,1,which.min)
    output$C01.LA.EW$crit.min <- 0*output$LA$crit.min
    output$C01.LA.EW$Neighb <- 0*output$LA$Neighb
    for (iK in 1:lK) {
      output$C01.LA.EW$crit.min[iK] <- crit[iK,ind.min[iK]]

      switch(ind.min[iK],
             #=1
             output$C01.LA.EW$Neighb[,,iK] <- output$LA$Neighb[,,iK],
             # =2
             output$C01.LA.EW$Neighb[,,iK] <- output$C01$Neighb[,,iK],
             # =3
             output$C01.LA.EW$Neighb[,,iK] <- output$EW$Neighb[,,iK],
             stop("Internal error: bad value ind.min")
             ) # fin switch
    }
    output$C01.LA.EW$G <- convGraph(output$C01.LA.EW$Neighb)
  }
# Suppress the last dimension of G and Neighb when equal to 1
  if (lK ==1)
    output <- simplifDim(output)

  
  return(output)
} # fin calcLarsNEW

simplifDim <- function(output) {
  # ---------------------------------------------------------------
  # FUNCTION  
  # Simplifier la structure des sorties:
  # enlever la derniere dimension de G et Neighb si elle
  # est egale a 1
  # INPUT
  #   output: liste contenant tous les resultats finaux
  # OUTPUT
  #   output: all the components G and Neighb have one dimension less
  # INPUT CONDITION
  #   The last dimension of G and Neighb is equal to 1
  # CALLED BY
  # calcLarsNEW calcLarsNEWQE
  # ---------------------------------------------------------------
for (meth in names(output)) {
  output[[meth]][["G"]] <- as.matrix(output[[meth]][["G"]][,,1])
  output[[meth]][["Neighb"]] <- as.matrix(output[[meth]][["Neighb"]][,,1])
}
return(output)
}

  

