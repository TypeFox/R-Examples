# VERSION 2009
# Copyright@INRA-2009
# Programs for the QE family


selectQE <- function(X, 
                  dmax=min(3,nrow(X)-3,ncol(X)-1),
                  K=2.5,   
                  min.ev=10**(-8),
                  max.iter=10**6,
                  max.nG=10**8,
                  max.size= 10**8,
                  verbose=FALSE) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Main function for the QE family of the package GGMselect
  #   End-user function.
  #   Estimate the graph of partial correlation from data
  # INPUT
  # Entrees identiques a celles de selectFast:
  #   X : n x p matrix (data set)
  #   dmax : scalar or p dimensional vector
  #         (maximum degree of the nodes of the graph)
  #         positive integers <=  n-3, p-1
  #         Default value: min(c(3,n-3,p-1))
  #   K : scalar or vector (tuning parameter)
  #   min.ev : positive real number (for inversion)
  # Entrees specifiques a selectQE:
  #   max.iter: integer scalar = maximum number of stepwise
  #             iterations
  #   max.nG: integer scalar= maximum number of graphs
  #     per collection. Stepwise if more.
  #    When the size of Mod (matrix of the models) is greater
  #    than this value,  a stepwise procedure is run.
  #   max.size: integer scalar = maximum threshold of the memory
  #          required to calculate SCR. The calculations
  #          are trigged only if the total dimension of the
  #          structures SCR and SCRmin do not exceed this value.
  #          Execution is abandoned, otherwise.
  #          (see calcSCRQE)
  # Other:
  #   verbose: logical, TRUE if intermediary output is required
  #         (to trace the current process in real time)
  # OUTPUT
  #   A list with components:
  #     G, Neighb, crit.min: see selectFast
  # SEEALSO
  #   selectFast selectMyFam
  # --------------------------------------------------------
  # PS:
  # Each slice of G can be an input for the "ggm" package
  # Example: require(ggm); drawGraph(output$QE$G[,,1])
  # ---------------------------------------------------------------
  # Liberer la memoire on exit
  # car celle-ci n'est pas liberee en cas d'interruption
  # de l'execution par l'utilisateur
  on.exit(gc(verbose=FALSE), add=TRUE)
  
  # Verifier les arguments qui sont communs aux autres methodes
  # (see selectFast) et calculer Dmax, le degre de chaque noeud
  res <- verifyArg( X, dmax, K, min.ev, max.iter)
  # Deplier la liste retournee= extract.named(res)
  n <- res$n; p <- res$p;  Dmax <- res$Dmax

  # Rescaling
  X <- scale(X,center=TRUE,scale=FALSE)
  XNorm <- scale(X,center=FALSE)/sqrt(n-1)

  # Calculate the penalty
  pen <- penalty(p, n,  Dmax, K)

  # Main calculations
  retour <- calcLarsNEWQE(X, XNorm, Dmax,  pen,
                      min.ev,
                      max.nG, max.iter, max.size,  verbose)
  return(retour$QE)
} # fin selectQE




calcLarsNEWQE <- function(X, XNorm, Dmax, pen,
                        min.ev,
                        max.nG, max.iter,
                        max.size,  verbose) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Compute the estimated graph for QE family
  # INPUT
  #   X :  n x p data  matrix
  #   XNorm: n x p data matrix after rescaling (columns of norm 1)
  #   Dmax : p dimensional vector
  #     (Dmax[j] = maximum degree of node j)
  #   pen: max(Dmax) x lK array (penalty)
  #   min.ev, max.nG, max.size, verbose: see selectQE
  # OUTPUT
  #   See selectQE
  # CALLED BY
  #   selectQE
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
  Neighb <- array(0,c(p, Dmaxmax, lK))
  dimnames(Neighb) <- list(as.character(1:p),
                           NULL,  dimnames(pen)[[2]])
  
  if (verbose==TRUE)
        cat("*** Run family QE\n")
  scr <- scr.init
  crit.min <- crit.min.init

  # Call to calcSCRQE;
  # Calculate SCR (residual sums of squares for models
  # with each one of the Dmax dimensions),
  # and Mat.Chap (adjacency  matrix of the current graph for
  # each value of K) and mod.Chap (the chosen models)
  SCRQE.out <- calcSCRQE(X, Dmax, Dmaxmax, min.ev,
                                      max.size,
                                      scr.init, pen, lK,   
                                      verbose)

  # Call to calcGrChapSymQE;
  # Calculate Mat.Chap (adjacency  matrix of the "directed" graph for
  # each value of K) and the minimum criteria
  retSymQE <-  calcGrChapSymQE(SCRQE.out$Mat.Chap,
                                  SCRQE.out$SCR,
                                  scr.init,
                                  pen, Dmax,n, p,
                                  max.nG, max.iter,  verbose)
  output$QE$G <-  retSymQE$Mat.Chap
  output$QE$crit.min <- retSymQE$crit.min

  # Calculer les voisins sous forme d'un tableau d'indices
  # a partir de la matrice d'adjacence
  output$QE$Neighb <- GtoNeighb(output$QE$G, Dmaxmax)

      # Suppress the last dimension of G and neighb when equal to 1
  if (lK ==1)
    output <- simplifDim(output)

  if (verbose==TRUE)
      cat("*** End family QE\n")
  return(output)
} # fin calcLarsNEWQE

verifTaille  <- function( p, Dmax, verbose) {
    #-------------------------------------------------------------
  # FUNCTION
  #  Calculate the total number of elements in SCR
  # INPUT
  #   p : number of observations
  #   Dmax :  integer vector with dimension p. Maximum degrees
  #   verbose: flag for printing results
  # OUTPUT 
  #    number of elements in SCR 
  # CALLED BY
  #  calcSCRQE
  #-------------------------------------------------------------
  dimSCR <- 0
  for (j in 1:p) {
    for (d in 1:Dmax[j]) {
      dimSCR <- dimSCR + choose(p-1,d)
    } # fin d
  } # fin j
  if (verbose)  
    cat("Number of RSS calculated:", dimSCR, "\n")
  return(dimSCR)
}


  


calcSCRQE <- function(X, Dmax, Dmaxmax,  min.ev,
                      max.size, NormX, pen, lK,
                      verbose)  {
  #-------------------------------------------------------------
  # FUNCTION
  #  Compute SCR (residual sums of squares),
  #          Mat.Chap (current adjacency  matrix), and
  #          mod.Chap (the chosen models)
  # INPUT
  #   X : matrix nxp of observations
  #   Dmax :  integer vector with dimension p. Maximum degrees
  #   Dmaxmax : integer scalar = max(Dmax)
  #   min.ev : positive real number (for inversion)
  #   max.size,  etc, ..., verbose: see selectQE
  # OUTPUT : list
  #   SCR: list of p components.
  #      Each component j is a list of Dmax[j] components.
  #       Each component d is a vector of length nMod[d]
  #   SCR[[j]][[d]] : residual sums of squares for models with
  #                          dimension d
  #   Mat.Chap: array (lK, p, p). Adjacency  matrix of
  #    the "directed" graph for K=K[iK]
  #     (0 = no edge, 1 = an edge)
  #   mod.Chap: list of lK components.
  #      Each component is a list of p components.
  #       Each component j is a vector of length Dmax[j]
  #   mod.Chap[[iK]][[j]] the chosen models for variable j
  # NOTE: A VOIR: on ne se sert pas de  mod.Chap
  # CALLED BY
  #  calcLarsNEWQE
  #-------------------------------------------------------------
  p <- dim(X)[2]
  n <- nrow(X)
  if (!is.null(max.size)) {
    # verifier que la taille de SCR n'est pas trop grosse
    sumMem <- verifTaille(p,  Dmax, verbose)
    if (sumMem > max.size) 
      stop("the number of RSS to be calculated, ",
                 sumMem,
               ", is greater than the maximum bound, max.size, i.e ", max.size)
  } # fin (!is.null(max.size))

  
  
  nDmaxmax <- n*Dmaxmax
  DmaxmaxDmaxmax <-  Dmaxmax*Dmaxmax
    
  lesMod <- list(NULL)
  nMod <- rep(0,Dmaxmax)

  for (d in 1:Dmaxmax) {
      nMod[d] <- choose(p-1,d)
      lesMod[[d]] <- t(combinations(p-1,d)) 
    }
    # Structures allocation for the C programme
  
  SCR <- list(NULL)
  SCRminmod <- list(NULL)
  for (d in 1:Dmaxmax) {
    SCRminmod[[d]] <- rep(0, nrow(lesMod[[d]]))
      }
  SCRminSCR <- rep(0.0, Dmaxmax)
  
  for (j in 1:p) {
      SCR[[j]] <- list(NULL)
        for (d in 1:Dmax[j]) {
        SCR[[j]][[d]] <- rep(0.0,nMod[d])
        }
    }


#  // Pour le calcul de modChap et matChap
  mod.Chap <- list(NULL)
  Mat.Chap <- array(0,c(lK,p,p))
  dimnames(Mat.Chap) <- list(dimnames(pen)[[2]],
                                as.character(1:p),
                             as.character(1:p))
  for (iK in 1:lK) {
      mod.Chap[[iK]] <- list(NULL)
      for (j in 1:p) {
        mod.Chap[[iK]][[j]] <- rep(-1,Dmax[j])
      }
    }
  

# The input and working structures for GGMcalcSCRQE
  z <- list(  X=as.double(X),
              Dmax= as.integer(Dmax), nMod=as.integer(nMod),
              lesMod=lesMod,
              n= as.integer(n), p=as.integer(p),
              minvp=  as.double(min.ev),
              M= as.double(array(0, DmaxmaxDmaxmax)),
              W1= as.double(array(0, nDmaxmax)),
              W2= as.double(array(0, DmaxmaxDmaxmax)),
              W3= as.double(array(0, DmaxmaxDmaxmax)),
              svdMd= as.double(array(0, Dmaxmax)),
              vu= as.double(array(0, DmaxmaxDmaxmax)),
              svdMv= as.double(array(0, DmaxmaxDmaxmax)),
              iwork= as.integer(array(0, Dmaxmax)),
              xvals= as.double(array(0, DmaxmaxDmaxmax)),
              W4= as.double(array(0, n)),
              r1= as.double(array(0, nDmaxmax)),
              Proj= as.double(array(0.0, n)),
              workj= as.double(array(0, n)),
              work= as.double(array(0, nDmaxmax)))

# The output list of GGMcalcSCRQE
  z$SCRout <- list(SCR=SCR,Mat.Chap= Mat.Chap, mod.Chap=mod.Chap)


  
    # ENORME BOUCLE (complexite = p**(Dmax+1))
    SCRQE.out <- .Call("GGMcalcSCRQE", NormX,
                       pen , nrow(pen),
                      lK, SCRminSCR, SCRminmod, z)
  
# Pour des raisons pratiques, mod.Chap avait ete initialise a -1:
# il faut oter ces -1
  SCRQE.out$mod.Chap <- rapply(SCRQE.out$mod.Chap,
                               function(X) {
    if (all(X==-1)) return(0) else return(X[X!=-1])
  }, how = "replace" )
                       
     # Garbage collection
  rm(Mat.Chap)
  rm(mod.Chap)
  rm(z)
  gc(verbose=FALSE)

  return(SCRQE.out)
  } # fin calcSCRQE



calcGrChapSymQE <- function( MatG, SCR, NormX, pen,  Dmax, n, p, 
                            max.nG=10**8, max.iter=10**6,
                            verbose=FALSE){
  #-------------------------------------------------------------
  # FUNCTION
  # Compute Mat.Chap (adjacency  matrix of the "directed" graph )
  # from the current matrix MatG 
  # and crit.min (minimum criterions)
  # INPUT
  #   MatG :  output from calcSCRQE
  #   SCR :  output from calcSCRQE
  #   NormX : vector px1
  #   pen : output from penalty 
  #   Dmax : vector with dimension p. Maximum degrees 
  #   n,p : integers
  #   max.nG,...   verbose : see selectQE
  # OUTPUT
  #  List with components:
  #   Mat.Chap :
  #       array ( p, p, lK)
  #        each slice iK (for iK=1 to lK) is the adjacency
  #         matrix of the  "directed" graph for K=K[iK]
  #     (0 = no edge, 1 = an edge)
  #   crit.min : vector with dimension lK. 
  # CALLED BY
  #  calcLarsNEWQE
  #-------------------------------------------------------------
  lK <- dim(pen)[2]
  crit.out <- rep(0,lK)
# Initialisation de la liste en argument des programmes C
  z <- list(n=as.integer(n),
      p=as.integer(length(NormX)),
      penrows= as.integer(nrow(pen)),
      pen=as.double(pen),
      NormX=as.double(NormX),
      Dmax=as.integer(Dmax),
      SCR=SCR,
      G = as.integer(rep  (0,p*p)),
      ind=as.integer(rep  (0,p)),
      scrG=list(scr=rep(0.0,p), d=as.integer(rep(0,p))),
      sumcrit=as.double(0.0))
  W2 <- W3 <- NULL
  Mat.Chap <- MatG*0
  if (verbose)
          cat("***    Running loop GrSymQE from iK = 1 to", lK, "\n")
  
  for (iK in 1:lK) {
    SW <- 0 # 1 if stepwise should be done
    lMod <- 0 # lgueur max de Mod (pour les messages, si verbose)
    z$iK <- as.integer(iK)
    MatG.et <- MatG[iK,,]*t(MatG[iK,,])
    MatG.ou <- pmin(MatG[iK,,]+t(MatG[iK,,]),1)
    Gdiff <- MatG.ou-MatG.et
    if (sum(Gdiff)==0) {
      Mat.Chap[iK,,] <- MatG.et
      z$G <- as.integer(MatG.et)
      crit.out[iK] <- .Call("GGMscrgcritQE", z)$sumcrit
    }    else {
      vect <- Gdiff[lower.tri(Gdiff)]
      ind1 <- (1:length(vect))[vect==1]
      d1 <- length(ind1)
      crit.min <- rep(Inf,d1)
      mod.min <- list(NULL)
      mm <- 0
      Mod <- NULL

      if (verbose)
          cat("***    Running loop GrSymQE for iK =", iK,
              "and d = 1 to", d1, "\n")
      
# dStop : compteur du nombre de pas effectues
      dStop <- 0
      for (d in 1:d1) {
        if (verbose)
          cat("***    GrSymQE: iK =", iK, "d =", d, "\n")

# Le try permet de recuperer les eventuels pbes de place memoire
        # et d'enchainer sur du stepwise s'il y en a
        go <- tryCatch({
         # calcul de Mod
        if (mm==0) {
          rm(Mod); gc(verbose=FALSE) # faire de la place en memoire
          Mod <- t(combinations(d1,d,ind1))
        } else {
          # Quand mm!=0, on est sur que cr existe:
          # c'est la sortie de calcCritminQE
          if (cr$nModTG > 0) {
            if (dim(Mod)[2]==cr$nModTG) break # sortie de la boucle d
            Mod <- as.matrix(Mod[,-cr$ModTG, drop=F])
          }
          j1 <- Mod[d-1,]<ind1[d1]
          sj1 <- sum(j1)
          if (sj1==0) break
# s'arreter dans le cas sj1=1 et d=2          
          if ((sj1==1)  && (d==2)) break
          if ((sj1>0)  && (dim(Mod)[1]==1)) {
            Mod <- t(combinations(sj1,2,Mod[,j1]))
            sj1 <- 0
          }
          if ((sj1>0)  && (dim(Mod)[1]>1)) {
# si Mod est un vecteur le transformer en matrice avec 1 colonne
            Mod <- as.matrix(Mod[,j1])
            # Calcul du nbre de colonnes de la matrice Mod
            # en sortie car il faut l'allouer avant l'appel au C
            # On n'utilise pas  apply car cela amene a
            # dupliquer la matrice: c'est trop gros
             dd1 <- d-1
             ncolModOut <- 0
             for (aa in 1:dim(Mod)[2]) {
               for (bb in 1: length(ind1)) {
                 if (ind1[bb] >Mod[dd1,aa]) ncolModOut <- ncolModOut + 1
               }
             }
            ModOut <- matrix(0, nrow=d, ncol=ncolModOut)
            Mod<- .C("GGMloopGrSymQE",
               as.integer(Mod), as.integer(d),
               as.integer(d-1), as.integer(ncol(Mod)), 
	    as.integer(d), as.integer(ncolModOut),
		    as.integer(length(ind1)),
               as.integer(ind1), ModOut=as.integer(ModOut))$ModOut
            rm(ModOut); gc(verbose=FALSE)# faire de la place en memoire
            Mod<-matrix(Mod, nrow=d)
          } # fin ((sj1>0)  && (dim(Mod)[1]==1))
        }  # fin  (mm==1)
# Fin du calcul de Mod
        # calcCritminQE calculates min (minimum value of the
        # criteria over all models in Mod) and
        # argmin (models for which the minimum is attained)
        # and ModTG (all the models)
        cr <- calcCritminQE(MatG.et,Mod, z)
      }, error=function(e){cat("\ncalcGrSymQE:\n");print(e);
                           cat("Careful: Iterations stop with iK=", iK,
                               "and d=",d,
              "because of memory problems when allocating matrix with",
                               d, "rows and",
          ncolModOut, "columns. Stepwise procedure\n")})

        # Des pbes memoire sont survenus: on continue sur du stepwise
if (!is.list(go)) {
  warning(" *** Memory problem: => stepwise is run for iK=",
          iK)
  SW <- 1
  break
  #Sortie de la boucle d
}

        dStop <- dStop+1
        if ((mm==0)&&(cr$nModTG != 0)) mm <- 1
        crit.min[d] <- cr$min
        mod.min[[d]] <- cr$argmin
        # lMod: lgueur max de Mod (pour les messages, si verbose)
        lMod <- max(lMod,length(Mod))

        if (length(Mod)>max.nG) {
        warning(" *** Run stepwise procedure for iK=",
                      iK,
              " because the number of graphs in the collection =",
                      length(Mod), " (>max.nG=", max.nG,")")
          SW <- 1
          break
        }

      } # fin d
      rm(Mod); gc(verbose=FALSE)
      crit.chap <- min(crit.min)

      z$G <- as.integer(MatG.et)
      crit.et <- .Call("GGMscrgcritQE", z)$sumcrit
  # calcul de Mat.Chap et crit.out : comparaison Graphe.et avec le
#      minimum entre et/ou                  
      if (crit.et < crit.chap) {
        Mat.Chap[iK,,] <- MatG.et
        crit.out[iK] <- crit.et
      }      else {
        d <- which.min(crit.min)
        vect.chap <- rep(0,length(vect))
        vect.chap[ mod.min[[d]] ] <- 1
        MatChap <- Mat.Chap[iK,,]
        MatChap[lower.tri(MatChap)] <- vect.chap
        Mat.Chap[iK,,] <- MatG.et+MatChap+t(MatChap)
        crit.out[iK] <- crit.chap
      }

      # Stepwise procedure
      if (SW==1) {
        if (verbose==1)
          cat(" *** Stepwise procedure for iK =",
                      iK, "\n")
          if (is.null(W2)) {
          # Les tableaux de travail n'ont pas encore ete alloues
          W2=rep(0, p*p)
          W3=rep(0, p*p)
        }
        # Calcul du graphe correspondant au dernier graphe calcule
#       pour d=dStop
        vect.chap <- rep(0,length(vect))
        vect.chap[ mod.min[[dStop]] ] <- 1
        MatChap <- 0*Mat.Chap[iK,,]
        MatChap[lower.tri(MatChap)] <- vect.chap
        Mat.Stop <- MatG.et+MatChap+t(MatChap)
        crit.Stop <- crit.min[dStop]
        retsw <-
          calcSW(Mat.Stop, MatG.ou,crit.Stop,
                   max.iter, z, W2, W3)
                # comparaison de retsw$critmin avec crit.out[iK]
        if (retsw$critmin < crit.out[iK]) {
          Mat.Chap[iK,,] <- retsw$MatChap
          crit.out[iK] <- retsw$critmin
        }
      }       else {
        if (verbose==1)
          cat(" *** No stepwise procedure for iK =",
            iK,
          " (the number of graphs in the collection =", lMod,")\n")
      }
          
    } # fin (sum(Gdiff)!=0)
  } # fin iK
  return(list(Mat.Chap=aperm(Mat.Chap, c(3,2,1)),
              crit.min= crit.out))
} # fin calcGrChapSymQE 



calcCritminQE <- function(MatGetiK, Mod, z) {
  #-------------------------------------------------------------
  # FUNCTION
  # Compute the criterion and update minimum if needed
  # INPUT
  # MatGetiK : adjacency pxp matrix of graph "LA" for a value iK
  # Mod : matrix with d rows containing models with dimension d
  # INPUT/OUTPUT
  # z: list (see calcGrSymChap)
  # OUTPUT list with components
  # min : minimum value of the criteria over all models in Mod
  # argmin : vector with dimension d.
  #        Models for which the minimum is attained
  # ModTG: integer vector of the models
  #       with dimension=the number of models
  # nModTG: integer scalar. Effective length of ModTG
  # CALLED BY
  #   calcGrChapSymQE
  #-------------------------------------------------------------

  ll <- dim(Mod)[2]   # nombre de modeles
  dd <- dim(Mod)[1]
  nModTG <- 0
  critmin <- Inf;
  
  z$ll <-  as.integer(ll)
  z$dd <-  as.integer(dd)
  z$Mod  <- as.integer(Mod)
  z$critmin <- as.double(critmin)
  z$critargmin <- as.integer(rep(0,dd))
  z$ModTG <- as.integer(rep(0,ll))
  z$nModTG <- as.integer(nModTG)
  z$MatGetiK <- as.integer(MatGetiK)

  result <- .Call("GGMcritminQE", z)
 # en sortie, le composant G contient la derniere mat G calculee
  if (result$nModTG ==0) {
    result$ModTG <- NULL
  }  else {
    result$ModTG <- result$ModTG[1:result$nModTG]
  }
  
  return(
         list(min=result$critmin, argmin=result$critargmin,
              nModTG=result$nModTG, ModTG=result$ModTG))
} # fin fonction
      
  
calcSW <- function(Mat.min ,Mat.max, crit.min,
                   max.iter, z, W2, W3)
  {
  #-------------------------------------------------------------
  # FUNCTION
  # Stepwise procedure
  # INPUT
  # Mat.min : adjacency pxp matrix of the current graph
  # Mat.max : pxp matrix
  # crit.min : minimum criterion
  # max.iter: see selectQE
  # z: list with the other required input    
  # WORKING
  #   W2, W3 : working array, of length >= (p*(p+1))/2 -p
  # OUTPUT
  #   Mat.Chap: adjacency pxp matrix
  #   crit.min: scalar
  # CALLED BY
  #   calcGrChapSymQE
  #-------------------------------------------------------------
    
    p <- nrow(Mat.max)
    Mat.Chap <- Mat.min
    
    # APPEL C
    ret <-.Call("GGMbcSW",
          as.integer(Mat.max), as.integer(Mat.min),
          as.integer(p), as.integer(max.iter),
                list( MatChap=as.integer(Mat.Chap),
                     critmin=as.double(crit.min)),
          z,
          as.integer(W2), as.integer(W3))
    return(ret)
  }
 



GtoNeighb <- function(AdjG, Dmaxmax) {
  # ---------------------------------------------------------------
  # FUNCTION
  #   Convert  adjacency matrices into graphs
  # INPUT
  #   AdjG: array p x p x length(K)
  #   Dmaxmax: second dimension of Neighb
  # OUTPUT
  #   Neighb:  array p x max(dmax) x length(K)
  # CALLED BY
  # calcLarsNEWQE
  # ---------------------------------------------------------------
  lK <- dim(AdjG)[3]
  p <- dim(AdjG)[1]
  Neighb <- array(0,c(p, Dmaxmax, lK))
  for (iK in 1:lK) {
    for (a in 1:p) {
      liste <- grep(1, AdjG[a,, iK])
      if (length(liste) > 0)
        Neighb[a,1:length(liste), iK] <- liste
    }
  }
  
  dimnames(Neighb) <- list(dimnames(AdjG)[[1]],
                           NULL,
                         dimnames(AdjG)[[3]])
return(Neighb)
} # fin GtoNeighb



