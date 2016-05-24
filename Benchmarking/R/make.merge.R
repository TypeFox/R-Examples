# $Id: make.merge.R 140 2015-05-15 21:48:02Z B002961 $


make.merge <- function(grp, nFirm=NULL, X=NULL, names=NULL)  {
   # Opstiller aggregeringsmatrix for at danne grupperne grp ud fra X.
   # Hvad der skal merges skal angives som indeks i en liste af arrays
   # hvor hvert array er indeks for de enheder der skal indgaa i en given
   # gruppe
   if ( class(grp) == "factor" )  {
      # print("Faktor")
      g <- nlevels(grp)
      K <- Kg <- length(grp)
      Kn <- -1
   } else if ( class(grp)=="list" && class(grp[[1]])=="character" ) {
      # print("Liste af navne")
      g <- length(grp)
      Kn <- K <- length(names) 
      Kg <- K
   } else {
     # print("Liste af numre")
     g <- length(grp)
     Kg <- -1
   }
   if ( !is.null(nFirm) && class(nFirm)!="numeric" && 
                                class(nFirm)!="integer" )
      stop("The argument nFirm must be numeric or integer")
   if ( !is.null(X) && class(X)!="matrix" )
      stop("The argument X must be a matrix")
   # print(g)
   if ( Kg == -1 & is.null(X) & is.null(nFirm) ) {
      stop("Either X or nFirm must be in the call to merge.matrix or grp must be a factor")
   }
   Kx <- -1
   if ( !is.null(X) ) {
      K <- Kx <- dim(X)[1]
   }
   if ( !is.null(nFirm) )
      K <- nFirm
   if ( !is.null(names) )
      Kn <- length(names) 
   if ( !is.null(nFirm) && !is.null(X) && Kx != K )
      stop("nFirm must be the number of rows in X")
   if (Kg!=-1 && !is.null(nFirm) && Kg!=K )
      stop("nFirm must be the length of the facotr grp")
   if (Kg!=-1 && !is.null(X) && Kg!=Kx )
      stop("The length of the factor grp must be the number of rows in X")
   if ( !is.null(names) && K>0 && K!=Kn )
      stop("The length of names must be the number of firms")
   if ( class(grp)=="list" && class(grp[[1]])=="character" && Kn <= 0)
      stop("When grp is a list of names for mergers the argument names must also be supplied")
   if ( K < 0 && Kn > 0 )
      K <- Kn

   Mer <- matrix(0, nrow=g, ncol=K)
   if ( class(grp) == "factor" )  {
      for ( i in 1:g )  {  # Saet 1-taller soejler for dem der skal merges
     	    Mer[i,as.numeric(grp)==i] <- 1 
      }
   } else if ( class(grp)=="list" && class(grp[[1]])=="character")  {
      for ( i in 1:g )  {
         Mer[i,which(names %in% grp[[i]])] <- 1
      }
   } else {
      for ( i in 1:g )  {  # Saet 1-taller soejler for dem der skal merges
     	    Mer[i,grp[[i]]] <- 1 
      }
   }
   if ( !is.null(names(grp)) )
      rownames(Mer) <- names(grp)
   if ( !is.null(names) )
      colnames(Mer) <- names
   return(Mer)    # returnerer merge matrix
   # X %*% Mer    # returnerer merge input/output data
}

