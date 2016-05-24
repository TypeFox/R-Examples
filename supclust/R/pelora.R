## Pelora - a supervised algorithm for grouping predictor variables
pelora <- function(x, y, u = NULL, noc = 10, lambda = 1/32, flip = "pm",
                   standardize = TRUE, trace = 1)
  {
    ## Bedeutung der Input-Variablen
    ## -----------------------------
    ## x            Expressionmatrix im Format (n x p)
    ## y            Binärer Responsevektor, codiert mit 0 und 1
    ## u            Klinische Variablen im Format (n x m)
    ## noc          Anzahl Variablen, die ins Modell aufgenommen werden
    ## lambda       Reskalierter Penalty-Parameter, sollte in [0,1] sein
    ## flip         Sign-Flipping in der x-Matrix: ständig, einmalig, nicht
    ## standardize  Variablen-Standardisierung in der x-Matrix: ja oder nein?
    ## trace        Run-Control Output: 0 nichts, 1 moderat, 2 ausführlich


    ## Check input
    if(!is.matrix(x) || !is.numeric(x))
        stop("'x' must be a numeric matrix (e.g. gene expressions)")
    n  <- nrow(x)
    if(!is.numeric(y) || length(y) != n || any(y != 0 & y != 1))
        stop("'y' must be a numeric vector of length n = ", nrow(x),
             " with only 0/1 entries")
    yvals <- 0:1 # the y-values, aka "class labels"

    ## Sign-Flip
    if (flip=="pm")
      {
        x      <- cbind(x, -x)
        signs  <- c(rep(0,ncol(x)))
      }
    if (flip=="cor")
      {
        sgnChg <- sign.change(x,y)
        x      <- sgnChg$x
        signs  <- sgnChg$signs
      }
    if (flip=="none")
      {
        signs  <- NULL
      }

    ## Standardisierung der x-Variablen
    if (standardize)
      {
        stndz  <- standardize.genes(x)
        x      <- stndz$x
        means  <- stndz$means
        sdevs  <- stndz$sdevs
        if (any(sdevs==0))
          stop ("There are predictor variables with st. dev. = 0")
      }
    else
      {
        means <- NULL
        sdevs <- NULL
      }

    ## Setting up the clinical variables
    if((have.u <- !is.null(u))) { ## with or without "clinical variables"?
        if(!is.matrix(u) || !is.numeric(u))
            stop("'u' must be a numeric matrix (e.g. clinical variables)")
        m <- ncol(u)
        E <- cbind(x, u)
    } else {
        m <- 0
        E <- x
    }

    ## Initialisierung
    g         <- ncol(E) # if flip %in% c("cor","none") = px + m
    p         <- noc+1
    px        <- ncol(x)
    X         <- matrix(0, n, p)
    X[,1]     <- 1
    P         <- matrix(0, p, p) ##= diag(apply(X,2,var)) * lambda*n
    theta0    <- log(mean(y)/(1-mean(y)))
    theta     <- c(theta0, rep(0, noc))
    prob      <- rep(1/(1 + exp(-theta0)), n)## == 1/(1 + exp(-X %*% theta))
    W         <- diag(prob*(1-prob))

    ## Earlier options that are now held fixed:
    penflag   <- 0 ## for the standard penalty, !=0 for the nonstandard penalty
    blockflag <- 1 ## for not allowing a gene to enter a cluster more than once
    critflag  <- 0 ## for the log-likelihood criterion, 1 for the L2 criterion
    valiflag  <- 1 ## for validation of genes via backdeletion, 0 to do without

    lSize     <- 2*g*p

    ## Aufruf der C-Funktion
    res <- .C(R_clusterer,
              E = 	as.double(E),
              X = 	as.double(X),
              W = 	as.double(W),
              P = 	as.double(P),
              y = 	as.double(y),
              prob = 	as.double(prob),
              theta =   as.double(theta),
              lambda=   as.double(lambda*n),
              n = 	as.integer(n),
              g = 	as.integer(g),
              m = 	as.integer(m),
              p = 	as.integer(p),
              critflag = as.integer(critflag),
              penflag  = as.integer(penflag),
              blockflag= as.integer(blockflag),
              valiflag = as.integer(valiflag),
              traceflag= as.integer(trace),
              ## Output:
              genliste = 	integer(lSize),
              kriterium = 	double (lSize),
              DUP = FALSE)[c("genliste", "kriterium")]

    ## Auswertung, bilden der Liste mit Genen und Mittelwerten
    i         <- 1
    genliste  <- res$genliste
    kriterium <- res$kriterium
    genes     <- integer()
    cCrit     <- double()
    alle      <- list()
    crit      <- list()
    while (genliste[i] != 0) {

        if (genliste[i] == -1) {
            i <- i+1
            genes <- c(genes, genliste [i])
            cCrit <- c(cCrit, kriterium[i])
            i <- i+1
        }
        if (genliste[i] == -2) {
            kickout <- genliste[i+1]
            iKick <- which(genes == kickout)
            genes <- genes[-iKick]
            cCrit <- cCrit[-iKick]
            i <- i+2
        }
        if (genliste[i] == -4) { # cluster end
            alle    <- c(alle, list(genes))
            crit    <- c(crit, list(cCrit))
            genes <- integer()
            cCrit <- double()
            i <- i+1
        }
    }
    noclu <- length(alle)

    var.type <- factor(sapply(alle, function(j) any(j > px)),
                       levels = c(FALSE, TRUE),
                       labels = c("Cluster", "Clinical"))
    values <- sapply(alle, function(ids) rowMeans(E[, ids, drop=FALSE]))

    ## FIXME: The samp.names should be attached to 'values' / kept
    ## -----  where available;  values should get "predictor names" !!
    dnE <- dimnames(E)
    if(is.null(sNames <- dnE[[1]]))
        sNames <- names(y)
    if(is.null(sNames)) sNames <- as.character(1:n)

    ## Output
    out <- list(genes = alle, values = values, y = y, yvals = yvals,
                lambda = lambda, noc = noc, px = px, flip = flip,
                var.type = var.type, crit = crit, signs = signs,
                samp.names = sNames, gene.names = dnE[[2]], call = match.call())
    class(out) <- "pelora"
    out
  }


## Returns the fitted values
fitted.pelora <- function(object, ...)
  {
    ## Fitted values
    out <- object$values
    if (is.null(vNames <- colnames(out)))
      vNames      <- paste("Predictor", 1:object$noc)
    dimnames(out) <- list(object$samp.names, vNames)
    out
  }


## A short overview of what has been found by Pelora
print.pelora <- function(x, digits = getOption("digits"), details = FALSE, ...)
  {
    ## Preliminaries
    noc      <- x$noc
    fin.crit <- sapply(1:noc, function(i) {cr <- x$crit[[i]]; cr[length(cr)]})
    nGenes   <- unlist(lapply(x$genes,length))
    smDig    <- max(2, digits - 4)
    cCrit    <- format(round(fin.crit, smDig), nsmall = smDig, digits = digits)
    cG       <- format(nGenes)
    cI       <- format(1:noc)

    ## General overview
    cat("\nPelora called with lambda = ", x$lambda, ",", sep="")
    isClust <- x$var.type == "Cluster"
    if (allClust <- all(isClust)) ##  predictors are clusters only
      cat(" ", noc, " cluster", if(noc > 1)"s"," fitted\n\n", sep="")
    else ## predictors are both clusters *and* clinical variables
      cat("\n", (nC <- sum(isClust)), " cluster", if(nC>1)"s", " and ",
          noc - nC, " clinical variable", if ((noc-nC)>1)"s",
          " fitted\n\n", sep="")

    ## Information about each of the predictors
    for (i in 1:noc)
      {
	gic <- x$genes[[i]]
        hic <- x$genes[[i]]
        if(!is.null(x$signs) && any(x$signs == 0))
          hic[hic>(x$px/2)] <- hic[hic>(x$px/2)]-(x$px/2)
	if (isClust[i])
          {
	    ## the i-th predictor is a gene cluster
	    ng <- length(gic) # == nGenes[i]
	    if(allClust)
              cat("Cluster", cI[i], ": Contains ")
	    else
              cat("Predictor", cI[i], ": Cluster with ")
	    cat(cG[i], " gene", if(ng != 1)"s",
		", final criterion ", cCrit[i],"\n", sep="")

	    if(details)
              {
		## printing all the genes
		cj  <- format(1:ng)
		cgi <- format(hic)
		for (j in 1:ng)
                  {
                    cat("Entry", cj[j], ": Gene", cgi[j])
                    if (!is.null(x$signs))
                      {
                        if (x$signs[gic[j]] == 0)
                          {
                            cat(if (gic[j]>(x$px/2))
                                " (flipped)" else "          ")
                          }
                        else
                          {
                            cat(if (x$signs[gic[j]] == -1)
                                " (flipped)" else "          ")
                          }
                        if (!is.null(x$gene.names))
                          cat(" : Name", x$gene.names[gic[j]])
                        cat("\n")
                      }
                  }
              }
          }
	else
          { ## x$var.type[i] == "Clinical"
	    cat("Predictor ", cI[i], " : Clinical variable ",
                gic[1]-x$px,
		if (!is.null(x$gene.names))
		paste(" named '", x$gene.names[gic[1]], "'"),
		", final criterion ", cCrit[i], "\n", sep="")
          }
	if(details) cat("\n")
      } ## for(j )
    cat("\n")
    invisible(x)
  }

## Yields printed output about the clustering in more detail
summary.pelora <- function(object, digits = getOption("digits"), ...)
  {
    print(object, details = TRUE, digits = digits, ...)
  }


## Plots a 2-dimensional projection of the first 2 Pelora-clusters
plot.pelora <- function(x, main = "2-Dimensional Projection Pelora's output",
                        xlab = NULL, ylab = NULL, col = seq(x$yvals), ...)
  {
    if(x$noc <= 1)
      {
        if(is.null(xlab))
          xlab <- if (x$var.type[1] == "Cluster")
            "Mean Expression of Pelora's Cluster 1"
          else "Expression of a Clinical Variable"

        if(is.null(ylab))
          ylab <- "Class label"

        plot(x$values[,1], x$y, type="n", xlab=xlab, ylab=ylab, main=main, ...)
        text(x$values[,1], x$y, x$y, col = col[1 + x$y])
        return()
      }

    if(is.null(xlab))
    xlab <-
      if (x$var.type[1] == "Cluster")
        "Mean Expression of Pelora's Cluster 1"
      else "Expression of a Clinical Variable"

    if(is.null(ylab))
    ylab <-
        if (x$var.type[2] == "Cluster")
            "Mean Expression of Pelora's Cluster 2"
        else "Expression of a Clinical Variable"

    plot(x$values[,1], x$values[,2], type = "n",
         xlab = xlab, ylab = ylab, main = main, ...)
    text(x$values[,1], x$values[,2], x$y, col = col[1 + x$y])
    invisible()
}


## Returns the coefficients of the penalized logistic regression classifier
coef.pelora <- function(object, ...)
  {
    clusnames  <- character(object$noc)
    for (i in 1:object$noc) clusnames[i] <- paste("Predictor", i)
    out        <- ridge.coef(fitted(object), object$y, object$lambda)
    names(out) <- c("Intercept", clusnames)
    out
  }


## Predictions with Pelora's clusters
predict.pelora <- function(object, newdata = NULL, newclin = NULL,
                           type = c("fitted", "probs", "class"),
                           noc = object$noc, ...)
  {
    ## Checking the input
    type <- match.arg(type)
    stopifnot(length(noc) >= 1, length(object$noc) == 1, object$noc >= 1)
    if (max(noc) > object$noc)
      stop("You cannot predict with more predictors than you have fitted")

    ## Return fitted values, probabilities or class labels for the training data
    if (is.null(newdata))
      {
        X   <- cbind(rep(1,length(object$y)),fitted(object))

        if (length(noc)==1) {
          if(noc==object$noc)
          {
            prvec           <- matrix(1/(1+exp(-c(X%*%coef(object)))), ncol = 1)
          } else { ## length(noc)==1 && noc < object$noc
            koeff <- ridge.coef(object$val[, 1:noc], object$y, object$lambda)
            prvec <- matrix(c((1/(1+exp(-(X[, 1:(noc+1)]%*%koeff))))), ncol = 1)
          }
            dimnames(prvec) <- list(object$samp.names, paste(noc, "Predictors"))
            return(switch(type,
                          fitted = fitted(object),
                          probs  = prvec,
                          class  = (prvec > 0.5)*1))
        }
        if (length(noc)>1 & max(noc)<=object$noc)
          {
            prmt <- NULL
            for (i in 1:length(noc)) {
              koeff <- ridge.coef(object$val[,1:(noc[i])],object$y,object$lamb)
              prvec <- matrix(c((1/(1+exp(-(X[,1:(noc[i]+1)]%*%koeff))))),nc=1)
              dimnames(prvec) <- list(object$samp.n,paste(noc[i],"Predictors"))
              prmt  <- cbind(prmt, prvec)
            }
            return(switch(type,
                          fitted = fitted(object)[, noc, drop = FALSE],
                          probs  = prmt,
                          class  = (prmt>0.5)*1))
          }
      }
    else ## Returning fitted values, probabilities or class labels for 'newdata'
      {
        ## Customizing newdata according to the choice of the flipping method
        if (object$flip=="pm") newdata <- cbind(newdata, -newdata)

        ## Check if new clinical variables are provided too
        if (is.null(newclin) && any(object$var.type=="Clinical"))
          stop("You also need to provide new clinical variables")

        ## Check the dimensions of the new data
        if (ncol(newdata)!=object$px)
          stop("The new data need to have the same number of ",
               "predictor variables as the training data")

        ## Flip the signs of the new data
        if (object$flip=="cor") newdata <- t(t(newdata)*object$signs)

        ## Standardization of the new data
        if (!is.null(object$means))
          {
            for (i in 1:(ncol(newdata)))
              {
                newdata[,i] <- (newdata[,i]-object$means[i])/object$sdevs[i]
              }
          }

        ## Merge expression data and clinical variables (if available)
        if (!is.null(newclin)) newdata <- cbind(newdata, newclin)

        ## Determine the fitted values
        Xt <- cbind(rep(1,nrow(newdata)))
        for (j in 1:object$noc) {
          Xt <- cbind(Xt, rowMeans(newdata[,object$genes[[j]], drop = FALSE]))
        }

        ## Naming the predictors
        sampnames <- 1:nrow(newdata)
        clusnames <- paste("Predictor", 1:object$noc)
        dimnames(Xt) <- list(sampnames, c("Intercept", clusnames))

        if (length(noc)==1) {
          eta <-
            if(noc == object$noc) {
              Xt %*% coef(object)
            } else { ## length(noc) == 1  &&  noc < object$noc
              koeff <- ridge.coef(object$val[,1:noc],object$y,object$lambda)
              Xt[,1:(noc+1)] %*% koeff
            }
          prvec <- matrix(1 /(1 + exp(- drop(eta))), ncol = 1,
                          dimnames = list(1:nrow(newdata), paste(noc, "Predictors")))
          switch(type,
                 fitted = Xt[,2:ncol(Xt)],
                 probs  = prvec,
                 class  = (prvec > 0.5)*1)
        }
        else { ##  (length(noc) > 1 && max(noc) <= object$noc)
          prmt <- NULL
          for (i in 1:length(noc)) {
            koeff <- ridge.coef(object$val[,1:(noc[i])],object$y,object$lamb)
            prvec <- matrix(1/ (1 + exp(-Xt[,1:(noc[i]+1)] %*% koeff)),
                            ncol = 1,
                            dimnames = list(1:nrow(newdata),paste(noc[i],"Predictors")))
            prmt  <- cbind(prmt, prvec)
          }
          switch(type,
                 fitted = Xt[,2:ncol(Xt)],
                 probs  = prmt,
                 class  = (prmt>0.5)*1)
        }
      }
  }












