## -*- mode: R; ess-indent-level: 2 -*-
## $Id: sca.R,v 1.17 2003/07/16 17:11:15 maechler Exp $

sca <- function(S, b = if(interactive) 5, d = 0,
                qmin = if(interactive) 0 else 5,
                corblocks = if(interactive) 0 else 0.3,
		criterion = c("csv", "blp"),
		cluster = c("median","single","complete"),
		withinblock = TRUE, invertsigns = FALSE,
		interactive = dev.interactive())
{
  ## Interactive program for Simple Component Analysis (SCA)

  ## S: correlation (or variance-covariance) matrix of data to be analyzed
  ##
  ## b: number of block-components initially proposed
  ## d: number of difference-components initially proposed
  ##
  ## qmin>0:  take over parameter d
  ##	      Select number of difference-components such that
  ##	      number of components at least equal to qmin
  ##
  ## corblocks>0: take over parameter b
  ##		 Select number of block-components such that
  ##		 "between-blocks" correlations smaller than corblocks
  ##
  ## criterion="csv" for "corrected sum of variances"
  ## criterion="blp" for "best linear predictor"
  ##
  ## withinblock=TRUE if difference-components necessarily "within-block"
  ## withinblock=FALSE if difference-components not necessarily "within-block"
  ##
  ## cluster="single" for "single linkage clustering"
  ## cluster="median" for "median linkage clustering"
  ## cluster="complete" for "complete linkage clustering"
  ##
  ## invertsigns=FALSE: the sign of the variables are not inverted
  ## invertsigns=TRUE: the sign of some variables may be inverted in order to
  ##		avoid negative correlations
  ## interactive = TRUE: interactive program active
  ## interactive =FALSE: computation without user interaction

  S <- as.matrix(S)
  dS <- dim(S)
  p <- dS[2]
  ## check the correlation matrix
  if(dS[1] != p) stop("First argument `S' must be a square matrix")
  if(!is.logical(tst <- all.equal.numeric(S, t(S), 100 * .Machine$double.eps))
     || !tst)
      stop("First argument `S' must be a symmetric matrix")
  eS <- eigen(S,symmetric = TRUE)
### MM: FIXME pos.def. test
  if(eS$value[p] < (-0.0001))
      warning("argument \"S\" is not a definite positive matrix")
  invertind <- NULL
  if(invertsigns && isdiffcomp(eS$vec[,1])) {
    p.plus <- sum(iPos <- eS$vec[,1] >= 0)
    pminus <- sum(!iPos)
    invertind <- (1:p)[if(p.plus >= pminus) !iPos else iPos]
    S[,invertind] <- -S[,invertind]
    S[invertind,] <- -S[invertind,]
    warning("some rows and columns in `S' have been inverted to avoid negative values")
  }

  S0 <- S

  ## check values of parameters

  cluster <- match.arg(cluster)
  criterion <- match.arg(criterion)

  if(corblocks > 1) corblocks <- 1
  qmin <- round(qmin)
  if(qmin > p) qmin <- p

  ## create labels

  labvar <- dimnames(S)[[1]]
  if(is.null(labvar)) labvar <- dimnames(S)[[2]]
  if(is.null(labvar)) labvar <- paste("V",1:p,sep = "")
  if(!is.null(invertind))
      labvar[invertind] <- paste("-",labvar[invertind],sep = "")

  ##---------------Calculate first b blocks and first d differences-------------

  ## We keep everything (and index by `I') such that
  ##  quick "back" and "forth" can be used interactively.
  I <- 1
  phase <- as.integer(1)
  P <- list(sortmatrix(S,diag(p)))
  dimnames(P[[I]]) <- list(labvar,paste("B",1:p,sep = ""))
  crit <- list(NA)
  nextpc <- list(NA)

  if(corblocks > 0) { ## determine b := #{block comp.} from correlations:
    maxcor <- maxmatrix(abs(corcomp(S,P[[I]])))
    while(I < p && maxcor$val > corblocks) {
      phase[I+1] <- 1
      P[[I+1]] <- agglomblock(S,P[[I]],cluster)
      crit[[I+1]] <- NA
      nextpc[[I+1]] <- NA
      I <- I+1
      maxcor <- maxmatrix(abs(corcomp(S,P[[I]])))
    }
    b <- p+1-I
  }
  else { ## check and use the argument `b'

    b <- round(b)
    if(b < 1) {
      warning("at least one block component required, setting b = 1");
      b <- 1
    }
    else if(b > p) b <- p

    while(I <= (p-b)) { ## ie. (b < p) & !(corblocks > 0)
        phase[I+1] <- 1
        P[[I+1]] <- agglomblock(S,P[[I]],cluster)
        crit[[I+1]] <- NA
        nextpc[[I+1]] <- NA
        I <- I+1
    }
  }

  if(qmin > 0)
      d <- qmin-b
  else { ## check and use the argument `d' :
    d <- round(d)
    if(b+d > p) ## should give a warning - at least if corblocks = 0 (?)
        d <- p-b
  }

  if(d > 0) {
    phase[I+1] <- 2
    P[[I+1]] <- P[[I]]
    crit[[I+1]] <- NA
    nextpc[[I+1]] <- NA
    I <- I+1
    for(i in 1:d) {
      phase[I+1] <- 2
      tmp <- nextdiff(S,P[[I]],withinblock,criterion)
      P[[I+1]] <- tmp[[1]]
      crit[[I+1]] <- NA
      nextpc[[I+1]] <- tmp[[2]]
      I <- I+1
    }
  }

  if(!interactive) {
    phase[I+1] <- 4 # directly jump to phase 4
    P[[I+1]] <- P[[I]]
    crit[[I+1]] <- allcrit(S,P[[I]],criterion,sortP = TRUE)
    I <- I+1
  }
  else { ##  if(interactive)

    clickvalid <- 0.1
    clickentries <- 0.025
    Imax <- I
    cat("Use your graphics window for interaction with sca() from now on!\n")
  }

  ## phase=1: Definition of block-components
  ## phase=2: Definition of difference-components
  ## phase=3: Deletion of difference-components
  ## phase=4: Final system (non-interactive)

  while(phase[I] < 4) { ##-------- interactive algorithm -----------------------

    ##---------------Initialize the present step--------------------------------

    ## calculate values of criteria if the system has not yet been printed:

    if(is.na(crit[[I]][[1]][1]))
	crit[[I]] <- allcrit(S,P[[I]],criterion,sortP = TRUE)

    ## initilaize the screen for this step

    Ptmp <- P[[I]]
    qI <- ncol(Ptmp)
    crittmp <- crit[[I]]
    okcomp <- rep(TRUE, qI)
    message <- ""
    warn <- ""
    if(I > 1 && phase[I] == (phase[I-1]+1) && !eqmatrix(P[[I]], P[[I-1]])) {
      warn <- "Warning: components have been reordered !"
    }
    r0 <- redrawMatrix(phase[I],I,Imax,Ptmp,crittmp,nextpc[[I]],
		       okcomp,message,warn)
    warn <- ""

    ## new   = T : the bottom left  button on the screen has been activated
    ## valid = T : the bottom right button on the screen has been activated
    ## back  = T : the  top   left  button on the screen has been activated
    ## forth = T : the  top   right button on the screen has been activated
    new <- valid <- back <- forth <- FALSE

    ##---------------Treat the present step-------------------------------------

    while(!new && !valid && !back && !forth) {

      message <- ""
      replot <- FALSE
      r <- unlist(locator(1)) ## << User Input  c(x, y)

      j <- which.min(abs(r[1]-r0[[1]]))
      i <- which.min(abs(r[2]-r0[[2]]))

      ##---------------Change the blocks----------------------------------------

      if(phase[I] == 1 && I > 1 && Ptmp[i,j] == 0 &&
	 abs(r[1]-r0[[1]][j]) <= clickentries &&
	 abs(r[2]-r0[[2]][i]) <= clickentries)
      {
	replot <- TRUE
	Ptmp[i,j] <- 1
	Ptmp[i,-j] <- 0
	okcomp <- isblockempty(Ptmp)
	if(all(okcomp)) {
	  message <- ""
	  crittmp <- allcrit(S,Ptmp,criterion,sortP = TRUE)
	} else message <- "some block is empty !"
      }

      ##---------------Change the differences-----------------------------------

      if(phase[I] == 2 && j > b &&
	 abs(r[1]-r0[[1]][j]) <= clickentries &&
	 abs(r[2]-r0[[2]][i]) <= clickentries)
      {
	replot <- TRUE
	Ptmp[i,j] <- sign(Ptmp[i,j])-1
	if(Ptmp[i,j] == (-2)) Ptmp[i,j] <- 1
	Ptmp[,j] <- simpvector(Ptmp[,j])
	okcomp[j] <- isdiffcomp(Ptmp[,j]) && (b == 1 || !withinblock ||
					      iswithinblock(Ptmp[,j],Ptmp[,1:b]))
	if(!withinblock && sum(okcomp) < qI)
	    message <- "this is not a difference !"
	if(withinblock && sum(okcomp) < qI)
	    message <- "this is not a (within-block) difference !"
	if(sum(okcomp) == qI) {
	  message <- ""
	  crittmp <- allcrit(S,Ptmp,criterion,sortP = TRUE)
	}
      }

      ##---------------Delete difference-components-----------------------------

      if(phase[I] == 3 && j > b &&
	     abs((r[1]-r0[[1]])[j]) <= clickentries &&
	 min(abs((r[2]-r0[[2]]))) <= clickentries)
      {
	new <- TRUE
	Ptmp <- Ptmp[, -j , drop = FALSE]
        qI <- ncol(Ptmp)
	## not needed: if(qI == 1) dimnames(Ptmp)[[2]] <- "B1"
	crittmp <- allcrit(S,Ptmp,criterion,sortP = TRUE)

        ## MM FIXME: shouldn't we break out in this case anyway?
        ##    -----  otherwise can happen that
      }

      ##---------------Select the next step------------------------------------

      if(sum(okcomp) == qI &&
         ((phase[I] == 1 && qI > 1) ||
          (phase[I] == 2 && qI < p)))
      {
        new <- r["x"] < 0 && all(abs(c(0,0) - r) <= clickvalid)
      }

      if(sum(okcomp) == qI)
	valid <- r["x"] > 1 && all(abs(c(1,0) - r) <= clickvalid)

      if(I > 1)
	back  <- r["x"] < 0 && all(abs(c(0,1) - r) <= clickvalid)

      if(I < Imax)
	forth <- r["x"] > 1 && all(abs(c(1,1) - r) <= clickvalid)

      if(phase[I] < 3 && new) message <- "calculation in progress ..."
      if(message != "") replot <- TRUE
      if(replot)
	  r0 <- redrawMatrix(phase[I],I,Imax,Ptmp,crittmp,nextpc[[I]],
			     okcomp,message,warn)
    } ## while(.. present step ..)

    ##---------------Prepare the next step when new=TRUE------------------------

    if(new) {
      if(phase[I] == 1) {
	if(eqmatrix(Ptmp,P[[I]])) {
	  J <- I+1
	}
	else {
	  phase[I+1] <- 1
	  P[[I+1]] <- Ptmp
	  crit[[I+1]] <- crittmp
	  nextpc[[I+1]] <- NA
	  J <- I+2
	}
	phase[J] <- 1
	P[[J]] <- agglomblock(S,Ptmp,cluster)
	crit[[J]] <- allcrit(S,P[[J]],criterion,sortP = TRUE)
	nextpc[[J]] <- NA
      }
      else if(phase[I] == 2) {
	if(eqmatrix(Ptmp,P[[I]])) {
	  J <- I+1
	} else { ## !eqmatrix(Ptmp,P[[I]])
	  phase[I+1] <- 2
	  P[[I+1]] <- Ptmp
	  crit[[I+1]] <- crittmp
	  nextpc[[I+1]] <- nextpc[[I]]
	  J <- I+2
	}
	phase[J] <- 2
	tmp <- nextdiff(S,Ptmp,withinblock,criterion)
	P[[J]] <- tmp[[1]]
	crit[[J]] <- allcrit(S,P[[J]],criterion,sortP = TRUE)
	nextpc[[J]] <- tmp[[2]]
      }
      else if(phase[I] == 3) {
	J <- I+1
	phase[J] <- 3
	P[[J]] <- Ptmp
	crit[[J]] <- crittmp
	nextpc[[J]] <- NA
      }
    } ## if(new)

    ##---------------Prepare the next step when valid=TRUE----------------------

    if(valid) {
      if(phase[I] < 3) {
	if(eqmatrix(Ptmp,P[[I]])) {
	  J <- I+1
	} else {
	  phase[I+1] <- phase[I]
	  P[[I+1]] <- Ptmp
	  crit[[I+1]] <- crittmp
	  nextpc[[I+1]] <- nextpc[[I]]
	  J <- I+2
	}
	phase[J] <- phase[I]+1
	P[[J]] <- sortmatrix(S,Ptmp)
	crit[[J]] <- allcrit(S,P[[J]], criterion, sortP = TRUE)
	nextpc[[J]] <- NA
      }
      else if(phase[I] == 3) {
	J <- I+1
	phase[J] <- 4
	P[[J]] <- Ptmp
	crit[[J]] <- crittmp
	nextpc[[J]] <- NA
      }

    }## if(valid)

    ##---------------Prepare the next step when back=TRUE or forth=TRUE---------

    if(back)
      J <- if(eqmatrix(Ptmp,P[[I]])) I - 1 else I
    if(forth) J <- I+1

    ##---------------End of the step--------------------------------------------

    if(valid && phase[I] == 1) b <- dim(P[[I]])[2]
    if(new || valid) Imax <- J
    I <- J

  } ## while( phase[] < 4 )

  ##---------------End of the process-------------------------------------------

  if(interactive) {
    okcomp <- rep(TRUE,dim(P[[I]])[2])
    message <- ""
    r0 <- redrawMatrix(phase[I],I,Imax,P[[I]],crit[[I]],nextpc[[I]],
		       okcomp,message,warn)
  }

  dimnames(S0) <- list(labvar,labvar)
  cNam <- dimnames(P[[I]])[[2]]
  dimnames(crit[[I]]$corsc) <- list(cNam,cNam)
  crit[[I]]$maxcor$row <- cNam[crit[[I]]$maxcor$row]
  crit[[I]]$maxcor$col <- cNam[crit[[I]]$maxcor$col]

  r <- list(simplemat = P[[I]],
	    loadings = normmatrix(P[[I]]),
	    allcrit  = crit[[I]],
	    nblock   = b,
	    ndiff    = dim(P[[I]])[2] - b,
	    criterion = criterion,
	    cluster  = cluster,
	    withinblock = withinblock,
	    invertsigns = invertsigns,
	    vardata  = S0)
  class(r) <- "simpcomp"
  r
} ## sca()

##----------------------------------------------------------------------------
##-------------------PRINTING SCREEN PROCEDURE--------------------------------
##----------------------------------------------------------------------------

redrawMatrix <- function(phase, I, Imax, P, crit, nextpc,
			 okcomp, message, warn)
{
  ## print the screen for a given system

  if(0 > (phase <- as.integer(phase)) || phase > 4)
      stop("`phase' must be in 1:4")
  if(is.null(dP <- dim(P))) stop("`P' must be a component matrix")
  p <- dP[1]
  q <- dP[2]
  if(q > p) stop("q > p : should never happen")
  if(length(okcomp) != q)
      warning("** okcomp[] has wrong length ", length(okcomp),
              " instead of ", q)

  ## determine  b := #{block components} <= q
  if(phase == 1)
    b <- q
  else {
    b <- 1
    while(b < q && okcomp[b+1] && !isdiffcomp(P[,b+1]))
      b <- b+1
  }

  dn <- dimnames(P)

  ## control the size of the loadings and the variables labels
  maxvar <- 20
  cex <- min(1, maxvar/p)
  maxchar <- 8
  cexlab <- cex* min(1, maxchar/ nchar(dn[[1]]))

  ## colors for block-, difference-, and error - components and commands
  col <- c(b = 2, d = 4, e = 6, c = 4)

  ## create message
  if(message == "")
    message <- paste(q,"components:    SCA =", percent(crit$cumsc,1),
                     "    PCA =", percent(crit$cumpc,1),
                     "    Opt =", percent(crit$opt,1),
                     "    Max cor =", round(crit$maxcor$val,2),
                     "(", dn[[2]][crit$maxcor$row],
                     "-", dn[[2]][crit$maxcor$col],")")

  ## coordinates strictly inside (0,1) (will be returned)
  xcord <- seq(0,1,len = q+2)[-c(1,q+2)] # length q
  ycord <- seq(1,0,len = p+2)[-c(1,p+2)] # length p

  ## plot the labels

  plot(0,1,xlim = c(-0.05,1.05), ylim = c(0,1),
       type = "n", axes = FALSE, xlab = "", ylab = "")

  par(adj = 1)
  text(rep(0,p),ycord,dn[[1]],cex = cexlab)
  par(adj = 0.5)
  text(xcord,1,dn[[2]],cex = cex)

  ## plot the system

  for(j in 1:q) {
    if(okcomp[j]) {
      par(adj = 1)
      text(rep(xcord[j],p), ycord, P[,j],
	   col = col[if(j <= b) "b" else "d"], cex = cex)
      par(adj = 0.5)
      text(rep(xcord[j],p),rep(0,p), paste(percent(crit$varsc[j],1)),
	   cex = cex)
    }
    else { ## !okcomp[j]
      par(adj = 1)
      text(rep(xcord[j],p),ycord,P[,j],col = col["e"],cex = cex)
    }
  }
  if(phase == 2 && !is.na(nextpc[1])) {
    par(adj = 1)
    text(rep(1,p),ycord, round(100 * nextpc), cex = cex)
  }

  ## plot the buttons

  par(adj = 1)
  if(phase < 4 && I > 1)
      text(0,1,"BACK",cex = 1,col = col["c"])
  if(phase == 1 && sum(okcomp) == q && q > 1)
      text(0,0,"GROUP",cex = 1,col = col["c"])
  else if(phase == 2 && sum(okcomp) == q && q < p)
      text(0,0,"NEXT",cex = 1,col = col["c"])

  par(adj = 0)
  if(phase < 4 && I < Imax)
      text(1,1,"FORW",cex = 1,col = col["c"])
  if(phase < 4 && sum(okcomp) == q)
      text(1,0,"VALID",cex = 1,col = col["c"])

  ## create titles

  if(phase <= 3)
      tit1 <- paste("Stage ",phase," (Step ",I,"):", sep="")
  tit1 <-
      switch(phase,
	     paste(tit1,"Definition of block-components"),	## 1
	     paste(tit1,"Definition of difference-components"),	## 2
	     paste(tit1,"Deletion of difference-components"),	## 3
	     "Simple components - Final system")	## 4

  cl1 <- "Click on 0s to modify blocks"
  cl2 <- paste("Click on entries of ",dn[[2]][b+1]," to change signs",sep = "")
  cl3 <- paste("Click on entries of ",dn[[2]][b+1],"-",dn[[2]][q]," to change signs",sep = "")
  cl4 <- paste("Click on entries of ",dn[[2]][b+1]," to delete component",sep = "")
  cl5 <- paste("Click on entries of ",dn[[2]][b+1],"-",dn[[2]][q]," to delete component",sep = "")
  cl6 <- "Click on GROUP to agglomerate blocks"
  cl7 <- paste("Click on NEXT to calculate D",q+1,sep = "")

  cl11 <- "There is only one block"
  cl12 <- "Maximum number of components attained"
  cl13 <- "No difference-components to delete"

  tit2 <-
      switch(phase,
	 { ## phase 1:
	   cl8 <- "Click on VALID to go to Stage 2"
	   if(q > 1 && q == p)	paste(cl6,"-",cl8)
	   else if(q > 1 && q < p)  paste(cl1,"-",cl6,"-",cl8)
	   else if(q == 1)  paste(cl11,"-",cl8)
	   else warning("q (= ncol(P)) is invalid in phase 1")
	 },

	 { ## phase 2 :
	   cl9 <- "Click on VALID to go to Stage 3"
	   if(q > b+1) {
	     if(q < p)	paste(cl3,"-",cl7,"-",cl9)
	     else ## q == p
		 paste(cl12,"-",cl3,"-",cl9)
	   }
	   else if(q == b+1) {
	     if(q < p)	paste(cl2,"-",cl7,"-",cl9)
	     else ## q == p
		 paste(cl12,"-",cl2,"-",cl9)
	   }
	   else if(q == b) {
	     if(q < p)	paste(cl7,"-",cl9)
	     else ## q == p
		 paste(cl12,"-",cl9)
	   }
	   ## else q < b : impossible
	 },

	 { ## phase 3 :
	   cl10 <- "Click on VALID to end the program"
	   if(q >  b+1)	 paste(cl5,"-",cl10)
	   else if(q == b+1)  paste(cl4,"-",cl10)
	   else if(q == b)  paste(cl13,"-",cl10)
	   ## else q < b : impossible
	 },
	 { ## phase 4 :
	   ""
	 })

  ## print titles and messages
  par(adj = 0.5)

  title(main = paste(tit1,"\n",tit2),	cex = 0.7)
  ## FIXME: use mtext() instead
  title(sub = message, mgp = c(1,1,0),	cex = 1)
  title(sub = warn,    mgp = c(2.5,1,0),cex = 1)

  list(x= xcord, y=ycord)
} ## redrawMatrix()

##----------------------------------------------------------------------------
##-------------------SCA PROCEDURES-------------------------------------------
##----------------------------------------------------------------------------

agglomblock <- function(S, P, cluster = c("median","single","complete"))
{
## agglomerate the two block-components according to cluster

  S <- as.matrix(S)
  P <- as.matrix(P)
  p <- dim(P)[1]
  b <- dim(P)[2]
  if(b == 1) return(P)
  cluster <- match.arg(cluster)

  R <- cov2cor(S)
  g <- vector("list", b)
  for(i in 1:b) g[[i]] <- (1:p)[P[,i] == 1]
  cormax <- -Inf
  for(i in 1:(b-1)) {
    for(j in (i+1):b) {
      newcor <-
	  switch(cluster,
		 "single"  = max   (R[g[[i]],g[[j]]]),
		 "median"  = median(R[g[[i]],g[[j]]]),
		 "complete"= min   (R[g[[i]],g[[j]]]))
      if(newcor >= cormax) {
	cormax <- newcor
	i0 <- i
	j0 <- j
       }
    }
  }
  P[,i0] <- P[,i0]+P[,j0]
  P <- sortmatrix(S,P[,-j0])
  dimnames(P)[[2]] <- paste("B",1:(b-1),sep = "")
  return(P)
}

##----------------------------------------------------------------------------

nextdiff <- function(S,P, withinblock, criterion = c("csv", "blp"))
{
  ## calculate the next simple difference-component
  S <- as.matrix(S)
  P <- as.matrix(P)
  criterion <- match.arg(criterion)
  zcomp <- firstpcres(S,P)
  ## Eigen vectors are only unique up to  * (+- 1)
  ## make sure the component is "unique": first (non-zero) entry positive:
  k <- 1; while(zcomp[k] == 0) k <- k + 1
  if(zcomp[k] < 0) zcomp <- -zcomp

  q <- ncol(P)

  if(!withinblock) {
    scompmax <- shrinkdiff(zcomp,S,P,criterion)$scompmax
    nextpc <- zcomp
  }
  else { ## withinblock :
    k <- 1
    critmax <- 0
    while(k <= q && !isdiffcomp(P[,k])) {
      not0 <- P[,k] != 0
      if(sum(not0) >= 2) {
	zc0 <- zcomp
	zc0[!not0] <- 0
	sd.z <- shrinkdiff(zc0, S,P, criterion)
	if(sd.z$critmax >= critmax) {
	  critmax <- sd.z$critmax
	  scompmax <- sd.z$scompmax
	  nextpc <- zc0
	}
      }
      k <- k+1
    }
  }

  P <- cbind(P,scompmax)
  dimnames(P)[[2]][q+1] <- paste("D",q+1,sep = "")

  list(P=P, nextpc=nextpc)
}

##----------------------------------------------------------------------------

shrinkdiff <- function(zcomp,S,P, criterion) {
## shrink a component towards a simple difference-component

  S <- as.matrix(S)
  P <- as.matrix(P)
  scompmax <- simpvector(sign(zcomp))
  critmax <- quickcrit(scompmax,S,P,criterion)
  abs.z <- abs(zcomp)
  cut <- sort(abs.z[abs.z > 0])
  if(length(cut) > 2) {
    for(i in 2:(length(cut)-1)) {
      newscomp <- simpvector(sign(zcomp*(abs.z >= cut[i])))
      newcrit <- quickcrit(newscomp,S,P,criterion)
      if(isdiffcomp(newscomp) && newcrit >= critmax) {
	critmax <- newcrit
	scompmax <- newscomp
      }
    }
  }
  list(scompmax=scompmax, critmax=critmax)
}

##----------------------------------------------------------------------------
##-------------------MATRICES AND VECTORS PROCEDURES--------------------------
##----------------------------------------------------------------------------

simpvector <- function(x) {
## simplify the vector x

  x <- sign(x)
  npos <- sum(x > 0)
  nneg <- sum(x < 0)
  if(nneg > 0) x[x > 0] <- nneg
  if(npos > 0) x[x < 0] <- -npos
  if(nneg > 0 || npos > 0) {
    while(sum(x%%2) == 0) x <- x/2
    while(sum(x%%3) == 0) x <- x/3
    while(sum(x%%5) == 0) x <- x/5
    while(sum(x%%7) == 0) x <- x/7
    if(nneg > 1) while(sum(x%%nneg) == 0) x <- x/nneg
    if(npos > 1) while(sum(x%%npos) == 0) x <- x/npos
  }
  return(x)
}

##----------------------------------------------------------------------------

eqmatrix <- function(A,B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  all(dim(A) == dim(B)) && all(A == B)
}

##----------------------------------------------------------------------------

isblockempty <- function(P)
  apply(as.matrix(P), 2, function(v) max(v) > 0)

##----------------------------------------------------------------------------

isdiffcomp <- function(x) { prod(range(x)) < 0 }

##----------------------------------------------------------------------------

iswithinblock <- function(d,P) {
  P <- as.matrix(P)
  b <- ncol(P)
  ## FIXME: faster using partial sort() :
  return(sort(crossprod(abs(d), abs(P)))[b-1] == 0)
}


##----------------------------------------------------------------------------
##-------------------COMPUTATIONAL PROCEDURES---------------------------------
##----------------------------------------------------------------------------

allcrit <- function(S,P, criterion, sortP = TRUE)
{
## calculate all kinds of criteria

  S <- as.matrix(S)
  P <- normmatrix(P)
  q <- ncol(P)
  tr.S <- sum(diag(S))
  varpc <- eigen(S,symmetric = TRUE)[[1]][1:q] / tr.S
  varsc <- diag(covcomp(S,P)) / tr.S
  cumpc <- sum(varpc)## == cumsum(varpc)[q]
  cumsc <- sccrit(S,P,criterion,sortP)[q]
  opt <- cumsc/cumpc
  corsc <- corcomp(S,P)
  maxcor <- maxmatrix(abs(corsc))
  list(varpc=varpc, varsc=varsc,
       cumpc=cumpc, cumsc=cumsc,
       opt=opt, corsc=corsc, maxcor=maxcor)
}

##----------------------------------------------------------------------------

sccrit <- function(S,P, criterion, sortP = TRUE) {
## return value of criterion for components P on S (cumulative)

  S <- as.matrix(S)
  if(sortP) P <- sortmatrix(S,P)
  P <- normmatrix(P)
  q <- ncol(P)
  varsc <- numeric(q)
  for(k in 1:q) {
    if(k > 1) {
      Pk <- P[,1:(k-1)]
      Ak <- crossprod(Pk, S)
      Ak <- S - crossprod(Ak, solve(Ak %*% Pk, Ak))
    } else Ak <- S
    pk <- P[,k]
    varsc[k] <-
	switch(criterion,
	       "csv" = crossprod(pk,Ak)%*% pk,
	       "blp" = (crossprod(pk, Ak %*% Ak) %*% pk) / (crossprod(pk,Ak)%*% pk))
  }

  return(cumsum(varsc) / sum(diag(S)))
}

##----------------------------------------------------------------------------

quickcrit <- function(newcomp,S,P, criterion) {
## return additional contribution of a new component to the system P on S

  S <- as.matrix(S)
  Pk <- normmatrix(P)
  Ak <- crossprod(Pk, S)
  Ak <- S - crossprod(Ak, solve(Ak %*% Pk, Ak))

  pk <- normmatrix(newcomp)
  switch(criterion,
	 "csv" = crossprod(pk, Ak)%*% pk,
	 "blp" = (crossprod(pk, Ak %*% Ak) %*% pk) / (crossprod(pk, Ak)%*% pk)
         )
}

##----------------------------------------------------------------------------

firstpcres <- function(S,P) {
## return the first principal component of residuals of S given components P

  S <- as.matrix(S)
  P <- normmatrix(P)
  A <- crossprod(P, S)
  A <- S - crossprod(A, solve(A %*% P, A))

  pc <- eigen(A, symmetric = TRUE)[[2]][,1]
  return(pc)
}

corcomp <- function(S,P) {
## return the correlation matrix of the components P on S

  P <- normmatrix(P)
  A <- crossprod(P, as.matrix(S)) %*% P ## the covariance
  return(cov2cor(A))
}

##----------------------------------------------------------------------------

## FIXME: efficiency : this is only used as  diag(covcomp(S,P)) !!!
covcomp <- function(S,P) {
## return the variance-covariance matrix of the components P on S

  P <- normmatrix(P)
  A <- crossprod(P, as.matrix(S)) %*% P
  return(A)
}

##----------------------------------------------------------------------------

normmatrix <- function(P) {
## normalize the matrix P

  P <- as.matrix(P)
  sweep(P,2,apply(P,2,function(v) { sqrt(sum(v^2)) }),"/")
}

##----------------------------------------------------------------------------
##-------------------SORTING PROCEDURES---------------------------------------
##----------------------------------------------------------------------------

maxmatrix <- function(R) {
## return position and value of the largest element of a correlation matrix R
## (without taking into account the diagonal elements)

  R <- as.matrix(R)
  p <- dim(R)[1]
  if(p <= 1) {
    row <- NA
    col <- NA
    val <- NA
  } else { ## (p > 1)
    i <- sort.list(R)[p*(p-1)]
    row <- i%%p+p*(i%%p == 0)
    col <- 1+i%/%(p+1/p)
    val <- R[row,col]
    if(row > col) { ## swap
      tmp <- col; col <- row; row <- tmp
    }
  }
  list(row=row, col=col, val=val)
}

##----------------------------------------------------------------------------

sortmatrix <- function(S,P) {
## sort a matrix P by decreasing variances of components
## (where the block-components are first)

  S <- as.matrix(S)
  P <- as.matrix(P)
  q <- ncol(P)
  diffcomp <- apply(P,2,isdiffcomp)
  iblock <- (1:q)[!diffcomp]
  i.diff <- (1:q)[diffcomp]
  d <- diag(covcomp(S,P))
  P[,] <- P[,c(iblock[sort.list(-d[iblock])],
	       i.diff[sort.list(-d[i.diff])])]
  return(P)
}

##----------------------------------------------------------------------------
##-------------------PRINTING CRITERIA PROCEDURE------------------------------
##----------------------------------------------------------------------------

percent <-  {
  if(is.R())## work around R's format(*, nsmall = . ) bug:
      function(p, d = 0, sep=" ")
          paste(formatC(round(100 * p, d), format="f", digits= d), "%",
                sep = sep)
  else
      function(p, d = 0, sep=" ")
          paste(format(round(100 * p, d), nsmall= d), "%", sep = sep)
}

print.simpcomp <- function(x, ndec = 2, ...) {

  ## print method for the result list of the sca() function above

  catn <- function(...) cat(..., "\n")

  labcrit <-
      switch(x$criterion,
             "csv" = "corrected sum of variances",
             "blp" = "best linear predictor")

  catn("------------------------------------------------------------")
  catn("Simple Component Analysis")
  catn("------------------------------------------------------------")
  catn("Optimality criterion        :",labcrit)
  catn("Clustering procedure        :",paste(x$cluster,"linkage"))
  catn("Within-block differences    :",x$withinblock)
  catn("Possible invertion of signs :",x$invertsigns)
  catn("Number of block-components  :",x$nblock)
  catn("Number of diff.-components  :",x$ndiff)
  catn("------------------------------------------------------------")
  catn("Simple matrix:")
  print(x$simplemat)

  AC <- x$allcrit
  catn("------------------------------------------------------------")
  catn("Variance principal components:\n", percent(AC$varpc,ndec))
  catn("Variance simple components   :\n", percent(AC$varsc,ndec))
  catn("------------------------------------------------------------")
  catn("Extracted variability PCA:", percent(AC$cumpc,ndec))
  catn("Extracted variability SCA:", percent(AC$cumsc,ndec))
  catn("Optimality SCA           :", percent(AC$opt,ndec))
  catn("------------------------------------------------------------")
  catn("Correlations simple components:")
  print(round(AC$corsc,ndec))
  catn("------------------------------------------------------------")
  cat("Max (abs) correlation: ")
  catn(round(AC$maxcor$val,ndec), "(",AC$maxcor$row,"-",AC$maxcor$col,")")
  catn("------------------------------------------------------------")
  invisible(x)
}
