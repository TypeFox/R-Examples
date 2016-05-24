#
eigenmap <- function(x,opt.coord=NA,weighting=Wf.sqrd,boundaries,wpar,select=.Machine$double.eps^0.5) {
#
  if(!is.function(weighting)) stop("Parameter 'weighting' must be a weighting function.")
  if(!is.numeric(x)) stop("Parameter 'x' must be numeric!")
  if(inherits(x,"dist")) {
    D <- as.matrix(x)
    if(all(is.na(opt.coord))) {
      opt.coord <- as.matrix(x=1:nrow(D))
      rownames(opt.coord) <- rownames(D)
    } else {
      opt.coord <- as.matrix(opt.coord)
      if(nrow(opt.coord) != nrow(D)) {
        stop(paste("You provided",nrow(opt.coord),"optional coordinates to reference", nrow(D),"observations!"))
      } else rownames(opt.coord) <- rownames(D)
    }
  } else {
    opt.coord <- as.matrix(x)
    x <- dist(x,method="euclidean")
    D <- as.matrix(x)
    rownames(opt.coord) <- rownames(D)
  }
  wformals <- names(formals(weighting))
  if(any(wformals == "boundaries") && missing(boundaries)) {
    boundaries <- c(0,max(hclust(x,method="single")$height))
    warning("No boundaries given, they were set to (",boundaries[1L],"; ",boundaries[2L],")")
  }
  if(any(wformals == "boundaries")) {
    if(any(wformals == "wpar") && !missing(wpar))
      W <- weighting(D,boundaries,wpar)
    else
      W <- weighting(D,boundaries)
  } else {
    if(any(wformals == "wpar") && !missing(wpar))
      W <- weighting(D,wpar)
    else
      W <- weighting(D)
  }
  diag(W) <- 0
  n <- nrow(W)
  colWbar <- matrix(colMeans(W),1L,n)
  MW <- mean(W)
  term <- matrix(1,n,1L) %*% colWbar
  O <- W - term - t(term) + MW
  eigO <- eigen(O)
  variables <- eigO$vectors[,abs(eigO$values) >= select]
  rownames(variables) <- rownames(D) ; colnames(variables) <- paste("dbMEM",1L:ncol(variables),sep="")
  return(structure(list(
    coordinates=opt.coord,
    D=D,weighting=weighting,
    boundaries=if(!missing(boundaries)) boundaries else NULL,
    wpar=if(!missing(wpar)) wpar else NULL,
    W=W,colWbar=colWbar,MW=MW,
    lambda=eigO$values[abs(eigO$values) >= select],
    U=variables),class="eigenmap"))
}
#
Wf.sqrd <- function(D) return(-0.5 * D)
#
Wf.binary <- function(D,boundaries) {
  b <- matrix(0, nrow(D), ncol(D))
  b[D > boundaries[1L] & D <= boundaries[2L]] <- 1
  return(b)
}
#
Wf.PCNM <- function(D,boundaries) {
  b <- Wf.binary(D, boundaries)
  a <- 1 - (D / (4 * boundaries[2L]))**2
  W <- b * a
  return(W)
}
#
Wf.Drayf1 <- function(D,boundaries) {
  b <- Wf.binary(D, boundaries)
  a <- 1 - D / max(D)
  W <- b * a
  return(W)
}
#
Wf.Drayf2 <- function(D,boundaries,wpar=1) {
  b <- Wf.binary(D, boundaries)
  a <- 1 - (D / max(D))^wpar
  W <- b * a
  return(W)
}
#
Wf.Drayf3 <- function(D,boundaries,wpar=1) {
  b <- Wf.binary(D, boundaries)
  a <- 1 / D^wpar
  W <- b * a
  return(W)
}
#
eigenmap.score <- function(object,target) {
  if(!is.matrix(target)) attr(target,"dim") <- c(1L,length(target))
  n <- ncol(object$D)
  if(n != ncol(target))
    stop("Mismatch distances to target: ",n," but ",ncol(object$D)," expected.")
  ntgt <- nrow(target)
  #
  if(!is.null(object$boundaries)) {
    if(!is.null(object$wpar))
      Wnk <- object$weighting(target,object$boundaries,object$wpar)
    else
      Wnk <- object$weighting(target,object$boundaries)
  } else {
    if(!is.null(object$wpar))
      Wnk <- object$weighting(target,object$wpar)
    else
      Wnk <- object$weighting(target)
  }
  Wnk[target==0] <- 0
  scores <- (Wnk - matrix(rowMeans(Wnk),ntgt,1L) %*% matrix(1,1L,n) -
             matrix(1,ntgt,1L) %*% object$colWbar + object$MW) %*%
             object$U %*% diag(object$lambda**(-1))
  colnames(scores) <- colnames(object$U)
  return(scores)
}
#
print.eigenmap <- function(x,...) {
  cat(paste("\nMoran's eigenvector map containing",length(x$lambda),"basis functions.\n"))
  cat(paste("Functions span",nrow(x$U),"observations.\n\n"))
  cat(paste("Eigenvalues:\n"))
  print.default(x$lambda)
  cat("\n")
  return(invisible(NULL))
}
#
plot.eigenmap <- function(x,...) {
  if (ncol(x$coordinates) > 2) {
    warning(paste(ncol(x$coordinates),"dimensions were provided but only the first 2 were used for plotting."))
  }
  cat("Left-click on the graphical display to see further variables or right-click (Mac: esc) to terminate plotting.\n")
  for (i in 1:length(x$lambda)) {
    layout(matrix(c(1,2,2),1,3))
    plot(y=x$lambda,x=1:length(x$lambda),ylab=expression(lambda),xlab="Order",type="b")
    title("Eigenvalues diagram")
    points(y=x$lambda[i],x=i,pch=21,bg="black")
    if (ncol(x$coordinates) == 1) plot(y=x$U[,i],x=x$coordinates[,1],ylab="value",xlab="Location",type="l")
    if (ncol(x$coordinates) > 1) {
      plot(y=x$coordinates[,2],x=x$coordinates[,1],
           asp=1,ylab="Location (y)",xlab="Location (x)",type="p",pch=3)
      gcol <- grey((sign(x$U[,i])+1)/2)
      gsize <- 3*abs(x$U[,i])/max(abs(x$U[,i]))
      points(y=x$coordinates[,2],x=x$coordinates[,1], pch = 21, bg=gcol, cex = gsize)
    }
    title(paste("Variable display:",colnames(x$U)[i]))
    ttt <- locator(1)
    if (is.null(ttt)) break
  }
  return(invisible(NULL))
}
#
