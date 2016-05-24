gnesfa <- function(dudi, Focus, Reference, centering = c("single", "twice"),
                   scannf=TRUE, nfFirst=2, nfLast=0)
  {
    call <- match.call()
    centering <- match.arg(centering)
    if (!inherits(dudi, "dudi"))
      stop("dudi should be of class \"dudi\"")
    X <- dudi$tab
    if (missing(Reference))
      Reference <- dudi$lw
    if (missing(Focus))
      Focus <- dudi$lw
    if (!is.vector(Focus))
      stop("Focus should be a vector")
    if (!is.vector(Reference))
      stop("Reference should be a vector")
    if (length(Focus)!=length(Reference))
      stop("Reference and Focus of unequal length")
    if (all(Focus==Reference))
      stop("Identical Focus and Reference distributions")
    Focus <- Focus/sum(Focus)
    Reference <- Reference/sum(Reference)

    ## centering
    me <- apply(X,2,function(x) sum(x*Reference))
    Z1 <- as.matrix(sweep(X,2,me, "-"))

    if (centering == "twice") {
      m <- apply(Z1,2,function(x) sum(x*Focus))
      b <- m/sqrt(sum(m^2))
      Z <- Z1%*%(diag(ncol(Z1))-as.matrix(outer(b,b)))
    } else {
      Z <- Z1
    }

    ## First PCA
    rank <- dudi$rank
    pc1 <- dudi.pca(Z, center=FALSE, scale=FALSE,
                    row.w=Reference, nf=rank, scannf=FALSE)
    c1 <- as.matrix(pc1$c1)
    L <- sweep(Z%*%c1,2,sqrt(pc1$eig),"/")

    ## Second PCA
    pc2 <- dudi.pca(L, center=FALSE, scale=FALSE,
                    row.w=Focus, scannf=FALSE, nf=pc1$rank)
    if (scannf) {
      barplot(pc2$eig, main="Eigenvalues")
      cat("Select the first number of axes (>=1): ")
      nfFirst <- as.integer(readLines(n = 1))
      barplot(1/pc2$eig, main="1/Eigenvalues")
      cat("Select the second number of axes (>=0): ")
      nfLast <- as.integer(readLines(n = 1))
    }
    if (nfFirst <= 1)
        nfFirst <- 1
    if (nfLast <= 0)
        nfLast <- 0

    if (pc2$rank!=pc1$rank) {
      pc2 <- dudi.pca(L, center=FALSE, scale=FALSE,
                      row.w=Focus, scannf=FALSE, nf=pc2$rank)
    }

    ## Results
    res <- list()
    res$tab <- Z
    foo <- function(x) {
      if (x>0) {
        return((pc2$rank-x+1):pc2$rank)
      } else {
        return(NULL)
      }
    }
    keep <- unique(c(1:nfFirst, foo(nfLast)))
    namco <- paste("Axis", keep, sep=".")
    namli <- paste("Component", keep, sep=".")

    pc2$c1 <- pc2$c1[,keep]
    co <- as.matrix(sweep(c1,2,sqrt(pc1$eig),"/"))%*%as.matrix(pc2$c1)
    c1 <- do.call("cbind",
                  lapply(1:ncol(co), function(i) co[,i]/sqrt(sum(co[,i]^2))))
    res$l1 <- as.data.frame(Z%*%co)
    res$li <- as.data.frame(Z%*%c1)
    res$co <- as.data.frame(c1)
    row.names(res$co) <- names(dudi$tab)
    row.names(res$li) <- row.names(dudi$tab)
    row.names(res$l1) <- row.names(dudi$tab)
    names(res$co) <- namco
    names(res$li) <- namli
    names(res$l1) <- namli
    res$Reference <- Reference
    res$Focus <- Focus
    res$eig <- pc2$eig
    res$nfFirst <- nfFirst
    res$nfLast <- nfLast
    res$call <- call
    res$centering <- centering
    if (centering=="twice") {
      res$mar <- m
      res$nmar <- b
      res$projmar <- apply(Z1, 1, function(x) sum(x*b))
    }
    res$cor <- cor(dudi$tab, res$li)
    class(res) <- "gnesfa"
    return(res)
  }

print.gnesfa <- function (x, ...)
{
  if (!inherits(x, "gnesfa"))
    stop("Object of class 'gnesfa' expected")
  cat("GNESFA")
  cat("\n$call: ")
  print(x$call)
  cat("\n$centering: ")
  cat(x$centering)
  cat("\neigenvalues: ")
  l0 <- length(x$eig)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...")
  cat("\n$nfFirst:", x$nfFirst, "first axes saved")
  cat("\n$nfLast:", x$nfLast, "last axes saved")
  cat("\n")
  cat("\n")
  sumry <- array("", c(3, 4), list(1:3, c("vector", "length",
                                          "mode", "content")))
  sumry[1, ] <- c("$Reference", length(x$Reference),
                  mode(x$Reference), "Weighting matrix of reference distribution")
  sumry[2, ] <- c("$Focus", length(x$Focus), mode(x$Focus), "Weighting matrix of focus distribution")
  sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values of specialization")
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow",
                                          "ncol", "content")))
  sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
  sumry[3, ] <- c("$l1", nrow(x$l1), ncol(x$l1), "row coordinates (variance weighted by $Reference =1)")
  sumry[4, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
  sumry[5, ] <- c("$cor", nrow(x$co), ncol(x$co), "correlation between variables and axes")
  class(sumry) <- "table"
  print(sumry)
  if (x$centering=="twice") {
    cat("\nother elements: ")
    cat(names(x)[(length(x)-2):(length(x))], "\n")
  }
}



