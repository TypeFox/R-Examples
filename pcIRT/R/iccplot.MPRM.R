#'@rdname iccplot
#'@method iccplot MPRM
#'@export

iccplot.MPRM <- function(object, items="all", ...){

  if(nrow(object$itempar) > 3){stop("only up to three dimensions can be plotted!")}

  pp1 <- seq(-5,5, by=0.5)
  pp2 <- seq(-5,5, by=0.5)

if(is.numeric(items)){
  y <- sapply(items, function(l){

        pp_comb <- expand.grid(theta1 = pp1, theta2 = pp2)
        m1 <- exp(pp_comb[,1] - object$itempar[1,l])
        m2 <- exp(pp_comb[,2] - object$itempar[2,l])

        n <- rowSums(cbind(m1,m2, exp(0)))
        respc1 <- m1/n
        respc2 <- m2/n
        r2 <- cbind(respc1, respc2)
        respc3 <- 1-rowSums(r2)

        mat <- cbind(respc1, respc2, respc3)
        colnames(mat) <- paste0("cat", 1:3)
        par(ask=TRUE)
        z <- sapply(1:ncol(mat), function(p1){
          zmat <- matrix(mat[,p1], ncol=length(pp1))
          persp(pp1, pp2, zmat, xlab="theta 1", ylab="theta 2", zlab="response probability", theta=35, phi=25, ticktype="detailed", main=paste0(colnames(object$itempar)[l]," ", colnames(mat)[p1]))
        })
        par(ask=FALSE)
  })
} else if(items == "all"){
  y <- sapply(1:ncol(object$itempar), function(l){

    pp_comb <- expand.grid(theta1 = pp1, theta2 = pp2)
    m1 <- exp(pp_comb[,1] - object$itempar[1,l])
    m2 <- exp(pp_comb[,2] - object$itempar[2,l])

    n <- rowSums(cbind(m1,m2, exp(0)))
    respc1 <- m1/n
    respc2 <- m2/n
    r2 <- cbind(respc1, respc2)
    respc3 <- 1-rowSums(r2)

    mat <- cbind(respc1, respc2, respc3)
    colnames(mat) <- paste0("cat", 1:3)
    par(ask=TRUE)
    z <- sapply(1:ncol(mat), function(p1){
      zmat <- matrix(mat[,p1], ncol=length(pp1))
      persp(pp1, pp2, zmat, xlab="theta 1", ylab="theta 2", zlab="response probability", theta=35, phi=25, ticktype="detailed", main=paste0(colnames(object$itempar)[l]," ", colnames(mat)[p1]))
    })
    par(ask=FALSE)
  })
} else {stop("Items must be a numeric vector to choose a subset of items or must be 'all' to choose all items")
}
}

