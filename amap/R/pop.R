pop <- function(x,fmbvr=TRUE,triabs=TRUE,allsol=TRUE)
  {
    couts <- as.matrix(x)

    n <- as.integer(nrow(couts))

    ysave <-  as.integer(matrix(0,nrow=n,ncol=n))
    renum <- y <- ysave

    
    
    bornth <- z0 <- z <- as.double(0)
    
    res <- .Fortran("pnkfmb",
                    as.integer(fmbvr),
                    as.integer(triabs),
                    as.integer(allsol),
                    n = as.integer(n),
                    couts = as.double(couts),
                    ysave = ysave,
                    y = ysave,
                    renum= renum,
                    bornth = bornth,
                    nbcl0 = as.integer(0),
                    z0 = z0 ,
                    nbcl = as.integer(0),
                    z = z,
                    nbemp = as.integer(0),
                    nbdep = as.integer(0),
                    nbsol = as.integer(0),
                    nap = as.integer(0),
                    PACKAGE="amap")


    class(res) <- "pop"
    return(res)
    
  }


print.pop <- function(x,...)
  {
    i <- 1:x$n
    classes <- x$y[i+(i-1)*x$n]
    cat("Upper bound     (half cost)   :",x$bornth,'\n')
    cat("Final partition (half cost)   :",x$z,'\n')
    cat("Number of classes             :",x$nbcl,"\n")
    cat("Forward move count            :",x$nbemp,"\n")
    cat("Backward move count           :",x$nbdep,"\n")
    cat("Constraints evaluations count :",x$nap,"\n")
    cat("Number of local optima        :",x$nbsol,"\n\n")
    print(data.frame(Individual=i,class=classes))
  }

