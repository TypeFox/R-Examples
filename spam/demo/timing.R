# This is file ../spam/demo/timing.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     








# We construct a few large matrices and we compare how much faster (slower)
# we are compared to the full matrix analysis.
# Since all the calculation are also done with full matrices, we do not
# exagerate with the sizes.

set.seed(14)

compare <- function(expr1,expr2,tag=NULL)
  {
    if( !is.null(tag)) cat( "Comparing: ", tag, fill=TRUE)
    print(data.frame(full=system.time( expr1, TRUE)[1:3],
                     sparse=system.time( expr2, TRUE)[1:3],
                     row.names=c("user","system","elapsed")))
  }

xn <- 1000
xm <- 1200

# first start with a full matrix.
fmat1 <- matrix(rnorm(xn*xm),xn,xm)
smat1 <- as.spam(fmat1)

compare(fmat2 <- t(fmat1), smat2 <- t(smat1), "Transpose")

compare(ffmat <- fmat1 %*% fmat2,
        ssmat <- smat1 %*% smat2, "multiplication")

compare( solve(ffmat),  solve(ssmat), "solving")


compare(rbind(fmat1,fmat1),rbind(smat1,smat1))
compare(cbind(fmat1,fmat1),cbind(smat1,smat1))





# now create a sparse matrix.
fmat1[fmat1<3] <- 0
smat1 <- as.spam(fmat1)



compare(fmat2 <- t(fmat1), smat2 <- t(smat1), "Transpose")

compare(ffmat <- fmat1 %*% fmat2,
        ssmat <- smat1 %*% smat2, "multiplication")

compare(ffmat <- ffmat + diag(xn),
        ssmat <- ssmat + diag.spam(xn), "add identity")

compare(ffmat <- 1:xn %d+% ffmat,
        ssmat <- 1:xn %d+% ssmat, "add identity quicker")

compare( solve(ffmat),  solve(ssmat), "solving")

summary(ssmat)


# compare a few cbind/rbinds

compare(rbind(fmat1,fmat1),rbind(smat1,smat1))
compare(cbind(fmat1,fmat1),cbind(smat1,smat1))


