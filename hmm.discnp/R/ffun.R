ffun <- local({

fixNA <- function(x,Rho) {
  if(any(is.na(x))) {
    rn <- c(rownames(Rho),"missing")
    Rho  <- rbind(Rho,rep(1,ncol(Rho)))
    rownames(Rho) <- rn
    x[is.na(x)] <- "missing"
  }
  Rho[x,1:ncol(Rho)]
}

function(y,Rho)
{
#
# Function ffun to calculate f(x) = Pr(Y=x | the model parameters)
# for each entry of each vector in the list y, for each value of the
# state k, k = 1, ..., K.  The returned result, fy, is a matrix whose
# ***rows*** correspond to the states and whose ***columns***
# correspond to the observations y.
#

fy <- lapply(y,fixNA,Rho=Rho)
fy <- do.call(rbind,fy)
t(fy)
}
})
