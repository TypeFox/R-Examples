##  SEQIMPUTE is a sequential imputation method.
##
##  The seqimpute method is described in:
##     S. Verboven, K. Vanden Branden and P. Goos (2007)
##     "Sequential Imputation for missing values", Computational Biology and
##     Chemistry, 31: 320-327.
##
##
##  a) If the data set is complete, just return it without doing anything.
##  b) The algorithm needs at least p complete cases - if there are no or less than p
##      complete cases, find a way to impute them in some other way.
##
##  library(rrcov)
##  data(phosphor)
##  x <- phosphor[,1:2]
##  x[10,2] <- NA
##  x[15,1] <- NA
##  y <- impSeq(x)
##  y[10,2]     # 41.02312
##  y[15,1]     # 16.60647
##  cbind(phosphor[,1:2],x,y)
##

impSeq <- function(x){

    if(is.data.frame(x))
    {
        x <- data.matrix(x)
    }else if(!is.matrix(x))
    {
        x <- matrix(x, length(x), 1, dimnames = list(names(x), deparse(substitute(x))))
    }
    xcall <- match.call()

    n <- nrow(x)
    p <- ncol(x)
    isnanx = is.na(x) + 0
    risnanx = apply(isnanx,1,sum) # observations with missing values

    if(length(which(risnanx > 0)) == 0)         # no missing data - return x
        return(x)


    ## sort according to percentage of missing values (so that first the observations
    ## with smallest number of missing values are handled)

    ##[sortx,Ix] = sort(risnanx);
    sortx <- sort.int(risnanx, index.return=TRUE)
    sorth <- sort.int(sortx$ix, index.return=TRUE)
    x = x[sortx$ix,]

##isnanx = isnanx(Ix,:);
##risnanx = sortx; %observations with missings
    isnanx = is.na(x) + 0
    risnanx = sortx$x #observations with missing values

    complobs = which(risnanx == 0)
    misobs = which(risnanx != 0)
    nmisobs = length(misobs)
    ncomplobs = length(complobs)

    ## start sequential imputation
    for(inn in 1:nmisobs)
    {
        if(inn == 1)
        {
            covx = cov(x[complobs,])
            mx = colMeans(x[complobs,])
        }else #updating formula
        {
            mxo = mx
            mx = ((ncomplobs-1)*mx + x[misobs[inn-1],])/ncomplobs
            covx = (ncomplobs-2)/(ncomplobs-1)*covx + 1/(ncomplobs-1) *
                    as.matrix(x[misobs[inn-1],] - mx) %*% t(as.matrix(x[misobs[inn-1],] - mx)) +
                    as.matrix(mxo - mx) %*% t(as.matrix(mxo - mx))
        }

        if(p >= length(complobs))
        {
            icovx = solve(covx + 0.01*diag(p))
        }else
        {
            icovx = solve(covx)
        }

        mvar = as.logical(isnanx[misobs[inn],])
        xo = x[misobs[inn],!mvar]

        # estimate missing part of x
        x[misobs[inn], mvar] =
            mx[mvar] - solve(icovx[mvar,mvar]) %*% icovx[mvar,!mvar] %*% as.matrix(xo - mx[!mvar])

        complobs = c(complobs, misobs[inn])
        ncomplobs = ncomplobs + 1
    }

    x[sorth$ix,]
}
