#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[getDesignMatrix.R] by DSB Don 17/06/2010 14:57 (CEST)>
##
## Description:
## Extract the design matrix of the first element (== model) of a BayesMfp object.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   center the (non-intercept) columns
## 04/09/2008   argument center in getFpTransforms is now used, correct centering of uc columns
## 05/09/2008   use xCentered attribute from BayesMfp object
## 09/11/2008   add parameter to control if fixed columns (intercept) should be in the
##              return matrix
## 17/06/2010   do the centering optionally and only at the end and give the
##              column means back as an attribute of the design matrix. This is
##              important for the prediction of new data points, where exactly
##              the same shifts must be used in the construction of the new
##              design matrix! So we cannot just center the new design matrix
##              with its own column means but we must use the column means of
##              the old design matrix.
#####################################################################################

getDesignMatrix <- function (x, # a valid BayesMfp-Object of length 1 (otherwise only first element
                                # recognized)
                             fixedColumns=TRUE, # return the fixed columns
                                        # inside the matrix or not?
                             center=TRUE # do the centering?
                             )
{
    full <- attr (x, "x")
    
    inds <- attr (x, "indices")
    powers <- x[1][[1]]$powers
    ucSet <- x[1][[1]]$ucTerms

    nFix <- length (inds$fixed)         
    stopifnot(identical(nFix, 1L))

    ucColInds <- inds$uc %in% ucSet
    nUc <- sum (ucColInds)

    nFp <- length (unlist (powers))

    ## reserve space for return matrix
    nColumns <-
        if(fixedColumns)
            nFix + nUc + nFp
        else
            nUc + nFp
    
    ret <- matrix (nrow = nrow (full),
                   ncol = nColumns)
    retColnames <- character (ncol (ret))
    col <- 0L                           # invariant: already col columns written

    if(fixedColumns)
    {
        ## fixed columns
        new <- full[, inds$fixed, drop = FALSE]
        newInds <- col + seq_along (inds$fixed)
        
        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)
        
        col <- col + nFix
    }
    
    ## fp part
    for (i in seq_along (inds$bfp)){
        pi <- powers[[i]]
        if (len <- length (pi)) {       # if there is at least one power
            new <- getFpTransforms (full[, inds$bfp[i], drop = FALSE], pi, center=FALSE)
            newInds <- col + seq_along (pi)

            ret[, newInds] <- new
            retColnames[newInds] <- colnames (new)

            col <- col + len
        }
    }

    ## uc part
    if (length (ucSet)){
        new <- full[, ucColInds, drop = FALSE]
        newInds <- col + seq_len (nUc)

        ret[, newInds] <- new
        retColnames[newInds] <- colnames (new)

        col <- col + nUc
    }

    rownames (ret) <- rownames (full)
    colnames (ret) <- retColnames

    ## only once center if that was wished
    shifts <-
        if(center)
        {
            colMeans(ret)
        }
        else
        {
            ## no shifting at all
            rep.int(0, times=nColumns)
        }

    ## intercept column is not shifted in any case
    if(fixedColumns)
        shifts[1L] <- 0

    ## return the (optionally centered) design matrix and the used shifts
    return(structure(sweep(ret, 2, shifts),
                     shifts=shifts))
}



