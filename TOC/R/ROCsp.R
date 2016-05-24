.ROCsp <- function(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE) {

if(!is.null(nthres) & !is.null(thres)) stop("Enter nthres OR thres as input, not both at the same time")

# extract cell values from the boolean and index maps 
boolval <- getValues(boolean)
indval <- getValues(index)

# extract cell values from the mask map if given
if(!is.null(mask)) mask <- getValues(mask)

# calculate basic roc (toc) table
tocd <- .ROCnosp(indval, boolval, mask=mask, nthres=nthres, thres=thres, NAval=NAval, progress=progress, 
                 ones.bool=NULL, zeros.bool=NULL)

return(tocd) 
}
