### rescale function

## dataFlag - flag, whether the the data is plotted or not
## rescaleFlag

# function returns the list of rescaling arguments to be passed on the
# corresponding diagnostic function

setMethod("rescaleFunction", signature(L2Fam="ANY"),
   function(L2Fam, ...){
      scaleList <- list(scaleX = substitute(FALSE)
                        ,scaleY = substitute(FALSE))
    return(scaleList)}
)
