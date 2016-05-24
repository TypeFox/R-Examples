##
##  f i n d . R  Finds indices of nonzero elements
##


finds <- function(v) 
    which( if (is.logical(v)) v else v != 0 )
