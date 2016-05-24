## identify_xo : identify crossover. 2nd version for speed by C coding
#
# Input : data matrix.
# Return : xo information matrix.
# xomat[1,i] : observation index.
# xomat[2,i] : xo strating marker. (left)
# xomat[3,i] : xo ending marker.  (right)
# ob_ind : xo matrix index for each observation.
#
# ** have get_n_xo() as a sub function
# return : number of crossover happened

identify_xo <- function(sdat)
{
    n_ind <- nrow(sdat)
    n_pos <- ncol(sdat)

    if(!is.matrix(sdat) || n_ind==0 || n_pos==0)
        stop("Input must be a matrix with at least one row and one column.")

    n_xo <- get_n_xo(sdat)

    output <- .C("R_identify_xo",
                 as.integer(sdat),
                 as.integer(n_ind),
                 as.integer(n_pos),
                 as.integer(n_xo),
                 left = as.integer(rep(0,n_xo)),
                 right = as.integer(rep(0,n_xo)),
                 ind_id = as.integer(rep(0,n_xo)),
                 ob_ind = as.integer(rep(0,n_ind)),
                 PACKAGE="xoi" )
    xomat <- rbind(output$ind_id, output$left, output$right)
    return(list(xomat=xomat, ob_ind=output$ob_ind))
}
