## kwak_get_n_xoi.R

get_n_xo <- function(sdat)
{
    n_ind <- nrow(sdat)
    n_pos <- ncol(sdat)

    if(!is.matrix(sdat) || n_ind==0 || n_pos==0)
        stop("Input must be a matrix with at least one row and one column.")

    output <- .C("R_get_N_xo",
                 as.integer(n_ind),
                 as.integer(n_pos),
                 as.integer(sdat),
                 n_xo=as.integer(0),
                 PACKAGE="xoi" )
    output$n_xo
}
