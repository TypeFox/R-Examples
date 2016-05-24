##' Convert an IRT response matrix in the wide format to a data frame
##' in the long format. In the wide format each row corresponds to a
##' respondent and each column to an item.
##' 
##' @title Convert a wide-format matrix to long format
##' @param mat an integer IRT response matrix (i.e. a wide format)
##' @return a long-format data frame
tolong <- function(mat) {
    stopifnot(is.matrix(mat))
    nr <- nrow(mat)
    nc <- ncol(mat)
    data.frame(value = as.vector(mat),
               item = gl(nc, nr,
               labels = sprintf(paste("i%0", nchar(nc), "d", sep=''), seq_len(nc))),
               person = gl(nr, 1, length=length(mat),
               labels = sprintf(paste("p%0", nchar(nr), "d", sep=''), seq_len(nr))))
}

dendrify <- function(mat, cmx) {
    ff <- factor((m1 <- tolong(mat))[["value"]])
    stopifnot(is.matrix(cmx),
              (nr <- nrow(cmx)) == length(levels(ff)),
              (nc <- ncol(cmx)) < nr)
    subset(within(data.frame(value = as.vector(cmx[as.integer(ff), ]),
                             item = rep(m1$item, nc),
                             person = rep(m1$person, nc),
                             node = gl(nc, nrow(m1),
                             labels = sprintf(paste("node%0", nchar(nc), "d", sep=''),
                             seq_len(nc)))),
                  sub <- item:node), !is.na(value))
}

