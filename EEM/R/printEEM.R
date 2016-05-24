#' Print EEM
#' 
#' Print EEM
#' 
#' @param x EEM class object
#' @param ... arguments for \code{print} function
#' 
#' @examples
#' data(applejuice)
#' print(applejuice)
#' 
#' @export
print.EEM <- function(x, ...) {
    dimension <- sapply(x, function(y) {DIM <- dim(y)
                                        paste(DIM[1], DIM[2], sep = "x") })
    EX_range <- sapply(x, function(y) {RANGE <- range(as.numeric(colnames(y))) 
                                        paste(RANGE[1], RANGE[2], sep = "~")})
    EM_range <- sapply(x, function(y) {RANGE <- range(as.numeric(rownames(y)))
                                        paste(RANGE[1], RANGE[2], sep = "~")})
    
    # assign sName to Null if 
    if (!is.null(names(x))) sName <- names(x) else sName <- rep("none", length(x))
    
    # assemble table
    Table <- data.frame( Sample_name = sName,       
                        Dimension = dimension,
                        Excitation_range = EX_range,
                        Emission_range = EM_range, 
                        row.names = NULL)
    print(Table)
}