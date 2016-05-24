#' Fold EEM matrix into a list
#' 
#' Fold EEM matrix into a list 
#' 
#' @param EEM_uf Unfolded EEM matrix where columns are wavelength condition and rows are samples.
#' It should have corresponding column names (formatted as EX###EM###) and row names. 
#' @param name optional for data.frame input to specify the sample names
#' @param ... arguments for other methods
#' 
#' @return EEM a list containing EEM/EEM data
#' 
#' @examples
#' data(applejuice)
#' applejuice_uf <- unfold(applejuice) # unfold list into matrix
#' applejuice_uf_norm <- normalize(applejuice_uf) # normalize matrix
#' drawEEM(fold(applejuice_uf_norm), 1) # visualize normalized EEM
#' 
#' @export
#' 
#' @importFrom reshape2 acast
fold <- function(EEM_uf, ...) UseMethod("fold")

#' @rdname fold
#' @export
fold.matrix <- function(EEM_uf, ...){
    
  # information from EEM_uf
  sName <- rownames(EEM_uf)
  N <- dim(EEM_uf)[1]
  var <- colnames(EEM_uf)
  ex <- getEX(var)
  em <- getEM(var)
  
  # add data into list
  EEM <- list()
  for (i in 1:N){
    data <- EEM_uf[i,]
    dataFrame <- data.frame(x = ex, y = em, z = as.numeric(data))
    EEM[[i]] <- as.matrix(acast(dataFrame, y ~ x, value.var = "z"))
  }
    
  # return 
  names(EEM) <- sName
  class(EEM) <- "EEM"
  return(EEM)
}

#' @describeIn fold fold unfolded data.frame
#' @export
fold.data.frame <- function(EEM_uf, name = NULL, ...){
    
    # turn into matrix
    EEM_uf <- as.matrix(EEM_uf)
    if (!is.null(name)) rownames(EEM_uf) <- name
    
    # use fold.matrix
    fold.matrix(EEM_uf)
}

#' @rdname fold
#' @export
fold.numeric <- function(EEM_uf, ...){
    
    # information from EEM_uf
    var <- names(EEM_uf)
    ex <- getEX(var)
    em <- getEM(var)
    
    # add data into list
    EEM <- list()
    data <- EEM_uf
    dataFrame <- data.frame(x = ex, y = em, z = as.numeric(data))
    EEM[[1]] <- as.matrix(acast(dataFrame, y ~ x, value.var = "z"))
    
    # return 
    class(EEM) <- "EEM"
    return(EEM)
}

#' @export
fold.tbl_df <- function(EEM_uf, name = NULL, ...){
    EEM_uf <- as.data.frame(EEM_uf)
    fold.data.frame(EEM_uf, name = name)
}

#' @export
fold.tbl <- function(EEM_uf, name = NULL, ...){
    EEM_uf <- as.data.frame(EEM_uf)
    fold.data.frame(EEM_uf, name = name)
}