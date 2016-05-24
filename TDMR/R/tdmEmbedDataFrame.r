######################################################################################
##' Embed columns in a data frame.
##'
##' Embeds the columns named \code{cols} in the data frame \code{x} into a space of dimension \code{dimension}.
##'
##' @param x The data frame containing the columns to embed.
##' @param cols [NULL] A vector with a list of names of the columns to embed. If NULL, embed all columns of x.
##' @param dimension [1] The additional dimensions to generate when embedding.
##' @param take [1] Take from these dimensions only every take'th column.
##' @return The data frame, augmented with embedded columns, shortended by \code{dimension} rows.
######################################################################################
tdmEmbedDataFrame <- function(x, cols = NULL, dimension = 1, take = 1) {
  d <- dimension + 1
  if(d <= 1 || d >= nrow(x))
    stop("dimension must be greater than 0 and smaller than nrow(x)");
  if(take <= 0 || take > dimension)
    stop("parameter take must be greater than 0 and not larger than dimension");
    
  columns <- if (missing(cols)) names(x) else cols

  # build a list of matrices, one for each element A,B,... in columns, where the 1st matrix 
  # has column names A.P1, A.P2, ..., the 2nd matrix has column names B.P1, B.P2, ... and so on. 
  embeddings <- lapply(columns, function(column) {
                                  embedded <- as.matrix(embed(x[[column]], d)[,-1])
                                  colnames(embedded) <- paste(column, ".P", 1:ncol(embedded), sep = "")
                                  embedded
                                });
                                
  # select only every take'th column from each of the matrices (does nothing if take==1)
  embeddings <- lapply(embeddings, function(n) { 
                                  tcol <- (1:floor(dimension/take))*take;
                                  em <- as.matrix(n[,tcol]); 
                                  colnames(em) <- colnames(n)[tcol]; 
                                  em; 
                                  });
                                
  data.frame(list(x[-1:-(d-1),], embeddings)); 
}


