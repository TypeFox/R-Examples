#' @encoding UTF-8
#' @title Untable
#'
#' @description Method for recreate the data.frame out of a contingency table, i.e., converts from summarized data to long.
#' @param x The table object as a data.frame, table, or, matrix.
#' @param freq The column with count values.
#' @param rownames Row names to add to the data.frame.
#' @param \dots Extra parameters.
#'
#' @examples
#' if (interactive()) {
#' gss <- data.frame(
#' expand.grid(sex=c("female", "male"),
#' party=c("dem", "indep", "rep")),
#' count=c(279,165,73,47,225,191))
#'
#' print(gss)
#'
#' # Then expand it:
#' GSS <- untable(gss, freq="count")
#' head(GSS)
#' }
#' @export
`untable` <- function(x, ...){
  UseMethod("untable")
}
NULL


#' @rdname untable
#' @export
`untable.data.frame` <- function(x, freq = "Freq", rownames = NULL, ...){

  if(all(is.na(match(freq, names(x)))))
    stop(gettextf("Frequency column %s does not exist!", freq))

  res <- x[untable(x[,freq], type="as.numeric")[,], -grep(freq, names(x))]
  rownames(res) <- rownames

  return(res)
}
NULL



#' @rdname untable
#' @param dimnames Set dimnames of an object if require.
#' @param type The type of variable. If NULL, ordered factor is returned.
#' @param colnames Column names to add to the data.frame.
#' @export
`untable.default` <- function(x, dimnames=NULL, type = NULL, rownames = NULL, colnames = NULL, ...) {
  # coerce to table, such as also be able to handle vectors
  x <- as.table(x)
  if(!is.null(dimnames)) dimnames(x) <- dimnames
  if(is.null(dimnames) && identical(type, "as.numeric")) dimnames(x) <- list(seq_along(x))
  # set a title for the table if it does not have one
  # if(is.null(names(dimnames(x)))) names(dimnames(x)) <- ""
  # if(length(dim(x))==1 && names(dimnames(x))=="") names(dimnames(x)) <- "Var1"
  # replaced 26.3.2013
  for( i in 1:length(dimnames(x)) )
    if (is.null(names(dimnames(x)[i])) || names(dimnames(x)[i]) == "")
      if (length(dimnames(x)) == 1) names(dimnames(x)) <- gettextf("Var%s", i)
      else names(dimnames(x)[i]) <- gettextf("Var%s", i)

      res <- as.data.frame(expand.grid(dimnames(x))[rep(1:prod(dim(x)), as.vector(x)),])
      rownames(res) <- NULL
      if(!all(names(dimnames(x))=="")) colnames(res) <- names(dimnames(x))

      # return ordered factors, if wanted...
      if(is.null(type)) type <- "as.factor"
      # recycle type:
      if(length(type) < ncol(res)) type <- rep(type, length.out=ncol(res))

      for(i in 1:ncol(res)){
        if(type[i]=="as.numeric"){
          res[,i] <- as.numeric(as.character(res[,i]))
        } else {
          res[,i] <- eval(parse(text = gettextf("%s(res[,i])", type[i])))
        }
      }

      # overwrite the dimnames, if requested
      if(!is.null(rownames)) rownames(res) <- rownames
      if(!is.null(colnames)) colnames(res) <- colnames

      return(res)
}
NULL

