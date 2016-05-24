#' Multiple code chunks
#'
#' This function creates multiple code chunks from a function and along arguments marked with a star (*).
#' Each of these special arguments is a list.
#' The nth code chunk will use the nth element of each marked list (recycled if necessary) as argument.
#' 
#' @param FUN the function to use for the chunks.
#' @param arg.names a character vector giving the argument names of \code{FUN} to set.
#' @param arg.values a vector giving the values or object names to assign to
#' each argument given with \code{arg.names} (they must match in order).
#' Object names must be backquoted or quoted.
#' Lists names marked with a star (e.g. \code{"*L1"} for \code{L1})
#' indicate their elements will be used sequentially in chunks. 
#' @param type the type of chunk to produce. Can be "\code{block}", "\code{inline}" or "\code{none}".
#' @param echo logical indicating whether to include R source code in the result.
#' @param warning logical indicating whether to print warnings in the result.
#' @param error logical indicating whether to stop on errors.
#' @param message logical indicating whether to print messages in the result.
#' @param fig.width,fig.height numeric value setting width and height of the plots in inches.
#' @param fig.align character string setting the alignment of the plots. Can be "\code{left}", "\code{right}" and "\code{center}".
#' @param options a character string to specify the knitr options.
#' This will overwrite the options set with the other arguments.
#' 
#' @return a character vector of R code chunks which can be evaluated by \pkg{knitr}.
#' 
#' @export
chunkerize <- function(FUN, arg.names, arg.values, type = "block",
                       echo = FALSE, warning = FALSE, error = FALSE, message = TRUE,
                       fig.width = 4, fig.height = 4, fig.align = "center",
                       options = NULL){
  
  if(!is.character(FUN)){
    FUN <- deparse(substitute(FUN))
  }
  
  type <- match.arg(type, c("block", "inline", "none"))
  
  arg.star <- sapply(arg.values, function(x) substr(x, 1, 1) == "*" & nchar(x) > 1)
  arg.values[arg.star] <- substr(arg.values[arg.star], 2, 100000L)

  star.len <- sapply(arg.values[arg.star], function(x) length(eval(parse(text = x))))
  
  if(max(star.len) != min(star.len)){
    stop("Rolling lists must have same lengths.")
  }
  
  arg.values.mat <- matrix(arg.values, nrow = max(star.len), ncol = length(arg.values), byrow = TRUE)
  idx <- paste("[[", seq(1:max(star.len)), "]]", sep = "")

  arg.values.mat[, which(arg.star)] <- paste(arg.values.mat[, which(arg.star)], idx, sep = "")

  arg.comp <- paste(arg.names, " = ", t(arg.values.mat), sep = "")
  arg.comp <- matrix(arg.comp, nrow = max(star.len), ncol = length(arg.values), byrow = TRUE)
  arg.comp <- apply(arg.comp, 1, paste, sep = "", collapse = ", ")
  
  res <- paste(FUN, "(", arg.comp, ")", sep = "")
  
  if(type == "inline"){
    res <- paste("`r ", res, "`", sep = "")
  }
  
  if(type == "block"){
    if(is.null(options)){
      options <- paste("echo=", echo,
                       ", warning=", warning,
                       ", error=", error,
                       ", message=", message,
                       ", fig.width=", fig.width,
                       ", fig.height=", fig.height,
                       ", fig.align='", fig.align, "'", sep = "")
    }
    res <- paste("\n```{r ", options, "}\n\t", res, "\n```\n", sep = "")
  }
  
  return(res)
}

