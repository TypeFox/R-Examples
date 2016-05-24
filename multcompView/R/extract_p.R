#' Extracts p-values
#' 
#' For a given object it will look for the column named "p adj" or "difference"
#' and extract its value mantaining its names
#' 
#' 
#' @aliases extract_p extract_p.default extract_p.TukeyHSD extract_p.mc
#' @param x A object that has p-values or logical values.
#' @return A named vector with p-values or logical values.
#' @author Luciano Selzer
#' @seealso \code{\link{multcompLetters}} \code{\link{multcompTs}}
#' @keywords manip array
#' @export
#' @examples
#' 
#' experiment <- data.frame(treatments = gl(11, 20, labels = c("dtl", "ctrl", "treat1", 
#'               "treat2", "treatA2", "treatB", "treatB2",
#'               "treatC", "treatD", "treatA1", "treatX")),
#'               y = c(rnorm(20, 10, 5), rnorm(20, 20, 5), rnorm(20, 22, 5), rnorm(20, 24, 5),
#'                rnorm(20, 35, 5), rnorm(20, 37, 5), rnorm(20, 40, 5), rnorm(20, 43, 5),
#'                rnorm(20, 45, 5), rnorm(20, 60, 5), rnorm(20, 60, 5)))
#' exp_tukey <- TukeyHSD(exp_aov <- aov(y  ~ treatments, data = experiment))
#' 
#' extract_p(exp_tukey)
#' 
#' require(pgirmess)
#' extract_p(kruskalmc(y ~ treatments, data = experiment))
#' 
#' 
"extract_p" <- function(x) {
  UseMethod("extract_p")
}

#' @export
#' @describeIn extract_p
#' 
extract_p.default <- function(x){
  ans <- x[ ,"p adj"]
  #To be sure that names are kept
  names(ans) <- rownames(x)
  ans
}

#' @export
#' @describeIn extract_p

extract_p.TukeyHSD <- function(x) {
   x <- lapply(x, extract_p.default)
   x
}

#' @export
#' @describeIn extract_p
#' 

extract_p.mc <- function(x){
  ans <- x[["dif.com"]][, "difference"]
  names(ans) <- rownames(x[["dif.com"]])
  ans
}

