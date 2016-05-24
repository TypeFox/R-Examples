#' Convert an object to an mpoly
#' 
#' mpoly is the most basic function used to create objects of class
#' mpoly.
#' 
#' @param x an object of class lm
#' @param ... additional arguments to pass to methods
#' @return the object formated as a mpoly object.
#' @author David Kahle \email{david.kahle@@gmail.com}
#' @seealso \code{\link{mp}}
#' @export
#' @examples
#' 
#' library(plyr)
#' 
#' n <- 101
#' s <- seq(-5, 5, length.out = n)
#' 
#' # one dimensional case
#' df <- data.frame(x = s)
#' df <- mutate(df, y = -x^2 + 2*x - 3 + rnorm(n, 0, 2))
#' with(df, plot(x, y))
#' mod <- lm(y ~ x + I(x^2), data = df)
#' (p <- as.mpoly(mod))
#' f <- as.function(p)
#' lines(s, f(s), col = "red")
#' 
#' 
#' 
#' # one dimensional case with ggplot2
#' library(ggplot2); theme_set(theme_bw())
#' 
#' qplot(x, y, data = df) 
#' qplot(x, y, data = df) +
#'   stat_function(fun = f, colour = "red")
#' 
#' 
#' # two dimensional case with ggplot2
#' 
#' df <- expand.grid(x = s, y = s)
#' df <- mutate(df, z = x^2 - y^2 + 2 * x*y + rnorm(n^2, 0, 3))
#' qplot(x, y, data = df, geom = "raster", fill = z)
#' mod <- lm(z ~ x + y + I(x^2) + I(y^2) + I(x*y), data = df)
#' p <- as.mpoly(mod)
#' f <- as.function(p)
#' df$fit <- apply(df[,c("x","y")], 1, f)
#' qplot(x, y, data = df, geom = "raster", fill = fit)
#' qplot(x, y, data = df, geom = "raster", fill = z - fit) # residuals
#' 
#' 
#' 
#' 
#' 
#' 
#' 
as.mpoly <- function(x, ...) UseMethod("as.mpoly")



#' @export  
as.mpoly.default <- function(x, ...) 
  stop("object not supported.  see ?as.mpoly for details.")




#' @export  
as.mpoly.lm <- function(x, ...){
  coefs <- coefficients(x)
  coef_names <- names(coefs)
  coef_names[coef_names == "(Intercept)"] <- 1
  I_ndcs <- which(str_detect(coef_names, "I([0-9a-zA-Z]*)"))
  if(length(I_ndcs) > 0){
    coef_names[I_ndcs] <- sapply(as.list(coef_names[I_ndcs]), 
      function(s) {
          str_sub(s, 3, -2)
      }
    )
  }
  coef_names  <- str_replace_all(coef_names, " \\* ", " ")
  mp_str <- paste(coefs, coef_names, sep = " ", collapse = " + ")
  mp(mp_str)
}





#' @export  
as.mpoly.polynomial <- function(x, indeterminate = "x", ...){
  as.mpoly.numeric(unclass(x), indeterminate)
}





#' @export  
as.mpoly.numeric <- function(x, indeterminate = "x", ...){
  n <- length(x)
  
  ## make list and populate
  p <- list()
  p[[1]] <- c(coef = x[1])
  if(n > 1) for(deg in 1:(n-1)){
    v <- c(deg, x[deg+1])
    names(v) <- c(indeterminate, "coef")
    p[[deg+1]] <- v
  }
  
  ## clean out zeros
  p <- Filter(function(v) v[["coef"]] != 0, p)
  
  ## clean in case everything is removed
  if(length(p) == 0) p <- list(c(coef = 0))
  
  ## class and out
  class(p) <- "mpoly"
  p
}


