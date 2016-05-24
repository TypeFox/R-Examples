#' Standard measures of complexity for strings
#' 
#' Functions to compute different measures of complexity for strings: Entropy, Second-Order Entropy, and Change Complexity
#' 
#' @usage entropy(string)
#' 
#' entropy2(string)
#' 
#' change_complexity(string)
#' 
#' @param string \code{character} vector containing the to be analyzed strings (can contain multiple strings for the entropy measures).
#' 
#' @return \code{numeric}, the complexity of the string. For \code{entropy} and \code{entropy2} of the same length as \code{string}. \code{change_complexity} currently only works with inputs of length 1.
#' 
#' @details For users who need advanced functions, a comprehensive package computing various versions of entropy estimators is available \pkg{entropy}. For users who just need first and second-order entropy and which to apply them to short string, the \pkg{acss} package provides two functions: \code{entropy} (first-order entropy) and \code{entropy2} second-order entropy.
#' 
#' Change complexity (\code{change_complexity}) assesses cognitive complexity or the subjective perception of complexity of a binary string. It has been comprehensively defined by Aksentijevic and Gibson (2012). Although the algorithm will work with any number of symbols up to 10, the rationale of Change Complexity only applies to binary strings. 
#' 
#' @references Aksentijevic & Gibson (2012). Complexity equals change. \emph{Cognitive Systems Research}, 15-17, 1-16.
#' 
#' @name entropy
#' @aliases entropy entropy2 change_complexity
#' @export entropy entropy2 change_complexity
#' 
#' @example examples/examples.standard_complexity.R
#' 
#' 

######### Entropy
# returns the entropy of a given string.
# there is a package already : http://cran.r-project.org/web/packages/entropy/entropy.pdf
# but it is too complicated for our purpose.
entropy <- function(string){
  check_string(string)
  splitted <- strsplit(string,"")
  y <- lapply(splitted, function(x) as.vector(table(x)))
  l <- nchar(string)
  y <- mapply(function(x, y) x/y, y, l, SIMPLIFY=FALSE)
  names(y) <- string
  vapply(y, function(x) -sum(x*log2(x)), 0)
}

########## Second order entropy
# There are different ways to compute second order entropy. Because we have small strings, I prefer using a sliding window of 2 symbols.
entropy2 <- function(string){
  check_string(string)
  l <- nchar(string)
  if (any(l < 2)) stop("length of strings need to be > 1.")
  splitted <- strsplit(string,"")
  new.string <- lapply(splitted, function(x) paste0(x[-length(x)], x[-1]))
  y <- lapply(new.string, function(x) as.vector(table(x)))
  l <- vapply(new.string, length, 0)
  y <- mapply(function(x, y) x/y, y, l, SIMPLIFY=FALSE)
  names(y) <- string
  vapply(y, function(x) -sum(x*log2(x)), 0)  
}


# ref : Aksentijevic & Gibson (2012) Complexity equals change, Cognitive Systems Research, 15-17, 1-16
## change complexity for a binary string s

change_complexity <- function(string) {
  check_string(string)
  vapply(string, .change_complexity, 0)
}

.change_complexity <- function(string){
  #browser()
  l <- nchar(string)
  if (any(l<3)){stop("length of strings need to be > 2.")}
  #splitted <- strsplit(string,"") 
  
  s=strsplit(string,"")[[1]]
  m=matrix(rep(0,(l-1)^2),ncol=l-1)
  
  for (j in 1:(l-1)){
    if  (s[j]!=s[j+1]) {m[1,j]=1}
  }
  
  for (i in 2:(l-1)){
    for (j in 1:(l-i)){
      if  (m[(i-1),j]!=m[(i-1),j+1]) {m[i,j]=1}
    }
  }
  
  # print(m)
  
  p=rep(0,l-1)
  w=rep(0,l-1)
  for (i in 1:l-1){p[i]=sum(m[,i])
                   w[i]=1/(l-i)}
  C=sum(p*w)
  names(C) <- string
  return(C)
}

