#
#Wolfgang Konen
#Cologne University of Applied Science
#fnArchive.R
#
# Wrapper fnArchiveF(x) for a function fn(x) returning res=c(obj,g1,...,gn).
# Extra feature: all x are stored in (ncall x length(x))-matrix soluArchive 
# and all res are stored in (ncall x length(res))-matrix funcArchive.
# ncall is the number of calls to fnArchiveF since the last resetSoluArchive().



#' construct wrapper \code{fnArchiveF} for a function \code{fn}
#' 
#' @param fn  function with vector argument x, should return a vector (obj,g1,g2,...,gn)
#' 
#' @return \code{fnArchiveF}, a function which is a wrapper for \code{fn} and stores as a side effect the argument x and 
#'        the return value of \code{fn} in soluArchive and funcArchive, resp. These archives can be 
#'        retrieved at any time with \cr 
#'            \code{(environment(fnArchiveF))$getSoluArchive()}  and \cr
#'            \code{(environment(fnArchiveF))$getFuncArchive()}  
#' 
#' @keywords internal 
#' 
fnArchiveFactory <- function(fn) {
  ffx <- fn; 
  soluArchive <- NULL;
  funcArchive <- NULL;
  #' reset soluArchive and funcArchive
  #' @keywords internal 
  resetSoluArchive <- function() {
    soluArchive <<- NULL
    funcArchive <<- NULL
  }
  
  getSoluArchive <- function() { return(soluArchive) }
  getFuncArchive <- function() { return(funcArchive) }
  
  #' execute \code{res=fn(x)}, store x in soluArchive and res in funcArchive
  #' 
  #' @param x
  #' @keywords internal 
  fnArchiveF <- function(x) {
    res <- ffx(x);
    if (!is.null(soluArchive)) testit::assert("Dimension mismatch soluArchive <--> x",length(x)==ncol(soluArchive))
    if (!is.null(funcArchive)) testit::assert("Dimension mismatch funcArchive <--> res",length(res)==ncol(funcArchive))
    soluArchive <<- rbind(soluArchive,x)
    funcArchive <<- rbind(funcArchive,res)
    return (res)  
  }
  
  return(fnArchiveF);
}

 
# --- only for debug / example ---
#
testFnArchive <- function() {
  fn <- function(x) { x^2+3; }

  fnArchiveF <- fnArchiveFactory(fn);
  
  fnArchiveF(3);
  fnArchiveF(4);
  
  I <- (environment(fnArchiveF))$getSoluArchive()
  R <- (environment(fnArchiveF))$getFuncArchive()
  print(I)
  print(R)
} 

#testFnArchive()
#
# --- only for debug / example ---

