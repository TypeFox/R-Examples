##' Concatenate functions
##'
##' @param ... functions
##' @author David Hajage
##' @keywords internal
funs2fun <- function(...) {
  fnames <- as.character(match.call()[-1])
  fs <- list(...)
  n <- length(fs)
  function(x, ...) {
    results <- NULL
    args <- list(...)
    namesargs <- names(args)
    for (i in 1:n) {
      func <- match.fun(fs[[i]])
      forms <- formals(func) # Pour min et max (et les autres
                             # primitives), il faudrait mettre
                             # 'formals(args(func))'. Le problème est
                             # que min et max retourne le minimum de
                             # tout ce qui n'est pas 'na.rm', donc si
                             # je met un autre argument (genre probs =
                             # 1/3), min et max prennent en compte sa
                             # valeur, d'où surprises... Je préfère
                             # laisser comme ça.
      namesforms <- names(forms)
      if (all(namesforms != "...")) {
        finalargs <- c(list(x = x), args[namesargs %in% namesforms])
      } else {
        finalargs <- c(list(x = x), args)
      }
      tmp <- do.call(func, finalargs)
      names(tmp) <- paste(fnames[i], names(tmp))
      results <- c(results, tmp)
    }
    results
  }
}
