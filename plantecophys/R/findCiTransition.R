findCiTransition <- function(photofun, ...){
  O <- function(Ci, photofun, ...){
    x <- photofun(Ci=Ci, ...)
    Ac <- x$Ac - x$Rd
    Aj <- x$Aj - x$Rd
    (Ac - Aj)^2
  }
  o <- optimize(O, c(25,1500), photofun=photofun, ...)
  return(o$minimum)
}
