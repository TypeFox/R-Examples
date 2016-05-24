#########################################################
# Class to store the matrix which describes the polynomial
#########################################################
setClass("PCEdesign", slots=c(.Data="matrix"))
#########################################################
# print method
# INPUT
# x : return of Structure()
print.PCEdesign <- function (x, all=FALSE, ...) {
  # oter le terme constant
  planx <- x[-1,, drop=FALSE]
  degree <- max(planx)
  if (all) {
    cat("Polynomial expression:\n")
  # inverser les colonnes
#  planx <- planx[, seq(ncol(planx), 1, -1), drop=FALSE]
    descr <- paste(rownames(x),collapse=" + ")

    cat(descr,"\n" )
    # Number of monomials per degree
    so <- apply(planx, 1, sum)
    for (i in 1:degree) {
        cat("Number of monomials of degree ", i,": ",
            length(so[so==i]), "\n", sep="")
    } # fin i
  } # fin all
  
  cat("Total number of monomials: ", nrow(planx), "\n")
  cat("Number of factors: ", ncol(planx), "\n")
  cat("Polynomial degree: ", degree, "\n")
  return(invisible())
} # end print.PCEdesign

#########################################################
# show method
# INPUT
#  return of Structure()
show.PCEdesign  <- function(object){
  print.PCEdesign(object)
    return(invisible())
} # end show.PCEdesign


# --------------------------------------
setMethod("show", signature(object="PCEdesign"),
          definition=show.PCEdesign)

                
  

         
