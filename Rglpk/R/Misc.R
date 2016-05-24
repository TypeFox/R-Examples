## Miscellaneous functions used in this package

## Get the right representation of integer vector.
## This function takes an index vector and returns a vector of logicals
## with length 'n'.

#glp_integers <-
#function(x, n)
#{
#  if(!all(x <= n))
#stop("Indices must not exceed the number of objective coefficients.")
#  out <- logical(n)
#  out[x] <- TRUE
#  out
#}

print.MP_data_from_file <- function(x, ...){
  if(!inherits(x, "MP_data_from_file"))
     stop("'x' must be of class 'MP_data_from_file'")
  if(attr(x, "n_integer_vars") > 0L){
    writeLines(paste("A mixed integer linear program with",
                     attr(x, "n_objective_vars"), "objective variables,"))
    writeLines(paste(attr(x, "n_integer_vars"), "are integer and",
                     attr(x, "n_binary_vars"), "of which are binary variables."))
    writeLines(paste("This problem has", attr(x, "n_constraints"),
                     "constraints with", attr(x, "n_nonzeros"),
                     "non-zero values in the constraint matrix."))
  } else{
    writeLines(paste("A linear program with", attr(x, "n_objective_vars"), "objective variables."))
    writeLines(paste("This problem has", attr(x, "n_constraints"),
                     "constraints with", attr(x, "n_nonzeros"),
                     "non-zero values in the constraint matrix."))
  }
}
