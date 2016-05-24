### S3 method for conjugate transpose, Hermitian transpose, or
### Hermitian conjugate.

### This is a bad idea.
# t.complex <- function(x){
#   if(!is.complex(x)){
#     stop("x is not in complex.")
#   }
#   if(!is.matrix(x)){
#     x <- as.matrix(x)
#   }
#   Conj(t.default(x))
# } # End of t.complex().

H <- function(x){
  if(!is.complex(x)){
    stop("x is not in complex.")
  }
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  Conj(t.default(x))
} # End of H().
