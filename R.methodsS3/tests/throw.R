library("R.methodsS3")

message("TESTING: throw()...")

rbern <- function(n=1, prob=1/2) {
  if (prob < 0 || prob > 1)
    throw("Argument 'prob' is out of range: ", prob)
  rbinom(n=n, size=1, prob=prob)
}

rbern(10, 0.4)
# [1] 0 1 0 0 0 1 0 0 1 0
tryCatch({
  rbern(10, 10*0.4)
}, error=function(ex) {
  print(ex)
})

message("TESTING: throw()...DONE")
