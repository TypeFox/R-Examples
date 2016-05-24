randEnv <- local({
  env <- new.env(parent=emptyenv())
  env[["halton.bases"]] <- c(2,3,5,7,11,13,17,19,23,29,31,37)
  env[["halton.perm"]] <- 0:127;
  env
})

# Get the next number from a <dimension>-dimensional
# Halton sequence
halton.getNext <- function(dimension)
{
  n <- randEnv[["halton.state"]] %/% dimension;
  b <- randEnv[["halton.bases"]][(randEnv[["halton.state"]] %% dimension) + 1];
  perm <- randEnv[["halton.perm"]];
  x <- 0;
  e <- -1
  while(n>0){
    x <- x + perm[n%%b+1]*b^e;
    e <- e - 1;
    n <- n %/% b;
  }
  randEnv[["halton.state"]] <- randEnv[["halton.state"]] + 1;
  return(x);
}

# Sample <n> sets of <dimension>-dimensional values
# from a Halton sequence with a random starting point.
halton.sample <- function(n,dimension)
{
  # random initialization
  randEnv[["halton.state"]] <- 0
  return(sapply(1:(n*dimension),function(x)halton.getNext(dimension)))
}
