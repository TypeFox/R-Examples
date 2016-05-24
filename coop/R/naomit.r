naomit_mat <- function(x)
{
  .Call(R_fast_naomit, x)
}

# naomit_coo <- function(a, i, j)
# {
#   .Call(R_naomit_coo, a, i, j)
# }

### TODO unify
naomit <- function(x)
{
  if (is.numeric(x) && is.matrix(x))
    naomit_mat(x)
  else
    stop("unsupported type in naomit")
}
