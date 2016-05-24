co_sparse <- function(n, a, i, j, index, type, use)
{
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.integer(i))
    storage.mode(i) <- "integer"
  if (!is.integer(j))
    storage.mode(j) <- "integer"
  
  use <- check_use(use)
  if (use == "everything")
  {}
  else if (use == "all.obs")
  {
    if (anyNA(a))
      stop("missing observations in covar/pcor/cosine")
  }
  ### TODO
  # else if (use == "complete.obs")
  # {
  #   if (anyNA(x))
  #   {
  #     out <- naomit_coo(a, i, j)
  #     a <- out[[1]]
  #     i <- out[[2]]
  #     j <- out[[3]]
  #   }
  # }
  else
    stop("unsupported 'use' method")
  
  .Call(R_co_sparse, as.integer(n), a, i, j, as.integer(index), as.integer(type))
}



csc_to_coo <- function(row_ind, col_ptr)
{
  .Call(R_csc_to_coo, row_ind, col_ptr)
}
