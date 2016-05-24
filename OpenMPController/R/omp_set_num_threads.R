omp_set_num_threads <- function(n)
  {
    .C("omp_set_num_threads_ptr", as.integer(n))
  }
