myBsample <- function(val, original_sample_size, n.boot, prob = NULL) {
      # val is a vector to resample the sample size is original_sample_size*n.boot
      if (length(val) == 0) {
            res <- NULL
      } else if (length(val) == 1) {
            # problems with weighted sample if only one term
            res <- matrix(rep(val, original_sample_size * n.boot), n.boot,
                          original_sample_size, byrow = TRUE)
      } else if (length(val) > 1) {
            res <- matrix(sample(val, original_sample_size * n.boot,
                                 replace = TRUE, prob = prob), n.boot,
                          original_sample_size, byrow = TRUE)
      }
      res
}
