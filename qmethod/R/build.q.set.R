build.q.set <- function(q.concourse, q.sample, q.distribution) {
  q.sample <- as.character(q.sample) #  just to be safe

  # Validate input =============================================================
  if (!is.matrix(q.concourse)) {
  	stop("The input specified for q.concourse is not a matrix.")
  }
  if (!is.vector(q.distribution)) {
    stop("The input specified for q.distribution is not a matrix.")
  }
  if (!is.vector(q.sample)) {
   stop("The input specified for q.sample is not a vector.")
  }
  if (length(q.sample) != sum(q.distribution)) { #  test if sums are equal
    stop(
      paste(
        "There are",
        length(q.sample),
        "items in your q-sample, but",
        sum(q.distribution),
        "entries expected in the q-distribution",
        sep=" "
      )
    )
  }
  missing.in.concourse <- !q.sample %in% rownames(q.concourse)
  if (any(missing.in.concourse)) {  # if any missing, stop
    stop(
      paste(
         "There are item handles in your sample not defined in the concourse:",
        q.sample[missing.in.concourse],
        sep=" "
      )
    )
  }

  # Subset the concourse =================================================
  q.set <- q.concourse[q.sample,]  # only add sampled rows from concourse
  q.set <- as.matrix(q.set)
  message(paste("Build a q.set of", nrow(q.set), "items."))
	return(q.set)
}
