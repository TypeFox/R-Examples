CTM_registry <- list(CTM_VEM.fit = c("VEM", "CTM_VEM", "CTM_VEM.fit"))

CTM <- function(x, k, method = "VEM", control = NULL, model = NULL, ...) {
  if (is(x, "DocumentTermMatrix")) {
    if (!any(attr(x, "weighting") %in% c("term frequency", "tf"))) {
      stop("The DocumentTermMatrix needs to have a term frequency weighting")
    }
  } else if (!is(x, "simple_triplet_matrix")) {
    x <- as.simple_triplet_matrix(x)
  }
  if (!all.equal(x$v, as.integer(x$v)))
    stop("Input matrix needs to contain integer entries")
  if (!all(row_sums(x) > 0))
    stop("Each row of the input matrix needs to contain at least one non-zero entry")
  mycall <- match.call()
  
  if (!is.null(model)) {
    x <- match_terms(x, model)
    k <- model@k
  }

  if (as.integer(k) != k || as.integer(k) < 2) stop("'k' needs to be an integer of at least 2")  

  if(missing(method) && !missing(model))
    method <- paste(class(model), "fit", sep = ".")
  if(!is.function(method)) {
    MATCH <- which(sapply(CTM_registry, function(x) length(grep(tolower(method), tolower(x)))) > 0)
    if (!length(MATCH) == 1)
      stop("'method' not specified correctly")
    method <- get(names(CTM_registry)[MATCH])
  }

  method(x, k, control, model, mycall, ...)
}

CTM_VEM.fit <- function(x, k, control = NULL, model = NULL, call, ...) {
  control <- as(control, "CTM_VEMcontrol")
  if (length(control@seed) != control@nstart)
    stop(paste("Need ", control@nstart, " seeds", sep = ""))
  if (control@initialize == "random") control@initialize <- "rand"
  else if (control@initialize == "seeded") control@initialize <- "seed" 
  else if (control@initialize == "model") {
    if (!is(model, "CTM")) stop("Need a model of class 'CTM' for initialization")
  }
  if (is(model, "CTM")) control@initialize <- "model"
  result_dir <- path.expand(paste(control@prefix, "-ctm", sep = ""))
  if (control@save) dir.create(result_dir, showWarnings = FALSE)
  if (!control@estimate.beta) control@em@iter.max <- -1L

  obj <- vector("list", control@nstart)
  for (i in seq_len(control@nstart)) {
    control_i <- control
    control_i@seed <- control@seed[i]
    obj[[i]] <- .Call("rctm",
                      ## simple_triplet_matrix
                      as.integer(x$i),
                      as.integer(x$j),
                      as.numeric(x$v),
                      as.integer(x$nrow),
                      as.integer(x$ncol),                 
                      ## CTMcontrol
                      control_i,
                      ## number of topics
                      as.integer(k),
                      ## directory for output files
                      result_dir,
                      ## initial model
                      model,
                      PACKAGE = "topicmodels")
    obj[[i]]@gamma <- cbind(exp(obj[[i]]@gamma), 1)
    obj[[i]]@gamma <- obj[[i]]@gamma/rowSums(obj[[i]]@gamma)
    obj[[i]] <- new(class(obj[[i]]), obj[[i]], call = call, control = control_i,
                    documents = x$dimnames[[1]], terms = x$dimnames[[2]], n = as.integer(sum(x$v)))
  }
  if (control@best) {
    MAX <- which.max(sapply(obj, logLik))
    if (length(MAX)) {
      obj <- obj[[MAX]]
    } else warning("no finite likelihood")
  }
  obj

}
