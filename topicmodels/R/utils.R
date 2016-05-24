setMethod("posterior", signature(object = "TopicModel", newdata = "missing"),
function(object, newdata, ...) {
  terms <- exp(object@beta)
  dimnames(terms) <- list(seq_len(object@k), object@terms)
  topics <- object@gamma
  dimnames(topics) <- list(object@documents, seq_len(object@k))
  list(terms = terms,
       topics = topics)
})

setMethod("posterior", signature(object = "TopicModel", newdata = "ANY"),
function(object, newdata, control = list(), ...) {
  if (!is(newdata, "simple_triplet_matrix"))  {
    newdata <- as.simple_triplet_matrix(newdata)
  }
  CLASS <- strsplit(class(object), "_")[[1]]
  control <- as(control, paste(class(object), "control", sep = ""))
  control@estimate.beta <- FALSE
  terms <- exp(object@beta)
  dimnames(terms) <- list(seq_len(object@k), object@terms)
  topics <- get(CLASS[1])(newdata, method = CLASS[2], 
         model = object, control = control)@gamma
  dimnames(topics) <- list(rownames(newdata), seq_len(object@k))
  list(terms = terms,
       topics = topics)
})

setGeneric("terms")
setGeneric("topics", function(x, ...) standardGeneric("topics"))

setMethod("topics", signature(x = "TopicModel"), function(x, k, threshold, ...) 
  most_likely(x, "topics", k, threshold, ...))

setMethod("terms", signature(x = "TopicModel"), function(x, k, threshold, ...) 
  most_likely(x, "terms", k, threshold, ...))

most_likely <- function(x, which = c("terms", "topics"), k, threshold, ...)
{
  which <- match.arg(which)
  labels <- if (which == "terms") x@terms else seq_len(x@k)
  post <- posterior(x)[[which]]
  if (missing(k)) {
    k <- ifelse(missing(threshold), 1, ncol(post))
  } else {
    k <- min(ncol(post), k)
  }
  if (!missing(threshold)) {
    most <- sapply(seq_len(nrow(post)), function(i) {
      index <- which(post[i,] > threshold)
      labels[index[index %in% order(post[i,], decreasing = TRUE)[seq_len(k)]]]
    }, ...)
  } else {
    most <- sapply(seq_len(nrow(post)), function(i) 
                    labels[order(post[i,], decreasing = TRUE)[seq_len(k)]], ...)
  }
  if (is(most, "matrix")) {
    colnames(most) <- if (which == "terms") paste("Topic", seq_len(ncol(most))) else x@documents
  } else {
    names(most) <- if (which == "terms") paste("Topic", seq_along(most)) else x@documents
  }
  return(most)
}

setGeneric("get_df", function(object, ...) standardGeneric("get_df"))

setMethod("get_df", signature(object="LDA_Gibbs"),
function(object, ...) length(object@beta))

setMethod("get_df", signature(object="LDA_VEM"),
function(object, ...) 
  as.integer(object@control@estimate.alpha) + length(object@beta))

setMethod("get_df", signature(object="CTM_VEM"),
function(object, ...) 
  (object@k - 1) * (object@k/2 + 1) + length(object@beta))

setMethod("logLik", signature(object="TopicModel"),
function(object, ...) {
  val <- sum(object@loglikelihood)
  attr(val, "df") <- get_df(object)
  attr(val, "nobs") <- object@Dim[1]
  class(val) <- "logLik"
  val
})

setMethod("logLik", signature(object="Gibbs_list"),
function(object, ...) sapply(object@fitted, logLik))

distHellinger <- function(x, y, ...) UseMethod("distHellinger")

distHellinger.default <- function(x, y, ...) 
{
  if (missing(y)) {
    x <- sqrt(x) 
    z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    for (k in seq_len(nrow(x)-1)) {
      z[k,-seq_len(k)] <- z[-seq_len(k),k] <- 
        sqrt(1/2 * colSums((t(x[-seq_len(k),,drop=FALSE]) - x[k,])^2))
    }
  } else {
    if(ncol(x) != ncol(y))
      stop("'x' and 'y' must have the same number of columns")
    x <- sqrt(x); y <- sqrt(y)
    z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (k in seq_len(nrow(y))) z[,k] <- sqrt(1/2 * colSums((t(x) - y[k,])^2))
  }
  z
}

distHellinger.simple_triplet_matrix <- function(x, y, ...) 
{
  if (!missing(y))
    stop("if x is a 'simple_triplet_matrix' y is not allowed to be specified.")
  x$v <- sqrt(x$v)
  z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  for (k in seq_len(nrow(x)-1)) {
    x1 <- x[-seq_len(k),]
    xk <- x[k,]
    xk$i <- rep(seq_len(nrow(x1)), each = length(xk$i))
    xk$j <- rep(xk$j, nrow(x1))
    xk$v <- rep(xk$v, nrow(x1))
    xk$nrow <- nrow(x1)
    z[k,-seq_len(k)] <- z[-seq_len(k),k] <- 
      sqrt(1/2 * row_sums((xk - x1)^2))
  }
  z
}

get_terms <- function(...) terms(...)
get_topics <- function(...) topics(...)

ldaformat2dtm <- function(documents, vocab, omit_empty = TRUE) {
  stm <- simple_triplet_matrix(i = rep(seq_along(documents), sapply(documents, ncol)),
                               j = as.integer(unlist(lapply(documents, "[", 1, )) + 1L),
                               v = as.integer(unlist(lapply(documents, "[", 2, ))),
                               nrow = length(documents),
                               ncol = length(vocab),
                               dimnames = list(names(documents), vocab))
  dtm <- as.DocumentTermMatrix(stm, weightTf)
  if (omit_empty)
    dtm <- dtm[row_sums(dtm) > 0,]
  dtm
}

dtm2ldaformat <- function(x, omit_empty = TRUE) {
  split.matrix <- 
    function (x, f, drop = FALSE, ...) 
      lapply(split(seq_len(ncol(x)), f, drop = drop, ...),
             function(ind) x[,ind, drop = FALSE])

  documents <- vector(mode = "list", length = nrow(x))
  names(documents) <- rownames(x)
  documents[row_sums(x) > 0] <- split(rbind(as.integer(x$j) - 1L, as.integer(x$v)), as.integer(x$i))
  if (omit_empty)
    documents[row_sums(x) == 0] <- NULL
  else 
    documents[row_sums(x) == 0] <- rep(list(matrix(integer(), ncol = 0, nrow = 2)), sum(row_sums(x) == 0))
  list(documents = documents,
       vocab = colnames(x))
}

setOldClass("DocumentTermMatrix")
setOldClass("simple_triplet_matrix")

setGeneric("perplexity", function(object, newdata, ...) standardGeneric("perplexity"))

setMethod("perplexity", signature(object = "VEM", newdata = "missing"), function(object, newdata, ...)  
  exp(-as.numeric(logLik(object))/object@n))

setMethod("perplexity", signature(object = "ANY", newdata = "matrix"), function(object, newdata, ...) 
          perplexity(object, as.simple_triplet_matrix(newdata), ...))

setMethod("perplexity", signature(object = "ANY", newdata = "DocumentTermMatrix"), function(object, newdata, ...)  {
  class(newdata) <- "simple_triplet_matrix"
  perplexity(object, newdata, ...)
})  

setMethod("perplexity", signature(object = "VEM", newdata = "simple_triplet_matrix"), function(object, newdata, control, ...) {
  CLASS <- strsplit(class(object), "_")[[1]]
  if (missing(control)) {
    control <- object@control
  } else {
    control <- as(control, paste(class(object), "control", sep = ""))
  }
  control@estimate.beta <- FALSE
  control@nstart <- 1L
  object_inf <- get(CLASS[1])(newdata, method = CLASS[2], model = object, control = control)
  perplexity(object_inf, ...)
})

setMethod("perplexity", signature(object = "Gibbs", newdata = "simple_triplet_matrix"), function(object, newdata, control, use_theta = TRUE, estimate_theta = TRUE, ...) {
  if (use_theta) {
    if (estimate_theta) {
      CLASS <- strsplit(class(object), "_")[[1]]
      if (missing(control)) {
        control <- object@control
      } else {
        control <- as(control, paste(class(object), "control", sep = ""))
      }
      control@estimate.beta <- FALSE
      control@nstart <- 1L
      object <- get(CLASS[1])(newdata, method = CLASS[2], model = object, control = control)
    }
    return(exp(-sum(log(colSums(exp(object@beta[,newdata$j]) * t(object@gamma)[,newdata$i])) * newdata$v)/sum(newdata$v)))
  } else {
    return(exp(-sum(log(colSums(exp(object@beta[,newdata$j])/object@k)) * newdata$v)/sum(newdata$v)))
  }
})

setMethod("perplexity", signature(object = "Gibbs_list", newdata = "simple_triplet_matrix"), function(object, newdata, control, use_theta = TRUE, estimate_theta = TRUE, ...) {
  if (use_theta) {
    if (estimate_theta) {
      CLASS <- strsplit(class(object@fitted[[1]]), "_")[[1]]
      if (missing(control)) {
        control <- lapply(object@fitted, slot, "control")
      } else {
        control <- rep(list(as(control, paste(class(object@fitted[[1]]), "control", sep = ""))), length(object@fitted))
      }
      object@fitted <- lapply(seq_along(object@fitted), function(i) {
        control[[i]]@estimate.beta <- FALSE
        control[[i]]@nstart <- 1L
        control[[i]]@iter <- control[[i]]@thin
        control[[i]]@best <- TRUE
        get(CLASS[1])(newdata, method = CLASS[2], 
                      model = object@fitted[[i]], control = control[[i]])
      })
    } else if (nrow(newdata) != nrow(object@fitted[[1]]@gamma)) stop("newdata needs to have the same number of documents")
    logs <- sapply(object@fitted, function(z)
                   log(colSums(exp(z@beta[,newdata$j]) * t(z@gamma)[,newdata$i])))
  } else {
    logs <- sapply(object@fitted, function(z)
                   log(colSums(exp(z@beta[,newdata$j])/z@k)))
  }
  exp(-sum(log(rowMeans(exp(logs))) * newdata$v)/sum(newdata$v))
})

setMethod("perplexity", signature(object = "list", newdata = "missing"), function(object, newdata, ...) {
  if (any(!sapply(object, inherits, "VEM"))) stop("if newdata is missing only VEM objects can be used")
  exp(-mean(sapply(object, logLik))/object[[1]]@n)
}) 

setMethod("perplexity", signature(object = "list", newdata = "simple_triplet_matrix"), function(object, newdata, ...) {
  perplexities <- sapply(object, perplexity, newdata = newdata, ...)
  exp(mean(log(perplexities)))
})
