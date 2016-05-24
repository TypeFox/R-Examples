#' @rdname rfPermute
#' 
#' @importFrom stats na.fail model.response model.frame terms reformulate
#' @export rfPermute.formula
#' @export
#' 
rfPermute.formula <- function(formula, data = NULL, ..., subset, 
                              na.action = na.fail, nrep = 100) {
  if (!inherits(formula, "formula")) stop("method is only for formula objects")
  m <- match.call(expand.dots = FALSE)
  if (any(c("xtest", "ytest") %in% names(m))) {
    stop("xtest/ytest not supported through the formula interface")
  }
  
  # extract formula terms
  names(m)[2] <- "formula"
  if (is.matrix(eval(m$data, parent.frame()))) m$data <- as.data.frame(data)
  m$... <- NULL
  m$nrep <- NULL
  m$na.action <- na.action
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  y <- model.response(m)
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  m <- model.frame(terms(reformulate(attributes(Terms)$term.labels)), data.frame(m))
  for (i in seq(along = ncol(m))) {
    if (is.ordered(m[[i]])) m[[i]] <- as.numeric(m[[i]])
  }
  
  # run rfPermute
  rf.call <- match.call()
  rf.call[[1]] <- as.name("rfPermute")
  names(rf.call)[2:3] <- c("x", "y")
  rf.call$x <- m
  rf.call$y <- y
  rf.call$subset <- rf.call$na.action <- NULL
  rf.call[-1] <- lapply(rf.call[-1], eval, envir = parent.frame())
  rf <- eval(rf.call)
  
  # reconstitute original randomForest call
  rf.call <- match.call()
  rf.call[[1]] <- as.name("randomForest")
  rf$call <- rf.call
  rf$call$nrep <- NULL
  rf$terms <- Terms
  if (!is.null(attr(m, "na.action"))) rf$na.action <- attr(m, "na.action")
  
  class(rf) <- if(nrep > 0) {
    c("rfPermute", "randomForest.formula", "randomForest")
  } else {
    c("randomForest.formula", "randomForest")
  }
  
  return(rf)
}
