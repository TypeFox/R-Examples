transformedTrain <- function(obj,
                             target=obj$data$event,
                             select
                             ) {

  ## check that obj inherits from gam or glm
  if (!inherits(obj, "glm"))
    stop("use only with \"gam\" or \"glm\" objects")

  nbVar <- dim(obj$data)[1]
  if (inherits(target,"integer")) {
    if (!identical(sort(unique(target)),as.integer(c(0,1)))) {
      ## target is interpreted as a vector of indices
      newTarget <- integer(nbVar)
      newTarget[target] <- 1
      target <- newTarget
      rm(newTarget)
    }
  } else {
    if (identical(sort(unique(target)),c(0,1))) {
      ## target is interpreted as a vector of indices
      newTarget <- integer(nbVar)
      newTarget[target] <- 1
      target <- newTarget
      rm(newTarget)
    }
  } ## End of conditional on inherits(target,"integer")
  ## target is now a binary vector
  
  if (!missing(select)) {
    if (class(select)[1] == "character") {
      myCommand <- paste("with(obj$data,",select,")")
      select <- eval(parse(text=myCommand))
    } ## End of conditional on class(select)[1] == "character"
  } else {
    select <- !logical(dim(obj$data)[1])
  } ## End of consitional on !missing(subset) 

  if (!is.null(obj$na.action)) {
    subset <- - as.integer(obj$na.action)
  } else {
    subset <- seq(dim(obj$data)[1])
  }

  target <- target[subset][select[subset]]
  Lambda <- cumsum(obj$fitted.values[select[subset]])
  tST <- as.spikeTrain(Lambda[target==1])
  class(tST) <- c("transformedTrain",class(tST))
  tST
  
}
