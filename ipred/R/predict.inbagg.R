predict.inbagg <- function(object, newdata, ...) {
  if(!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if(any(names(object$W) %in% names(newdata))) newdata <- newdata[!(names(newdata) %in% names(object$W))]
  NBAGG <- length(object$mtrees)
  N <- nrow(newdata)
  classes <- levels(object$y)
  vote <- matrix(0, nrow=N, ncol=length(classes))
  for(i in 1:NBAGG) {
    intermed <- object$mtrees[[i]]$bfct(newdata)
#    XX <- data.frame(newdata, intermed)
    if(!is.null(object$mtrees[[i]]$btree$fixed.function)) {
      names(intermed) <- sub(".[0-9]$", "", names(intermed))
      XX <- data.frame(newdata, intermed)
#      names(XX)[(ncol(XX)-ncol(intermed)+1):ncol(XX)] <- sub(".[0-9]$", "", names(intermed))
      res <- object$mtrees[[i]]$btree$fixed.function(XX)
    } else {
      XX <- data.frame(newdata, intermed)
      if(is.null(object$mtrees[[i]]$btree$predict)) {
        res <- try(predict(object$mtrees[[i]]$btree$model, newdata = XX, ...))
      } else {
        res <- try(object$mtrees[[i]]$btree$predict(object$mtrees[[i]]$btree$model, newdata = XX, ...))
      }
    }
    res <- cbind(1:N, res)
    vote[res] <- vote[res] +1 
  }

  RET <- factor(classes[apply(vote, 1, uwhich.max)])
  RET
}
