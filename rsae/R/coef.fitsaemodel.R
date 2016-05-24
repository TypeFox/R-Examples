coef.fitsaemodel <-
function(object, type="both", ...){
   model <- attr(object, "saemodel")
   theta <- object$theta
   beta <- object$beta
   if (type == "both"){
      raneff <- t(as.matrix(object$theta))
      colnames(raneff) <- c("ResidualVar", "AreaVar")
      rownames(raneff) <- "raneff"
      fixeff <- t(as.matrix(object$beta))
      colnames(fixeff) <- colnames(model$X)
      rownames(fixeff) <- "fixeff"
      res <- list(fixeff=fixeff, raneff=raneff)
      print(fixeff)
      cat("\n")
      print(raneff)
   }
   if (type == "raneff"){
      res <- t(as.matrix(object$theta))
      colnames(res) <- c("ResidualVar", "AreaVar")
      rownames(res) <- "raneff"
      print(res)
   }
   if (type == "fixeff"){
      res <- t(as.matrix(object$beta))
      colnames(res) <- colnames(model$X)
      rownames(res) <- "fixeff"
      print(res)
   }
   invisible(res)
}

