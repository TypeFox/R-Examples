rdigitsBFOS <- function(n, eta=0.25, alpha=NULL, silent=FALSE) {
  stopifnot(eta>0.01 && eta <= 0.5)
  z <- rep(c(1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,
             0,1,0,1,1,0,1,0,1,1,1,1,0,1,1,1,1,1,0,1,0,0,1,0,1,1,1,1,1,1,1,1,1,
             1,1,0,1,1),n)
    if (is.null(alpha)) {#default: generate with predetermined eta
    alpha <- (-0.00014393510482067046 + eta*(0.3743903928013429 +
              eta*(0.003068677095495819 +
              eta*0.003068677095495819)))
    } else { #generate with input alpha, ignore input eta
      stopifnot(alpha>0 && alpha <= 0.2)
      eta <- (-0.00018468194429807232 + alpha*(2.7140277221213913 -
                alpha*(0.876012944101314 + alpha*2.51929256972748)))
    }
    m <- matrix((z + ifelse(runif(length(z)) > alpha, 0, 1))%%2,
                byrow=TRUE, ncol=7)
    colnames(m) <- paste0("x", 1:7)
    ans<-cbind(as.data.frame.matrix(m), digit=rep(factor(0:9),n))
    title <- paste("Electronic Digit Recognition Problem, alpha =",
                   round(alpha, 5),
                   "BayesRate =", round(eta, 5))
    if(!silent) {
      cat(title, fill=TRUE)
    }
    attr(ans, "title") <- title
    ans
}
