predictPTLparams <- function(d,PTLmodel){
        beta <- max(c(0.01,predict.lm(PTLmodel$beta,newdata=data.frame(dists=d))))
        alpha <- max(c(0,predict.lm(PTLmodel$alpha,newdata=data.frame(dists=d))))
        gamma <- min(c(0,predict.lm(PTLmodel$gamma,newdata=data.frame(dists=d))))
        list(alpha=alpha,beta=beta,gamma=gamma)
}

