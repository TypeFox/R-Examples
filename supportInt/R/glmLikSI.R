#' @import ProfileLikelihood
#' @import splines
#' @import stats
#' @export
glmLikSI <- function(glm.mod, level){
  if(class(glm.mod)[1]!="glm") {stop("First argument must be of class glm")}
  na <- colnames(glm.mod$model)
  cis <- confint(glm.mod, level=.9999)
  sis <- matrix(nrow=length(glm.mod$coef)-1, ncol=2)
  rownames(sis) <- names(glm.mod$coef)[-1]
  for(i in 2:length(na)){
    if(length(na)==2){form <- as.formula(paste(na[1], "~", 1))
    } else {
      form <- as.formula(paste(na[1], "~", paste(na[-c(1,i)], collapse="+")))
    }
    pl <- profilelike.glm(form, glm.mod$model, profile.theta=na[i], 
                          family=glm.mod$family, lo.theta=cis[i,1], hi.theta=cis[i,2])
    nk <- 8
    fl <- 0
    while(fl==0){
      sp <- ns(pl$theta,nk)
      mod <- lm(pl$profile.lik.norm~sp)
      if(1-summary(mod)$r.squared<.01) {
        fl <- fl+1
      } else {nk <- nk+2}
    }
    kn <- as.numeric(sort(c(attr(sp, "knots"), attr(sp, "Boundary.knots"))))
    p.lik <- function(z) sapply(z, function(a) c(1, predict(sp,a))%*%mod$coef-1/level)
    sis[i-1,] <- c(uniroot(p.lik, c(cis[i,1], glm.mod$coef[i]))$root, uniroot(p.lik, c(glm.mod$coef[i], cis[i,2]))$root)
  }
  colnames(sis) <- c(paste("low 1/", level, sep=""), paste("upp 1/", level, sep=""))
  return(sis)
}