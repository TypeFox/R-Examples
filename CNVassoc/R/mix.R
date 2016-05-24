mix <-
function (mixdat, method, num.class, mean.ini, sigma.ini, pi.ini, var.equal)
{
    if (length(num.class) > 1) {
        G.opt <- Mclust(mixdat, G = num.class)$G
    }
    else {
        G.opt <- num.class
    }
    if (method == 1) {
        if (missing(mean.ini))
          mean.ini <- seq(min(mixdat), max(mixdat), len = G.opt)
        if (missing(sigma.ini))
          sigma.ini <- rep(0.1, G.opt)
        if (missing(pi.ini))
          pi.ini <- rep(1 / G.opt, G.opt)
        ans <- mixdist::mix(mixgroup(mixdat), mixparam(mean.ini, sigma.ini, pi.ini), emsteps = 1)
        ans$method <- "function mix {package: mixdist}"
        ans$G <- G.opt
    }
    if (method == 2) {
        if (missing(var.equal))
          ans <- Mclust(mixdat, G = G.opt)
        else 
          if (var.equal)
            ans <- Mclust(mixdat, G = G.opt, modelNames = "E")
          else 
            ans <- Mclust(mixdat, G = G.opt, modelNames = "V")
        ans$parameter$mu <- ans$parameters$mean
        ans$parameter$sigma <- sqrt(ans$parameters$variance$sigmasq)
        if (ans$modelName == "E")
          ans$parameter$sigma <- rep(ans$parameter$sigma, ans$G)
        ans$parameter$pi <- ans$parameter$pro
        ans$P <- NULL
        ans$method <- "function Mclust {package: mclust}"
    }
    if (method == 3) {
        if (missing(mean.ini))
          mean.ini <- seq(min(mixdat), max(mixdat), len = G.opt)
        if (missing(sigma.ini))
          sigma.ini<- rep(0.1, G.opt)
        if (missing(pi.ini))
          pi.ini <- rep(1 / G.opt, G.opt)          
        ans.var.equal <- EMmixt(mixdat, mean.ini, sigma.ini, pi.ini, var.equal = TRUE)
        ans.var.inequal <- EMmixt(mixdat, mean.ini, sigma.ini, pi.ini, var.equal = FALSE)
        if (missing(var.equal))
          if (ans.var.inequal$bic > ans.var.equal$bic)
            ans <- ans.var.equal
          else 
            ans <- ans.var.inequal
        else
          if (var.equal)
            ans <- ans.var.equal
          else 
            ans <- ans.var.inequal
        ans$G <- G.opt
        ans$parameter <- list()
        ans$parameter$mu <- ans$mu.surrog
        ans$parameter$sigma <- ans$sd.surrog
        ans$parameter$pi <- ans$pi
        ans$P <- NULL
        ans$method <- "function EMmixt"        
    }
    ans
}
