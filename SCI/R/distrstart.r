lmom.start <- function(x,distr=c("gamma","genlog","gev","gumbel",
                             "lnorm","norm","pe3","weibull"),...){
### estimates parameters of distributions using L-moments
###
### x : data vector
### distr : A character string "name" naming a distribution for which
### the corresponding density function ('dname'), the
### corresponding distribution function ('pname') and the
### quantile function ('qname') must be defined (see for example 'GammaDist'
    distr <- match.arg(distr)
    x.lmom <- try(lmom.ub(x),TRUE)
    if(distr == "gamma"){
        ppar <- try(pargam(x.lmom,checklmom=FALSE),TRUE)
        ppar <-  if(class(ppar)!="try-error"){
            ppar$para
        } else {
            c(NA,NA)
        }
        names(ppar) <- c("shape","rate")        
    } else if(distr == "genlog"){
        ppar <- try({
            k.est <- -x.lmom$TAU3
            a.est <- x.lmom$L2/(gamma(1+k.est)*gamma(1-k.est))
            e.est <- x.lmom$L1+(x.lmom$L2-a.est)/k.est
            c(shape=k.est,scale=a.est,location=e.est)
        },TRUE)
          if(class(ppar)=="try-error")
              ppar <- c(shape = NA, scale = NA, location=NA)           
    } else if(distr == "gev"){
        if(requireNamespace("evd", quietly = TRUE)){
            ppar <- try(pargev(x.lmom,checklmom=FALSE),TRUE)
            ppar <-  if(class(ppar)!="try-error"){
                ppar$para <- c(loc=ppar$para["xi"],
                               scale=ppar$para["alpha"],
                               shape=-ppar$para["kappa"]
                               )            
            } else {
                c(NA,NA,NA)
            }
            names(ppar) <- c("loc","scale","shape")
        } else {
            stop("package 'evd' needed fro 'distr=gev'")
        }
    } else if(distr == "gumbel"){
        if(requireNamespace("evd", quietly = TRUE)){
            ppar <- try(pargum(x.lmom,checklmom=FALSE),TRUE)
            ppar <-  if(class(ppar)!="try-error"){
                ppar$para
            } else {
                c(NA,NA)
            }
            names(ppar) <- c("loc","scale")
        } else {
            stop("package 'evd' needed fro 'distr=gumbel'")
        }
    } else if(distr == "lnorm"){
        ppar <- try(parln3(x.lmom,checklmom=FALSE),TRUE)
        ppar <-  if(class(ppar)!="try-error"){
            ppar$para[c("mulog","sigmalog")]
        } else {
            c(NA,NA)
        }
        names(ppar) <- c("meanlog","sdlog")
    } else if(distr == "norm"){
        ppar <- try(parnor(x.lmom,checklmom=FALSE),TRUE)
        ppar <-  if(class(ppar)!="try-error"){
            ppar$para
        } else {
            c(NA,NA)
        }
        names(ppar) <- c("mean","sd") 
    } else if(distr == "pe3"){
        ppar <- try({
            estimate <-  parpe3(x.lmom,checklmom=FALSE)
            if(estimate$para[[3]] > 1.8)
		estimate$para[[3]] <- 1.75
            if(estimate$para[[3]] < -1.8)
		estimate$para[[3]]<- -1.75
            c(shape=estimate$para[[3]],scale=estimate$para[[2]],location=estimate$para[[1]])
        },TRUE)
        if(class(ppar)=="try-error")
            ppar <- c(shape=NA, scale=NA,location=NA)
    } else if(distr == "weibull"){
        ppar <- try(parwei(x.lmom,checklmom=FALSE),TRUE)
        ppar <-  if(class(ppar)!="try-error"){
            ppar$para[c("delta","beta")]
        } else {
            c(NA,NA)
        }
        names(ppar) <- c("shape","scale")
    }
    ppar <- as.list(ppar)
    return(ppar)
}


mom.start <- function(x,distr=c("gamma","gumbel","logis","lnorm","norm",
                            "weibull"),...){
### estimates parameters of distributions using moments
###
### x : data vector
### distr : A character string "name" naming a distribution for which
### the corresponding density function ('dname'), the
### corresponding distribution function ('pname') and the
### quantile function ('qname') must be defined (see for example 'GammaDist'
    distr <- match.arg(distr)
    if(distr %in% c("norm","lnorm","gamma","logis")){
        ppar <- try(fitdistrplus::mmedist(x,distr),TRUE)
        ppar <- if(class(ppar)!="try-error"){
            ppar$estimate
        } else {
            switch(distr,
                   norm=c(mean = NA, sd= NA ),
                   lnorm=c(meanlog = NA, sdlog = NA),
                   gamma=c(shape = NA, rate = NA),
                   logis=c(location = NA, scale = NA))
        }
    } else if(distr == "gumbel"){
        if(requireNamespace("evd", quietly = TRUE)){
            ppar <- try({
                x <- x[x>0]
                m <- mean(x)
                sd.data <- sd(x) 
                scale <- (sd.data*(6)^0.5)/pi
                location <- m-0.5772*scale
                c(loc = location,scale = scale)
            },TRUE)
            if(class(ppar)=="try-error")
                ppar <- c(location = NA,scale = NA)
        } else {
            stop("package 'evd' needed fro 'distr=gumbel'")
        }
    } else if(distr == "weibull"){
        ## adapted from 'mledist' in package 'fitdistrplus'
        ppar <- try({
            m <- mean(log(x))
            v <- var(log(x))
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)
            c(shape = shape, scale = scale)
        },TRUE)
        if(class(ppar)=="try-error")
            ppar <- c(shape = NA, scale = NA)
    }
    ppar <- as.list(ppar)
    return(ppar)
}



dist.start <- function(x,distr,...){
### estimates starting values for distributions
###
### first 'lmom.start' is called. In case of failure of 'lmom.start'
### an attempt is undertaken to use 'mom.start'
###
### x : data vector
### distr : A character string "name" naming a distribution for which
### the corresponding density function ('dname'), the
### corresponding distribution function ('pname') and the
### quantile function ('qname') must be defined (see for example 'GammaDist'
    par1 <- try(suppressWarnings(lmom.start(x=x,distr=distr,...)),silent=TRUE)
    if(class(par1)!="try-error"&all(is.finite(unlist(par1)))){
        return(par1)
    } else {
        par2 <- try(suppressWarnings(mom.start(x=x,distr=distr,...)),silent=TRUE)
        if(class(par2)=="try-error"){
            if(class(par1)!="try-error"){
                return(par1)
            } else {
                stop("parameters of distribution",distr,"could not be identifyed")
            }
        } else {
            return(par2)
        }
    }
}
 
