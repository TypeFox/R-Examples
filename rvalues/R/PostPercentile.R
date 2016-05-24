PostPercentile <- function(object)  {
   
    ### object - an object of class "rvalues"

    if(object$aux$prior=="conjugate")  {
        if(object$aux$family=="gaussian")  {
           XX <- object$aux$unsorted$MLE
           samp.var <- object$aux$unsorted$SE^2
           prior.mean <- object$aux$hypers[1]
           prior.var <- object$aux$hypers[2]
    
           denom <- sqrt((samp.var + prior.var)*(prior.var^2 + 2*samp.var*prior.var))
           num <- (prior.var + 2*samp.var)*prior.mean + prior.var*XX
           pep <- pnorm(num/denom,lower.tail=FALSE)
        }
        else if(object$aux$family=="poisson") {
           XX <- object$aux$unsorted$xx
           eta <- object$aux$unsorted$eta
           
           ### For poisson-gamma, the hyperparameters correspond
           ### to the rate parameter version of the gamma
           
           prior.shape <- object$aux$hypers[1]
           prior.rate <- object$aux$hypers[2]
           
           num <- prior.shape*(eta + prior.rate)
           den <- prior.rate*(prior.shape + XX)
           
           pep <- pf(num/den,df1=2*(prior.shape + XX),df2=2*prior.shape)
        }
        else if(object$aux$family=="binomial") {
           
           ff <- function(u,xx,nn,alpha,beta) {
              ans <- pbeta(u,shape1=alpha+xx,shape2=beta + nn - xx)*dbeta(u,shape1=alpha,shape2=beta)
              return(ans)
           }
           x <- object$aux$unsorted$xx
           m <- object$aux$unsorted$nn
           a <- object$aux$hypers[1]
           b <- object$aux$hypers[2]
           
           pep <- rep(0,length(x))
           for(i in seq_len(length(x))) {
              pep[i] <- integrate(ff,lower=0,upper=1,xx=x[i],nn=m[i],alpha=a,beta=b)$value
           }
           
           #pep <- 1 - (object$aux$V%*%delta.alpha)
        }
    }
    else if(object$aux$prior=="nonparametric")  {
         ## object$aux$prior is null, for example if using the Valpha() function
         
         ### get distances between points on the alpha grid.
         delta.alpha <- c(object$aux$alpha.grid[1],diff(object$aux$alpha.grid))

         ### Compute estimate posterior percentiles.
         ### This relies on the fact that 
         ###  pep = 1 - \int_0^1 V_\alpha(D_i) d\alpha
         
         pep <- 1 - (object$aux$V%*%delta.alpha)
    
    }
    return(pep)
}

