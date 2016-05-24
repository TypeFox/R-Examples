"powercalc" <- function(cross,n,effect,sigma2,env.var,gen.var,
                        thresh=3,sel.frac=1,theta=0,bio.reps=1)
  {
    # if error variance missing, calculate it
    if(missing(sigma2))
      {
        if((missing(env.var))|(missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
        sigma2 <- error.var(cross,env.var,gen.var,bio.reps)
      }
    else
      {
        if((!missing(env.var))|(!missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
      }

    
    if(cross=="bc")
      ans <- powercalc.bc(n,effect,sigma2,thresh,sel.frac,theta)
    else if(cross=="f2")
      ans <- powercalc.f2(n,effect,sigma2,thresh,sel.frac,theta)      
    else if(cross=="ri")
      ans <- powercalc.ri(n,effect,sigma2,thresh)      
    else
      stop("Unknown cross ", cross, ".")

    ans <- c(ans, 100*prop.var(cross,effect,sigma2))
    names(ans) <- c("power","percent.var.explained")
    t(ans)
  }

"detectable" <- function(cross,n,effect=NULL,
                         sigma2,env.var,gen.var,power=0.8,thresh=3,
                         sel.frac=1,theta=0,bio.reps=1) 
  {
    # argument check
    if(missing(sigma2))
      {
        if((missing(env.var))|(missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
        sigma2 <- error.var(cross,env.var,gen.var,bio.reps)
      }
    else
      {
        if((!missing(env.var))|(!missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
      }

    if(cross=="bc")
      ans <- detectable.bc(n,sigma2,power,thresh,sel.frac,theta)
    else if(cross=="f2")
      {
        if(is.null(effect))
          {
            ans <- detectable.f2(n,effect="add",sigma2,power,thresh,
                                 sel.frac,theta)
          }
        else
          {
            ans <- detectable.f2(n,effect=effect,sigma2,power,thresh,
                                 sel.frac,theta)
          }
      }
    else if(cross=="ri")
      ans <- detectable.ri(n,sigma2,power,thresh)
    else
      stop("Unknown cross ", cross, ".")

    ans <- c(ans, 100*prop.var(cross,ans,sigma2))

    if( cross=="f2" )
          names(ans) <- c("additive.effect","dominance.effect",
                          "percent.var.explained")
    else
      names(ans) <- c("effect","percent.var.explained")

    t(ans) 
  }

"samplesize" <- function(cross,effect,sigma2,env.var,gen.var,power=0.8,
                         thresh=3,sel.frac=1,theta=0,bio.reps=1) 
  {
    # argument check
    if(missing(sigma2))
      {
        if((missing(env.var))|(missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
        sigma2 <- error.var(cross,env.var,gen.var,bio.reps)
      }
    else
      {
        if((!missing(env.var))|(!missing(gen.var)))
          stop("Need either sigma2 or both env.var and gen.var.")
      }

    if(cross=="bc")
      ans <- samplesize.bc(effect,sigma2,power,thresh,sel.frac,theta)
    else if(cross=="f2")
      ans <- samplesize.f2(effect,sigma2,power,thresh,sel.frac,theta)
    else if(cross=="ri")
      ans <- samplesize.ri(effect,sigma2,power,thresh)
    else
      stop("Unknown cross ", cross, ".")

    ans <- c(ans, 100*prop.var(cross,effect,sigma2))
    names(ans) <- c("sample.size","percent.var.explained")

    t(ans)
  }



"powercalc.bc" <- function(n,effect,sigma2,thresh=3,sel.frac=1,theta=0)
  {
    delta <- (effect/2)/sqrt(sigma2)
    # effective sample size
    m <- n * info.bc(sel.frac,theta)
    # non-centrality parameter
    ncp <- m*delta^2
    if( m<30 )
      {
        if( (sel.frac<1) | (theta>0) )
          {
            warning("Approximation not reliable as effective sample size < 30.")
          }
      }
    # threshold in 2*loglikelihood units
    T <- 2*log(10)*thresh
    # power using non-central chi-square
    pow <- 1-pchisq(T,df=1,ncp=ncp)
    return(pow)
  }


"detectable.bc" <- function (n, sigma2,
                             power = 0.8, thresh = 3, sel.frac = 1, theta = 0) 
{
  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
    effect <- uniroot(function(x) {
      powercalc.bc(n, x, sigma2, thresh, sel.frac, theta) -
        power}, interval = c(0,30*sqrt(sigma2/n)))$root
    effect
}

"samplesize.bc" <- function (effect, sigma2, power = 0.8, thresh = 3,
                             sel.frac = 1, theta = 0)
{

  # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.bc(2^m, effect, sigma2, thresh,
                                         sel.frac, theta))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.bc(n, effect, sigma2, thresh, sel.frac, theta) - power
    }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))

}


"powercalc.f2" <-  function (n, effect, sigma2, thresh = 3, sel.frac = 1,
                             theta = 0) 
{
  if(length(effect)!=2)
    stop("Incorrect effect specification; need vector of size 2.")
  
  # get info per individual
  iii <- info.f2(sel.frac, theta)
  # additive and dominance components
  a <- effect[1]/sqrt(sigma2)
  d <- effect[2]/sqrt(sigma2)
  # calculate non-centrality parameter
  ncp <- n * ( iii$add*a^2/2 + iii$dom*d^2/4 )
  m <- n*min(iii$add,iii$dom)
  # if effective sample size not big enough, stop
  if (m < 30)
    {
      if( (sel.frac<1) | (theta>0) )
        {
          warning("Approximation not reliable as effective sample size < 30.")
        }
    }
  # calculate threshold in chi-square scale
  T <- 2 * log(10) * thresh
  # calculate power
  pow <- 1 - pchisq(T, df = 2, ncp = ncp)
  pow

}


"detectable.f2" <- function (n, effect="add", sigma2, power = 0.8, thresh = 3,
                             sel.frac = 1, theta = 0)
{
  # effect
  if(effect=="add")
    {
      a <- 1
      d <- 0
      effect <- c(a,d)
    }
  else if(effect=="dom")
    {
      a <- 1
      d <- 1
      effect <- c(a,d)
    }
  else if( is.numeric(effect) && (length(effect)==2) )
    {
      effect <- effect
    }
  else
    {
      stop("Cannot understand effect argument.")
    }

  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
  del <- uniroot(function(x) {
    powercalc.f2(n, x*effect, sigma2, thresh, sel.frac, theta) - power  },
                 interval = c(0,30*sqrt(sigma2/n)))$root

    # decide what to return depending on delta flag
  return(del*effect)
}

"samplesize.f2" <- function (effect, sigma2, power = 0.8, thresh = 3,
                             sel.frac = 1, theta = 0)
{

  if(length(effect)!=2)
    {
      warning("Assuming additive effect.")
      effect <- c(effect[1],0)
    }
  # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.f2(2^m, effect, sigma2, thresh,
                                         sel.frac, theta))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.f2(n, effect, sigma2, thresh, sel.frac, theta) - 
          power }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))

}



"powercalc.ri" <- function(n,effect,sigma2,thresh=3)
{
  powercalc.bc(n,effect*2,sigma2,thresh,sel.frac=1,theta=0)
}

"detectable.ri" <- function (n, sigma2, power = 0.8, thresh = 3) 
{
  # proportion of variance explained for given sample size,
  # power, threshold, selection fraction, and size of marker interval
    effect <- uniroot(function(x) {
      powercalc.ri(n, x, sigma2, thresh) -  power},
                     interval = c(0,30*sqrt(sigma2/n)))$root
    effect
}


"samplesize.ri" <- function (effect,sigma2,power = 0.8, thresh = 3)
{

   # find an interval for the sample size to search
  # search in powers of 2
  p <- 0
  m <- 0
  while(p<power)
    {
      m <- m+1
      p <- suppressWarnings(powercalc.ri(2^m, effect, sigma2, thresh))
    }

  # refine the number solving the power equation
    n <- uniroot(function(n) {
        powercalc.ri(n, effect, sigma2, thresh) - power
    }, interval = c(2^(m-1), 2^m))$root

  # return the nearest largest integer
    return(ceiling(n))
}



"error.var" <- function(cross,env.var=1,gen.var=0,bio.reps=1)
  {
    # get the genetic variance multiplier 
    if(cross=="bc")
      CC <- 1/4
    else if (cross=="f2")
      CC <- 1/2
    else if (cross=="ri")
      CC <- 1
    else
      stop("Cross type not recognized.")

    env.var/bio.reps + CC*gen.var
    
  }



# function to convert the genotype means to additive and dominant
# effects segegating in a cross
"gmeans2effect" <- function(cross,means)
{
aa <- means[1]
ab <- means[2]
bb <- means[3]

if( cross == "f2" )
  {
    a <- (aa-bb)/2
    d <- ab - (bb+aa)/2
    effect <- c(a,d)
  }
else if( cross == "bc" )
  {
    effect <- c(aa-ab,ab-bb)
  }
else if (cross == "ri" )
  {
    effect <- (aa-bb)/2
  }

effect
}

"prop.var" <- function(cross,effect,sigma2)
  {
    if(cross=="bc")
      ans <- effect^2/4
    else if (cross=="f2")
      ans <- effect[1]^2/2 + effect[2]^2/4
    else if (cross=="ri")
      ans <- effect^2
    else
      stop("Cross type not recognized.")

    ans <- ans/(ans+sigma2)
    ans
  }
