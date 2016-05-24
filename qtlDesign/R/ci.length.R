ci.length <- function(cross,n,effect,p=0.95,sigma2=1,
                      env.var,gen.var,bio.reps=1)
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
    {
      drift <- effect^2/sigma2
    }
  else if(cross=="f2")
    {
      drift <- (2*effect[1]^2+effect[2]^2)/sigma2
    }
  else if(cross=="ri-self")
    {
      drift <- 2*effect^2/sigma2
    }
  else if(cross=="ri-sib")
    {
      drift <- 4*effect^2/sigma2
    }
  else
    stop("Unknown cross ", cross, ".")

  ans <- 400*qexp(p)/(n*drift)
  ans
}
