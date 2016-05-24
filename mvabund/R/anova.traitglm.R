anova.traitglm <- function(object, ..., nBoot=99, resamp="pit.trap", test="LR", block = NULL, show.time="all", bootID=NULL)
{
    if("manyglm" %in% class(object) == FALSE)
      stop("Sorry, the anova function only works for traitglm objects fitted using manyglm.")
  
      n.sites    = dim(object$L)[1]
    n.spp      = dim(object$L)[2]
    if(is.null(bootID))
    {
      if(is.null(block)) #use block argument if provided, else make one up
        blockID = 1:n.sites
      else
      {
        if(is.factor(block)==FALSE)
          stop("block argument must be a factor")
        blockID = block
      }
      block      = factor(rep(blockID,n.spp))
    }
    else
      block=NULL

    dots <- list(...)
    ndots <- length(dots)
    if (ndots==0)
    {
      env.times.trait = object
      object$call$get.fourth = FALSE
      env.plus.trait = eval(object$call)

      #exception handling - don't call anova if there are no fourth corner terms in the model
      if(length(env.times.trait$fourth.corner)==0)
         stop("Sorry, your trait model has no fourth corner terms in it so you can't call anova without specifying an alternate model too")
      an = anova.manyglm(env.plus.trait, env.times.trait, nBoot=nBoot, resamp=resamp, test=test, block=block, show.time=show.time, bootID=bootID)
    }
    else
    {
      example = object
      an = anova.manyglm(example, ..., nBoot=nBoot, resamp=resamp, test=test, block=block, show.time=show.time, bootID=bootID)
    }
    
    #get rid of multivariate terms that don't apply, since it is a univariate fit:
    an$p.uni="none"
    an$cor.type="I"
    an$uni.p[is.numeric(an$uni.p)] <- NA
    if(ndots==0)
      row.names(an$table)=c("Main effects only", "env:trait (fourth corner)")
    return(an)
}
