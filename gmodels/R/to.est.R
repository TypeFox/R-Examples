# to.est.R
# return a vector for cm in estimable()
# Randy Johnson
# Laboratory of Genomic Diversity at NCI-Frederick

.to.est <- function(obj, params)
{
  ## if('lme' %in% class(obj) | 'mer' %in% class(obj))
  ##   {
  ##     eff.obj <- fixef(obj)
  ##   }
  ## else
  if('geese' %in% class(obj))
    {
      eff.obj <- obj$beta
    }
  else
    {
      eff.obj <- coef(obj)
    }

  if(is.null(obj))
    stop("Error obtaining model coefficients")

  est <- rep(0, length(eff.obj))
  names(est) <- names(eff.obj)

  if(!missing(params))
    {
      if(is.null(names(params)))
        if(length(params)==length(est))
          names(params) <- names(est)
        else
          stop("'param' has no names and does not match number of coefficients of model. Unable to construct coefficient vector")
      else
        {
          matches <- names(params) %in% names(est)
          if(!(all(matches)))
             stop(
                  '\n\t',
                  'Invalid parameter name(s): ',
                  paste(names(params)[!matches], collapse=', '),
                  '\n\t',
                  'Valid names are: ',
                  paste(names(est), collapse=', ')
                  )
        }

      if(is.list(params))
        est[names(params)] <- unlist(params)
      else
        est[names(params)] <- params
    }

  return(est)
}
