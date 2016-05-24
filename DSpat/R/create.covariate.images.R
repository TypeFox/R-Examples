create.covariate.images=function(covariates, varnames)
#######################################################################
# Creates a list of covariates images from a dataframe of covariates
# defined on a grid for the study area.  Images are created for
# variables contained in varnames and the values of the covariates are
# in covariates
#
#  Arguments:
#
#   covariates  - covariate dataframe (see DSpat for structure)
#   varnames    - names of variables contains in covariates
#
#  Value: list of covariate images for covariates in formula derived
#  from covariates.  Elements in list are named using covariate names.
#
#  Jeff Laake
#  21 April 2008
#######################################################################
{
# Using the variables specified in formula construct a list of covariate images
  covariate.im=vector("list",length=length(varnames))
  names(covariate.im)=varnames
  y.pg = unique(covariates$y)
  x.pg = unique(covariates$x)
  for (cname in varnames)
    covariate.im[[cname]]=im(covariates[,cname],xcol=x.pg,yrow=y.pg)
  return(covariate.im)
}
