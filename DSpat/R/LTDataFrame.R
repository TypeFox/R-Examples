LTDataFrame=function(study.area, lines, lines.psp, int.formula, det.formula,
                            covariates, Q.lt)
################################################################################
# Creates a dataframe of covariates for data and dummy points of the quadrature
#
# Arguments:
#
#  study.area   - owin object that defines study area
#  lines        - data frame of lines with structure as shown in quadscheme.lt
#  lines.psp    - psp class with added list elements width and label
#  int.formula  - model formula for intensity process
#  det.formula  - model formula for detection process
#  covariates   - covariate dataframe (see DSpat for structure)
#  Q.lt         - lt quadscheme of class quad
#
# Value:
#  dataframe with covariates at observations followed  by (rows) of covariates
#        at dummy points
#  list of covariate images
#
# Devin Johnson & Jeff Laake
# 4 April 2008
################################################################################
{
# Make sure all variable names are included in covariates
# or in lines. If not, stop
  varnames=all.vars(int.formula)
  if("distance" %in% varnames)
    stop("\n distance should not be included in intensity formula\n")
  varnames=varnames[!varnames %in% c("x","y")]
  covnames=names(covariates)
  if(any(names(lines)%in% varnames[!varnames %in%covnames]))
    stop(paste("\n Covariates in lines cannot be used in intensity formula\n"))
  if(any(!varnames %in% covnames))
    stop(paste("\n Following variables specified in intensity formula not in covariates\n",
        paste(varnames[!varnames %in% covnames],collapse=" ")))
#
# Construct the detection formula as an interaction with -distance^2/2
#
  detvars=all.vars(as.formula(paste("~",det.formula,sep="")))
  if("distance" %in% detvars)
    stop("\n distance should not be included in detection formula because it specifies interaction with distance\n")
  if(det.formula!= ~-1)
  {
     if(det.formula==~1)
       det.formula="I(-distance^2/2)"
     else
        det.formula=paste("(",as.character(det.formula)[2],"):I(-distance^2/2)",sep="")
     detvars=all.vars(as.formula(paste("~",det.formula,sep="")))
     detvars=detvars[!detvars %in% c("distance","x","y")]
     varnames=c(varnames,detvars[detvars %in% covnames])
     xvars=detvars[!detvars %in% covnames]
     if(any(!xvars %in% names(lines)))
       stop(paste("\n Following variables specified in formula not in covariates or in lines\n",
          paste(varnames[!detvars %in% names(lines)],sep=" ")))
     formula=formula(paste("~",as.character(int.formula)[2],"+",det.formula,collapse=""))
     varnames=unique(c(varnames,xvars))
  }
  else
     formula=int.formula
# Turn spatial covariates into spatstat images
  covariate.im <- create.covariate.images(covariates,varnames)
# Extract elements
  obs.ppp = Q.lt$data
  dummy.ppp = Q.lt$dummy
  obs.labels=obs.ppp$label
  dummy.labels=dummy.ppp$label
  ends = lines.psp$ends
  labels=lines.psp$label
  number.lines=length(labels)
# Get covariates for observations contained in cov.im
  covnms = names(covariate.im)
  cov.obs=function(label, obs.ppp, lines.psp, covariates)
  {
    pts = obs.ppp[obs.ppp$label==label]
    ends = lines.psp$ends[lines.psp$label==label,]
    distance = dist2line(pts,ends)$distance
    if(pts$n>0)
    {
       covdata =  lapply(covnms,
                  function(nm, covariates){assign(nm, covariates[[nm]][pts])},
                                 covariates=covariates)
       names(covdata) = covnms
       covdata = as.data.frame(covdata)
       if(dim(covdata)[1]==0)
          return(data.frame(cbind(distance,label=label)))
       else
          return(data.frame(cbind(distance,covdata,label=label)))
    }
    else
      return(NULL)
  }
  obs.df =  do.call('rbind', lapply(labels,cov.obs,
                   obs.ppp=obs.ppp,lines.psp=lines.psp,covariates=covariate.im))
# Merge with covariates from lines used in formula; these should only be
# covariates that affect detection probability because they are no known for the
# entire area of innference
  if(any(names(lines)%in% detvars))
     obs.df=merge(obs.df,subset(lines,
               select=c("label",names(lines)[names(lines)%in% detvars])))
# Get covariates for dummy points contained in cov.im
  cov.dummy=function(label, dummy.ppp, lines.psp, covariates)
  {
       pts = dummy.ppp[dummy.ppp$label==label]
       pts = pts[!is.na(pts$x)&!is.na(pts$y),]
       ends = lines.psp$ends[lines.psp$label==label,]
       distance = dist2line(pts,ends)$distance
       covdata =  lapply(covnms,
                          function(nm, covariates){assign(nm, covariates[[nm]][pts])},
                                covariates=covariates)
       names(covdata) = covnms
       covdata = as.data.frame(covdata)
       if(dim(covdata)[1]==0)
          return(data.frame(cbind(distance,label=label)))
       else
          return(data.frame(cbind(distance,covdata,label=label)))
  }
  dummy.df =  do.call('rbind',lapply(labels,cov.dummy,
                dummy.ppp=dummy.ppp,lines.psp=lines.psp,covariates=covariate.im))
  if(any(names(lines)%in% varnames))
    dummy.df=merge(dummy.df,subset(lines,
             select=c("label",names(lines)[names(lines)%in% varnames])))
  return(list(formula=formula,cov.df=rbind(obs.df, dummy.df),
                 covariate.im=covariate.im))
}


