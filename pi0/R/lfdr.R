lfdr=ppee=function(object, ...)UseMethod('lfdr')
lfdr.default=#ppee.default=
function(object, ...)
{
    idx=c(grep('^lfdr$', names(object), ignore.case=TRUE), grep('^ppee$', names(object), ignore.case=TRUE))
    if(length(idx)>0) return (object[[idx[idx]]])
    idx=c(grep('^lfdr$', names(attributes(object)), ignore.case=TRUE), grep('^ppee$', names(attributes(object)), ignore.case=TRUE))
    if(length(idx)>0) return (attr(object, names(attributes(object))[idx[idx]]))
    NA
}

lfdr.parncpF=
lfdr.nparncpF=
function(object, ...)
{
	if(any(is.na(object))) return (NA)
    pmin(pmax(object$pi0*df(object$data$tstat, object$data$df1, object$data$df2)/fitted(object), 0), 1)
}
lfdr.parncpt=#ppee.parncpt=
lfdr.nparncpt=#ppee.nparncpt=
lfdr.discTMix=#ppee.discTMix=
function(object, ...)
{
    if(any(is.na(object))) return (NA)
    pmin(pmax(object$pi0*dt(object$data$tstat, object$data$df)/fitted(object), 0), 1)
}
lfdr.sparncpt=#ppee.sparncpt=
function(object, ...)
{
    pmin(pmax(object$pi0*dt(object$parfit$data$tstat, object$parfit$data$df)/fitted(object), 0), 1)
}
lfdr.sparncpF=#ppee.sparncpt=
function(object, ...)
{
    pmin(pmax(object$pi0*df(object$parfit$data$Fstat, object$parfit$data$df1, object$parfit$data$df2)/fitted(object), 0), 1)
}

lfdr.nparncpp=#ppee.nparncpp=
function(object, ...)
{
    object$LFDR
}

lfdr.CBUM=#ppee.CBUM=
function(object, ...)
{
    attr(object, 'lfdr')
}

lfdr.znormix=#ppee.znormix=
function(object, ...)
{
    attr(object, 'lfdr')
}

lfdr.convest=#ppee.convest=
function(object, ...)
{
    attr(object, 'lfdr')
}


