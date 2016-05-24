svyPVquantile <-
function(by, svydat, pvs, quantile, interval.type="quantile", colN=FALSE, ...)
{

ergquant <- pv_quantile(by=by, svydat=svydat, pvs=pvs, quantile=quantile,  interval.type=interval.type, ...)
pm       <- mergeALL(ergquant)  



### um die ordnung der factors gleich zu lassen (vor allem wichtig bezogen auf grafiken) wird hier nochmal umgeordnet so wie es im datensatz ?blich ist

mybys <- all.vars(by)
# facordall <- mapply(function(x,number) factor(pm[,number], levels=levels(svydat$variables[[x]])), x=mybys, number=1:length(mybys),SIMPLIFY=FALSE)
# 
# facordallDF <- data.frame(facordall)

facordallDF <- fALL(mybys,pm, svydat)

pm[,1:length(mybys)] <- facordallDF


if(colN)
{
  colnames(pm)[1:length(mybys)] <- c(mybys) 
  
}




# if(addcountry)
# {  
#   pm  <- data.frame("Country"=unique(svydat$variables$CNTRYID), pm) 
# }

pm  
}
