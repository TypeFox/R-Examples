svyPVlevel <-
#function(by, svydat, pvs, CATDEF=NA, addcountry=TRUE)
function(by, svydat, pvs, CATDEF, levlab=NA, right=TRUE, colN=FALSE)
# CATEDEF = similar to the 'breaks' argument in the cut() function.
##1 so if you want to create level1 between 0 and 100, level2 between 100 and 300 and level3 between 300 and 500 the input has to look like:
## c(0,100,300,500)
  # which creates half opened intervals (0,100] ; (100,300] ; (300,500] which means that the interval is closed on the left an open on the right. so 0 is NOT included in the first interval, but 100 is. 100 is not included in the second interval but 300 is. etc ...
# so it is better to give the following CATDEF statement if you want to create levels like mentioned in ##1
  ## c(-0.1,100,300,501) --> or something similar
## colN --> add colnames
  
{

  
if(any(is.na(levlab)))
{
CATDEFlab <- paste0("level",1:(length(CATDEF)-1))
  
} else {
  
  #if(length(grep("\\.",levlab)) > 0) {"Check levels. No '.' allowed!"}
  
  CATDEFlab  <- levlab
}
  

PVPRS <- svydat$variables[, pvs]


LEVPVs <- data.frame(lapply(PVPRS, function(x) cut(x,CATDEF,CATDEFlab, right=right)))
colnames(LEVPVs) <- paste0(colnames(PVPRS),"_cat")

svydat2 <- 0
addALL <- paste0("svydat2 <- update(svydat, ",paste0(colnames(LEVPVs),"=LEVPVs$",colnames(LEVPVs),collapse=" , "),")")


# create svydat2 which includes the categorized plausible values
eval(parse(text=addALL))  


ergperc <- pv_perc(by=by, svydat=svydat2, pvcat=colnames(LEVPVs))
#ergperc


pm  <- mergeALL(ergperc)


### um die ordnung der factors gleich zu lassen (vor allem wichtig bezogen auf grafiken) wird hier nochmal umgeordnet so wie es im datensatz ?blich ist

mybys <- all.vars(by)
# facordall <- mapply(function(x,number) factor(pm[,number], levels=levels(svydat2$variables[[x]])), x=mybys, number=1:length(mybys),SIMPLIFY=FALSE)
# 
# facordallDF <- data.frame(facordall)

facordallDF <- fALL(mybys,pm, svydat)


pm[,1:length(mybys)] <- facordallDF 

# Beschriftung
if(colN)
{
colnames(pm)[1:(length(mybys)+1)] <- c(mybys,"levels") 

}

# if(addcountry)
# {  
#   pm  <- data.frame("Country"=unique(svydat$variables$CNTRYID), pm) 
# }

return(pm)  
  
}
