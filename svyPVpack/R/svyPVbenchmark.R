svyPVbenchmark <-
function(by, svydat, pvs, BENCH=NA, colN=FALSE)
# put in a Benachmark - a this benchmark (e.g.: 300) the ratio of people below and at/above this benchmark is identified and reported - as well as the SE of this ratio
  # additionally the sum of weights and the number of observations for each group is reported.

# pvs = a vertor which contains the column-names of the plausible value variables  
{
  

if(is.na(BENCH)) stop("Please submit a benchmark!")



BENCH1 <- paste0("x >=",BENCH)

PVPRS <- svydat$variables[, pvs]


catPVs <- data.frame(lapply(PVPRS,function(x) factor(ifelse(eval(parse(text=BENCH1)),paste0(">=",BENCH),paste0("<",BENCH)))))
colnames(catPVs) <- paste0(colnames(PVPRS),"_cat")

svydat2 <- 0
addALL <- paste0("svydat2 <- update(svydat, ",paste0(colnames(catPVs),"=catPVs$",colnames(catPVs),collapse=" , "),")")


# create svydat2 which includes the categorized plausible values
eval(parse(text=addALL))  
  
  

ergperc <- pv_perc(by=by, svydat=svydat2, pvcat=colnames(catPVs))
pm      <- mergeALL(ergperc)



### um die ordnung der factors gleich zu lassen (vor allem wichtig bezogen auf grafiken) wird hier nochmal umgeordnet so wie es im datensatz ?blich ist

mybys <- all.vars(by)
# facordall <- mapply(function(x,number) factor(pm[,number], levels=levels(svydat2$variables[[x]])), x=mybys, number=1:length(mybys),SIMPLIFY=FALSE)
# 
# facordallDF <- data.frame(facordall)

browser()

facordallDF <- fALL(mybys,pm, svydat)


browser()

pm[,1:length(mybys)] <- facordallDF 

# Beschriftung
if(colN)
{
  colnames(pm)[1:(length(mybys)+1)] <- c(mybys,"benchmark") 
  
}




# if(addcountry)
# {  
#   pm  <- data.frame("Country"=unique(svydat$variables$CNTRYID), pm) 
# }
  
return(pm)  
  
}
