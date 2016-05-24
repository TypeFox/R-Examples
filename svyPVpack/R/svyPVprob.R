svyPVprob <-
function(by, svydat, pvs=NULL, colN=FALSE)
{
# if plausible values (or any other column names) are provided by the argument "pvs" the weighted ratio of group members is computed only for the observations where all pvs variables are NOT NA!  "pvs" must contain valid names of variables which exist in the "svydat" survey-design!

# check input

checkds <- data.frame(svydat$variables[,all.vars(by)])
only1gr <- all(sapply(checkds,function(A)length(unique(A))) == 1)  
if(only1gr) stop("by variable must contain more th an 1 category.")


# compute group ratios
cellp   <- opv_perc(by=by, svydat=svydat, pvs=pvs)


#########################################################
########### additional information ######################
######################################################### 
ADC <- additional_comp(by=by,svydat=svydat)

colnames(ADC$Ncases) <- c(paste0("Group",1:length(all.vars(by))),"Number.of.cases")
colnames(ADC$Sumweights) <- c(paste0("Group",1:length(all.vars(by))),"Sum.of.weights")


pmV <- data.frame(ADC$Ncases, "Sum.of.weights"=ADC$Sumweights[,length(all.vars(by))+1])


# merge the outcomes reasonably
pm  <- merge(pmV,cellp,sort=FALSE)



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

return(pm)
  
}
