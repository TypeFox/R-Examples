

calibratePackageSurvey = function(grossSample,modelformula,popTotals=NULL,response=NULL,
          popData=NULL,samplingWeights=NULL,calfun="linear",totalReturn=0,returnwGross=FALSE,...)
{
  if(returnwGross) calWeights1 = NULL
  modelformula =as.formula(modelformula)

  lm_model = lm(modelformula,data=grossSample)
  xFromModel = model.matrix(lm_model)
  rFromModel = model.frame(lm_model)[1]

  if(is.null(popTotals))
  {
    if(is.null(popData)) tp = getTotal(grossSample,lm_model,samplingWeights)   # Ny
	else tp = getTotal(popData,lm_model)

	if(totalReturn==2) return(tp)
    popTotals =  tp$N*tp$colSum/tp$colN
  }
  else
  {
    popTotals=setTotal(popTotals,lm_model)
  }

  if(totalReturn==1) return(popTotals)

  if(sum(is.na(popTotals))>0)
  {
    # not netDesign but same name
    netDesign = svydesign(ids=~1,data=data.frame(rFromModel,xFromModel),weights=samplingWeights)
    col1= !is.na(popTotals)
    col1[1] = FALSE     # Antar Intercept nr 1 og denne tas bort
    varNames1 = colnames(xFromModel)[col1] # Variabler der totaler finnes
    modelformula1 = update(modelformula,paste("~ ",paste(varNames1,collapse="+")))
    calWeights1 = weights(calibrate(netDesign,modelformula1 ,popTotals[!is.na(popTotals)],calfun=calfun,...))
    popTotals[is.na(popTotals)] =  colSums(calWeights1*xFromModel[,is.na(popTotals)])
  }

  if(!(totalReturn==0)) return(popTotals)

  netDesign = svydesign(ids=~1,data=grossSample[rFromModel==1,],weights=samplingWeights[rFromModel==1])

  calWeights = rep(0,dim(grossSample)[1])
  calWeights[rFromModel==1] = weights(calibrate(netDesign,modelformula ,popTotals,calfun=calfun,...))
  if(returnwGross){
    return(list(w=calWeights,wGross=calWeights1))
  }
  calWeights
}
