#==============================================================================
# File: generateStartValue.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Generate the starting values to fit the model in step2.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Generate the starting values, step 2
#------------------------------------------------------------------------------
.generateStartValue = function(model, delta, w, startValueControl)
{
  Q = function(theta)
      {
        diff = delta - .thetaToDelta(model, theta)
        return(sum(diff*diff*w))
      }

  positive = grep("<", paramNames(model))

  initialPopulation = length(paramNames(model))*120
  if( !is.null(startValueControl[["initPopulation"]]) )
    initialPopulation = startValueControl[["initPopulation"]]

  startVals = lapply(1:initialPopulation,
                     function(dummy)
                     {
                       startVal = rnorm(length(paramNames(model)))
                       startVal[positive] = abs(startVal[positive])
                       return(list(startVal=startVal, Q=Q(startVal)))
                     })

  nChildren         = length(paramNames(model))*2
  nGenerations      = 20
  selection1        = 15
  selection2        = 15

  if( !is.null(startValueControl[["nChildren"]]) )
    nChildren = startValueControl[["nChildren"]]

  if( !is.null(startValueControl[["nGenerations"]]) )
    nGenerations = startValueControl[["nGenerations"]]

  if( !is.null(startValueControl[["selection1"]]) )
    selection1 = startValueControl[["selection1"]]

  if( !is.null(startValueControl[["selection2"]]) )
    selection2 = startValueControl[["selection2"]]

  QS = sapply(startVals, function(x) return(x$Q))
  startVals = startVals[order(QS)[1:selection1]]

  Delta = 0 != .funDeriv(fun = function(theta) .thetaToDelta(model,theta),
                        startVals[[1]]$startVal )

  for( gen in 1:nGenerations )
  {
    newGen = list()
    keepTrack = rep(0,length(startVals))

    for( i in 1:length(startVals) )
    {
      deltaTemp = sign(.thetaToDelta(model, startVals[[i]]$startVal)) != sign(delta)
      keepTrack[i] = length(deltaTemp)-sum(deltaTemp)
      deltaTemp = (Delta%*%deltaTemp)
      gammai = which(deltaTemp>=1)
      gammai = gammai[!(gammai %in% positive)]
      stVali = startVals[[i]]$startVal
      stValiGammai = stVali[gammai]

      for( child in 1:nChildren )
      {		
        newgenTmp = stVali
        newgenTmp[gammai] = stValiGammai * (rbinom(length(gammai), 1, 0.5)-.5)*2
        newgenTmp = newgenTmp+rnorm(length(stVali), sd=.01)	
        deltaTemp = sign(.thetaToDelta(model, newgenTmp)) != sign(delta)
        newGen[[child]] = list(startVal=newgenTmp, 
                               Q=Q(newgenTmp)+ sum(deltaTemp)*1000)
      }

      QS = sapply(newGen, function(x) return(x$Q))
      startVals[i] = newGen[order(QS)[1]]
    }

    #cat(gen, keepTrack, "\n")
    tbl = as.matrix(table(keepTrack))
    tmp = as.integer(rownames(tbl))
    maxrow = which(tmp==max(tmp))

    if( tbl[maxrow]/sum(tbl) > .85 )
      break
  }

  QS = sapply(startVals, function(x) return(x$Q))
  startVals = startVals[order(QS)[1:selection2]]

  return(startVals)
}
