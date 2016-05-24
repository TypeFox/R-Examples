#==============================================================================
# File: strum.R
#
# Author: Nathan Morris
#         Yeunjoo Song
#
# Notes: Main driver of strum analysis.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# The main analysis given StrumModel and input data
#------------------------------------------------------------------------------
strum = function(myStrumModel, myStrumData,
                 step1OptimControl = list(maxit=5500, fnscale=-10),
                 startValueControl = list(initPopulation=NULL, nChildren=NULL, nGenerations=NULL, selection1=NULL, selection2=NULL),
                 step2OptimControl = list(maxit=5000, reltol=.Machine$double.eps),
                 ibdMarkers = NULL)
{
  # 1. Check input parameters
  #---------------------------
  if( class(myStrumModel) != "strumModel" )
    stop("Wrong class input! Need an object of 'strumModel' class.")

  if( class(myStrumData) != "strumData" )
    stop("Wrong class input! Need an object of 'strumData' class.")

  isPedigree = (dataType(myStrumData) == "Pedigree")

  # 2. Check ascertainment
  #------------------------
  probands = list()

  if( isPedigree )
  {
    aName = ascertainment(myStrumModel)

    if( class(aName) == "character" )
    {
      if( !(aName %in% names(dataVals(myStrumData))) )
      {
        stop(paste("Model contains ascertainment but no column of that name exist in data. ",
                   "Either set ascertainment in the model to NULL",
                   "or use a data set with the required proband column.",
                   sep = " "))
      } else
      {
        aValues = dataVals(myStrumData)[,aName, drop=FALSE]
        if( any(is.na(aValues)) )
        {
          stop(paste("Data contains the missing values(s) in the column '",
                     aName, "'. ",
                     "Complete data is required to run a strum model with ",
                     "ascertainment. Please check your data!",
                     sep = ""))
        }

        aValues = split(aValues, dataVals(myStrumData)$family)

        probands = lapply(aValues, function(pk) return(pk==1))

        pCount = unlist(lapply(aValues, sum))
        if( any(pCount > 1) )
          warning("Pedigree(s) with more than 1 proband exist!")
        if( any(pCount== 0) )
          warning("Pedigree(s) with no proband exist!")
      }

    } else if( !is.null(aName) )
    {
      warning("Model contains a wrongly specified ascertainment!  Ignoring...")
    }
  }

  # 3. Check random effect
  #------------------------
  allRE = allRandomEffects(myStrumModel)

  if( any(c("a","p","c") %in% allRE) & !isPedigree )
    stop("Can't fit a genetic, pologenic or common environment effect without pedigree structures.")

  cat("\n Start STRUM analysis ...\n")

  # 3. Prepare analysis data
  #--------------------------
  y     = .getAnalysisY( myStrumModel, myStrumData, probands)
  x     = .getAnalysisX( myStrumModel, myStrumData, probands)
  vcPEC = .getAnalysisVC(myStrumModel, myStrumData, allRE)

  myFittedModel = list()

  if( any(allRE == "a") )
  {
    iMarkers = NULL
    
    if( is.null(ibdMarkers) )
    {
      iMarkers = names(ibdMatrix(ibd(myStrumData)))
    } else
    {
      for( m in 1:length(ibdMarkers) )
      {
        mName = ibdMarkers[m]

        if( !(mName %in% names(ibdMatrix(ibd(myStrumData)))) )
          warning(paste("IBD marker \"", mName, "\" doesn't exist in the IBD file!  Ignoring...",
                        sep=""))
        else
          iMarkers = c(iMarkers, mName)
      }

      if( is.null(iMarkers) )
      {
        stop(paste("No valid IBD marker available in 'ibdMarkers'!  ",
                   "Will use all markers in 'myStrumData' object.",
                   sep=""))

        iMarkers = names(ibdMatrix(ibd(myStrumData)))
      }
    }

    if( is.null(iMarkers) )
    {
      stop(paste("No IBD marker available in 'myStrumData' object!  Need an object of ",
                 "strumData class constructed with the IBD information of the marker(s) ",
                 "to be included in the model as additive genetic component(s).",
                 sep=""))
    }

    for( m in 1:length(iMarkers) )
    {
      mName = iMarkers[m]

      cat(paste("\n  ", mName, ":\n", sep=""))
      vca = ibdMatrix(ibd(myStrumData))[[mName]]

      vcAll = mapply(function(peci, vcai)
                     {
                       peci[["a"]] = vcai
                       return(peci)
                     },
                     vcPEC, vca,
                     SIMPLIFY = FALSE)

      vcPro = list()
      if( length(probands) > 0 )
        vcPro = mapply(.filterMissingVC, vcAll, probands, probands, SIMPLIFY=FALSE) 

      vc = list(vcAll = vcAll, vcPro = vcPro)
      filtered = .getValidAnalysisData(y, x, vc)

      mFittedModel = .fitModel(myStrumModel, filtered$y, filtered$x, filtered$vc, step1OptimControl, startValueControl, step2OptimControl)

      myFittedModel[[mName]] = mFittedModel
    }
  } else
  {
    vcPro = list()
    if( length(probands) > 0 )
      vcPro = mapply(.filterMissingVC, vcPEC, probands, probands, SIMPLIFY=FALSE) 

    vc = list(vcAll = vcPEC, vcPro = vcPro)
    filtered = .getValidAnalysisData(y, x, vc)

    myFittedModel = .fitModel(myStrumModel, filtered$y, filtered$x, filtered$vc, step1OptimControl, startValueControl, step2OptimControl)
  }

  cat("\nAnalysis completed!\n\n")

  return(myFittedModel)
}

#------------------------------------------------------------------------------
# Model analysis data Y from StrumData given StrumModel
# - Extract Y
#------------------------------------------------------------------------------
.getAnalysisY = function(myStrumModel, myStrumData, probands)
{
  yNames = myStrumModel@varList$name[myStrumModel@varList$inY==TRUE]

  # Error check if names not there
  if( !length(yNames) )
    stop("No y names exist in StrumModel")

  yAll = tryCatch(dataVals(myStrumData)[,yNames, drop=FALSE],
                  error=function(e)
                        {          
                          stop(paste("Some of variables expected from the model ",
                                     "are not in the data.  Please check the ",        
                                     "variable names in your model and data!",          
                               sep=""))          
                        })          

  yAll = split(yAll, dataVals(myStrumData)$family)
  yAll = lapply(yAll, data.matrix)
  
  yPro = list()
  if( length(probands) > 0 )
    yPro = mapply(.filterMissing, yAll, probands, TRUE, SIMPLIFY=FALSE) 

  return(list(yAll = yAll, yPro = yPro))
}

#------------------------------------------------------------------------------
# Model analysis data X from StrumData given StrumModel
# - Extract X from covariates
#------------------------------------------------------------------------------
.getAnalysisX = function(myStrumModel, myStrumData, probands)
{
  xNames = myStrumModel@varList$name[myStrumModel@varList$covariate==TRUE]

  xAll = list()

  # If names not there, then intercept only
  if( !length(xNames) )
  {
    xAll = rep(1, nrow(dataVals(myStrumData)))
    xAll = split(xAll, dataVals(myStrumData)$family)
    xAll = lapply(xAll, data.matrix)
  } else
  {
    xAll = tryCatch(dataVals(myStrumData)[,xNames, drop=FALSE],
                    error=function(e)
                          {          
                            stop(paste("Some of covariates expected from the model ",
                                       "are not in the data.  Please check the ",        
                                       "covariate names in your model and data!",          
                                 sep=""))          
                          })          

    na_default = options("na.action")$na.action
    options(na.action='na.pass')
    xAll = lapply(split(xAll, dataVals(myStrumData)$family),
                  function(xi)
                  {
                    xi = data.matrix(xi)
                    model.matrix(~xi)
                  })
    options(na.action=na_default)
  }
  
  xPro = list()
  if( length(probands) > 0 )
    xPro = mapply(.filterMissing, xAll, probands, TRUE, SIMPLIFY=FALSE) 

  return(list(xAll = xAll, xPro = xPro))
}

#------------------------------------------------------------------------------
# Model analysis data VC from StrumData given StrumModel
# - Use allRandomEffects to extract VC values
#------------------------------------------------------------------------------
.getAnalysisVC = function(myStrumModel, myStrumData, allRE)
{
  vcPEC = list()

  if( any(allRE == "p") )
    vcPEC = lapply(phi(myStrumData), function(p) return(list(p=p)))

  if( any(allRE == "e") )
  {
    vce = lapply(split(dataVals(myStrumData), dataVals(myStrumData)$family),
                 function(p) return(list(e=diag(nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcei) return(c(vci, vcei)),
                     vcPEC, vce,
                     SIMPLIFY = FALSE)
    else
      vcPEC = vce
  }

  if( any(allRE == "c") )
  {
    vcc = lapply(phi(myStrumData), 
                 function(p) return(list(c=matrix(rep(1,nrow(p)*nrow(p)),nrow=nrow(p)))))

    if( length(vcPEC) )
      vcPEC = mapply(function(vci, vcci) return(c(vci, vcci)),
                     vcPEC, vcc,
                     SIMPLIFY = FALSE)

    else
      vcPEC = vcc
  }

  return(vcPEC)
}

#------------------------------------------------------------------------------
# Get valid analysis data only - take out empty ped after filtering by X & Y
#------------------------------------------------------------------------------
.getValidAnalysisData = function(y, x, vc)
{
  missingX = lapply(x$xAll, .findMissing, isX=TRUE)
  missingY = lapply(y$yAll, function(yk) return(rowSums(is.na(yk))!= ncol(yk)))

  missingXY = mapply("&", missingX, missingY, SIMPLIFY=FALSE)

  xa  = mapply(.filterMissing,   x$xAll,   missingXY, TRUE,      SIMPLIFY=FALSE)
  ya  = mapply(.filterMissing,   y$yAll,   missingXY, TRUE,      SIMPLIFY=FALSE)
  vca = mapply(.filterMissingVC, vc$vcAll, missingXY, missingXY, SIMPLIFY=FALSE)

  xan = Filter(function(x) nrow(x) > 0, xa)
  ped_names = names(xan)
  yan  = ya[ped_names]
  vcan = vca[ped_names]

  ypn = list()
  xpn = list()
  vcpn = list()

  if( length(y$yPro) > 0 )
  {
    missingXp = lapply(x$xPro, .findMissing, isX=TRUE)
    missingYp = lapply(y$yPro, function(yk) return(rowSums(is.na(yk))!= ncol(yk)))
    
    missingXYp = mapply("&", missingXp, missingYp, SIMPLIFY=FALSE)
    
    xp  = mapply(.filterMissing,   x$xPro,   missingXYp, TRUE,       SIMPLIFY=FALSE)
    yp  = mapply(.filterMissing,   y$yPro,   missingXYp, TRUE,       SIMPLIFY=FALSE)
    vcp = mapply(.filterMissingVC, vc$vcPro, missingXYp, missingXYp, SIMPLIFY=FALSE)
    
    ypn  = yp[ped_names]
    xpn  = xp[ped_names]
    vcpn = vcp[ped_names]
  }


  return(list(y  = list(yAll=yan,yPro=ypn),
              x  = list(xAll=xan,xPro=xpn),
              vc = list(vcAll=vcan,vcPro=vcpn)))
}
