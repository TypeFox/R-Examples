#==============================================================================
# File: simulateData.R
#
# Author: Nathan Morris
#
# Notes:
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Simulate strumData object with given information.
#------------------------------------------------------------------------------
simulateStrumData = function(simModel, inData = NULL, N = NULL)
{
  if( class(simModel) != "strumSimModel" )
    stop("'simModel' must be of class strumSimModel")

  # Do error checking and marker simulation
  #  depending on what kind of inData was passed in
  #-------------------------------------------------------------------
  if( class(inData) == "NULL" )
  {
    if( is.null(N) )
      stop("If inData is not specified, N must be specified!")
    else if( N <= 0 )
      stop("N must be a positive integer!")

    if( sum(varList(simModel)$covariate) != 0 )
      stop(paste("Model contains covariates but inData=NULL. ",
                 "Either add an additional formula to change the",
                 "covariate into an endogenous (i.e. dependent) variable,",
                 "or pass in a data set with the required covariates.",
                 sep = " "))

    #fam = rep(1, N)
    id  = seq(1:N)
    fm  = rep(0, N)
    x   = data.frame(cbind(id, id, fm, fm))
    names(x) = c("family", "id", "father", "mother")

    retObj = createStrumData(x, "RawData")
    isPedigree = FALSE
	
  } else if( class(inData) == "data.frame" )
  {
    if( !length(inData) )
      stop("inData must be non-empty data.frame!")

    # Create new strumData object
    #-----------------------------
    if( .checkPed(inData) == 0 )
    {
      retObj = createStrumData(inData, "Pedigree")
      isPedigree = TRUE
    } else 
    {
      retObj = createStrumData(inData, "RawData")
      isPedigree = FALSE		
    }

    if( is.null(N) ) 
      N = nrow(dataVals(retObj))
    else if( is.null(simModel@ascertainment) )
    {
      warning("Ignoring N becuase data has been input but ascertainment has not been specified.")
      N = nrow(dataVals(retObj))
    }

  } else if( class(inData) == "strumData" )
  {
    if( !length(dataVals(inData)) )
      stop("inData strumData object must contain non-empty data!")

    if( dataType(inData) == "RawData" )
      isPedigree = FALSE
    else if( dataType(inData) == "Pedigree" )
      isPedigree = TRUE
    else
      stop(paste("inData is a strumData object of type ",
		 dataType(inData),
                 ". It must be of type 'RawData' or 'Pedigree'.",
                 sep = ""))

    retObj = inData

    if( is.null(N) ) 
      N = nrow(dataVals(retObj))
    else if( is.null(simModel@ascertainment) )
    {
      warning("Ignoring N becuase data has been input but ascertainment has not been specified.")
      N = nrow(dataVals(retObj))
    }

  } else
  {
    stop("inData must be a data.frame, strumData object or NULL.")
  }

  inDFVarNames = names(dataVals(retObj))

  if( any(c("p","c") %in% allRandomEffects(simModel)) & !isPedigree )
    warning("Simulating a polygenic or common environment effect without pedigree structures.")

  if( !is.null(markerInfo(simModel)) )
    inDFVarNames = c(inDFVarNames, rownames(haplotypes(markerInfo(simModel))))

  # Calculate covariance matrix
  #-----------------------------
  theta = numeric()
  L     = simModel@L(theta)
  B     = simModel@B(theta)
  E     = lapply(simModel@E, function(f) return(f(theta)))
  Z     = lapply(simModel@Z, function(f) return(f(theta)))
  Gs    = simModel@Gs(theta)
  Gm    = simModel@Gm(theta)
  H     = solve(diag(nrow(B))-B)
  LH    = L %*% H
  C     = Gm + LH %*% Gs
  VC    = lapply(simModel@allRandomEffects,
                 function(v) return(LH %*% Z[[v]] %*% t(LH) + E[[v]]))

  names(VC) = names(Z)

  # also dataVals(inData) is a data.frame containing the input phenotypes/covariates etc.
  #--------------------------------------------------------------------------------------
  myVarList = varList(simModel)

  if( nrow(dataVals(retObj)) != 0 )
  {
    covPhenos = colnames(C)

    if( length(covPhenos) > 0 )
    {
      # Make sure that the needed covariates are present in dataVals(retObj) 
      # (note that marker names are checked to exist when the model object is created)
      #--------------------------------------------------------------------------------
      if( !all(covPhenos %in% c(inDFVarNames, "[intercept]")) )
      {
        badones = covPhenos[!(covPhenos %in% c(inDFVarNames,"[intercept]"))]
        stop(paste("The model contains covariates that are not present in the 'inData' parameter.\n",
                   "Either add an additional formula to change the covariate into an endogenous\n",
                   "(i.e. dependent) variable, or pass in a data set with the required covariates.\n",
                   "Covariates missing from inData are:",
                   paste(badones, collapse = ", ")))
      }
    }
  }

  if( isPedigree )
  {
    if( "p" %in% allRandomEffects(simModel) )
      phisqrt = lapply(phi(retObj), .sqrtMat)
    else
      phisqrt = NULL

    if( is.null(ascertainment(simModel)) )
    {
      # simulate markers
      #------------------
      if( !is.null(markerInfo(simModel)) )
        retObj = .simulateMarker(markerInfo(simModel), retObj)

      X = as.matrix(cbind("[intercept]"=rep(1,nrow(dataVals(retObj))), dataVals(retObj)[,covPhenos[-1]]))
      myMean = tcrossprod(X, C)

      simdata = .simulateVC(VC, myMean, N, phiSqrt=phisqrt, pedPhi=phi(retObj))
      simdata = .simulateMissing(simdata, traitMissingRate(simModel))
      simdata = cbind(simdata, dataVals(retObj))
      retObj@dataVals = simdata
    } else
    {
      simdata = matrix(0,0,0)
      phiMat  = NULL

      needIBD = !is.null(ibd(retObj))
      if( !is.null(markerInfo(simModel)) )
        needIBD = !is.null(ibd(retObj)) | returnIBD(markerInfo(simModel))

      if( needIBD )
        ibdMat = NULL

      maxFamilyIndex = 0
      shortN = nrow(dataVals(retObj))

      counter = 1
      while( nrow(simdata) < N && counter < (shortN*1000) )
      {  
        # keep simulating until samplesize is greater than N
        #----------------------------------------------------
        if( !is.null(markerInfo(simModel)) )
          simObj = .simulateMarker(markerInfo(simModel), retObj)
        else
          simObj = retObj

        X = as.matrix(cbind("[intercept]"=rep(1,nrow(dataVals(retObj))), dataVals(simObj)[,covPhenos[-1]]))
        myMean = tcrossprod(X, C)

        simdataTmp = .simulateVC(VC, myMean, N=shortN, phiSqrt=phisqrt, pedPhi=phi(simObj))
        simdataTmp = data.frame(simdataTmp, dataVals(simObj))
        simdataTmp$family = simdataTmp$family + maxFamilyIndex
        simObj@dataVals = simdataTmp
        fname = as.character(as.numeric(names(phi(simObj))) + maxFamilyIndex)
        names(phi(simObj)) = fname

        for( m in as.character(ibdMarker(ibd(simObj))$Marker) )
          names(simObj@ibd@ibdMatrix[[m]]) = fname

        probandExist = FALSE

        ascInfo = lapply(split(simdataTmp, dataVals(retObj)$family),
                         function(thisFam)
                         {
                           tmp = simModel@ascertainment(thisFam)

                           if( is.list(tmp) )
                           {
                             if( length(tmp) != 2 )
                             {
                               stop(paste("Return value(s) from ascertainment function of the ",
                                          "simulation model is not correct!  The first component ",
                                          "should be a TRUE/FALSE value indicating the ascertainment ",
                                          "of the pedigree.  The second component should be a vector ",
                                          "of TRUE/FALSE, the same length as the size of pedigree, ",
                                          "stating the proband status of each member of the pedigree.", 
                                          sep=""))
                             }

                             if( length(tmp[[1]]) != 1 || !is.logical(tmp[[1]]) || nrow(thisFam) != length(tmp[[2]]) )
                             {
                               stop(paste("Return value(s) from ascertainment function of the ",
                                          "simulation model is not correct!  The first component ",
                                          "should be a TRUE/FALSE value indicating the ascertainment ",
                                          "of the pedigree.  The second component should be a vector ",
                                          "of TRUE/FALSE, the same length as the size of pedigree, ",
                                          "stating the proband status of each member of the pedigree.", 
                                          sep=""))
                             }

                             probandExist <<- TRUE
                             tmp[["keep"]] = rep(tmp[[1]], nrow(thisFam))
                             names(tmp) = c("aStatus", "proband", "keep")

                             return(tmp)
                           }

                           if( !is.logical(tmp) )
                             stop(paste("Return value(s) from ascertainment function of the ",
                                        "simulation model is not correct!  It has to be a TRUE/FALSE ",
                                        "value indicating the ascertainment of the pedigree.",
                                        sep=""))

                           if( length(tmp) == nrow(thisFam) )
                           {
                             probandExist <<- TRUE
                             return(list(aStatus = any(tmp),
                                         proband = tmp,
                                         keep    = rep(any(tmp), nrow(thisFam))))
                           } else if( length(tmp) == 1 )
                           {
                             return(list(aStatus = tmp,
                                         proband = rep(FALSE, nrow(thisFam)),
                                         keep    = rep(tmp, nrow(thisFam))))
                           }
			 })

        if( probandExist )
        {
          pStatus = lapply(ascInfo, function(aInfo) return(aInfo$proband))

          pCount = unlist(lapply(pStatus, sum))
          if( any(pCount > 1) )
            warning("Pedigree(s) with more than 1 proband exist!")

          proband = as.numeric(unlist(pStatus))
          simdataTmp = cbind(simdataTmp, proband)
        }

        keep = unlist(lapply(ascInfo, function(aInfo) return(aInfo$keep)))

        if( sum(keep) > 0 )
        {
          simObj = simObj[keep]

          if( nrow(simdata) == 0 )
          {
            simdata = simdataTmp[keep,]
            phiMat  = phi(simObj)

            if( needIBD )
              ibdMat = ibd(simObj)

          } else
          {
            simdata = rbind(simdata, simdataTmp[keep,])
            phiMat = c(phiMat, phi(simObj))

            if( needIBD )
              ibdMat@ibdMatrix = mapply(function(a,b) return(c(a,b)),
                                        ibdMatrix(ibdMat),
                                        ibdMatrix(ibd(simObj)),
                                        SIMPLIFY=FALSE)
          }

          maxFamilyIndex = max(simdataTmp$family)
        }
      
        counter = counter + 1
      }

      if( counter >= shortN*1000 )
        stop("Ascertainment probability too close to 0, no probands exist with 1000 replicates!")

      simdata = .simulateMissing(simdata, traitMissingRate(simModel))

      retObj@dataVals = simdata
      retObj@phi = phiMat

      if( needIBD )
        retObj@ibd = ibdMat

      if( nrow(simdata) > N )
      {	
        lastFam = simdata$family[N]
        lastIndividualPos = max(which(lastFam==simdata$family))
        retObj = retObj[1:lastIndividualPos,]
      }
    }
  } else
  {
    #myVar = VC[[1]]

    #if( length(VC) > 1 )
    #  for( i in 2:length(VC) )
    #    myVar = myVar + VC[[i]]
    if( "p" %in% allRandomEffects(simModel) )
      phi = lapply(split(dataVals(retObj), dataVals(retObj)$family),
                   function(p) return(diag(nrow(p))))
    else
      phi = NULL

    if( is.null(ascertainment(simModel)) )
    {
      # simulate markers
      #------------------
      if( !is.null(markerInfo(simModel)) )
        retObj = .simulateMarker(markerInfo(simModel), retObj)

      X = as.matrix(cbind("[intercept]"=rep(1,nrow(dataVals(retObj))), dataVals(retObj)[,covPhenos[-1]]))
      myMean = tcrossprod(X, C)

      # simulate the phenotypes
      #-------------------------
      #sqrtVar = .sqrtMat(myVar, method="c")
      #simdata = myMean + matrix(rnorm(N*nrow(sqrtVar)),N,nrow(sqrtVar)) %*% sqrtVar
      simdata = .simulateVC(VC, myMean, N, phiSqrt=phi, pedPhi=phi)
      simdata = .simulateMissing(simdata, traitMissingRate(simModel))
      simdata = cbind(simdata, dataVals(retObj))
        
      retObj@dataVals = simdata
    } else
    {
      stop("Ascertainment simulation has only been implemented for raw data!")
    }
  }

  .printInfoLine("Simulating strumData", "Done", 52, 0)

  return(retObj)
}

#------------------------------------------------------------------------------
# Simulate missing
#------------------------------------------------------------------------------
.simulateMissing = function(simdata, trtMissingRate)
{
  if( all(trtMissingRate == 0) )
    return(simdata)

  for( t in 1:ncol(simdata) )
  {
    mPos = sample.int(length(simdata[,t]), length(simdata[,t])*trtMissingRate[t], replace=T)
    simdata[mPos,t] = NA
  }

  return(simdata)
}

#------------------------------------------------------------------------------
# Simulate variance components
#------------------------------------------------------------------------------
.simulateVC = function(VC, myMean, N, phiSqrt = NULL, pedPhi =  NULL)
{
  retVal = myMean

  if( "p" %in% names(VC) )
  {
    sqrtVar = .sqrtMat(VC[["p"]])

    polyComponent = lapply(phiSqrt,
                           function(kSqrtPhi)
                           {
                             P = matrix(rnorm(nrow(kSqrtPhi)*nrow(sqrtVar)),nrow=nrow(kSqrtPhi))
                             return(kSqrtPhi %*% P %*% sqrtVar)
                           })

    polyComponent = do.call(rbind,polyComponent)
    retVal = retVal + polyComponent
  }

  if( "e" %in% names(VC) )
  {
    sqrtVar = .sqrtMat(VC[["e"]])
    retVal  = retVal + matrix(rnorm(N*nrow(sqrtVar)),N,nrow(sqrtVar)) %*% sqrtVar
  }

  if( "c" %in% names(VC) )
  {
    sqrtVar = .sqrtMat(VC[["c"]])
    nPeds   = length(pedPhi)
    tmp     = unlist(lapply(1:nPeds, 
                            function(k) return(rep(k,nrow(pedPhi[[k]])))))

    commonEnvComponent = matrix(rnorm(nPeds*nrow(sqrtVar)),nPeds,nrow(sqrtVar)) %*% sqrtVar

    commonEnvComponent = commonEnvComponent[tmp,]
    retVal = retVal + commonEnvComponent
  }

  return(retVal)
}

#------------------------------------------------------------------------------
# Square-root of a matrix
#------------------------------------------------------------------------------
.sqrtMat = function(S, method = "e")
{
  if( method == "e" )
  {
    e = eigen(S)
    return(e$vec %*% diag(sapply(e$val, function(x) sqrt(max(0, x)))) %*% t(e$vec))
  } else
  {
    return(chol(S))
  }
}

#------------------------------------------------------------------------------
# Performs the genetic part of the simulation.
#------------------------------------------------------------------------------
.simulateMarker = function(mInfo, sData)
{
  if( is.null(mInfo) )
    return(sData)

  if( dataType(sData) == "Pedigree" )
  {
    status = .checkPed(dataVals(sData))  #Check if a pedigree
    if( status )
      stop("SimulateMarker: Data is not a pedigree.\n")

    pedData = .orderPedigrees(data.frame(sDataTempRowId=1:nrow(dataVals(sData)),
                             dataVals(sData)))

    pos = .getParentPos(pedData)

    MFpos = cbind(absoluteMother=unlist(lapply(pos,
                                               function(posForFamily)
                                               {
                                                 return(posForFamily$absoluteMother)
                                               })),
                  absoluteFather=unlist(lapply(pos,
                                               function(posForFamily)
                                               {
                                                 return(posForFamily$absoluteFather)
                                               })))

    if( mInfo@returnIBD )
    {
      dummyMarkerInfo = .findDummyMarkerInfo(mInfo)
      dummyHaps = list()
      founderHapCount = 0
    }
  } else if( dataType(sData) == "RawData" )
  {
    if( returnIBD(mInfo) )
    {
      warning("IBD information can only be returned for a pedigree")
      mInfo@returnIBD = FALSE
    }

    MFpos = matrix(NA,nrow(dataVals(sData)),2)
  } else
    stop("Cannot simulate marker data when data type is not either 'Pedigree' or 'RawData'")

  ret = as.data.frame(matrix(0.0,nrow(MFpos), nrow(markerFacts(mInfo))))
  haps = list()
  hap1 = list()
  hap2 = list()

  for( i in 1:nrow(MFpos) )
  {
    #is haplotype1 a founder?
    if( is.na(MFpos[i,1]) )
    {
      hap1 = .simulateHaplotypes(mInfo)

      if( returnIBD(mInfo) )
      {
        founderHapCount = founderHapCount+1
        hap1$dummy = rep(founderHapCount, dummyMarkerInfo$dummyMarkerCountAll)
      }
    } else
    {
      if( returnIBD(mInfo) )
      {
        hap1 = .simulateHaplotypes(mInfo,
                                   offspring = TRUE,
                                   parentHaplotypes = haps[[MFpos[i,1]]],
                                   parentDummyHaplotypes = dummyHaps[[MFpos[i,1]]],
                                   dummyMarkerInfo = dummyMarkerInfo)
      } else
      {
        hap1 = .simulateHaplotypes(mInfo,
                                   offspring = TRUE,
                                   parentHaplotypes = haps[[MFpos[i,1]]])
      }
    }

    #is haplotype2 a founder?
    if( is.na(MFpos[i,2]) )
    {
      hap2 = .simulateHaplotypes(mInfo)

      if( returnIBD(mInfo) )
      {
        founderHapCount = founderHapCount+1
        hap2$dummy = rep(founderHapCount, dummyMarkerInfo$dummyMarkerCountAll)
      }
    } else
    {
      if( returnIBD(mInfo) )
      {
        hap2 = .simulateHaplotypes(mInfo,
                                   offspring = TRUE,
                                   parentHaplotypes = haps[[MFpos[i,2]]],
                                   parentDummyHaplotypes = dummyHaps[[MFpos[i,2]]],
                                   dummyMarkerInfo = dummyMarkerInfo)
      } else
      {
        hap2 = .simulateHaplotypes(mInfo,
                                   offspring = TRUE,
                                   parentHaplotypes = haps[[MFpos[i,2]]])
      }
    }

    haps[[i]] = cbind(hap1$main, hap2$main)

    if( returnIBD(mInfo) )
      dummyHaps[[i]] = cbind(hap1$dummy, hap2$dummy)

    ret[i,] = mInfo@coding[1+apply(haps[[i]],1,sum)] #code the two haplotypes
  }

  colnames(ret) = rownames(mInfo@haplotypes)

  if( returnIBD(mInfo) )
  {
    numbDummyMarkers = nrow(dummyHaps[[1]])
    dummyMarkerNames = unlist(dummyMarkerInfo$markerNames)
    dummyMarkerIndex = 1:numbDummyMarkers
    
    names(dummyMarkerIndex) = dummyMarkerNames
    
    familyIndex = tapply(1:nrow(pedData), pedData$family, function(x) return(x) )

    IBDList = lapply(dummyMarkerIndex,
                     function(thisMarkerIndex)
                     {
                       lapply(familyIndex,
                              function(thisFamilyIndex)
                              {
                                thisFamilyDummyHaps = dummyHaps[thisFamilyIndex]
                                thisFamilyDummyHaps1 = sapply(thisFamilyDummyHaps,function(x) return(x[thisMarkerIndex,1]))
                                thisFamilyDummyHaps2 = sapply(thisFamilyDummyHaps,function(x) return(x[thisMarkerIndex,2]))

                                ret = outer(thisFamilyDummyHaps1,thisFamilyDummyHaps1,"==")+
                                      outer(thisFamilyDummyHaps1,thisFamilyDummyHaps2,"==")+
                                      outer(thisFamilyDummyHaps2,thisFamilyDummyHaps1,"==")+
                                      outer(thisFamilyDummyHaps2,thisFamilyDummyHaps2,"==")
                                ret = ret/2
                              })
                     })

    return(new("strumData",
               dataVals = cbind(dataVals(sData), ret[order(pedData$sDataTempRowId),]),
               dataType = "Pedigree",
               phi = phi(sData),
               ibd = new("strumIBD",
                         ibdMarker = data.frame(Marker=dummyMarkerNames,
                                                Position=unlist(dummyMarkerInfo$markerPos),
                                                Chromosome=unlist(dummyMarkerInfo$Chrom)),
                         ibdMatrix = IBDList)))
  } else
  {
    if( sData@dataType == "Pedigree" )
      return(new("strumData",
                 dataVals = cbind(dataVals(sData), ret[order(pedData$sDataTempRowId),]),
                 dataType = "Pedigree",
                 phi = phi(sData)))
    else
      return(new("strumData", 
                 dataVals = cbind(dataVals(sData), ret),
                 dataType = "RawData"))
  }
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
.findDummyMarkerInfo = function(mInfo)
{
  retChrom            = list()
  retPos              = list()
  retMarkerNames      = list()
  dummyMarkerCountAll = 0

  for( chrom in unique(mInfo@markerFacts$chrom) )
  {
    chromRowPos      = which(markerFacts(mInfo)$chrom==chrom)
    chromMapPos      = markerFacts(mInfo)$mapPos[chromRowPos]
    minChromMapPos   = min(chromMapPos)
    maxChromMapPos   = max(chromMapPos)
    chromLength      = maxChromMapPos - minChromMapPos
    dummyMarkerCount = floor(chromLength/mInfo@intervalIBD)+1

    dummyMarkerCountAll = dummyMarkerCountAll+dummyMarkerCount

    thisPos = (0:(dummyMarkerCount-1))*intervalIBD(mInfo)+minChromMapPos

    retChrom[[as.character(chrom)]]       = rep(chrom, dummyMarkerCount)
    retPos[[as.character(chrom)]]         = thisPos
    retMarkerNames[[as.character(chrom)]] = paste("Chr",chrom,"Pos",thisPos,sep="")
  }

  return(list(Chrom               = retChrom,
              markerPos           = retPos,
              markerNames         = retMarkerNames,
              dummyMarkerCountAll = dummyMarkerCountAll))	
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
.computeIntervals = function(recombMapPos, markerMapPos)
{	
  numbRecombs = length(recombMapPos)
  numbMarkers = length(markerMapPos)

  recomVec = .Call("computeInterval", numbMarkers, recombMapPos, markerMapPos)

  dim(recomVec) = c(2, numbRecombs)
  if( recomVec[2,1] == 0 )
    recomVec[1,1] = 0 #little fix it

  return(recomVec)
}

#------------------------------------------------------------------------------
# Performs the genetic part of the simulation.
#  If offspring equals false,
#    a random mosaic of haplotypes from the marker model haplotypes is created.
#  If offspring is true,
#    a random mosaic from the parents haplotypes is created
#------------------------------------------------------------------------------
.simulateHaplotypes = function(
                       mInfo,
                       offspring = FALSE,
                       parentHaplotypes = character(),
                       parentDummyHaplotypes = character(),
                       dummyMarkerInfo
                     )
{
  numbMarkers = length(markerFacts(mInfo)$markerName)

  retMain  = numeric()
  retDummy = numeric()

  for( chrom in unique(markerFacts(mInfo)$chrom) )
  {
    chromRowPos    = which(markerFacts(mInfo)$chrom==chrom)
    chromMapPos    = markerFacts(mInfo)$mapPos[chromRowPos]
    minChromMapPos = min(chromMapPos)
    maxChromMapPos = max(chromMapPos)
    chromLength    = maxChromMapPos - minChromMapPos

    if( offspring )
      numbRecombs = rpois(1, 0.01 * chromLength)
    else
      numbRecombs = rpois(1, populationRecombRate(mInfo) * chromLength)
  }

  if( numbRecombs != 0 )
  {
    recombMapPos = minChromMapPos + runif(numbRecombs) * chromLength
    recombMapPos = sort(recombMapPos)

    # add extra recomb after end of the chromosome
    #  to get rid of "if" statements in the C function
    #-------------------------------------------------------
    recombMapPos    = c(recombMapPos, maxChromMapPos+1)
    recombIntervals = .computeIntervals(recombMapPos, markerFacts(mInfo)$mapPos)	

    if( offspring )
    {
      startHap = sample.int(n=2, 1)
      whichHap = 1+(1+(-1)^(startHap+1:ncol(recombIntervals)))/2

      recombIntervals = rbind(recombIntervals, whichHap)

      hapChrom = apply(recombIntervals, 2,
                       function(i) return(parentHaplotypes[i[1]:i[2],i[3]]))

      # only return this if IBD is supposed to be returned
      #  and if the individual is not a founder
      #----------------------------------------------------
      if( returnIBD(mInfo) )
      {
        whichHapDummy = 1+(1+(-1)^(startHap+1:ncol(recombIntervals)))/2
        recombIntervalsDummy = .computeIntervals(recombMapPos, dummyMarkerInfo$markerPos[[as.character(chrom)]])
        recombIntervalsDummy = rbind(recombIntervalsDummy, whichHapDummy)

        hapChromDummy = apply(recombIntervalsDummy, 2,
                              function(i) return(parentDummyHaplotypes[i[1]:i[2],i[3]]))
      }
    } else
    {
      recombIntervals = rbind(recombIntervals,
                              sample.int(n=ncol(haplotypes(mInfo)),
                              size = ncol(recombIntervals),
                              replace = TRUE))

      hapChrom = apply(recombIntervals, 2,
                       function(i) return(haplotypes(mInfo)[i[1]:i[2],i[3]]))
    }
  } else
  {
    if( offspring )
    {
      tmp = sample.int(n=2, 1)
      hapChrom = parentHaplotypes[,tmp]

      # only return this if IBD is supposed to be returned
      #  and if the individual is not a founder
      #----------------------------------------------------
      if( returnIBD(mInfo) )
        hapChromDummy = parentDummyHaplotypes[,tmp]

    } else
      hapChrom = haplotypes(mInfo)[,sample.int(n=ncol(haplotypes(mInfo)), size = 1, replace = TRUE)]
  }

  retMain = c(retMain, unlist(hapChrom))

  if( returnIBD(mInfo) && offspring )
    retDummy = c(retDummy, unlist(hapChromDummy))

  return(list(main=retMain, dummy=retDummy))			
}
