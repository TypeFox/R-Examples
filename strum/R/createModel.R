#==============================================================================
# File: createModel.R
#
# Author: Nathan Morris
#
# Notes: Functions to create a model in strum analysis.
#
# History: Initial implementation
#          Revision - yes Jan 2013
#==============================================================================

#------------------------------------------------------------------------------
# Create strumModel object
#------------------------------------------------------------------------------
createStrumModel = function(
                     formulas,
                     ascertainment = NULL,
                     defaultError = '<p,e>',
                     assumeExogCovariate = TRUE,
                     fixLoadingToOne = TRUE
                   )
{
  retObj = .createModel(formulas = formulas,
                        parseAsSimModel = FALSE,
                        defaultError = defaultError,
                        assumeExogCovariate = assumeExogCovariate,
                        fixLoadingToOne = fixLoadingToOne)

  retObj@ascertainment = ascertainment

  .printInfoLine("Creating strumModel", "Done", 52, 0)

  return(retObj)
}

#------------------------------------------------------------------------------
# Create simModel object
#------------------------------------------------------------------------------
createSimModel = function(
                   formulas,
                   tMissingRate = c(0),
                   markerInfo = NULL,
                   ascertainment = NULL,
                   defaultError = '<p,e>',
                   assumeExogCovariate = TRUE
                 )
{
  retObj = .createModel(formulas = formulas,
                        parseAsSimModel = TRUE,
                        defaultError = defaultError,
                        assumeExogCovariate = assumeExogCovariate,
                        fixLoadingToOne = FALSE)

  retObj@ascertainment = ascertainment
  retObj@markerInfo    = markerInfo

  if( any(tMissingRate < 0) )
  {
    warning(paste("\nThe missing rate of trait should be greater than or equal to 0. ",
                  "Setting to 0!\n", sep=""))

    tMissingRate[which(tMissingRate < 0)] = 0
  } else if( any(tMissingRate >= 1) )
  {
    warning(paste("\nAssuming the value of tMissingRate as ",
                  "missing proportion!\n", sep=""))

    tMissingRate[which(tMissingRate >= 1)] = tMissingRate[which(tMissingRate >= 1)]/100
  }

  varY = varList(retObj)$name[retObj@varList$inY]
  if( length(tMissingRate) == 1 )
  {
    mRate = rep(tMissingRate, length(varY))
  } else if( length(tMissingRate) != length(varY) )
  {
    warning(paste("\nThe number of missing rates specified is not the same ",
                  "as the number of observed trait variables. ",
                  "Setting to the first value for all trait variables!\n", sep=""))

    mRate = rep(tMissingRate[1], length(varY))    
  } else
  {
    mRate = tMissingRate
  }

  retObj@traitMissingRate = mRate

  .printInfoLine("Creating strumSimModel", "Done", 52, 0)

  return(retObj)
}

#------------------------------------------------------------------------------
# Create a strum model object - simModel or strumModel
#   if assumeExogCovariate = TRUE, then it is assumed that
#     any observed, quantitative exogenous traits are covariates.
#   Otherwise, they are treated as latent & dependent variables.
#------------------------------------------------------------------------------
.createModel = function(
                 formulas,
                 parseAsSimModel = FALSE,
                 defaultError = '<p,e>',
                 assumeExogCovariate = TRUE,
                 fixLoadingToOne = TRUE
               )
{
  if( class(formulas) != "character" | length(formulas) != 1 | !nzchar(formulas) )
    stop("formulas must be a character string.")

  formulasOrig = formulas

  re = .getRandomError(defaultError, paste("defaultError=",defaultError))
  if( length(re) == 0 )
    defaultError = '<e>'

  formulas = .parseFormulas(formulas, parseAsSimModel, defaultError, fixLoadingToOne)
  varLists = .buildVarList(formulas, parseAsSimModel, assumeExogCovariate)

  varList             = varLists$varList
  formulasMeasurement = varLists$formulasMeasurement
  formulasTrimmed     = formulas$formulasTrimmed
  formulasCausal      = formulas$formulasCausal
  formulasConstraint  = formulas$formulasConstraint
  ordinalVars         = formulas$ordinalVars
  indicatorVars       = formulas$indicatorVars

  eitaNames = varList$name[varList$inEita]
  yNames    = varList$name[varList$inY]
  covNames  = c("[intercept]",varList$name[varList$covariate])
  exoNames  = varList$name[varList$exogenous]
  exoNames  = exoNames[!is.na(exoNames)]
  
  randomEfectNames = unique(c(unique(unlist(lapply(formulasCausal, function(x) x$re))),
                              unique(unlist(lapply(formulasMeasurement, function(x) x$re)))))

  randomEfectNames = sort(randomEfectNames, decreasing = TRUE)

  #  create the function by first creating prototypes with NA and filling them in
  #-------------------------------------------------------------------------------

  # prototypes
  #------------
  ptLambda  = matrix(0,nrow=length(yNames),    ncol=length(eitaNames), dimnames=list(yNames,eitaNames))
  ptBeta    = matrix(0,nrow=length(eitaNames), ncol=length(eitaNames), dimnames=list(eitaNames,eitaNames))
  ptGammaS  = matrix(0,nrow=length(eitaNames), ncol=length(covNames),  dimnames=list(eitaNames,covNames))
  ptGammaM  = matrix(0,nrow=length(yNames),    ncol=length(covNames),  dimnames=list(yNames,covNames))
  ptE       = lapply(randomEfectNames,
                     function(x)
                     {
                       tmp = matrix(0,nrow=length(yNames),ncol=length(yNames), dimnames=list(yNames,yNames))
                       return(tmp)
                     })
  ptZ       = lapply(randomEfectNames,
                     function(x)
                     {
                       tmp = matrix(0,nrow=length(eitaNames),ncol=length(eitaNames), dimnames=list(eitaNames,eitaNames))
                       return(tmp)
                     })

  names(ptE) = randomEfectNames
  names(ptZ) = randomEfectNames

  defVal1 = function(p1,p2)
            {
              if( parseAsSimModel )
              {
                if( p1 == p2 )
                  return(1)
                else
                  return(0)
              } else
                  return(NA)
            }

  dfVC = strsplit(defaultError,split="[<]|[,]|[>]")[[1]]
  dfVC = dfVC[nchar(dfVC)!=0]
  exoPos = match(exoNames,eitaNames)

  for( k1 in exoPos )
    for( k2 in exoPos )
    {
      for( j in dfVC )
        ptZ[[j]][k1,k2] = defVal1(k1,k2)
    }

  # First add in information from the measurement and causal formulas
  #-------------------------------------------------------------------
  defVal2 = function()
            {
              if( parseAsSimModel )
                return(1)
              else
                return(NA)
            }

  # exogenous variables correlated by defualt
  #-------------------------------------------

  for( eq in formulasMeasurement )
  {
    lhsPos = match(eq$lhs, eitaNames)

    for( myi in 1:length(eq$rhsVar) )
    {
      myvar    = eq$rhsVar[myi]
      mycoef   = eq$rhsCoef[myi]
      myvarPos = match(myvar,yNames)
      ptLambda[myvarPos,lhsPos] = tryCatch(as.numeric(mycoef),
                                           warning=function(x) stop(paste("Loading value '",mycoef,"' most be numeric.",sep="")))
      for(j in eq$re)
        ptE[[j]][myvarPos,myvarPos] = defVal2()

      ptGammaM[myvarPos,"[intercept]"] = eq$intercept
    }
  }

  for( eq in formulasCausal )
  {
    isCovariateEq = eq$rhsVar %in% covNames

    if( eq$lhs %in% indicatorVars )
    {
      lhsPos = match(eq$lhs, yNames)

      if( length(eq$rhsVar) > 0 )
      {
        for( myi in 1:length(eq$rhsVar) )
        {
          myvar     = eq$rhsVar[myi]
          mycoef    = eq$rhsCoef[myi]
          myvarPos = match(myvar,yNames)

          if( isCovariateEq[myi] )
            ptGammaM[lhsPos,myvarPos] = tryCatch(as.numeric(mycoef),
                                                 warning=function(x) stop(paste("Coeficient '",mycoef,"' most be numeric.",sep="")))
          else
            stop("A model indicator variable is beeing adjusted by a non-covariate.")
        }
      }

      for( j in eq$re )
        ptE[[j]][lhsPos,lhsPos] = defVal2()

      ptGammaM[lhsPos,"[intercept]"] = eq$intercept

    } else
    {
      lhsPos = match(eq$lhs, eitaNames)
      if( length(eq$rhsVar) > 0 )
      {
        for( myi in 1:length(eq$rhsVar) )
        {
          myvar  = eq$rhsVar[myi]
          mycoef = eq$rhsCoef[myi]

          if( isCovariateEq[myi] )
          {
            myvarPos = match(myvar,colnames(ptGammaS))
            ptGammaS[lhsPos,myvarPos] = tryCatch(as.numeric(mycoef),
                                                 warning=function(x) stop(paste("Coeficient '",mycoef,"' most be numeric.",sep="")))
        
          } else
          {
            myvarPos = match(myvar,colnames(ptBeta))
            ptBeta[lhsPos,myvarPos] = tryCatch(as.numeric(mycoef),
                                             warning=function(x) stop(paste("Coeficient '",mycoef,"' most be numeric.",sep="")))
          }
        }
      }

      for( j in eq$re )
        ptZ[[j]][lhsPos,lhsPos] = defVal2()

      ptGammaS[lhsPos,"[intercept]"] = eq$intercept
    }
  }

  # Next add information from the value constraints (as opposed to equality constrains)
  #-------------------------------------------------------------------------------------
  valueConstraints           = sapply(formulasConstraint,function(eq) length(eq$rhs))==1
  ordinalPropertyConstraints = sapply(formulasConstraint,function(eq) eq$lhs[1])=="orinalProperties"

  for( s in formulasConstraint[valueConstraints] )
  {
    cmd  = s$lhs[1] 
    args = s$lhs[2:length(s$lhs)]

    if( cmd == "cov" )
    {
      if( length(args) < 2 )
        stop("Error in constraint equation: cov must include at least two arguments.")

      if( length(args) == 2 )
        vc = randomEfectNames
      else
      {
        vc = args[3:length(args)]
        tmp = vc[!(vc %in% randomEfectNames)]
        if( length(tmp) )
          stop(paste("Error in constraint equation: random effects '", tmp[1],
                     "' is not recognizable.", sep=""))
      }

      eitaMatch      = match(args[1:2],eitaNames)
      yMatch         = match(args[1:2],yNames)
      covariateMatch = match(args[1:2],covNames)

      if( any(!is.na(covariateMatch)) )
      {
        covariateMatch = covariateMatch[!is.na(covariateMatch)]
        tmp = ((args[1:2])[covariateMatch])[1]
        stop(paste("Error in constraint equation: covariance involving '", tmp,
                   "' cannot exist because it is a covariate. Consider changing the ",
                   "'assumeExogCovariate' option.", sep=""))
      } else if( all(!is.na(eitaMatch)) )
      {		
        for( thisj in vc )
        {
          ptZ[[thisj]][eitaMatch[1],eitaMatch[2]] = s$rhs
          ptZ[[thisj]][eitaMatch[2],eitaMatch[1]] = s$rhs
        }
      } else if( all(!is.na(yMatch)) )
      {
        for( j in vc )
        {
          ptE[[j]][yMatch[1],yMatch[2]] = s$rhs
          ptE[[j]][yMatch[2],yMatch[1]] = s$rhs
        }
      } else
        stop(paste("Error in constraint equation: covariance between '", args[1],
                   "' and '", args[2],"' is not allowed.",sep=""))
    } else if( cmd == "var" )
    {
      if( length(args) < 1 )
        stop("Error in constraint equation: var must include at least one argument.")

      if( length(args) == 1 )
        vc = randomEfectNames
      else
      {
        vc = args[2:length(args)]
        tmp = vc[!(vc %in% randomEfectNames)]
        if( length(tmp) )
          stop(paste("Error in constraint equation: random effects '", tmp[1],
                     "' is not recognizable.",sep=""))
      }

      eitaMatch      = match(args[1],eitaNames)
      yMatch         = match(args[1],yNames)
      covariateMatch = match(args[1],covNames)

      if( !is.na(covariateMatch) )
      {
        stop(paste("Error in constraint equation: variance involving '", args[1],
                   "' cannot exist because it is a covariate. Consider changing the ",
                   "'assumeExogCovariate' option.",sep=""))
      } else if( !is.na(eitaMatch) )
      {		
        for( j in vc )
          ptZ[[j]][eitaMatch,eitaMatch] = s$rhs
      } else if( !is.na(yMatch) )
      {
        for( j in vc )
          ptZ[[j]][yMatch,yMatch] = s$rhs
      } else
          stop(paste("Error in constraint equation: varaiance of '", args[1],
                     "' is not allowed because variable does not exist.",sep=""))
    } else if( cmd == "coef" )
    {
      if( length(args) != 2 )
        stop("Error in constraint equation: coef must have exactly two arguments.")

      eitaMatch      = match(args[1:2],eitaNames); isEitaMatch=!is.na(eitaMatch)
      yMatch         = match(args[1:2],yNames); isYMatch=!is.na(yMatch)
      covariateMatch = match(args[1:2],covNames); isCovariateMatch=!is.na(covariateMatch)
      
      if( all(isEitaMatch) )
      {
        ptBeta[eitaMatch[1],eitaMatch[2]] = s$rhs
      } else if( any(isYMatch) & any(isEitaMatch) )
      {
        tmp = .getIndex(yMatch,eitaMatch,isYMatch,isEitaMatch)
        ptLambda[tmp[1],tmp[2]] = s$rhs
      } else if( any(isYMatch) & any(isCovariateMatch) )
      {
        tmp = .getIndex(yMatch,covariateMatch,isYMatch,isCovariateMatch)
        ptGammaM[tmp[1],tmp[2]] = s$rhs
      } else if( any(isEitaMatch) & any(isCovariateMatch) )
      {
        tmp = .getIndex(eitaMatch,covariateMatch,isEitaMatch,isCovariateMatch)
        ptGammaS[tmp[1],tmp[2]] = s$rhs
      } else
        stop(paste("Error in constraint equation: coef between '", args[1],
                   "' and '", args[2],"' is not allowed.",sep=""))
    } else if( cmd == "ordinalProperties" )
    {
      stop("Error in constraint equation: ordinal properties not specified correctly.")
    } else
      stop(paste("Error in constraint equation: constraint term '", cmd,
                 "' is not recognizable.",sep=""))
  }

  buildFunction = function(protNum, measurement=FALSE)
  {
    paramCount = sum(is.na(protNum))
    if( parseAsSimModel )
      proto = paste(deparse(protNum),collapse ="")
    else
      proto = paste("matrix(",paste(deparse(as.numeric(protNum)),collapse=""),",",nrow(protNum),",",ncol(protNum),")",sep="")

    for( myi in 1:paramCount )
    {
      proto = sub("NA",paste("thB[",myi+curPos,"]",sep=""),proto)
      proto = sub("_real_","",proto)
    }

    proto = paste("function(thB) return(",proto,")")
    proto = parse(text=proto)
    curPos <<- curPos+paramCount #note incremented outside function

    if( measurement )
    {
      tmp = outer(rownames(protNum),colnames(protNum), function(x1,x2) paste(x2,x1,sep="=~"))
      tmp = tmp[is.na(protNum)]
      thBNames <<- c(thBNames,paste(tmp))
    } else
    {
      tmp = outer(rownames(protNum),colnames(protNum), function(x1,x2) paste(x1,x2,sep="~"))
      tmp = tmp[is.na(protNum)]
      thBNames <<- c(thBNames,paste(tmp))
    }

    return(eval(proto))	
  }

  buildVarFunction = function(protNum, vcname)
  {
    NACount = sum(is.na(protNum))
    if( parseAsSimModel )
      proto = paste(deparse(protNum),collapse ="")
    else
      proto = paste("matrix(",paste(deparse(as.numeric(protNum)),collapse=""),",",nrow(protNum),",",ncol(protNum),")",sep="")		

    tmp = matrix(0,nrow(protNum),nrow(protNum))
    tmp[is.na(protNum)] = NA
    paramCount = 1
    myindex = rep(0,NACount)
    keepName = rep(TRUE,NACount)
    myna = 1
    mynames = outer(rownames(protNum),colnames(protNum), function(x1,x2) paste(x1,x2,sep="~~"))

    for( i in 1:nrow(protNum) )
      for( j in 1:nrow(protNum) )
      {
        if( is.na(tmp[j,i]) )
        {
          if( i > j )
          {
            myindex[myna] = tmp[i,j]
            keepName[myna] = FALSE

            if( mynames[i,j] > mynames[j,i] )
              mynames[i,j] = mynames[j,i]
          } else
          {
            tmp[j,i] = paramCount
            myindex[myna] = paramCount
            paramCount = paramCount+1
          }

          myna = myna+1
        }
      }

    paramCount = paramCount-1 #adjust for the one overcount

    for( myi in 1:NACount )
    {
      proto = sub("NA",paste("thB[",myindex[myi]+curPos,"]",sep=""),proto)
      proto = sub("_real_","",proto)
    }

    proto = paste("function(thB) return(",proto,")")
    proto = parse(text=proto)
    curPos <<- curPos+paramCount #note incremented outside function

    mynames = mynames[is.na(protNum)]
    mynames = mynames[keepName]
    if( paramCount > 0 )
      thBNames <<- c(thBNames,paste(mynames,"<",vcname,">",sep=""))#note incremented outside function

    return(eval(proto))	
  }

  curPos = 0
  thBNames = character()

  Lambda = buildFunction(ptLambda,T)
  Beta   = buildFunction(ptBeta)
  GammaS = buildFunction(ptGammaS)
  GammaM = buildFunction(ptGammaM)
  E      = mapply(buildVarFunction,ptE, names(ptE))
  Z      = mapply(buildVarFunction,ptZ, names(ptZ))

  thToThB = .buildThToThB(formulasConstraint, valueConstraints,
                          ordinalPropertyConstraints, randomEfectNames, thBNames)

  if( parseAsSimModel )
    retObj = new("strumSimModel",
                 varList          = varList,
                 formulas         = formulasTrimmed,
                 allRandomEffects = randomEfectNames,
                 paramNames       = thBNames,
                 ascertainment    = "ANY",
                 L                = Lambda,
                 B                = Beta,
                 E                = E,
                 Z                = Z,
                 Gs               = GammaS,
                 Gm               = GammaM,
                 thToThB          = thToThB)
  else
    retObj = new("strumModel",
                 varList          = varList,
                 formulas         = formulasTrimmed,
                 allRandomEffects = randomEfectNames,
                 paramNames       = thBNames,
                 ascertainment    = "ANY",
                 L                = Lambda,
                 B                = Beta,
                 E                = E,
                 Z                = Z,
                 Gs               = GammaS,
                 Gm               = GammaM,
                 thToThB          = thToThB)

  # return it
  #-----------
  return(retObj)
}

#------------------------------------------------------------------------------
# Parse out formulas
#------------------------------------------------------------------------------
.parseFormulas = function(formulas, parseAsSimModel, defaultError, fixLoadingToOne)
{
  # remove end of line
  #--------------------
  formulas = strsplit(formulas, split="\n")[[1]]

  formulas0 = lapply(formulas, function(s) return(.trimSpace(s)))
  formulas0 = formulas0[nzchar(formulas0)]
  formulas0 = lapply(formulas0, function(s) return(paste("  ", s,"\n",sep="")))
  formulas0 = paste(formulas0, collapse = '')

  # remove comments
  #-----------------
  formulas = sapply(strsplit(formulas,split="#"), function(x) return(.trimSpace(x[1])))

  # remove empty lines
  #--------------------
  formulas = formulas[nzchar(formulas)]
  formulas = formulas[!is.na(formulas)]

  if( length(formulas) == 0 )
    stop("No valid eqautions exist in the formulas!")

  # find the measurement formulas first
  #-------------------------------------
  formulasTmp         = strsplit(formulas, split="=~")
  formulasMeasurement = formulasTmp[sapply(formulasTmp, function(x) return(length(x)==2))]
  formulasRest        = unlist(formulasTmp[sapply(formulasTmp, function(x) return(length(x)==1))])

  if( (length(formulasMeasurement)+length(formulasRest)) != length(formulasTmp) ) 
    stop("Formulas can contain only one '=~' per line!")

  formulasCausal     = list()
  formulasConstraint = list()

  if( length(formulasRest) != 0 )
  {
    # find the causal formulas
    #--------------------------
    formulasTmp        = strsplit(formulasRest, split="~")
    formulasCausal     = formulasTmp[sapply(formulasTmp, function(x) return(length(x)==2))]
    formulasConstraint = unlist(formulasTmp[sapply(formulasTmp, function(x) return(length(x)==1))])

    if( (length(formulasCausal)+length(formulasConstraint)) != length(formulasTmp) ) 
      stop("Formulas can contain only one '~' per line!")
  
    if( length(formulasConstraint) != 0 )
    {
      # find the constraint equations
      #-------------------------------
      formulasConstraint = strsplit(formulasConstraint, split="=")
      tmp = sapply(formulasConstraint, function(x) return(length(x)==1))
      if( any(tmp) )
        stop("All non empty lines in the formula must contain '=~' or '~' or '='!" )

      tmp = sapply(formulasConstraint, function(x) return(length(x)==3))
      if( any(tmp) )
        stop("Formulas can contain only one '=' per line!") #TODO: expand instead of error out
    }
  }

  # parse out the constraint equations
  #------------------------------------
  formulasConstraint = lapply(formulasConstraint,
                              function(x)
                              {
                                x = strsplit(x , split="[(]|[)]|[,]")
                                x = lapply(x, function(y)
                                              {
                                                y = sapply(y, .trimSpace)
                                                y = y[nchar(y)!=0 & !is.na(y)]
                                                return(as.vector(y))
                                              })

                                names(x) = c("lhs","rhs")
                                if( length(x$rhs) == 1 )
                                {
                                  if( x$rhs == "NA" )
                                  {
                                    if( parseAsSimModel )
                                      stop("Error in constraint equation: constriant equations cannot be set to 'NA' for the purpose of simulation.")
                                    x$rhs = NA
                                  } else
                                    x$rhs = tryCatch(as.numeric(x$rhs),
                                                     warning=function(x) stop(paste("Constraint value",x$rhs,"most be numeric.")))
                                }

                                return(x)
                              })

  ordinalVars = sapply(formulasConstraint,
                       function(eq)
                       {
                         if(eq$lhs[1] %in% c("orinalProperties","op"))
                           return(eq$lhs[2])
                         else
                           return(NA)
                       })

  ordinalVars = ordinalVars[!is.na(ordinalVars)]
  
  # parse out the random effects and terms
  #----------------------------------------
  formulasCausal = lapply(formulasCausal,
                          .getFormulaStuff,
                          parseAsSimModel  = parseAsSimModel,
                          defualt1         = FALSE,
                          ordinalVars      = ordinalVars,
                          defualtIntercept = FALSE,
                          defaultError     = defaultError)

  formulasMeasurement = lapply(formulasMeasurement,
                               .getFormulaStuff,
                               parseAsSimModel  = parseAsSimModel,
                               defualt1         = fixLoadingToOne,
                               ordinalVars      = ordinalVars,
                               defualtIntercept = TRUE,
                               defaultError     = defaultError)

  names(formulasCausal)      = sapply(formulasCausal, function(x) x$lhs)
  names(formulasMeasurement) = sapply(formulasMeasurement, function(x) x$lhs)

  indicatorVars  = unique(unlist(lapply(formulasMeasurement, function(x) x$rhsVar)))

  return(list(formulasTrimmed     = formulas0,
              formulasMeasurement = formulasMeasurement,
              formulasCausal      = formulasCausal,
              formulasConstraint  = formulasConstraint,
              ordinalVars         = ordinalVars,
              indicatorVars       = indicatorVars))
}

#------------------------------------------------------------------------------
# Fill out varList
#   If assumeExogCovariate = TRUE, then it is assumed that
#     any observed, quantitative exogenous traits are covariates.
#   Otherwise, they are treated as latent & dependent variables.
#------------------------------------------------------------------------------
.buildVarList = function(formulas, parseAsSimModel, assumeExogCovariate)
{
  formulasConstraint  = formulas$formulasConstraint
  formulasCausal      = formulas$formulasCausal
  formulasMeasurement = formulas$formulasMeasurement
  ordinalVars         = formulas$ordinalVars
  indicatorVars       = formulas$indicatorVars

  latentVars     = c()
  rhsVars        = c()
  endogenousVars = c()
  
  if( length(formulasMeasurement) )
    latentVars = as.vector(sapply(formulasMeasurement, function(x) x$lhs))
  
  if( length(formulasCausal) )
  {
    rhsVars        = unique(unlist(lapply(formulasCausal, function(x) x$rhsVar)))
    endogenousVars = as.vector(sapply(formulasCausal, function(x) x$lhs))
  }

  allVars = unique(c(latentVars, endogenousVars, rhsVars, indicatorVars))
  varList = data.frame(name=allVars, stringsAsFactors=FALSE)

  #hmm... so what about variables that dont appear in the measurement model?
  #-------------------------------------------------------------------------
  varList$obs = !(varList$name %in% latentVars)
  isIndicator = varList$name %in% indicatorVars
  isEndog     = varList$name %in% endogenousVars
  isOrdinal   = varList$name %in% ordinalVars

  strangeVarsOrdinal = varList$name[!isIndicator & varList$obs &  isOrdinal]
  strangeVarsEndog   = varList$name[!isIndicator & varList$obs &  isEndog & !isOrdinal] #These may be added to the measurement model with no error
  strangeVarsExog    = varList$name[!isIndicator & varList$obs & !isEndog & !isOrdinal] #These may be treated as covariates
 
  #assume these variables are observed
  #------------------------------------
  for( vname in strangeVarsEndog )
  {
    #add measurement equation with a loading of 1
    if( !formulasCausal[[vname]]$interceptSet )
      formulasCausal[[vname]]$intercept=NA

    intercept = 0
    if( !parseAsSimModel )
      intercept = as.numeric(NA)

    formulasMeasurement[[vname]] = list(lhs = vname,
                                        rhsVars = vname,
                                        rhsCoef = "1",
                                        re = character(),
                                        intercept = intercept, #0)
                                        interceptSet = FALSE)
  }

  for( vname in strangeVarsOrdinal )
  {
    #add measurement equation with a loading of 1
    formulasMeasurement[[vname]] = list(lhs = vname,
                                        rhsVars = vname,
                                        rhsCoef = "1",
                                        re = character(),
                                        intercept = 0)
  }

  if( !assumeExogCovariate )
  {
    for( vname in strangeVarsExog )
    {
      #add measurement equation with a loading of 1
      intercept = 0
      if( !parseAsSimModel )
        intercept = as.numeric(NA)

      formulasMeasurement[[vname]] = list(lhs = vname, 
                                          rhsVars = vname,
                                          rhsCoef = "1",
                                          re = character(),
                                          intercept = intercept, #0)
                                          interceptSet = FALSE)
    }

    varList$covariate = FALSE

  } else
  {
    if( length(strangeVarsOrdinal) > 0 )
      stop(paste("Latent class analysis not implemented. Variable",
                 strangeVarsOrdinal[1],
                 "is defined as latent but also ordinal."))

    varList$covariate = varList$name %in% strangeVarsExog
  }
  
  eitaNames      = lapply(formulasMeasurement, function(x) x$lhs)
  varList$inEita = varList$name %in% eitaNames
  varList$inY    = varList$obs & !varList$covariate
  
  # variables that are in eita but never appear on the lhs are exogenous
  #----------------------------------------------------------------------
  varList$exogenous = !(varList$name %in% endogenousVars)
  varList$exogenous[!varList$inEita] = NA

  return(list(varList = varList,
              formulasMeasurement = formulasMeasurement))
}

#------------------------------------------------------------------------------
# Get the components from a formula
#------------------------------------------------------------------------------
.getFormulaStuff = function(x, parseAsSimModel, defualt1,
                            ordinalVars, defualtIntercept, defaultError)
{
  # parse out the error terms first
  #---------------------------------
  lhs = .trimSpace(x[1])
  rhstmp = .trimSpace(x[2])
  search1 = regexpr("[+]([[:space:]]*)<(.*)>",rhstmp) #search for + space < something >

  if( search1 > 0 )
  {
    rhs = substr(rhstmp,1,search1-1)
    re  = substr(rhstmp,search1+1,search1+attr(search1,"match.length"))
  } else
  {
    search2 = regexpr("<(.*)>",rhstmp)
    if( search2 > 0 )
    {
      rhs = substr(rhstmp,1,search2-1)
      re  = substr(rhstmp,search2,search2+attr(search2,"match.length"))
    } else
    {
      rhs = rhstmp
      re = defaultError
    }
  }

  re = .getRandomError(re, rhstmp)

  if( length(re) == 0 )
    re = .getRandomError(defaultError, paste("defaultError=",defaultError))

  # parse out other terms next
  #----------------------------
  rhs = strsplit(rhs, split="+", fixed=TRUE)[[1]]

  rhs = lapply(rhs,
               function(term)
               {
                 myterm = strsplit(term, split="*", fixed=TRUE)[[1]]

                 if( length(myterm) == 1 )
                   myterm = c(NA,.doSimpleCheck(myterm[1]))
                 else if( length(myterm) == 2 )
                   myterm = c(.doSimpleCheck(myterm[1]),.doSimpleCheck(myterm[2]))
                 else
                   stop(paste("Error in expression:",term))
               })

  rhsCoef = as.vector(sapply(rhs,
                             function(x)
                             {
                               mycoef = x[1]
                               if( !is.na(mycoef) )
                                 mycoef = tryCatch(as.numeric(mycoef),
                                                   warning = function(x) stop(paste("Cefficient value '",mycoef,"' most be numeric.", sep="")))
                               else if(parseAsSimModel)
                                 mycoef = 1
                               else
                                 mycoef = as.numeric(NA)
                               return(mycoef)
                             }))

  rhsVar = sapply(rhs, function(x) return(x[2]))
  suppressWarnings(varIsNumeric<-!is.na(as.numeric(rhsVar)))

  isIntercept = varIsNumeric | rhsVar=="NA"

  if( parseAsSimModel )
    intercept = 0
  else if(lhs %in% ordinalVars)
    intercept = 0
  else if(defualtIntercept)
    intercept = as.numeric(NA)
  else
    intercept = 0

  interceptSet = FALSE
  if( sum(isIntercept) == 1 )
  {
    intercept = as.numeric(rhsVar[isIntercept])
    if( lhs %in% ordinalVars )
      if( intercept != 0 )
        stop("Ordinal variables must have an intercept of 0.")

    rhsVar = rhsVar[!isIntercept]
    rhsCoef = rhsCoef[!isIntercept]
    interceptSet = TRUE #this is a dumb workaround to allow more defualt stuff to be implemented later
  } else if( sum(isIntercept) > 1 )
    stop("Formula found which defines the intercept twice.")

          
  if( all(is.na(rhsCoef)) & defualt1 )
    rhsCoef[1] = 1

  if( length(rhsCoef) > 0 )
    names(rhsCoef) = paste(lhs,rhsVar,sep="~")  

  return(list(lhs          = lhs,
              rhsCoef      = rhsCoef,
              rhsVar       = rhsVar,
              intercept    = intercept,
              interceptSet = interceptSet,
              re           = re))
}

#------------------------------------------------------------------------------
# Function for some simple parsing and checking
#------------------------------------------------------------------------------
.doSimpleCheck = function(z)
{
  tmp = .trimSpace(z)

  if( !grepl('[[\\t\\r\\n\\v\\f\\[\\^!#$%&()*+,/:;<=>?@`{|}~-]]', tmp) )
    return(tmp)
  else
    stop(paste("Error in model expression:",tmp))
}

#------------------------------------------------------------------------------
# Parse random errors
#------------------------------------------------------------------------------
.getRandomError = function(re, formula)
{
  re = strsplit(re,split="<")[[1]]

  if( length(re) != 2 )
    stop(paste("Error in the statement:",formula))

  re = strsplit(re[2],split=">")[[1]][1]
  re = strsplit(re,split=",|[+]")[[1]]

  return(re)
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
.getIndex = function(i1, i2, found1, found2)
{
  if( found1[1] & found2[2] )
    return(c(i1[1], i2[2]))
  else if( found1[2] & found2[1] )
    return(c(i1[2], i2[1]))
  else 
    stop("Error in constrain of a coeficient.")
}

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
.buildThToThB = function(formulasConstraint, valueConstraints,
                         ordinalPropertyConstraints, randomEfectNames, thBNames)
{
  equalityClass = list()

  if( sum(!valueConstraints & !ordinalPropertyConstraints) > 0 )
  {
    for( s in formulasConstraint[!valueConstraints & !ordinalPropertyConstraints] )
    {
      myrhs = .processEqualConstrSide(s$rhs, randomEfectNames)
      mylhs = .processEqualConstrSide(s$lhs, randomEfectNames)

      if( length(myrhs) != length(mylhs) )
        stop("Error in constrain: number of variance components being constrained is not the same.")

      for( myi in 1:length(myrhs) )
      {
        m1 = match(c(myrhs[myi],mylhs[myi]), thBNames)

        if( any(is.na(m1)) )
          stop(paste("Error in constraint: could not find parameter '",
                     ((c(myrhs,mylhs))[is.na(m1)])[1],
                     "'. You may need to free the parameters first,",
                     " and then constrain them to be equal.",sep=""))

        if( length(equalityClass) > 0 )
        {
          grp1 = which(sapply(equalityClass, function(cl) return(m1[1] %in% cl)))
          grp2 = which(sapply(equalityClass, function(cl) return(m1[1] %in% cl)))
        } else
        {
          grp1 = logical()
          grp2 = logical()
        }

        if( length(grp1)==1 & length(grp2)==1 )
        {
          #merge the two groups
          equalityClass[[grp1]] = c(equalityClass[[grp1]],equalityClass[[grp2]])
          equalityClass = equalityClass[-grp2]
        } else if( length(grp1)==0 & length(grp2)==0 )
        {
          #new group
          equalityClass[[length(equalityClass)+1]] = m1
        } else if( length(grp1)==1 & length(grp2)==0 )
        {
          #add element to grp
          equalityClass[[grp1]] = c(equalityClass[[grp1]],m1[2])
        } else if( length(grp1)==0 & length(grp2)==1 )
        {
          #add element to grp
          equalityClass[[grp2]] = c(equalityClass[[grp2]],m1[1])
        } else
          stop("Error in equality constraint.")
      }
    }

    equalityClass = lapply(equalityClass, sort)

    tmp = rep(0, length(thBNames))

    for( cl in equalityClass )
      for( pos in cl[-1] )
        tmp[pos] = cl[1]

    mycount = 1
    for( i in 1:length(tmp) )
    {
      if( tmp[i] == 0 )
      {
        tmp[i] = mycount
        mycount = mycount+1
      } else
        tmp[i]=tmp[tmp[i]]
    }

    thToThB = tmp
  } else
    thToThB = 1:length(thBNames)

  return(thToThB)
}

#------------------------------------------------------------------------------
#  apply the equality constraints 
#------------------------------------------------------------------------------
.processEqualConstrSide = function(mys, randomEfectNames)
{
  cmd = mys[1] 
  args = mys[2:length(mys)]

  if( cmd == "cov" )
  {
    if( length(args) < 2 )
      stop("Error in constraint equation: cov must include at least two arguments.")
    if( length(args) == 2 )
      vc = randomEfectNames
    else
    {
      vc = args[3:length(args)]
      tmp = vc[!(vc %in% randomEfectNames)]
      if( length(tmp) )
        stop(paste("Error in constraint equation: random effects '", tmp[1],
                   "' is not recognizable.",sep=""))
    }

    retVal = paste(args[1],"~~",args[2],"<",vc,">",sep="")
  } else if( cmd == "var" )
  {
    if( length(args) < 1 )
      stop("Error in constraint equation: var must include at least one argument.")
    if( length(args) == 1 )
      vc = randomEfectNames
    else
    {
      vc = args[2:length(args)]
      tmp = vc[!(vc %in% randomEfectNames)]
      if( length(tmp) )
        stop(paste("Error in constraint equation: random effects '", tmp[1],
                   "' is not recognizable.",sep=""))
    }

    retVal = paste(args[1],"~~",args[1],"<",vc,">",sep="")
  } else if( cmd == "coef" )
  {
    if( length(args) != 2 )
      stop("Error in constraint equation: coef must have exactly two arguments.")

    retVal = paste(args[1],"~",args[2],sep="")
  } else if( cmd == "ordinalProperties" )
  {
    stop(paste("Error in constraint equation: ordinalProperties not specified correctly")) 
  } else
    stop(paste("Error in constraint equation: constraint term '", cmd,
               "' is not recognizable.",sep=""))
  return(retVal)
}
