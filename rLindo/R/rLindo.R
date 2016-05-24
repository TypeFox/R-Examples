#rLindo.R
#The R interface to LINDO API 8.0.
#This file includes all R interface wrapper functions for LINDO API C functions.
#Copyright (C) 2013 LINDO Systems.


#************************************************************##
# Structure Creation and Deletion Routines (5)                #
#*************************************************************#
rLScreateEnv <- function(){

    ans <- .Call("rcLScreateEnv", PACKAGE = "rLindo")
    return(ans)
}

rLScreateModel <- function(env){

    ans <- .Call("rcLScreateModel", PACKAGE = "rLindo", env)
    return(ans)
}

rLSdeleteEnv <- function(env){

    ans <- .Call("rcLSdeleteEnv", PACKAGE = "rLindo", env)
    return(ans)
}

rLSdeleteModel <- function(model){

    ans <- .Call("rcLSdeleteModel", PACKAGE = "rLindo", model)
    return(ans)
}

rLScopyParam <- function(smodel, tmodel, nSolverType){

    ans <- .Call("rcLScopyParam", PACKAGE = "rLindo", smodel,tmodel,as.integer(nSolverType))
    return(ans)
}

#*************************************************************#
# Model I-O Routines (18)                                     #
#*************************************************************#
rLSreadMPSFile <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSreadMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSwriteMPSFile <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSwriteMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSreadLINDOFile <- function(model, pszFname){

    ans <- .Call("rcLSreadLINDOFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteLINDOFile <- function(model, pszFname){

    ans <- .Call("rcLSwriteLINDOFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSreadLINDOStream <- function(model, pszStream, nStreamLen){

    ans <- .Call("rcLSreadLINDOStream", PACKAGE = "rLindo", 
                 model,
                 as.character(pszStream),
                 as.integer(nStreamLen))
    return(ans)
}

rLSwriteLINGOFile <- function(model, pszFname){

    ans <- .Call("rcLSwriteLINGOFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteDualMPSFile <- function(model, pszFname, nFormat, nObjSense){

    ans <- .Call("rcLSwriteDualMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat),
                 as.integer(nObjSense))
    return(ans)
}

rLSwriteSolution <- function(model, pszFname){

    ans <- .Call("rcLSwriteSolution", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteSolutionOfType <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSwriteSolutionOfType", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSwriteIIS <- function(model, pszFname){

    ans <- .Call("rcLSwriteIIS", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteIUS <- function(model, pszFname){

    ans <- .Call("rcLSwriteIUS", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSreadMPIFile <- function(model, pszFname){

    ans <- .Call("rcLSreadMPIFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteMPIFile <- function(model, pszFname){

    ans <- .Call("rcLSwriteMPIFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSwriteWithSetsAndSC <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSwriteWithSetsAndSC", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSreadBasis <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSreadBasis", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSwriteBasis <- function(model, pszFname, nFormat){

    ans <- .Call("rcLSwriteBasis", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSreadLPFile <- function(model, pszFname){

    ans <- .Call("rcLSreadLPFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSreadLPStream <- function(model, pszStream, nStreamLen){

    ans <- .Call("rcLSreadLPStream", PACKAGE = "rLindo", 
                 model,
                 as.character(pszStream),
                 as.integer(nStreamLen))
    return(ans)
}

rLSsetPrintLogNull <- function(model){

    ans <- .Call("rcLSsetPrintLogNull", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

#*************************************************************#
# Error Handling Routines (3)                                 #
#*************************************************************#
rLSgetErrorMessage <- function(env, errorcode){

    ans <- .Call("rcLSgetErrorMessage", PACKAGE = "rLindo", 
                 env,
                 as.integer(errorcode))
    return(ans)
}

rLSgetFileError <- function(model){

    ans <- .Call("rcLSgetFileError", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetErrorRowIndex <- function(model){

    ans <- .Call("rcLSgetErrorRowIndex", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

#*************************************************************#
# Routines for Setting and Retrieving Parameter Values (17)   #
#*************************************************************#
rLSsetModelDouParameter <- function(model,nParameter,dValue){

    ans <- .Call("rcLSsetModelDouParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter),
                 as.numeric(dValue))
    return(ans)
}

rLSgetModelDouParameter <- function(model,nParameter){

    ans <- .Call("rcLSgetModelDouParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter))
    return(ans)
}

rLSsetModelIntParameter <- function(model,nParameter,nValue){

    ans <- .Call("rcLSsetModelIntParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter),
                 as.integer(nValue))
    return(ans)
}

rLSgetModelIntParameter <- function(model,nParameter){

    ans <- .Call("rcLSgetModelIntParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter))
    return(ans)
}

rLSsetEnvDouParameter <- function(env,nParameter,dValue){

    ans <- .Call("rcLSsetEnvDouParameter", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParameter),
                 as.numeric(dValue))
    return(ans)
}

rLSgetEnvDouParameter <- function(env,nParameter){

    ans <- .Call("rcLSgetEnvDouParameter", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParameter))
    return(ans)
}

rLSsetEnvIntParameter <- function(env,nParameter,nValue){

    ans <- .Call("rcLSsetEnvIntParameter", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParameter),
                 as.integer(nValue))
    return(ans)
}

rLSgetEnvIntParameter <- function(env,nParameter){

    ans <- .Call("rcLSgetEnvIntParameter", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParameter))
    return(ans)
}

rLSreadModelParameter <- function(model,pszFname){

    ans <- .Call("rcLSreadModelParameter", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSreadEnvParameter <- function(env,pszFname){

    ans <- .Call("rcLSreadEnvParameter", PACKAGE = "rLindo", 
                 env,
                 as.character(pszFname))
    return(ans)
}

rLSwriteModelParameter <- function(model,pszFname){

    ans <- .Call("rcLSwriteModelParameter", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))
    return(ans)
}

rLSgetIntParameterRange <- function(model,nParameter){

    ans <- .Call("rcLSgetIntParameterRange", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter))
    return(ans)
}

rLSgetDouParameterRange <- function(model,nParameter){

    ans <- .Call("rcLSgetDouParameterRange", PACKAGE = "rLindo", 
                 model,
                 as.integer(nParameter))
    return(ans)
}

rLSgetParamShortDesc <- function(env,nParam){

    ans <- .Call("rcLSgetParamShortDesc", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParam))
    return(ans)
}

rLSgetParamLongDesc <- function(env,nParam){

    ans <- .Call("rcLSgetParamLongDesc", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParam))
    return(ans)
}

rLSgetParamMacroName <- function(env,nParam){

    ans <- .Call("rcLSgetParamMacroName", PACKAGE = "rLindo", 
                 env,
                 as.integer(nParam))
    return(ans)
}

rLSgetParamMacroID <- function(env,szParam){

    ans <- .Call("rcLSgetParamMacroID", PACKAGE = "rLindo", 
                 env,
                 as.character(szParam))
    return(ans)
}
    
#*************************************************************#
# Model Loading Routines (19)                                 #
#*************************************************************#
rLSloadLPData <- function(model,
                          nCons,
                          nVars,
                          nObjSense,
                          dObjConst,
                          padC,
                          padB,
                          pszConTypes,
                          nAnnz,
                          paiAcols,
                          panAcols = NULL,
                          padAcoef,
                          paiArows,
                          padL = NULL,
                          padU = NULL)
{

    if (is.null(panAcols)) 
    {
        rpanAcols <- as.null(panAcols)
    }
    else 
    {
        rpanAcols <- as.integer(panAcols)
    }
    if (is.null(padL)) 
    {
        rpadL<- as.null(padL)
    }
    else 
    {
        rpadL <- as.numeric(padL)
    }
    if (is.null(padU)) 
    {
        rpadU <- as.null(padU)
    }
    else 
    {
        rpadU <- as.numeric(padU)
    }
    ans <- .Call("rcLSloadLPData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(nVars),
                 as.integer(nObjSense),
                 as.numeric(dObjConst),
                 as.numeric(padC),
                 as.numeric(padB),
                 as.character(pszConTypes),
                 as.integer(nAnnz),
                 as.integer(paiAcols),
                 rpanAcols,
                 as.numeric(padAcoef),
                 as.integer(paiArows),
                 rpadL,
                 rpadU)
    return(ans)
}

rLSloadQCData <- function(model,
                          nQCnnz,
                          paiQCrows,
                          paiQCcols1,
                          paiQCcols2,
                          padQCcoef)                          
{

    ans <- .Call("rcLSloadQCData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQCnnz),
                 as.integer(paiQCrows),
                 as.integer(paiQCcols1),
                 as.integer(paiQCcols2),
                 as.numeric(padQCcoef))         
    return(ans)
}

rLSloadConeData <- function(model,
                            nCone,
                            pszConeTypes,
                            paiConebegcone,
                            paiConecols)                          
{

    ans <- .Call("rcLSloadConeData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCone),
                 as.character(pszConeTypes),
                 as.integer(paiConebegcone),
                 as.integer(paiConecols))         
    return(ans)
}

rLSloadSETSData <- function(model,
                            nSETS,
                            pszSETStype,
                            paiCARDnum,
                            paiSETSbegcol,
                            paiSETScols)                            
{

    ans <- .Call("rcLSloadSETSData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSETS),
                 as.character(pszSETStype),
                 as.integer(paiCARDnum),    
                 as.integer(paiSETSbegcol),
                 as.integer(paiSETScols))         
    return(ans)
}

rLSloadSemiContData <- function(model,
                                nSCVars,
                                paiVars,
                                padL,
                                padU)                          
{

    ans <- .Call("rcLSloadSemiContData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSCVars),
                 as.integer(paiVars),
                 as.numeric(padL),
                 as.numeric(padU))         
    return(ans)
}

rLSloadVarType <- function(model,
                           spszVarTypes)                          
{

    ans <- .Call("rcLSloadVarType", PACKAGE = "rLindo", 
                 model,
                 as.character(spszVarTypes))         
    return(ans)
}

rLSloadNameData <- function(model,
                            pszTitle = NULL,
                            pszObjName = NULL,
                            pszRhsName = NULL,
                            pszRngName = NULL,
                            pszBndname = NULL,
                            paszConNames = NULL,
                            paszVarNames = NULL,
                            paszConeNames = NULL)                          
{
    if (is.null(pszTitle)) 
    {
        rpszTitle <- as.null(pszTitle)
    }
    else 
    {
        rpszTitle <- as.character(pszTitle)
    }
    if (is.null(pszObjName)) 
    {
        rpszObjName <- as.null(pszObjName)
    }
    else 
    {
        rpszObjName <- as.character(pszObjName)
    }
    if (is.null(pszRhsName)) 
    {
        rpszRhsName <- as.null(pszRhsName)
    }
    else 
    {
        rpszRhsName <- as.character(pszRhsName)
    }
    if (is.null(pszRngName)) 
    {
        rpszRngName <- as.null(pszRngName)
    }
    else 
    {
        rpszRngName <- as.character(pszRngName)
    }
    if (is.null(pszBndname)) 
    {
        rpszBndname <- as.null(pszBndname)
    }
    else 
    {
        rpszBndname <- as.character(pszBndname)
    }
    if (is.null(paszConNames)) 
    {
        rpaszConNames <- as.null(paszConNames)
    }
    else 
    {
        rpaszConNames <- as.character(paszConNames)
    }
    if (is.null(paszVarNames)) 
    {
        rpaszVarNames <- as.null(paszVarNames)
    }
    else 
    {
        rpaszVarNames <- as.character(paszVarNames)
    }
    if (is.null(paszConeNames)) 
    {
        rpaszConeNames <- as.null(paszConeNames)
    }
    else 
    {
        rpaszConeNames <- as.character(paszConeNames)
    }
    ans <- .Call("rcLSloadNameData", PACKAGE = "rLindo", 
                 model,
                 rpszTitle,
                 rpszObjName,
                 rpszRhsName,
                 rpszRngName,
                 rpszBndname,
                 rpaszConNames,
                 rpaszVarNames,
                 rpaszConeNames)         
    return(ans)
}

rLSloadNLPData <- function(model,
                           paiNLPcols,
                           panNLPcols,
                           padNLPcoef = NULL,
                           paiNLProws,
                           nNLPobj,
                           paiNLPobj,
                           padNLPobj = NULL)                          
{
    if (is.null(padNLPcoef)) 
    {
        rpadNLPcoef <- as.null(padNLPcoef)
    }
    else 
    {
        rpadNLPcoef <- as.numeric(padNLPcoef)
    }
    if (is.null(padNLPobj)) 
    {
        rpadNLPobj <- as.null(padNLPobj)
    }
    else 
    {
        rpadNLPobj <- as.character(padNLPobj)
    } 

    ans <- .Call("rcLSloadNLPData", PACKAGE = "rLindo", 
                 model,
                 as.integer(paiNLPcols),
                 as.integer(panNLPcols),
                 rpadNLPcoef,
                 as.integer(paiNLProws),
                 as.integer(nNLPobj),
                 as.integer(paiNLPobj),
                 rpadNLPobj)         
    return(ans)
}

rLSloadInstruct <- function(model,
                            nCons,
                            nObjs,
                            nVars,
                            nNumbers,
                            panObjSense,
                            pszConType,
                            pszVarType = NULL,
                            panInstruct,
                            nInstruct,
                            paiVars = NULL,
                            padNumVal,
                            padVarVal,
                            paiObjBeg,
                            panObjLen,
                            paiConBeg,
                            panConLen,
                            padLB = NULL,
                            padUB = NULL)                          
{
    if (is.null(pszVarType)) 
    {
        rpszVarType <- as.null(pszVarType)
    }
    else 
    {
        rpszVarType <- as.character(pszVarType)
    }
    if (is.null(paiVars)) 
    {
        rpaiVars <- as.null(paiVars)
    }
    else 
    {
        rpaiVars <- as.integer(paiVars)
    }
    if (is.null(padLB)) 
    {
        rpadLB <- as.null(padLB)
    }
    else 
    {
        rpadLB <- as.numeric(padLB)
    }
    if (is.null(padUB)) 
    {
        rpadUB <- as.null(padUB)
    }
    else 
    {
        rpadUB <- as.numeric(padUB)
    } 

    ans <- .Call("rcLSloadInstruct", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(nObjs),
                 as.integer(nVars),
                 as.integer(nNumbers),
                 as.integer(panObjSense),
                 as.character(pszConType),
                 rpszVarType,
                 as.integer(panInstruct),
                 as.integer(nInstruct),
                 rpaiVars,
                 as.numeric(padNumVal),
                 as.numeric(padVarVal),
                 as.integer(paiObjBeg),
                 as.integer(panObjLen),
                 as.integer(paiConBeg),
                 as.integer(panConLen),
                 rpadLB,
                 rpadUB)         
    return(ans)
}

rLSaddInstruct <- function(model,
                           nCons,
                           nObjs,
                           nVars,
                           nNumbers,
                           panObjSense,
                           pszConType,
                           pszVarType = NULL,
                           panInstruct,
                           nInstruct,
                           paiCons = NULL,
                           padNumVal,
                           padVarVal,
                           paiObjBeg,
                           panObjLen,
                           paiConBeg,
                           panConLen,
                           padLB = NULL,
                           padUB = NULL)                          
{
    if (is.null(pszVarType)) 
    {
        rpszVarType <- as.null(pszVarType)
    }
    else 
    {
        rpszVarType <- as.character(pszVarType)
    }
    if (is.null(paiCons)) 
    {
        rpaiCons <- as.null(paiCons)
    }
    else 
    {
        rpaiCons <- as.integer(paiCons)
    }
    if (is.null(padLB)) 
    {
        rpadLB <- as.null(padLB)
    }
    else 
    {
        rpadLB <- as.numeric(padLB)
    }
    if (is.null(padUB)) 
    {
        rpadUB <- as.null(padUB)
    }
    else 
    {
        rpadUB <- as.numeric(padUB)
    } 

    ans <- .Call("rcLSaddInstruct", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(nObjs),
                 as.integer(nVars),
                 as.integer(nNumbers),
                 as.integer(panObjSense),
                 as.character(pszConType),
                 rpszVarType,
                 as.integer(panInstruct),
                 as.integer(nInstruct),
                 rpaiCons,
                 as.numeric(padNumVal),
                 as.numeric(padVarVal),
                 as.integer(paiObjBeg),
                 as.integer(panObjLen),
                 as.integer(paiConBeg),
                 as.integer(panConLen),
                 rpadLB,
                 rpadUB)         
    return(ans)
}

rLSloadStringData <- function(model,
                              nStrings,
                              paszStringData)                          
{

    ans <- .Call("rcLSloadStringData", PACKAGE = "rLindo", 
                 model,
                 as.integer(nStrings),
                 as.character(paszStringData))         
    return(ans)
}

rLSloadString <- function(model,
                          pszString)                          
{

    ans <- .Call("rcLSloadString", PACKAGE = "rLindo", 
                 model,
                 as.character(pszString))         
    return(ans)
}

rLSdeleteStringData <- function(model)                          
{

    ans <- .Call("rcLSdeleteStringData", PACKAGE = "rLindo", 
                 model)         
    return(ans)
}

rLSdeleteString <- function(model)                          
{

    ans <- .Call("rcLSdeleteString", PACKAGE = "rLindo", 
                 model)         
    return(ans)
}

rLSgetStringValue <- function(model,
                              iString)                          
{

    ans <- .Call("rcLSgetStringValue", PACKAGE = "rLindo", 
                 model,
                 as.integer(iString))         
    return(ans)
}

rLSgetConstraintProperty <- function(model,
                                     ndxCons)                          
{

    ans <- .Call("rcLSgetConstraintProperty", PACKAGE = "rLindo", 
                 model,
                 as.integer(ndxCons))         
    return(ans)
}

rLSsetConstraintProperty <- function(model,
                                     ndxCons,
                                     nConptype)                          
{

    ans <- .Call("rcLSsetConstraintProperty", PACKAGE = "rLindo", 
                 model,
                 as.integer(ndxCons),
                 as.integer(nConptype))         
    return(ans)
}

rLSloadMultiStartSolution <- function(model,
                                      nIndex)                          
{

    ans <- .Call("rcLSloadMultiStartSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(nIndex))         
    return(ans)
}

rLSloadGASolution <- function(model,
                              nIndex)                          
{

    ans <- .Call("rcLSloadGASolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(nIndex))         
    return(ans)
}

#**************************************************************#
# Solver Initialization Routines (9)                           #
#**************************************************************#
rLSloadBasis <- function(model,
                         panCstatus,
                         panRstatus)                          
{

    ans <- .Call("rcLSloadBasis", PACKAGE = "rLindo", 
                 model,
                 as.integer(panCstatus),
                 as.integer(panRstatus))         
    return(ans)
}

rLSloadVarPriorities <- function(model,
                                 panCprior)                          
{

    ans <- .Call("rcLSloadVarPriorities", PACKAGE = "rLindo", 
                 model,
                 as.integer(panCprior))         
    return(ans)
}

rLSreadVarPriorities <- function(model,
                                 pszFname)                          
{

    ans <- .Call("rcLSreadVarPriorities", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))         
    return(ans)
}

rLSloadVarStartPoint <- function(model,
                                 padPrimal)                          
{

    ans <- .Call("rcLSloadVarStartPoint", PACKAGE = "rLindo", 
                 model,
                 as.numeric(padPrimal))         
    return(ans)
}

rLSloadVarStartPointPartial <- function(model,
                                        nCols,
                                        paiCols,
                                        padPrimal)                          
{

    ans <- .Call("rcLSloadVarStartPointPartial", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCols),
                 as.integer(paiCols),
                 as.numeric(padPrimal))         
    return(ans)
}

rLSloadMIPVarStartPoint <- function(model,
                                    padPrimal)                          
{

    ans <- .Call("rcLSloadMIPVarStartPoint", PACKAGE = "rLindo", 
                 model,
                 as.numeric(padPrimal))         
    return(ans)
}

rLSloadMIPVarStartPointPartial <- function(model,
                                           nCols,
                                           paiCols,
                                           paiPrimal)                          
{

    ans <- .Call("rcLSloadMIPVarStartPointPartial", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCols),
                 as.integer(paiCols),
                 as.integer(paiPrimal))         
    return(ans)
}

rLSreadVarStartPoint <- function(model,
                                 pszFname)                          
{

    ans <- .Call("rcLSreadVarStartPoint", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFname))         
    return(ans)
}

rLSloadBlockStructure <- function(model,
                                  nBlock,
								  panRblock,
								  panCblock,
								  nType)                          
{

    ans <- .Call("rcLSloadBlockStructure", PACKAGE = "rLindo", 
                 model,
                 as.integer(nBlock),
                 as.integer(panRblock),
                 as.integer(panCblock),
                 as.integer(nType))         
    return(ans)
}

#**************************************************************#
# Optimization Routines (6)                                    #
#**************************************************************#
rLSoptimize <- function(model,
                        nMethod)
{

    ans <- .Call("rcLSoptimize", PACKAGE = "rLindo", 
                 model,
                 as.integer(nMethod))
    return(ans)
}

rLSsolveMIP <- function(model)
{
    tryCatch
    (
         ans <- .Call("rcLSsolveMIP", PACKAGE = "rLindo", 
                      model)
     )

     return(ans)
}

rLSsolveGOP <- function(model)
{

    ans <- .Call("rcLSsolveGOP", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSoptimizeQP <- function(model)
{

    ans <- .Call("rcLSoptimizeQP", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLScheckConvexity <- function(model)
{

    ans <- .Call("rcLScheckConvexity", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSsolveSBD <- function(model,
                        nStages,
					    panRowStage,
					    panColStage)
{

    ans <- .Call("rcLSsolveSBD", PACKAGE = "rLindo", 
                 model,
                 as.integer(nStages),
                 as.integer(panRowStage),
                 as.integer(panColStage))
    return(ans)
}

#**************************************************************#
# Solution Query Routines (15)                                 #
#**************************************************************#
rLSgetIInfo <- function(model,
                        nQuery)
{

    ans <- .Call("rcLSgetIInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery))
    return(ans)
}

rLSgetDInfo <- function(model,
                        nQuery)
{

    ans <- .Call("rcLSgetDInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery))
    return(ans)
}

rLSgetPrimalSolution <- function(model)
{

    ans <- .Call("rcLSgetPrimalSolution", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetDualSolution <- function(model)
{

    ans <- .Call("rcLSgetDualSolution", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetReducedCosts <- function(model)
{

    ans <- .Call("rcLSgetReducedCosts", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetReducedCostsCone <- function(model)
{

    ans <- .Call("rcLSgetReducedCostsCone", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetSlacks <- function(model)
{

    ans <- .Call("rcLSgetSlacks", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetBasis <- function(model)
{

    ans <- .Call("rcLSgetBasis", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetSolution <- function(model,
                           nWhich)
{

    ans <- .Call("rcLSgetSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(nWhich))
    return(ans)
}

rLSgetMIPPrimalSolution <- function(model)
{

    ans <- .Call("rcLSgetMIPPrimalSolution", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPDualSolution <- function(model)
{

    ans <- .Call("rcLSgetMIPDualSolution", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPReducedCosts <- function(model)
{

    ans <- .Call("rcLSgetMIPReducedCosts", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPSlacks <- function(model)
{

    ans <- .Call("rcLSgetMIPSlacks", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPBasis <- function(model)
{

    ans <- .Call("rcLSgetMIPBasis", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetNextBestMIPSol <- function(model)
{

    ans <- .Call("rcLSgetNextBestMIPSol", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

#**************************************************************#
#    Model Query Routines (29)                                 #
#**************************************************************#
rLSgetLPData <- function(model)
{

    ans <- .Call("rcLSgetLPData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetQCData <- function(model)
{

    ans <- .Call("rcLSgetQCData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetQCDatai <- function(model,
                          iCon)
{

    ans <- .Call("rcLSgetQCDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCon))
    return(ans)
}

rLSgetVarType <- function(model)
{

    ans <- .Call("rcLSgetVarType", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetVarStartPoint <- function(model)
{

    ans <- .Call("rcLSgetVarStartPoint", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetVarStartPointPartial <- function(model)
{

    ans <- .Call("rcLSgetVarStartPointPartial", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPVarStartPointPartial <- function(model)
{

    ans <- .Call("rcLSgetMIPVarStartPointPartial", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetMIPVarStartPoint <- function(model)
{

    ans <- .Call("rcLSgetMIPVarStartPoint", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetSETSData <- function(model)
{

    ans <- .Call("rcLSgetSETSData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetSETSDatai <- function(model,
                            iSet)
{

    ans <- .Call("rcLSgetSETSDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iSet))
    return(ans)
}

rLSgetSemiContData <- function(model)
{

    ans <- .Call("rcLSgetSemiContData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetLPVariableDataj <- function(model,
                                  iVar)
{

    ans <- .Call("rcLSgetLPVariableDataj", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar))
    return(ans)
}

rLSgetVariableNamej <- function(model,
                                iVar)
{

    ans <- .Call("rcLSgetVariableNamej", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar))
    return(ans)
}

rLSgetVariableIndex <- function(model,
                                pszVarName)
{

    ans <- .Call("rcLSgetVariableIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszVarName))
    return(ans)
}

rLSgetConstraintNamei <- function(model,
                                  iCon)
{

    ans <- .Call("rcLSgetConstraintNamei", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCon))
    return(ans)
}

rLSgetConstraintIndex <- function(model,
                                  pszConName)
{

    ans <- .Call("rcLSgetConstraintIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszConName))
    return(ans)
}

rLSgetConstraintDatai <- function(model,
                                  iCon)
{

    ans <- .Call("rcLSgetConstraintDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCon))
    return(ans)
}

rLSgetLPConstraintDatai <- function(model,
                                    iCon)
{

    ans <- .Call("rLSgetLPConstraintDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCon))
    return(ans)
}

rLSgetConeNamei <- function(model,
                            iCone)
{

    ans <- .Call("rcLSgetConeNamei", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCone))
    return(ans)
}

rLSgetConeIndex <- function(model,
                            pszConeName)
{

    ans <- .Call("rcLSgetConeIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszConeName))
    return(ans)
}

rLSgetConeDatai <- function(model,
                            iCone)
{

    ans <- .Call("rcLSgetConeDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCone))
    return(ans)
}

rLSgetNLPData <- function(model)
{

    ans <- .Call("rcLSgetNLPData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetNLPConstraintDatai <- function(model,
                                     iCon)
{

    ans <- .Call("rcLSgetNLPConstraintDatai", PACKAGE = "rLindo", 
                 model,
                 as.integer(iCon))
    return(ans)
}

rLSgetNLPVariableDataj <- function(model,
                                   iVar)
{

    ans <- .Call("rcLSgetNLPVariableDataj", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar))
    return(ans)
}

rLSgetNLPObjectiveData <- function(model)
{

    ans <- .Call("rcLSgetNLPObjectiveData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetDualModel <- function(model,
                            dualmodel)
{

    ans <- .Call("rcLSgetDualModel", PACKAGE = "rLindo", 
                 model,
                 dualmodel)
    return(ans)
}

rLScalinfeasMIPsolution <- function(model,
                                    padPrimalMipsol = NULL)
{
    if (is.null(padPrimalMipsol)) 
    {
        rpadPrimalMipsol <- as.null(padPrimalMipsol)
    }
    else 
    {
        rpadPrimalMipsol <- as.numeric(padPrimalMipsol)
    }

    ans <- .Call("rcLScalinfeasMIPsolution", PACKAGE = "rLindo", 
                 model,
                 rpadPrimalMipsol)
    return(ans)
}

rLSgetRoundMIPsolution <- function(model,
                                   padPrimal = NULL,
                                   iUseOpti)
{
	if (is.null(padPrimal)) 
    {
        rpadPrimal <- as.null(padPrimal)
    }
    else 
    {
        rpadPrimal <- as.numeric(padPrimal)
    }
    
    ans <- .Call("rcLSgetRoundMIPsolution", PACKAGE = "rLindo", 
                 model,
                 rpadPrimal,
                 as.integer(iUseOpti))
    return(ans)
}

rLSgetRangeData <- function(model)
{

    ans <- .Call("rcLSgetRangeData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

#**************************************************************#
# Model Modification Routines (25)                             #
#**************************************************************#
rLSaddConstraints <- function(model,
                              nNumaddcons,
                              pszConTypes,
                              paszConNames = NULL,
                              paiArows,
                              padAcoef,
                              paiAcols,
                              padB)
{
    if (is.null(paszConNames)) 
    {
        rpaszConNames <- as.null(paszConNames)
    }
    else 
    {
        rpaszConNames <- as.character(paszConNames)
    }

    ans <- .Call("rcLSaddConstraints", PACKAGE = "rLindo", 
                 model,
                 as.integer(nNumaddcons),
                 as.character(pszConTypes),
                 rpaszConNames,
                 as.integer(paiArows),
                 as.numeric(padAcoef),
                 as.integer(paiAcols),
                 as.numeric(padB))
    return(ans)
}

rLSaddVariables <- function(model,
                            nNumaddvars,
                            pszVarTypes,
                            paszVarNames = NULL,
                            paiAcols,
                            panAcols = NULL,
                            padAcoef,
                            paiArows,
                            padC,
                            padL = NULL,
                            padU = NULL)
{
    if (is.null(paszVarNames)) 
    {
        rpaszVarNames <- as.null(paszVarNames)
    }
    else 
    {
        rpaszVarNames <- as.character(paszVarNames)
    }

    if (is.null(panAcols)) 
    {
        rpanAcols <- as.null(panAcols)
    }
    else 
    {
        rpanAcols <- as.integer(panAcols)
    }

    if (is.null(padL)) 
    {
        rpadL <- as.null(padL)
    }
    else 
    {
        rpadL <- as.numeric(padL)
    }

    if (is.null(padU)) 
    {
        rpadU <- as.null(padU)
    }
    else 
    {
        rpadU <- as.numeric(padU)
    }

    ans <- .Call("rcLSaddVariables", PACKAGE = "rLindo", 
                 model,
                 as.integer(nNumaddvars),
                 as.character(pszVarTypes),
                 rpaszVarNames,
                 as.integer(paiAcols),
                 rpanAcols,
                 as.numeric(padAcoef),
                 as.integer(paiArows),
                 as.numeric(padC),
                 rpadL,
                 rpadU)
    return(ans)
}

rLSaddCones <- function(model,
                        nCone,
                        pszConeTypes,
                        paszConenames = NULL,
                        paiConebegcol,
                        paiConecols)
{
    if (is.null(paszConenames)) 
    {
        rpaszConenames <- as.null(paszConenames)
    }
    else 
    {
        rpaszConenames <- as.character(paszConenames)
    }

    ans <- .Call("rcLSaddCones", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCone),
                 as.character(pszConeTypes),
                 rpaszConenames,
                 as.integer(paiConebegcol),
                 as.integer(paiConecols))
    return(ans)
}

rLSaddSETS <- function(model,
                       nSETS,
                       pszSETStype,
                       paiCARDnum,
                       paiSETSbegcol,
                       paiSETScols)
{
    ans <- .Call("rcLSaddSETS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSETS),
                 as.character(pszSETStype),
                 as.integer(paiCARDnum),
                 as.integer(paiSETSbegcol),
                 as.integer(paiSETScols))
    return(ans)
}

rLSaddQCterms <- function(model,
                          nQCnonzeros,
                          paiQCconndx,
                          paiQCvarndx1,
                          paiQCvarndx2,
                          padQCcoef)
{
    ans <- .Call("rcLSaddQCterms", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQCnonzeros),
                 as.integer(paiQCconndx),
                 as.integer(paiQCvarndx1),
                 as.integer(paiQCvarndx2),
                 as.numeric(padQCcoef))
    return(ans)
}

rLSdeleteConstraints <- function(model,
                                 nCons,
                                 paiCons)
{
    ans <- .Call("rcLSdeleteConstraints", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(paiCons))
    return(ans)
}

rLSdeleteCones <- function(model,
                           nCones,
                           paiCones)
{
    ans <- .Call("rcLSdeleteCones", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCones),
                 as.integer(paiCones))
    return(ans)
}

rLSdeleteSETS <- function(model,
                          nSETS,
                          paiSETS)
{
    ans <- .Call("rcLSdeleteSETS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSETS),
                 as.integer(paiSETS))
    return(ans)
}

rLSdeleteSemiContVars <- function(model,
                                  nSCVars,
                                  paiSCVars)
{
    ans <- .Call("rcLSdeleteSemiContVars", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSCVars),
                 as.integer(paiSCVars))
    return(ans)
}

rLSdeleteVariables <- function(model,
                               nVars,
                               paiVars)
{
    ans <- .Call("rcLSdeleteVariables", PACKAGE = "rLindo", 
                 model,
                 as.integer(nVars),
                 as.integer(paiVars))
    return(ans)
}

rLSdeleteQCterms <- function(model,
                             nCons,
                             paiCons)
{
    ans <- .Call("rcLSdeleteQCterms", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(paiCons))
    return(ans)
}

rLSdeleteAj <- function(model,
                        iVar1,
                        nRows,
                        paiRows)
{
    ans <- .Call("rcLSdeleteAj", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar1),
                 as.integer(nRows),
                 as.integer(paiRows))
    return(ans)
}

rLSmodifyLowerBounds <- function(model,
                                 nVars,
                                 paiVars,
                                 padL)
{
    ans <- .Call("rcLSmodifyLowerBounds", PACKAGE = "rLindo", 
                 model,
                 as.integer(nVars),
                 as.integer(paiVars),
                 as.numeric(padL))
    return(ans)
}

rLSmodifyUpperBounds <- function(model,
                                 nVars,
                                 paiVars,
                                 padU)
{
    ans <- .Call("rcLSmodifyLowerBounds", PACKAGE = "rLindo", 
                 model,
                 as.integer(nVars),
                 as.integer(paiVars),
                 as.numeric(padU))
    return(ans)
}

rLSmodifyRHS <- function(model,
                         nCons,
                         paiCons,
                         padB)
{
    ans <- .Call("rcLSmodifyRHS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(paiCons),
                 as.numeric(padB))
    return(ans)
}

rLSmodifyObjective <- function(model,
                               nVars,
                               paiVars,
                               padC)
{
    ans <- .Call("rcLSmodifyObjective", PACKAGE = "rLindo", 
                 model,
                 as.integer(nVars),
                 as.integer(paiVars),
                 as.numeric(padC))
    return(ans)
}

rLSmodifyAj <- function(model,
                        iVar1,
                        nRows,
                        paiRows,
                        padAj)
{
    ans <- .Call("rcLSmodifyAj", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar1),
                 as.integer(nRows),
                 as.integer(paiRows),
                 as.numeric(padAj))
    return(ans)
}

rLSmodifyCone <- function(model,
                          cConeType,
                          iConeNum,
                          iConeNnz,
                          paiConeCols)
{
    ans <- .Call("rcLSmodifyCone", PACKAGE = "rLindo", 
                 model,
                 as.character(cConeType),
                 as.integer(iConeNum),
                 as.integer(iConeNnz),
                 as.integer(paiConeCols))
    return(ans)
}

rLSmodifySET <- function(model,
                         cSETtype,
                         iSETnum,
                         iSETnnz,
                         paiSETcols)
{
    ans <- .Call("rcLSmodifySET", PACKAGE = "rLindo", 
                 model,
                 as.character(cSETtype),
                 as.integer(iSETnum),
                 as.integer(iSETnnz),
                 as.integer(paiSETcols))
    return(ans)
}

rLSmodifySemiContVars <- function(model,
                                  nSCVars,
                                  paiSCVars,
                                  padL,
                                  padU)
{
    ans <- .Call("rcLSmodifySemiContVars", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSCVars),
                 as.integer(paiSCVars),
                 as.numeric(padL),
                 as.numeric(padU))
    return(ans)
}

rLSmodifyConstraintType <- function(model,
                                    nCons,
                                    paiCons,
                                    pszConTypes)
{
    ans <- .Call("rcLSmodifyConstraintType", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCons),
                 as.integer(paiCons),
                 as.character(pszConTypes))
    return(ans)
}

rLSmodifyVariableType <- function(model,
                                  nVars,
                                  paiVars,
                                  pszVarTypes)
{
    ans <- .Call("rcLSmodifyVariableType", PACKAGE = "rLindo", 
                 model,
                 as.integer(nVars),
                 as.integer(paiVars),
                 as.character(pszVarTypes))
    return(ans)
}

rLSaddNLPAj <- function(model,
                        iVar1,
                        nRows,
                        paiRows,
                        padAj)
{
    ans <- .Call("rcLSaddNLPAj", PACKAGE = "rLindo", 
                 model,
                 as.integer(iVar1),
                 as.integer(nRows),
                 as.integer(paiRows),
                 as.numeric(padAj))
    return(ans)
}

rLSaddNLPobj <- function(model,
                         nCols,
                         paiCols,
                         padColj)
{
    ans <- .Call("rcLSaddNLPobj", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCols),
                 as.integer(paiCols),
                 as.numeric(padColj))
    return(ans)
}

rLSdeleteNLPobj <- function(model,
                            nCols,
                            paiCols)
{
    ans <- .Call("rcLSdeleteNLPobj", PACKAGE = "rLindo", 
                 model,
                 as.integer(nCols),
                 as.integer(paiCols))
    return(ans)
}

#**************************************************************#
# Model & Solution Analysis Routines (10)                      #
#**************************************************************#
rLSgetConstraintRanges <- function(model)
{
    ans <- .Call("rcLSgetConstraintRanges", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetObjectiveRanges <- function(model)
{
    ans <- .Call("rcLSgetObjectiveRanges", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetBoundRanges <- function(model)
{
    ans <- .Call("rcLSgetBoundRanges", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetBestBounds <- function(model)
{
    ans <- .Call("rcLSgetBestBounds", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfindIIS <- function(model,
                       nLevel)
{
    ans <- .Call("rcLSfindIIS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nLevel))
    return(ans)
}

rLSfindIUS <- function(model,
                       nLevel)
{
    ans <- .Call("rcLSfindIUS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nLevel))
    return(ans)
}

rLSfindBlockStructure <- function(model,
                                  nBlock,
                                  nType)
{
    ans <- .Call("rcLSfindBlockStructure", PACKAGE = "rLindo", 
                 model,
                 as.integer(nBlock),
                 as.integer(nType))
    return(ans)
}

rLSgetIIS <- function(model)
{
    ans <- .Call("rcLSgetIIS", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetIUS <- function(model)
{
    ans <- .Call("rcLSgetIUS", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetBlockStructure <- function(model)
{
    ans <- .Call("rcLSgetBlockStructure", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

#**************************************************************#
# Memory Related Routines(9)                                   #
#**************************************************************#
rLSfreeSolverMemory <- function(model)
{
    ans <- .Call("rcLSfreeSolverMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfreeHashMemory <- function(model)
{
    ans <- .Call("rcLSfreeHashMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfreeSolutionMemory <- function(model)
{
    ans <- .Call("rcLSfreeSolutionMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfreeMIPSolutionMemory <- function(model)
{
    ans <- .Call("rcLSfreeMIPSolutionMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfreeGOPSolutionMemory <- function(model)
{
    ans <- .Call("rcLSfreeGOPSolutionMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSsetProbAllocSizes <- function(model,
                                 n_vars_alloc,
                                 n_cons_alloc,
                                 n_QC_alloc,
                                 n_Annz_alloc,
                                 n_Qnnz_alloc,
                                 n_NLPnnz_alloc)
{
    ans <- .Call("rcLSsetProbAllocSizes", PACKAGE = "rLindo", 
                 model,
                 as.integer(n_vars_alloc),
                 as.integer(n_cons_alloc),
                 as.integer(n_QC_alloc),
                 as.integer(n_Annz_alloc),
                 as.integer(n_Qnnz_alloc),
                 as.integer(n_NLPnnz_alloc))
    return(ans)
}

rLSsetProbNameAllocSizes <- function(model,
                                     n_varname_alloc,
                                     n_rowname_alloc)
{
    ans <- .Call("rcLSsetProbNameAllocSizes", PACKAGE = "rLindo", 
                 model,
                 as.integer(n_varname_alloc),
                 as.integer(n_rowname_alloc))
    return(ans)
}

rLSaddEmptySpacesAcolumns <- function(model,
                                      paiColnnz)
{
    ans <- .Call("rcLSaddEmptySpacesAcolumns", PACKAGE = "rLindo", 
                 model,
                 as.integer(paiColnnz))
    return(ans)
}

rLSaddEmptySpacesNLPAcolumns <- function(model,
                                         paiColnnz)
{
    ans <- .Call("rcLSaddEmptySpacesNLPAcolumns", PACKAGE = "rLindo", 
                 model,
                 as.integer(paiColnnz))
    return(ans)
}

#**************************************************************#
# Stochastic Programming Interface (88)                        #
#**************************************************************#
rLSwriteDeteqMPSFile <- function(model,
                                 pszFilename,
                                 nFormat,
                                 iType)
{
    ans <- .Call("rcLSwriteDeteqMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFilename),
                 as.integer(nFormat),
                 as.integer(iType))
    return(ans)
}

rLSwriteDeteqLINDOFile <- function(model,
                                   pszFilename,
                                   iType)
{
    ans <- .Call("rcLSwriteDeteqLINDOFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszFilename),
                 as.integer(iType))
    return(ans)
}

rLSwriteSMPSFile <- function(model,
                             pszCorefile,
                             pszTimefile,
                             pszStocfile,
                             nMPStype)
{
    ans <- .Call("rcLSwriteSMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszCorefile),
                 as.character(pszTimefile),
                 as.character(pszStocfile),
                 as.integer(nMPStype))
    return(ans)
}

rLSreadSMPSFile <- function(model,
                            pszCorefile,
                            pszTimefile,
                            pszStocfile,
                            nMPStype)
{
    ans <- .Call("rcLSreadSMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszCorefile),
                 as.character(pszTimefile),
                 as.character(pszStocfile),
                 as.integer(nMPStype))
    return(ans)
}

rLSwriteSMPIFile <- function(model,
                             pszCorefile,
                             pszTimefile,
                             pszStocfile)
{
    ans <- .Call("rcLSwriteSMPIFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszCorefile),
                 as.character(pszTimefile),
                 as.character(pszStocfile))
    return(ans)
}

rLSreadSMPIFile <- function(model,
                            pszCorefile,
                            pszTimefile,
                            pszStocfile)
{
    ans <- .Call("rcLSreadSMPIFile", PACKAGE = "rLindo", 
                 model,
                 as.character(pszCorefile),
                 as.character(pszTimefile),
                 as.character(pszStocfile))
    return(ans)
}

rLSwriteScenarioSolutionFile <- function(model,
                                         jScenario,
                                         pszFname)
{
    ans <- .Call("rcLSwriteScenarioSolutionFile", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.character(pszFname))
    return(ans)
}

rLSwriteNodeSolutionFile <- function(model,
                                     jScenario,
                                     iStage,
                                     pszFname)
{
    ans <- .Call("rcLSwriteNodeSolutionFile", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iStage),
                 as.character(pszFname))
    return(ans)
}

rLSwriteScenarioMPIFile <- function(model,
                                    jScenario,
                                    pszFname)
{
    ans <- .Call("rcLSwriteScenarioMPIFile", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.character(pszFname))
    return(ans)
}

rLSwriteScenarioMPSFile <- function(model,
                                    jScenario,
                                    pszFname,
                                    nFormat)
{
    ans <- .Call("rcLSwriteScenarioMPSFile", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.character(pszFname),
                 as.integer(nFormat))
    return(ans)
}

rLSwriteScenarioLINDOFile <- function(model,
                                      jScenario,
                                      pszFname)
{
    ans <- .Call("rcLSwriteScenarioLINDOFile", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.character(pszFname))
    return(ans)
}

rLSsetModelStocDouParameter <- function(model,
                                        iPar,
                                        dVal)
{
    ans <- .Call("rcLSsetModelStocDouParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(iPar),
                 as.numeric(dVal))
    return(ans)
}

rLSgetModelStocDouParameter <- function(model,
                                        iPar)
{
    ans <- .Call("rcLSgetModelStocDouParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(iPar))
    return(ans)
}

rLSsetModelStocIntParameter <- function(model,
                                        iPar,
                                        iVal)
{
    ans <- .Call("rcLSsetModelStocIntParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(iPar),
                 as.integer(iVal))
    return(ans)
}

rLSgetModelStocIntParameter <- function(model,
                                        iPar)
{
    ans <- .Call("rcLSgetModelStocIntParameter", PACKAGE = "rLindo", 
                 model,
                 as.integer(iPar))
    return(ans)
}

rLSgetScenarioIndex <- function(model,
                                pszName)
{
    ans <- .Call("rcLSgetScenarioIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszName))
    return(ans)
}

rLSgetStageIndex <- function(model,
                             pszName)
{
    ans <- .Call("rcLSgetStageIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszName))
    return(ans)
}

rLSgetStocParIndex <- function(model,
                               pszName)
{
    ans <- .Call("rcLSgetStocParIndex", PACKAGE = "rLindo", 
                 model,
                 as.character(pszName))
    return(ans)
}

rLSgetStocParName <- function(model,
                              nIndex)
{
    ans <- .Call("rcLSgetStocParName", PACKAGE = "rLindo", 
                 model,
                 as.integer(nIndex))
    return(ans)
}

rLSgetScenarioName <- function(model,
                               nIndex)
{
    ans <- .Call("rcLSgetScenarioName", PACKAGE = "rLindo", 
                 model,
                 as.integer(nIndex))
    return(ans)
}

rLSgetStageName <- function(model,
                            nIndex)
{
    ans <- .Call("rcLSgetStageName", PACKAGE = "rLindo", 
                 model,
                 as.integer(nIndex))
    return(ans)
}

rLSgetStocIInfo <- function(model,
                            nQuery,
                            nParam)
{
    ans <- .Call("rcLSgetStocIInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nParam))
    return(ans)
}

rLSgetStocDInfo <- function(model,
                            nQuery,
                            nParam)
{
    ans <- .Call("rcLSgetStocDInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nParam))
    return(ans)
}

rLSgetStocSInfo <- function(model,
                            nQuery,
                            nParam)
{
    ans <- .Call("rcLSgetStocDInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nParam))
    return(ans)
}

rLSgetStocCCPIInfo <- function(model,
                               nQuery,
                               nScenarioIndex,
                               nCPPIndex)
{
    ans <- .Call("rcLSgetStocCCPIInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nScenarioIndex),
                 as.integer(nCPPIndex))
    return(ans)
}

rLSgetStocCCPDInfo <- function(model,
                               nQuery,
                               nScenarioIndex,
                               nCPPIndex)
{
    ans <- .Call("rcLSgetStocCCPDInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nScenarioIndex),
                 as.integer(nCPPIndex))
    return(ans)
}

rLSgetStocCCPSInfo <- function(model,
                               nQuery,
                               nScenarioIndex,
                               nCPPIndex)
{
    ans <- .Call("rcLSgetStocCCPSInfo", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(nScenarioIndex),
                 as.integer(nCPPIndex))
    return(ans)
}

rLSloadSampleSizes <- function(model,
                               panSampleSize)
{
    ans <- .Call("rcLSloadSampleSizes", PACKAGE = "rLindo", 
                 model,
                 as.integer(panSampleSize))
    return(ans)
}

rLSloadConstraintStages <- function(model,
                                    panStage)
{
    ans <- .Call("rcLSloadConstraintStages", PACKAGE = "rLindo", 
                 model,
                 as.integer(panStage))
    return(ans)
}

rLSloadVariableStages <- function(model,
                                  panStage)
{
    ans <- .Call("rcLSloadVariableStages", PACKAGE = "rLindo", 
                 model,
                 as.integer(panStage))
    return(ans)
}

rLSloadStageData <- function(model,
                             numStages,
                             panRstage,
                             panCstage)
{
    ans <- .Call("rcLSloadStageData", PACKAGE = "rLindo", 
                 model,
                 as.integer(numStages),
                 as.integer(panRstage),
                 as.integer(panCstage))
    return(ans)
}

rLSloadStocParData <- function(model,
                               panSparStage,
                               padSparValue)
{
    ans <- .Call("rcLSloadStocParData", PACKAGE = "rLindo", 
                 model,
                 as.integer(panSparStage),
                 as.numeric(padSparValue))
    return(ans)
}

rLSloadStocParNames <- function(model,
                                nSvars,
                                paszSVarNames)
{
    ans <- .Call("rcLSloadStocParNames", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSvars),
                 as.character(paszSVarNames))
    return(ans)
}

rLSgetDeteqModel <- function(model,
                             iDeqType)
{
    ans <- .Call("rcLSgetDeteqModel", PACKAGE = "rLindo", 
                 model,
                 as.integer(iDeqType))
    return(ans)
}

rLSaggregateStages <- function(model,
                               panScheme,
                               nLength)
{
    ans <- .Call("rcLSaggregateStages", PACKAGE = "rLindo", 
                 model,
                 as.integer(panScheme),
                 as.integer(nLength))
    return(ans)
}

rLSgetStageAggScheme <- function(model)
{
    ans <- .Call("rcLSgetStageAggScheme", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSdeduceStages <- function(model,
                            nMaxStage,
                            panRowStagesIn,
                            panColStagesIn,
                            panSparStage)
{
    ans <- .Call("rcLSdeduceStages", PACKAGE = "rLindo", 
                 model,
                 as.integer(nMaxStage),
                 as.integer(panRowStagesIn),
                 as.integer(panColStagesIn),
                 as.integer(panSparStage))
    return(ans)
}

rLSsolveSP <- function(model)
{
    ans <- .Call("rcLSsolveSP", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSsolveHS <- function(model,
                       nSearchMethod)
{
    ans <- .Call("rcLSsolveHS", PACKAGE = "rLindo", 
                 model,
                 as.integer(nSearchMethod))
    return(ans)
}

rLSgetScenarioObjective <- function(model,
                                    jScenario)
{
    ans <- .Call("rcLSgetScenarioObjective", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetNodePrimalSolution <- function(model,
                                     jScenario,
                                     iStage)
{
    ans <- .Call("rcLSgetNodePrimalSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iStage))
    return(ans)
}

rLSgetNodeDualSolution <- function(model,
                                   jScenario,
                                   iStage)
{
    ans <- .Call("rcLSgetNodeDualSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iStage))
    return(ans)
}

rLSgetNodeReducedCost <- function(model,
                                  jScenario,
                                  iStage)
{
    ans <- .Call("rcLSgetNodeReducedCost", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iStage))
    return(ans)
}

rLSgetNodeSlacks <- function(model,
                             jScenario,
                             iStage)
{
    ans <- .Call("rcLSgetNodeSlacks", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iStage))
    return(ans)
}

rLSgetScenarioPrimalSolution <- function(model,
                                         jScenario)
{
    ans <- .Call("rcLSgetScenarioPrimalSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetScenarioReducedCost <- function(model,
                                      jScenario)
{
    ans <- .Call("rcLSgetScenarioReducedCost", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetScenarioDualSolution <- function(model,
                                       jScenario)
{
    ans <- .Call("rcLSgetScenarioDualSolution", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetScenarioSlacks <- function(model,
                                 jScenario)
{
    ans <- .Call("rcLSgetScenarioSlacks", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetNodeListByScenario <- function(model,
                                     jScenario)
{
    ans <- .Call("rcLSgetNodeListByScenario", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetProbabilityByScenario <- function(model,
                                        jScenario)
{
    ans <- .Call("rcLSgetProbabilityByScenario", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetProbabilityByNode <- function(model,
                                    iNode)
{
    ans <- .Call("rcLSgetProbabilityByNode", PACKAGE = "rLindo", 
                 model,
                 as.integer(iNode))
    return(ans)
}

rLSgetStocParData <- function(model)
{
    ans <- .Call("rcLSgetStocParData", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSaddDiscreteBlocks <- function(model,
                                 iStage,
                                 nRealzBlock,
                                 padProb,
                                 pakStart,
                                 paiRows,
                                 paiCols,
                                 paiStvs,
                                 padVals,
                                 nModifyRule)
{
    ans <- .Call("rcLSaddDiscreteBlocks", PACKAGE = "rLindo", 
                 model,
                 as.integer(iStage),
                 as.integer(nRealzBlock),
                 as.numeric(padProb),
                 as.integer(pakStart),
                 as.integer(paiRows),
                 as.integer(paiCols),
                 as.integer(paiStvs),
                 as.numeric(padVals),
                 as.integer(nModifyRule))
    return(ans)
}

rLSaddScenario <- function(model,
                           jScenario,
                           iParentScen,
                           iStage,
                           dProb,
                           nElems,
                           paiRows,
                           paiCols,
                           paiStvs,
                           padVals,
                           nModifyRule)
{
    ans <- .Call("rcLSaddScenario", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario),
                 as.integer(iParentScen),
                 as.integer(iStage),
                 as.numeric(dProb),
                 as.integer(nElems),
                 as.integer(paiRows),
                 as.integer(paiCols),
                 as.integer(paiStvs),
                 as.numeric(padVals),
                 as.integer(nModifyRule))
    return(ans)
}

rLSaddDiscreteIndep <- function(model,
                                iRow,
                                jCol,
                                iStv,
                                nRealizations,
                                padProbs,
                                padVals,
                                nModifyRule)
{
    ans <- .Call("rcLSaddDiscreteIndep", PACKAGE = "rLindo", 
                 model,
                 as.integer(iRow),
                 as.integer(jCol),
                 as.integer(iStv),
                 as.integer(nRealizations),
                 as.numeric(padProbs),
                 as.numeric(padVals),
                 as.integer(nModifyRule))
    return(ans)
}

rLSaddParamDistIndep <- function(model,
                                 iRow,
                                 jCol,
                                 iStv,
                                 nDistType,
                                 nParams,
                                 padParams,
                                 iModifyRule)
{
    ans <- .Call("rcLSaddParamDistIndep", PACKAGE = "rLindo", 
                 model,
                 as.integer(iRow),
                 as.integer(jCol),
                 as.integer(iStv),
                 as.integer(nDistType),
                 as.integer(nParams),
                 as.numeric(padParams),
                 as.integer(iModifyRule))
    return(ans)
}

rLSaddChanceConstraint <- function(model,
                                   iSense,
                                   nCons,
                                   paiCons,
                                   dPrLevel,
                                   dObjWeight)
{
    ans <- .Call("rcLSaddChanceConstraint", PACKAGE = "rLindo", 
                 model,
                 as.integer(iSense),
                 as.integer(nCons),
                 as.integer(paiCons),
                 as.numeric(dPrLevel),
                 as.integer(dObjWeight))
    return(ans)
}

rLSsetNumStages <- function(model,
                            numStages)
{
    ans <- .Call("rcLSsetNumStages", PACKAGE = "rLindo", 
                 model,
                 as.integer(numStages))
    return(ans)
}

rLSgetStocParOutcomes <- function(model,
                                  jScenario)
{
    ans <- .Call("rcLSgetStocParOutcomes", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSloadCorrelationMatrix <- function(model,
                                     nDim,
                                     nCorrType,
                                     nQCnnz,
                                     paiQCcols1,
                                     paiQCcols2,
                                     padQCcoef)
{
    ans <- .Call("rcLSloadCorrelationMatrix", PACKAGE = "rLindo", 
                 model,
                 as.integer(nDim),
                 as.integer(nCorrType),
                 as.integer(nQCnnz),
                 as.integer(paiQCcols1),
                 as.integer(paiQCcols2),
                 as.numeric(padQCcoef))
    return(ans)
}

rLSgetCorrelationMatrix <- function(model,
                                    iFlag,
                                    nCorrType)
{
    ans <- .Call("rcLSgetCorrelationMatrix", PACKAGE = "rLindo", 
                 model,
                 as.integer(iFlag),
                 as.integer(nCorrType))
    return(ans)
}

rLSgetStocParSample <- function(model,
                                iStv,
                                iRow,
                                jCol)
{
    ans <- .Call("rcLSgetStocParSample", PACKAGE = "rLindo", 
                 model,
                 as.integer(iStv),
                 as.integer(iRow),
                 as.integer(jCol))
    return(ans)
}

rLSgetDiscreteBlocks <- function(model,
                                 iEvent)
{
    ans <- .Call("rcLSgetDiscreteBlocks", PACKAGE = "rLindo", 
                 model,
                 as.integer(iEvent))
    return(ans)
}

rLSgetDiscreteBlockOutcomes <- function(model,
                                        iEvent,
                                        iRealz)
{
    ans <- .Call("rcLSgetDiscreteBlockOutcomes", PACKAGE = "rLindo", 
                 model,
                 as.integer(iEvent),
                 as.integer(iRealz))
    return(ans)
}

rLSgetDiscreteIndep <- function(model,
                                iEvent)
{
    ans <- .Call("rcLSgetDiscreteIndep", PACKAGE = "rLindo", 
                 model,
                 as.integer(iEvent))
    return(ans)
}

rLSgetParamDistIndep <- function(model,
                                 iEvent)
{
    ans <- .Call("rcLSgetParamDistIndep", PACKAGE = "rLindo", 
                 model,
                 as.integer(iEvent))
    return(ans)
}

rLSgetScenario <- function(model,
                           jScenario)
{
    ans <- .Call("rcLSgetScenario", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSgetChanceConstraint <- function(model,
                                   iChance)
{
    ans <- .Call("rcLSgetChanceConstraint", PACKAGE = "rLindo", 
                 model,
                 as.integer(iChance))
    return(ans)
}

rLSgetSampleSizes <- function(model)
{
    ans <- .Call("rcLSgetSampleSizes", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetConstraintStages <- function(model)
{
    ans <- .Call("rcLSgetConstraintStages", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetVariableStages <- function(model)
{
    ans <- .Call("rcLSgetVariableStages", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetStocRowIndices <- function(model)
{
    ans <- .Call("rcLSgetStocRowIndices", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSsetStocParRG <- function(model,
                            iStv,
                            iRow,
                            jCol,
                            pRG)
{
    ans <- .Call("rcLSsetStocParRG", PACKAGE = "rLindo", 
                 model,
                 as.integer(iStv),
                 as.integer(iRow),
                 as.integer(jCol),
                 pRG)
    return(ans)
}

rLSgetScenarioModel <- function(model,
                                jScenario)
{
    ans <- .Call("rcLSgetScenarioModel", PACKAGE = "rLindo", 
                 model,
                 as.integer(jScenario))
    return(ans)
}

rLSfreeStocMemory <- function(model)
{
    ans <- .Call("rcLSfreeStocMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSfreeStocHashMemory <- function(model)
{
    ans <- .Call("rcLSfreeStocHashMemory", PACKAGE = "rLindo", 
                 model)
    return(ans)
}

rLSgetModelStocParameterInt <- function(model,
                                        nQuery)
{
    ans <- .Call("rcLSgetModelStocParameterInt", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery))
    return(ans)
}

rLSgetModelStocParameterDou <- function(model,
                                        nQuery)
{
    ans <- .Call("rcLSgetModelStocParameterDou", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery))
    return(ans)
}

rLSgetModelStocParameterChar <- function(model,
                                         nQuery)
{
    ans <- .Call("rcLSgetModelStocParameterChar", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery))
    return(ans)
}

rLSsetModelStocParameterInt <- function(model,
                                        nQuery,
                                        pnResult)
{
    ans <- .Call("rcLSsetModelStocParameterInt", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.integer(pnResult))
    return(ans)
}

rLSsetModelStocParameterDou <- function(model,
                                        nQuery,
                                        pdResult)
{
    ans <- .Call("rcLSsetModelStocParameterInt", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.numeric(pdResult))
    return(ans)
}

rLSsetModelStocParameterChar <- function(model,
                                         nQuery,
                                         pacResult)
{
    ans <- .Call("rcLSsetModelStocParameterChar", PACKAGE = "rLindo", 
                 model,
                 as.integer(nQuery),
                 as.character(pacResult))
    return(ans)
}

rLSgetEnvStocParameterInt <- function(env,
                                      nQuery)
{
    ans <- .Call("rcLSgetEnvStocParameterInt", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery))
    return(ans)
}

rLSgetEnvStocParameterDou <- function(env,
                                      nQuery)
{
    ans <- .Call("rcLSgetEnvStocParameterDou", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery))
    return(ans)
}

rLSgetEnvStocParameterChar <- function(env,
                                       nQuery)
{
    ans <- .Call("rcLSgetEnvStocParameterChar", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery))
    return(ans)
}

rLSsetEnvStocParameterInt <- function(env,
                                      nQuery,
                                      pnResult)
{
    ans <- .Call("rcLSsetEnvStocParameterInt", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery),
                 as.integer(pnResult))
    return(ans)
}

rLSsetEnvStocParameterDou <- function(env,
                                      nQuery,
                                      pdResult)
{
    ans <- .Call("rcLSsetEnvStocParameterInt", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery),
                 as.numeric(pdResult))
    return(ans)
}

rLSsetEnvStocParameterChar <- function(env,
                                       nQuery,
                                       pacResult)
{
    ans <- .Call("rcLSsetEnvStocParameterChar", PACKAGE = "rLindo", 
                 env,
                 as.integer(nQuery),
                 as.character(pacResult))
    return(ans)
}

#**************************************************************#
# Statistical Calculations Interface (15)                      #
#**************************************************************#
rLSsampCreate <- function(env,
                          nDistType)
{
    ans <- .Call("rcLSsampCreate", PACKAGE = "rLindo", 
                 env,
                 as.integer(nDistType))
    return(ans)
}

rLSsampDelete <- function(sample)
{
    ans <- .Call("rcLSsampDelete", PACKAGE = "rLindo", 
                 sample)
    return(ans)
}

rLSsampSetDistrParam <- function(sample,
                                 nIndex,
                                 dValue)
{
    ans <- .Call("rcLSsampSetDistrParam", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nIndex),
                 as.numeric(dValue))
    return(ans)
}

rLSsampGetDistrParam <- function(sample,
                                 nIndex)
{
    ans <- .Call("rcLSsampGetDistrParam", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nIndex))
    return(ans)
}

rLSsampEvalDistr <- function(sample,
                             nFuncType,
                             dXval)
{
    ans <- .Call("rcLSsampEvalDistr", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nFuncType),
                 as.numeric(dXval))
    return(ans)
}

rLSsampEvalUserDistr <- function(sample,
                                 nFuncType,
                                 padXval,
                                 nDim)
{
    ans <- .Call("rcLSsampEvalUserDistr", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nFuncType),
                 as.numeric(padXval),
                 as.integer(nDim))
    return(ans)
}

rLSsampSetRG <- function(sample,
                         RG)
{
    ans <- .Call("rcLSsampSetRG", PACKAGE = "rLindo", 
                 sample,
                 RG)
    return(ans)
}

rLSsampGenerate <- function(sample,
                            nMethod,
                            nSize)
{
    ans <- .Call("rcLSsampGenerate", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nMethod),
                 as.integer(nSize))
    return(ans)
}

rLSsampGetPoints <- function(sample)
{
    ans <- .Call("rcLSsampGetPoints", PACKAGE = "rLindo", 
                 sample)
    return(ans)
}

rLSsampLoadPoints <- function(sample,
                              nSampSize,
                              padXval)
{
    ans <- .Call("rcLSsampLoadPoints", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nSampSize),
                 as.numeric(padXval))
    return(ans)
}

rLSsampGetCIPoints <- function(sample)
{
    ans <- .Call("rcLSsampGetCIPoints", PACKAGE = "rLindo", 
                 sample)
    return(ans)
}

rLSsampLoadDiscretePdfTable <- function(sample,
                                        nLen,
                                        padProb,
                                        padVals)
{
    ans <- .Call("rcLSsampLoadDiscretePdfTable", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nLen),
                 as.numeric(padProb),
                 as.numeric(padVals))
    return(ans)
}

rLSsampGetDiscretePdfTable <- function(sample)
{
    ans <- .Call("rcLSsampGetDiscretePdfTable", PACKAGE = "rLindo", 
                 sample)
    return(ans)
}

rLSsampGetIInfo <- function(sample,
                            nQuery)
{
    ans <- .Call("rcLSsampGetIInfo", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nQuery))
    return(ans)
}

rLSsampGetDInfo <- function(sample,
                            nQuery)
{
    ans <- .Call("rcLSsampGetDInfo", PACKAGE = "rLindo", 
                 sample,
                 as.integer(nQuery))
    return(ans)
}

#**************************************************************#
# Random Number Generation Interface (12)                      #
#**************************************************************#
rLScreateRG <- function(env,
                        nMethod)
{
    ans <- .Call("rcLScreateRG", PACKAGE = "rLindo", 
                 env,
                 as.integer(nMethod))
    return(ans)
}

rLScreateRGMT <- function(env,
                          nMethod)
{
    ans <- .Call("rcLScreateRGMT", PACKAGE = "rLindo", 
                 env,
                 as.integer(nMethod))
    return(ans)
}

rLSgetDoubleRV <- function(rg)
{
    ans <- .Call("rcLSgetDoubleRV", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

rLSgetInt32RV <- function(rg,
                          iLow,
                          iHigh)
{
    ans <- .Call("rcLSgetInt32RV", PACKAGE = "rLindo", 
                 rg,
                 as.integer(iLow),
                 as.integer(iHigh))
    return(ans)
}

rLSsetRGSeed <- function(rg,
                         nSeed)
{
    ans <- .Call("rcLSsetRGSeed", PACKAGE = "rLindo", 
                 rg,
                 as.integer(nSeed))
    return(ans)
}

rLSdisposeRG <- function(rg)
{
    ans <- .Call("rcLSdisposeRG", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

rLSsetDistrParamRG <- function(rg,
                               iParam,
                               dParam)
{
    ans <- .Call("rcLSsetDistrParamRG", PACKAGE = "rLindo", 
                 rg,
                 as.integer(iParam),
                 as.numeric(dParam))
    return(ans)
}

rLSsetDistrRG <- function(rg,
                          nDistType)
{
    ans <- .Call("rcLSsetDistrRG", PACKAGE = "rLindo", 
                 rg,
                 as.integer(nDistType))
    return(ans)
}

rLSgetDistrRV <- function(rg)
{
    ans <- .Call("rcLSgetDistrRV", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

rLSgetInitSeed <- function(rg)
{
    ans <- .Call("rcLSgetInitSeed", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

rLSgetRGNumThreads <- function(rg)
{
    ans <- .Call("rcLSgetRGNumThreads", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

rLSfillRGBuffer <- function(rg)
{
    ans <- .Call("rcLSfillRGBuffer", PACKAGE = "rLindo", 
                 rg)
    return(ans)
}

#**************************************************************#
# Sprint Interface (1)                                         #
#**************************************************************#
rLSsolveFileLP <- function(model,
                           szFileNameMPS,
                           szFileNameSol,
                           nNoOfColsEvaluatedPerSet,
                           nNoOfColsSelectedPerSet,
                           nTimeLimitSec)
{
    ans <- .Call("rcLSsolveFileLP", PACKAGE = "rLindo", 
                 model,
                 as.character(szFileNameMPS),
                 as.character(szFileNameSol),
                 as.integer(nNoOfColsEvaluatedPerSet),
                 as.integer(nNoOfColsSelectedPerSet),
                 as.integer(nTimeLimitSec))
    return(ans)
}

#**************************************************************#
# Branch and price (1)                                         #
#**************************************************************#
rLSsolveMipBnp <- function(model,
                           nBlock,
                           pszFname)
{
    ans <- .Call("rcLSsolveMipBnp", PACKAGE = "rLindo", 
                 model,
                 as.integer(nBlock),
                 as.character(pszFname))
    return(ans)
}