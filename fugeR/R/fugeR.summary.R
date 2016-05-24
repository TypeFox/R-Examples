######################################################################################
# fugeR.summary
#
#'   Summarize a fuzzy system.
#'   
#'   Show the text description of a fuzzy system in a human readable form.
#'   
#'   @param fuzzySystem [NULL] The fuzzy system to show.
#'
#'   @examples
#'   ##
#'   ##
#'   \dontrun{
#'      fis <- fugeR.run (
#'                  In,
#'                  Out,
#'                  generation=100,
#'                  population=200,
#'                  elitism=40,
#'                  verbose=TRUE,
#'                  threshold=0.5,
#'                  sensiW=1.0,
#'                  speciW=1.0,
#'                  accuW=0.0,
#'                  rmseW=1.0,
#'                  maxRules=10,
#'                  maxVarPerRule=2,
#'                  labelsMf=2
#'              )
#'      
#'      fugeR.summary(fis)
#'   }
#'
#' @seealso  \code{\link{fugeR.run}}
#' 
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.summary <-
function(fuzzySystem=NULL) {
    #should perfrom more check
#     if(fuzzySystem == NULL){
#         stop("Please provide a valid fugeR fuzzy system")
#     }
    
  InNames <- fuzzySystem$inNames
  OutNames <- fuzzySystem$outNames
  
  #antecedent part
  fuzzySystemIn <- list(
    varIdMat = matrix(as.integer(fuzzySystem$inputVarIds), fuzzySystem$nbRule - 1),
    mfIdMat  = matrix(as.integer(fuzzySystem$inputMfIds), fuzzySystem$nbRule - 1),
    mfValMat = fuzzySystem$inputMfs
  )
  
  #consequent part
  fuzzySystemOut <- list(
    varIdMat = matrix(as.integer(fuzzySystem$outputVarIds), fuzzySystem$nbRule - 1),
    mfIdMat  = matrix(as.integer(fuzzySystem$outputMfIds), fuzzySystem$nbRule - 1),
    mfValMat = fuzzySystem$outputMfs
  )
    
  fugeR.nbMaxVarInPerRule <- fuzzySystem$nbMaxIn
  fugeR.nbInputSet <- fuzzySystem$nbInMf
  
  #Find which var are used by the system
  nbVar <- length(InNames)
  nbRule <- fuzzySystem$nbRule
  lstVarUsed <- fuzzySystemIn$varIdMat[fuzzySystemIn$varIdMat %in% 1:nbVar]
  lstVarUsed <- unique(lstVarUsed)
  lstMf <- list()
  lstValue <- list()
  sapply(lstVarUsed, function(x)  {
                        lstMf[[x]] <<- c(NA)
                      } )
  
  lstRule <- list()
  
  
  #------------- CONSEQUENT PART -----------#
  fugeR.nbVarOut <- fuzzySystem$nbOut
  fugeR.nbOutputSet <- fuzzySystem$nbOutMf
  
  #Find var
  lstInRule <- list()
  sapply(1:(nbRule-1), function(x) {
    lstInRule[[x]] <<- fuzzySystemOut$varIdMat[x,]
    lstInRule[[x]] <<- lstInRule[[x]] %in% 1:fugeR.nbVarOut
  } )
  #Find mf ids
  lstMfId <- list()
  sapply(1:(nbRule-1), function(x) {     
    lstMfId[[x]] <<- fuzzySystemOut$mfIdMat[x,]
  } )
  #Find membership functions
  lstMfOut <- list()
  sapply(1:fugeR.nbVarOut, function(x)  {
        lstMfOut[[x]] <<- sort(fuzzySystemOut$mfValMat[x,])
        lstMfOut[[x]] <<- fuzzySystem$minOut[x] + (fuzzySystem$intervalOut[x] * lstMfOut[[x]])
  } )

  
  #----------- THE TEXTUAL REPRESENTATION -----------#
  txtVarDefinition <- c('INPUT VARIABLES : \n')
  txtRules <- c()

  #Rule to compute
  #Minus 1, because the default rule
  nbRuleToCompute <- 1:(nbRule-1)
  cat('FUZZY SYSTEM :\n')
  for(i in nbRuleToCompute) {

    for(j in 1:fugeR.nbMaxVarInPerRule) {
      idVar <- fuzzySystemIn$varIdMat[i,j]
      #if var not used... we skip to next var  
      if(!(idVar %in% lstVarUsed)) {
        next
      }
      #otherwise we fuzzify values corresponding to this var
      
      #if is equal to NA it means we have to define the membership functions
      if(is.na(lstMf[[idVar]][1])) {
        lstMf[[idVar]] <- sort(fuzzySystemIn$mfValMat[i, ((j*fugeR.nbInputSet) - (fugeR.nbInputSet-1)):(j*fugeR.nbInputSet)])
        lstMf[[idVar]] <-  fuzzySystem$minIn[idVar] + (fuzzySystem$intervalIn[idVar] * lstMf[[idVar]])
        
        txtVarDefinition <- c(txtVarDefinition, paste(InNames[idVar], ' : (', sep=''))
        first <- TRUE
        for(mf in lstMf[[idVar]]){
          if(!first){
            txtVarDefinition <- c(txtVarDefinition, paste(', ',sep=''))
          }
          txtVarDefinition <- c(txtVarDefinition, paste(mf,sep=''))
          first <- FALSE
        }
        txtVarDefinition <- c(txtVarDefinition, paste(')\n',sep=''))
      }
      
      #Get the mf list
      mfList <- as.vector(lstMf[[idVar]])
      idMf = fuzzySystemIn$mfIdMat[i,j]
      
      txtRules <- c(txtRules,  paste(InNames[idVar],' is MF_',idMf,sep=''))
    }
    
    cat('RULE (', i, ') : ',sep='')
    cat(txtRules, sep=' AND ')
    cat(' THEN ')
    
    for(h in 1:fugeR.nbVarOut) {
      if (lstInRule[[i]][h]) {
        cat(OutNames[h],' is MF_', lstMfId[[i]][h], sep='')
        cat(' ')
      }
    }
    
    cat('\n')
    txtRules <- c()
    #lstRule[[i]] <- minValue
  }
  
  cat('DEFAULT RULE (', nbRule, ') : ', sep='')
  for(h in 1:fugeR.nbVarOut) {
    cat(OutNames[h],' is MF_', as.integer(fuzzySystem$defautMfIds[h]), sep='')
    #cat(OutNames[h],' is MF_', lstMfId[[nbRule]][h], sep='')
    cat(' ')
  }
  cat('\n\n')
  cat(txtVarDefinition,sep='')
  
  cat('\n')
  cat('OUTPUT VARIABLES')
  cat('\n')
  for(h in 1:fugeR.nbVarOut) {
    cat(OutNames[h],' : (',sep='')
    cat(lstMfOut[[h]], sep=', ')
    cat(')\n')
  }

}
