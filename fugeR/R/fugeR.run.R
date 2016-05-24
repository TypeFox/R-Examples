######################################################################################
# fugeR.run
# TO DO: add support for more than 2 membership function for input variables of the fuzzy system.
#
#'   R based Fuzzy logic evolutionary algorithm
#'   
#'   A machine learning algorithm for fuzzy system. \cr
#'   This function use a genetic algorithm in order to construct a fuzzy system able
#'   to fit the values given as \code{labels}. The \code{data} and \code{labels} are used
#'   has learning data.
#'
#'   This is a fuzzy system evolutionnary algotithm. A genetic algorithm is used to find a
#'   fuzzy system able to fit the the data given as labels. \cr
#'   The genetic algorithm generate a a random \code{population} of fuzzy system.
#'   At each \code{generation} all the fuzzy system are tested with the data.
#'   Their predictions are then compared with the labels and a "performance" is given at each system.
#'   The top best system (\code{elitsm} are taken without modification for the next generation.
#'   The population is then used to generate the population for the next generation using crossover
#'   and mutation. At the end of the process (at the last generation) the fuzzy system that obtained
#'   the best performance is returned.
#'   
#'    
#'   @param data            [NULL] Data frame to be used for training (only numeric values are supported).
#'   @param labels          [NULL] Labels of \code{data} (only numeric values are supported).
#'   @param maxRules        [4] Maximum number of rule.
#'   @param maxVarPerRule   [3] Maximum number of input variable per rule.
#'   @param labelsMf        [2] Number of singleton for output variable membership function.
#'   @param population      [200] The population size.
#'   @param elitism         [NA] The number of chromosomes that are kept into the next generation.
#'                               By default is about 20\% of the population size.
#'   @param mutation        [0.01] The chance that a gene in the chromosome mutates.
#'   @param generation      [100] The number of generation made by the genetic algorithm.
#'   @param sensiW          [1.0] The weight of the sensitivity in the fitness function.
#'   @param speciW          [1.0] The weight of the specificity in the fitness function.
#'   @param accuW           [0.0] The weight of the accuracy in the fitness function.
#'   @param threshold       [0.5] The threshold to apply in order to calculate sensitivity, specificity and accuracy.
#'   @param rmseW           [0.2] The weight of the "root mean square error" between labels and
#'                                values predicted by the fuzzy system.
#'   @param verbose         [FALSE] If true the algorithm will be more verbose. By default False.
#' 
#'   @return  \code{fis}, A list containing the logs of the evolution, the peformances of the best system and
#'                        its description.
#'                        \item{inputVarIds}{The IDs of the variable used in the fuzzy system}
#'                        \item{inputMfIds}{The IDs of the membership function used by each variable in the fuzzy system}
#'                        \item{inputMfs}{The values used to caclculate the membership functions}
#'                        \item{outputVarIds}{The IDs of each output variable of each rule}
#'                        \item{outputMfIds}{The IDs of the membership function used by each output variable}
#'                        \item{outputMfs}{the value used to compute the membership functions of the output variables}
#'                        \item{fitness}{The fitness value reached by the best fuzzy system}
#'                        \item{mse}{The Mean Square Error of the best fuzzy system}
#'                        \item{rmse}{The Root Mean Square Error between labels and the prediction made by the best fuzzy system}
#'                        \item{accu}{The accuracy of the prediction made by the best fuzzy system
#'                                    (only if a threshold different of NA was given as argument)}
#'                        \item{sensi}{The sensitivity of the prediction made by the best fuzzy system
#'                                    (only if a threshold different of NA was given as argument)}
#'                        \item{speci}{The specificity of the prediction made by the best fuzzy system
#'                                    (only if a threshold different of NA was given as argument)}
#'                        \item{evo}{A list containing the evolution logs}
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
#'   }
#'
#' @seealso  \code{\link{fugeR.sfRun}} \code{\link{fugeR.predict}}\code{\link{fugeR.summary}}
#'           \code{\link{fugeR.save}} \code{\link{fugeR.load}}
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.run <-
function(   data=NULL,
            labels=NULL,
            maxRules=4,
            maxVarPerRule=3,
            labelsMf=2,
            population=200,
            elitism=NA, #about 20% by default
            mutation=0.01,
            generation=100,
            sensiW=1.0,
            speciW=1.0,
            accuW=0.0,
            threshold=0.5,
            rmseW=0.2,
            verbose=FALSE
            )
{
    #Should be an argument when more than 2 mf will be implemented
    dataMf <- 2
    #System size parameters
    fugeR.nbRule <- maxRules #MIN VALUE IS 2!!! (Default Rule inclusive)
    fugeR.nbMaxVarInPerRule <- maxVarPerRule #MIN is 1.. no max
    fugeR.nbVarOut <- ncol(labels) #MIN is 1... no Max
    fugeR.nbInputSet <- dataMf # only 2 is supported for the moment
    fugeR.nbOutputSet <- labelsMf #Min is 1, no max
    
    #SAVE BEST SYSTEM DURING EVOLUTION
    fugeR.ACTUAL_VALUES <- data
    fugeR.TO_PREDICT <- labels
    
    #Fill fitness parameter list
    #if the vector has not the good length
    #then the last value is taken as default value

    if( length(sensiW) != fugeR.nbVarOut) {
      stop('SensiW vector should have a length equal to the number of column in "labels"')
    }
    if( length(speciW) != fugeR.nbVarOut ) {
      stop('SensiW vector should have a length equal to the number of column in "labels"')
    }
    if( length(accuW) != fugeR.nbVarOut ) {
      stop('AccuW vector should have a length equal to the number of column in "labels"')
    }
    if( length(threshold) != fugeR.nbVarOut ) {
      stop('Threshold vector should have a length equal to the number of column in "labels"')
    }
    if( length(rmseW) != fugeR.nbVarOut ) {
      stop('RmseW vector should have a length equal to the number of column in "labels"')
    }
    
    #System size parameters
    assign("fugeR.nbRule"           , maxRules      , envir=fugeRglobal)
    assign("fugeR.nbMaxVarInPerRule", maxVarPerRule , envir=fugeRglobal)
    assign("fugeR.nbVarOut"         , fugeR.nbVarOut, envir=fugeRglobal)
    assign("fugeR.nbInputSet"       , dataMf        , envir=fugeRglobal)
    assign("fugeR.nbOutputSet"      , labelsMf      , envir=fugeRglobal)
    
    #Fitness parameter
    fugeR.lstFitParameter <- list(
      SENSI_W = c(sensiW),
      SPECI_W = c(speciW),
      ACCU_W  = c(accuW),
      THRESH  = c(threshold),
      RMSE_W  = c(rmseW)
    )
    #Fitness parameters
    assign("fugeR.lstFitParameter", fugeR.lstFitParameter, envir=fugeRglobal)
    
    #Data and labels
    assign("fugeR.ACTUAL_VALUES", data, envir=fugeRglobal)
    assign("fugeR.TO_PREDICT", labels, envir=fugeRglobal )
    
    #Best System
    assign("fugeR.BEST_FIT", NA, envir=fugeRglobal)
    assign("fugeR.BEST_SYS", NA, envir=fugeRglobal)
    
    #################### CONSTRUCT THE CHROMOSOME ####################
    #-----------------------------------------------------------------
    chromosomeInputMin  <- rep(0.0, (fugeR.nbMaxVarInPerRule*(2+(fugeR.nbInputSet))))
    chromosomeInputMax  <- rep(0.0, (fugeR.nbMaxVarInPerRule*(2+(fugeR.nbInputSet))))
    chromosomeOutputMin <- rep(0.0, (fugeR.nbVarOut*2))
    chromosomeOutputMax <- rep(0.0, (fugeR.nbVarOut*2))
    #Define input Var ID
    colToSelect <- seq(1, length(chromosomeInputMin), fugeR.nbInputSet+2)
    chromosomeInputMin[colToSelect] <- 1
    chromosomeInputMax[colToSelect] <- ncol(fugeR.ACTUAL_VALUES) + 0.99 # 0 DON'T CARE ADDED
    #Define input MF ID
    colToSelect <- seq(2, length(chromosomeInputMin), fugeR.nbInputSet+2)
    chromosomeInputMin[colToSelect] <- 1
    chromosomeInputMax[colToSelect] <- fugeR.nbInputSet + 0.99
    #Define input MF
    seqMf <- rep(1:fugeR.nbInputSet, fugeR.nbMaxVarInPerRule)
    toAdd <- rep(seq(2, length(chromosomeInputMin), 2+fugeR.nbInputSet), each=fugeR.nbInputSet)
    colToSelect <- seqMf + toAdd
    chromosomeInputMin[colToSelect] <- 0.0
    chromosomeInputMax[colToSelect] <- 1.0
    #Define output Var ID
    colToSelect <- seq(1, length(chromosomeOutputMin), 2)
    chromosomeOutputMin[colToSelect] <- 1
    chromosomeOutputMax[colToSelect] <- fugeR.nbVarOut + 0.99 #+ 5  #5 DON'T CARE ADDED
    #Define output MF ID
    colToSelect <- seq(2, length(chromosomeOutputMin), 2)
    chromosomeOutputMin[colToSelect] <- 1
    chromosomeOutputMax[colToSelect] <- fugeR.nbOutputSet + 0.99

    #Put the input and output chromosome part together
    min <- rep(c(chromosomeInputMin, chromosomeOutputMin), fugeR.nbRule - 1)
    max <- rep(c(chromosomeInputMax, chromosomeOutputMax), fugeR.nbRule - 1)
    
    #Default Rule code MF of all the ouput variable
    chromosomeDefaultMin <- rep(0.0, (fugeR.nbVarOut*(1 + fugeR.nbOutputSet)))
    chromosomeDefaultMax <- rep(0.0, (fugeR.nbVarOut*(1 + fugeR.nbOutputSet)))
    #Define output MF ID
    colToSelect <- seq(1, length(chromosomeDefaultMin), 1 + fugeR.nbOutputSet)
    chromosomeDefaultMin[colToSelect] <- 1
    chromosomeDefaultMax[colToSelect] <- fugeR.nbOutputSet + 0.99
    #Define output var MF
    seqMf <- rep(1:fugeR.nbOutputSet, fugeR.nbVarOut)
    toAdd <- rep(seq(1, length(chromosomeDefaultMin), 1+fugeR.nbOutputSet), each=fugeR.nbOutputSet)
    colToSelect <- seqMf + toAdd
    chromosomeDefaultMin[colToSelect] <- 0.0
    chromosomeDefaultMax[colToSelect] <- 1.0
    
    #Add default rule to the genome
    min <- c(min, chromosomeDefaultMin)
    max <- c(max, chromosomeDefaultMax)
    
    #--------------------------------------------------------------

    
    if(is.na(elitism)) {
        elitism <- as.integer(population * 0.20)
    }

    if(is.na(mutation)) {
        mutation <- 0.01
    }
    
    #Compute (max - min) interval and min for every input and output
    #variable as we need it to find the real value coded by the MFs
    #because the MFs code a relative value include in 0 and 1
    fugeR.minIn <- sapply(data, function(x){min(x, na.rm=TRUE)})
    fugeR.maxIn <- sapply(data, function(x){max(x, na.rm=TRUE)})
    
    fugeR.intervalIn <- fugeR.maxIn - fugeR.minIn
    assign("fugeR.intervalIn", fugeR.intervalIn, envir=fugeRglobal)
    assign("fugeR.minIn", fugeR.minIn, envir=fugeRglobal)
    
    fugeR.minOut <- sapply(labels, min)
    fugeR.maxOut <- sapply(labels, max)
    
    fugeR.intervalOut <- fugeR.maxOut - fugeR.minOut
    assign("fugeR.intervalOut", fugeR.intervalOut, envir=fugeRglobal)
    assign("fugeR.minOut", fugeR.minOut, envir=fugeRglobal)

    #pass the verbose option
    assign("fugeR.verbose", verbose, envir=fugeRglobal)

    #################### LAUNCH THE EVOLUTION ####################
    #-------------------------------------------------------------- 
    evoLog <- fugeR.evo(  min,
                          max,
                          monitorFunc=NULL,
                          evalFunc=fugeR.evaluate,
                          verbose=verbose,
                          mutationChance=mutation,
                          popSize=population,
                          elitism=elitism,
                          iters=generation )
    #--------------------------------------------------------------
    
    toReturn <- fugeRglobal$fugeR.BEST_SYS
    toReturn$inNames = names(data)
    toReturn$outNames = names(labels)

    rm(list = ls(envir=fugeRglobal), envir=fugeRglobal)

    return(toReturn)
}
