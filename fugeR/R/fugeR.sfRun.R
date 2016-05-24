######################################################################################
# fugeR.sfRun
#
#'   The parallel version of fugeR.run using snowfall.
#'   
#'   The parallel version of fugeR.run. Will launch \code{fugeR.run} a number of times
#'   given as argument. This function use \code{\link{snowfall}} package in order to take benefit
#'   of mutli-core computers.
#'
#'   \code{fugeR.sfRun} sould be used when you want to repeat an experience many times. \cr
#'   This is usefull when you are searching the good parameters (maxRules, macVarPerRule) for a problem.
#'   \code{fugeR.sfRun} will launch \code{fugeR.run} and test the obtained system. It automatically resamples the data
#'   using bootstrapping method. \cr
#'   For example if the argument \code{rep} has the value 1000 and the number of sample in data is 100.
#'   FugeR.sfRun resample the data with replacement with the size of the resample equal to 100
#'   (the size of the original data set) this constitute the training set, the samples that were not picked are taken
#'   to create the validation set. FugeR.run is then called with the training set and the obtained fuzzy systems
#'   is tested on the validation set. If \code{rep} value was 1000, this operation is repeated 1000 times. \cr
#'   FugeR.sfRun saves every systems in the directory specified by \code{path} and return a resume of the performance
#'   obtained by each system on their training and validation set.
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
#'   @param path            [NULL] THe path where to save the fuzzy systems.
#'   @param rep             [300] Number fuzzy system to find.
#'   @param parallel        [TRUE] Logical value indicating if the function can run in parralel
#'   @param cpus            [1] number of cpus that can be used
#'         
#'   @return  \code{res}, A data.frame of size \code{rep} containing the performance of each fuzzy system on
#'                        training and validation set.
#'
#'   @examples
#'   ##
#'   ##
#'   \dontrun{
#'      expResume <- fugeR.sfRun (
#'                      In,
#'                      Out,
#'                      generation=100,
#'                      population=200,
#'                      elitism=40,
#'                      verbose=TRUE,
#'                      threshold=0.5,
#'                      sensiW=1.0,
#'                      speciW=1.0,
#'                      accuW=0.0,
#'                      rmseW=1.0,
#'                      maxRules=10,
#'                      maxVarPerRule=2,
#'                      labelsMf=2,
#'                      path=\'./exp\',
#'                      rep=100,
#'                      parallel=TRUE,
#'                      cpus=2
#'                  )
#'   }
#'
#' @seealso  \code{\link{fugeR.run}} \code{\link{fugeR.predict}}\code{\link{fugeR.summary}}
#'           \code{\link{fugeR.save}} \code{\link{fugeR.load}}
#' @author Alexandre Bujard, HEIG-VD, Jul'2012
#'
#' @export
######################################################################################
fugeR.sfRun <-
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
            verbose=FALSE,
            path=NULL,
            rep=300,
            parallel=FALSE,
            cpus=1
            )
{
    #--- Initializate Parallel Computing ---#
    sfInit(parallel=parallel,cpus=cpus)
    
    if(parallel) 
    { 
        sfLibrary("fugeR",character.only=TRUE); 
    }
    
    nbExp <- 1:rep
    
    learning <- vector('list', rep)
    validation <- vector('list', rep)
    
    for( i in nbExp) {
        learning[[i]] <- sample(nrow(data), size=nrow(data), replace=T)
        all <- (1:nrow(data))
        validation[[i]] <- all[-learning[[i]]]
    }
    
    
    sfExport(list=list(
        'data',
        'labels',
        'maxRules',
        'maxVarPerRule',
        'labelsMf',
        'population',
        'elitism',
        'mutation',
        'generation',
        'sensiW',
        'speciW',
        'accuW',
        'threshold',
        'rmseW'))
    
	cat('\n', 'sfRun: START', '\n')
    res <- sfLapply(
        nbExp, function(x){ 
            args <- list (
                data[learning[[x]],,drop=FALSE],
                labels[learning[[x]],,drop=FALSE],
                maxRules,
                maxVarPerRule,
                labelsMf,
                population,
                elitism,
                mutation,
                generation,
                sensiW,
                speciW,
                accuW,
                threshold,
                rmseW,
                verbose=FALSE
            )
            do.call('fugeR.run', args) }
    )
    cat('\n', 'sfRUN: END', '\n')
    #Stop cluster
    sfStop()
    
    #Valid all the system on the validation set
    #the cases not taken for the training set
    #Fitness parameter
    lstFitParameter <- list(
        SENSI_W = c(sensiW),
        SPECI_W = c(speciW),
        ACCU_W  = c(accuW),
        THRESH  = c(threshold),
        RMSE_W  = c(rmseW)
    )
    #Fitness
    valRes <- vector('list',rep)
    for(i in nbExp) {
        prediction <- fugeR.predict(fuzzySystem=res[[i]], dataset=data[validation[[i]],,drop=FALSE])
        valRes[[i]] <- fugeR.evalFitness(prediction, labels[validation[[i]],,drop=FALSE], res[[i]], lstFitParameter)
    }
    
	cat('\n', 'sfRUN: Save systems', '\n')
	#Save the fuzzy systems
    for(i in nbExp) {
        fugeR.save(res[[i]],file.path(path,paste('sys-',i,'.R',sep='')))
    }
    
    #Save all the training set 
    dput(learning, file.path(path,paste('learningSet.R',sep='')) )
    
    idSys <- c()
    sensiL <- c()
    speciL <- c()
    accuL <- c()
    rmseL <- c()
    
    sensiV <- c()
    speciV <- c()
    accuV <- c()
    rmseV <- c()
    
    for(i in nbExp) {
        idSys <- c(idSys, paste('sys-', i,sep=''))
        sensiL <- c(sensiL, res[[i]]$sensi)
        speciL <- c(speciL, res[[i]]$speci)
        accuL <- c(accuL, res[[i]]$accu)
        rmseL <- c(rmseL, res[[i]]$rmse)
        
        sensiV <- c(sensiV, valRes[[i]]$sensi)
        speciV <- c(speciV, valRes[[i]]$speci)
        accuV <- c(accuV, valRes[[i]]$accu)
        rmseV <- c(rmseV, valRes[[i]]$rmse)
    }
    
    resume <- data.frame(idSys=idSys,
                         sensiLearning=sensiL,
                         speciLearning=speciL,
                         accuLearning=accuL,
                         rmseLearning=rmseL,
                         sensiValidation=sensiV,
                         speciValidation=speciV,
                         accuValidation=accuV,
                         rmseValidation=rmseV)
    
	#TO DO
    #should return a data.frame with the result of each run
    #on training and validation
    #res <- data.frame()
    return(resume)
}
