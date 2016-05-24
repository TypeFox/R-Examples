#' Estimate the connectivity matrix of a causal graph
#'
#' @description Estimates the connectivity matrix of a directed causal graph, 
#' using various possible methods. Supported methods at the moment are 
#' backShift, bivariateANM, bivariateCAM, CAM, hiddenICP, ICP, GES, GIES, LINGAM, 
#' PC, regression and RFCI.
#' @param X A (nxp)-data matrix with n observations of p variables.
#' @param environment An optional vector of length n, where the entry for 
#' observation i is an index for the environment in which observation i took 
#' place (Simplest case: entries \code{1} for observational data and entries
#'  \code{2} for interventional data of unspecified type. Encoding for observational
#'  data can we changed with \code{indexObservationalData}). Is required for 
#'  methods \code{ICP}, \code{hiddenICP} and \code{backShift}.
#' @param interventions A optional list of length n. The entry for observation
#'  i is a numeric vector that specifies the variables on which interventions 
#'  happened for observation i (a scalar if an intervention happened on just 
#'  one variable and \code{numeric(0)} if no intervention occured for this 
#'  observation). Is used for methods \code{gies} and \code{CAM} and will 
#'  generate the vector \code{environment} if the latter is set to \code{NULL}.
#'  (However, this might generate too many different environments for some data 
#'  sets, so a hand-picked vector \code{environment} is preferable). Is also used 
#'  for \code{ICP} and \code{hiddenICP} to exclude interventions on the target 
#'  variable of interest.
#' @param parentsOf The variables for which we would like to estimate the 
#' parents. Default are all variables.
#' @param method A string that specfies the method to use. The methods 
#' \code{pc} (PC-algorithm), \code{LINGAM} (LINGAM), \code{ges} 
#' (Greedy equivalence search), \code{gies} (Greedy interventional equivalence 
#' search) and \code{rfci} (Really fast causal inference) are imported from the 
#' package "pcalg" and are documented there in more detail, including the 
#' additional options that can be supplied via \code{setOptions}. The method 
#' \code{CAM} (Causal additive models) is documented in the package "CAM" and 
#' the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} 
#' (Invariant causal prediction with hidden variables) are from the package 
#' "InvariantCausalPrediction". The method \code{backShift} comes from the 
#' package "backShift". Finally, the methods \code{bivariateANM} and 
#' \code{bivariateCAM} are for now implemented internally but will hopefully 
#' be part of another package at some point in the near future.
#' @param alpha The level at which tests are done. This leads to confidence 
#' intervals for \code{ICP} and \code{hiddenICP} and is used internally for 
#' \code{pc} and \code{rfci}.
#' @param variableSelMat An optional logical matrix of dimension (pxp). An 
#' entry \code{TRUE} for entry (i,j) says that variable i should be considered 
#' as a potential parent for variable j and vice versa for \code{FALSE}. If the 
#' default value of \code{NULL} is used, all variables will be considered, but 
#' this can be very slow, especially for methods \code{pc}, \code{ges}, 
#' \code{gies}, \code{rfci} and \code{CAM}.
#' @param excludeTargetInterventions When looking for parents of variable k 
#' in 1,...,p, set to \code{TRUE} if observations where an intervention on 
#' variable k occured should be excluded. Default is \code{TRUE}.
#' @param onlyObservationalData If set to \code{TRUE}, only observational data 
#' is used. It will take the index in \code{environment} specified by 
#' \code{indexObservationalData}. If \code{environment} is \code{NULL}, all 
#' observations are used. Default is \code{FALSE}.
#' @param indexObservationalData Index in \code{environment} that encodes 
#' observational data. Default is \code{1}.
#' @param returnAsList If set to \code{TRUE}, will return a list, where entry 
#' k is a list containing the estimated parents of variable k. The option 
#' \code{directed} will be ignored if set to \code{TRUE}. Default is 
#' \code{FALSE}.
#' @param pointConf If \code{TRUE}, numerical estimates will be returned if 
#' possible. For methods \code{ICP} and \code{hiddenICP}, these are the values 
#' in the individual confidence intervals (at chosen level \code{alpha}) that 
#' are closest to 0; for other methods these are point estimates. Some methods 
#' do not return numerical point estimates; for these the output will remain 
#' binary 0/1 (no-edge/edge). Default is \code{FALSE}.
#' @param setOptions A list that can take method-specific options; see the 
#' individual documentations of the methods for more options and their 
#' possible values.
#' @param directed If \code{TRUE}, an edge will be returned if and only if an 
#' edge has been detected to be directed (ie entry will be set to 0 for entry 
#' (j,k) if both j->k and k-> j are estimated). Ignored if not the whole graph 
#' is estimated or if \code{returnAsList} is \code{TRUE}.
#' @param verbose If \code{TRUE}, detailed output is provided.
#'
#' @return If option \code{returnAsList} is \code{FALSE}, a sparse matrix, 
#' where a 0 entry in position (j,k) corresponds to an estimate of "no edge" 
#' \code{j} -> \code{k}, while an entry 1 corresponds to an 
#' estimated egde. If option \code{pointConf} is \code{TRUE}, the 1 entries 
#' will be replaced by numerical values that are either point estimates of the 
#' causal coefficients or confidence bounds (see above). 
#' If option \code{returnAsList} is \code{TRUE}, a list will be returned. 
#' The k-th entry in the list is the numeric vector with the indices of the 
#' estimated parents of node \code{k}. 
#' 
#' @author Christina Heinze \email{heinze@@stat.math.ethz.ch}, 
#'  Nicolai Meinshausen \email{meinshausen@@stat.math.ethz.ch}
#' 
#' @seealso \code{\link{getParentsStable}} for stability selection-based 
#' estimation of the causal graph.
#' 
#' @examples
#' ## load the backShift package for data generation and plotting functionality
#' if(!requireNamespace("backShift", quietly = TRUE))
#'  stop("The package 'backShift' is needed for the examples to 
#'  work. Please install it.", call. = FALSE)
#' 
#' require(backShift)
#' 
#' # Simulate data with connectivity matrix A with assumptions 
#' # 1) hidden variables present
#' # 2) precise location of interventions is assumed unknown
#' # 3) different environments can be distinguished
#' 
#' ## simulate data
#' myseed <- 1
#' 
#' # sample size n
#' n <- 10000
#' 
#' # p=3 predictor variables and connectivity matrix A
#' p <- 3
#' labels <- c("1", "2", "3")
#' A <- diag(p)*0
#' A[1,2] <- 0.8
#' A[2,3] <- 0.8
#' A[3,1] <- -0.4  
#' 
#' # divide data in 10 different environments
#' G <- 10
#' 
#' # simulate
#' simResult <- simulateInterventions(n, p, A, G, intervMultiplier = 3, 
#'              noiseMult = 1, nonGauss = TRUE, hiddenVars = TRUE, 
#'              knownInterventions = FALSE, fracVarInt = NULL, simulateObs = TRUE, 
#'              seed = myseed)
#' X <- simResult$X
#' environment <- simResult$environment
#' 
#' ## apply all  methods given in vector 'methods'
#' ## (using all data pooled for pc/LINGAM/rfci/ges -- can be changed with option 
#' ## 'onlyObservationalData=TRUE')
#' 
#' methods <- c("backShift", "LINGAM") #c("pc", "rfci", "ges") 
#' 
#' # select whether you want to run stability selection
#' stability <- FALSE
#' 
#' # arrange graphical output into a rectangular grid
#' sq <- ceiling(sqrt(length(methods)+1))
#' par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
#' 
#' ## plot and print true graph
#' cat("\n true graph is  ------  \n" )
#' print(A)
#' plotGraphEdgeAttr(A, plotStabSelec = FALSE, labels = labels, thres.point = 0, 
#'  main = "TRUE GRAPH")
#' 
#' ## loop over all methods and compute and print/plot estimate
#' for (method in methods){
#'   cat("\n result for method", method,"  ------  \n" )
#'  
#'   if(!stability){
#'     # Option 1): use this estimator as a point estimate
#'     Ahat <- getParents(X, environment, method=method, alpha=0.1, pointConf = TRUE)
#'   }else{
#'     # Option 2): use a stability selection based estimator
#'     # with expected number of false positives bounded by EV=2
#'     Ahat <- getParentsStable(X, environment, EV=2, method=method, alpha=0.1)
#'   }
#'  
#'   # print and plot estimate (point estimate thresholded if numerical estimates
#'   # are returned)
#'   print(Ahat)
#'   if(!stability)
#'     plotGraphEdgeAttr(Ahat, plotStabSelec = FALSE, labels = labels,
#'      thres.point = 0.05,
#'      main=paste("POINT ESTIMATE FOR METHOD\n", toupper(method)))
#'   else
#'     plotGraphEdgeAttr(Ahat, plotStabSelec = TRUE, labels = labels, 
#'      thres.point = 0, main = paste("STABILITY SELECTION 
#'      ESTIMATE\n FOR METHOD", toupper(method)))
#'  }
#'  
#' @keywords Causality, Graph estimations
#'  
getParents <- function(X, environment = NULL, interventions = NULL, 
                       parentsOf = 1:ncol(X),
                       method= c("ICP", "hiddenICP", "backShift", "pc", 
                                 "LINGAM", "ges", "gies", "CAM", "rfci",
                                 "regression", "bivariateANM", 
                                 "bivariateCAM")[1],  
                       alpha = 0.1, variableSelMat = NULL,
                       excludeTargetInterventions = TRUE, 
                       onlyObservationalData = FALSE, 
                       indexObservationalData = 1,
                       returnAsList=FALSE, pointConf = FALSE, 
                       setOptions = list(), directed=TRUE, verbose = FALSE){

    # check whether method is supported and dependencies are installed
    checkDependencies(method)
  
    # check validity of other input arguments
    if(is.data.frame(X)) X <- as.matrix(X)
    if(!is.matrix(X)) stop("'X' needs to be a matrix")
    if(!all(as.numeric(parentsOf) %in% (1:ncol(X)))) 
      stop("'parentsOf' needs to be a subset of 1:ncol(X)")
    if(!is.list(interventions) & !is.null(interventions)) 
      stop("'interventions' needs to be a list or NULL")
    if(length(interventions)!=nrow(X) & !is.null(interventions)) 
      stop("'interventions' needs to have as many entries as there are rows 
           in 'X' (or be 'NULL')")
    if(is.null(environment) &is.null(interventions) & method %in% 
       c("hiddenICP","ICP","backShift","gies") ) 
      stop(paste("'environment' and 'interventions' cannot 
                 both be 'NULL' for method", method))
    if(is.null(environment) & !is.null(interventions)){
        environment <- match(interventions, unique(interventions))
        if((lu <- length(unique(environment)))>50) 
          warning(paste(
                  "'environment' was set to NULL and has been created via \n 
                  '> environment <- match(interventions,unique(interventions))'\n 
                  but this results in", lu," different environments",
                  "(unique intervention combinations);\n 
                  very likely better to define a smaller number of environments 
                  using subject knowledge about the experiment by grouping 
                  various intervention targets into a single environment"))
    }
    if(length(environment) != nrow(X) & !is.null(environment)) 
      stop("'environment' needs to have the same length as there 
           are rows in 'X' (or be 'NULL')")
    if(alpha<0) stop("alpha needs to be positive")
    if(!is.null(variableSelMat)){
        if(!is.logical(variableSelMat)) 
          stop("'variableSelMat' needs to be a matrix with boolean entries")
        if(nrow(variableSelMat)!=ncol(X)) 
          stop("'variableSelMat' needs to have as many rows as there 
               are variables (columns of 'X')")
    }
    
    # find unique settings
    uniqueSettings <- unique(environment)
    
    # select observations to use
    # use observational data only or 
    # use data from different environments/interventions
    if(onlyObservationalData & length(uniqueSettings) > 1){ 
      # use only observational data
      
      if(!is.null(indexObservationalData)){
        # if indexObservationalData is given, extract observational data
        sel <- which(environment %in% indexObservationalData)
        warning(paste("Will use only environment", 
                      indexObservationalData,
                      "(= observational data?) among the", 
                      length(uniqueSettings),
                      "given distinct environments for method", 
                      method, "(", length(sel),"observations)"))
      
      }else{
        # if indexObservationalData is set to NULL, use all data points
        sel <- 1:nrow(X)
        
        warning(paste("Will use all observations for method", method, 
                      "assuming all data points are observational data (", 
                      length(sel),"observations)"))
      }
      
      X <- X[sel,]
      environment <- environment[sel]
      interventions <- interventions[sel]
    }
    
    ## gather estimated parents of each "parentsOf" node in an element of a list
    result <- list()
    for (k in 1:length(parentsOf)){
        result[[k]] <- numeric(0)
        attr(result[[k]],"parentsOf") <- parentsOf[k]
    }
    
    # run method
    switch(method,
           "ICP" = {
             result <- runICP(X, environment, interventions, parentsOf, alpha, 
                              variableSelMat, excludeTargetInterventions, 
                              pointConf, setOptions, verbose, result)
           },
           
           "hiddenICP" = {
             result <- runHiddenICP(X, environment, interventions, parentsOf, 
                                    alpha, variableSelMat, 
                                    excludeTargetInterventions, pointConf, 
                                    setOptions, verbose, result)
            },
           
            "backShift" = {
              result <- runBackShift(X, environment, parentsOf, variableSelMat, 
                                     pointConf, setOptions, verbose, result)
            },
           
            "regression" = {
              result <- runRegression(X, parentsOf, variableSelMat, pointConf, 
                                      setOptions, verbose, result)
            },
           
            "gies" = {
              result <- runGIES(X, interventions, parentsOf, variableSelMat, 
                                setOptions, directed, verbose, result)
            },
            
            "ges" = {
              result <- runGES(X, parentsOf, variableSelMat, setOptions, 
                               directed, verbose, result)
            },
            
            "pc" = {
              result <- runPC(X, parentsOf, alpha, variableSelMat, setOptions, 
                              directed, verbose, result)
            },
           
            "rfci" = {
              result <- runRFCI(X, parentsOf, alpha, variableSelMat, setOptions, 
                                directed, verbose, result)
            },
           
            "LINGAM" = {
              result <- runLINGAM(X, parentsOf, pointConf, setOptions, directed, 
                                  verbose, result)
            },
           
           "CAM" = {
              result <- runCAM(X, interventions, parentsOf, variableSelMat, 
                               setOptions, directed, verbose, result)
            },
           
            "bivariateCAM" = {
              result <- runBivariateCAM(X, parentsOf, variableSelMat, pointConf,
                                        verbose, result)
            },
           
            "bivariateANM" = {
              result <- runBivariateANM(X, parentsOf, variableSelMat, pointConf, 
                                        verbose, result)
            },
           
           {
               warning(paste("method ", method," not implemented"))
           }
           )
    
    # prepare output
    if(returnAsList){
        out <- result
    }else{
        rowind <- unlist(result)
        colind <- numeric(length(rowind))
        x <- unlist(lapply(result, function(x) attr(x,"coefficients")))
        
        if(is.null(x)) x <- 1
        
        cc <- 0
        
        for (k in 1:length(result)){
            norep <- length(result[[k]])
            if(norep>0){
                colind[ cc+(1:norep)] <- rep(k,norep)
                cc <- cc+norep
            }
        }
        
        resmat <- sparseMatrix(i=rowind,
                               j=colind,
                               x=x,
                               dims=c(ncol(X), length(parentsOf)))
        
        colnames(resmat) <- parentsOf

        out <- resmat
    }
    
    rownames(out) <- 
      if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
    colnames(out) <- rownames(out)[parentsOf]
    
    return(out)
}