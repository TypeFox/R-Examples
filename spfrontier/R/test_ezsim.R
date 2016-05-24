spfrontier.dgp <- function(){
    cat("")
    .mu <- NULL
    if (exists("mu")) .mu<-get("mu")
    .rhoU <- NULL
    if (exists("rhoU")) .rhoU<-get("rhoU")
    .rhoV <- NULL
    if (exists("rhoV")) .rhoV<-get("rhoV")
    .rhoY <- NULL
    if (exists("rhoY")) .rhoY<-get("rhoY")
    .control <- list()
    if (exists("control")) .control<-get("control")
    .beta1 <- NULL
    if (exists("beta1")) .beta1<-get("beta1")
    .beta0 <- NULL
    if (exists("beta0")) .beta0<-get("beta0")
    .beta2 <- NULL
    if (exists("beta2")) .beta2<-get("beta2")
    .sigmaU <- NULL
    if (exists("sigmaU")) .sigmaU<-get("sigmaU")
    .sigmaV <- NULL
    if (exists("sigmaV")) .sigmaV<-get("sigmaV")
    .parDef <- NULL
    if (exists("parDef")) .parDef<-get("parDef")
    .loggingLevel <- NULL
    if (exists("loggingLevel")) .loggingLevel<-get("loggingLevel")
    .inefficiency <- NULL
    if (exists("inefficiency")) .inefficiency<-get("inefficiency")
    .n <- NULL
    if (exists("n")) .n<-get("n")
    
    if (!is.null(.mu)){
        if (abs(.mu)>.sigmaU){
            cat("DGP: Truncated normal mean (mu) is higher than standard deviation, which can lead to non-skewed residuals")
        }
    }
    
    formula <- as.formula("y ~ X1 + X2")
    beta<-c(.beta1,.beta2)
    k <- length(beta)
    x1 <- runif(.n, 1,10)
    x2 <- runif(.n, 1,10)
    X <- cbind(log(x1), log(x2))
    #X <- matrix(rnorm(.n*k,0,.sigmaX),.n, k)
    
    
    W_v <- NULL
    SpW2 <- diag(.n)
    if (!is.null(.rhoV)){
        W_v <- rowStdrt(genW(.n,type="queen"))
        SpW2 <- solve(diag(.n)-.rhoV*W_v)
    }
    v <-  SpW2%*%rmvnorm(1,mean = rep(0, .n),sigma = .sigmaV^2*diag(.n))[1,]
    #if (!is.null(.rhoV)){
    #    print(coef(lm(v~Wv-1, data=data.frame(v, Wv = W_v%*%v))))
    #}
    
    W_u <- NULL
    muval <- 0
    
    if (!is.null(.mu))
        muval <- .mu
    SpW3 <- diag(.n)
    if (!is.null(.rhoU)){
        W_u <- rowStdrt(genW(.n,type="queen"))
        SpW3 <- solve(diag(.n)-.rhoU*W_u)
    }
    u <- SpW3%*%rtmvnorm(1,mean = rep(muval, .n),sigma = .sigmaU^2*diag(.n),algorithm="gibbs", lower=rep(0, .n))[1,]
    #u <- rtmvnorm(1,mean = rep(muval, .n),sigma = .sigmaU^2*SpW3*t(SpW3),algorithm="gibbs", lower=rep(0, .n))[1,]
    #u <- SpW3%*%abs(rnorm(.n, mean = muval, sd = .sigmaU))
    
    #if (!is.null(.rhoU)){
    #   print(coef(lm(u~Wu-1, data=data.frame(u, Wu = W_u%*%u))))
    #}
    sk <- skewness(v-u)
    if (sk>=0){
        cat("DGP: Skewness of generated residuals is non-negative = ",sk)
    }
    
    #plot(density(v - u))
    y <- .beta0 + X %*% beta + v - u
    
    W_y <- NULL
    if (!is.null(.rhoY)){
        W_y <- rowStdrt(genW(.n,type="queen"))
        SpW <- solve(diag(.n)-.rhoY*W_y)
        y <- SpW%*%y
    }
    dat <- data.frame(y,X)
    colnames(dat) <-c('y',paste("X", seq(k), sep = ""))
    tv <- evalFunctionOnParameterDef(.parDef,spfrontier.true.value)
    #plot(x, y, main="D (spatial u)")
    #xf<-seq(1,10,by=0.01)
    #yf<-.beta0+.beta1*log(xf)+.beta2*log(xf)^2
    #lines(xf, yf,col="red")
    if (!is.null(.control$replaceWyWv) && .control$replaceWyWv) W_y <- W_v
    if (!is.null(.control$replaceWyWu) && .control$replaceWyWu) W_y <- W_u
    if (!is.null(.control$replaceWvWu) && .control$replaceWvWu) W_v <- W_u
    
    if (!is.null(.control$ignoreWy) && .control$ignoreWy) W_y <- NULL
    if (!is.null(.control$ignoreWv) && .control$ignoreWv) W_v <- NULL
    if (!is.null(.control$ignoreWu) && .control$ignoreWu) W_u <- NULL
    m <- round(runif(1)*1000)
    print(paste("Mark = ",m))
    save(W_y, W_u, W_v, data,formula, file=paste("test",Sys.getpid(),"-",m,".rData",sep=""))
    result <- list(formula=formula, data=dat,W_y=W_y,W_v=W_v,W_u=W_u, tv=tv,y=y, X=X,
                   loggingLevel=.loggingLevel,inefficiency=.inefficiency,control=.control)
    return(result)
}

spfrontier.estimator <- function(d){
    run <- 0 
    n <- nrow(d$data)
    nam <- paste(".run",n)
    if (exists(nam, envir=.spfrontierEnv)) run<-get(nam, envir=.spfrontierEnv)
    run <- run + 1
    assign(nam, run, envir=.spfrontierEnv)
    
    message(paste("Start [pid=",Sys.getpid()," n =",n,", run =",run,"]------------------------->",rnorm(1)))
    initialValues <- NULL
    logl<-logLikelihood(formula=d$formula, data=d$data,
                        W_y = d$W_y, W_v = d$W_v,W_u = d$W_u,
                        inefficiency = d$inefficiency,
                        values=d$tv)
    message("Log-likelihood (true DGP) (pid=",Sys.getpid(),") = ", logl)
    if (!is.null(d$control) && d$control$true.initial){
        message("Using true values as initial (pid=",Sys.getpid(),"):", (logl>Infin))
        if(logl>Infin){
            initialValues <- d$tv
        }
    }
    modelEstimates <- spfrontier(d$formula,d$data,W_y=d$W_y,W_v=d$W_v,W_u=d$W_u,
                                 logging = d$loggingLevel,inefficiency=d$inefficiency,onlyCoef=T,
                                 control=list(),initialValues=initialValues)
    if (status(modelEstimates) > 0){ 
        fake = rep(1000,length(d$tv))
        names(fake) <- names(d$tv)
        out <- fake #Livehack for ezsim to exclude failure results later
    }else{
        coef <- coefficients(modelEstimates)
        out <- c(coef$beta,coef$rhoY, coef$sigmaV, coef$sigmaU, coef$rhoV, coef$rhoU, coef$mu)
        
    }
    message(paste("<-------------------------End pid=",Sys.getpid()))
    #Sys.sleep(0.1)
    flush.console()
    return(out)
}


#' @title True value for simulation
#'
#' @description
#' \code{spfrontier.true.value} returns true parameter values for a simulation process
#' 
#' @details
#' The \code{spfrontier.true.value} function should notbe used directly, it is exported for supporting \code{\link{ezsim}}
#' 
#' @rdname simulation

spfrontier.true.value <- function(){
    .mu <- NULL
    if (exists("mu")) .mu<-get("mu")
    .rhoU <- NULL
    if (exists("rhoU")) .rhoU<-get("rhoU")
    .rhoV <- NULL
    if (exists("rhoV")) .rhoV<-get("rhoV")
    .rhoY <- NULL
    if (exists("rhoY")) .rhoY<-get("rhoY")
    .control <- list()
    if (exists("control")) .control<-get("control")
    .beta1 <- NULL
    if (exists("beta1")) .beta1<-get("beta1")
    .beta0 <- NULL
    if (exists("beta0")) .beta0<-get("beta0")
    .beta2 <- NULL
    if (exists("beta2")) .beta2<-get("beta2")
    .sigmaU <- NULL
    if (exists("sigmaU")) .sigmaU<-get("sigmaU")
    .sigmaV <- NULL
    if (exists("sigmaV")) .sigmaV<-get("sigmaV")
    .parDef <- NULL
    if (exists("parDef")) .parDef<-get("parDef")
    .loggingLevel <- NULL
    if (exists("loggingLevel")) .loggingLevel<-get("loggingLevel")
    .inefficiency <- NULL
    if (exists("inefficiency")) .inefficiency<-get("inefficiency")
    .n <- NULL
    if (exists("n")) .n<-get("n")
    
    tv <- c(.beta0, .beta1, .beta2)
    tvNames <- c("beta0","beta1","beta2")
    if (is.null(.control$ignoreWy) || !.control$ignoreWy){
        if(!is.null(.rhoY)){
            tv <- c(tv, .rhoY)
            tvNames <- c(tvNames, "rhoY")
        }
    }
    tv <- c(tv, .sigmaV, .sigmaU)
    tvNames <- c(tvNames, "sigmaV","sigmaU")
    if (is.null(.control$ignoreWv) || !.control$ignoreWv){
        if(!is.null(.rhoV)){
            tv <- c(tv, .rhoV)
            tvNames <- c(tvNames, "rhoV")
        }
    }
    
    if (is.null(.control$ignoreWu) || !.control$ignoreWu){
        if(!is.null(.rhoU)){
            tv <- c(tv, .rhoU)
            tvNames <- c(tvNames, "rhoU")
        }
    }
    if(!is.null(.mu)){
        tv <- c(tv, .mu)
        tvNames <- c(tvNames, "mu")
    }
    names(tv) <- tvNames
    return(tv)
}




#' @title Spatial stochastic frontier model simulation tests
#'
#' @description
#' \code{ezsimspfrontier} tests estimators of a spatial stochastic frontier model with different parameters
#' 
#' @details
#' The \code{ezsimspfrontier} function executes multiple calls of the \code{spfrontier} estimator on a simulated data set, 
#' generated on the base of provided parameters. The resulting estimates can be analysed for biasedness, efficiency, etc.
#' 
#' 
#' @param runs a number of simulated samples 
#' @param params a set with parameters to be used in simulation.   
#' @param inefficiency sets the distribution for inefficiency error component. Possible values are 'half-normal' (for half-normal distribution) and 'truncated' (for truncated normal distribution). 
#' By default set to 'half-normal'. See references for explanations
#' @param logging an optional level of logging. Possible values are 'quiet','warn','info','debug'. 
#' By default set to quiet.
#' @param control an optional list of control parameters for simulation process. Currently the procedure supports:\cr
#' ignoreWy (TRUE/FALSE) - the spatial contiguity matrix for a dependent variable is not provided to  \code{\link{spfrontier}} estimator (but used in DGP) 
#' ignoreWv (TRUE/FALSE) - the spatial contiguity matrix for a symmetric error term is not provided to  \code{\link{spfrontier}} estimator (but used in DGP) 
#' ignoreWu (TRUE/FALSE) - the spatial contiguity matrix for a inefficiency error term is not provided to  \code{\link{spfrontier}} estimator (but used in DGP) 
#' parallel (TRUE/FALSE) - whether to use parallel computer
#' seed  - a state for random number generation in R. If NULL (default), the initial state is random. See \code{\link{set.seed}} for details.
#' auto_save  - saves intermediate results to files. See \code{\link{ezsim}} for details.
#' 
#' @keywords spatial stochastic frontier, simulation
#' @export
#' @seealso 
#' \code{\link{ezsim}}
#' @rdname simulation
#' 
#' @examples
#' params000 <- list(n=c(50, 100),beta0=5,
#'                  beta1=10,
#'                  beta2=1,
#'                  sigmaV=0.5, 
#'                  sigmaU=2.5)
#' ctrl <- list(seed=999, cores=1)
#' res000 <- ezsimspfrontier(2, params = params000,  
#'                  inefficiency = "half-normal",
#'                  logging = "info", 
#'                  control=ctrl)
#' summary(res000)

ezsimspfrontier <- function(runs,
                            params,
                            inefficiency = "half-normal",
                            logging = "info",
                            control = list()){
    start <- Sys.time()
    con <- list(auto_save=0,seed = NULL,ignoreWy=F,ignoreWv=F,ignoreWu=F,replaceWyWu=F,replaceWyWv=F, cores = detectCores()-1, outfile="spfrontier_simulations.log",true.initial=F)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    
    #rm(list = ls(envir=.SpfrontierEnv), envir=.SpfrontierEnv)
    .SpfrontierEnv <- new.env(hash = TRUE)
    
    parDef <- createParDef(selection = params, banker = list(loggingLevel=logging,inefficiency=inefficiency,control=con,parDef=createParDef(selection = params, banker=list(control=con))))
    if (!is.null(con$seed)) {
        message("Predefined random seed: ", con$seed)
        set.seed(con$seed)
    }
    cl <- NULL
    if (con$cores>1){
        message("Parallel computing is in action, number of workers: ",con$cores)
        cl <- makeCluster(con$cores, outfile=con$outfile)
        clusterApply(cl, 1:con$cores, function(x) {
            set.seed(con$seed+x)
            message("Predefined seed: ", con$seed+x)
            #devtools::load_all(".")
            require(spfrontier)
        })
    }
    ezsim_spfrontier<-ezsim(
        m             = runs,
        run           = TRUE,
        run_test      = FALSE,
        parameter_def = parDef,
        display_name = c(Intercept='beta[0]',X1='beta[1]',X2='beta[2]',beta0='beta[0]',beta1='beta[1]',beta2='beta[2]',sigmaV='sigma[v]',sigmaU='sigma[u]',sigmaX='sigma[X]',rhoV='rho[v]',rhoU='rho[u]',rhoY='rho[Y]',mu='mu'),
        dgp           = spfrontier.dgp,
        estimator     = spfrontier.estimator,
        true_value    = spfrontier.true.value,
        auto_save   = con$auto_save,
        use_core = con$cores,
        cluster = cl
    )
    if(!is.null(cl)){
        stopCluster(cl)
    }
    ezsim_spfrontier <- clearFakes(ezsim_spfrontier)
    message("Executed in ", (Sys.time()-start))
    return(ezsim_spfrontier)
}


clearFakes = function(ezsim_ob){
    parSets = length(ezsim_ob$simulation_result)
    results = list()
    for (j in 1:parSets){
        results[[j]] = data.frame()
        runs = length(ezsim_ob$simulation_result[[j]])
        for (i in 1:runs){
            if (ezsim_ob$simulation_result[[j]][[runs-i+1]][1]==1000){
                ezsim_ob$simulation_result[[j]][[runs-i+1]] = NULL
            }else{
                results[[j]] = rbind(results[[j]],ezsim_ob$simulation_result[[j]][[runs-i+1]])
            }
        }
    }
    ezsim_ob = createSimulationTable(ezsim_ob)
    ezsim_ob$results = results
    return(ezsim_ob)
}