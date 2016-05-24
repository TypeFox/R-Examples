#' This function captures output files from MLwiN for estimation in
#' WinBUGS/OpenBUGS.
#' 
#' This function allows R to call WinBUGS using the output files from MLwiN.
#' This function uses functionalities in the \code{\link[rbugs]{rbugs}}
#' package.
#' 
#' @param D A vector specifying the type of distribution used in the model.
#' @param levID A character (vector) specifying the level ID(s).
#' @param datafile A file name where the BUGS data file will be saved in
#' .txt format.
#' @param initfiles A list of file names where the BUGS initial values will
#' be saved in .txt format.
#' @param modelfile A file name where the BUGS model will be saved in .txt
#' format.
#' @param bugEst A file name where the estimates from BUGS will be stored in
#' .txt format.
#' @param fact A list of objects used to specify factor analysis. See `Details'
#' below.
#' @param addmore A vector of strings specifying additional coefficients to be
#' monitored.
#' @param n.chains The number of chains to be monitored.
#' @param n.iter The number of iterations for each chain
#' @param n.burnin The length of burn-in for each chain
#' @param n.thin Thinning rate
#' @param debug A logical value indicating whether (\code{TRUE}) or not
#' (\code{FALSE}; the default) to close the BUGS window after completion of the
#' model run
#' @param bugs The full name (including the path) of the BUGS executable
#' @param bugsWorkingDir A directory where all the intermediate files are to be
#' stored; defaults to \code{tempdir()}.
#' @param OpenBugs If \code{TRUE}, OpenBUGS is used, if \code{FALSE} (the
#' default) WinBUGS is used.
#' @param cleanBugsWorkingDir If \code{TRUE}, the generated files will be
#' removed from the \code{bugsWorkingDir}; defaults to \code{FALSE}.
#' @param seed An integer specifying the random seed.
#'
#' @details
#' A list of objects to specify factor analysis, as used in the
#' argument \code{fact}:
#' \itemize{
#' \item \code{nfact}: specifies the number of factors;
#' \item \code{lev.fact}: Specifies the level/classification for the random part of
#' the factor for each factor;
#' \item \code{nfactcor}: specifies the number of
#' correlated factors;
#' \item \code{factcor}: a vector specifying the correlated
#' factors: the first element corresponds to the first factor number, the
#' second to the second factor number, the third element corresponds to the
#' starting value for the covariance and the fourth element to whether this
#' covariance is constrained
#' (\code{1}) or not (\code{0}). If more than one pair of factors is correlated,
#' then repeat this sequence for each pair.
#' \item \code{loading}: a matrix specifying the
#' starting values for the factor loadings and the starting value of the factor
#' variance. Each row corresponds to a factor.
#' \item \code{constr}: a matrix
#' specifying indicators of whether the factor loadings and the factor variance
#' are constrained (\code{1}) or not (\code{0}).
#' }
#'
#' @return Returns an \code{\link[coda]{mcmc}} object.
#'
#' @note This function is derived from \code{\link[rbugs]{rbugs}} (written by
#' Jun Yan and Marcos Prates).
#'
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso \code{\link{runMLwiN}},\code{\link[rbugs]{rbugs}}
#' @export
mlwin2bugs <- function(D,levID, datafile, initfiles, modelfile, bugEst, fact, addmore, n.chains, n.iter, n.burnin, n.thin, debug=FALSE, bugs,
                       bugsWorkingDir=tempdir(), OpenBugs = FALSE, cleanBugsWorkingDir = FALSE, seed = NULL){
  
  rbugs2 <- function (data.file, inits.files, paramSet, model, bugEst, fact, n.chains = 1, n.iter = 2000,
                      n.burnin = floor(n.iter/2), n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/1000)), dic = FALSE, debug = FALSE,
                      bugs = system("which OpenBUGS", TRUE), bugsWorkingDir, OpenBugs = TRUE,
                      cleanBugsWorkingDir = FALSE, genFilesOnly = FALSE, verbose = FALSE,
                      seed = NULL) {
    Windows = FALSE
    os.type <- .Platform$OS.type
    if (length(bugsWorkingDir) == 0)
      stop("It is required to specify the bugsWorkingDir")
    if (os.type == "windows") {
      if (!file.exists(bugs))
        stop(paste("BUGS executable", bugs, "does not exists."))
      Windows = TRUE
    }
    else if (os.type == "unix") {
      if (length(bugs) == 0)
        bugs <- system("which OpenBUGS", TRUE)
      if (length(bugs) == 0)
        stop(paste("BUGS executable", bugs, "does not exists."))
    }
    else warning("This function has not been tested on mac-os.")
    
    if (is.null(bugsWorkingDir)) {
      bugsWorkingDir <- tempfile("bugsWorkingDir")
      if (!file.exists(bugsWorkingDir))
        dir.create(bugsWorkingDir)
      on.exit(if (cleanBugsWorkingDir) unlink(bugsWorkingDir,
                                              TRUE))
    }
    workingDir <- bugsWorkingDir
    
    if (!file.exists(model)) stop("Model file doesn't exist.")
    model.file <- file.path(workingDir, "model.txt")
    file.copy(model, model.file, overwrite = TRUE)
    
    
    script.file <- paste(workingDir, "script.txt", sep = "/")
    rbugs::genBugsScript(paramSet, n.chains, n.iter, n.burnin, n.thin,
                  dic, model.file, data.file, inits.files, bugsWorkingDir,
                  script.file, debug, OpenBugs, Windows, seed)
    
    trLbr <- function(unix) {
      lines <- readLines(unix)
      writeLines(lines, unix, sep="\r\n")
    }
    if (OpenBugs) {
      trLbr(model.file)
      trLbr(data.file)
      for (i in inits.files) trLbr(i)
    }
    if (genFilesOnly) {
      cat("Files are generated in", workingDir, "\n")
      return(TRUE)
    }
    rbugs::runBugs(bugs, script.file, n.chains, workingDir, OpenBugs,
            Windows, verbose)
    all <- rbugs::getBugsOutput(n.chains, workingDir, OpenBugs)  ##Modified by Marcos
    all$n.iter = n.iter                                   ##Modified by Marcos
    all$n.burnin = n.burnin                               ##Modified by Marcos
    all$n.thin = n.thin                                   ##Modified by Marcos
    
    #### functions to get the output
    getCodaFileNames <- function(n.chains, workingDir, OpenBugs) {
      CODA <- if (OpenBugs) "codaCODA" else "coda"
      INDEX <- if (OpenBugs) "index" else "Index"
      CHAIN <- if (OpenBugs) "chain" else NULL
      coda  <- file.path(workingDir, CODA)
      codaFiles <- paste(coda, CHAIN, 1:n.chains, ".txt", sep="")
      codaIndexFile <- paste(coda, INDEX, ".txt", sep="")
      list(codaFiles=codaFiles, codaIndexFile=codaIndexFile)
    }       
    if(cleanBugsWorkingDir) { ##Modified by Marcos
      fnames <- getCodaFileNames(n.chains, workingDir, OpenBugs) ##Modified by Marcos
      coda.files <- c(fnames$codaIndexFile, fnames$codaFiles) ##Modified by Marcos
      log.file <- file.path(workingDir, "log.txt") ##Modified by Marcos
      file.remove(c(model.file, log.file, data.file, ##Modified by Marcos
                    inits.files[1], script.file, coda.files)) ##Modified by Marcos
    }
    
    all
  }
  environment(rbugs2)<-environment(rbugs)
  
  nlev= length(levID)
  if(nlev==1){
    parameters=c("beta","va0","sigma")
  }else{
    if (D[1]=="Multivariate Normal"||D[1]=="Multinomial"||D[1]=="Mixed"){
      ux=sapply(3:nlev, function(x) paste("u",x,sep=""))
      sigma2=sapply(2:nlev, function(x) paste("sigma2.u",x,sep=""))
      if (!is.null(fact)){
        nfact=fact$nfact
        loadname=NULL
        sigma2.fact.name=factname=rep(0,nfact)
        for (i in 1:nfact){
          loadname=c(loadname,sapply(1:(ncol(fact$loading)-1), function(x) paste("load",i,".",x,sep="")))
          factname[i]=paste("fact",i,sep="")
          sigma2.fact.name[i]=paste("sigma2.fact",i,sep="")
        }
        parameters=c("beta",ux, sigma2, loadname,factname,sigma2.fact.name)
      }else{
        parameters=c("beta",ux, sigma2)
      }
    }else{
      if(D[1]=="t"){
        ux=sapply(2:nlev, function(x) paste("u",x,sep=""))
        sigma2=sapply(2:nlev, function(x) paste("sigma2.u",x,sep=""))
        parameters=c("beta",ux, sigma2, "va0", "sigma","sigma2","df")
      }else{
        ux=sapply(2:nlev, function(x) paste("u",x,sep=""))
        sigma2=sapply(2:nlev, function(x) paste("sigma2.u",x,sep=""))
        parameters=c("beta",ux, sigma2, "va0", "sigma","sigma2")
      }
    }
  }
  if (!is.null(addmore)) parameters=c(parameters,addmore)
  chains.bugs=rbugs2(data.file = datafile, inits.files = initfiles,
                     paramSet=parameters, model=modelfile, bugEst=bugEst, n.chains = n.chains, n.iter = n.iter, n.burnin=n.burnin, n.thin=n.thin, debug=debug, bugs=bugs,
                     bugsWorkingDir=bugsWorkingDir, OpenBugs=OpenBugs, cleanBugsWorkingDir=cleanBugsWorkingDir, seed = seed)
  chains.bugs.mcmc=rbugs::rbugs2coda(chains.bugs,burnin=1,n.thin)
  chains.bugs.mcmc
}
