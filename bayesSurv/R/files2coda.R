#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       files2coda.R                        ####
####                                                 ####
#### FUNCTIONS:  files2coda                          ####
#########################################################
##
## 11/04/2004: rewritten such that it returns directly mcmc object
##             with chains for all parameters
##
## 01/03/2013: minor changes
##

### ======================================
### files2coda
### ======================================
files2coda <- function(files,
                       data.frames,
                       variant = 1,
                       dir = getwd(),
                       start = 1,
                       end,
                       thin = 1,
                       header = TRUE,
                       chain)
{
  filesindir <- dir(dir)    ## character vector with available files

  if (missing(files)) misfiles <- TRUE
  else                misfiles <- FALSE

  if (missing(data.frames)) misdatfr <- TRUE
  else                      misdatfr <- FALSE
  
  ## Appropriate files, if not given
  if (misfiles && misdatfr){
    if (variant == 1)
      files <- paste(c("mixmoment", "beta", "b", "D", "Y", "r", "otherp", "u", "MHinfo", "MHinfob", "loglik"), ".sim", sep = "")
    else
      stop ("Not yet implemented variant of bayessurvreg function.")
  }
  else{
    if (!misdatfr && misfiles) files <- character(0)
  }  

  ## Indeces of iterations
  iters <- NULL
  if (sum(!is.na(match(filesindir, "iteration.sim")))){
    help <- read.table(file = paste(dir, "/", "iteration.sim", sep = ""), header = header)
    if (!dim(help)[1]) stop("Incorrect 'header' parameter.")
    iters <- as.numeric(help[,1])
    if (length(iters) == 0) iters <- NULL
  }

  ## Start and end for mcmc function of coda library
  if (!is.null(iters)){
    if (start > length(iters)) stop("start is not compatible with iteration.sim file supplied.")
    mcstart <- iters[start]
    if (missing(end))
      end <- length(iters)
    else{
      if (end > length(iters))
        stop("end is not compatible with iteration.sim file supplied.")
      else
        if (end < start) stop("start and end are not compatible.")
    }      
    mcend <- iters[end]
  }
  else{
    mcstart <- start
    if (missing(end))
      end <- numeric(0)
    else
      if (end < start) stop("start and end are not compatible.")
    mcend <- end
  }    

  ## Create a big matrix with all sampled values
  tmc <- numeric(0)

  comp <- 1
  if (length(files) > 0){
    for (i in 1:length(files)){
      if (sum(!is.na(match(filesindir, files[i])))){

        if (files[i] == "mixture.sim"){
          help <- read.table(file = paste(dir, "/", files[i], sep = ""), header = header, nrows = 1)
          if (!dim(help)[1]){
            stop(files[i], ": Incorrect 'header' parameter.", sep="")            
          }  
          kncol <- dim(help)[2]
          help <- scan(file = paste(dir, "/", files[i], sep = ""), skip = 1*header)
          wanna <- seq(1, length(help), by = kncol)
          help <- data.frame(k = help[wanna])
        }
        else{        
          help <- read.table(file = paste(dir, "/", files[i], sep = ""), header = header)
          if (!dim(help)[1]){
            ### This usually happens for those files where only a subset of iterations was stored (like Y, r, u, ...)
            ### * skip reading of this file
            next
          }  
        }          

        if (start > nrow(help)) stop("start is not compatible with data.")
        if (!length(end))
          end <- nrow(help)
        else{
         if (end > nrow(help))
            stop("end is not compatible with data.")
          else
            if (end < start) stop("start and end are not compatible.")
        }        
          
        wanna <- seq(start, end, by = thin)
        if (ncol(help) == 1){
          cname <- colnames(help)
          help <- as.data.frame(help[wanna, ])
          colnames(help) <- cname                
        }
        else{
          help <- help[wanna, ]
        }
        if (comp == 1) tmc <- help
        else           tmc <- cbind(tmc, help)        
        comp <- comp + 1
      }  
    }
  }

  if (!misdatfr){
    for (i in 1:length(data.frames)){
        if (missing(chain)) help <- get(data.frames[i])
        else                help <- get(data.frames[i])[[chain]]
        if (!dim(help)[1]) stop("Incorrect data.frame supplied.")

        if (start > nrow(help)) stop("start is not compatible with a data.frame.")
        if (!length(end))
          end <- nrow(help)
        else{
          if (end > nrow(help))
            stop("end is not compatible with a data.frame.")
          else
            if (end < start) stop("start and end are not compatible.")
        }
        wanna <- seq(start, end, by = thin)
        if (ncol(help) == 1){
          cname <- colnames(help)
          help <- as.data.frame(help[wanna, ])
          colnames(help) <- cname                
        }
        else{
          help <- help[wanna, ]
        }
        if (comp == 1) tmc <- help
        else           tmc <- cbind(tmc, help)        
        comp <- comp + 1
    }
  }

  ## Create an mcmc object
  if (length(mcend)) mc <- mcmc(tmc, start = mcstart, end = mcend, thin = thin)
  else               mc <- mcmc(tmc, start = mcstart, thin = thin)
  
  return(mc)  
}  


