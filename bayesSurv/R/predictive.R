####################################################
#### AUTHOR:     Arnost Komarek                 ####
####             (2004)                         ####
####                                            ####
#### FILE:       predictive                     ####
####                                            ####
#### FUNCTIONS:  predictive                     ####
####################################################

### ======================================
### predictive
### ======================================
predictive <- function(
     formula,
     random,
     time0 = 0,
     data = parent.frame(),
     grid,
     type = "mixture",
     subset,
     na.action = na.fail,
     quantile = c(0, 0.025, 0.5, 0.975, 1),                       
     skip = 0,
     by = 1,
     last.iter,
     nwrite,
     only.aver = FALSE,
     predict = list(Et=TRUE, t=FALSE, Surv=TRUE, hazard=FALSE, cum.hazard=FALSE),
     store = list(Et=TRUE, t = FALSE, Surv = FALSE, hazard = FALSE, cum.hazard=FALSE),
     Eb0.depend.mix = FALSE,
     dir = getwd(),
     toler.chol = 1e-10,
     toler.qr = 1e-10)
{
   thispackage = "bayesSurv"
   #thispackage = NULL
  
   transform = function(t){log(t)}
   dtransform = function(t){1/t}

   typeError<- pmatch(type, table = c("mixture", "spline", "polya.tree"), nomatch = 0) - 1
   if (typeError == -1 || typeError >= 2) stop("Unknown or not yet implemented error type.")
   
   control <- predictive.control(predict, store, only.aver, quantile)
   predict <- control$predict
   store <- control$store
   only.aver <- control$only.aver
   quantile <- control$quantile
   predictShch <- predict$Surv || predict$hazard || predict$cum.hazard

   ## Starting time for the survival model
   ## ====================================
   if (missing(time0)) time0 <- 0
   if (time0 < 0) stop("time0 must be non-negative.")

   ## Extract all the design information from the function call
   ## ==========================================================
   m <- match.call(expand.dots = FALSE)
   des <- bayessurvreg.design(m, formula, random, data, transform, dtransform)

   ## Check whether needed files are available
   ## and whether at least first row has correct number of elements
   ## ==============================================================
   filesindir <- dir(dir)    ## character vector with available files
   if (!length(filesindir)) stop("Empty directory with simulated values?")   
   if (sum(!is.na(match(filesindir, "mweight.sim")))){
     mix <- read.table(paste(dir, "/mweight.sim", sep = ""), nrows = 1)     
     kmax <- length(mix)
   }
   else
     stop("File with simulated values of mixture weights not found.")

   if (sum(!is.na(match(filesindir, "mmean.sim")))){
     mix <- read.table(paste(dir, "/mmean.sim", sep = ""), nrows = 1)
     kmax2 <- length(mix)
     if (kmax != kmax2) stop("Different kmax indicated by files mweight.sim and mmean.sim.")
   }
   else
     stop("File with simulated values of mixture means not found.")
   
   if (sum(!is.na(match(filesindir, "mvariance.sim")))){
     mix <- read.table(paste(dir, "/mvariance.sim", sep = ""), nrows = 1)     
     kmax2 <- length(mix)
     if (kmax != kmax2) stop("Different kmax indicated by files mweight.sim and mvariance.sim.")
   }
   else
     stop("File with simulated values of mixture variances not found.")
      
   nbeta <- des$nfixed + des$nrandom - des$randomInt
   if (nbeta){
     if (sum(!is.na(match(filesindir, "beta.sim")))){
       beta <- scan(file = paste(dir, "/beta.sim", sep = ""), nlines = 1, skip = 1)
       if (length(beta) != nbeta) stop("Incorrect 'beta.sim' file supplied.")
     }
     else
       stop("File with simulated values of regression parameters not found.")
   }

   if (des$nrandom){
     nD <- 0.5*(des$nrandom*(1 + des$nrandom))
     if (sum(!is.na(match(filesindir, "D.sim")))){
       D <- scan(file = paste(dir, "/D.sim", sep = ""), nlines = 1, skip = 1)     
       if (length(D) != 1 + nD) stop("Incorrect 'D.sim' file supplied.")
     }
     else
       stop("File with simulated values of a covariance matrix of random effects not found.")
   }

   ## nsimul, skip, by
   ## ================   
   if (sum(!is.na(match(filesindir, "mixmoment.sim")))){
     mix <- read.table(paste(dir, "/mixmoment.sim", sep = ""), header = TRUE)
     if (missing(last.iter)) M <- dim(mix)[1]
     else{
       M <- last.iter
       if (last.iter > dim(mix)[1]) M <- dim(mix)[1]
       if (last.iter <= 0)          M <- dim(mix)[1]
     }      
   }
   else
     stop("File mixmoment.sim not found.")     
    
   if (missing(skip)) skip <- 0
   else{
     if (skip > M) stop("You ask to skip more rows from the file than available.")
     if (is.na(skip) || skip < 0) skip <- 0
   }
   if (missing(by)) by <- 1
   else{
     if (is.na(by) || by <= 0) by <- 1
   }    

   lvalue <- 1 + (M - skip - 1) %/% by
   if (missing(nwrite)) nwrite <- lvalue
   if (nwrite > lvalue) nwrite <- lvalue
   
   ## Grids
   ## ======  
   if (predictShch){
     if (is.list(grid)){
       if (length(grid) != des$n) stop("Incorrect 'grid' parameter supplied.")
       ngrid <- sapply(grid, length)
       gridall <- unlist(grid)
     }
     else{
       ngrid <- rep(length(grid), des$n)
       gridall <- rep(grid, des$n)
     }
     if (!is.numeric(ngrid)) stop("Incorrect 'grid' parameter supplied.")
     if (sum(ngrid <= 0)) stop("Incorrect 'grid' parameter supplied.")
     if (sum(gridall <= 0.0)) stop("All grid values must be positive.")
   }
   cumngrid <- c(0, cumsum(ngrid))
  
   ## Create files to store predictive values
   ## =======================================
   filesindir <- dir(dir)    ## character vector with files in dir
      
   write.headers <- function(filename, predict0, dir, filesindir, n, gridall, cumngrid, obs = TRUE, label, namesX)
   {
     lf <- nchar(filename)
     remove <- substring(filesindir, first = 1, last = lf)
     remove <- pmatch(remove, table = filename, nomatch = 0, duplicates.ok = TRUE)
     if (any(remove)) file.remove(paste(dir, "/", filesindir[remove == 1], sep = ""))
     if (predict0){
       if (obs){
         for (i in 1:n){
           sink(paste(dir, "/", filename, i, ".sim", sep = ""), append = FALSE)
           cat(paste(gridall[(cumngrid[i] + 1):cumngrid[i+1]], sep = ""), "\n", sep = "  "); sink()
         }
       }
       else{
           sink(paste(dir, "/", filename, ".sim", sep = ""), append = FALSE)
           cat(paste(label, namesX, sep = ""), "\n", sep = "     "); sink()
       } 
     }          
   }

   write.headers("predET", store$Et, dir, filesindir, des$n, gridall, cumngrid, FALSE, "ET", des$rnamesX)
   write.headers("quantET", predict$Et, dir, filesindir, des$n, gridall, cumngrid, FALSE, "ET", des$rnamesX)
   write.headers("predT", store$t, dir, filesindir, des$n, gridall, cumngrid, FALSE, "T", des$rnamesX)
   write.headers("quantT", predict$t, dir, filesindir, des$n, gridall, cumngrid, FALSE, "T", des$rnamesX)
   write.headers("predS", store$Surv, dir, filesindir, des$n, gridall, cumngrid)
   write.headers("predhazard", store$hazard, dir, filesindir, des$n, gridall, cumngrid)
   write.headers("predcumhazard", store$cum.hazard, dir, filesindir, des$n, gridall, cumngrid)
   write.headers("quantS", predict$Surv, dir, filesindir, des$n, gridall, cumngrid)
   write.headers("quanthazard", predict$hazard, dir, filesindir, des$n, gridall, cumngrid)
   write.headers("quantcumhazard", predict$cum.hazard, dir, filesindir, des$n, gridall, cumngrid)


   ## Sample
   ## ========
   dims <- c(des$n, des$ncluster, des$nwithin, des$nY, des$nX, des$nfixed, des$nrandom, 1*des$randomInt, nwrite)
   dims2 <- c(length(quantile) ,ngrid)
   predictV <- c(predict$Et, predict$t, predict$Surv, predict$hazard, predict$cum.hazard)
   storeV <- c(store$Et, store$t, store$Surv, store$hazard, store$cum.hazard)
   nsimul <- c(M, 0, nwrite)       ## (nthin is ignored)
   tolers <- c(toler.chol, toler.qr)
   prior.pari <- c(kmax, 0, 1*Eb0.depend.mix)
   
   cat("Simulation started on                       ", date(), "\n", sep = "")
   fit <- .C("predictive", as.integer(typeError),
                           as.character(dir),
                           as.integer(dims),
                           as.integer(dims2),
                           X = as.double(des$X),
                           indb = as.integer(des$indb),
                           quant = as.double(quantile),
                           grid = as.double(gridall),
                           prior.pari = as.integer(prior.pari),
                           prior.pard = as.double(time0),
                           nsimul = as.integer(nsimul),
                           skip = as.integer(skip),
                           by = as.integer(by),
                           only.aver = as.integer(only.aver),
                           predict = as.integer(predictV),
                           store = as.integer(storeV),
                           tolers = as.double(tolers),
                           err = integer(1),
             PACKAGE = thispackage)

   if (fit$err != 0) warning ("Something went wrong during the simulation.")
   cat("Simulation finished on                      ", date(), "\n", sep = "")

   return(fit$err)      
}
