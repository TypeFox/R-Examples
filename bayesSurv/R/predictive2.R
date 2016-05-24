#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2004)                              ####
####                                                 ####
#### FILE:       predictive2.R                       ####
####                                                 ####
#### FUNCTIONS:  predictive2                         ####
#########################################################

### ======================================
### predictive2
### ======================================
## 23/01/2005: start working on it
## 07/02/2005: G-spline random intercept allowed
##
predictive2 <- function(
     formula,
     random,
     obs.dim,
     time0,
     data = parent.frame(),
     grid,
     na.action = na.fail,
     Gspline,
     quantile = c(0, 0.025, 0.5, 0.975, 1),
     skip = 0,
     by = 1,
     last.iter,
     nwrite,
     only.aver = TRUE,
     predict = list(density=FALSE, Surv=TRUE, hazard=FALSE, cum.hazard=FALSE),
     dir = getwd(),
     extens = "",
     extens.random = "_b",
     version = 0
)                        
{
   thispackage = "bayesSurv"
   #thispackage = NULL
   
   transform = function(t){log(t)}
   dtransform = function(t){1/t}

   ## Extract all the design information from the function call
   ## ==========================================================
   m <- match.call(expand.dots = FALSE)
   des <- bayessurvreg.design(m, formula, random, data, transform, dtransform)
   dims <- c(des$n, des$ncluster, des$nwithin)

   ## Check some other parameters and create vectors for C++ function call
   ## =====================================================================
   if (missing(obs.dim)) obs.dim <- numeric(0)   
   if (missing(time0)) time0 <- numeric(0)
   if (missing(Gspline)) stop("Gspline parameter must be given")
   control   <- predictive2.control(predict, only.aver, quantile, obs.dim, time0, Gspline, des$n)
   predict   <- control$predict
   only.aver <- control$only.aver
   quantile  <- control$quantile
   obs.dim   <- control$obs.dim
   time0     <- control$time0
   Gspline   <- control$Gspline

   if (length(quantile) == 0) only.aver <- TRUE
   n.quantile <- ifelse(only.aver, 0, length(quantile))
   
   ## Onset/event? (Needed only by version = 32)
   ## ==========================================
   if (extens == "") Onset <- 1
   else{
     if (extens == "_2") Onset <- 0
     else                warning("Somewhat strange value of the argument extens")
   }  
   
   ## C++ parameters to hold 'betaGamma' object
   ## ==========================================
   obj.Beta <- bayessurvreg2.priorBeta(prior.beta=list(mean.prior=rep(0, des$nX), var.prior=rep(1, des$nX)),
                                       init=list(beta=rep(0, des$nX)),
                                       design=des)
   ## Check whether needed files are available
   ## * further,  whether at least first row has correct number of elements
   ## * and determine the MC sample size
   ## =======================================================================
   filesindir <- dir(dir)
   if (!length(filesindir)) stop("Empty directory with simulated values?")
   
   if (obj.Beta$parmI["nbeta"]){
     if (sum(!is.na(match(filesindir, paste("beta", extens, ".sim", sep=""))))){
       bbeta <- read.table(paste(dir, "/beta", extens, ".sim", sep = ""), nrows = 1)
       if (length(bbeta) != obj.Beta$parmI["nbeta"])
         stop(paste("Your formula indicates that there are ", obj.Beta$parmI["nbeta"], " regression parameters in the model however the file beta.sim contains ", length(bbeta), " columns", sep=""))
     }
     else
       stop("File beta.sim not found.")
   }

   if (sum(!is.na(match(filesindir, paste("mweight", extens, ".sim", sep=""))))){
     mix <- read.table(paste(dir, "/mweight", extens, ".sim", sep = ""), nrows = 1)
     k.max <- length(mix)
     if (k.max != Gspline$total.length) stop("Different total_length of the G-spline indicated by the file mweight.sim and Gspline parameter of this function")
   }
   else
     stop("File with simulated values of mixture weights not found.")

   if (sum(!is.na(match(filesindir, paste("mmean", extens, ".sim", sep=""))))){
     mix <- read.table(paste(dir, "/mmean", extens, ".sim", sep = ""), nrows = 1)
     kmax2 <- length(mix)/Gspline$dim
     if (k.max != kmax2) stop("Different total_length of the G-spline indicated by files mweight.sim and mmean.sim.")
   }
   else
     stop("File with simulated values of mixture means indeces  not found.")
   
   if (sum(!is.na(match(filesindir, paste("gspline", extens, ".sim", sep=""))))){
     mix <- read.table(paste(dir, "/gspline", extens, ".sim", sep = ""), nrows = 1)     
     lmix <- length(mix)
     if (lmix != 5*Gspline$dim) stop(paste("You indicate that dimension is ", Gspline$dim, " so that file gspline.sim must have ", 5*Gspline$dim, " columns", sep=""))
   }
   else{
     stop("File with simulated values of gamma/sigma/delta/intercept/scale not found.")
   }     

   if (sum(!is.na(match(filesindir, paste("mixmoment", extens, ".sim", sep=""))))){
     mix <- read.table(paste(dir, "/mixmoment", extens, ".sim", sep = ""), header = TRUE)
     if (missing(last.iter)) M <- dim(mix)[1]
     else{
       M <- last.iter
       if (last.iter > dim(mix)[1]) M <- dim(mix)[1]
       if (last.iter <= 0)          M <- dim(mix)[1]
     }      
   }
   else{
     stop("File mixmoment.sim not found.")
   }   

   ## C++ parameters to hold possible 'b' object
   ## ===========================================
   if (des$nrandom){
     if (version == 3){
       if (sum(!is.na(match(filesindir, paste("mweight", extens.random, ".sim", sep=""))))){
         mix.b <- read.table(paste(dir, "/mweight", extens.random, ".sim", sep = ""), nrows = 1)
         k.max.b <- length(mix.b)
       }
       else{
         stop("File with simulated values of random intercept mixture weights not found.")
       }
       if (sum(!is.na(match(filesindir, paste("mmean", extens.random, ".sim", sep=""))))){
         mix.b <- read.table(paste(dir, "/mmean", extens.random, ".sim", sep = ""), nrows = 1)
         kmax2.b <- length(mix.b)
         if (k.max.b != kmax2.b) stop("Different total_length of the G-spline indicated by files mweight_b.sim and mmean_b.sim.")
       }
       else{
         stop("File with simulated values of random intercept mixture means indeces not found.")
       }
       if (sum(!is.na(match(filesindir, paste("gspline", extens.random, ".sim", sep=""))))){
         mix.b <- read.table(paste(dir, "/gspline", extens.random, ".sim", sep = ""), nrows = 1)     
         lmix.b <- length(mix.b)
         if (lmix.b != 5) stop(paste("You indicate that dimension is ", 1, " so that file gspline_b.sim must have ", 5, " columns", sep=""))
       }
       else{
         stop("File with simulated values of random intercept gamma/sigma/delta/intercept/scale not found.")
       }
       if (sum(!is.na(match(filesindir, paste("mixmoment", extens.random, ".sim", sep=""))))){
         mix.b<- read.table(paste(dir, "/mixmoment", extens.random, ".sim", sep = ""), header = TRUE)
         if (missing(last.iter)) M.b <- dim(mix.b)[1]
         else{
           M.b <- last.iter
           if (last.iter > dim(mix)[1]) M.b <- dim(mix.b)[1]
           if (last.iter <= 0)          M.b <- dim(mix.b)[1]
         }      
       }
       else{
         stop("File mixmoment_b.sim not found.")
       }
       K.b <- (k.max.b - 1)/2

       b.GsplI <- c(1, k.max.b)
       prior.b <- list(specification=2, K=K.b, izero=0, neighbor.system="uniCAR", order=3, equal.lambda=TRUE,
                       prior.lambda="fixed", prior.gamma="fixed", prior.intercept="fixed", prior.sigma="fixed", prior.scale="fixed",
                       c4delta=1.5)
       mcmc.par.b <- list(type.update.a.b="slice", k.overrelax.a.b=1, k.overrelax.sigma.b=1, k.overrelax.scale.b=1)
       init.b <- list(b=rep(0, des$ncluster), lambda.b=1, sigma.b=0.3, gamma.b=0, scale.b=1, intercept.b=0, a.b=numeric(0))
       obj.b <- bayessurvreg3.priorb(prior.b=prior.b, init=init.b, design=des, mcmc.par=mcmc.par.b)
     }    ## end of if (version == 3)     
     else{
       if (version == 32){
         lD <- 3
         if (sum(!is.na(match(filesindir, paste("D", ".sim", sep=""))))){
           DD <- read.table(paste(dir, "/D", ".sim", sep = ""), nrows = 1)
           if (length(DD) != lD + 1){
             stop(paste("File D.sim should contain ", lD+1, " columns, however it has ", length(DD), " columns", sep=""))
           }             
         }
         else{
           stop("File D.sim not found.")
         }                  
       }
       else{
         lD <- 0.5*(des$nrandom*(1 + des$nrandom))
         if (sum(!is.na(match(filesindir, paste("D", extens, ".sim", sep=""))))){
           DD <- read.table(paste(dir, "/D", extens, ".sim", sep = ""), nrows = 1)
           if (length(DD) != lD + 1){
             stop(paste("File D", extens, ".sim should contain ", lD+1, " columns, however it has ", length(DD), " columns", sep=""))
           }             
         }
         else{
           stop(paste("File D", extens, ".sim not found.", sep=""))
         }         
       }  
       
       b.GsplI <- c(1, 0)       
       if (version == 32){
         obj.b <- list()
         obj.b$bparmI <- c(0, 1, des$ncluster, des$nwithin)
         obj.b$bparmD <- rep(0, des$ncluster)
       }
       else{
         prior.b <- list(prior.D = "inv.wishart", df.D = des$nrandom+2,
                         scale.D = diag(des$nrandom)[lower.tri(diag(des$nrandom), diag = TRUE)])
         obj.b <- bayessurvreg2.priorb(prior.b = prior.b, init = list(beta=rep(0, des$nX)), design=des)            
       }         
     } 
   }    ## end of if(des$nrandom)
   else{
     b.GsplI <- c(1, 0)     
     prior.b <- list()
     obj.b <- bayessurvreg2.priorb(prior.b = prior.b, init = list(beta=rep(0, des$nX)), design=des)
     M.b <- M
   }    ## end of else(des$nrandom)
   names(b.GsplI) <- c("dim", "total.length")
   

   ## Grids
   ## =====
   tt0 <- time0[obs.dim]
   if (missing(grid)) stop("grid must be supplied")   
   if (is.list(grid)){
     stop("Different grid values for each combination of covariates not yet implemented")
     simple.grid <- FALSE
     if (length(grid) != des$n) stop("Incorrect 'grid' parameter supplied.")
     ngrid <- sapply(grid, length)     
     gridall <- unlist(grid)
   }
   else{
     simple.grid <- TRUE
     ngrid <- rep(length(grid), des$n)
     gridall <- rep(grid, des$n)
   }
   sum.ngrid <- sum(ngrid)
   
   gridall <- gridall - rep(tt0, ngrid)   
   if (!is.numeric(ngrid)) stop("Incorrect 'grid' parameter supplied.")
   if (sum(ngrid <= 0)) stop("Incorrect 'grid' parameter supplied)")
   if (sum(gridall <= 0.0)) stop("All grid values must be higher than corresponding time0.")
   if (sum(is.na(gridall))) stop("NA's appeared in grid")   
   log.gridall <- log(gridall)
   if (sum(is.na(log.gridall))) stop("NA's appeared in log(grid)")    

   if (version == 3){
     if (M != M.b) stop("Different sample sized indicated by mixmoment.sim and mixmoment_b.sim files")
   }     

   if (missing(skip)) skip <- 0
   else{
     if (skip > M) stop("You ask to skip more rows from the file than available.")
     if (skip < 0) skip <- 0
   }
   if (missing(by)) by <- 1
   else{
     if (by <= 0) by <- 1
   }    

   ## Needed space for various quantities
   ## ===================================
   M.now <- 1 + (M-skip-1) %/% by
   lvalue      <- ifelse(only.aver, sum.ngrid, sum.ngrid * M.now)
   lquantvalue <- ifelse(only.aver, 1, sum.ngrid * n.quantile)

   if (predict$density){
     lv.density <- lvalue
     lvquant.density <- lquantvalue
   }
   else{
     lv.density <- if (predict$hazard) sum.ngrid else 1
     lvquant.density <-  1
   }

   if (predict$Surv){
     lv.Surv <- lvalue
     lvquant.Surv <- lquantvalue
   }
   else{
     lv.Surv <- if (predict$hazard || predict$cum.hazard) sum.ngrid else 1
     lvquant.Surv <-  1
   }

   if (predict$hazard){
     lv.hazard <- lvalue
     lvquant.hazard <- lquantvalue     
   }
   else{
     lv.hazard <- 1
     lvquant.hazard <- 1
   }     

   if (predict$cum.hazard){
     lv.cum.hazard <- lvalue
     lvquant.cum.hazard <- lquantvalue     
   }
   else{
     lv.cum.hazard <- 1
     lvquant.cum.hazard <- 1
   }     

   if (missing(nwrite)) nwrite <- M
   
   ## Simulation
   ## ===========
   mcmc <- .C("predictive_GS", av.dens          = double(ifelse(predict$density, sum.ngrid, 1)),
                               av.Surv          = double(ifelse(predict$Surv, sum.ngrid, 1)),
                               av.hazard        = double(ifelse(predict$hazard, sum.ngrid, 1)),
                               av.cum.hazard    = double(ifelse(predict$cum.hazard, sum.ngrid, 1)),
                               val.dens         = double(lv.density),
                               val.Surv         = double(lv.Surv),
                               val.hazard       = double(lv.hazard),
                               val.cum.hazard   = double(lv.cum.hazard),
                               quant.dens       = double(lvquant.density),
                               quant.Surv       = double(lvquant.Surv),
                               quant.hazard     = double(lvquant.hazard),
                               quant.cum.hazard = double(lvquant.cum.hazard),              
                               dims             = as.integer(dims),
                               X                = as.double(if(des$nX) t(des$X) else 0),              
                               obs.dims         = as.integer(obs.dim - 1),
                               M.now            = as.integer(M.now),
                               dir              = as.character(dir),
                               extens           = as.character(extens),
                               extens.random    = as.character(extens.random),
                               GsplI            = as.integer(c(Gspline$dim, Gspline$total.length, Gspline$K)),
                               objBetaI         = as.integer(obj.Beta$parmI),
                               objBetaD         = as.double(obj.Beta$parmD),
                               objbI            = as.integer(obj.b$bparmI),
                               objbD            = as.double(obj.b$bparmD),
                               b.GsplI          = as.integer(b.GsplI),
                               grid             = as.double(gridall),
                               log.grid         = as.double(log.gridall),
                               ngrid            = as.integer(ngrid),
                               probs            = as.double(if(only.aver) 0.5 else quantile),
                               nquant           = as.integer(n.quantile),
                               only.aver        = as.integer(only.aver),
                               predict          = as.integer(c(predict$density, predict$Surv, predict$hazard, predict$cum.hazard)),
                               M                = as.integer(M),
                               skip             = as.integer(skip),
                               by               = as.integer(by),
                               nwrite           = as.integer(nwrite),
                               version          = as.integer(version),
                               Onset            = as.integer(Onset),
                               err              = integer(1),
              PACKAGE = thispackage)

   if (mcmc$err) stop("No results produced, something is wrong.")   

   toreturn <- list()
   if (simple.grid){
     toreturn$grid <- grid              ## this starts again with time0 and we give it back to the user
     if (predict$density)     toreturn$density    <- matrix(mcmc$av.dens, nrow=des$n, byrow=TRUE)
     if (predict$Surv)        toreturn$Surv       <- matrix(mcmc$av.Surv, nrow=des$n, byrow=TRUE)
     if (predict$hazard)      toreturn$hazard     <- matrix(mcmc$av.hazard, nrow=des$n, byrow=TRUE)
     if (predict$cum.hazard)  toreturn$cum.hazard <- matrix(mcmc$av.cum.hazard, nrow=des$n, byrow=TRUE)
     if (!only.aver){
       l.obs <- ngrid[1] * n.quantile
       if (predict$density){
         toreturn$quant.density <- list()
         for (i in 1:des$n){
           toreturn$quant.density[[i]] <- matrix(mcmc$quant.dens[((i-1)*l.obs+1):(i*l.obs)], byrow=TRUE, ncol=ngrid[1])
           rownames(toreturn$quant.density[[i]]) <- paste(quantile*100, "%", sep="")
         }
       }
       if (predict$Surv){
         toreturn$quant.Surv <- list()
         for (i in 1:des$n){
           toreturn$quant.Surv[[i]] <- matrix(mcmc$quant.Surv[((i-1)*l.obs+1):(i*l.obs)], byrow=TRUE, ncol=ngrid[1])
           rownames(toreturn$quant.Surv[[i]]) <- paste(quantile*100, "%", sep="")          
         }                    
       }
       if (predict$hazard){
         toreturn$quant.hazard <- list()
         for (i in 1:des$n){
           toreturn$quant.hazard[[i]] <- matrix(mcmc$quant.hazard[((i-1)*l.obs+1):(i*l.obs)], byrow=TRUE, ncol=ngrid[1])
           rownames(toreturn$quant.hazard[[i]]) <- paste(quantile*100, "%", sep="")          
         }                             
       }
       if (predict$cum.hazard){
         toreturn$quant.cum.hazard <- list()
         for (i in 1:des$n){
           toreturn$quant.cum.hazard[[i]] <- matrix(mcmc$quant.cum.hazard[((i-1)*l.obs+1):(i*l.obs)], byrow=TRUE, ncol=ngrid[1])
           rownames(toreturn$quant.cum.hazard[[i]]) <- paste(quantile*100, "%", sep="")          
         }                                      
       }       
     }
   }     
   else{
     stop("Different grid values for each combination of covariates not yet implemented")
   }

   attr(toreturn, "sample.size") <- mcmc$M.now
   return(toreturn)  
}
