########## R script: gamm-slice ##########

# Fit a generalised linear mixed model using
# slice sampling. The model is two smooth
# functions with an binary offset.

# Last changed: 06 JAN 2015

# Set flag for type of likelihood

gSlc <- function(formula,data = NULL,random = NULL,family,control = gSlc.control())
{
   
   nIter <- control$nIter
   nBurnin <- control$nBurnin
   nThin <- control$nThin
   fEPV <- control$fixedEffPriorVar  
   sdPS <- control$sdPriorScale
   method <- "stepout"

   par(mfrow=c(1,1))
   if ((nBurnin<100)|(nIter<100))
      stop("currently only working for chains longer than 100")

   if ((100*round(nBurnin/100)!=nBurnin)|(100*round(nIter/100)!=nIter))
      warning("chain lengths not multiples of 100; truncation may occur.")

   blankPlot <- function()
      plot(0,0,type="n",xlim=c(0,1),ylim=c(0,1),xaxt ="n",yaxt="n",
           xlab="",ylab="",bty="n")

   ### Present introduction window

   blankPlot()
   text(0.5,0.5,"Welcome to gammSlice!",cex=2,col="navy")

   # Using the interpret.gammSlice to obtain the parameter's value

   formula.infor <- interpret.gSlc(formula)

   # Find the number of smooth functions and linear functions

   vars.pred <- formula.infor$varlist
   vars.lin <- formula.infor$lvar
   vars.smth <- formula.infor$svar
   var.resp <- formula.infor$response   
   num.basis <- formula.infor$nbas
   
   num.smooth.funcs <- length(vars.smth)
   num.linear.funcs <- length(vars.lin)
   num.vars <- length(vars.pred)
   
  
   # Find the random factor name
   if (!is.null(random)) {
	  random.factor <- TRUE
       if (is.list(random)) {
          r.names <- names(random)
       random.vars <- c(unlist(lapply(random, function(x) all.vars(formula(x)))), 
            r.names)
       } else  {random.vars <- NULL}
   } else {random.vars <- NULL; random.factor <- FALSE}

   if (length(random.vars) > 1) {stop("There are more than 1 identity number.")}

   # Match the data and the formula

   if (is.null(data)) { 
       if (exists(var.resp) == TRUE) {
           y <- get(var.resp)
           y <- array(y,dim = c(length(y),1))
           colnames(y) <- var.resp
       } else {stop("There is no response variable.")}
       
       if (num.vars > 0) {
         x.data <- array(dim=c(length(y),0))

         for (i in 1: num.vars ) {
              if(exists(vars.pred[i]) == TRUE) {
                 temp.obj.i <- array(get(vars.pred[i]),dim = c(length(y),1))
              	 
                 if (is.character(temp.obj.i)) {
              
                     x.cate <- array(NA,dim = c(length(y),0))

                     names.catevars <- sort(unique(temp.obj.i))
                    
                     for(j in 2:length(names.catevars)) {
                        
                         x.cate <- cbind(x.cate, as.numeric(temp.obj.i == names.catevars[j]))
                         vars.lin <- c(vars.lin,paste("(",vars.pred[i],")",names.catevars[j]  ,"Vs",names.catevars[1],sep = ""))

                     }
                     x.data <- cbind(x.data, x.cate)
                     vars.lin <- vars.lin[which(vars.lin != vars.pred[i])]

                 } else {x.data <- cbind(x.data,temp.obj.i) }

               
              } else {stop("There is missing predictor.")}
         }
         colnames(x.data) <- c(vars.lin,vars.smth)
         data <- cbind(y,x.data)
       } else {data <- cbind(y)}

       if (random.factor) {
           if (exists(random.vars) == TRUE) {
              idnumcol <- array(get(random.vars), dim = c(length(y),1))
              colnames(idnumcol) <- random.vars   
              data <- cbind(idnumcol, data) 
           } else {stop("No random factor id list is provided.")}
       }
  
   } 

   # Update variables list in case of category variables
   	
   
   vars.lin  -> formula.infor$lvar
   vars.smth -> formula.infor$svar
   c(vars.lin,vars.smth) -> vars.pred -> formula.infor$varlist
   
   num.smooth.funcs <- length(vars.smth)
   num.linear.funcs <- length(vars.lin)
   num.vars <- length(vars.pred)

   if (!is.null(data)) { 
     incl.pred  <- setequal(intersect(vars.pred, colnames(data)), formula.infor$varlist)
     incl.resp  <- is.element(var.resp,colnames(data))
     if(!incl.pred) { stop("The variables names in the data set and the formula are not matched.")}
     if(!incl.resp) { stop("The response name in the data set and the formula are not matched.")}

     if (random.factor) {
        
        if (length(random.vars) == 1) {
           if (!is.element(random.vars,colnames(data))) {
               stop("The random name is not in the data label.")
           }
        } else {
           stop("The id list must be provided because there is a random intercept.")
        }
     }
   } 

   # Derive X from the data frame of X and the formula

   num.obs <- nrow(data)  
   data.original <- data
   
   if (num.vars > 0) {
      X.original <- as.matrix(data[,formula.infor$varlist])

   # Scale data set  X.original

      X.min <-  apply(X.original, 2, min)
      X.max <-  apply(X.original, 2, max)
      X.range <- X.max - X.min
      X <- t((t(X.original) - X.min)/X.range)

      data[,vars.pred] <- X
   } 

   if (random.factor) {

   ## Reorder the data set
   
   data <- data[order(data[,random.vars]),]

   ## Find the numbers of observations and repeated measurements
   
   obsid <- c(data[1:num.obs,random.vars])
   num.gps <- max(obsid)
   } else { num.gps <- nrow(data) 
   }
 
   ## Derive X from the data frame of X and the formula
   
   X <- cbind(rep(1,num.obs))
  
   if (num.vars > 0){

     X <- cbind(rep(1,num.obs), as.matrix(data[,formula.infor$varlist]))
   }

   y <- as.matrix(data[,formula.infor$response])

   ## Build the default sparse matrix Z01

   Z01 <- array(dim = c(nrow(X),0))

   if (random.factor) { 
 
      # Find the number of observations of groups 
      # which have index less than or equal gp.index 
      # for 1\leq gp.index \leq num.gps

      up.bk.id <- rep(0,num.gps)
      low.bk.id <- rep(0,num.gps)
  
      num.in.gp <- rep(0,num.gps)
      up.bk.id[1]<- num.in.gp[1] <- sum(data[,random.vars] == 1)
      low.bk.id[1] <- 1

      for (gp.index in 2 : num.gps){  
         num.in.gp[gp.index] <- sum(data[,random.vars] == gp.index)
 
         up.bk.id[gp.index] <- up.bk.id[gp.index -1]+ num.in.gp[gp.index]
         low.bk.id[gp.index] <- up.bk.id[gp.index-1] + 1
      }

      # Build the sparse matrix Z01
 
      Z01 <- matrix(1,num.in.gp[1],1)
      for (gp.index in 2: num.gps) {
         Z01 <- rbind(Z01,matrix(0,num.in.gp[gp.index],gp.index-1))
         Z01 <- cbind(Z01,c(rep(0,up.bk.id[gp.index-1]) ,rep(1,num.in.gp[gp.index])))
      }
   } else {up.bk.id <- low.bk.id <- c(1:num.gps)
          num.in.gp <- rep(1,num.gps) 
   }

   ## Build the default matrix of Zspline
   Zspline <- array(dim = c(nrow(X),0))


   nrebVal <- num.smooth.funcs 
   rbeiVal <- rep(0,nrebVal)
   
   if (num.smooth.funcs > 0) {
    
      # Set knots 
  
      num.knots <- c()

      for (i in 1: num.smooth.funcs) {
           
          if (num.basis[i] == -1) {
              spred.i  <- data[,formula.infor$svar[i]]
              nk.i <- min(20,as.integer(0.25*length(unique(spred.i))) + 1)   
              num.knots <- c(num.knots, nk.i)
          } else {
                 num.knots <- c(num.knots,num.basis[i])
          }
     }
      

      # Using information from the formula and data to
      # construct Zspline matrix and accompanying 
      # variance component block structure.

      for (i in 1: num.smooth.funcs){
 
          knots.i <- seq(0,1,length=(num.knots[i]+2))[-c(1,(num.knots[i]+2))]
 
          knots.i <- quantile(unique(data[,formula.infor$svar[i]]),knots.i)

          svd.Omega.i <- svd(abs(outer(knots.i,knots.i,"-"))^3)
     
          matrix.sqrt.Omega.i <- t(svd.Omega.i$v%*%(t(svd.Omega.i$u)*sqrt(svd.Omega.i$d)))

          data.svar.i <- c(data[,formula.infor$svar[i]])

          Z.spline.i <- t(solve(matrix.sqrt.Omega.i,t(abs(outer(data.svar.i, knots.i,"-")^3))))

          Zspline <- cbind(Zspline,Z.spline.i)
          rbeiVal[i] <- sum(num.knots[1:i])
      }   
   }

   ## Build the matrix Z

   if (random.factor) {
   Z <- cbind(Z01,Zspline)
   rbeiVal <- c(num.gps,rbeiVal + num.gps)
   nrebVal <- nrebVal + 1
   } else {
   Z <- Zspline
   }

   ## Do fitting via slice sampling.

   nc.X <- ncol(X)
   nc.Z <- ncol(Z)
   
   
   ## Form the C matrix and C^Ty vector.

   C.mat <- cbind(X,Z)
   nc.C <- ncol(C.mat)

   CTy <- rep(0,nc.C)

   for (j in 1:nc.C)
   {CTy[j] = t(C.mat[,j])%*%y
   } 

   S.u <- rep(sdPS,nrebVal)

   ## Information list

   datainfor <- list(y=y,Cmat = C.mat, CTy = CTy, 
                     numobs = num.obs, ncC = nc.C, ncX = nc.X,
                     upbkid = up.bk.id, lowbkid = low.bk.id,
                     numgps = num.gps, nrebVal = nrebVal,
                     rbeiVal = rbeiVal, Su = S.u)

   ## Set priors for parameters

   sigsq.beta <- rep(fEPV,nc.X)

   ## Set starting values
   
   lennuVal <- nc.C * nBurnin
   lenssuVal <- nrebVal * nBurnin

   nulast.Ini <- rep(0,nc.C)
   sigsq.u.Ini <- c()
   sigsq.nu.Ini <- sigsq.beta

   if(random.factor) {
      sigsq.u.Ini <-  c(sigsq.u.Ini,1)
      sigsq.nu.Ini <- c(sigsq.nu.Ini,rep(1,num.gps))
   } 
 
   if (num.smooth.funcs > 0) {
      for (i in 1: num.smooth.funcs ) {
           sigsq.u.Ini <- c(sigsq.u.Ini,1) 
           sigsq.nu.Ini <- c(sigsq.nu.Ini,rep(10,num.knots[i]))   
      }
   }  

   sum.dists.Ini <- rep(0,nc.C)
   Cnu.rest.Ini <- rep(0,num.obs)
   nucurr.Ini <- rep(1,nc.C)   
    
   pripara.Ini <- list(sumdists = sum.dists.Ini, 
         sigsqulast = sigsq.u.Ini, sigsqnu = sigsq.nu.Ini, 
         nulast = nulast.Ini, nucurr = nucurr.Ini, 
         Cnurest = Cnu.rest.Ini)

   ## Call the appropriate routine

   if (method == "NR") method <- "nr"
   if (method == "stepout") method <- "sp"

   if (random.factor == T ) {rand <- 1} 
   if (random.factor == F ) {rand <- 0}

   if (family == "binomial") {dist <- 0}
   if (family == "poisson")  {dist <- 1}

   routindex <- as.integer(2*dist + rand)

   blankPlot()
   text(0.5,0.7,"Initialisation completed.",cex=2,col="green3")
   text(0.5,0.3,"Burnin commencing.......",cex=2,col="green3")
   
   BIncSize <- round(nBurnin/100) 
   priparaB <- pripara.Ini

   for (j in 1:100) {    
        
        resultBj <- gSlcMCMC(routindex,datainfor, BIncSize,
                            priparaB)
        priparaB <- resultBj$paratoSave  
        if (j == 1) {
           blankPlot()
           text(0.5,0.7,"Percentage of burnin completed is",cex=1,col="red")
        }
        
        rect(0.2,0,0.8,0.5, col = "white",density = 1000, border = "white")

        text(0.5,0.2,j,cex=7,col="red")  
   }     

   priparaI <- resultBj$paratoSave
  
   IteIncSize <- round(nIter/100)
   rnIter <- IteIncSize * 100

   nu <- array(dim = c(nc.C,0))
   sigsq.u <- array(dim = c(nrebVal,0))
   timeTaken <- 0

   blankPlot()
   text(0.5,0.7,"Percentage of iterations completed is",cex=1,col="blue")
   for (k in 1:100) {
       
        resultIk <- gSlcMCMC(routindex,datainfor, IteIncSize + 1,
                            priparaI)
        MCMCIk <- resultIk$MCMCres
        nuIk <- array(MCMCIk$nu[,2:(IteIncSize+1)],dim = c(nc.C, IteIncSize))
        nu <- cbind(nu,nuIk)

        sigsquIk <- array(MCMCIk$sigsqu[,2:(IteIncSize+1)], dim = c(nrebVal,IteIncSize))
        sigsq.u <- cbind(sigsq.u,sigsquIk)
        
        timeTaken <- timeTaken + MCMCIk$timeTaken

        priparaI <- resultIk$paratoSave  
       
        rect(0.2,0,0.8,0.5, col = "white",density = 1000, border = "white")

        text(0.5,0.2,k,cex=7,col="blue") 
       
   }

   if (nThin > nIter) {
       stop("The number of Thin is smaller than the number of iteration.")
   }
   nu <- matrix(nu[col(nu)%% nThin == 0],nrow = nrow(nu), ncol = ncol(nu)%/% nThin)

   nu.beta <- matrix(nu[1:nc.X,], nrow = nc.X, ncol = ncol(nu))
   
   if (nrebVal > 0) {
       sigsq.u <- matrix(sigsq.u[col(sigsq.u)%% nThin == 0], nrow = nrow(sigsq.u), 
                         ncol = ncol(sigsq.u) %/% nThin)
   }
    
   ## Build the default matrix of edfMCMC,
   ## the below part is to compute edfMCMC.

   edfMCMC <- array(dim = c(num.smooth.funcs, ncol(sigsq.u)))
   nu.spline <- array(dim = c(0,ncol(nu)))

   
   if (num.smooth.funcs >0) {
      if (random.factor == F) { 
         nu.spline <- nu[(nc.X+1):(nc.X+ sum(num.knots)),]
      } else {
      nu.spline <- nu[(nc.X+num.gps+1):(nc.X+num.gps+sum(num.knots)),]
      }

   ## First set the indicator matrices for each smooth component
      
      e.linear.part <- array(0, dim = c(num.smooth.funcs,length(vars.pred)))
      colnames(e.linear.part) <- vars.pred
    
      e.smooth.part <- array(0,dim = c(num.smooth.funcs, sum(num.knots)))
   
      up.index.part <- 0
      low.index.part <- 0
      
      for (i in 1: num.smooth.funcs) {
           e.linear.part[i,num.linear.funcs + i] <- 1
           low.index.part <- up.index.part + 1
           up.index.part <- up.index.part + num.knots[i]
           e.smooth.part[i, low.index.part : up.index.part] <- rep(1,num.knots[i]) 
      }
      e.vector <- cbind(e.linear.part,e.smooth.part)

   ## Reduced forms of Cmat and nuMCMC

      nuMCMC <- nu
      if(random.factor) {
      sigsq.MCMC <- matrix(sigsq.u[2:(num.smooth.funcs+1),], nrow = num.smooth.funcs, ncol = ncol(sigsq.u))
      } else {
      sigsq.MCMC <- matrix(sigsq.u[1:num.smooth.funcs,],nrow = num.smooth.funcs, ncol = ncol(sigsq.u))  
      }
      Red.Cmat <- cbind(X,Zspline)
      
      if (!random.factor) {
          Red.nuMCMC <- nuMCMC
      } else {
        nuMCMC.linear <- matrix(nuMCMC[1:nc.X,], nrow = nc.X, ncol = ncol(nuMCMC))
        nuMCMC.smooth <- matrix(nuMCMC[(nc.X+num.gps+1):(nc.X+sum(num.knots) + num.gps), ],
                                nrow = sum(num.knots), ncol = ncol(nuMCMC) )
        Red.nuMCMC <- rbind(nuMCMC.linear,nuMCMC.smooth)
      }
      
      etaMCMC <- Red.Cmat%*% Red.nuMCMC

   ## Compute the vector w

      if(family =="binomial") wVecMCMC <- exp(etaMCMC)/((1+exp(etaMCMC))^2)
      if(family =="poisson")  wVecMCMC <- exp(etaMCMC)

   ## Compute the effective degrees of freedom MCMC samples:

      for (k in 1 : num.smooth.funcs) { 
          for (i in 1: ncol(edfMCMC)) {
              CTWC <- crossprod(Red.Cmat*wVecMCMC[,i],Red.Cmat) 
              lamVec <- rep(0,num.vars+1)

              for(j in 1: num.smooth.funcs) {
                 lamVec <- c(lamVec, rep((1/sigsq.MCMC[j,i]),num.knots[j]) )
              }
              e.vector.k  <- c(0,as.matrix(e.vector[k,]))
              edfMCMC[k,i] <- sum(diag(diag(e.vector.k)%*% solve(CTWC + diag(lamVec),CTWC)))
          }
      }
   }

     
   ## Set default value of matrix nu.linear.beta
   nu.linear.beta <- array(dim = c(0,ncol(nu.beta)))

   if (num.linear.funcs > 0) {
       nu.linear.beta <- matrix(nu.beta[2: (1+num.linear.funcs),], num.linear.funcs, ncol(nu.beta))
   }
   nu.all <- rbind(nu.beta,nu.spline)

   if (num.linear.funcs > 0) {
       nu.linear.beta <- nu.linear.beta/X.range[1:num.linear.funcs]
   }

   if (num.smooth.funcs >0) {
       summ.part <- rbind(nu.linear.beta, edfMCMC)
   } else { 
         num.knots <- NULL
         if (num.linear.funcs > 0) {
            linear.beta <- apply(nu.linear.beta,1,mean)
            beta0 <- nu[1,] - linear.beta %*% X.min 
            summ.part <- rbind(beta0,nu.linear.beta)
         } else {summ.part <-  beta0 <- nu[1,]
                X.min <- NULL
                X.max <- NULL
                X.range <- NULL
         }
   }
   
   if(random.factor) {
      sigsqU <- matrix(sigsq.u[1,],nrow = 1, ncol = ncol(sigsq.u))
      summ.part <- rbind(summ.part,sigsqU)
      u <- nu[(nc.X + 1) : (nc.X + num.gps),]
   } else {u <- NULL}
   
   # Obtain coefficients specific to this setting

   ret <- list(nu= nu.all, beta = nu.beta, u = u, sigmaSquared = sigsq.u, 
               summ = summ.part, scaleData = data,nBasis = num.knots, 
               formulaInfor = formula.infor,family = family, timeTaken = timeTaken, 
               randomFactor = random.factor, 
               Xmin = X.min, Xmax = X.max, Xrange = X.range)
 
   class(ret) <- "gSlc"
   return(ret)

}


########## End of gamm-slice ##########