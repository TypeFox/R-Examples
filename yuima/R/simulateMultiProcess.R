#
setMethod("simulate", "yuima.multimodel",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL, method="euler",
                   hurst, methodfGn="WoodChan",
                   sampling, subsampling,
                   #Initial = 0, Terminal = 1, n = 100, delta,
                   #	grid, random = FALSE, sdelta=as.numeric(NULL),
                   #	sgrid=as.numeric(NULL), interpolation="none"
                   ...){

            tmpsamp <- NULL
            if(missing(sampling)){
              tmpsamp <- setSampling(...)
              #		 tmpsamp <- setSampling(Initial = Initial, Terminal = Terminal, n = n,
              #				delta = delta, grid = grid, random = random, sdelta=sdelta,
              #				sgrid=sgrid, interpolation=interpolation)
            } else {
              tmpsamp <- sampling
            }
            tmpyuima <- setYuima(model=object, sampling=tmpsamp)

            out <- simulate(tmpyuima, nsim=nsim, seed=seed, xinit=xinit,
                            true.parameter=true.parameter,
                            space.discretized=space.discretized,
                            increment.W=increment.W, increment.L=increment.L,
                            method=method,
                            hurst=hurst,methodfGn=methodfGn, subsampling=subsampling)
            return(out)
          })

 aux.simulate.multimodel<-function(object, nsim, seed, xinit, true.parameter,
                               space.discretized, increment.W, increment.L,method,
                               hurst,methodfGn,
                               sampling, subsampling,
                               Initial, Terminal, n, delta,
                               grid, random, sdelta,
                               sgrid, interpolation){


   ##:: errors checks

   ##:1: error on yuima model
   yuima <- object

   if(missing(yuima)){
     yuima.warn("yuima object is missing.")
     return(NULL)
   }

   tmphurst<-yuima@model@hurst

   if(!missing(hurst)){
     yuima@model@hurst=hurst
   }

   if (is.na(yuima@model@hurst)){
     yuima.warn("Specify the hurst parameter.")
     return(NULL)
   }

   tmpsamp <- NULL
   if(is.null(yuima@sampling)){
     if(missing(sampling)){
       tmpsamp <- setSampling(Initial = Initial,
                              Terminal = Terminal, n = n, delta = delta,
                              grid = grid, random = random, sdelta=sdelta,
                              sgrid=sgrid, interpolation=interpolation)
     } else {
       tmpsamp <- sampling
     }
   } else {
     tmpsamp <- yuima@sampling
   }

   yuima@sampling <- tmpsamp

   sdeModel <- yuima@model
   Terminal <- yuima@sampling@Terminal[1]
   n <- yuima@sampling@n[1]
   r.size <- sdeModel@noise.number
   d.size <- sdeModel@equation.number

   ##:2: error on xinit
   if(missing(xinit)){
     xinit <- sdeModel@xinit
   } else {
     if(length(xinit) != d.size){
       if(length(xinit)==1){
         xinit <- rep(xinit, d.size)
       } else {
         yuima.warn("Dimension of xinit variables missmatch.")
         return(NULL)
       }
     }
   }

   xinit <- as.expression(xinit)  # force xinit to be an expression


   par.len <- length(sdeModel@parameter@all)

   if(missing(true.parameter) & par.len>0){
     true.parameter <- vector(par.len, mode="list")
     for(i in 1:par.len)
       true.parameter[[i]] <- 0
     names(true.parameter) <-   sdeModel@parameter@all
   }

   yuimaEnv <- new.env()

   if(par.len>0){
     for(i in 1:par.len){
       pars <- sdeModel@parameter@all[i]

       for(j in 1:length(true.parameter)){
         if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
           assign(sdeModel@parameter@all[i], true.parameter[[j]],yuimaEnv)
         }
       }
       #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
     }
   }


   if(space.discretized){
     if(r.size>1){
       warning("Space-discretized EM cannot be used for multi-dimentional models. Use standard method.")
       space.discretized <- FALSE
     }
     if(length(sdeModel@jump.coeff)){
       warning("Space-discretized EM is for only Wiener Proc. Use standard method.")
       space.discretized <- FALSE
     }
   }

   ##:: Error check for increment specified version.
   if(!missing(increment.W) & !is.null(increment.W)){
     if(space.discretized == TRUE){
       yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
       return(NULL)
     }else if(dim(increment.W)[1] != r.size){
       yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
       return(NULL)
     }else if(dim(increment.W)[2] != n){
       yuima.warn("Length of increment's column must be same as sampling@n[1].")
       return(NULL)
     }
   }

   ##:: Error check for increment specified version.
   if(!missing(increment.L) & !is.null(increment.L)){
     if(space.discretized == TRUE){
       yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
       return(NULL)
     }else if(dim(increment.L)[1] != length(yuima@model@jump.coeff[[1]]) ){ #r.size){
       yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
       return(NULL)
     }else if(dim(increment.L)[2] != n){
       yuima.warn("Length of increment's column must be same as sampling@n[1].")
       return(NULL)
     }
   }


   yuimaEnv$dL <- increment.L


   if(space.discretized){
     ##:: using Space-discretized Euler-Maruyama method
     yuima@data <- space.discretized(xinit, yuima, yuimaEnv)

     yuima@model@hurst<-tmphurst
     return(yuima)
   }


   ##:: using Euler-Maruyama method
   delta <- Terminal/n

   if(missing(increment.W) | is.null(increment.W)){

     if( sdeModel@hurst!=0.5 ){

       grid<-sampling2grid(yuima@sampling)
       isregular<-yuima@sampling@regular

       if((!isregular) || (methodfGn=="Cholesky")){
         dW<-CholeskyfGn(grid, sdeModel@hurst,r.size)
         yuima.warn("Cholesky method for simulating fGn has been used.")
       } else {
         dW<-WoodChanfGn(grid, sdeModel@hurst,r.size)
       }

     } else {

       delta<-Terminal/n
       if(!is.Poisson(sdeModel)){ # if pure CP no need to setup dW
         dW <- rnorm(n * r.size, 0, sqrt(delta))
         dW <- matrix(dW, ncol=n, nrow=r.size,byrow=TRUE)
       } else {
         dW <- matrix(0,ncol=n,nrow=1)  # maybe to be fixed
       }
     }

   } else {
     dW <- increment.W
   }

#   if(is.Poisson(sdeModel)){
#      yuima@data <- simCP(xinit, yuima, yuimaEnv)
#    } else {
#      yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
#    }

   if(is(sdeModel,"yuima.multimodel")){
     if(length(sdeModel@measure.type)==1){
        if(sdeModel@measure.type=="CP"){
          intens <- as.character(sdeModel@measure$intensity)
          dens <- as.character(sdeModel@measure$df$expr)
          dumCP <- setPoisson(intensity = intens, df = dens,
            dimension = length(sdeModel@jump.coeff[[1]]))
          dummSamp <- yuima@sampling
          samp <- setSampling(Initial = dummSamp@Initial,
                              Terminal = dummSamp@Terminal,
                              n = dummSamp@n)
         traj <- simulate(object = dumCP,
                           sampling = samp,
                           true.parameter = true.parameter)
         Incr.levy <- diff(traj@data@zoo.data[[1]])
         if(length(traj@data@zoo.data)>1){
           for(i in c(2:length(traj@data@zoo.data))){
               Incr.levy<-cbind(Incr.levy,diff(traj@data@zoo.data[[i]]))
            }
          }


        }else{
          dummSamp <- yuima@sampling
          samp <- setSampling(Initial = dummSamp@Initial,
                              Terminal = dummSamp@Terminal,
                              n = dummSamp@n)
          xinitCode <- yuima@model@xinit
          dimJumpCoeff <- length(yuima@model@jump.coeff[[1]])
          dumjumpCoeff <- matrix(as.character(diag(rep(1,dimJumpCoeff))),dimJumpCoeff,dimJumpCoeff)
          Dumsolve.variable<-paste0("MyLevyDum",c(1:dimJumpCoeff))
          LevyMod <- setMultiModel(drift=rep("0",dimJumpCoeff),
                              diffusion = NULL,
                              jump.coeff = dumjumpCoeff,
                              df = as.character(sdeModel@measure$df$expr),
                              measure.type = sdeModel@measure.type,
                              solve.variable = Dumsolve.variable)
          yuimaLevy <- setYuima(model=LevyMod, sampling = samp)
          yuimaLevy@model@dimension <- dimJumpCoeff

          traj<- simCode(xinit=xinitCode,yuima = yuimaLevy, env=yuimaEnv)

          Incr.levy <- diff(traj@data@zoo.data[[1]])
          if(length(traj@data@zoo.data)>1){
            for(i in c(2:length(traj@data@zoo.data))){
              Incr.levy<-cbind(Incr.levy,diff(traj@data@zoo.data[[i]]))
            }
          }

          #yuima.stop("code multivariate Levy will be implemented as soon as possible")
        }
     }else{
       intens <- as.character(sdeModel@measure$intensity)
       dens <- as.character(sdeModel@measure$df$expr)
       # If we consider independence between CP and the Other Levy
       # we have:
       numbLev <- length(sdeModel@measure.type)
       posCPindex <- c(1:numbLev)[sdeModel@measure.type%in%"CP"]
       CPmeasureComp <- paste0(dens,"[,c(",toString(posCPindex),")]")
       intens <- as.character(sdeModel@measure$intensity)
       dumCP <- setPoisson(intensity = intens, df = CPmeasureComp,
         dimension = length(posCPindex))

       # Simulation CP part
       dummSamp <- yuima@sampling
       samp <- setSampling(Initial = dummSamp@Initial,
                           Terminal = unique(dummSamp@Terminal),
                           n = unique(dummSamp@n))
       trajCP <- simulate(object = dumCP, sampling = samp,
          true.parameter = true.parameter)


       dimJumpCoeff <- length(yuima@model@measure.type)
       dumjumpCoeff <- matrix(as.character(diag(rep(1,dimJumpCoeff))),dimJumpCoeff,dimJumpCoeff)
       Dumsolve.variable <- paste0("MyLevyDum",c(1:dimJumpCoeff))
       dummy.measure.code <- as.character(sdeModel@measure$df$expr)
       LevyMod <- setMultiModel(drift=rep("0",dimJumpCoeff),
                                diffusion = NULL,
                                jump.coeff = dumjumpCoeff,
                                df = dummy.measure.code,
                                measure.type = "code",
                                solve.variable = Dumsolve.variable)
       yuimaLevy <- setYuima(model=LevyMod, sampling = samp)
       yuimaLevy@model@dimension <- dimJumpCoeff

       trajcode<- simCode(xinit=rep("0",length=dimJumpCoeff),
          yuima = yuimaLevy, env=yuimaEnv)

       countCP <- 0
       countcode <- 0
       if(yuima@model@measure.type[1]=="CP"){
          Incr.levy <- as.matrix(as.numeric(diff(trajCP@data@zoo.data[[1]])))
          countcode <- countcode+1
       }else{
         if(yuima@model@measure.type[1]=="code"){
            Incr.levy <- as.matrix(as.numeric(diff(trajcode@data@zoo.data[[1]])))
            countCP <- countCP+1
         }
       }

       if(length(yuima@model@measure.type)>1){
         for(i in c(2:length(yuima@model@measure.type))){
           if(yuima@model@measure.type[i]=="CP"){
             Incr.levy<-cbind(Incr.levy,as.numeric(diff(trajCP@data@zoo.data[[(i-countCP)]])))
             countcode <- countcode+1
           }else{
            if(yuima@model@measure.type[i]=="code"){
              Incr.levy <- cbind(Incr.levy,as.numeric(diff(trajcode@data@zoo.data[[i]])))
              countCP <- countCP+1
            }
           }
          }
       }


        # yuima.stop("Levy with CP and/or code")
      }
      assign("dL",t(Incr.levy),envir=yuimaEnv)
      sim <- Multi.Euler(xinit,yuima,dW,env=yuimaEnv)

   }
#   yuima@data@zoo.data <- NULL
   yuima@data@zoo.data<-as.list(numeric(length=length(sim@zoo.data)))

   for(i in 1:length(yuima@data@zoo.data)){
     yuima@data@zoo.data[[i]]<-sim@zoo.data[[i]]
     index(yuima@data@zoo.data[[i]]) <- yuima@sampling@grid[[1]]
   }## to be fixed
   yuima@data@original.data <- sim@original.data
   yuima@model@xinit <- xinit
   yuima@model@hurst <-tmphurst

   if(missing(subsampling))
     return(yuima)
   subsampling(yuima, subsampling)

 }


 simCode <- function(xinit,yuima,env){


   sdeModel<-yuima@model

   modelstate <- sdeModel@solve.variable
   modeltime <- sdeModel@time.variable
   Terminal <- yuima@sampling@Terminal[1]
   Initial <- yuima@sampling@Initial[1]
   dimension <- yuima@model@dimension
   dummy.val <- numeric(dimension)
   if(length(xinit) !=  dimension)
     xinit <- rep(xinit,  dimension)[1:dimension]
   if(length(unique(as.character(xinit)))==1 &&
      is.numeric(tryCatch(eval(xinit[1],envir=env),error=function(...) FALSE))){
     dX_dummy<-xinit[1]
     dummy.val<-eval(dX_dummy, envir=env)
     if(length(dummy.val)==1){
       dummy.val<-rep(dummy.val,dimension)
     }
     for(i in 1:length(modelstate)){
       assign(modelstate[i],dummy.val[i] ,envir=env)
     }

   } else {
     for(i in 1:dimension){
       dummy.val[i] <- eval(xinit[i], envir=env)
     }

   }

   ###    Simulation of CP using Lewis' method

   ##:: Levy
   JP <- eval(sdeModel@jump.coeff[[1]], envir=env)
   mu.size <- length(JP)
   #  print(str(JP))

   #assign(sdeModel@measure$intensity, env) ## intensity param
#    .CPintensity <- function(.t) {
#      assign(modeltime, .t, envir=env)
#      eval(sdeModel@measure$intensity, envir=env)
#    }


   dummyList<-as.list(env)

   lgth.meas<-length(yuima@model@parameter@measure)
   if(lgth.meas>1){
     for(i in c(2:lgth.meas)){
       idx.dummy<-yuima@model@parameter@measure[i]
       assign(idx.dummy,as.numeric(dummyList[idx.dummy]))
     }
   }


   # we use Lewis' acceptance/rejection method

   #if(grep("^[dexp|dnorm|dgamma|dconst]", sdeModel@measure$df$expr)){
   ##:: e.g. dnorm(z,1,1) -> rnorm(mu.size*N_sharp,1,1)
   #gsub("^d(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", sdeModel@measure$df$expr, perl=TRUE)
   #} else{
   #stop("Sorry. CP only supports dconst, dexp, dnorm and dgamma yet.")
   #}
   dumStringMeas <- toString(sdeModel@measure$df$expr)
   dumStringMeas1 <- substr(x=dumStringMeas, start=2,stop=nchar(x = dumStringMeas))
   dumStringMeas2 <- paste0("r",dumStringMeas1)
   tmpMeas2 <- strsplit(x=dumStringMeas2,split="")
   posMeas2 <- match("(" , tmpMeas2[[1]])[1]
   dumStringMeas3 <- substr(x=dumStringMeas2, start=1,stop=(posMeas2-1))
   a<-deparse(args(eval(parse(text=dumStringMeas3))))[1]
   b<-gsub("^function (.+?)","(",a)
   b1 <- substr(x=b,start =1, stop=(nchar(b)-1))
   FinalMeasRandn<-paste0(dumStringMeas3,b1)
   dummyvarMaes <- all.vars(parse(text=FinalMeasRandn))
   posDum<- match(c(sdeModel@jump.variable,sdeModel@parameter@measure),dummyvarMaes)
   if(length(posDum)+1!=length(dummyvarMaes)){
     yuima.stop("too input variables in the random number function")
   }
   deltaVar <- dummyvarMaes[-posDum]
#    ell <- optimize(f=.CPintensity, interval=c(Initial, Terminal), maximum = TRUE)$objective
#    ellMax <- ell * 1.01
   F <- suppressWarnings(parse(text=gsub("^r(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", parse(text=FinalMeasRandn), perl=TRUE)))

   F.env <- new.env(parent=env)
   N_sharp <- unique(yuima@sampling@n)
   TrueDelta <- unique(yuima@sampling@delta)
   assign(deltaVar, TrueDelta, envir=F.env)
   assign("mu.size", mu.size, envir=F.env)
   assign("N_sharp", N_sharp, envir=F.env)
   randJ <- eval(F, envir=F.env)  ## this expression is evaluated in the F.env


   randJ <- rbind(dummy.val, randJ)
   RANDLevy <- apply(randJ,2,cumsum)
   tsX <- zoo(x=RANDLevy)
   yuimaData <- setYuima(data=setData(tsX, delta=TrueDelta))
   #yuimaData <- subsampling(yuimaData, sampling=yuima@sampling)
   return(yuimaData)
 }



 Multi.Euler<-function(xinit,yuima,dW,env){

   sdeModel<-yuima@model

   modelstate <- sdeModel@solve.variable
   modeltime <- sdeModel@time.variable
   V0 <- sdeModel@drift
   V <- sdeModel@diffusion
   r.size <- sdeModel@noise.number
   d.size <- sdeModel@equation.number
   Terminal <- yuima@sampling@Terminal[1]
   n <- yuima@sampling@n[1]
   dL <- env$dL

   #	dX <- xinit

   # 06/11 xinit is an expression: the structure is equal to that of V0
   if(length(unique(as.character(xinit)))==1 &&
      is.numeric(tryCatch(eval(xinit[1],env),error=function(...) FALSE))){
     dX_dummy<-xinit[1]
     dummy.val<-eval(dX_dummy, env)
     if(length(dummy.val)==1){dummy.val<-rep(dummy.val,length(xinit))}
     for(i in 1:length(modelstate)){
       assign(modelstate[i],dummy.val[i] ,env)
     }
     dX<-vector(mode="numeric",length(dX_dummy))

     for(i in 1:length(xinit)){
       dX[i] <- dummy.val[i]
     }
   }else{
     dX_dummy <- xinit
     if(length(modelstate)==length(dX_dummy)){
       for(i in 1:length(modelstate)) {
         if(is.numeric(tryCatch(eval(dX_dummy[i],env),error=function(...) FALSE))){
           assign(modelstate[i], eval(dX_dummy[i], env),env)
         }else{
           assign(modelstate[i], 0, env)
         }
       }
     }else{
       yuima.warn("the number of model states do not match the number of initial conditions")
       return(NULL)
     }

     # 06/11 we need a initial variable for X_0
     dX<-vector(mode="numeric",length(dX_dummy))

     for(i in 1:length(dX_dummy)){
       dX[i] <- eval(dX_dummy[i], env)
     }
   }
   ##:: set time step
   delta <- Terminal/n

   ##:: check if DRIFT and/or DIFFUSION has values
   has.drift <- sum(as.character(sdeModel@drift) != "(0)")
   var.in.diff <- is.logical(any(match(unlist(lapply(sdeModel@diffusion, all.vars)), sdeModel@state.variable)))
   #print(is.Poisson(sdeModel))

   ##:: function to calculate coefficients of dW(including drift term)
   ##:: common used in Wiener and CP
   p.b <- function(t, X=numeric(d.size)){
     ##:: assign names of variables
     for(i in 1:length(modelstate)){
       assign(modelstate[i], X[i], env)
     }
     assign(modeltime, t, env)
     ##:: solve diffusion term
     if(has.drift){
       tmp <- matrix(0, d.size, r.size+1)
       for(i in 1:d.size){
         tmp[i,1] <- eval(V0[i], env)
         for(j in 1:r.size){
           tmp[i,j+1] <- eval(V[[i]][j],env)
         }
       }
     } else {  ##:: no drift term (faster)
       tmp <- matrix(0, d.size, r.size)
       if(!is.Poisson(sdeModel)){ # we do not need to evaluate diffusion
         for(i in 1:d.size){
           for(j in 1:r.size){
             tmp[i,j] <- eval(V[[i]][j],env)
           } # for j
         } # foh i
       } # !is.Poisson
     } # else
     return(tmp)
   }

   X_mat <- matrix(0, d.size, n+1)
   X_mat[,1] <- dX

   if(has.drift){  ##:: consider drift term to be one of the diffusion term(dW=1)
     dW <- rbind( rep(1, n)*delta , dW)
   }


   if(!length(sdeModel@measure.type)){ ##:: Wiener Proc
     ##:: using Euler-Maruyama method

     if(var.in.diff & (!is.Poisson(sdeModel))){  ##:: diffusions have state variables and it is not Poisson
       ##:: calcurate difference eq.
       for( i in 1:n){
         dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i]
         X_mat[,i+1] <- dX
       }
     }else{  ##:: diffusions have no state variables (not use p.b(). faster)
       sde.tics <- seq(0, Terminal, length=(n+1))
       sde.tics <- sde.tics[2:length(sde.tics)]

       X_mat[, 1] <- dX

       ##:: assign names of variables
       for(i in 1:length(modelstate)){
         assign(modelstate[i], dX[i])
       }
       assign(modeltime, sde.tics)
       t.size <- length(sde.tics)

       ##:: solve diffusion term
       if(has.drift){
         pbdata <- matrix(0, d.size*(r.size+1), t.size)
         for(i in 1:d.size){
           pbdata[(i-1)*(r.size+1)+1, ] <- eval(V0[i], env)
           for(j in 1:r.size){
             pbdata[(i-1)*(r.size+1)+j+1, ] <- eval(V[[i]][j], env)
           }
         }
         dim(pbdata)<-(c(r.size+1, d.size*t.size))
       }else{
         pbdata <- matrix(0, d.size*r.size, t.size)
         if(!is.Poisson(sdeModel)){
           for(i in 1:d.size){
             for(j in 1:r.size){
               pbdata[(i-1)*r.size+j, ] <- eval(V[[i]][j], env)
             } # for j
           } # for i
         } # !is.Poisson
         dim(pbdata)<-(c(r.size, d.size*t.size))
       } # else

       pbdata <- t(pbdata)

       ##:: calcurate difference eq.
       for( i in 1:n){
         if(!is.Poisson(sdeModel))
           dX <- dX + pbdata[((i-1)*d.size+1):(i*d.size), ] %*% dW[, i]
         X_mat[, i+1] <- dX
       }
     }
     tsX <- ts(data=t(X_mat), deltat=delta , start=0)

   }else{ ##:: Levy
     JP <- sdeModel@jump.coeff
     mu.size <- length(JP)
     # cat("\n Levy\n")
     ##:: function to solve c(x,z)
     p.b.j <- function(t, X=numeric(d.size)){
       for(i in 1:length(modelstate)){
         assign(modelstate[i], X[i], env)
       }
       assign(modeltime, t, env)
       #      tmp <- numeric(d.size)
       j.size <- length(JP[[1]])
       tmp <- matrix(0, mu.size, j.size)
       # cat("\n inside\n")
       #print(dim(tmp))
       for(i in 1:mu.size){
         for(j in 1:j.size){
           tmp[i,j] <- eval(JP[[i]][j],env)
         }
         #        tmp[i] <-  eval(JP[i], env)
       }
       return(tmp)
     }
     #  print(ls(env))

     ### WHY I AM DOING THIS?
     #    tmp <- matrix(0, d.size, r.size)
     #
     #for(i in 1:d.size){
     #        for(j in 1:r.size){
     #            cat("\n here\n")
     #            tmp[i,j] <- eval(V[[i]][j],env)
     #        } # for j
     #    }
     ###

#      if(sdeModel@measure.type == "CP" ){ ##:: Compound-Poisson type
#
#        ##:: delete 2010/09/13 for simulate func bug fix by s.h
#        ## eta0 <- eval(sdeModel@measure$intensity)
#
#        ##:: add 2010/09/13 for simulate func bug fix by s.h
#        eta0 <- eval(sdeModel@measure$intensity, env) ## intensity param
#
#        ##:: get lambda from nu()
#        tmp.expr <- function(my.x){
#          assign(sdeModel@jump.variable,my.x)
#          return(eval(sdeModel@measure$df$expr))
#        }
#        #lambda <- integrate(sdeModel@measure$df$func, 0, Inf)$value * eta0
#        #lambda <- integrate(tmp.expr, 0, Inf)$value * eta0 ##bug:2013/10/28
#
#        dummyList<-as.list(env)
#        #   print(str(dummyList))
#        #print(str(idx.dummy))
#        lgth.meas<-length(yuima@model@parameter@measure)
#        if(lgth.meas>1){
#          for(i in c(2:lgth.meas)){
#            idx.dummy<-yuima@model@parameter@measure[i]
#            #print(i)
#            #print(yuima@model@parameter@measure[i])
#            assign(idx.dummy,as.numeric(dummyList[idx.dummy]))
#          }
#        }
#
#
#        #lambda <- integrate(tmp.expr, -Inf, Inf)$value * eta0
#
#        ##:: lambda = nu() (p6)
#        N_sharp <- rpois(1,Terminal*eta0)	##:: Po(Ne)
#        if(N_sharp == 0){
#          JAMP <- FALSE
#        }else{
#          JAMP <- TRUE
#          Uj <- sort( runif(N_sharp, 0, Terminal) )
#          ij <- NULL
#          for(i in 1:length(Uj)){
#            Min <- min(which(c(1:n)*delta > Uj[i]))
#            ij <- c(ij, Min)
#          }
#        }
#
#        ##:: make expression to create iid rand J
#        if(grep("^[dexp|dnorm|dgamma|dconst]", sdeModel@measure$df$expr)){
#          ##:: e.g. dnorm(z,1,1) -> rnorm(mu.size*N_sharp,1,1)
#          F <- suppressWarnings(parse(text=gsub("^d(.+?)\\(.+?,", "r\\1(mu.size*N_sharp,", sdeModel@measure$df$expr, perl=TRUE)))
#        }else{
#          stop("Sorry. CP only supports dconst, dexp, dnorm and dgamma yet.")
#        }
#
#        ##:: delete 2010/09/13 for simulate func bug fix by s.h
#        ## randJ <- eval(F)  ## this expression is evaluated locally not in the yuimaEnv
#
#        ##:: add 2010/09/13 for simulate func bug fix by s.h
#        F.env <- new.env(parent=env)
#        assign("mu.size", mu.size, envir=F.env)
#        assign("N_sharp", N_sharp, envir=F.env)
#
#        randJ <- eval(F, F.env)  ## this expression is evaluated in the F.env
#
#        j <- 1
#        for(i in 1:n){
#          if(JAMP==FALSE || sum(i==ij)==0){
#            Pi <- 0
#          }else{
#            if(is.null(dL)){
#              J <- eta0*randJ[j]/lambda
#              j <- j+1
#              ##cat(paste(J,"\n"))
#              ##Pi <- zeta(dX, J)
#              assign(sdeModel@jump.variable, J, env)
#
#              if(sdeModel@J.flag){
#                J <- 1
#              }
#
#              Pi <- p.b.j(t=i*delta,X=dX) * J
#            }else{# we add this part since we allow the user to specify the increment of CP LM 05/02/2015
#              Pi <- p.b.j(t=i*delta,X=dX) %*% dL[, i]
#            }
#            ##Pi <- p.b.j(t=i*delta, X=dX)
#          }
#          dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] + Pi
#          X_mat[, i+1] <- dX
#        }
#        tsX <- ts(data=t(X_mat), deltat=delta, start=0)
#        ##::end CP
#      }else if(sdeModel@measure.type=="code"){  ##:: code type
#        ##:: Jump terms
#        code <- suppressWarnings(sub("^(.+?)\\(.+", "\\1", sdeModel@measure$df$expr, perl=TRUE))
#        args <- unlist(strsplit(suppressWarnings(sub("^.+?\\((.+)\\)", "\\1", sdeModel@measure$df$expr, perl=TRUE)), ","))
#        #print(args)
#        dZ <- switch(code,
#                     rNIG=paste("rNIG(n, ", args[2], ", ", args[3], ", ", args[4], "*delta, ", args[5], "*delta, ", args[6],")"),
#                     rIG=paste("rIG(n,", args[2], "*delta, ", args[3], ")"),
#                     rgamma=paste("rgamma(n, ", args[2], "*delta, ", args[3], ")"),
#                     rbgamma=paste("rbgamma(n, ", args[2], "*delta, ", args[3], ", ", args[4], "*delta, ", args[5], ")"),
#                     ##                   rngamma=paste("rngamma(n, ", args[2], "*delta, ", args[3], ", ", args[4], ", ", args[5], "*delta, ", args[6], ")"),
#                     rngamma=paste("rngamma(n, ", args[2], "*delta, ", args[3], ", ", args[4], ", ", args[5], "*delta,", args[6],")"),
#                     ##                   rstable=paste("rstable(n, ", args[2], ", ", args[3], ", ", args[4], ", ", args[5], ", ", args[6], ")")
#                     rstable=paste("rstable(n, ", args[2], ", ", args[3], ", ", args[4], "*delta^(1/",args[2],"), ", args[5], "*delta)")
#        )
#        dummyList<-as.list(env)
#        #print(str(dummyList))
#        lgth.meas<-length(yuima@model@parameter@measure)
#        #print(lgth.meas)
#        if(lgth.meas!=0){
#          for(i in c(1:lgth.meas)){
#            #print(i)
#            #print(yuima@model@parameter@measure[i])
#            idx.dummy<-yuima@model@parameter@measure[i]
#            #print(str(idx.dummy))
#            assign(idx.dummy,dummyList[[idx.dummy]])
#            #print(str(idx.dummy))
#            #print(str(dummyList[[idx.dummy]]))
#            #print(get(idx.dummy))
#          }
#        }

#        if(is.null(dZ)){  ##:: "otherwise"
#          cat(paste("Code \"", code, "\" not supported yet.\n", sep=""))
#          return(NULL)
#        }
       # if(!is.null(dL))
         dZ <- dL
       # else
         # dZ <- eval(parse(text=dZ))
       ##:: calcurate difference eq.
       #print(str(dZ))
#        if(is.null(dim(dZ)))
#          dZ <- matrix(dZ,nrow=1)
       # print(dim(dZ))
       #  print(str(sdeModel@jump.variable))
       for(i in 1:n){
         assign(sdeModel@jump.variable, dZ[,i], env)

         if(sdeModel@J.flag){
           dZ[,i] <- 1
         }
         #           cat("\np.b.j call\n")
         tmp.j <- p.b.j(t=i*delta, X=dX)
         #print(str(tmp.j))
         #cat("\np.b.j cback and dZ\n")
         # print(str(dZ[,i]))
         # print(sum(dim(tmp.j)))

         #print(str(tmp.j))
         #print(str(p.b(t = i * delta, X = dX) %*% dW[, i]))
         dX <- dX + p.b(t=i*delta, X=dX) %*% dW[, i] +tmp.j %*% dZ[,i]
         X_mat[, i+1] <- dX
       }
       tsX <- ts(data=t(X_mat), deltat=delta, start=0)
       ##::end code
#      }else{
#        cat(paste("Type \"", sdeModel@measure.type, "\" not supported yet.\n", sep=""))
#        return(NULL)
#      }
   }##::end levy
   yuimaData <- setData(original.data=tsX)
   return(yuimaData)


 }
