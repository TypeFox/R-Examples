## We have splitted the simulate function into blocks to allow for future
## methods to be added. S.M.I. & A.B.
## Interface to simulate() changed to match the S3 generic function in the
## package stats
## added an environment to let R find the proper values defined in the main
## body of the function, which in turn calls different simulation methods
## All new simulation methods should look into the yuimaEnv for local variables
## when they need to "eval" R expressions

##:: function simulate
##:: solves SDE and returns result
subsampling <- function(x,y) return(x)

setGeneric("simulate",
           function(object, nsim=1, seed=NULL, xinit, true.parameter, space.discretized=FALSE,
                    increment.W=NULL, increment.L=NULL, method="euler", hurst, methodfGn="WoodChan",
                    sampling=sampling, subsampling=subsampling, ...
                    #		Initial = 0, Terminal = 1, n = 100, delta,
                    #		grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL),
                    #		sgrid=as.numeric(NULL), interpolation="none"
           )
             standardGeneric("simulate")
)


setMethod("simulate", "yuima.model",
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


# We rewrite the code of simulate method. We build a new internal function aux.simulate containing
# the old code. This function starts directly if the model is an object of yuima.model-class
# or yuima.carma-class. If the model is an object of class yuima.cogarch, simulate method runs
# the internal function aux.simulateCogarch that generates first a path of the underlying levy and then
# uses this path to construct the trajectories of the Cogarch model by calling the aux.simulate function.

setMethod("simulate", "yuima",
          function(object, nsim=1, seed=NULL, xinit, true.parameter,
                   space.discretized=FALSE, increment.W=NULL, increment.L=NULL,
                   method="euler",
                   hurst,methodfGn="WoodChan",
                   sampling, subsampling,
                   Initial = 0, Terminal = 1, n = 100, delta,
                   grid = as.numeric(NULL), random = FALSE, sdelta=as.numeric(NULL),
                   sgrid=as.numeric(NULL), interpolation="none"){

            if(is(object@model,"yuima.cogarch")){
              res<-aux.simulateCogarch(object, nsim, seed, xinit, true.parameter,
                           space.discretized, increment.W, increment.L, method,
                           hurst,methodfGn,
                           sampling, subsampling,
                           Initial, Terminal, n, delta,
                           grid, random, sdelta,
                           sgrid, interpolation)
            }else{
              if(is(object@model,"yuima.multimodel")){
                res <- aux.simulate.multimodel(object, nsim, seed, xinit, true.parameter,
                                               space.discretized, increment.W, increment.L,method,
                                               hurst,methodfGn,
                                               sampling, subsampling,
                                               Initial, Terminal, n, delta,
                                               grid, random, sdelta,
                                               sgrid, interpolation)
              }else{
                res<-aux.simulate(object, nsim, seed, xinit, true.parameter,
                                     space.discretized, increment.W, increment.L,method,
                                     hurst,methodfGn,
                                     sampling, subsampling,
                                     Initial, Terminal, n, delta,
                                     grid, random, sdelta,
                                     sgrid, interpolation)
              }
            }
            return(res)

#             ##:: errors checks
#
#             ##:1: error on yuima model
#             yuima <- object
#
#             if(missing(yuima)){
#               yuima.warn("yuima object is missing.")
#               return(NULL)
#             }
#
#             tmphurst<-yuima@model@hurst
#
#             if(!missing(hurst)){
#               yuima@model@hurst=hurst
#             }
#
#             if (is.na(yuima@model@hurst)){
#               yuima.warn("Specify the hurst parameter.")
#               return(NULL)
#             }
#
#             tmpsamp <- NULL
#             if(is.null(yuima@sampling)){
#               if(missing(sampling)){
#                 tmpsamp <- setSampling(Initial = Initial,
#                                        Terminal = Terminal, n = n, delta = delta,
#                                        grid = grid, random = random, sdelta=sdelta,
#                                        sgrid=sgrid, interpolation=interpolation)
#               } else {
#                 tmpsamp <- sampling
#               }
#             } else {
#               tmpsamp <- yuima@sampling
#             }
#
#             yuima@sampling <- tmpsamp
#
#             sdeModel <- yuima@model
#             Terminal <- yuima@sampling@Terminal[1]
#             n <- yuima@sampling@n[1]
#             r.size <- sdeModel@noise.number
#             d.size <- sdeModel@equation.number
#
#             ##:2: error on xinit
#             if(missing(xinit)){
#               xinit <- sdeModel@xinit
#             } else {
#               if(length(xinit) != d.size){
#                if(length(xinit)==1){
#                  xinit <- rep(xinit, d.size)
#                } else {
#                 yuima.warn("Dimension of xinit variables missmatch.")
#                 return(NULL)
#                }
#               }
#             }
#
#             xinit <- as.expression(xinit)  # force xinit to be an expression
#
#
#             par.len <- length(sdeModel@parameter@all)
#
#             if(missing(true.parameter) & par.len>0){
#               true.parameter <- vector(par.len, mode="list")
#               for(i in 1:par.len)
#                 true.parameter[[i]] <- 0
#               names(true.parameter) <-   sdeModel@parameter@all
#             }
#
#             yuimaEnv <- new.env()
#
#             if(par.len>0){
#               for(i in 1:par.len){
#                 pars <- sdeModel@parameter@all[i]
#
#                 for(j in 1:length(true.parameter)){
#                   if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
#                     assign(sdeModel@parameter@all[i], true.parameter[[j]],yuimaEnv)
#                   }
#                 }
#                 #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
#               }
#             }
#
#
#             if(space.discretized){
#               if(r.size>1){
#                 warning("Space-discretized EM cannot be used for multi-dimentional models. Use standard method.")
#                 space.discretized <- FALSE
#               }
#               if(length(sdeModel@jump.coeff)){
#                 warning("Space-discretized EM is for only Wiener Proc. Use standard method.")
#                 space.discretized <- FALSE
#               }
#             }
#
#             ##:: Error check for increment specified version.
#             if(!missing(increment.W) & !is.null(increment.W)){
#               if(space.discretized == TRUE){
#                 yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
#                 return(NULL)
#               }else if(dim(increment.W)[1] != r.size){
#                 yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
#                 return(NULL)
#               }else if(dim(increment.W)[2] != n){
#                 yuima.warn("Length of increment's column must be same as sampling@n[1].")
#                 return(NULL)
#               }
#             }
#
#             ##:: Error check for increment specified version.
#             if(!missing(increment.L) & !is.null(increment.L)){
#               if(space.discretized == TRUE){
#                 yuima.warn("Parameter increment must be invalid if space.discretized=TRUE.")
#                 return(NULL)
#               }else if(dim(increment.L)[1] != length(yuima@model@jump.coeff[[1]]) ){ #r.size){
#                 yuima.warn("Length of increment's row must be same as yuima@model@noise.number.")
#                 return(NULL)
#               }else if(dim(increment.L)[2] != n){
#                 yuima.warn("Length of increment's column must be same as sampling@n[1].")
#                 return(NULL)
#               }
#             }
#
#
#             yuimaEnv$dL <- increment.L
#
#
#             if(space.discretized){
#               ##:: using Space-discretized Euler-Maruyama method
#               yuima@data <- space.discretized(xinit, yuima, yuimaEnv)
#
#               yuima@model@hurst<-tmphurst
#               return(yuima)
#             }
#
#
#             ##:: using Euler-Maruyama method
#             delta <- Terminal/n
#
#             if(missing(increment.W) | is.null(increment.W)){
#
#               if( sdeModel@hurst!=0.5 ){
#
#                 grid<-sampling2grid(yuima@sampling)
#                 isregular<-yuima@sampling@regular
#
#                 if((!isregular) || (methodfGn=="Cholesky")){
#                   dW<-CholeskyfGn(grid, sdeModel@hurst,r.size)
#                   yuima.warn("Cholesky method for simulating fGn has been used.")
#                 } else {
#                   dW<-WoodChanfGn(grid, sdeModel@hurst,r.size)
#                 }
#
#               } else {
#
#                 delta<-Terminal/n
#                 if(!is.Poisson(sdeModel)){ # if pure CP no need to setup dW
#                  dW <- rnorm(n * r.size, 0, sqrt(delta))
#                  dW <- matrix(dW, ncol=n, nrow=r.size,byrow=TRUE)
#                 } else {
#                     dW <- matrix(0,ncol=n,nrow=1)  # maybe to be fixed
#                 }
#               }
#
#             } else {
#               dW <- increment.W
#             }
#
#             if(is.Poisson(sdeModel)){
#                 yuima@data <- simCP(xinit, yuima, yuimaEnv)
#             } else {
#                 yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
#             }
#
#             for(i in 1:length(yuima@data@zoo.data))
#               index(yuima@data@zoo.data[[i]]) <- yuima@sampling@grid[[1]]  ## to be fixed
#             yuima@model@xinit <- xinit
#             yuima@model@hurst <-tmphurst
#
#             if(missing(subsampling))
#               return(yuima)
#             subsampling(yuima, subsampling)
#
           }
          )

aux.simulate<-function(object, nsim, seed, xinit, true.parameter,
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

  if(is.Poisson(sdeModel)){
    yuima@data <- simCP(xinit, yuima, yuimaEnv)
  } else {
    yuima@data <- euler(xinit, yuima, dW, yuimaEnv)
  }

  for(i in 1:length(yuima@data@zoo.data))
    index(yuima@data@zoo.data[[i]]) <- yuima@sampling@grid[[1]]  ## to be fixed
  yuima@model@xinit <- xinit
  yuima@model@hurst <-tmphurst

  if(missing(subsampling))
    return(yuima)
  subsampling(yuima, subsampling)

}

aux.simulateCogarch<-function(object, nsim, seed, xinit, true.parameter,
                              space.discretized, increment.W, increment.L, method,
                              hurst,methodfGn,
                              sampling, subsampling,
                              Initial, Terminal, n, delta,
                              grid, random, sdelta,
                              sgrid, interpolation){
  yuimaCogarch<-object
  model<-yuimaCogarch@model
  info<-model@info
  samp <- yuimaCogarch@sampling
  if(model@measure.type=="CP" && !is.null(increment.L)){
    method="euler"
  }

  if(method=="euler"||(method=="mixed" && model@measure.type=="code")){
     if(length(increment.L)==0){
          aux.Noise<-setModel(drift="0",
                              diffusion="0",
                              jump.coeff="1",
                              measure=info@measure,
                              measure.type=info@measure.type)


    #    aux.samp<-setSampling(Initial = samp@Initial, Terminal = samp@Terminal[1], n = samp@n[1], delta = samp@delta,
    #                         grid=samp@grid, random = samp@random, sdelta=samp@sdelta,
    #                         sgrid=samp@sgrid, interpolation=samp@interpolation )

        aux.samp<-setSampling(Initial = samp@Initial,
                              Terminal = samp@Terminal[1],
                              n = samp@n[1])
        auxModel<-setYuima(model=aux.Noise, sampling= aux.samp)

      if(length(model@parameter@measure)==0){
        aux.incr2<-aux.simulate(object=auxModel, nsim=nsim, seed=seed,
                               space.discretized=space.discretized, increment.W=increment.W,
                               increment.L=increment.L,
                               hurst=0.5,methodfGn=methodfGn)
      }else{
        aux.incr2<-aux.simulate(object=auxModel, nsim=nsim, seed=seed,
                                true.parameter = true.parameter[model@parameter@measure],
                                space.discretized=space.discretized, increment.W=increment.W,
                                increment.L=increment.L,
                                hurst=0.5,methodfGn=methodfGn)
      }
      increment<-diff(as.numeric(get.zoo.data(aux.incr2)[[1]]))
     } else{increment<-increment.L}
    # Using the simulated increment for generating the quadratic variation
    # As first step we compute it in a crude way. A more fine approach is based on
    # the mpv function.
    quadratVariation <- increment^2
    incr.L <- t(matrix(c(increment,quadratVariation),ncol=2))
    incr.W <- matrix(0, nrow=1,ncol=length(increment))
    # Simulate trajectories Cogarch
  }
  d.size <- model@equation.number
  if(missing(xinit)){
    xinit <- yuimaCogarch@model@xinit
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
  if(method=="euler"){
#   result <- aux.simulate(object=yuimaCogarch, nsim=nsim, seed=seed, xinit=xinit,
#                          true.parameter = true.parameter,
#                          space.discretized = space.discretized,increment.W =incr.W,
#                          increment.L=incr.L, method=method,
#                     hurst=hurst,methodfGn=methodfGn,
#                     sampling=sampling, subsampling=subsampling,
#                     Initial=Initial, Terminal=Terminal, n=n, delta=delta,
#                     grid=grid, random=random, sdelta=sdelta,
#                     sgrid=sgrid, interpolation=interpolation)



    Terminal <- samp@Terminal[1]
    n <- samp@n[1]
    Delta <- Terminal/n
    name.ar <- paste0(info@ar.par,c(1:info@q))
    name.ma <- paste0(info@ma.par,c(1:info@p))
    name.loc <- info@loc.par
    name.param <- names(true.parameter)
    parms <- as.numeric(true.parameter)
    names(parms)<-name.param
    value.ar <- parms[name.ar]
    value.ma <- parms[name.ma]
    value.a0 <- parms[name.loc]
    AMatrix <- MatrixA(value.ar)
    avect<-evect<-matrix(0,info@q,1)
    evect[info@q,] <- 1
    avect[c(1,info@p),1] <- value.ma
    Indent<-diag(info@q)
  # Inputs: incr.L
    tavect<-t(avect)

    ncolsim <- (info@q+2)
    sim <- matrix(0,n+1,ncolsim)

    par.len <- length(model@parameter@all)
    if(missing(true.parameter) & par.len>0){
      true.parameter <- vector(par.len, mode="list")
      for(i in 1:par.len)
        true.parameter[[i]] <- 0
        names(true.parameter) <-   model@parameter@all
      }

      yuimaEnv <- new.env()

      if(par.len>0){
        for(i in 1:par.len){
          pars <- model@parameter@all[i]
          for(j in 1:length(true.parameter)){
            if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
              assign(model@parameter@all[i], true.parameter[[j]], yuimaEnv)
            }
          }
        #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
        }
      }

      for(i in c(1:ncolsim)){
        sim[1,i] <- eval(xinit[i], yuimaEnv)
      }

      for(t in c(2:n)){
          #sim[t,3:ncolsim] <- value.a0*expm(AMatrix*Delta)%*%evect*incr.L[2,t-1]+
          #  expm(AMatrix*Delta)%*%(Indent+evect%*%tavect*incr.L[2,t-1])%*%sim[t-1,3:ncolsim]
          #         sim[t,2]<-value.a0+tavect%*%sim[t,3:ncolsim]
          #         sim[t,1]<-sim[t-1,1]+sqrt(sim[t,2])*incr.L[1,t]
          #        sim[t,3:ncolsim]<-expm(AMatrix*Delta)%*%sim[t-1,3:ncolsim]+expm(AMatrix)%*%evect*sim[t-1,2]*incr.L[2,t]
          sim[t,2]<-value.a0+tavect%*%sim[t-1,3:ncolsim]
          sim[t,3:ncolsim]<-sim[t-1,3:ncolsim]+(AMatrix*Delta)%*%sim[t-1,3:ncolsim]+evect*sim[t-1,2]*incr.L[2,t]
          sim[t,1]<-sim[t-1,1]+sqrt(sim[t,2])*incr.L[1,t]
      }
      X <- ts(sim[-(samp@n[1]+1),])
      Data <- setData(X,delta = Delta)
      result <- setYuima(data=Data,model=yuimaCogarch@model, sampling=yuimaCogarch@sampling)











  }else{
    Terminal <- samp@Terminal[1]
    n <- samp@n[1]
    Delta <- Terminal/n
    name.ar <- paste0(info@ar.par,c(1:info@q))
    name.ma <- paste0(info@ma.par,c(1:info@p))
    name.loc <- info@loc.par
    name.param <- names(true.parameter)
    parms <- as.numeric(true.parameter)
    names(parms)<-name.param
    value.ar <- parms[name.ar]
    value.ma <- parms[name.ma]
    value.a0 <- parms[name.loc]
    AMatrix <- MatrixA(value.ar)
    avect<-evect<-matrix(0,info@q,1)
    evect[info@q,] <- 1
    avect[c(1,info@p),1] <- value.ma
    Indent<-diag(info@q)
    # Inputs: incr.L
    tavect<-t(avect)

    ncolsim <- (info@q+2)
    sim <- matrix(0,n+1,ncolsim)

    par.len <- length(model@parameter@all)
    if(missing(true.parameter) & par.len>0){
      true.parameter <- vector(par.len, mode="list")
      for(i in 1:par.len)
        true.parameter[[i]] <- 0
      names(true.parameter) <-   model@parameter@all
    }

    yuimaEnv <- new.env()

    if(par.len>0){
      for(i in 1:par.len){
        pars <- model@parameter@all[i]
        for(j in 1:length(true.parameter)){
          if( is.na(match(pars, names(true.parameter)[j]))!=TRUE){
            assign(model@parameter@all[i], true.parameter[[j]], yuimaEnv)
          }
        }
        #assign(sdeModel@parameter@all[i], true.parameter[[i]], yuimaEnv)
      }
    }

    for(i in c(1:ncolsim)){
      sim[1,i] <- eval(xinit[i], yuimaEnv)
    }

    if(yuimaCogarch@model@measure.type=="code"){
            for(t in c(2:n)){

#         sim[t,2]<-value.a0+tavect%*%sim[t,3:ncolsim]
#         sim[t,1]<-sim[t-1,1]+sqrt(sim[t,2])*incr.L[1,t]
#        sim[t,3:ncolsim]<-expm(AMatrix*Delta)%*%sim[t-1,3:ncolsim]+expm(AMatrix)%*%evect*sim[t-1,2]*incr.L[2,t]
#        sim[t,3:ncolsim]<-sim[t-1,3:ncolsim]+AMatrix*Delta%*%sim[t-1,3:ncolsim]+evect*sim[t-1,2]*incr.L[2,t-1]
        sim[t,2]<-value.a0+tavect%*%sim[t-1,3:ncolsim]
        sim[t,3:ncolsim] <- value.a0*expm(AMatrix*Delta)%*%evect*incr.L[2,t]+
          expm(AMatrix*Delta)%*%(Indent+evect%*%tavect*incr.L[2,t])%*%sim[t-1,3:ncolsim]
        sim[t,1]<-sim[t-1,1]+sqrt(sim[t,2])*incr.L[1,t]

      }
      X <- ts(sim[-(samp@n[1]+1),])
      Data <- setData(X,delta = Delta)
      result <- setYuima(data=Data,model=yuimaCogarch@model, sampling=yuimaCogarch@sampling)
    return(result)
    }else{
        lambda <- eval(model@measure$intensity, yuimaEnv)


        #Simulating jump times
      #intensity <- lambda*Delta
      intensity<-lambda
      jump_time<-numeric()
      jump_time[1] <- rexp(1, rate = intensity)
      # In yuima this part is evaluated using function eval
      Time <-numeric()
      Time[1] <- jump_time[1]
      j <- 1
      numb_jum<-numeric()
#       for (i in c(1:n) ){
#         numb_jum[i]<-0
#         while(Time[j]<i){
#           numb_jum[i]<-numb_jum[i]+1
#           jump_time[j+1]<-rexp(1,rate=intensity)
#           Time[j+1]<-Time[j]+jump_time[j+1]
#           j<-j+1
#         }
#       }

      while(Time[j] < Terminal){
        jump_time[j+1]<-rexp(1,rate=intensity)
        Time[j+1]<-Time[j]+jump_time[j+1]
        j<-j+1
      }

      total_NumbJ <- j
      # Counting the number of jumps
#       N<-matrix(1,n,1)
#       N[1,1]<-numb_jum[1]
#       for(i in c(2:n)){
#         N[i,1]=N[i-1,1]+numb_jum[i]
#       }
      # Simulating the driving process
      F <- suppressWarnings(parse(text=gsub("^d(.+?)\\(.+?,", "r\\1(total_NumbJ,", model@measure$df$expr, perl=TRUE)))
      assign("total_NumbJ",total_NumbJ, envir=yuimaEnv)
      dL<-eval(F, envir=yuimaEnv)
      #dL<-rnorm(total_NumbJ,mean=0,sd=1)
#       L<-matrix(1,total_NumbJ,1)
#       L[1]<-dL[1]
#       for(j in c(2:total_NumbJ)){
#         L[j]<-L[j-1] + dL[j]
#       }
      # Computing the processes V and Y at jump
      V<-matrix(1,total_NumbJ,1)
      Y<-matrix(1,info@q,total_NumbJ)
      Y[,1]<-matrix(sim[1,c(3:(3+info@q-1))],info@q,1) #Starting point for unobservable State Space Process.

      V[1,]<-value.a0+sum(tavect*Y[,1])
      G<-matrix(1, total_NumbJ,1)
      G[1]<-0

      for(j in c(2:total_NumbJ)){
        Y[,j]<-as.numeric(expm(AMatrix*jump_time[j])%*%Y[,j-1])+(V[j-1,])*evect*dL[j]^2
        V[j,]<-value.a0+sum(tavect*Y[,j])
        #       }
#       # Computing the process G at jump time
#
#       for(j in c(2:total_NumbJ)){
        G[j]<-G[j-1]+sqrt(V[j-1])*dL[j]
      }


        res<-approx(x=c(0,Time), y = c(0,G),
                    xout=seq(0,Terminal, by=Terminal/n),
                    method = "constant")
        sim[,1]<-res$y
        i<-1
        for(j in 1:length(Time)){
          while (i*Delta < Time[j] && i <= n){
            sim[i+1,3:ncolsim]<-expm(AMatrix*(Time[j]-i*Delta))%*%Y[,j]
            sim[i+1,2]<-value.a0+as.numeric(tavect%*%sim[i,3:ncolsim])
            i<-i+1

          }
        }


#       # Realizations observed at integer times
#       i<-1
#       while(N[i]==0){
#         i <- i+1
#       }
# #       G_obs<-numeric()
#       # L_obs<-numeric()
# #       V_obs<-numeric()
# #       Y_obs<-matrix(0,info@q,)
#       sim[c(1:(i-1)),1]<-0
#       sim[c(i:n),1]<-G[N[c(i:n)]]
# #       L_obs[c(1:(i-1))]<-0
# #       L_obs[c(i:n)]<-L[N[c(i:n)]]
#       for(j in c(1:(i-1))){
#         sim[j,3:ncolsim]<-as.numeric(Y[,j])
#         sim[j,2]<-value.a0+tavect%*%expm(AMatrix*j)%*%matrix(1,info@q,1)#Starting point for unobservable State Space Process
#       }
#       for(j in c(i:n)){
#         sim[j,3:ncolsim]<-as.numeric(Y[,N[j]])
#         sim[j,2]<-value.a0+as.numeric(tavect%*%expm(AMatrix*(j-Time[N[j]]))%*%Y[,N[j]])
#       }
    }
  X <- ts(sim[-1,])
  Data <- setData(X,delta = Delta)
  result <- setYuima(data=Data,model=yuimaCogarch@model, sampling=yuimaCogarch@sampling)
  }
  return(result)
}

# Simulate method for an object of class cogarch.gmm.incr

setMethod("simulate","cogarch.gmm.incr",
          function(object, nsim=1, seed=NULL, xinit,  ...){

              out <-aux.simulategmm(object=object, nsim=nsim, seed=seed, xinit=xinit, ...)
#             out <- simulate(object = model, nsim = nsim, seed=seed, xinit=xinit,
#                      sampling = samp,
#                      method = "euler",
#                      increment.L = t(as.matrix(c(0,Incr.L))),
#                      true.parameter = true.parameter,
#                      )

             return(out)
          }
          )

aux.simulategmm<-function(object, nsim=1, seed=NULL, xinit, ...){
  Time<-index(object@Incr.Lev)
  Incr.L<-coredata(object@Incr.Lev)

  model <- object@model
  EndT <- Time[length(Time)]
  numb <- (length(Incr.L)+1)
  valpar<-coef(object)

  idx <- na.omit(match(names(valpar),model@parameter@xinit))
  solnam <- model@parameter@xinit[-idx]
  solval <- as.numeric(Diagnostic.Cogarch(object, display=FALSE)$meanStateVariable)

#  solval <-50.33

  names(solval) <- solnam

  true.parameter <- as.list(c(valpar,solval))

  samp <- setSampling(Initial = 0, Terminal = EndT, n = numb)
  out <- simulate(object = model, nsim = nsim, seed=seed, xinit=xinit,
                  sampling = samp,
                  method = "euler",
                  increment.L = t(as.matrix(c(0,Incr.L))),
                  true.parameter = true.parameter
  )
  return(out)
}
