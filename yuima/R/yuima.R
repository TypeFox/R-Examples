yuima.stop <- function(x) 
stop(sprintf("\nYUIMA: %s\n", x))

yuima.warn <- function(x) 
   warning(sprintf("\nYUIMA: %s\n", x))

# 22/11/2013
# We introduce a new utility yuima.simplify that allows us to simplify 
# the expressions in the drift, diffusion and jump terms.

# yuima.Simplify modified from the original code  Simplify.R 
# by Andrew Clausen <clausen@econ.upenn.edu> in 2007. 
# http://economics.sas.upenn.edu/~clausen/computing/Simplify.R

# This isn't a serious attempt at simplification code.  It just does some
# obvious things like 0 + x => x.  It was written to support Deriv.R.

yuima.Simplify <- function(expr, yuima.env){
    
    ### 
    
    
    Simplify_ <- function(expr)
    {
        if (is.symbol(expr)) {
            expr
        } else if (is.language(expr) && is.symbol(expr[[1]])) {
            # is there a rule in the table?
            sym.name <- as.character(expr[[1]])
            if (class(try(Simplify.rule <-
            get(sym.name, envir=yuima.env,
            inherits=FALSE), silent=TRUE))
            != "try-error")
            return(Simplify.rule(expr))
        }
        expr
    }
    
    Simplify.function <- function(f, x=names(formals(f)), env=parent.frame())
    {
        stopifnot(is.function(f))
        as.function(c(as.list(formals(f)),
        Simplify_(body(f))),
        envir=env)
    }
    
    `Simplify.+` <- function(expr)
    {
        if (length(expr) == 2)
        {
            if (is.numeric(expr[[2]]))
            return(+expr[[2]])
            return(expr)
        }
        
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            b
        } else if (is.numeric(b) && all(b == 0)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a + b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.-` <- function(expr)
    {
        if (length(expr) == 2)
        {
            if (is.numeric(expr[[2]]))
            return(-expr[[2]])
            return(expr)
        }
        
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            -b
        } else if (is.numeric(b) && all(b == 0)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a - b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.(` <- function(expr)
    expr[[2]]
    
    `Simplify.*` <- function(expr)
    {
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            0
        } else if (is.numeric(b) && all(b == 0)) {
            0
        } else if (is.numeric(a) && all(a == 1)) {
            b
        } else if (is.numeric(b) && all(b == 1)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a * b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.^` <- function(expr)
    {
        a <- Simplify_(expr[[2]])
        b <- Simplify_(expr[[3]])
        
        if (is.numeric(a) && all(a == 0)) {
            0
        } else if (is.numeric(b) && all(b == 0)) {
            1
        } else if (is.numeric(a) && all(a == 1)) {
            1
        } else if (is.numeric(b) && all(b == 1)) {
            a
        } else if (is.numeric(a) && is.numeric(b)) {
            a ^ b
        } else {
            expr[[2]] <- a
            expr[[3]] <- b
            expr
        }
    }
    
    `Simplify.c` <- function(expr)
    {
        args <- expr[-1]
        args.simplified <- lapply(args, Simplify_)
        if (all(lapply(args.simplified, is.numeric))) {
            as.numeric(args.simplified)
        } else {
            for (i in 1:length(args))
            expr[[i + 1]] <- args.simplified[[i]]
            expr
        }
    }
    
    
    assign("+", `Simplify.+`, envir=yuima.env)
    assign("-", `Simplify.-`, envir=yuima.env)
    assign("*", `Simplify.*`, envir=yuima.env)
    assign("(", `Simplify.(`, envir=yuima.env)
    assign("c", `Simplify.c`, envir=yuima.env)
    assign("^", `Simplify.^`, envir=yuima.env)
    

    ###
    
    
    
    
    
    
    
    
  as.expression(Simplify_(expr[[1]]))

}


## Constructor and Initializer of class 'yuima'

# we convert objects to "zoo" internally
# we should change it later to more flexible classes

setMethod("initialize", "yuima",
          function(.Object, data=NULL, model=NULL, sampling=NULL, characteristic=NULL, functional=NULL){  
            eqn <- NULL

            if(!is.null(data)){
              .Object@data <- data
              eqn <- dim(data)
              if(is.null(sampling))
               sampling <- setSampling(grid=list(index(get.zoo.data(data)[[1]])))
            }
            
            if(!is.null(model)){
              if(!is.null(eqn)){
                if(eqn!=model@equation.number){
                  yuima.warn("Model's equation number missmatch.")
                  return(NULL)
                }
              }else{
                eqn <- model@equation.number
              }
              .Object@model <- model
            }
            
            if(!is.null(sampling)){
              if(!is.null(eqn)){
                if(eqn!=length(sampling@Terminal)){
                  if(length(sampling@Terminal)==1){
                    sampling@Terminal <- rep(sampling@Terminal, eqn)
                    sampling@n <- rep(sampling@n, eqn)
                  }else{
                    yuima.warn("Sampling's equation number missmatch.")
                    return(NULL)
                  }
                }
              }else{
                eqn <- length(sampling@Terminal)
              }
              .Object@sampling <- sampling
            }
            
            if(!is.null(characteristic)){
              if(!is.null(eqn)){
                if(eqn!=characteristic@equation.number){
                  yuima.warn("Characteristic's equation number missmatch.")
                  return(NULL)                  
                }
              }
              .Object@characteristic <- characteristic
            }else if(!is.null(eqn)){
              characteristic <- new("yuima.characteristic", equation.number=eqn, time.scale=1)
              .Object@characteristic <- characteristic
            }
            
			if(!is.null(functional)) .Object@functional <- functional
			
            return(.Object)
          })

# setter
setYuima <-
  function(data=NULL, model=NULL, sampling=NULL, characteristic=NULL, functional=NULL){
    if(is.CARMA(model)&& !is.null(data)){
      if(dim(data@original.data)[2]==1){
        dum.matr<-matrix(0,length(data@original.data),
                         (model@info@p+1))
        dum.matr[,1]<-as.numeric(data@original.data)
        data<-setData(zoo(x=dum.matr, order.by=time(data@zoo.data[[1]])))  
      
      }       
    }
    if(is.COGARCH(model)&& !is.null(data)){
      if(dim(data@original.data)[2]==1){
#         data<-setData(zoo(x=matrix(as.numeric(data@original.data),length(data@original.data),
#                                    (model@info@p+1)), order.by=time(data@zoo.data[[1]])))  
        dum.matr<-matrix(0,length(data@original.data),
                         (model@info@q+2))
        dum.matr[,1]<-as.numeric(data@original.data)
        data<-setData(zoo(x=dum.matr, order.by=time(data@zoo.data[[1]])))  
        
        
      }       
    }
    
    # LM 25/04/15 
    return(new("yuima", data=data, model=model, sampling=sampling, characteristic=characteristic,functional=functional))
  }




setMethod("show", "yuima.functional",
function(object){
    str(object)
} )

setMethod("show", "yuima.sampling",
function(object){
    str(object)
} )


setMethod("show", "yuima.data",
function(object){
    show(setYuima(data=object))
    
} )


setMethod("show", "yuima.model",
function(object){
    show(setYuima(model=object))
    
} )

setMethod("show", "yuima",
function(object){
    
    myenv <- new.env()
    mod <- object@model
    has.drift <- FALSE
    has.diff <- FALSE
    has.fbm <- FALSE
    has.levy <- FALSE
    is.wienerdiff <- FALSE
    is.fracdiff <- FALSE
    is.jumpdiff <- FALSE
    is.carma <- FALSE
    is.cogarch <- FALSE
    is.poisson <- is.Poisson(mod)

    if(length(mod@drift)>0) has.drift <- TRUE
    if(length(mod@diffusion)>0) has.diff <- TRUE
    if(length(mod@jump.coeff)>0) has.levy <- TRUE
    if(!is.null(mod@hurst)) has.fbm <- TRUE
    
    if( has.drift | has.diff ) is.wienerdiff <- TRUE
    if( has.fbm  ) is.fracdiff <- TRUE
    if( has.levy ) is.jumpdiff <- TRUE
    ldif <- 0
    if(length(mod@diffusion)>0)
     ldif <- length(mod@diffusion[[1]])
    if(ldif==1 & (length(mod@diffusion)==0)){
     if( as.character(mod@diffusion[[1]]) == "(0)" ){
      has.diff <- FALSE
      is.wienerdiff <- FALSE
      is.fracdiff <- FALSE
     }
    }
    if( class(mod) == "yuima.carma")
     is.carma <- TRUE
    if( class(mod) == "yuima.cogarch"){
      is.cogarch <- TRUE
      is.wienerdiff <- FALSE
      is.fracdiff <- FALSE
    }
      
    if( is.wienerdiff | is.fracdiff | is.jumpdiff  ){
        if( is.wienerdiff & ! is.carma & !is.poisson & !is.cogarch){
         cat("\nDiffusion process")
         if( is.fracdiff ){
             if(!is.na(mod@hurst)){
                 if(mod@hurst!=0.5){
                  cat(sprintf(" with Hurst index:%.2f", mod@hurst))
                 }
             } else {
                 cat(" with unknown Hurst index")
             }
         }
        }
        if(is.carma)
          cat(sprintf("\nCarma process p=%d, q=%d", mod@info@p, mod@info@q))
        if(is.cogarch)
          cat(sprintf("\nCogarch process p=%d, q=%d", mod@info@p, mod@info@q))
        if(is.poisson)
          cat("\nCompound Poisson process")
          
        if( (is.jumpdiff & ! is.cogarch) ){
            if( (is.wienerdiff | is.carma) & !is.poisson ){
                cat(" with Levy jumps")
            } else {
                if(!is.poisson)
                cat("Levy jump process")
            }
        }else{
          cat(" with Levy jumps")
        }
        
        cat(sprintf("\nNumber of equations: %d", mod@equation.number))
        if((is.wienerdiff | is.fracdiff) & !is.poisson)
         cat(sprintf("\nNumber of Wiener noises: %d", mod@noise.number))
        if(is.jumpdiff)
         cat(sprintf("\nNumber of Levy noises: %d", 1))
        if(is.cogarch)
          cat(sprintf("\nNumber of quadratic variation: %d", 1))
         if(length(mod@parameter@all)>0){
             cat(sprintf("\nParametric model with %d parameters",length(mod@parameter@all)))
         }
    }
    
    if(length(object@data@original.data)>0){
        n.series <- 1
        if(!is.null(dim(object@data@original.data))){
            n.series <- dim(object@data@original.data)[2]
            n.length <- dim(object@data@original.data)[1]
        } else {
            n.length <- length(object@data@original.data)
        }
        
        cat(sprintf("\n\nNumber of original time series: %d\nlength = %d, time range [%s ; %s]", n.series, n.length, min(time(object@data@original.data)), max(time(object@data@original.data))))
    }
    if(length(object@data@zoo.data)>0){
        n.series <- length(object@data@zoo.data)
        n.length <- unlist(lapply(object@data@zoo.data, length))
        t.min <- unlist(lapply(object@data@zoo.data, function(u) as.character(round(time(u)[which.min(time(u))],3))))
        t.max <- unlist(lapply(object@data@zoo.data, function(u) as.character(round(time(u)[which.max(time(u))],3))))
        
        delta <- NULL
        is.max.delta <- rep("", n.series)
        have.max.delta <- FALSE
        for(i in 1:n.series){
            tmp <- length(table(round(diff(time(object@data@zoo.data[[i]])),5)))
            if(tmp>1){
             tmp <- max(diff(time(object@data@zoo.data[[i]])), na.rm=TRUE)
             is.max.delta[i] <- "*"
             have.max.delta <- TRUE
             #tmp <- NULL
            } else {
             tmp <- diff(time(object@data@zoo.data[[i]]))[1]
            }
            if(is.null(tmp)){
                delta <- c(delta, NA)
            } else {
                delta <- c(delta, tmp)
            }
        }
        
        
        cat(sprintf("\n\nNumber of zoo time series: %d\n", n.series))
        tmp <- data.frame(length=n.length, time.min = t.min, time.max =t.max, delta=delta)
        if(have.max.delta)
         tmp <- data.frame(tmp, note=is.max.delta)
        nm <- names(object@data@zoo.data)
        if(is.null(nm)){
         rownames(tmp) <- sprintf("Series %d",1:n.series)
        } else {
         rownames(tmp) <- nm
        }
        print(tmp)
        if(have.max.delta)
         cat("================\n* : maximal mesh")
    }
    
})

