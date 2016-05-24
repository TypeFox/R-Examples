.constructArg.list <- function(fun,matchCall, onlyFormal=FALSE, debug =FALSE){
    dprint <- function(...) if(debug) print(...)
    form0 <- form <- formals(fun)
    form0["..."] <- NULL
    dprint("----------------------")
    es.call.0 <- list(as.list(matchCall)[[1]])
    if(! "..." %in% names(form) && onlyFormal){
       match.i <- 2
       while(match.i<= length(matchCall)){
           nm0 <- names(matchCall)[match.i]
           if(!nm0 %in% names(form0)){
              es.call.0[[nm0]] <- if(is.call(matchCall[[match.i]])){
                       eval(matchCall[[match.i]])
                       }else matchCall[[match.i]]
              matchCall[[match.i]] <- NULL
           }else match.i <- match.i + 1
       }
    }else{es.call.0 <- matchCall}
    dprint("----------------------")
    dprint(matchCall)
    dprint("----------------------")
    for(form.i in seq(along=form0)) {
        if(!is.null(names(form0)[form.i])){
           dprint(form.i)
           dprint(names(form0)[form.i])
           dprint(form0[[form.i]])
           dprint("----------------------")
           nam0 <- names(form0)[form.i]
           if(!nam0 %in% names(matchCall) &&
              !is.symbol(form0[[form.i]])
              ){
                  dprint("Try assign")
                  dprint(paste(nam0,": ",sep=""))
                  dprint(matchCall[[nam0]])
                  if(!is.null(form0[[form.i]])){
                      dprint(dum <- str(form0[[form.i]]))
                      xu <- if(is.call(form0[[form.i]])){
                               eval(form0[[form.i]])
                               }else form0[[form.i]]
                      if(!is.null(xu)) {dprint(xu)
                      matchCall[[nam0]] <- xu
                      dprint(matchCall[[nam0]])}
                  }else{
                      matchCall[nam0] <- NULL
                      dprint("NULL")
                  }
                  dprint("----------------------")

              }
           }
    }
    mc <- list(mc=matchCall)
    if(length(es.call.0)) mc$esc <- as.call(es.call.0)
    dprint(mc)
    return(mc)
}

.fix.in.defaults <- function(call.list, fun, withEval=TRUE){
 formals.fun <- formals(fun)
 k <- length(call.list)
 L <- length(formals.fun)
 if("..." %in% names(formals.fun)) L <- L-1
# cat("\nFORMALS\n");print(formals.fun); cat("\nCALL\n");print(call.list)
 for(i in 1:L){
     foP <- paste(formals.fun[[i]])
     if(length(formals.fun[[i]]))
     if(length(foP)||foP!=""){
        if(!names(formals.fun)[i] %in% names(call.list)&&!is.null(formals.fun[[i]])){
           k <- k + 1
           call.list[[k]] <- formals.fun[[i]]
           if(withEval && is.call(formals.fun[[i]])){
               xu <- eval(formals.fun[[i]])
               if(!is.null(xu)) call.list[[k]] <- xu
           }
           if(length(foP)){
              if(withEval && is.name(formals.fun[[i]])){
                  xu <- if(foP!="") get(foP) else NULL
                  if(!is.null(xu)) call.list[[k]] <- xu
              }
           }
           names(call.list)[k] <- names(formals.fun)[i]
        }
     }
  }
# cat("\nNEWLIST\n"); print(call.list)
 return(call.list)

}

.pretreat <- function(x, na.rm = TRUE){
    if(missing(x))
        stop("'x' is missing with no default")
    if(!is.numeric(x)){
        if(is.data.frame(x))
            x <- data.matrix(x)
        else
            x <- as.matrix(x)
        if(!is.matrix(x))
            stop("'x' has to be a numeric vector resp. a matrix or data.frame")
    }
    completecases <- complete.cases(x)
    if(na.rm) x <- na.omit(x)
    return(list(x=x,completecases=completecases))
}

.check.eps <- function(...){
   mc <- as.list(match.call(expand.dots=TRUE)[-1])

   eps <- eps.lower <- eps.upper <- NULL

   ine <- is.null(mc[["eps"]]) || is.symbol(mc[["eps"]])
   inl <- is.null(mc[["eps.lower"]]) || is.symbol(mc[["eps.lower"]])
   inu <- is.null(mc[["eps.upper"]]) || is.symbol(mc[["eps.upper"]])
   if(ine && inl && inu){
        eps.lower <- 0
        eps.upper <- 0.5
    }else if(ine){
        if(!inl && !inu)
            eps.lower <- mc[["eps.lower"]]
            eps.upper <- mc[["eps.upper"]]
        if(!inl && inu)
            eps.lower <- mc[["eps.lower"]]
            eps.upper <- 0.5
        if(inl && !inu)
            eps.lower <- 0
            eps.upper <- mc[["eps.upper"]]
        if(length(eps.lower) != 1 || length(eps.upper) != 1)
            stop("'eps.lower' and 'eps.upper' have to be of length 1")
        if(!is.numeric(eps.lower) || !is.numeric(eps.upper) || eps.lower >= eps.upper)
            stop("'eps.lower' < 'eps.upper' is not fulfilled")
        if((!inl &&(eps.lower < 0)) || (!inu&& (eps.upper > 0.5)))
            stop("'eps.lower' and 'eps.upper' have to be in [0, 0.5]")
    }else{
        eps <- mc[["eps"]]
        if(length(eps) != 1)
            stop("'eps' has to be of length 1")
        if(!ine&& (eps == 0))
            stop("'eps = 0'! => use functions 'mean' and 'sd' for estimation")
        if(!ine && ((eps < 0) || (eps > 0.5)))
            stop("'eps' has to be in (0, 0.5]")
    }

    x <- mc[["x"]]
    if(is.matrix(x))
        sqrtn <- sqrt(ncol(x))
    else
        sqrtn <- sqrt(length(x))
    erg <- list(sqn = sqrtn)
    erg[["e"]] <- eps
    erg[["lower"]] <- eps.lower
    erg[["upper"]] <- eps.upper
    return(erg)
}

.isOKsteps <- function(steps){
    if(!is.integer(steps))
        steps <- as.integer(steps)
    if(steps < 1){
        stop("'steps' has to be some positive integer value")
    }
    if(length(steps) != 1){
        stop("'steps' has to be of length 1")
    }
   return(invisible(NULL))
}

.isOKfsCor <- function(fsCor){
    if(fsCor <= 0)
        stop("'fsCor' has to be positive")
    if(length(fsCor) != 1)
        stop("'fsCor' has to be of length 1")
   return(invisible(NULL))
}


#.getROptICstart <- function(..., ..eval=FALSE){
#    mc <- match.call(expand.dots=TRUE)
#    eps <- mc$eps
#    dots <- mc$dots
#
#    if(is.null(eps$e)){
#        r.lower <- eps$sqn * eps$lower
#        r.upper <- eps$sqn * eps$upper
#        ICstart <- substitute(do.call(radiusMinimaxIC,
#                    c(list(L2Fam = mc$L2FamStart, neighbor = mc$neighbor,
#                                   risk = mc$risk,
#                                   loRad = r.lower, upRad = r.upper,
#                                   verbose = mc$verbose,
#                                   OptOrIter = mc$OptOrIter),dots)))
#        if(..eval) ICstart <- eval(ICstart)
#
#        if(!isTRUE(all.equal(mc$fsCor, 1, tol = 1e-3))){
#            neighbor@radius <- neighborRadius(ICstart)*mc$fsCor
#            infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
#            ICstart <- substitute(do.call(optIC, c(list( model = mc$infMod, risk = mc$risk,
#                               verbose = mc$verbose, OptOrIter = mc$OptOrIter),
#                               dots)))
#            if(..eval) ICstart <- eval(ICstart)
#       }
#    }else{
#        neighbor@radius <- eps$sqn*eps$e*mc$fsCor
#        infMod <- InfRobModel(center = mc$L2FamStart, neighbor = mc$neighbor)
#        ICstart <- substitute(do.call(optIC, c(list(model = mc$infMod, risk = mc$risk,
#                           verbose = mc$verbose, OptOrIter = mc$OptOrIter),
#                           dots)))
#        if(..eval) ICstart <- eval(ICstart)
#
#    }
#  return(ICstart)
#}

genkStepCtrl <- function(useLast = getRobAStBaseOption("kStepUseLast"),
                    withUpdateInKer = getRobAStBaseOption("withUpdateInKer"),
                    IC.UpdateInKer = getRobAStBaseOption("IC.UpdateInKer"),
                    withICList = getRobAStBaseOption("withICList"),
                    withPICList = getRobAStBaseOption("withPICList"),
                    scalename = "scale", withLogScale = TRUE,
                    withEvalAsVar = NULL){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,genkStepCtrl)
 return(es.list)
}
genstartCtrl<- function(initial.est = NULL, initial.est.ArgList = NULL,
                        startPar = NULL, distance = CvMDist, withMDE = NULL){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,genstartCtrl)
 return(es.list)
}
gennbCtrl <- function(neighbor = ContNeighborhood(),
                      eps, eps.lower, eps.upper){
  es.call <- match.call()
  es.list <- as.list(es.call[-1])
  es.list <- .fix.in.defaults(es.list,gennbCtrl)
 return(es.list)
}
