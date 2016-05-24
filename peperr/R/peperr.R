`peperr` <-
function(response, x, 
   indices=NULL,
   fit.fun, complexity=NULL, args.fit=NULL, args.complexity=NULL, 
   parallel=NULL, cpus=2, clustertype=NULL, clusterhosts=NULL,
   noclusterstart=FALSE, noclusterstop=FALSE,
   aggregation.fun=NULL, args.aggregation=NULL,
   load.list=extract.fun(list(fit.fun, complexity, aggregation.fun)), 
   load.vars=NULL, load.all=FALSE,
   trace=FALSE, debug=FALSE, peperr.lib.loc=NULL, 
   RNG=c("RNGstream", "SPRNG", "fixed", "none"), seed=NULL, 
   lb=FALSE, sr=FALSE, sr.name="default", sr.restore=FALSE)
{
  binary <- FALSE
   if(is.null(aggregation.fun)){
      if(is.Surv(response)){
         aggregation.fun <- aggregation.pmpec
      } else {
         if (is.vector(response)&& length(unique(response))<3){
            aggregation.fun <- aggregation.brier
         } else {
            stop("Please specify argument 'aggregation.fun' according to 'response'")  
         }
      }
   }

   if(is.null(indices)){
      if(nrow(x)<ncol(x)){
         indices <- resample.indices(n=nrow(x), method="sub632", sample.n=500)
      } else {
         indices <- resample.indices(n=nrow(x), method="boot", sample.n=500)
      }
   }

   if(!is.null(parallel)){
      sfInit(parallel=parallel, cpus=cpus, nostart=noclusterstart, 
         type=clustertype, socketHosts=clusterhosts)
      } else { 
      sfInit(nostart=noclusterstart)
   }
   sfLibrary(peperr, lib.loc=peperr.lib.loc)
   RNG <- match.arg(RNG)
   if (RNG=="RNGstream" && is.null(seed)){
      warning("You are using a parallel random number generator ('RNGstream') with its default seed. See Details of documentation for other options.")
   }
   if (RNG!="none"){
      if (RNG!="fixed"){
         if (!is.null(seed)){
            sfClusterSetupRNG(type=RNG, seed=seed)
         } else {
            sfClusterSetupRNG(type=RNG)
         }
      } else {
         if (length(seed)==1 || length(seed)==(length(indices$sample.index)+2)){
            set.seed(seed[1])
         } else {
            stop("Provide argument 'seed' that is integer or vector of length number of samples plus 2")
         }
      }
   }
   if(load.all==TRUE){
      sfExportAll()
      for (i in 1:length(.packages())){
         eval(call("sfLibrary", (.packages()[i]), character.only=TRUE))
      }
   } else {
      if (!is.null(load.list)){
         for (i in seq(along=load.list$packages)){
            try(eval("sfLibrary"(load.list$packages[i], 
               character.only=TRUE)), silent=!debug)
         }
         for (i in seq(along=load.list$functions)){
            try(eval(call("sfExport", load.list$functions[i], debug=debug)), silent=!debug)
         }
         for (i in seq(along=load.list$variables)){
            try(eval(call("sfExport", load.list$variables[i], debug=debug)), silent=!debug)
         }
      }
      if (!is.null(load.vars)){
         sfExport(names(load.vars))
      }
   }

   null.model <- NULL
   actual.data <- as.data.frame(x)
   if (is.Surv(response)){
      #xnames <- names(actual.data) 
      time <- response[,"time"]
      status <- response[, "status"]
      actual.data$time <- time
      actual.data$status <- status
      km.fit <- survfit(Surv(time, status)~1, data=actual.data)
      null.model <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
         response=response, x=x, model=km.fit), args.aggregation))
   } else {
      if (is.vector(response)&& length(unique(response))<3){
         binary <- TRUE
         actual.data$response <- response
         logreg.fit <- glm(formula=response~1, data=actual.data, family=binomial())
         null.model <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
            response=response, x=x, model=logreg.fit), args.aggregation))
      }
   }
   dimnames(null.model) <- NULL
   attr(null.model, "addattr") <- attributes(null.model)$addattr
   fullsample.attr <- attr(null.model, "addattr")

   if (is.Surv(response)){
   km.pred <- summary(object=km.fit, times=fullsample.attr)$surv 
   km.weight <- -1*diff(km.pred)
   }

   fullsample.index <- resample.indices(n=nrow(x), method="no")
   sample.index <- indices$sample.index
   not.in.sample <- indices$not.in.sample
   sample.index.full <- indices$sample.index
   not.in.sample.full <- indices$not.in.sample
   sample.n <- length(sample.index)+1
   sample.index.full[[sample.n]] <- fullsample.index$sample.index[[1]]
   not.in.sample.full[[sample.n]] <- fullsample.index$not.in.sample[[1]]

   pP. <- ls(pattern="predictProb.", name=".GlobalEnv")
   try(eval(call("sfExport", pP., debug=debug)), silent=!debug)
   PLL. <- ls(pattern="PLL.", name=".GlobalEnv")
   try(eval(call("sfExport", PLL., debug=debug)), silent=!debug)
   
   if (environmentName(environment(fit.fun))!="peperr"){sfExport("fit.fun")}
   if (environmentName(environment(aggregation.fun))!="peperr"){sfExport("aggregation.fun")}
   if (environmentName(environment(complexity))!="peperr"){sfExport("complexity")}
   sfExport("response", "x", "sample.index.full", "sample.n", "not.in.sample.full", 
      "args.fit", "args.complexity", "args.aggregation", "binary", "RNG", "seed")
   try(sfExport("km.weight"), silent=!debug)

   if(trace){message("Evaluation on slaves starts now")}

   sample.fun <- function(actual.sample){
       if (RNG=="fixed"){
         if (length(seed)==1){
           set.seed(seed+actual.sample)
         } else {
           set.seed(seed[actual.sample+1])
         } 
      }

      if (trace && actual.sample<sample.n){
         cat("Sample run", actual.sample, "of", (sample.n-1), "\n")
      }
      if (!is.Surv(response)&&!is.matrix(response)){
         response <- matrix(response, ncol=1)
      }
      if (is.function(complexity)){
         sample.complexity <- do.call("complexity",
            c(list(response=response[unique(sample.index.full[[actual.sample]]),],
            x=x[unique(sample.index.full[[actual.sample]]),, drop=FALSE], full.data=actual.data), args.complexity))
      } else {
         if (is.vector(complexity)){
            sample.complexity <- complexity
         } else sample.complexity <- 0
      }

      if (is.Surv(response)){
         lipec.oob <- c()
         lipec.oob.null <- c()
         pll.oob <- c()
         #pll.oob.null <- c()

      }
      actual.error <- c()
      sample.fit.list <- list()
      if (is.list(sample.complexity)){
         for (i in 1:length(sample.complexity[[1]])){
            list.sample.complexity <- lapply(sample.complexity, function(arg) arg[i])
            sample.fit <- do.call("fit.fun", c(list(response=response[unique(sample.index.full[[actual.sample]]),],
               x=x[unique(sample.index.full[[actual.sample]]),, drop=FALSE], 
               cplx=list.sample.complexity), args.fit))
            sample.fit.list[[i]] <- sample.fit
         #class(sample.fit) <- c(class(sample.fit), "peperrinterinternal")

            actual.error.i <-  do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",
               response=response[not.in.sample.full[[actual.sample]],], 
               x=x[not.in.sample.full[[actual.sample]],, drop=FALSE], model=sample.fit, 
               cplx=list.sample.complexity, fullsample.attr=fullsample.attr), 
               args.aggregation))
            actual.error <- rbind(actual.error, actual.error.i)
 
            if (is.Surv(response)){
               km.fit <- survival::survfit(Surv(time, status)~1,
                  data=actual.data[unique(sample.index.full[[actual.sample]]),])
               km.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
                  response=response[unique(not.in.sample.full[[actual.sample]]),],
                  x=x[unique(not.in.sample.full[[actual.sample]]),, drop=FALSE], model=km.fit, 
                  fullsample.attr=fullsample.attr), args.aggregation))
               lipec.oob.null.i  <- sum(km.apparent[1:(length(km.weight))]*km.weight, na.rm=TRUE)
               lipec.oob.null <- rbind(lipec.oob.null, lipec.oob.null.i)
               lipec.oob.i <- sum(actual.error.i[1:(length(km.weight))]*km.weight, na.rm=TRUE)
               lipec.oob <- rbind(lipec.oob, lipec.oob.i)
               if (exists(paste("PLL.", class(sample.fit), sep=""))){
# 		   pll.oob.null.i <- PLL(object=km.fit, newdata=x[not.in.sample.full[[actual.sample]],, drop=FALSE],
#                      newtime=time[not.in.sample.full[[actual.sample]]], 
#                      newstatus=status[not.in.sample.full[[actual.sample]]], complexity=list.sample.complexity)
#                   pll.oob.null <- rbind(pll.oob.null, pll.oob.null.i)
                  pll.oob.i <- try(PLL(object=sample.fit, newdata=x[not.in.sample.full[[actual.sample]],, drop=FALSE],
                     newtime=time[not.in.sample.full[[actual.sample]]], 
                     newstatus=status[not.in.sample.full[[actual.sample]]], complexity=list.sample.complexity))
                  pll.oob <- rbind(pll.oob, pll.oob.i)
               }
            } else {
               if(binary){
                  logreg.fit <- glm(formula=response~1, data=actual.data[sample.index.full[[actual.sample]],],
                     family=binomial())
                  logreg.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",  
                     response=response[unique(sample.index.full[[actual.sample]]),],
                     x=x[unique(sample.index.full[[actual.sample]]),, drop=FALSE], model=logreg.fit), 
                     args.aggregation))
                  }
               }
            }
         } else {
            if (is.vector(sample.complexity)){
                for (i in 1:length(sample.complexity)){
                   sample.fit <- do.call("fit.fun",
                      c(list(response=response[unique(sample.index.full[[actual.sample]]),],
                      x=x[unique(sample.index.full[[actual.sample]]),, drop=FALSE], 
                     cplx=sample.complexity[i]), args.fit))
                     sample.fit.list[[i]] <- sample.fit
         #class(sample.fit) <- c(class(sample.fit), "peperrinterinternal")

                   actual.error.i <-  do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",
                      response=response[not.in.sample.full[[actual.sample]],], 
                      x=x[not.in.sample.full[[actual.sample]],,drop=FALSE], model=sample.fit, 
                      cplx=sample.complexity[i], fullsample.attr=fullsample.attr), 
                      args.aggregation))
                   actual.error <- rbind(actual.error, actual.error.i)

                   if (is.Surv(response)){
                     km.fit <- survival::survfit(Surv(time, status)~1,
                        data=actual.data[sample.index.full[[actual.sample]],])
cat("Kaplan-Meier, step", actual.sample)
                     km.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
                        response=response[unique(not.in.sample.full[[actual.sample]]),],
                        x=x[unique(not.in.sample.full[[actual.sample]]),, drop=FALSE], model=km.fit,
                        fullsample.attr=fullsample.attr), args.aggregation))
                     lipec.oob.null.i  <- sum(km.apparent[1:(length(km.weight))]*km.weight, na.rm=TRUE)
                     lipec.oob.null <- rbind(lipec.oob.null, lipec.oob.null.i)
                     lipec.oob.i <- sum(actual.error.i[1:(length(km.weight))]*km.weight, na.rm=TRUE)
                     lipec.oob <- rbind(lipec.oob, lipec.oob.i)
                     if (exists(paste("PLL.", class(sample.fit), sep=""))){
#  			 pll.oob.null.i <- PLL(object=km.fit, newdata=x[not.in.sample.full[[actual.sample]],, drop=FALSE],
#                      newtime=time[not.in.sample.full[[actual.sample]]], 
#                      newstatus=status[not.in.sample.full[[actual.sample]]], complexity=list.sample.complexity)
#                         pll.oob.null <- rbind(pll.oob.null, pll.oob.null.i)
                        pll.oob.i <- try(PLL(object=sample.fit, newdata=x[not.in.sample.full[[actual.sample]],, drop=FALSE],
                           newtime=time[not.in.sample.full[[actual.sample]]], 
                           newstatus=status[not.in.sample.full[[actual.sample]]], complexity=sample.complexity[i]))
                        pll.oob <- rbind(pll.oob, pll.oob.i)
                     }
                  } else {
                     if (binary){
                     logreg.fit <- glm(formula=response~1, data=actual.data[sample.index.full[[actual.sample]],],
                        family=binomial())
                     logreg.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",  
                        response=response[unique(sample.index.full[[actual.sample]]),],
                        x=x[unique(sample.index.full[[actual.sample]]),, drop=FALSE], 
                        model=logreg.fit), args.aggregation))
                     }
                     }
                  }
              }
          }
      dimnames(actual.error) <- NULL
      
      if (actual.sample<sample.n) sample.fit.list <- NULL

      if (is.Surv(response)){
         dimnames(lipec.oob) <- NULL
         dimnames(lipec.oob.null) <- NULL
         dimnames(pll.oob) <- NULL
  #      dimnames(pll.oob.null) <- NULL
         out <- list(actual.error=actual.error, sample.complexity=sample.complexity, 
            sample.fit=sample.fit.list,
            lipec.oob=lipec.oob, lipec.oob.null=lipec.oob.null, pll.oob=pll.oob, #pll.oob.null=pll.oob.null, 
km.apparent=km.apparent)
      } else {
         if(binary){
            out <- list(actual.error=actual.error, sample.complexity=sample.complexity,
                      sample.fit=sample.fit.list, logreg.apparent=logreg.apparent)
         } else {
            out <- list(actual.error=actual.error, sample.complexity=sample.complexity,
               sample.fit=sample.fit.list)
         }
      }
   out
}
   if (lb){
      sample.error.list <- sfClusterApplyLB(as.list(1:sample.n), sample.fun)
   } else {
      if(sr){
      sample.error.list <- sfClusterApplySR(as.list(1:sample.n), sample.fun, 
                             name=sr.name, restore=sr.restore)
      } else {
      sample.error.list <- sfLapply(as.list(1:sample.n), sample.fun)
      }
   }
   Stop <- FALSE
   for(i in 1:length(sample.error.list)){
      if (class(sample.error.list[[i]])=="try-error"){ 
         if (i != sample.n){
            cat("Error in run", i, ":", sample.error.list[[i]], "\n")
         } else {
            cat("Error in full sample run:", sample.error.list[[i]], "\n")
         }
         Stop <- TRUE
      }
   }
   if (Stop) stop("Error(s) occurred in (bootstrap) sample run(s), see above")

   sample.error.full <- lapply(sample.error.list, function(arg) arg$actual.error)
   sample.complexity.full <- unlist(lapply(sample.error.list, function(arg) arg$sample.complexity))
   sample.fit.full <- lapply(sample.error.list, function(arg) arg$sample.fit)
   if (is.Surv(response)){
   sample.lipec.full <- lapply(sample.error.list, function(arg) arg$lipec.oob)
   sample.lipec.null.full <- lapply(sample.error.list, function(arg) arg$lipec.oob.null)
   sample.pll.full <- lapply(sample.error.list, function(arg) arg$pll.oob)
   #sample.pll.null.full <- lapply(sample.error.list, function(arg) arg$pll.oob.null)
   sample.km.full <- lapply(sample.error.list, function(arg) arg$km.apparent)
   attr(null.model, "addattr") <- NULL
   } else {
   if (binary){
   sample.lrm.full <- lapply(sample.error.list, function(arg) arg$logreg.fit) # not used at the time
   sample.null.model.full <- lapply(sample.error.list, function(arg) arg$logreg.apparent)
   }
   }

   if (is.function(complexity)){
      if (is.list(sample.error.list[[1]]$sample.complexity)){
         cplx <- sample.error.list[[sample.n]]$sample.complexity
      } else {
         cplx <- sample.complexity.full[[sample.n]]
      }
   } else {
      if (is.vector(complexity)){
         cplx <- complexity
      } else cplx <- 0
   }

   noinf.error <- c()

   if (is.list(cplx)){
      for (i in 1:length(cplx[[1]])){
           list.cplx <- lapply(cplx, function(arg) arg[i])
#           fullmodel <- do.call("fit.fun", c(list(response=response, x=x, 
#             cplx=list.cplx), args.fit))
#
         noinf.error.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="noinf", 
               response=response , x=x, model=sample.fit.full[[sample.n]][[i]], 
               cplx=list.cplx, fullsample.attr=fullsample.attr), args.aggregation))	
         noinf.error <- rbind(noinf.error, noinf.error.i)	
      }
   } else {
      if (is.vector(cplx)){
         for (i in 1:length(cplx)){

            noinf.error.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="noinf", 
                  response=response , x=x, model=sample.fit.full[[sample.n]][[i]], 
                  cplx=cplx[i], fullsample.attr=fullsample.attr), args.aggregation))	
            noinf.error <- rbind(noinf.error, noinf.error.i)	
         }
      }
   }
   if (is.null(fullsample.attr)) fullsample.attr <- attr(noinf.error.i, "addattr")
   dimnames(noinf.error) <- NULL

   sfStop(nostop=noclusterstop)

   if (is.vector(sample.complexity.full)){
      sample.complexity <- sample.complexity.full[-sample.n]
   }
   if (is.list(cplx)){
      sample.complexity <- sample.complexity.full[1:(length(sample.complexity.full)-length(cplx))]
   }

   if (!is.function(complexity) && !is.list(complexity)){
   sample.complexity <- unique(sample.complexity)
   }
  
   if (is.Surv(response)){
   output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
      selected.complexity=cplx, complexity=complexity, 
      response=response, full.model.fit=sample.fit.full[[sample.n]],
      full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error, 
      attribute=fullsample.attr,
      sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity,  
      sample.lipec=sample.lipec.full[1:(sample.n-1)], 
      sample.lipec.null=sample.lipec.null.full[1:(sample.n-1)],
      sample.pll=sample.pll.full[1:(sample.n-1)],
      null.model=null.model, null.model.fit=km.fit, 
      sample.null.model=sample.km.full[1:(sample.n-1)])
   } else {
      if (binary){
         output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
            selected.complexity=cplx, complexity=complexity, 
            response=response, full.model.fit=sample.fit.full[[sample.n]],
            full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error,
            attribute=fullsample.attr,
            sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity,
            null.model=null.model, null.model.fit=logreg.fit,
            sample.null.model=sample.null.model.full[1:(sample.n-1)])
      } else {
          output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
            selected.complexity=cplx, complexity=complexity, 
            response=response, full.model.fit=sample.fit.full[[sample.n]],
            full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error,
            attribute=fullsample.attr,
            sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity)
      }
   }
   class(output) <- "peperr"
   output
}


