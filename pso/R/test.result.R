setClass(Class="test.result",
         representation=representation(
           problem="test.problem",
           result="list",
           time="numeric"))

setMethod(f="show",
          signature="test.result",
          definition=function(object) {
            cat("TEST RESULTS\n")
            show(object@problem)
            r <- sapply(object@result,function(x) x$value)
            rate <- getSuccessRate(object)
            rate$feval <- c(0,rate$feval,object@problem@maxf)
            rate$rate <- c(0,rate$rate,1)
            rate <- sum((rate$rate[2:length(rate$rate)]+
                         rate$rate[1:(length(rate$rate)-1)])*
                        (rate$feval[2:length(rate$feval)]-
                         rate$feval[1:(length(rate$feval)-1)]))/
                           (2*object@problem@maxf)
            cat("Objective mean: ",mean(r),"\n",sep="")
            cat("Objective s.d.: ",sd(r),"\n",sep="")
            cat("Objective min: ",min(r),"\n",sep="")
            cat("Objective max.: ",max(r),"\n",sep="")
            cat("Success rate: ",100*sum(r<=object@problem@objective)/
                object@problem@ntest,"%\n",sep="")
            cat("Efficiency: ",100*rate,"%\n",sep="")
            cat("Timing: ",object@time[1]," s. (user), ",
                object@time[2]," s. (system), ",
                object@time[3]," s. (elapsed)\n",sep="")
            invisible(object)
          })

setGeneric("getSuccessRate", function(object) standardGeneric("getSuccessRate"))
setMethod("getSuccessRate", signature(object="test.result"),
          function(object) {
            r <- sapply(object@result,function(x) x$value)
            fn <- sapply(object@result,function(x)
                         as.integer(x$counts["function"]))
            z <- sort(fn,index.return=TRUE)
            r <- r[z$ix]<=object@problem@objective
            rate <- sapply(1:length(r),function(i) sum(r[1:i]))/
              object@problem@ntest
            r <- !duplicated(paste(rate,z$x))
            return(list(feval=z$x[r],rate=rate[r]))
          })

if (!isGeneric("plot")) {
  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}

setMethod("plot", signature(x="test.result", y="missing"),
          function(x, y, ...) {
            r <- getSuccessRate(x)
            args <- list(x=r$feval,
                         y=r$rate*100,
                         type="l",
                         xlab="function evaluations",
                         ylab="success rate (%)")
            args[names(list(...))] <- list(...)
            do.call(plot,args)
          })

if (!isGeneric("lines")) {
  setGeneric("lines", function(x, ...) standardGeneric("lines"))
}

setMethod("lines", signature(x="test.result"),
          function(x, ...) {
            r <- getSuccessRate(x)
            lines(x=r$feval, y=r$rate*100, ...)
          })

if (!isGeneric("points")) {
  setGeneric("points", function(x, ...) standardGeneric("points"))
}

setMethod("points", signature(x="test.result"),
          function(x, ...) {
            r <- getSuccessRate(x)
            points(x=r$feval, y=r$rate*100, ...)
          })
