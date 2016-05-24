#' @import methods

#' @rdname fit
#' @param object A cdtamodel object created by \link{cdtamodel} function.
#' @param data A data-frame with no missing values containg TP, TN, FP, FN, 'SID' and co-varaiables(if necessary).
#' @param SID A string indicating the name of the column with the study identifier.
#' @param chains A positive numeric value specifying the number of chains, default is 3.
#' @param iter A positive numeric value specifying the number of iterations per chain. The default is 6000.
#' @param warmup A positive numeric value (<iter) specifying the number of iterations to be discarded(burn-in/warm-up). The default is 1000.
#' @param thin A positive numeric value specifying the interval in which the samples are stored. The default is 10.
#' @param cores A positive numeric values specifying the number of cores to use to execute parallel sampling. When the hardware has more at least 4 cores,
#' the default is 3 cores and otherwise 1 core.
#' @param ... Other optional parameters as specified in \link[rstan]{stan}.
setGeneric(name="fit",
          function(object, ...){ standardGeneric("fit") })


#' Fit cdtamodel object.
#' @rdname fit
#' @method fit cdtamodel
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
setMethod("fit", signature = "cdtamodel",
          function(object,
                   data,
                   SID,
                   cores=3,
                   chains=3,
                   iter=6000,
                   warmup=1000,
                   thin=10,
                   ...){
              fit.cdtamodel(cdtamodel=object,
                            data=data,
                            SID=SID,
                            cores=cores,
                            chains=chains,
                            iter=iter,
                            warmup=warmup,
                            thin=thin,
                            ...)
          })


#' Print a cdtamodel object.
#' @rdname show
#' @param object A cdtamodel object returned by \link{cdtamodel} function.
#' @method show cdtamodel
#' @keywords internal
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
setMethod("show", signature = "cdtamodel",
          function(object){
              cat(text=object@modelcode)
          })



#' Print a cdtafit object.
#' @rdname show
#' @method show cdtafit
#' @keywords internal
setMethod("show", signature = "cdtafit",
          function(object){
              print.cdtafit(object)
})




#' Forestplots for cdtafit objects.
#' @rdname plot
#' @param x A cdtafit object from \link{fit}.
#' @param graph An optional numeric value indicating which forest to plot(s) to graph. Valid values are:0 - for no graph, 1 - yielding a forest plot of the
#' sensitivity and specificity with a 95 percent exact confidence intervals, 2 - yielding a forest plot of the posterior study-specific sensitivity and specificity
#' and the marginal mean sensitivy and specificity and their corresponding 95 percent credible intervals, 3 - yielding a combination of 1 and 2 in one plot, and NULL(default) - yielding plots of
#' 1, 2 and 3.
#' @param title.1 An optional string indicating the title of graph 1.
#' @param title.2 An optional string indicating the title of graph 2.
#' @param title.3 An optional string indicating the title of graph 3.
#' @param width An optional numeric value to adjust the dogding position. The default is 0.2.
#' @param shape.1 An optional numeric value(0-255) indicating the symbol to plot in graph 1. The defualt is 19 which is a solid circle. See \link[graphics]{points} for more details.
#' @param size.1 An optional positive numeric value indicating the size of symbols in graph 1. The defualt is 2.5.
#' @param shape.2 An optional numeric value(0-255) indicating the symbol to plot in graph 2. The defualt is 8 which is a star. See \link[graphics]{points} for more details.
#' @param size.2 An optional positive numeric value indicating the size of symbols in graph 2. The defualt is 2.5.
#' @param shape.O An optional numeric value(0-255) indicating the symbol representing the posterior marginal mean in graph 2. The defualt is 19 which is a solid circle. See \link[graphics]{points} for more details.
#' @param size.O An optional numeric value indicating the size of symbols representing the posterior marginal means in graph 2.
#' @param cols.1 An optional string vector specifying colours of shapes in graph 1.
#' @param cols.2 An optional string vector specifying colours of shapes in graph 2.
#' @param digits An optional positive value to control the number of digits to print when printing numeric values. The default is 3.
#' @param ... Extra optional arguments. See \link[rstan]{stan} options.
#' @method plot cdtafit
 #'@importFrom graphics plot
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}

setMethod("plot", signature(x="cdtafit", y="missing"),
          function(x,
                   title.1=NULL,
                   title.2=NULL,
                   title.3=NULL,
                   graph=NULL,
                   width=0.2,
                   shape.1=19,
                   size.1=2.5,
                   shape.2=8,
                   size.2=2.5,
                   shape.O=9,
                   size.O=3.5,
                   cols.1=NULL,
                   cols.2=NULL,
                   digits=3,...){
              forestplot.cdtafit(x,
                             title.1=NULL,
                             title.2=NULL,
                             title.3=NULL,
                             graph=NULL,
                             width=0.2,
                             shape.1=19,
                             size.1=2.5,
                             shape.2=8,
                             size.2=2.5,
                             shape.O=9,
                             size.O=3.5,
                             cols.1=NULL,
                             cols.2=NULL,
                             digits=3,...)
})

#' Traceplots for cdtafit objects.
#' @param object A cdtafit object from \link{fit}
#' @param ... Extra optional arguments as defined in \link[rstan]{stan_trace}.
#' @rdname traceplot
setGeneric(name="traceplot", def=function(object, ...){standardGeneric("traceplot") })

#' @rdname  traceplot
#' @method traceplot cdtafit
#' @export
#' @author Victoria N Nyaga \email{victoria.nyaga@outlook.com}
setMethod("traceplot", signature = "cdtafit",
          function(object, ...){
              traceplot.cdtafit(x=object, ...)

})

