##' Class "stability.path"
##'
##' Class of object returned by the \code{stability} function, with
##' methods \code{print}, \code{show} and \code{plot}.
##'
##' @section Slots: \describe{
##'
##' \item{\code{probabilities}: }{a \code{Matrix} object containing the
##' estimated probabilities of selection along the path of solutions.}
##' \item{\code{penalty}: }{Object of class \code{"character"}
##' indicating the penalizer used.}
##' \item{\code{naive}: }{logical indicating whether rescaling of the
##' coefficients has been performed regarding the \eqn{\ell_2}{l2}-penalty.}
##' \item{\code{lambda1}: }{a vector with the levels of the first penalty.}
##' \item{\code{lambda2}: }{a scalar with the \eqn{\ell_2}{l2}-penalty level.}
##' \item{\code{folds}: }{a list that contains the folds used for each subsample.}
##' }
##' @aliases print,stability.path-method show,stability.path-method
##'
##' @docType class
##'
##' @keywords class
##'
##' @seealso See also \code{\link{plot,stability.path-method}}, and
##' \code{\link{stability}}.
##'
##' @name stability.path-class
##'
##' @rdname stability.path-class
##'
##' @exportClass stability.path
##' @exportMethod print
##' @exportMethod show
##'
setClass("stability.path",
  representation = representation(
     probabilities = "Matrix"   ,
     penalty       = "character",
     naive         = "logical"  ,
     lambda1       = "numeric"  ,
     lambda2       = "numeric"  ,
     folds         = "list"  )
)

setMethod("print", "stability.path", definition =
   function(x, ...) {
     if (x@naive) {
       cat("Stability path for", x@penalty, "penalizer, no rescaling of the coefficients (naive).\n")
     } else {
       cat("Stability path for", x@penalty, "penalizer, coefficients rescaled by (1+lambda2).\n")
     }
     cat("- penalty parameter lambda1:", length(x@lambda1), "points from",
         format(max(x@lambda1), digits = 3),"to",
         format(min(x@lambda1), digits = 3),"\n")
     cat("- penalty parameter lambda2:", x@lambda2)
     cat("\n")
     invisible(x)
   }
)

setMethod("show", "stability.path", definition =
   function(object) {print(object)}
)

##' Plot method for \code{stability.path}.
##'
##' Produce a plot of the stability path obtained by stability
##' selection.
##'
##' @usage \\S4method{plot}{stability.path}(x, y, xvar = "lambda", annot=TRUE,
##'          main = paste("Stability path for ", slot(x, "penalty")," regularizer", sep=""),
##'          log.scale = TRUE,  labels = rep("unknown status",p), plot = TRUE,
##'          sel.mode = c("rank","PFER"), cutoff=0.75, PFER=2, nvar=floor(n/log(p)), ...)
##' @param x output of a \code{stability} run (must be of class
##' \code{stability.path}).
##' @param y used for S4 compatibility.
##' @param main main title. If none given, a somewhat appropriate
##' title is automatically generated.
##' @param xvar variable to plot on the X-axis: either \code{"lambda"}
##' (first penalty level) or \code{"fraction"} (fraction of the
##' penalty level applied tune by \eqn{\lambda_1}{lambda1}). Default
##' is \code{"lambda"}.
##' @param log.scale logical; indicates if a log-scale should be used
##' when \code{xvar="lambda"}. Default is \code{TRUE}.
##' @param plot logical; indicates if the graph should be
##' plotted. Default is \code{TRUE}. If \code{FALSE}, only the
##' \pkg{ggplot2} object is sent back.
##' @param sel.mode a character string, either \code{'rank'} or
##' \code{'PFER'}. In the first case, the selection is based on the
##' rank of total probabilties by variables along the path: the first
##' \code{nvar} variables are selected (see below). In the second
##' case, the PFER control is used as described in Meinshausen and
##' Buhlmannn's paper. Default is \code{'rank'}.
##' @param nvar number of variables selected (only relevant when
##' \code{sel.mode} equals \code{'rank'}. Default is \code{floor(n/log(p))}.
##' @param cutoff value of the cutoff probability (only relevant when
##' \code{sel.mode} equals \code{'PFER'}).
##' @param PFER value of the per-family error rate to control (only
##' relevant when \code{sel.mode} equals \code{'PFER'}).
##' @param labels an optional vector of labels for each variable in
##' the path (e.g., 'relevant'/'irrelevant'). See examples.
##' @param annot logical; should annotation be made on the graph
##' regarding controlled PFER (only relevant when \code{sel.mode}
##' equals \code{'PFER'})? Default is \code{TRUE}.
##' @param ... used for S4 compatibility.
##' @return a list with a \pkg{ggplot2} object which can be plotted
##' via the \code{print} method, and a vector of selected variables
##' corresponding to method of choice (\code{'rank'} or
##' \code{'PFER'})
##' 
##' @name plot,stability.path-method
##' @aliases plot,stability.path-method
##' @aliases plot.stability.path
##' @docType methods
##' @rdname plot.stability.path
##'
##'
##' @examples \dontrun{
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
##' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
##' diag(Sigma) <- 1
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' ## Build a vector of label for true nonzeros
##' labels <- rep("irrelevant", length(beta))
##' labels[beta != 0] <- c("relevant")
##' labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))
##'
##' ## Call to stability selection function, 200 subsampling
##' stab <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)
##' 
##' ## Build the plot an recover the selected variable
##' plot(stab, labels=labels)
##' plot(stab, xvar="fraction", labels=labels, sel.mode="PFER", cutoff=0.75, PFER=2)
##' }
##' @importFrom graphics plot
##' @exportMethod plot
##' @import ggplot2 scales grid
##' @export
setMethod("plot", "stability.path", definition =
   function(x, y, xvar = "lambda", annot=TRUE,
            main = paste("Stability path for ", slot(x, "penalty")," regularizer", sep=""),
            log.scale = TRUE,  labels = rep("unknown status",p), plot = TRUE,
            sel.mode = c("rank","PFER"), cutoff=0.75, PFER=2, nvar=floor(n/log(p)), ...) {

     p <- ncol(x@probabilities)
     n <- max(unlist(x@folds))
     sel.mode <- match.arg(sel.mode)
     if (length(x@lambda1) == 1)
       stop("Not available when length(lambda1) == 1")

     if(PFER <=0)
       stop("PFER should be at least equal to 1.")

     if(cutoff <0.5 | cutoff > 1)
       stop("The cutoff is supposed to be a probability in [.5,1] ...")

     nzeros <- which(colSums(x@probabilities) != 0)
     if(nvar>p)
       stop("The rank is supposed to be less than p ...")

     if (length(nzeros) == 0)
       stop("Nothing to plot: all probabilities are zero along the path.")

     prob  <- as.matrix(x@probabilities[, nzeros])
     rownames(prob) <- NULL

     selection <- rep("unselected",ncol(prob))
     if (sel.mode == "PFER") {
       ## estimate the average number of selected variables on the
       ## current path and pick the one controling the PFER at the
       ## desired level
       q    <- rowSums(prob >= cutoff)
       qLim <- sqrt(PFER * (2 * cutoff-1) * p)
       iq <- q <= qLim
       iq <- ifelse(which.min(iq) != 1,which.min(iq)-1,ifelse(sum(iq) == 0,1,length(iq)))
       selection[prob[iq, ] > cutoff] <- "selected"
     } else {
       selection[order(colSums(prob),decreasing=TRUE)[1:nvar]] <- "selected"
     }

     ## the x-axis variable
     xv <- switch(xvar,"fraction" = x@lambda1/max(x@lambda1), x@lambda1)
     if (xvar == "lambda")
       xv <- log10(xv)

     ## Build the data frame for ggploting
     data.coef <- melt(data.frame(xvar=xv, prob=prob),id="xvar")
     data.coef$selection <- factor(rep(selection, each=length(xv)))
     if (is.null(labels)) {
       data.coef$labels <- factor(rep(1:p, each=length(xv)))
     } else {
       if (sum(is.na(labels[nzeros]))>0 ) {
         labels <- NULL
         warning("The number of label is wrong, ignoring them.")
         data.coef$labels <- factor(rep(nzeros, each=length(xv)))
       } else {
         data.coef$labels <- factor(rep(labels[nzeros], each=length(xv)))
       }
     }
     colnames(data.coef) <- c("xvar","var","prob","selection","variables")

     ## Build the ggplot object
     d <- ggplot(data.coef,aes(x=xvar,y=prob, linetype=variables, colour=selection, group=var)) +
       geom_line(aes(x=xvar,y=prob)) +
         labs(x=switch(xvar,
                "fraction" = expression(lambda[1]/max[lambda[1]]),
                ifelse(log.scale,expression(log[10](lambda[1])),expression(lambda[1]))),
              y="selection probabilities") + ggtitle(main)
     d <- d + scale_x_reverse()
     if (is.null(labels)) {
       d <- d + theme(legend.position="none")
     } else {
       if (length(labels[nzeros]) != length(nzeros)) {
         d <- d + theme(legend.position="none")
       }
     }

     ## Manage the annotation for the selected variables (PFER mode)
     if (annot & sel.mode == "PFER") {
       d <- d + annotate("rect", xmin=xv[1], xmax=xv[iq], ymin=cutoff, ymax=1, alpha=0.15)
       d <- d + annotate("text", x=c(xv[length(xv)],xv[1],xv[iq]), y=c(cutoff,1,0), hjust=c(0,0,0.25), vjust=c(0,1.1,1.1),
                         label=c(paste("pi[thr]"),paste("PFER <=",PFER),paste("hat(q)==",round(q[iq],2))),
                         parse=TRUE, size=4, alpha=0.85)
       d <- d + geom_hline(yintercept=cutoff, linetype="dashed", alpha=0.35, size=.5)
       d <- d + geom_vline(xintercept=xv[iq], linetype="dashed", alpha=0.35, size=.5)
     }

     if (plot) {print(d)}
     invisible(list(ggplot.object=d, selected=nzeros[selection == "selected"]))

   }
)
