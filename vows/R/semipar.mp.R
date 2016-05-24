#' Massively parallel semiparametric regression
#' 
#' Fits a possibly very large number of semiparametric models by quadratically
#' penalized least squares.  The model may include a combination of parametric
#' terms, smooth terms, varying-coefficient terms, and simple random effect
#' structures.
#' 
#' The basic approach to massively parallel smoothing is described in Reiss et
#' al. (2014). Although simple mixed-effect models are available,
#' \code{\link{semipar.mix.mp}} is generally preferable for mixed models with a
#' single smooth term.
#' 
#' Each element of \code{list.all} corresponding to a \emph{nonparametric} term
#' of the model is a list with components \code{modmat}, \code{penmat},
#' \code{pen.order}, \code{start}, and \code{end}. For each \emph{parametric}
#' term, the same five components are included, plus \code{basis},
#' \code{argvals}, \code{effect}, \code{k}, and \code{norder}.
#' 
#' @param formula a formula object such as "\code{~ x1 + sf(x2) +sf(x2, effect
#' = x3)}" where \code{x1} is a linear (parametric) predictor, \code{x2} is a
#' predictor on which the responses depend smoothly, and \code{x3} is a
#' predictor whose effect is linear but varies smoothly with \code{x2} (i.e., a
#' varying-coefficient predictor).
#' @param Y an \eqn{n \times V} response matrix, where \eqn{V} is the number of
#' models fitted in parallel, e.g., voxels in neuroimaging applications.
#' @param lsp vector of candidate log tuning parameters (\eqn{log(\lambda)}).
#' @param data an optional data frame containing the variables in the model.
#' @param range.basis a numeric vector of length 2 defining the interval over
#' which the B-spline basis is created. If \code{NULL}, it will be set as the
#' range of the variable to be evaluated by the basis.
#' @param knots knot placement for the B-spline bases. The default,
#' \code{"quantile"}, gives knots at equally spaced quantiles of the data. The
#' alternative, \code{"equispaced"}, gives equally spaced knots.
#' @param rm.constr logical: should the constraints be removed for
#' varying-coefficient models?
#' @param random a formula or a matrix for random effects.
#' @param store.reml logical: should the pointwise REML criterion at each grid
#' point be included in the output?  \code{FALSE} by default, as this output
#' can be very large.
#' @param store.fitted logical: should the fitted values be included in the
#' output? \code{FALSE} by default.
#' @return An object of class \code{"semipar.mp"}, which is also of class
#' \code{"\link{qplsc.mp}"} but includes the following additional elements:
#' \item{where.sf, where.nsf}{vectors or scalars identifying where the smooth
#' and non-smooth terms, respectively, appear in the model formula.}
#' \item{list.all}{a list of lists, one for each term of the model; see
#' Details.} \item{formula}{model formula.} \item{Y}{response matrix.}
#' \item{lsp}{candidate values for the log smoothing parameter.}
#' \item{data}{the supplied data frame, if any.}
#' @author Yin-Hsiu Chen \email{enjoychen0701@@gmail.com} and Philip Reiss
#' \email{phil.reiss@@nyumc.org}
#' @references Reiss, P. T., Huang, L., Chen, Y.-H., Huo, L., Tarpey, T., and
#' Mennes, M. (2014). Massively parallel nonparametric regression, with an
#' application to developmental brain mapping. \emph{Journal of Computational
#' and Graphical Statistics}, \emph{Journal of Computational and Graphical
#' Statistics}, 23(1), 232--248.
#' @examples
#' 
#' n<-32
#' Ys <- matrix(0, n, 5)
#' for(i in 1:n) Ys[i,]<--2:2+rnorm(5, i^2, i^0.5)+sin(i)
#' x1 <- rnorm(n,0,5)
#' x2 <- 1:n+runif(n, 1, 20)
#' semipar.obj <- semipar.mp(~x1+sf(x2,k=10),Y=Ys,lsp=seq(5,50,,30))
#' @export
semipar.mp <- function(formula, Y, lsp, data = NULL, range.basis = NULL, knots = "quantile", rm.constr = FALSE, random = NULL, store.reml = FALSE, store.fitted = FALSE) {
    tf <- terms.formula(formula, specials="sf")
    trmstrings <- attr(tf, "term.labels")
    envir <- environment(formula)
    if (!is.null(range.basis)) envir$range.basis = range.basis
    envir$quantile = quantile
    where.sf <- attr(tf, "specials")$sf
    where.nsf <- which((1:length(trmstrings))%in%where.sf==FALSE)
    list.all <- vector("list", length(trmstrings))
    
    if (length(where.nsf)!=0) {
        for (i in 1:length(where.nsf)) {
            list.all[[where.nsf[i]]]$modmat <- if (!is.null(data)) as.matrix(eval(as.call(parse(text=trmstrings[where.nsf[i]]))[[1]],envir = data, enclos = envir)) else 
as.matrix(eval(as.call(parse(text=trmstrings[where.nsf[i]]))[[1]], envir = envir))        	
            ncol <- ncol(list.all[[where.nsf[i]]]$modmat)
            list.all[[where.nsf[i]]]$penmat <- matrix(0, ncol, ncol) 
            list.all[[where.nsf[i]]]$pen.order <- ncol
        }
    }
    if (length(where.sf)!=0) {
        if (is.null(data)) for (i in 1:length(where.sf)) {
                temptext = trmstrings[where.sf[i]]
                if (!is.null(range.basis)) temptext <- sub("\\)", ", range.basis = range.basis\\)", temptext) 
                list.all[[where.sf[i]]]<-eval(parse(text = temptext), envir = envir)
        }
        else for (i in 1:length(where.sf)) {
            formula.term <- sub("sf\\(", "", sub("\\)", "", trmstrings[where.sf[i]]))
            split.term <- as.vector(strsplit(formula.term, ",")[[1]])
            var.x <- split.term[1]
            text <- sub(var.x,paste("data$",var.x,sep=""),trmstrings[where.sf[i]])
            if (!is.na(stringr::str_locate(formula.term, ", effect")[1])) { 
    	        effect.ind <- which(is.na(stringr::str_locate(split.term, " effect = ")[,1])==FALSE)
                var.effect <- sub(" effect = ", "", split.term[effect.ind])
                text <- sub(var.effect,paste("data$",var.effect,sep=""),text)
            }
            list.all[[where.sf[i]]]<-eval(parse(text = text))
        }
    } 

    modmat <- do.call('cbind', lapply(list.all, function(xx) xx$modmat))
    subpen.dim <- lapply(list.all, function(xx) dim(xx$penmat)[1])
    penmat <- matrix(0, do.call('sum', subpen.dim), do.call('sum', subpen.dim))
    for (i in 1:length(subpen.dim)) {
        start.ind <- cumsum(subpen.dim)[i] - subpen.dim[[i]] + 1
        end.ind <- cumsum(subpen.dim)[i]
        penmat[start.ind:end.ind, start.ind:end.ind] <- list.all[[i]]$penmat
	    list.all[[i]]$start <- start.ind + 1*(attr(tf, "intercept")==1)
	    list.all[[i]]$end <- end.ind + 1*(attr(tf, "intercept")==1)
    }  
    if (!is.null(random)) {
        if (is.matrix(random)) {
            modmat <- cbind(modmat, random)
            penmat <- rbind(cbind(penmat, matrix(0, NCOL(penmat), NCOL(random))), cbind(matrix(0, NCOL(random), NCOL(penmat)), diag(1, NCOL(random))))
        }
        else {
            ran.string <- attr(terms.formula(random), "term.labels")
            ran.terms <- strsplit(ran.string, " \\| ")
            factor <- if (!is.null(data)) as.factor(eval(as.call(parse(text=ran.terms[[1]][2]))[[1]] ,envir = data, enclos = envir))  
                      else as.factor(eval(as.call(parse(text=ran.terms[[1]][2]))[[1]], envir = envir))

          modmat <- cbind(modmat, model.matrix(~ factor-1))
          penmat <- rbind(cbind(penmat, matrix(0, NCOL(penmat), length(levels(factor)))), cbind(matrix(0, length(levels(factor)), NCOL(penmat)), diag(1, length(levels(factor)))))
        }      
    }

    constr.list = NULL
    if (!rm.constr) {
        ind.constr <- do.call('c',lapply(list.all,function(xx) is.null(xx$effect) && !is.null(xx$basis)))
        if (sum(ind.constr)!=0) {
            constr.list <- vector("list", sum(ind.constr))
            for (i in 1:sum(ind.constr)) {
                start <- cumsum(subpen.dim)[which(ind.constr)[i]] - subpen.dim[[which(ind.constr)[i]]] + 1
                end <- cumsum(subpen.dim)[which(ind.constr)[i]]
                constr.list[[i]]$start <- if (attr(tf, "intercept")==1) start+1  else start
                constr.list[[i]]$end <- if (attr(tf, "intercept")==1) end+1  else end
                constr.list[[i]]$C <- matrix(colSums(modmat[,start:end]),1)            
            }
        }
    }
    nulldim <- do.call('sum', lapply(list.all, function(xx) xx$pen.order)) 
    if (attr(tf, "intercept")==1) {
        modmat <- cbind(1, modmat)
        penmat <- rbind(0, cbind(0, penmat))
        nulldim <- nulldim + 1
    }
    qplsc.obj <- qplsc.mp(Y = Y, modmat = modmat, penmat = penmat, constr.list = constr.list, lsp = lsp, nulldim = nulldim, store.reml = store.reml, store.fitted = store.fitted)
    qplsc.obj$where.sf <- where.sf; qplsc.obj$where.nsf <- where.nsf
    qplsc.obj$list.all <- list.all; qplsc.obj$formula <- formula
    qplsc.obj$Y <- Y; qplsc.obj$lsp <- lsp; qplsc.obj$data <- data
    class(qplsc.obj) = c("semipar.mp", "qplsc.mp")
    qplsc.obj
}
