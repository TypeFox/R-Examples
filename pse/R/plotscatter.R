#' @export
#' @rdname plots
plotscatter <-
	function (obj, res=NULL, index.data=NULL, index.res=NULL, add.lm=TRUE, ylab = NULL, ...) {
		if (class(obj)=="LHS" | class(obj)=="PLUE") {
			if (is.null(index.data)) index.data <- 1:get.ninputs(obj)
			if (is.null(index.res)) index.res <- 1:get.noutputs(obj)
			data=get.data(obj)
			results=get.results(obj)
			if(is.null(ylab)) ylab = obj$res.names[index.res]
			plotscatter(data[,index.data], results[,index.res], add.lm=add.lm, ylab = ylab, ...)
		} else {
			res <- as.data.frame(res)
			if (is.null(dim(obj))) {nplots <-1} else {nplots <- dim(obj)[2]}
			if (is.null(dim(res))) {nplots <-nplots*1} else {nplots <- nplots* dim(res)[2]}

			nl <- floor(sqrt(nplots))
			nc <- ceiling(nplots/nl)
			opar <- par(no.readonly=TRUE)
			par (mar=c(4,4,2,0.5), mfrow=c(nl, nc), pch='.')
			on.exit(par(opar))

			index.var <- 1
			index.res <- 1
			for (i in 1:nl) for (j in 1:nc) {
				this.ylab =names(res)[index.res]
				if (!is.null(ylab)) {
					if (length(ylab)==1) this.ylab=ylab
					else this.ylab=ylab[index.res]
				}
				
				oneplotscatter(res[,index.res],obj[, index.var], c(names(obj)[index.var],this.ylab), add.lm, ...)
				index.res <- index.res +1
				if (index.res > dim(res)[2]) {index.res <- 1; index.var <- index.var + 1}
				if (index.var > dim(obj)[2]) break;
			}
		}
	}

oneplotscatter <-
function(res, var, name, add.lm, ...) {
		plot(var,res, xlab=name[1], ylab= name[2], ...)
		if (add.lm) {
			l <- lm(res~var)
			if (!is.na(coefficients(l)[2])) 
					abline((lm(res~var)))
		}
}
