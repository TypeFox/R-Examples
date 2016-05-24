DanielPlot <- function(fit, ...){
    UseMethod("DanielPlot")
}
DanielPlot.design <- function(fit, ..., response=NULL){
    if (!"design" %in% class(fit)) 
        stop("DanielPlot.design works for obj from class design only.")
    di <- design.info(fit)
    if (is.null(di$response.names)) 
        stop("The design fit must have at least one response.")
    if (!(is.null(response))) 
      if (!response %in% di$response.names)
        stop("Requested response is not a response variable in fit.")
    if (!(length(grep("FrF2",di$type))>0 | 
           length(grep("pb",di$type))>0)) { 
           if (!(di$type=="full factorial" & all(di$nlevels==2)))
        stop("The design obj must be of a type containing FrF2 or pb.")
       }
    grad <- 1
    if (length(grep("pb",di$type)) > 0 & di$nfactors < di$nruns-1)
          warning("Effects plots for Plackett-Burman designs must be done with nruns-1 effects! The error effects are missing!")
    ## make sure there are as many effects as possible in the plots, redundant ones will not be shown
    if (length(grep("FrF2",di$type)) > 0 | (di$type=="full factorial" & all(di$nlevels==2))){
       grad <- 2
       hilf <- lm(fit, degree=grad)
       ncoef <- sum(!is.na(coef(hilf)))
       nhilf <- di$nrun
       if (!is.null(di$ncenter)) nhilf <- nhilf - di$ncenter
       while (ncoef < nhilf){
          grad <- grad+1
          hilf <- lm(fit, degree=grad)
          ncoef <- sum(!is.na(coef(hilf)))
          }
    }
    subtext <- ""
    if (length(grep("splitplot", di$type)) > 0){
        ## splitplot situation
        ## make sure to use different plot symbols for whole plot and split plot effects
        ## important to distinguish wbreps and bbreps
               mm <- model.matrix(hilf)
               coefs <- coef(hilf)[-1] ## only non-missing coefficients
               nms <- names(coefs)
               hilf <- apply(mm[,1+(1:di$nfac.WP),drop=FALSE],1,paste,collapse="")
               pchs <- rep("*", length(coefs))
               pchs[1:di$nfac.WP] <- "o"
               
               for (j in setdiff(1:(length(coefs)),1:di$nfac.WP)){
                  if (!length(table(paste(hilf,mm[,nms[j]],sep="")))>di$nWPs) pchs[j] <- "o"
                }
          subtext <- "WARNING: whole plot effects (marked by o) may have larger variation than split-plot effects"
          DanielPlot(lm(fit, degree=grad, response=response), pch=pchs, subtitle=subtext, ...)
          }
          else DanielPlot(lm(fit, degree=grad, response=response), ...)
}

DanielPlot.default <-
function (fit, code = FALSE, autolab = TRUE, alpha=0.05,
         faclab = NULL, 
         block = FALSE, datax = TRUE, 
         half = FALSE, pch = "*", 
         cex.fac = par("cex.lab"), cex.lab = par("cex.lab"), 
         cex.pch = par("cex"), cex.legend = par("cex.lab"), 
         main = NULL, subtitle=NULL, ...) 
{
   if (! ("lm" %in% class(fit) | "aov" %in% class(fit))) 
      stop("fit must be a linear model object (lm or aov), or a design of class design")
    ## transform into -1 and 1 coded model
    fit <- remodel(fit)$model
    ## check whether of appropriate type
    if (!check(fit))
     stop("This routine is applicable for 2-level factorial designs without partial aliasing only.")
   
    
    if (any(names(coef(fit)) == "(Intercept)")) {
        factor.effects <- 2 * coef(fit)[-1]
        }
    else {
        factor.effects <- 2 * coef(fit)
    }
    respnam <- colnames(fit$model)[attr(attr(fit$model,"terms"),"response")]
    names(factor.effects) <- attr(fit$terms, "term.labels")
    terms.ord <- attr(fit$terms, "order")[!is.na(factor.effects)] ##moved here
    factor.effects <- factor.effects[!is.na(factor.effects)]
    plotmain <- paste("Normal Plot for", respnam)
    if (autolab) {
        ### take simulated critical values from package DoE.base, if available
        crit <- ME.Lenth(factor.effects,alpha=alpha)$ME
        if (!code)
        faclab <- list(idx = which(crit<=abs(factor.effects)),
               lab = names(factor.effects)[which(crit<=abs(factor.effects))])
        plotmain <- paste(plotmain, ", ", "alpha=", alpha, sep="") 
        }
    if (half) {
        n <- length(factor.effects)
        ## qnorm(0.5 + ppoints(n, a=1/2)/2) in halfnormal
        tn <- list(x = qnorm(0.5 + ppoints(n, a=1/2)/2)[rank(abs(factor.effects))], 
            y = abs(factor.effects))
        xlab <- "half-normal scores"
        ylab <- "absolute effects"
        plotmain <- paste("Half", plotmain)
    }
    else {
        tn <- qqnorm(factor.effects, plot = FALSE)
        xlab <- "normal scores"
        ylab <- "effects"
    }
    names(tn$x) <- names(factor.effects)  ## moved here
    names(tn$y) <- names(factor.effects)  ## added
    if (datax) {
        tmp <- tn$x
        tn$x <- tn$y
        tn$y <- tmp
        tmp <- xlab
        xlab <- ylab
        ylab <- tmp
    }
    labx <- names(factor.effects)
    laby <- 1:length(tn$y)
    points.labels <- names(factor.effects)
    if (is.null(main)) main <- plotmain
    plot.default(tn, xlim = c(min(tn$x), max(tn$x) + diff(range(tn$x))/5), 
        pch = pch, xlab = xlab, ylab = ylab, cex=cex.pch, cex.lab = cex.lab, 
        mgp=c(2,1,0), main = main, ...)
    ## at the top below main title, mainly for warning in case of split-plot
    if (!is.null(subtitle)) mtext(subtitle)
    if (is.null(faclab)) {
        if (!code) {
            effect.code <- labx
        }
        else {
            max.order <- max(terms.ord)
            no.factors <- length(terms.ord[terms.ord == 1])
            factor.label <- attr(fit$terms, "term.labels")[terms.ord == 1]
            faclet <- c(LETTERS[-9],letters[-9])
            factor.code <- faclet[1:no.factors]
            if (block) 
                factor.code <- c("BK", factor.code)
            texto <- paste(factor.code[1], "=", factor.label[1])
            for (i in 2:no.factors) {
                texto <- paste(texto, ", ", factor.code[i], "=", 
                  factor.label[i])
            }
            mtext(side = 1, line = 3.5, texto, cex = cex.legend)
            get.sep <- function(string, max.order) {
                k <- max.order - 1
                get.sep <- rep(0, k)
                j <- 1
                for (i in 2:(nchar(string)-1)) {
                  if (substring(string, i, i) == ":") {
                    get.sep[j] <- i
                    if (j == k) 
                      break
                    j <- j + 1
                  }
                }
                get.sep
            }
            labeling <- function(string, get.sep, max.order, 
                factor.code, factor.label) {
                labeling <- ""
                sep <- get.sep(string, max.order)
                sep <- sep[sep > 0]
                n <- length(sep) + 1
                if (n > 1) {
                  sep <- c(0, sep, nchar(string) + 1)
                  for (i in 1:n) {
                    labeling <- paste(labeling, sep = "", factor.code[factor.label == 
                      substring(string, sep[i] + 1, sep[i + 1] - 
                        1)][1])
                  }
                }
                else labeling <- paste(labeling, sep = "", factor.code[factor.label == 
                  string][1])
                labeling
            }
            effect.code <- rep("", length(terms.ord))
            for (i in 1:length(terms.ord)) {
                effect.code[i] <- labeling(names(tn$x)[i], get.sep, 
                  max.order, factor.code, factor.label)
            }
        }
        if (autolab){ 
               faclab <- list(idx = which(crit<=abs(factor.effects)),
                  lab = effect.code[which(crit<=abs(factor.effects))])
               if (length(faclab$idx) > 0)
               text(as.data.frame(tn)[faclab$idx,], paste(" ", faclab$lab), cex = cex.fac, adj = 0, 
                  xpd = NA)
           }
        else
        text(tn, paste("   ", effect.code), cex = cex.fac, adj = 0, 
            xpd = NA)
    }
    else {
        if (!is.list(faclab)) 
            stop("* Argument 'faclab' has to be NULL or a list with idx and lab objects")
        if (length(faclab$lab)>0) text(tn$x[faclab$idx], tn$y[faclab$idx], labels = faclab$lab, 
            cex = cex.fac, adj = 0)
    }
    if (!length(pch)==length(factor.effects)) pchs <- rep(pch, length(factor.effects))
    if (code) aus <- cbind(as.data.frame(tn), no = 1:length(tn$x), effect=names(factor.effects), coded=effect.code, pchs=pch)
      else aus <- cbind(as.data.frame(tn), no = 1:length(tn$x), effect=names(factor.effects), pchs=pch)
    invisible(aus)
}

#removed August 15 2013
#halfnormal <- function(effects, labs, codes=labs, alpha=0.05, xlab="absolute effects", ...){
#    effects <- abs(effects)
#    labord <- order(effects)
#    effects <- sort(effects)
#    if (!identical(codes, labs)){ 
#         haupteff <- setdiff(1:length(labs), grep(":", labs)) 
#         legende <- paste(codes[haupteff], labs[haupteff], sep="=",collapse=", ")
#     }
#     else legende <- ""
#    n <- length(effects)
#    ui <- qnorm(0.5 + (0:(n-1)+0.5)/(2*n))  ## "+0.5" in parentheses added Aug 14 2013
#    codes <- paste(rep("  ",n), codes, sep="")
#    
#    crit <- LenthPlot(effects,alpha=alpha,plt=FALSE)["ME"]
#    nlab <- sum(effects>crit)
#    plot(effects, ui, ylab = "Half-normal scores", xlab = xlab, sub=legende, ...)
#    if (nlab < n) 
#        points(effects[1:(n - nlab)],ui[1:(n - nlab)])
#    text(effects[(n - nlab + 1):n], ui[(n - nlab + 1):n], codes[labord][(n - 
#        nlab + 1):n], adj=0, xpd=NA)
#}

