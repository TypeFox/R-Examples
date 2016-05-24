MEPlot <- function(obj, ...){
    UseMethod("MEPlot")
}
MEPlot.design <- function(obj, ..., response=NULL){
    if (!"design" %in% class(obj)) 
        stop("MEPlot.design works for obj from class design only.")
    di <- design.info(obj)
    if (is.null(di$response.names)) 
        stop("The design obj must have at least one response.")
    if (!(is.null(response))) 
      if (!response %in% di$response.names)
        stop("Requested response is not a response variable in fit.")
    if (!(length(grep("FrF2",di$type))>0 | 
           length(grep("pb",di$type))>0)){ 
           if (!(di$type=="full factorial" & all(di$nlevels==2)))
           stop("The design obj must be of a type containing FrF2 or pb.")
    }
    MEPlot(lm(obj, degree=1, response=response), ...)
}
MEPlot.default <-
function(obj, main=paste("Main effects plot for", respnam), pch=15, 
         cex.xax = par("cex.axis"), cex.yax = cex.xax, mgp.ylab = 4, 
         cex.title=1.5, cex.main=par("cex.main"), lwd=par("lwd"), 
         abbrev=3, select=NULL, ...){
   # main      overall title
   # pch       plot character
   # lwd       line width
   # cex.xax   cex for levels on x-axes 
   # cex.yax   cex for label and levels on y-axis
   # mgp.ylab  mgp entry y-axis label 
   #           (relative to invisible axis of left-most plot)
   # cex.title multiplier for cex.main for the overall title given in option main
   # abbrev  maximum number of characters used for levels on x-axes
   if (! ("lm" %in% class(obj) | "aov" %in% class(obj))) 
      stop("obj must be a linear model object (lm or aov), or a design of class design")
   obj <- remodel(obj)
   labs <- lapply(obj$labs,function(sp) substr(sp,1,abbrev))
   obj <- obj$model
   if (!check(obj)) 
     stop("This routine is applicable for 2-level factorial designs without partial aliasing only.")
   mod <- obj$model
   term.ord <- attr(terms(obj),"order")
   nmain <- length(which(term.ord==1))
   intcol <- attr(attr(mod,"terms"),"intercept")
   respnam <- colnames(mod)[attr(attr(mod,"terms"),"response")]
   ymean <- mean(mod[,respnam])
   mm <- model.matrix(obj)
   if (intcol > 0) mm <- mm[,-intcol]
   terms1 <- colnames(mm)[which(term.ord==1)]
   if (is.null(select)) select <- 1:nmain
   else {
       if (!is.numeric(select)) stop("select must be numeric")
       if (!all(floor(select)==select)) stop("select must contain integer numbers")
       if (any(select<1 | select>nmain)) stop("select must contain numbers betweeen 1 and ", nmain, " only")
   }
   predmat <- matrix(rep(0,2*nmain),2,nmain)
   colnames(predmat) <- terms1
   predmat <- predmat[,select]
   nmain <- length(select)
   terms1 <- terms1[select]
   labs <- labs[select]
   #addnam <- setdiff(colnames(obj$model), terms1)
   #names <- c(terms1,addnam)
   for (i in 1:nmain)
      predmat[,i] <- ymean+c(-1,1)*coef(obj)[terms1[i]]
   omfrow <- par("mfrow")   
   omar <- par("mar")
   ooma <- par("oma")
   ax <- pretty(c(min(predmat),max(predmat)))
   par(mfrow=c(1,nmain),mar=c(2, 1, 2, 1) + 0.1, oma=c(3,5,4,0.1))
   for (i in 1:nmain){
          ## plot effect without axis drawn (but y label shown for i==1)
          ## plot effect without axis drawn and without axis labels for i>1
          if (i==1)
          plot(c(-1,1),predmat[,i],main=terms1[i],xlab="",xpd=NA, ylab=respnam,type="b",
              xlim=c(-1.3,1.3),ylim=c(min(ax),max(ax)), axes=FALSE, cex=2, 
              cex.lab=cex.yax, cex.axis=1.5, pch=pch, mgp=c(mgp.ylab,1,0), 
              cex.main=cex.main, lwd=lwd)
          else plot(c(-1,1),predmat[,i],main=terms1[i],xlab="",ylab="",type="b",
              xlim=c(-1.3,1.3),ylim=c(min(ax),max(ax)), axes=FALSE, cex=2, 
              cex.lab=1.2, cex.axis=1.5, pch=pch, cex.main=cex.main, lwd=lwd)
              
          box(which="figure")
          abline(h=ymean,xpd=TRUE)   ## line for mean in all plots
          axis(1, at = c(-1,1), labels = labs[[i]], 
               cex.axis=cex.xax, xpd=NA,lwd=lwd) ## draw bottom axes; 
                                         ## annotation may extend
                                         ## into outer area
          if (i==1)
          axis(2, at = ax, labels = ax, cex.axis=cex.yax, outer=TRUE,lwd=lwd)
              ## draw left-hand-side axis into the outer area
              ## (axis label comes from first actual plot
              ## placement controlled by mgp.ylab
     }
   title(main, line=1.5, outer=TRUE, cex.main=cex.title*cex.main)
   par(mfrow=omfrow,mar=omar,oma=ooma)
   rownames(predmat) <- c("-","+")
   invisible(predmat)
}

