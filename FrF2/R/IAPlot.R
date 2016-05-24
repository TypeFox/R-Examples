IAPlot <- function(obj, ...){
    UseMethod("IAPlot")
}
IAPlot.design <- function(obj, ..., response=NULL){
    if (!"design" %in% class(obj)) 
        stop("IAPlot.design works for obj from class design only.")
    di <- design.info(obj)
    if (is.null(di$response.names)) 
        stop("The design obj must have at least one response.")
    if (!(is.null(response))) 
      if (!response %in% di$response.names)
        stop("Requested response is not a response variable in obj.")
    if (!(length(grep("FrF2",di$type))>0 | 
           length(grep("pb",di$type))>0)) { 
           if (!(di$type=="full factorial" & all(di$nlevels==2)))
        stop("The design obj must be of a type containing FrF2 or pb.")
       }
    IAPlot(lm(obj, degree=2, response=response), ...)
}

IAPlot.default <-
function(obj, main=paste("Interaction plot matrix for",respnam), pch=c(15,17), 
    cex.lab=par("cex.lab"), cex=par("cex"), cex.xax=par("cex.axis"), 
    cex.yax=cex.xax, cex.title = 1.5, lwd=par("lwd"), abbrev=4, select=NULL, 
    show.alias=FALSE, ...){
    # obj     a linear model
    # pch     plot characters used
    # cex     plot character size
    # cex.lab size of variable labels in diagonal panels
    # cex.xax size of annotation for x-axis
    # cex.yax size of annotation for y-axis
    # abbrev  number of characters for factor levels in diagonal panel
    # show.alias  show number of effect in each panel, 
    #         in order to allow immediate judgment which effects are aliased

   if (! ("lm" %in% class(obj) | "aov" %in% class(obj))) 
      stop("obj must be a linear model object (lm or aov), or a design of class design")
    
   ### base everything on model recoded to -1 and 1 numerics
   obj <- remodel(obj)
   ## extract labels for factor levels to desired length (abbrev)
   labs <- lapply(obj$labs,function(sp) substr(sp,1,abbrev))
   ## reduce obj list to its linear model component
   obj <- obj$model
   mod <- obj$model
   respnam <- colnames(mod)[attr(attr(mod,"terms"),"response")]
   ## for some reason, this does not work when using obj instead of mod
   ## it neither works if respnam is something like rnorm(16)
      ## how do I fix or at least report this ?
  
   if (!check(obj))
     stop("This routine is applicable for 2-level factorial designs without partial aliasing only.")
     
   ## prepare for simple access of important quantities
   term.ord <- attr(terms(obj),"order")
   mains <- which(term.ord==1)
   if (!any(term.ord==2)) stop("no 2-factor interaction terms in model")
   ints <- which(term.ord==2)
   nmain <- length(mains)
   nint <- length(ints)
   intcol <- attr(terms(obj),"intercept")
   termnames <- attr(terms(obj),"term.labels")
   
   if (!is.null(select)){
       ## take care of selecting only part of the interactions
       if (!is.numeric(select)) stop("select must be numeric")
       if (!all(floor(select)==select)) stop("select must contain integer numbers")
       if (any(select<1 | select>nmain)) stop("select must contain numbers betweeen 1 and ", nmain, " only")
       select <- unique(select)
       if (length(select)<2) stop("at least 2 effects must be selected for an interaction plot matrix")
       selnam <- termnames[select]
       labs <- labs[select]
       hilf <- strsplit(termnames[ints],":")
       selected <- sapply(hilf, function(.obj1) all(.obj1 %in% selnam))
       intselnam <- termnames[ints[selected]]
       if (!is.name(respnam)) formulasel <- paste(respnam, "~", paste(c(selnam,intselnam),collapse="+"))
       formulasel <- paste(respnam, "~", paste(c(selnam,intselnam),collapse="+"))
       obj <- lm(formulasel, data=mod)
       term.ord <- attr(terms(obj),"order")
       mains <- which(term.ord==1)
       ints <- which(term.ord==2)
       nmain <- length(mains)
       nint <- length(ints)
       intcol <- attr(terms(obj),"intercept")
   }
   mm <- model.matrix(obj)
   ## omit intercept, if present
   if (intcol > 0) mm <- mm[,-intcol]
   
   coefs <- coef(obj)
   if (intcol > 0) coefs <- coefs[-intcol] 

 
   terms1 <- colnames(mm)[mains]
   terms2 <- colnames(mm)[ints]
   ## used for easy prediction
   addnam <- setdiff(colnames(obj$model), terms1)
   names <- c(terms1, addnam)
   
   ## create matrix for plots
   predmat <- matrix(rep(0,4*nint),4,nint)
   varnums <- matrix(rep(0,2*nint),2,nint)
   colnames(predmat) <- terms2
   
   ## extract aliasing information
   ## for optional annotation of plots
   al <- aliases(obj)$aliases
   if (is.null(al)) exist <- names(coefs) 
   else exist <- sapply(al, function(sp) sp[1])
   alnum <- rep(0,length(coefs))
   names(alnum) <- names(coefs)
   
   ## calculate alias numbers for aliasing information
   ## and replace NA coefficients for aliased effects 
   ## with the appropriately signed replacement from "master"
   for (i in 1:length(coefs)){
      if (!is.na(coefs[i])) alnum[i] <- which(exist==names(coefs)[i])
      else{ 
          hilf <- grep(names(coefs)[i],al)
          coefs[i] <- coefs[exist[hilf]]
          alnum[i] <- hilf
          if (substr(al[[hilf]][grep(names(coefs[i]),al[[hilf]])],1,1)=="-")
             coefs[i] <- -coefs[i]
      }
   }
   ## re-include intercept coefficient
   if (intcol==1) coefs <- c(obj$coefficients[intcol], coefs)
   else if (intcol < length(coefs)) 
       coefs <- c(coefs[1:(intcol-1)],obj$coefficients[intcol])
   else coefs <- c(coefs, obj$coefficients[intcol])
   
   ## calculate model-based predicted values for all interactions
   for (i in 1:nint){
      dat <- matrix(rep(0,4*nmain),4,nmain)
      inam <- strsplit(terms2[i],":")[[1]]
      ml <- inam[1]
      mr <- inam[2]
      ## first in terms1
      m1 <- min(grep(ml,terms1,fixed=TRUE),grep(mr,terms1,fixed=TRUE))
      ## second in terms1
      m2 <- max(grep(ml,terms1,fixed=TRUE),grep(mr,terms1,fixed=TRUE))
      ## file in varnums matrix for later use
      varnums[,i] <- c(m1,m2)  
      dat[,m1] <- c(-1,1,-1,1)
      dat[,m2] <- c(-1,-1,1,1)
      colnames(dat) <- terms1
      ## allow for things like variables excluded by "-" etc.
      dat <- cbind(dat, matrix(rep(0,length(addnam)*4),4,length(addnam)))
      colnames(dat) <- names
      dat <- as.data.frame(dat)
      modmat <- model.matrix(lm(terms(obj),dat))
      predmat[,i] <- modmat %*% coefs
   }

   ## save old parameters for restoring them later
   omfrow <- par("mfrow")
   omar <- par("mar")
   ooma <- par("oma")
   ## determine axis limits and graphics parameters
   ax <- pretty(c(min(predmat),max(predmat)))
   par(mfcol=c(nmain,nmain), mar=c(1, 1, 1, 1) + 0.1, oma=c(3,5,4,0.1))
   # i is column index, j is row index
   for (i in 1:nmain){
     for (j in 1:nmain){
         ## interaction terms (regardless whether in model or not)
         ww <- c(paste(terms1[i],terms1[j],sep=":"), 
               paste(terms1[j],terms1[i],sep=":"))
         ## column of alnum corresponding to this ww 
         ## (if any, otherwise integer(0))
         hilf <- which(names(alnum) %in% ww)
         
         ## handle diagonal panels
          if (i==j){ plot(c(-1),min(ax)+5/6*(max(ax)-min(ax)),
                          ylim=c(min(ax),max(ax)),xlim=c(-1.1,1.1),
                          axes=FALSE,xlab="",ylab=respnam,
                          col="red", pch=pch[1], cex=cex)
                box(which="figure")
                points(c(-1),min(ax)+1/6*(max(ax)-min(ax)),pch=pch[2], cex=cex, lwd=lwd)
                text(c(-1,-1),min(ax)+c(5,1)*(max(ax)-min(ax))/6,
                      labs[[i]],pos=4, xpd=NA, cex=cex.lab)
                text(0,(max(ax)+min(ax))/2, terms1[i],adj=0.5,xpd=TRUE, 
                    cex=cex.lab, col="blue")
                ## axes for corner panels
                if (i==1) axis(2, at = ax, labels = ax, cex.axis=cex.yax, 
                      outer=TRUE, lwd=lwd)
                if (j==nmain) axis(1, at = c(-1,1), labels = labs[[j]], 
                    cex.axis=cex.xax, outer=TRUE, lwd=lwd)
                }
          if (i<j) {
              sp <- intfind(i,j,varnums)
              if (is.null(sp)) {
                    plot(c(-10,10),c(-1,1), axes=FALSE, xlab="",
                       xlim=c(-1.1,1.1),ylim=c(min(ax),max(ax)))
                       ## skip panel
                    box(which="figure")
                    }
              else {plot(c(-1,1),predmat[c(1,2),sp],ylab=respnam,xlab=terms1[j],
                     type="b", xlim=c(-1.1,1.1),ylim=c(min(ax),max(ax)), 
                     axes=FALSE, col="red",
                     lty=3, pch=pch[1], cex=cex, lwd=lwd)
                    box(which="figure")
                    lines(c(-1,1),predmat[c(3,4),sp],lwd=lwd)
                    points(c(-1,1),predmat[c(3,4),sp],pch=pch[2], cex=cex, lwd=lwd)
                    }       
              ## axes for border panels
                if (i==1) axis(2, at = ax, labels = ax, cex.axis=cex.yax, outer=TRUE, lwd=lwd)
                if (j==nmain) axis(1, at = c(-1,1), labels = labs[[i]], 
                    cex.axis=cex.xax, outer=TRUE, lwd=lwd)
          }
          if (i>j) {
               sp <- intfind(j,i,varnums)
              if (is.null(sp)) {
                    plot(c(-10,10),c(-1,1), axes=FALSE, xlab="",
                       xlim=c(-1.1,1.1),ylim=c(min(ax),max(ax)))## skip panel
                    box(which="figure")
                    }
              else {plot(c(-1,1),predmat[c(1,3),sp],ylab=respnam,xlab=terms1[j],
                       type="b", xlim=c(-1.1,1.1),ylim=c(min(ax),max(ax)), axes=FALSE, 
                       col="red", lty=3, pch=pch[1], cex=cex, lwd=lwd)
                    box(which="figure")
                    lines(c(-1,1),predmat[c(2,4),sp],lwd=lwd)
                    points(c(-1,1),predmat[c(2,4),sp],pch=pch[2], cex=cex, lwd=lwd)}
         }
         if (show.alias) numact <- if (length(hilf)>0)
             text(0,max(ax),alnum[hilf],xpd=TRUE,cex=cex.lab)
   }
   }
   title(main, line=1.5, outer=TRUE, cex.main=cex.title*par("cex.main"))
   par(mfrow=omfrow,mar=omar,oma=ooma)
   rownames(predmat) <- c("-:-","+:-","-:+","+:+")
   if (show.alias) attr(predmat, "aliasgroups") <- 
        lapply(unique(alnum), function(obj) names(alnum)[alnum==obj])
   invisible(predmat)
}

