cubePlot <- function(obj, eff1, eff2, eff3, main=paste("Cube plot for",respnam),
    cex.title=1.5, cex.lab=par("cex.lab"), cex.ax=par("cex.axis"), 
    cex.clab=1.2, size=0.3, round=NULL,
    abbrev=4,y.margin.add=-0.2, modeled=TRUE){
    # obj     a linear model
    # cex     plot character size
    # cex.lab size of variable labels in diagonal panels
    # cex.xax size of annotation for x-axis
    # cex.yax size of annotation for y-axis
    # abbrev  number of characters for factor levels in diagonal panel
    # show.alias  show effects with which this 3way-interaction is aliased 

   respnam <- gsub("$","",c(substitute(obj)),fixed=TRUE)
   sub <- paste("modeled =",modeled)
   
   if (!"lm" %in% class(obj)) {
   sub <- NULL
   hilf <- data.frame(obj, eff1, eff2, eff3)
   colnames(hilf) <- gsub("$","_",c(substitute(obj), substitute(eff1), substitute(eff2), substitute(eff3)),fixed=TRUE)
   respnam <- colnames(hilf[1])
   obj <- eval(parse(text=paste("lm(",respnam,"~.^3,hilf)")))
   #colnames(obj$model)[2:4] <- c(substitute(eff1),substitute(eff2),substitute(eff3))
   eff1<-colnames(obj$model)[2]
   eff2<-colnames(obj$model)[3]
   eff3<-colnames(obj$model)[4]
   }
   else respnam <- colnames(obj$model)[attr(obj$"terms","response")]
   xends <- yends <- zends <- c(-1,1)
   if (is.factor(eval(parse(text=paste("obj$model$",eff1,sep=""))))) 
       xends <- substr(levels(eval(parse(text=paste("obj$model$",eff1,sep="")))),1,abbrev)
   if (is.factor(eval(parse(text=paste("obj$model$",eff2,sep=""))))) 
       yends <- substr(levels(eval(parse(text=paste("obj$model$",eff2,sep="")))),1,abbrev)
   if (is.factor(eval(parse(text=paste("obj$model$",eff3,sep=""))))) 
       zends <- substr(levels(eval(parse(text=paste("obj$model$",eff3,sep="")))),1,abbrev)

   ### base everything on model recoded to -1 and 1 numerics
   obj <- remodel(obj)
   ## extract labels for factor levels to desired length (abbrev)
   labs <- lapply(obj$labs,function(sp) substr(sp,1,abbrev))
   ## reduce obj list to its linear model component
   obj <- obj$model
   mod <- obj$model
   respnam <- colnames(mod)[attr(attr(mod,"terms"),"response")]
   ## for some reason, this does not work when using obj instead of mod
   
   if (!check(obj))
     stop("This routine is applicable for 2-level factorial designs without partial aliasing only.")
   if (modeled) basis <- aggregate(predict(obj),list(eff1=eval(parse(text=paste("obj$model$",eff1,sep=""))),
        eff2=eval(parse(text=paste("obj$model$",eff2,sep=""))),
        eff3=eval(parse(text=paste("obj$model$",eff3,sep=""))) ),mean)
   else basis <- aggregate(eval(parse(text=paste("obj$model$",respnam,sep=""))),
        list(eff1=eval(parse(text=paste("obj$model$",eff1,sep=""))),
        eff2=eval(parse(text=paste("obj$model$",eff2,sep=""))),
        eff3=eval(parse(text=paste("obj$model$",eff3,sep=""))) ),mean)
   if (!is.null(round)) basis[,4] <- round(basis[,4],round)
 bild <- cubedraw(min=-1,max=1,ticks=TRUE,y.margin.add=y.margin.add,
             xlab=eff1,ylab=eff2,zlab=eff3,  
             xends=xends,yends=yends,zends=zends,
             cex.lab=cex.ax, cex.axlab=cex.lab, 
             cex.main=cex.title*par("cex.main"), main=main, sub=sub,
             mar=c(5,4,4,4)+0.1)
 cubecorners(bild,circles=1:8,size=size)
 cubelabel(bild,labs=basis[,4],pos=NULL,adj=c(0.4,0.4),cex=cex.clab)
}