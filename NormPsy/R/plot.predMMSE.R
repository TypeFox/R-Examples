plot.predMMSE <- function(x,legend.loc="topright",legend,add=FALSE,...)
{
  if (!inherits(x, "predMMSE")) stop("use only with \"predMMSE\" objects")
  if(is.na(as.logical(add))) stop("add should be TRUE or FALSE") 
  
   ng <- 1
   ndistr <- grep("MMSEdistr",colnames(x))
   if(length(ndistr))
   {
    ng <- length(ndistr)/3 
   }
   else
   {
    if(ncol(x)>2)
    {
     ng <- ncol(x)-1 
    }
   }
  
  if(missing(legend)) legend <- paste("class",1:ng,sep="") 
  
   dots <- list(...)

   if(length(list(...)$main)) 
   {
    main1 <- as.character(eval(match.call()$main))
    dots <- dots[setdiff(names(dots),"main")]
   }
    else main1 <- "Prediction of MMSE scores"
                   

   if(length(list(...)$type))    
   {
    type1 <- eval(match.call()$type)
    dots <- dots[-which(names(dots)=="type")]
   }
   else  type1 <- "l"

  
   if(length(list(...)$col))    
   {
    col1 <- eval(match.call()$col)
    col1 <- rep(col1,length.out=ng) 
    dots <- dots[-which(names(dots)=="col")]
   }
   else  col1 <- rainbow(ng)  
  

   if(length(list(...)$ylim)) 
   {
    ylim1 <- eval(match.call()$ylim)
    dots <- dots[setdiff(names(dots),"ylim")]
   }
   else ylim1 <- c(0,30)
   
   if(length(list(...)$xlab)) 
   {
    xlab1 <- as.character(eval(match.call()$xlab))
    dots <- dots[setdiff(names(dots),"xlab")]
   }
   else xlab1 <- "prediction time"

   if(length(list(...)$ylab)) 
   {
    ylab1 <- as.character(eval(match.call()$ylab))
    dots <- dots[setdiff(names(dots),"ylab")]
   }
   else ylab1 <- "MMSE"
 
   
   
   names.plot <- c("adj","ann","asp","axes","bg","bty","cex","cex.axis","cex.lab","cex.main","cex.sub","col","col.axis",
   "col.lab","col.main","col.sub","crt","err","family","fig","fin","font","font.axis","font.lab","font.main","font.sub",
   "frame.plot","lab","las","lend","lheight","ljoin","lmitre","lty","lwd","mai","main","mar","mex","mgp","mkh","oma",
   "omd","omi","pch","pin","plt","ps","pty","smo","srt","sub","tck","tcl","type","usr","xaxp","xaxs","xaxt","xlab",
   "xlim","xpd","yaxp","yaxs","yaxt","ylab","ylbias","ylim") 
   dots.plot <- dots[intersect(names(dots),names.plot)]   
  
   if(!isTRUE(add))
   {
    do.call("matplot",c(dots.plot,list(x=x[,rep(1,ng),drop=FALSE],y=x[,1+1:ng,drop=FALSE],type=type1,ylab=ylab1,xlab=xlab1,main=main1,ylim=ylim1,col=col1)))
   }
   else
   {
    do.call("matlines",c(dots.plot,list(x=x[,rep(1,ng),drop=FALSE],y=x[,1+1:ng,drop=FALSE],type=type1,col=col1)))
   }   
  
   if(length(ndistr))
   {
    if(length(list(...)$lty))  dots.plot <- dots.plot[setdiff(names(dots),"lty")] 
     
    do.call("matlines",c(dots.plot,list(x=x[,rep(1,ng),drop=FALSE],y=x[,1+ng+1:ng,drop=FALSE],col=col1,lty=2)))
    do.call("matlines",c(dots.plot,list(x=x[,rep(1,ng),drop=FALSE],y=x[,1+2*ng+1:ng,drop=FALSE],col=col1,lty=2))) 
   }
  
   if(!is.null(legend))
   {
    if(length(list(...)$box.lty))
    {
     box.lty1 <- as.integer(eval(match.call()$box.lty))
     dots <- dots[setdiff(names(dots),"box.lty")]
    }
    else box.lty1 <- 0

    if(length(list(...)$inset))
    {
     inset1 <- eval(match.call()$inset)
     dots <- dots[setdiff(names(dots),"inset")]
    }
    else inset1 <- c(0.02,0.02) 
     
    if(length(list(...)$lty))
    {
     lty1 <- eval(match.call()$lty)
     dots <- dots[setdiff(names(dots),"lty")]
    }
    else  lty1 <- 1
     
    names.legend <- c("fill","border","lty","lwd","pch","angle","density","bg","box.lwd",   
    "box.lty","box.col","pt.bg","cex","pt.cex","pt.lwd","xjust","yjust","x.intersp","y.intersp","adj","text.width",
    "text.col","text.font","merge","trace","plot","ncol","horiz","title","xpd","title.col","title.adj","seg.len") 

    dots.leg <- dots[intersect(names(dots),names.legend)]
    if(type1=="l" | type1=="b") dots.leg <- c(dots.leg,list(lty=lty1))
    if(!(type1 %in% c("l","b"))) dots.leg <- dots.leg[setdiff(names(dots),"lwd")]
     
    do.call("legend",c(dots.leg,list(x=legend.loc, legend=legend, box.lty=box.lty1, inset=inset1,col=col1)))     
   }
 
  return(invisible(NULL))
}












