  illustrateLLN <- function(Distr = Norm(), 
                            n = c(1,3,5,10,25,50,100,500,1000,10000), 
                            m = 50, step = 1, sleep = 0, 
                            withConf = TRUE, withCover = (length(n)<=12), 
                            withEline = TRUE,  withLegend = TRUE, 
                            CLTorCheb = "CLT", coverage = 0.95, ..., 
                            col.Eline = "blue", lwd.Eline = par("lwd"), 
                            lty.Eline = par("lty"), 
                            col.Conf = "red", lwd.Conf = par("lwd"), 
                            lty.Conf = 2, cex.Cover = 0.7, cex.legend = 0.8){
  
  Distrc <- match.call(call = sys.call(sys.parent(1)))$Distr
  if(is.null(Distrc)) Distrc <- "Norm()"
  else Distrc <- as.character(deparse(Distrc))

  if(!is.numeric(coverage) || length(coverage)>1 || any(coverage <= 0)
     || any(coverage >= 1) )
     stop("Argument 'coverage' must be a single number in (0,1)")
     
  dots <- match.call(call = sys.call(sys.parent(1)), 
                      expand.dots = FALSE)$"..."

  dots.for.lines <- dots[! names(dots) %in% c("col", "lwd", "lty")]
  dots.for.matplot <-   dots.for.lines[! names(dots.for.lines) %in% 
                           c("ylim", "xlim", "pch", "xlab", "ylab", "axes", 
                             "main")]
  dots.for.legend <-   dots.for.lines[! names(dots.for.lines) %in% 
                           c("legend", "main", "cex", "sub")]
  dots.for.text <-   dots[! names(dots) %in% c("cex", "main", "sub")]

  confType <- pmatch(CLTorCheb, c("CLT","Chebyshev"), nomatch = 1)
  
     col <- if (!hasArg(col)) par("col") else dots$col
     if (hasArg(col) && missing(col.Eline))
         col.Eline <- col
     if (hasArg(col) && missing(col.Conf))
         col.Conf <- col

     lwd <- if (!hasArg(lwd)) par("lwd") else dots$lwd
     if (hasArg(lwd) && missing(lwd.Eline))
         lwd.Eline <- lwd
     if (hasArg(lwd) && missing(lwd.Conf))
         lwd.Conf <- lwd

     lty <- if (!hasArg(lty)) 1 else dots$lty
     if (hasArg(lty) && missing(lty.Eline))
         lty.Eline <- lty
     if (hasArg(lty) && missing(lty.Conf))
         lty.Conf <- lty

     cex <- if (!hasArg(cex)) 1 else dots$cex
     if (hasArg(cex) && missing(cex.Cover))
         cex.Cover <- cex
     if (hasArg(cex) && missing(cex.legend))
         cex.legend <- cex

     if (hasArg(lty) && missing(lty.Eline))
         lty.Eline <- lty
     if (hasArg(lty) && missing(lty.Conf))
         lty.Conf <- lty

     pch <- if (!hasArg(pch)) 16 else dots$pch

   facConf <- switch(confType, qnorm((1+coverage)/2), (1-coverage)^(-.5))
   facName <- switch(confType, "CLT", "Chebyshev")
   #confidence bounds based on CLT
   legend.txt <- if (hasArg(legend)) dots$legend else                        
         gettextf("%s-based %3.0f%% (pointwise) confidence interval", 
                         facName, round(100*coverage,0))

  da <- matrix(NA,m,length(n))
  
  omar <- par(no.readonly = TRUE)
#  omar$cin <- omar$cra <- omar$csi <- omar$cxy <-  omar$din <- NULL
  on.exit(par(omar))
     ## getting the parameter

  slots <-  slotNames(param(Distr))
  slots <-  slots[slots != "name"]
  nrvalues <-  length(slots)
  if(nrvalues > 0){
        values <-  numeric(nrvalues)
    for(i in 1:nrvalues)
      values[i] <-  attributes(attributes(Distr)$param)[[slots[i]]]
    paramstring <-  paste(values, collapse = ", ")
    nparamstring <-  paste(slots, "=", values, collapse = ", ")
    qparamstring <- paste("(",paramstring,")",sep="")
  }
  else paramstring <- qparamstring <- nparamstring <- ""


  .mpresubs <- function(inx)
                 .presubs(inx, c("%C", "%D", "%N", "%P", "%Q", "%A",
                                         "%X"),
                       list(as.character(class(Distr)[1]),
                         as.character(date()),
                         nparamstring,
                         paramstring,
                         qparamstring,
                         Distrc,
                         expression(~~~~~bar(X)[n]==~~sum(X[i],i==1,n)/n)))

  xlab <- if (!hasArg(xlab)) gettext("Sample size n") else dots$xlab
  xlab <- .mpresubs(xlab)
  ylab <- if (!hasArg(ylab)) "Realisations of %X" else dots$ylab
  ylab <- .mpresubs(ylab)

  tit <- c("LLN: Convergence against E[X] -- case %C%Q",
           "called with %A")
  if ( is.na (E(Distr)))
       tit <- c("LLN: Non-Convergence against E[X] -- case %C%Q",
                "called with %A")
    
  main <- if ( hasArg(main)) dots$main else tit
  main <- .mpresubs(main)

  sub <- if ( hasArg(sub)) dots$sub else ""
  sub <- .mpresubs(sub)

  LLNin <- function(x, n, from, to = from){
  for(i in seq(length(n)))
      da[from:to,i] <<- 
         rowMeans(matrix(r(x)(n[i]*(to-from+1)),(to-from+1),n[i]))
  }

   ## location and scale:
   
   mE <-  if(!is.na(E(Distr)))   E(Distr) else median(Distr)
   msd <- if(!is.na(sd(Distr))) sd(Distr) else 2 * mad(Distr)
 
   Ns <- seq(length(n))
   nn <- if(is(Distr,"Cauchy")) n*0+1 else n


   
   for(j in seq(1, m, by = step)) 
     {LLNin(Distr, n, j, j+step-1)
      do.call(matplot, args = c(list(Ns, t(da), pch = pch, col=col,  
              axes = FALSE, ylim = q(Distr)(c(0.02,0.98)), 
              xlab = xlab, ylab = "",  main = main), dots.for.matplot ))
 
      title(ylab = ylab, line = 1.7)
      axis(1, at = Ns, labels = n)
      axis(2)
      if (withEline)
          do.call( abline, args= c(list(h = mE, col = col.Eline, lwd = lwd.Eline, 
                   lty = lty.Eline), dots.for.lines))
      if (withConf)
          do.call( matlines, args = c(list(
                   Ns, mE + facConf * msd/sqrt(nn)%o%c(-1,1), 
                   lty = lty.Conf, col = col.Conf, lwd = lwd.Conf
                                          ), dots.for.lines))
               
      coverage <- colMeans( t(t(da) <= (mE+facConf*msd*1/sqrt(nn)) & 
                          t(da) >= (mE-facConf*msd*1/sqrt(nn))) , 
                          na.rm= TRUE )
      if (withCover && withConf)
          do.call(mtext, args = c(list(at = c(0,Ns), 
                c(gettext("coverage"),sprintf("%3.2f",round(coverage,2))), 
                cex = cex.Cover), dots.for.text))
      if (withLegend && withConf)
          do.call(legend, args = c(list("bottomright", 
                legend = legend.txt, cex = cex.legend,
                col = col.Conf, lwd = lwd.Conf, lty = lty.Conf),  
                dots.for.legend))
                   
      Sys.sleep(sleep)
      }
  }
  

#------------------------------------
#### utility copied from package distr v.2.6  svn-rev 943
#------------------------------------
.presubs <- function(inp, frompat, topat){
### replaces in an expression or a string all frompat patterns to topat patterns

logic <- FALSE
inCx <- sapply(inp,
   function(inpx){
      inC <- deparse(inpx)
      l <- length(frompat)
      for(i in 1:l)
         { if (is.language(topat[[i]])){
               totxt <- deparse(topat[[i]])
               totxt <- gsub("expression\\(", "\", ", gsub("\\)$",", \"",totxt))
               if (length(grep(frompat[i],inC))) logic <<- TRUE
               inC <- gsub(frompat[i],totxt,inC)
           }else inC <- gsub(frompat[i], topat[[i]], inC)
         }
      return(inC)
    })
if(length(grep("expression",inCx))>0)
   inCx <- gsub("expression\\(", "", gsub("\\)$","",inCx))
if (length(inCx) > 1) {
   inCx <- paste(inCx, c(rep(",", length(inCx)-1), ""),
                 sep = "", collapse = "\"\\n\",")
   if ( any(as.logical(c(lapply(inp,is.language)))) | logic )
      inCx <- paste("expression(paste(", gsub("\\\\n"," ", inCx), "))", sep ="")
   else
      inCx <- paste("paste(",inCx,")", sep ="")
}else inCx <- paste("expression(paste(",inCx,"))",sep="")
outC <- eval(parse(text = eval(inCx)))
return(outC)
}
