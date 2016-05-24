"plot.pairwiseCI" <-
function(x, CIvert=NULL, CIlty = 1, CIlwd=1, CIcex=1, H0line=NULL, H0lty=1, H0lwd=1, main=NULL, ylab="", xlab="", ... )
{

old.par <- par(no.readonly=TRUE)

# # # Check the arguments:

if(is.null(CIvert) == FALSE && is.logical(CIvert) == FALSE ) 
 {stop("argument CIvert must be specified as single logical value: either TRUE or FALSE ")}

if(length(CIvert) > 1) 
 {stop("argument CIvert must be specified as single logical value: either TRUE or FALSE ")}

if(is.numeric(CIlty) == FALSE || length(CIlty) != 1)
 {stop("CIlty must be specified as single integer (as appropriate for argument lty in lines)")}

if(is.numeric(CIlwd) == FALSE || length(CIlty) != 1)
 {stop("CIlwd must be specified as single integer (as appropriate for argument lwd in lines)")}

if(is.numeric(CIcex) == FALSE || length(CIcex) != 1)
 {stop("CIcex must be specified as single integer (as appropriate for argument cex in par)")}

if(is.numeric(H0lty) == FALSE)
 {stop("H0lty must be specified as integer (as appropriate for argument lty in lines)")}

if(is.numeric(H0lwd) == FALSE)
 {stop("H0lwd must be specified as integer (as appropriate for argument lwd in lines)")}


byout <- x$byout
bynames <- x$bynames
method <- x$method
alternative <- x$alternative
conf.level <- x$conf.level


# # # shall the CI be plotted vertical or horizontal

if(is.null(CIvert)){CIvert <- FALSE}


# # # which H0line shall be plotted: find a default

if(is.null(H0line))
 {
  if( any(c("Param.diff", "Lognorm.diff", "HL.diff", "HD.diff", "Median.diff", "Prop.diff") == method )) { H0line <- 0 } 
   else { H0line <- 1 }
 }
 else
 {
  if(mode(H0line)!="numeric")
   {stop("Specify a numeric value to be plotted as H0line")}

  if( any(c("Param.ratio", "Lognorm.ratio", "HL.ratio", "HD.ratio", "Median.ratio", "Prop.ratio", "Prop.or") == method ) && any(H0line <= 0) ) 
   {stop("the ratio is defined positive, specify H0line as a value greater 0 ")}
 }

# # H0line, H0lty and H0lwd to equal length:

H0par <- cbind(H0line, H0lty, H0lwd)
H0line <- H0par[,1]
H0lty <- H0par[,2]
H0lwd <- H0par[,3]

# # # only one by-group:

if(length(byout)==1)

{

 LOWER <- byout[[1]]$lower
 UPPER <- byout[[1]]$upper
 ESTIMATE <- byout[[1]]$estimate
 COMPNAMES <- byout[[1]]$compnames
 nest <- length(ESTIMATE)
 num <- 1:nest # index vector for the number of comparisons
 

 # # # all vectors of equal length??

 if( any( c(length(LOWER), length(UPPER), length(ESTIMATE)) != length(COMPNAMES) ) )
  { stop("plot.pairwiseCI INTERNAL: length of any (LOWER, UPPER, ESTIMATE, COMPNAMES) is not the same") }

}else{

# # # more than one by-group:

  k <- length(bynames)
  
  upperk <- numeric()   
  lowerk <- numeric()
  estimatek <- numeric()
  compnamesk <- character()
  
   for(e in 1:k)
    {
    upperk <- c(upperk, byout[[e]]$upper)  
    lowerk <- c(lowerk, byout[[e]]$lower)     
    estimatek <- c(estimatek, byout[[e]]$estimate)     
    compnamesk <- c(compnamesk, paste(bynames[[e]], ":", byout[[e]]$compnames) ) 
    }
 
  LOWER <-lowerk
  UPPER <-upperk
  ESTIMATE <-  estimatek
  COMPNAMES <- compnamesk 

  nest <- length(estimatek)
  num <- 1:nest
}

 # # # find the PLOT margins:

 if(alternative=="two.sided")
  {
   lplot <- min(LOWER, H0line )  
   uplot <- max(UPPER, H0line )
  }

 if(alternative=="less")
  {
   lplot <- min(H0line, ESTIMATE) 
   uplot <- max(UPPER, H0line) 
  }

 if(alternative=="greater")
  {
   lplot <- min(H0line, LOWER) 
   uplot <- max(ESTIMATE, H0line) 
  }


 # # # define MAIN, SUB, YLAB,. . .

if(is.null(main)){main <- ""}

# adjustment of cex.axis?

dargs <- list(...)
if(is.null(dargs$cex.axis)){CEX.AXIS <- par("cex.axis")}else{CEX.AXIS <- dargs$cex.axis}


# # # plot with vertical CIs

if(CIvert)
 {

 plot.new() 

 # the default margin size in inches
  mymai <- par(no.readonly=TRUE)$mai

 # adjust margin under the x axis according to length of comparison names
  xwidth<- 1.3 * max(strwidth(COMPNAMES, units = "inches", cex = CEX.AXIS)) 

 if(mymai[1] < xwidth){mymai[1] <- xwidth}

 par(mai=mymai, new=TRUE)

plot(x = num, y = ESTIMATE, axes = FALSE, ylim = c(lplot, uplot), xlim=c(0.7, nest+0.3),
 type="p", pch=16, cex=CIcex,
 xlab=" ",
 ylab=ylab,
 main=main
 )

axis(side = 1, at = num, labels=COMPNAMES, las=2, ...)
axis(side=2, ...)
box()

abline(h=H0line, lty=H0lty, lwd=H0lwd)

CIlengthin <- 0.3*par("fin")[1]/nest

if(alternative=="two.sided"){arrows(x0 = num, x1=num, y0 = LOWER, y1=UPPER,length=CIlengthin, angle=90, code=3 ) }
if(alternative=="less"){arrows(x0 = num, x1=num, y0 = lplot, y1=UPPER,length=CIlengthin, angle=90, code=2 )}
if(alternative=="greater"){arrows(x0 = num, x1=num, y0 = LOWER, y1=uplot, length=CIlengthin, angle=90, code=1 )}
}


# # # plot with horizontal CIs

if(!CIvert)
 {

 plot.new()

# the default margin size in inches
   mymai <- par(no.readonly=TRUE)$mai

 # adjust margin left of the y axis according to length of comparison names
  ywidth<- 1.3 * max(strwidth(COMPNAMES, units = "inches", cex = CEX.AXIS)) 

 if(mymai[2] < ywidth){mymai[2] <- ywidth}
 
 par(mai=mymai, new=TRUE)

num<-rev(num)

plot(x = ESTIMATE, y = num, axes = FALSE, xlim = c(lplot, uplot), ylim=c(0.7, nest+0.3),
 type="p", pch=16, cex=CIcex,
 ylab=" ",
 xlab=xlab,
  main=main
 )

axis(side = 2, at = num, labels=COMPNAMES, las=2, ...)
axis(side=1, ...)
box()

CIlengthin <- 0.3*par("fin")[2]/nest
abline(v=H0line, lty=H0lty, lwd=H0lwd)
if(alternative=="two.sided"){arrows(y0 = num, y1=num, x0 = LOWER, x1=UPPER,length=CIlengthin, angle=90, code=3 ) }
if(alternative=="less"){arrows(y0 = num, y1=num, x0 = lplot, x1=UPPER,length=CIlengthin, angle=90, code=2 )}
if(alternative=="greater"){arrows(y0 = num, y1=num, x0 = LOWER, x1=uplot, length=CIlengthin, angle=90, code=1 )}

}

par(old.par)

}

