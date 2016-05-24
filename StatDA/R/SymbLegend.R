SymbLegend <-
function(X,Y,z,type="percentile",qutiles=c(0,0.05,0.25,0.75,0.95,1),q=NULL,symbtype="EDA",symbmagn=0.8,
    leg.position="topright",leg.title="", leg.title.cex=0.8, leg.round=2,leg.wid=4,leg.just="right",
    cex.scale=0.8,xf=9e3,logscale=TRUE,accentuate=FALSE)
{
# Plots symbols and Legend on a map
# either percentile symbols, or boxplot symbols
#
# X ... X-coordinates
# Y ... Y-coordinates
# z ... values on the coordinates
# type ... "percentile" for percentile legend
#      ... "boxplot" for boxplot legend
# qutiles ... considered quantiles if type="percentile" is used
# q ... if not NULL, provide manually data points where to break
# symbtype ... type of symbols to be used; "EDA","EDAacc","EDAacc2",
#              "EDAext","GSC","arbit"
# symbmagn ... magnification factor for symbols
# leg.position ... position of the legend, either character like "topright" or coordinates
# leg.title ... title for legend
# leg.title.cex ... cex for legend title
# leg.round ... round legend to specified digits "pretty"
# leg.wid ... width (space in numbers) for legend
# leg.just ... how to justify the legend
# cex.scale ... cex for text "log-scale" and for boxplot legend - only for type="boxplot"
# xf ... x-distance from boxplot to number for legend
# logscale ... if TRUE a log scale is used (for boxplot scale) and the log-boxplot is computed
# accentuate ... if TRUE, accentuated symb ls are used (here only EDA accentuated!)
# 

if (symbtype=="EDA"){ symb=EDA.symbol; ssize=EDA.size}
else if (symbtype=="EDAacc"){ symb=EDAacc.symbol; ssize=EDAacc.size}
else if (symbtype=="EDAacc2"){ symb=EDAacc2.symbol; ssize=EDAacc2.size}
else if (symbtype=="EDAext"){ symb=EDAext.symbol; ssize=EDAext.size}
else if (symbtype=="GSC"){ symb=GSC.symbol; ssize=GSC.size}
else if (symbtype=="arbit"){ symb=arbit.symbol; ssize=arbit.size}
else stop("Symbol type not known - select another or use default!")

if (accentuate) {symb=EDAacc3.symbol; ssize=EDAacc3.size} # for boxplot symbols

ssize=ssize*symbmagn

if (type=="percentile"){
  if (is.null(q)) q=quantile(z, qutiles) # compute quantiles of the variable
  # plot symbols:
  for (i in 1:(length(qutiles)-1)){
    points(X[z>=q[i] & z<q[i+1]], Y[z>=q[i] & z<q[i+1]], pch=symb[i], cex=ssize[i])
  }
#  points(X[z>=q[length(qutiles)]],Y[z>=q[length(qutiles)]],pch=symb[length(qutiles)],
#              cex=ssize[length(qutiles)])
  points(X[z>=q[length(qutiles)]],Y[z>=q[length(qutiles)]],pch=symb[length(symb)],
              cex=ssize[length(ssize)])
  #Legend
  leg=rep(NA,length(q)-1)
  leg[1]=paste("  ",roundpretty(q[1],leg.round),"-",format(roundpretty(q[1+1],leg.round),
         width=leg.wid,justify=leg.just))
  for (i in 2:(length(q)-1)){
    leg[i]=paste(">",roundpretty(q[i],leg.round),"-",format(roundpretty(q[i+1],leg.round),
         width=leg.wid,justify=leg.just))
  }
  legend(leg.position,pch=rev(symb),pt.cex=rev(ssize), legend=rev(leg), title=leg.title, cex=leg.title.cex)
  }
else { # boxplot legend:
  if (logscale) {
    boxinfo=boxplotlog(z,plot=FALSE)
  }
  else {
    boxinfo=boxplot(z,plot=FALSE)
  }
  ret <- boxplotlegend(X,Y,z,boxinfo,legend.title=leg.title,logscale=logscale,
       symb=symb,ssize=ssize,accentuate=accentuate,cex.legtit=leg.title.cex,cex.scale=cex.scale,xf=xf)
  symbnew=ret$symb
  ssizenew=ret$ssize
  if (is.null(q)) q=ret$q
  # plot symbols:
  for (i in 1:(length(symbnew))){
    points(X[z>=q[i] & z<q[i+1]], Y[z>=q[i] & z<q[i+1]], pch=symbnew[i], cex=ssizenew[i])
  }
}

invisible()
}

