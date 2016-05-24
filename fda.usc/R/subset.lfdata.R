################################################################################
################################################################################
plot.lfdata<-function(lfdata,ask=FALSE,color,...){
  mf=5
  nvar<-length(lfdata)
  if (nvar>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
  else{    mf<-switch(nvar,
                      "1"={c(1,1)},
                      "2"={c(1,2)},
                      "3"={c(1,3)},
                      "4"={c(2,2)})            
           par(mfrow =mf)                    }
  names1<-names2<-names<-lfdata[[1]][["names"]]
  names1$main<-"Multivariate Functional Data"
  #    tr<-paste("mode.tr",trim*100,"\u0025",sep="")         
  nam<-names(lfdata)
  
  if (is.null(nam)) nam<-1:nvar
  for (idat in 1:nvar) {
    data<-lfdata[[idat]]$data
    tt<-lfdata[[idat]]$argvals
    rtt<-lfdata[[idat]]$rangeval
    if (missing(color)) color2<-1
    else {
      if (is.list(color)) color2<-color[[idat]]
      else color2<-color
    }
    plot(lfdata[[idat]], col =  color2,lty=1, main =nam[idat],...)
  }
  
}
################################################################################
################################################################################
#plot.lfdata(ab,col=color,lwd=2)
#el color puede ser un alista para cada fdata o un vector el mismo para todos
#plot.lfdata(ab)
################################################################################
subset.lfdata<-function(x,subset,select,drop=TRUE,...){
  #if (any(class(x)!="lfdata")) stop("No list class object")
  nvar<-length(x)
  newx<-x
  for (i in 1:nvar){
    newx[[i]]<-subset.fdata(x[[i]],subset,select,drop,...)   
  } 
  return(invisible(newx))
}
################################################################################
################################################################################
subset.fdata<-function(x,subset,select,drop=TRUE,...){
  if (!is.fdata(x)) stop("No fdata class object")
  if (missing(select) & missing(subset)) stop("You must specify at least the argument 'subset' or 'select'")
  if (missing(subset)) subset<-!logical(nrow(x$data))
  range<-TRUE
  if (missing(select)) {
    range<-FALSE
    select<-1:ncol(x$data)
  }
  if (is.numeric(subset)) {
    subset2<-logical(nrow(x$data))
    subset2[subset]<-TRUE
    subset<-subset2    
  }
  newx<-x
  newx$data<-subset(x$data,subset,select,drop=drop,...)
  #print(select)
  newx$argvals<-newx$argvals[select]
  if (range) newx$rangeval<-c(min(newx$argvals),max(newx$argvals))
  return(invisible(newx))
}
################################################################################
################################################################################