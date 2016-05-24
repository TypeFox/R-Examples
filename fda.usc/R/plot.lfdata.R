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
    plot.fdata(lfdata[[idat]], col =  color2, main =nam[idat],...)
  }
  
}