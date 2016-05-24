wdConvert<-function(input,from="in",to="pt"){
    convert<-paste(from,to,sep="2")
    if (from==to) out<-input else {
    out<-switch(convert,
         in2pt=input*72,
         pt2in=input/72,
         in2cm=input*2.54,
         cm2in=input/2.54,
         cm2pt=input*28.34646,
         pt2cm=input/28.34646,
         stop("from or to dimension not recognized")
         )}
    return(out)
}
