slider.smooth.plot.ts<-function(x,...)
{ 
 x.name<-deparse(substitute(x))
 if(missing(x)||length(x)<2) return("Error: x must be a vector")
 y<-x; if(is.ts(x)) { x<-time(x) } else { x<-seq(x) }

 t.x<-x; x<-y
 args<-list(...)
 kind<-c("3RS3R", "3RSS", "3RSR", "3R", "3", "S")
 slider(obj.name="nts",obj.value=x) 
 slider(obj.name="kind",obj.value="3")
 slider(obj.name="history.kind",obj.value="")
 refresh<-function(...){
 # initialization
   choice<-slider(obj.name="kind"); print(choice)
   history.kind<-slider(obj.name="history.kind")
   history.kind<-c(history.kind,choice)
   slider(obj.name="history.kind",obj.value=history.kind)
   xx<-slider(obj.name="nts") 
   xx<-smooth(xx,kind=choice)   
   # plot(t.x,x,type="l",lty=3,bty="n",xlab="",ylab="")
   do.call("plot",c(alist(t.x,x,type="l",lty=3,bty="n"),args))
   title(paste(history.kind,collapse="->"))
   points(t.x,xx,type="l",col="red") # lty=1+which(choice==kind))
   slider(obj.name="nts",obj.value=xx) 
 }
 reset<-function(...){
   slider(obj.name="nts",obj.value=x) 
   slider(obj.name="kind",obj.value="3")
   slider(obj.name="history.kind",obj.value="")
   refresh()
 }
 f1<-function(...){slider(obj.name="kind",obj.value="3RS3R");refresh()}
 f2<-function(...){slider(obj.name="kind",obj.value="3RSS"); refresh()}
 f3<-function(...){slider(obj.name="kind",obj.value="3RSR"); refresh()}
 f4<-function(...){slider(obj.name="kind",obj.value="3R");   refresh()}
 f5<-function(...){slider(obj.name="kind",obj.value="3");    refresh()}
 f6<-function(...){slider(obj.name="kind",obj.value="S");    refresh()}
 refresh()
 slider(but.functions=c(f1,f2,f3,f4,f5,f6),
        but.names=kind,reset.function=reset)
 cat("select type of smoothing a la Tukey and look the result!\n")
}

