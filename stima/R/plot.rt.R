plot.rt <-
function(x,digits=2,...) {
     rtobj<-x
     mytext<-function(yval, dev, wt, ylevel, digits, n, use.n)
   {
      if (use.n) {
         paste(format(yval, digits), "\nn=", n, sep = "")
      }
      else {
         paste(format(yval, digits))
      }
   }
   # Determine number of categories toplevel.
   allrows<-as.numeric(rownames(rtobj$trunk))
   morethan1cat<-max(allrows)>10000
   ncat<-max(1,floor(max(allrows)/10000))
   # Reorder rows for rpart.
   for (i in 1:ncat)
   {
      if (morethan1cat)
      {  # Select rows needed.
         selrows<-allrows>i*10000&allrows<(i+1)*10000
         subtrunk<-rtobj$trunk[selrows,]
         rows<-allrows[selrows]-i*10000
         rownames(subtrunk)<-rows
      }
      else
      {
         subtrunk<-rtobj$trunk
         rows<-allrows
      }
      n1<-log2(rows)
      n2<-order(round(n1-floor(n1),digits=digits))
      trunk<-subtrunk[n2,]
      #windows()
      if (nrow(trunk)>1)
      {  #Create frame.
         rtframe<-data.frame(matrix(0,nrow(trunk),8))
         colnames(rtframe)<-c("var","n","wt","dev","yval",
                               "complexity","ncompete","nsurrogate")
         rownames(rtframe)<-rownames(trunk)

         varnr=as.numeric(rownames(trunk))
         vars=trunk[,1]
         vars[trunk[,6]!=""]<-"<leaf>"
         vars[trunk[,6]==""]<-trunk[varnr%%2==0,1]

         rtframe[,1]<-factor(vars)
         rtframe[,2]<-trunk[,4]
         rtframe[,3]<-trunk[,4]
         rtframe[,4]<-0
         rtframe[,5]<-round(as.numeric(trunk[,5]),digits=digits)

         #Create splits.
         rtsplits<-matrix(0,nrow=length(vars[varnr%%2==0]),ncol=5,
                           dimnames=list(vars[vars!="<leaf>"],
                           c("count","ncat","improve","index","adj")))
         rtsplits[,1]<-as.numeric(trunk[vars!="<leaf>",4])
         # -1 means <, 2 means =.
         rtsplits[,2]<- -1
         rtsplits[,4]<-round(as.numeric(trunk[varnr%%2==0,3]),digits=digits)

         rpObj1<-list(frame=rtframe,splits=rtsplits)
         rpObj1$functions$text<-mytext
         class(rpObj1)<-"rpart"
         plot(rpObj1,branch=0,margin=.10,uniform=TRUE)
         text(rpObj1,all=TRUE,use.n=TRUE,fancy=TRUE)
      }
      else
      {  # Create empty plot.
         #plot.new()
         plot(c(0,800),c(0,800),type="n",xlab="",ylab="",axes=FALSE)
         rect(350,600,550,700)
        text(450,670,round(as.numeric(trunk[1,5]),digits=digits),adj=0.5)
         text(450,630,paste("n=",trunk[1,4]),adj=0.5)
      }
      if (ncat!=1) title(main=paste(trunk[1,1],trunk[1,2],trunk[1,3]))
      par(ask=TRUE)
   }
}
