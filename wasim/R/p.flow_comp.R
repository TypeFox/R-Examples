`p.flow_comp` <-
function(data, xdata=1:NROW(data),
measured=NULL,data.names=data.types$beschreibung_en[csubset],csubset = c(10,29,6,8,9,5,11),
crain=7, ylab="flow [mm/h]",c.flow_com=NULL,
l.flow_com =NULL , interflow_correction=TRUE, interflow_row=11,
baseflow_row=8, legend.position= "right",...){

  if (is.null(dim(data)))              #handle the case when only one component is passed
    dim(data)=c(1,length(data))

  #Set linewidth if not user-selected
  if (is.null(l.flow_com)){
     #default set
     if(!any(csubset != c(10,29,6,8,9,5,11))){
            l.flow_com=c(2,1,2,1,1,1,1,1)
     } else {
        l.flow_com=c(rep(1, length(csubset)))
     }
  }
  
  #Set color if not user-selected
  if (is.null(c.flow_com)){
     #default set
     if(!any(csubset != c(10,29,6,8,9,5,11))){
            c.flow_com=c("black", "pink","red", "green", "blue",
                         "yellow", "darkgreen", "orange" )
     } else {
        c.flow_com=(1:length(csubset))+2
     }
  }
  if ((!is.null(measured)) && (length(c.flow_com)==length(csubset)))       #if no colour for measured is given, use "red"
    c.flow_com=c("red",c.flow_com)
  if ((!is.null(measured)) && (length(l.flow_com)==length(csubset)))       #if no colour for measured is given, use "red"
    l.flow_com=c(2,l.flow_com)
  if ((!is.null(measured)) && (length(data.names)==length(csubset)))       #if no label for measured data specified, use "observed"
    data.names=c("observed",data.names)
   
  maxy=1.2*quantile(c(measured,data[csubset,]),probs=0.999,na.rm=TRUE)      #get maximum value of measured and modelled data

   #correct interflow since wasim outputs the sum of interflow and baseflow
   if(interflow_correction){
      if(any(csubset ==interflow_row)){   #ist interflow in Ausgabe?
          data[interflow_row,] <- data[interflow_row,] - data[baseflow_row,]
          data.names[data.names=="Interflow_QBBpQBpQI"] <- "Interflow"
      }
   }

   plot(range(xdata),c(0,maxy),type="n", xlab="", ylab=ylab, ...)

  if (!is.null(measured))
    lines(xdata,measured,lty=1, lwd=l.flow_com[1], col=c.flow_com[1])  #plot observed

  j=1
  for (i in csubset){
          cat("range", data.names[j], ":", range(data[i,]), "\n")
          lines(xdata,data[i,], col=c.flow_com[j+!is.null(measured)], lwd=l.flow_com[j+!is.null(measured)])
          j=j+1
  }

  rain.max<-max(data[crain,], na.rm=TRUE)
  lines(xdata,data[i,], col=c.flow_com[j+!is.null(measured)], lwd=l.flow_com[j+!is.null(measured)])
    plot.window(c(1,length(xdata)),c(2*rain.max ,0))
    lines(1:length(xdata),data[crain,],  type="h")
    axis(4)
    mtext(side= 4, line=2, text="rain/mm" )

    try(line_sumrain(xdata=xdata, cum_sum = cumsum(data[crain,])/10, theMax=rain.max, col="grey"))

  legend(legend.position,legend= data.names,lty=1, lwd=l.flow_com, col=c.flow_com,bty="n", cex=0.9)


}

