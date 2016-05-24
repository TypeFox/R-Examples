fan <-
  function(data = NULL, data.type="simulations", style = "fan", type = "percentile",
           probs = if(type=="percentile") seq(0.01, 0.99, 0.01) else c(0.5, 0.8, 0.95), 
           start = 1, frequency = 1, anchor = NULL, anchor.time=NULL,
           fan.col = heat.colors, alpha = if (style == "spaghetti") 0.5 else 1, 
           n.fan = NULL,
           ln = if(length(probs)<10) probs else 
             probs[round(probs,2) %in% round(seq(0.1, 0.9, 0.1),2)],
           ln.col = if(style=="spaghetti") "gray" else NULL, 
           med.ln = if(type=="interval") TRUE else FALSE, 
           med.col= "orange",
           rlab = ln, rpos = 4, roffset = 0.1, rcex = 0.8, rcol = NULL, 
           llab = FALSE, lpos = 2, loffset = roffset, lcex = rcex, lcol = rcol, 
           upplab = "U", lowlab = "L", medlab=if(type == "interval") "M" else NULL,
           n.spag = 30, 
           space = if(style=="boxplot") 1/frequency else 0.9/frequency, 
           add = FALSE, ylim = range(data)*0.8,...){
      
    #probs[round(probs,2) %in% c(0.05,0.10,0.20,0.50,0.80,0.90,0.95)]
    if(add==TRUE)
      plot(data[,1], type="n", ylim=ylim, ...)
    
    ##
    ##check inputs
    ##
    if(!(data.type %in% c("values","simulations")))
      stop("data.type must be set to one of: values, simulations")
    if(!(style %in% c("fan","boxfan","spaghetti","boxplot")))
      stop("style must be set to one of: fan, boxfan, spaghetti or boxplot")
    
    ##
    ##data classes
    ##
    if(class(data)[1]=="mts"){
      start<-start(data)[1]
      frequency<-frequency(data)
    }
    
    ##
    ##create pp and tt
    ##
    if(style=="fan" | style=="boxfan" | style=="spaghetti"){
      if(!(type %in% c("percentile","interval")))
        stop("type must be set to one of: percentile or interval")
      #ensure p is okay
      p<-probs
      if(min(p)<0 | max(p)>100)
        stop("all probs must be between 0 and 1 (or 0 and 100)")
      if(max(p)>1)
        p<-p/100
      #make p symetrical
      if(type=="percentile")
        p<-c(p,1-p)
      if(type=="interval" & data.type=="simulations")
        p <- c(p + (1-p)/2, 1 - p - (1-p)/2)
      p<-round(p,5) #i dont know why, but you need this otherwise not unique
      p<-sort(unique(p))
      
      #work out quantiles
      if(data.type=="simulations"){
        pp<-as.matrix(data)
        pp<-apply(pp,2,quantile,probs=p)
      }
      if(data.type=="values"){
        pp<-as.matrix(data)
        if(type=="percentile" & length(p)!=nrow(pp))
          stop("probs must correspond to the nrows of data if data.type==values and type is percentile")
        if(type=="interval" & length(probs)!=2*nrow(pp)){
          p<-probs
          p<-c(p + (1-p)/2, 1 - p - (1-p)/2)
          p<-sort(p)
          p <- round(p,5)
        }
      }
      n<-nrow(pp)
      if(type=="interval"){
        rownames(pp)<-paste0(rep(c(lowlab,upplab),each=n/2), 200*abs(0.5-p)  ,"%")
      }
      
      #add ancohor
      if(!is.null(anchor)){
        pp<-cbind(rep(anchor,n),pp)
      }
      #transform pp
      pp<-t(pp)
      
      
      #ts info
      if(class(data)!="zoo"){
        ppts <- ts(pp, start = start, frequency = frequency)
        tt<-time(ppts)
        tt<-as.numeric(tt)
        if(!is.null(anchor)){
          ppts <- ts(pp, start = time(lag(ppts))[1], frequency = frequency)
          tt <- time(ppts)
          tt<-as.numeric(tt)
        }
      }
      if(class(data)=="zoo"){
        tt<-time(data)
        if(!is.null(anchor))
          tt<-c(anchor.time,tt)
      }
    }
    
    ##
    ##fan colours
    ##
    if(style=="fan" | style=="boxfan"){
      #plot polygons
      if(is.null(n.fan))
        fan.col<-fan.col(floor(n/2))
      if(!is.null(n.fan))
        fan.col<-fan.col(n.fan)
      fan.col<-adjustcolor(fan.col,alpha.f=alpha)
      if(is.null(ln.col))
        ln.col<-fan.col[1]
    }
    
    ##
    ##fan plot
    ##
    if(style=="fan"){
      fan.fill<-function(ts1, ts2, tt, fan.col="grey"){
        xx <- cbind(tt,rev(tt)) 
        yy <- cbind(as.vector(ts1),rev(as.vector(ts2)))
        polygon(xx,yy, col=fan.col, border=fan.col)
      }
      #plot(cpi, type = "l", xlim = c(y0-5, y0+3), ylim = c(-2, 7))
      #plot(NULL, type = "n", xlim = c(1, 945),  ylim = range(th.mcmc), ylab = "Theta")
      #plot(net, ylim=range(net-ips$net.ci, net+ips$net.ci), type = "n")
      for(i in 1:floor(n/2)){
        fan.fill(ts1=pp[,i],ts2=pp[,n-i+1],tt=tt, fan.col=fan.col[floor(n/2)+1-i])
      }
    }   
    
    ##
    ##plot box fans
    ##
    if(style=="boxfan"){
      #single time series to use for at=time 
      x<-pp[,1]
      for(i in 1:nrow(pp)){
        for(j in 1:floor(n/2)){
          rect(xleft=tt[i]-0.5*space, ybottom=pp[i,j], xright=tt[i]+0.5*space, ytop=pp[i,n-j+1], col=fan.col[floor(n/2)+1-j], border=fan.col[floor(n/2)+1-j])
        }
      }
    }
    
    ##
    ##contour lines
    ##
    if(style=="fan" | style=="boxfan"){
      #ln=seq(5,95,15); llab=seq(5,95,15); rlab=c(80,50,20); 
      #ensure rlab will evaluate to original ln rather than altered
      ln0<-ln
      #ensure ln is okay
      if(!is.null(ln0)){
        if(min(ln0)<0 | max(ln0)>100)
          stop("all ln must be between 0 and 1 (or 0 and 100)")
        if(max(ln0)>1)
          ln0<-ln0/100
        #default lines on available pi
        if(type=="interval"){
          ln0 <- c(ln0 + (1-ln0)/2, 1 - ln0 - (1-ln0)/2)
          ln0 <-sort(ln0)
        }
        ln0<-round(ln0,5)
        #plot lines
        if(style=="fan"){
          for(i in match(ln0, p))
            lines(x=tt, y=pp[,i], col=ln.col)
        }
        if(style=="boxfan"){
          for(i in 1:nrow(pp)){
            for(j in match(ln0, p)){
              lines(x=tt[i]+c(-0.5,0.5)*space, y=rep(pp[i,j],2), col=ln.col)
            }
          }
        }
        if(is.na(sum(match(ln0,p))))
          print("some lines not plotted as conflict with precentiles given in probs")
      }
    }
    
    ##
    ##labels
    ##
    if(style=="fan" | style=="boxfan"){
      #names will be plotted in text
      if(data.type=="values" & type=="percentile")
        colnames(pp)<-paste0(p*100, "%")
      #default right text on available deciles
      if(!is.null(rlab)){
        if(min(rlab)<0 | max(rlab)>100)
          stop("all rlab must be between 0 and 1 (or 0 and 100)")
        if(max(rlab)>1)
          rlab<-rlab/100
        if(type=="interval")
          rlab<-c(rlab + (1-rlab)/2, 1 - rlab - (1-rlab)/2)
        rlab<-sort(rlab)
        rlab<-round(rlab, 5)
        for(i in match(rlab, p)){
          if(style=="fan")
            text(tt[length(tt)], pp[nrow(pp),i], colnames(pp)[i], pos=rpos, offset=roffset, cex=rcex, col=rcol)
          if(style=="boxfan")
            text(tt[length(tt)]+0.5*space, pp[nrow(pp),i], colnames(pp)[i], pos=rpos, offset=roffset, cex=rcex, col=rcol)
        }
        if(is.na(sum(match(rlab,p))))
          print("some right labels not plotted as conflict with precentiles given in probs")
      }
      if(is.numeric(llab[1]) | llab[1]==TRUE){
        if(is.numeric(llab[1])){
          if(min(llab)<0 | max(llab)>100)
            stop("all llab must be between 0 and 1 (or 0 and 100)")
          if(max(llab)>1)
            llab<-llab/100
          if(type=="interval")
            llab <- c(llab + (1-llab)/2, 1 - llab - (1-llab)/2)
          llab<-sort(llab)
          llab<-round(llab, 5)
        }
        if(llab[1]==TRUE)
          llab<-rlab
        for(i in match(llab, p)){
          if(style=="fan")
            text(tt[1], pp[1,i],  colnames(pp)[i], pos=lpos, offset=loffset, cex=lcex, col=lcol)
          if(style=="boxfan")
            text(tt[1]-0.5*space, pp[1,i], colnames(pp)[i], pos=lpos, offset=loffset, cex=lcex, col=lcol)
        }
        if(is.na(sum(match(llab,p))))
          print("some left labels not plotted as conflict with precentiles given in probs")
      }
    }
    
    #add median line and labels
    if(style=="fan" | style=="boxfan"){
      if(med.ln==TRUE & data.type=="simulations"){
        pp<-data
        pm<-apply(pp,2,median)
        if(!is.null(anchor))
          pm<-c(anchor,pm)
        if(is.null(med.col))
          med.col<-ln.col
        if(style=="fan"){
          lines(x=tt, y=pm, col=med.col)
        }
        if(style=="boxfan"){
          for(i in 1:nrow(pp)){
            #lines(x=(i-1)/frequency+c(-0.5,0.5)*space, y=rep(pm[i],2), col=med.col)
            lines(x=tt[i]+c(-0.5,0.5)*space, y=rep(pm[i],2), col=med.col)
          }
        }
        if(!is.null(rlab) & style %in% c("fan","spaghetti"))
          text(tt[length(tt)], pm[length(pm)], medlab, pos=rpos, offset=roffset, cex=rcex, col=rcol)
        if(!is.null(rlab) & style=="boxfan")
          text(tt[length(tt)]+0.5*space, pm[length(pm)], medlab, pos=rpos, offset=roffset, cex=rcex, col=rcol)
        if(any(llab==TRUE,is.numeric(llab)) & style %in% c("fan","spaghetti"))
          text(tt[1], pm[1], medlab, pos=lpos, offset=loffset, cex=lcex, col=lcol)
        if(any(llab==TRUE,is.numeric(llab)) & style=="boxfan")
          text(tt[1]-0.5*space, pm[1], medlab, pos=lpos, offset=loffset, cex=lcex, col=lcol)
      }
    }
    
    ##
    ##plot spaghetti
    ##
    if(style=="spaghetti"){
      ps<-as.matrix(data)
      n<-nrow(ps)
      ps<-ps[sample(1:n,n.spag),]
      #add ancohor
      if(!is.null(anchor))
        ps<-cbind(rep(anchor,nrow(ps)),ps)
      spag.col<-adjustcolor(ln.col,alpha.f=alpha)
      for(i in 1:nrow(ps))
        lines(x=tt, y=ps[i,], col=spag.col)
    }
    
    ##
    ##box plots
    ##
    if(style=="boxplot"){
      if(data.type=="values")
        stop(print("data must be simulations"))
      pp<-data
      n<-ncol(pp)
      if(!is.null(anchor))
        print("anchor ignored for boxplots plots")
      #single time series to use for at=time 
      p<-ts(pp[1,], start=start, frequency=frequency)
      #plot(NULL, type = "n", xlim = c(1, 10),  ylim = range(pp), ylab = "Theta")
      for(i in 1:n)
        boxplot(pp[,i], add=TRUE, at=time(p)[i], boxwex=space, xaxt = "n", yaxt = "n",...)
    }
    box()
  }


