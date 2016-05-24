geoscaleBox<-function (data, ages,units = c("Age", "Epoch", "Period"), tick.scale = "myr", 
          boxes = "Age", abbrev, cex.age = 0.3, cex.ts = 0.4, cex.pt = 1, age.lim = NULL,  
          data.lim = NULL, box.width=1, user.scale, ts.col = TRUE, ts.width = 0.3, label, vers="ICS2015",
          no.axis=FALSE, notch=FALSE,log=FALSE,color,direction="horizontal",erotate,arotate,urotate,...) 
{
  if(missing(ages)){
    return(cat("ages must be provided for the geological position for each box."))
  }
  
  if(missing(label) | no.axis==TRUE){
  	label <- ""
  }
    
    if(all(direction != c("horizontal","vertical"))){
      return(cat("direction must be either 'horizontal' or 'vertical', here set to 'horizonal'."))
      direction <- "horizontal"
    }
  
    
   if(tick.scale == "User" & missing(user.scale)){
 	cat("user.scale is missing, tick.scale set to 'myr'")
  	 tick.scale <- "myr"
   }	         
   
    if(all(direction != c("horizontal","vertical"))){
      return(cat("direction must be either 'horizontal' or 'vertical', here set to 'horizonal'."))
      direction <- "horizontal"
    }
               
  ages <- as.numeric(ages)
  
  statistics <- boxplot(data,plot=FALSE)
  
  if(length(statistics$out) == 0){
    outlier <- FALSE
    if(missing(data.lim)){
     data.lim <- c(min(statistics$stats[1,],na.rm=TRUE),max(statistics$stats[5,],na.rm=TRUE))
    }    
  } else{
    outlier <- TRUE
     if(missing(data.lim)){
      data.lim <- c(min(statistics$stats[1,],na.rm=TRUE),max(statistics$out,na.rm=TRUE))
     }
  }
  
  if(missing(color)){
    color <- "white"
  }
  
  timescales <- NULL
   data(timescales,envir=environment())
    timescale<-timescales[[vers]]

 units <- tolower(units)
  units <- paste(toupper(substring(units,1,1)),substring(units,2),sep="")

  convert<-matrix(ncol=2,nrow=5)
   convert[,1] <- c("Eonothem","Erathem","System","Series","Stage")
   convert[,2] <- c("Eon","Era","Period","Epoch","Age")

   for(c in 1:length(units)){
    if (is.na(convert[match(units[c],convert[,1]),2]) == FALSE){
     units[c] <- convert[match(units[c],convert[,1]),2]
    }
   }

   units<-unique(units)

if(!missing(abbrev)){
	
 abbrev <- tolower(abbrev)
  abbrev <- paste(toupper(substring(abbrev,1,1)),substring(abbrev,2),sep="")	
   abbrev[abbrev == "User"] <- NA

 if(any(abbrev == "All")){
 	abbrev <- c("Eon","Era","Period","Epoch","Age")
 }

 abbrev <- paste(toupper(substring(abbrev,1,1)),substring(abbrev,2),sep="")

  for(c in 1:length(abbrev)){
    if (is.na(convert[match(abbrev[c],convert[,1]),2]) == FALSE){
     abbrev[c] <- convert[match(abbrev[c],convert[,1]),2]
    }
  }

  abbrev<-unique(abbrev)
  
   levels(timescale[,"Name"]) <- c(levels(timescale[,"Name"]),as.character(timescale[,"Abbrev"]))
  
  for(a in 1:length(abbrev)){
  	timescale[timescale[,"Type"] == abbrev[a],"Name"] <- as.character(timescale[timescale[,"Type"] == abbrev[a],"Abbrev"])
  }

}

   if (any(units == "User") && missing(user.scale) == TRUE) {
    print("You need to specify a file for the argument user.scale, option for 'User' is removed")
    units <- subset(units, units != "User")
   }
    tscale_data <- matrix(ncol = 3, nrow = 6)
     colnames(tscale_data) <- c("srt", "Depth", "size")
     rownames(tscale_data) <- c("Eon", "Era", "Period", "Epoch","Age", "User")
	   if(direction == "vertical"){
	   	tscale_data[, "srt"] <- c(90, 90, 90, 0, 0, 0) 
	   } else{
	    tscale_data[, "srt"] <- c(0, 0, 0, 90, 90, 90)
	   }	
	  
	    tscale_data[, "Depth"] <- c(0.5, 0.5, 0.5, 1, 1, 1)
        tscale_data[, "size"] <- c(1, 1, 1, 0.8, 0.7, 0.5)
        tscale <- timescale[order(timescale[, 1], decreasing = T), ]
         units <- rownames(tscale_data)[sort(match(units, rownames(tscale_data)),decreasing = T)]
                              
     ## ROTATING THE NAMES
  
  if(!missing(erotate) && !is.numeric(erotate)){
    return(cat("\n value for protate must be numeric."))
  }
  if(!missing(arotate) && !is.numeric(arotate)){
    return(cat("\n value for arotate must be numeric."))
  }
  if(!missing(urotate) && !is.numeric(urotate)){
    return(cat("\n value for urotate must be numeric."))
  }
  if(!missing(erotate)){
    tscale_data["Epoch","srt"] <- erotate
  }
  if(!missing(arotate)){
    tscale_data["Age","srt"] <- arotate
  }
  if(!missing(urotate)){
    tscale_data["User","srt"] <- urotate
  }

	if(tick.scale == "n"){
	 tick.scale <- "no"
	}
                             
  if(tick.scale != "no"){
  	if(tick.scale == "myr"){
  	 scale.ticks <- 1
  	 scale.ages <- 10
  	  ticks <- seq(0, 4600, scale.ticks)
  	  time <- seq(0, 4600, scale.ages)
  	   	tick.width <- c(1, 0.5, 0.5, 0.5, 0.5, 0.7, 0.5, 0.5, 0.5,0.5)
  	   	tick.col <- c("black", "grey", "grey", "grey", "grey")
  	} else if(any(tick.scale == c("User","Eon","Era","Period","Epoch","Age"))){
  		if(tick.scale == "User"){
  			time <- user.scale
  		} else time <- subset(timescale,timescale[,"Type"] == tick.scale)
  			time <- unique(c(time[,1],time[,2]))
  			 ticks <- time
  			  tick.width <- 1
  			  tick.col <- "black"
  	} else if(is.numeric(tick.scale)){
  		ticks <- seq(0,4600,tick.scale)
  		time <- seq(0,4600,tick.scale)
  		 tick.width <- 1
  		 tick.col <- "black"
  	}
  }                       
                             
  if(missing(age.lim)){
    age.lim <- as.numeric(range(ages))
  }

 # Plotting horizontally

 if(direction == "horizontal") {
  par(fig = c(0, 1, 0, (ts.width)))
  par(mar = c(1, 4, 0, 2))

    plot(as.numeric(ages),statistics$stats[3,],ylim=data.lim,xlim=age.lim,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

     timescale.names <- timescale
      timescale.names <- timescale.names[timescale.names[,"End"] < max(c(par()$usr[1],par()$usr[2])) & timescale.names[,"Start"] > min(c(par()$usr[1],par()$usr[2])),]
        timescale.names[timescale.names[,"End"] < min(c(par()$usr[1],par()$usr[2])),"End"] <- min(c(par()$usr[1],par()$usr[2]))
        timescale.names[timescale.names[,"Start"] > max(c(par()$usr[1],par()$usr[2])),"Start"] <- max(c(par()$usr[1],par()$usr[2]))
  
          timescale.names[,"Midpoint"] <- (timescale.names[,"Start"] + timescale.names[,"End"]) /2     
          timescale.names[,"Range"] <- timescale.names[,"Start"] - timescale.names[,"End"]     
  
            for(t in 1:length(timescale.names[,1])){
              if(timescale.names[t,"Range"] < timescale[rownames(timescale.names)[t],"Range"]*0.3){
                timescale.names[t,"Name"] <- NA 
              }
            }
  
      if(log == TRUE){
        ylower <- 0.0000001
      } else{
        ylower <- par()$usr[3]        
      }
        yupper <- par()$usr[4]
       
       vals <- tscale_data[units, "Depth"]
       vals <- c(vals, 0.8)
       vals <- cumsum(vals/sum(vals))
       vals <- c(par()$usr[4], par()$usr[4] - (vals * (par()$usr[4] - par()$usr[3])))
       val <- vals[length(vals) - 1] - vals[length(vals)]
  
        if (tick.scale != "n") {
         text(time, (vals[length(vals)] + val * 0.3), time, cex = cex.age,srt = 90)
          segments(ticks, (vals[length(vals) - 1]), ticks, (vals[length(vals)] + val * 0.75), lwd = tick.width)
        }
  
        for (t in 1:length(units)) {
         if (units[t] == "User") {
          tscale_n <- user.scale
         } else tscale_n <- subset(tscale, tscale[, "Type"] == units[t])
            if (ts.col == TRUE & units[t] != "User") {
             rect(tscale_n[, "Start"], vals[t], tscale_n[, "End"], 
             vals[t + 1], col = rgb(tscale_n[, "Col_R"], tscale_n[,"Col_G"], tscale_n[, "Col_B"], maxColorValue = 255))
            } else rect(tscale_n[, "Start"], vals[t], tscale_n[, "End"],vals[t + 1], col = "white")
            	
            if(units[t] == "User"){
       		 text(tscale_n[, "Midpoint"], (vals[t] + vals[t + 1])/2, tscale_n[, "Name"], cex = cex.ts * tscale_data[match(units[t],rownames(tscale_data)), "size"], srt = tscale_data[match(units[t],rownames(tscale_data)), "srt"])
            } else{
            
             text(timescale.names[timescale.names[,"Type"] == units[t], "Midpoint"], (vals[t] + vals[t + 1])/2, 
             timescale.names[timescale.names[,"Type"] == units[t], "Name"], cex = cex.ts * tscale_data[match(units[t],rownames(tscale_data)), "size"], srt = tscale_data[match(units[t], 
                                                                                                                         rownames(tscale_data)), "srt"])
         }
        } 
  
  par(fig = c(0, 1, (ts.width), 1), new = T)
  par(mar = c(0, 4, 2, 2))   

   if(log == TRUE){
    plot(as.numeric(ages),statistics$stats[3,],ylim=data.lim,xlim=age.lim,type="n",bty="n",yaxt="n",xaxt="n",ylab=label,xlab="",log="y")    
   } else{
    plot(as.numeric(ages),statistics$stats[3,],ylim=data.lim,xlim=age.lim,type="n",bty="n",yaxt="n",xaxt="n",ylab=label,xlab="")    
   }

   if (!missing(boxes)) {
    if (boxes == "User") {
      tscale_x <- user.scale
    } else tscale_x <- subset(tscale, tscale[, "Type"] == boxes)
        rect(tscale_x[, "Start"], ylower, tscale_x[, "End"], 
        yupper, col = c("grey90", "white"), border = NA)
    }

    if(notch == TRUE){
    
     segments(ages,statistics$stats[1,],ages,statistics$stats[5,],lty=2)
        
      poly_x <- matrix(ncol=length(statistics$stats[1,]),nrow=11)
      poly_y <- matrix(ncol=length(statistics$stats[1,]),nrow=11)  
       for(x in 1:length(statistics$stats[1,])){
        poly_x[c(1,2,4,5),x] <- ages[x] + box.width
        poly_x[3,x] <- ages[x]+(box.width*0.5)
        poly_x[8,x] <- ages[x]-(box.width*0.5)
        poly_x[c(6,7,9,10),x] <- ages[x] - box.width
    
        poly_y[1,] <- statistics$stats[2,]
        poly_y[10,] <- statistics$stats[2,]    
        poly_y[2,] <- statistics$conf[1,]
        poly_y[9,] <- statistics$conf[1,]      
        poly_y[8,] <- statistics$stats[3,]
        poly_y[3,] <- statistics$stats[3,]
        poly_y[4,] <- statistics$conf[2,]
        poly_y[7,] <- statistics$conf[2,]      
        poly_y[5,] <- statistics$stats[4,]
        poly_y[6,] <- statistics$stats[4,]
       }
     
      for(b in 1:length(statistics$stats[1,])){
        polygon(poly_x[,b],poly_y[,b],col=color)
    
          segments(ages[b]-(box.width*0.5),statistics$stats[3,b],ages[b]+(box.width*0.5),statistics$stats[3,b])
          segments(ages[b]-(box.width*0.3),statistics$stats[5,b],ages[b]+(box.width*0.3),statistics$stats[5,b])
          segments(ages[b]-(box.width*0.3),statistics$stats[1,b],ages[b]+(box.width*0.3),statistics$stats[1,b]) 
      }
  } else{

    segments(ages,statistics$stats[1,],ages,statistics$stats[5,],lty=2)
    
    for(b in 1:length(statistics$stats[1,])){
     rect(ages[b]-box.width,statistics$stats[2,b],ages[b]+box.width,statistics$stats[4,b],col=color)
     segments(ages[b]-box.width,statistics$stats[3,b],ages[b]+box.width,statistics$stats[3,b])
    
      segments(ages[b]-(box.width*0.3),statistics$stats[5,b],ages[b]+(box.width*0.3),statistics$stats[5,b])
      segments(ages[b]-(box.width*0.3),statistics$stats[1,b],ages[b]+(box.width*0.3),statistics$stats[1,b])
    }
  }
    if(outlier == TRUE){
      points(ages[statistics$group],statistics$out,cex=cex.pt,...)
    }
    
    if(no.axis == FALSE){
    	axis(2)
    }
    
    # Plotting vertically
    
    } else if(direction == "vertical"){
    	
       if(missing(age.lim)){
  	    age.lim <- rev(range(ages))
       }	
    	
      par(fig = c(0, ts.width, 0, 1))
      par(mar = c(4, 0, 2, 0))

    plot(statistics$stats[3,],as.numeric(ages),ylim=age.lim,xlim=data.lim,type="n",bty="n",xaxt="n",yaxt="n",xlab="",ylab="")

     timescale.names <- timescale
      timescale.names <- timescale.names[timescale.names[,"End"] < max(c(par()$usr[3],par()$usr[4])) & timescale.names[,"Start"] > min(c(par()$usr[3],par()$usr[4])),]
        timescale.names[timescale.names[,"End"] < min(c(par()$usr[3],par()$usr[4])),"End"] <- min(c(par()$usr[3],par()$usr[4]))
        timescale.names[timescale.names[,"Start"] > max(c(par()$usr[3],par()$usr[4])),"Start"] <- max(c(par()$usr[3],par()$usr[4]))
  
          timescale.names[,"Midpoint"] <- (timescale.names[,"Start"] + timescale.names[,"End"]) /2     
          timescale.names[,"Range"] <- timescale.names[,"Start"] - timescale.names[,"End"]     
  
            for(t in 1:length(timescale.names[,1])){
              if(timescale.names[t,"Range"] < timescale[rownames(timescale.names)[t],"Range"]*0.3){
                timescale.names[t,"Name"] <- NA 
              }
            }
  
      if(log == TRUE){
        ylower <- 0.0000001
      } else{
        ylower <- par()$usr[1]        
      }
        yupper <- par()$usr[2]
       
       vals <- tscale_data[units, "Depth"]
       vals <- c(vals, 0.8)
       vals <- cumsum(vals/sum(vals))
       vals <- c(par()$usr[2], par()$usr[2] - (vals * (par()$usr[2] - par()$usr[1])))
       val <- vals[length(vals) - 1] - vals[length(vals)]
  
        if (tick.scale != "n") {
         text((vals[length(vals)] + val * 0.3), time, time, cex = cex.age,srt = 0)
          segments((vals[length(vals) - 1]), ticks, (vals[length(vals)] + val * 0.75), ticks, lwd = tick.width)
        }
  
        for (t in 1:length(units)) {
         if (units[t] == "User") {
          tscale_n <- user.scale
         } else tscale_n <- subset(tscale, tscale[, "Type"] == units[t])
            if (ts.col == TRUE & units[t] != "User") {
             rect(vals[t], tscale_n[, "Start"],vals[t + 1],  tscale_n[, "End"], 
              col = rgb(tscale_n[, "Col_R"], tscale_n[,"Col_G"], tscale_n[, "Col_B"], maxColorValue = 255))
            } else rect(vals[t], tscale_n[, "Start"], vals[t + 1], tscale_n[, "End"], col = "white")
            	
            if(units[t] == "User"){
       		 text((vals[t] + vals[t + 1])/2,tscale_n[, "Midpoint"], tscale_n[, "Name"], cex = cex.ts * tscale_data[match(units[t],rownames(tscale_data)), "size"], srt = tscale_data[match(units[t],rownames(tscale_data)), "srt"])
            } else{
            
             text((vals[t] + vals[t + 1])/2,timescale.names[timescale.names[,"Type"] == units[t], "Midpoint"], 
             timescale.names[timescale.names[,"Type"] == units[t], "Name"], cex = cex.ts * tscale_data[match(units[t],rownames(tscale_data)), "size"], srt = tscale_data[match(units[t], rownames(tscale_data)), "srt"])
         }
        } 
  
  par(fig = c(ts.width, 1, 0, 1), new = T)
  par(mar = c(4, 0, 2, 1))   

   if(log == TRUE){
    plot(statistics$stats[3,],as.numeric(ages),ylim=age.lim,xlim=data.lim,type="n",bty="n",yaxt="n",xaxt="n",xlab=label,ylab="",log="y")    
   } else{
    plot(statistics$stats[3,],as.numeric(ages),ylim=age.lim,xlim=data.lim,type="n",bty="n",yaxt="n",xaxt="n",xlab=label,ylab="")    
   }

   if (!missing(boxes)) {
    if (boxes == "User") {
      tscale_x <- user.scale
    } else tscale_x <- subset(tscale, tscale[, "Type"] == boxes)
        rect(ylower, tscale_x[, "Start"], yupper, tscale_x[, "End"], col = c("grey90", "white"), border = NA)
    }

    if(notch == TRUE){
    
     segments(statistics$stats[1,],ages,statistics$stats[5,],ages,lty=2)
        
      poly_x <- matrix(ncol=length(statistics$stats[1,]),nrow=11)
      poly_y <- matrix(ncol=length(statistics$stats[1,]),nrow=11)  
       for(x in 1:length(statistics$stats[1,])){
        poly_y[c(1,2,4,5),x] <- ages[x] + box.width
        poly_y[3,x] <- ages[x]+(box.width*0.5)
        poly_y[8,x] <- ages[x]-(box.width*0.5)
        poly_y[c(6,7,9,10),x] <- ages[x] - box.width
    
        poly_x[1,] <- statistics$stats[2,]
        poly_x[10,] <- statistics$stats[2,]    
        poly_x[2,] <- statistics$conf[1,]
        poly_x[9,] <- statistics$conf[1,]      
        poly_x[8,] <- statistics$stats[3,]
        poly_x[3,] <- statistics$stats[3,]
        poly_x[4,] <- statistics$conf[2,]
        poly_x[7,] <- statistics$conf[2,]      
        poly_x[5,] <- statistics$stats[4,]
        poly_x[6,] <- statistics$stats[4,]
       }
     
      for(b in 1:length(statistics$stats[1,])){
        polygon(poly_x[,b],poly_y[,b],col=color)
    
          segments(statistics$stats[3,b],ages[b]-(box.width*0.5),statistics$stats[3,b],ages[b]+(box.width*0.5))
          segments(statistics$stats[5,b],ages[b]-(box.width*0.3),statistics$stats[5,b],ages[b]+(box.width*0.3))
          segments(statistics$stats[1,b],ages[b]-(box.width*0.3),statistics$stats[1,b],ages[b]+(box.width*0.3)) 
      }
  } else{

    segments(statistics$stats[1,],ages,statistics$stats[5,],ages,lty=2)
    
    for(b in 1:length(statistics$stats[1,])){
     rect(statistics$stats[2,b],ages[b]-box.width,statistics$stats[4,b],ages[b]+box.width,col=color)
     segments(statistics$stats[3,b],ages[b]-box.width,statistics$stats[3,b],ages[b]+box.width)
    
      segments(statistics$stats[5,b],ages[b]-(box.width*0.3),statistics$stats[5,b],ages[b]+(box.width*0.3))
      segments(statistics$stats[1,b],ages[b]-(box.width*0.3),statistics$stats[1,b],ages[b]+(box.width*0.3))
    }
  }
    if(outlier == TRUE){
      points(statistics$out,ages[statistics$group],cex=cex.pt,...)
    }
    
    if(no.axis == FALSE){
    	axis(1)
    }	
    	
    }
  }