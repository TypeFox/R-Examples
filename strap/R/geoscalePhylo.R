geoscalePhylo<-function(tree, ages, direction="rightwards",units=c("Period", "Epoch", "Age"), boxes="Age", tick.scale="myr", user.scale, cex.age=0.3, cex.ts=0.3, cex.tip=0.3, width=1, label.offset,ts.col=TRUE, vers="ICS2013", x.lim, quat.rm=FALSE,erotate,arotate,urotate,...){
  
  options <- as.list(match.call())
   if(any(names(options) == "type")){
     if(all(options$type != c("phylogram","cladogram","p","c"))){
       return(cat("type must be either 'phylogram' or 'cladogram'."))
     }
   }
  
   if(all(direction != c("rightwards","upwards"))){
      return(cat("direction must be either 'rightwards' or 'upwards', here set to 'rightwards'."))
    }
  
  if(is.null(tree$root.time)){     
    return(cat("\n tree$root.time is missing, check tree is time scaled."))
  } else {root.age <- tree$root.time}
  
  if(boxes == "User" && any(units != "User")){
    boxes <- "no"
  }
  
  if(all(boxes != units)){
    boxes <- "no"
  }
  
  if(tick.scale == "User" && all(units != "User")){
    tick.scale <- "myr"
  }
   
  if(missing(ages) == FALSE){
      ranges <- TRUE
  } else{
    ranges <- FALSE
  }
  
  if(missing(user.scale) & any(units == "User")){
    units <- units[units != "User"]
    cat("\n user.scale not provided, 'Other' removed from units.")
  }
  
  if(missing(ages) == FALSE){
    ages<-ages[tree$tip.label,]    
  }
  
  if(any(units == "User") & !missing(user.scale)){   
    Midpoint <- matrix(ncol=1,nrow=length(user.scale[,1]))
      Midpoint[,1] <- (user.scale[,"Start"] + user.scale[,"End"])/2
        user.scale <- cbind(user.scale,Midpoint)  
  }
  
  if(all(units != "Age") && boxes == "Age"){
    boxes <- "no"
  }
  
  units <- paste(toupper(substring(units,1,1)),substring(units,2),sep="")
    
  # Standardizing the names of temporal units
  
  units[units == "Eonothem"] <- "Eon"
  units[units == "Erathem"] <- "Era"
  units[units == "Series"] <- "Epoch"
  units[units == "System"] <- "Period"  
  units[units == "Stage"] <- "Age"
    units <- unique(units)
  
  boxes[boxes == "Eonothem"] <- "Eon"
  boxes[boxes == "Erathem"] <- "Era"
  boxes[boxes == "Series"] <- "Epoch"
  boxes[boxes == "System"] <- "Period"  
  boxes[boxes == "Stage"] <- "Age"
      
  if(length(units) == 1){
    ts.width=0.15
  } else if(length(units) == 2){
    ts.width=0.2
  } else if(length(units) >= 3){
    ts.width=0.25
  }

  if(ranges == TRUE && missing(ages) == FALSE){
     missing.tip.names <- setdiff(tree$tip.label,row.names(ages))
      if(length(missing.tip.names) > 0){            
        cat(paste("\n",missing.tip.names,"not present in ages file, ranges set to FALSE"))
        cat("\n ranges set to FALSE")
          ranges <- FALSE
      }    
  }
   
  timescales <- NULL
   data(timescales,envir=environment())
    timescale <- timescales[[vers]]
      if(quat.rm == TRUE){
       timescale[(timescale[,"Midpoint"] < 3),"Name"] <- NA
      }
  
  tscale.data<-matrix(ncol=3,nrow=6)
    colnames(tscale.data) <-c("srt","Depth","size")
    rownames(tscale.data) <-c("Eon","Era","Period","Epoch","Age","User")
      if(direction == "upwards"){
        tscale.data[,"srt"] <- c(90,90,90,0,0,0)
      } else tscale.data[,"srt"] <- c(0,0,0,90,90,90)
      
      tscale.data[,"Depth"] <- c(1,1,1,2,3.5,3.5)
      tscale.data[,"size"] <- c(1,1,1,0.8,0.8,0.8)

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
    tscale.data["Epoch","srt"] <- erotate
  }
  if(!missing(arotate)){
    tscale.data["Age","srt"] <- arotate
  }
  if(!missing(urotate)){
    tscale.data["User","srt"] <- urotate
  }
  
  ## GEOLOGICAL RANGES OF TAXA

  units<-rownames(tscale.data)[sort(match(units,rownames(tscale.data)),decreasing=T)] 
    
  if(!missing(x.lim) ){
    x.lim <- sort(root.age - x.lim)
  } else if(ranges == TRUE && !missing(ages) && missing(x.lim)){
    x.lim <- (root.age - min(ages)) + diff(range(ages))*0.05
  } else {
    x.lim <- NULL
  }
  
    timescale<-timescale[order(timescale[,1],decreasing=T),]
    	timescale.rescaled <- timescale
  		timescale.rescaled[,c("Start","End","Midpoint")] <- root.age - timescale[,c("Start","End","Midpoint")]

  first_la <- tree$root.time - dist.nodes(tree)[1,Ntip(tree)+1]
	
    if(ranges==TRUE && missing(ages) == FALSE){
      offset<-array(dim=length(tree$tip.label),data=1)
       offset.correction <- diff(range(ages)) * 0.01 
        taxon.ranges <- root.age - ages[,c("FAD","LAD")]
         if(first_la != ages[1,"LAD"]){
          if(!missing(label.offset)){
           offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+label.offset)
          } else {
           offset <- array(dim=length(ages[,"FAD"]),data=(ages[,"FAD"] - ages[,"LAD"])+offset.correction)
          }
         }
    } else if(!missing(label.offset)){
      offset <- label.offset 
    } else {
      offset = 1
    }

  ### ADDING A TIMESCALE

  if(tick.scale != "n" | tick.scale != "no"){
    if(tick.scale == "myr" | is.numeric(tick.scale)){    
      scale.ticks=1
      
      if(is.numeric(tick.scale)){
       scale.ages <- tick.scale
      } else {scale.ages=10}
        
        tick.position <- root.age - seq(0,4600,scale.ticks)
          age.name <- seq(0,4600,scale.ages)
          age.position <- root.age - age.name
            lwd<-c(1,0.5,0.5,0.5,0.5,0.7,0.5,0.5,0.5,0.5)
            col<-c("black","grey","grey","grey","grey")
    }

    if(tick.scale != "myr" & is.numeric(tick.scale) == FALSE){
      age.name<-subset(timescale,timescale[,"Type"] == tick.scale & timescale[,"Source"] == "ICS")
        age.name<-sort(unique(c(age.name[,"Start"],age.name[,"End"])))
        	age.position<- tick.position <- root.age - age.name
            lwd=1
              col="black" 
    } 
    
    if(tick.scale == "User"){     
      age.name <- sort(unique(c(user.scale[,"Start"],user.scale[,"End"])))
       age.position <- tick.position <- root.age - age.name
          lwd=1
            col="black"        
      }
   }
  
  # Plotting the tree vertically
  if(direction == "upwards"){
    
    # ADDING THE TIME SCALE
    par(fig=c(0,ts.width,0,1))
     par(mar=c(3,1,2,5))
    
        plot.phylo(tree,plot=FALSE,no.margin=T,y.lim=x.lim,direction="upwards",cex=cex.tip,...)
          
          timescale.rescaled.names <- timescale.rescaled
           timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[3],]
            timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[3],"Start"] <- par()$usr[3]
    
          if(min(timescale.rescaled.names[,"End"]) < par()$usr[4]){
            timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[4],]
             timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[4],"End"] <- par()$usr[4]
          }
                timescale.rescaled.names[,"Midpoint"] <-(timescale.rescaled.names[,"Start"] + timescale.rescaled.names[,"End"])/2
    
    
          unit.depths <- tscale.data[units,"Depth"]
            if(tick.scale == "n" | tick.scale == "no"){
              unit.depths <- c(unit.depths,0.5)          
            } else if(length(units) <= 3){
              unit.depths <- c(unit.depths,2)
            } else if(length(units) > 3) {
              unit.depths <- c(unit.depths,2)
            }

           unit.depths <- cumsum(unit.depths/sum(unit.depths))
            unit.depths<-c(par()$usr[2],par()$usr[2]-(unit.depths*(par()$usr[2]-par()$usr[1])))
             
         depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
        if(tick.scale != "n" && tick.scale != "no"){
          text((unit.depths[length(unit.depths)]+depth*0.3),age.position,age.name,cex=cex.age,srt=0)           
           segments((unit.depths[length(unit.depths)-1]),tick.position,(unit.depths[length(unit.depths)]+depth*0.75),tick.position,lwd=lwd,col=col)
        }
    
    for(t in 1:length(units)){
      
      if(units[t] == "User"){
        tscale<-user.scale
         tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
         tscale.names <- tscale
      } else {
        tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
        tscale.names<-subset(timescale.rescaled.names,timescale.rescaled.names[,"Type"] == units[t])
      }
          
      if(ts.col == TRUE & units[t] != "User"){
        rect(unit.depths[t],tscale[,"Start"],unit.depths[t+1],tscale[,"End"],col=rgb(tscale[,"Col_R"],tscale[,"Col_G"],tscale[,"Col_B"],maxColorValue=255))
      } else rect(unit.depths[t],tscale[,"Start"],unit.depths[t+1],tscale[,"End"],col="white")
      text((unit.depths[t] + unit.depths[t+1])/2,tscale.names[,"Midpoint"],tscale.names[,"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
    
    # ADDING THE PHYLOGENY
    par(fig=c(ts.width,1,0,1),new=T)
     par(mar=c(3,0,2,2))
    
      plot.phylo(tree,plot=FALSE,no.margin=T,y.lim=x.lim,direction="upwards",cex=cex.tip,...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    
      if (!missing(boxes) && boxes != "no" && boxes !="n"){
       if(boxes == "User"){
         tscale <- root.age - user.scale[,c("Start","End")]
       } else {
         tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)}
       
           rect(par()$usr[3],tscale[,"Start"],par()$usr[4],tscale[,"End"],col=c("grey90","white"),border=NA)
       }
        
      par(fig=c(ts.width,1,0,1),new=T)
        par(mar=c(3,0,2,2))
    
        plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,y.lim=x.lim,cex=cex.tip,direction="upwards",...)
          if (ranges == TRUE){
            par(lend=1); segments(lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"FAD"],lastPP$xx[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],col="black",lwd=width*2)
          }
    
    if(units[1] == "User"){
      segments(par()$usr[1],min(root.age - user.scale[,"Start"]),par()$usr[1],max(root.age - user.scale[,"End"]))
    } else {      
      segments(par()$usr[1],min(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"Start"]),par()$usr[1],max(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"End"]))
    }
       
  } else{
    # Plotting the tree horizontally
    
  # PLOT 1 - TIMESCALE
    par(fig=c(0,1,0,ts.width))
    par(mar=c(1,3,0,2))

	    plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim,direction="rightwards",cex=cex.tip,...)
    
        timescale.rescaled.names <- timescale.rescaled
         timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[1],]
          timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[1],"Start"] <- par()$usr[1]
    
        if(min(timescale.rescaled.names[,"End"]) < par()$usr[2]){
          timescale.rescaled.names <- timescale.rescaled.names[timescale.rescaled.names[,"Start"] < par()$usr[2],]
           timescale.rescaled.names[timescale.rescaled.names[,"End"] > par()$usr[2],"End"] <- par()$usr[2]
        }
            timescale.rescaled.names[,"Midpoint"] <-(timescale.rescaled.names[,"Start"] + timescale.rescaled.names[,"End"])/2
    
       unit.depths <- tscale.data[units,"Depth"]
        if(tick.scale == "n" & tick.scale == "no"){
          unit.depths <- c(unit.depths,0.5)          
        } else if(length(units) <= 3){
          unit.depths <- c(unit.depths,2)
        } else if(length(units) > 3) {
          unit.depths <- c(unit.depths,2)
        }
     
         unit.depths <- cumsum(unit.depths/sum(unit.depths))
          unit.depths<-c(par()$usr[4],par()$usr[4]-(unit.depths*(par()$usr[4]-par()$usr[3])))
    
       depth <- unit.depths[length(unit.depths) - 1] - unit.depths[length(unit.depths)]
    
    if(tick.scale != "n" && tick.scale != "no"){
      text(age.position,(unit.depths[length(unit.depths)]+depth*0.3),age.name, cex=cex.age,srt=90)
  	  segments(tick.position,(unit.depths[length(unit.depths)-1]),tick.position,(unit.depths[length(unit.depths)]+depth*0.6),lwd=lwd,col=col)
    }

    for(t in 1:length(units)){
 	
 	    if(units[t] == "User"){
 	    	tscale<-user.scale
  	    	tscale[,c("Start","End","Midpoint")] <- root.age - tscale[,c("Start","End","Midpoint")]
     	    	tscale.names <- tscale
      } else {
 	      tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == units[t])
 	      tscale.names<-subset(timescale.rescaled.names,timescale.rescaled.names[,"Type"] == units[t])
 	    }
 	 	
 	    if(ts.col == TRUE & units[t] != "User"){rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col=rgb(tscale[,"Col_R"],tscale[,"Col_G"],tscale[,"Col_B"],maxColorValue=255))
 		    } else rect(tscale[,"Start"],unit.depths[t],tscale[,"End"],unit.depths[t+1],col="white")
          text(tscale.names[,"Midpoint"],(unit.depths[t] + unit.depths[t+1])/2,tscale.names[,"Name"],cex=cex.ts*tscale.data[match(units[t],rownames(tscale.data)),"size"],srt=tscale.data[match(units[t],rownames(tscale.data)),"srt"])
    }
 
  ## PLOT 2: PHYLOGENY

    par(fig=c(0,1,ts.width,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,plot=FALSE,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
      if (!missing(boxes) && boxes != "no" && boxes !="n"){
        if(boxes == "User"){
          tscale <- root.age - user.scale[,c("Start","End")]
        } else {tscale<-subset(timescale.rescaled,timescale.rescaled[,"Type"] == boxes)} 
  	        rect(tscale[,"Start"],par()$usr[3],tscale[,"End"],par()$usr[4],col=c("grey90","white"),border=NA)}

    par(fig=c(0,1,ts.width,1),new=T)
    par(mar=c(0,3,2,2))
  
      plot.phylo(tree,label.offset=offset,edge.width=width,no.margin=T,x.lim=x.lim,cex=cex.tip,...)
        if (ranges == TRUE){
          par(lend=1); segments(taxon.ranges[,"FAD"],lastPP$yy[c(1:length(tree$tip.label))],taxon.ranges[,"LAD"],lastPP$yy[c(1:length(tree$tip.label))],col="black",lwd=width*2)
        }
  
  if(units[1] == "User"){
    segments(min(root.age - user.scale[,"Start"]),par()$usr[3],max(root.age - user.scale[,"End"]),par()$usr[3])
  } else {      
    segments(min(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"Start"]),par()$usr[3],max(timescale.rescaled[timescale.rescaled[,"Type"] == units[1],"End"]),par()$usr[3])
  }
  
  }
}