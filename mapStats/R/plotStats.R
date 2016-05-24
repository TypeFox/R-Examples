plotStats <- function(statmats,
                      map.file,
                      d.geo.var,
                      map.geo.var,
                      ngroups,
                      separate,
                      paletteName,
                      colorVec,
                      map.label,
                      map.label.names,
                      cex.label,
                      col.label,
                      titles,
                      cex.title,
                      wt.ind = FALSE,
                      wt.label,
                      var.pretty,
                      geo.pretty,
                      by.pretty,
                      sp_layout.pars,
                      plotbyvar=TRUE,
                      num.col,
                      ...) {
 

    #correction just in case: if there is only one variable, no by.variable
    if (ncol(statmats[[1]][[1]])==2) { by.pretty <- NULL }

    map.file@data$sort.id <- 1:nrow(map.file@data)
    for (s in names(statmats)) {

	#get rid of the last Freqs table
    
        statmats[[s]] <- statmats[[s]][-length(statmats[[s]])]
  
        #if no by.var , need to change the variable name for merging
        #loop over stat names for each variable so plotvars is same for all
        if (is.null(by.pretty)) {
          for (mn in names(statmats[[1]])) {
             colnames(statmats[[s]][[mn]])[2] <- "theanalysisvariable"
               }
     

        }
      
    }



    
    num_plots <- length(unlist(statmats))
    num_vars <- length(names(statmats))
    num_stats <- length(statmats[[1]])
    
    #list to contain trellis objects
    list_of_plots <- list()
    
    separate <- as.numeric(separate)
    #0 -> same split for all variables and statistics
    #1 -> separate split by statistic but same across variables
    #2 -> separate by variable but same for each statistic within variable
    #3 -> each one separate
    
    #style parameter
    
    #dots <- list(...)
    
    #if ("style" %in% names(dots)) { styleVec <- dots$style }
    #else { styleVec <- "quantile" }
    
    
    
    #set up plotting parameters

    usebrewer <- is.null(colorVec)

    #just in case both are null
    if( usebrewer==FALSE  & is.null(paletteName) ) { 
        paletteName <- "Reds" 
        usebrewer <- TRUE
       }

    #modify ngroups and palette as necessary
    if( usebrewer==FALSE) {
        if( ! is.list(colorVec) ) { colorVec <- list(colorVec) }
        
        ngroups <- unlist(lapply(colorVec, FUN=length))
        #expand lengths
        if (separate==3) { colorVec <- rep_len(colorVec, length.out=num_plots) ; ngroups <- rep_len(ngroups, length.out=num_plots) }
        else if (separate==0) { colorVec <- rep_len(colorVec[[1]], length.out=num_plots) ; ngroups <- rep_len(ngroups[1], length.out=num_plots) }
        else if (separate==2) { colorVec <- rep(rep_len(colorVec, length.out=num_vars), each=num_stats) 
                                ngroups <- rep(rep_len(ngroups, length.out=num_vars), each=num_stats)
        }                
        else if (separate==1) { colorVec <- rep(rep_len(colorVec, length.out=num_stats), times=num_vars) 
                                ngroups <- rep(rep_len(ngroups, length.out=num_stats), times=num_vars)
        }               
    
      }    

     else {
       
       
       if (separate==3) {      paletteName <- rep_len(paletteName, length.out=num_plots) ; ngroups <- rep_len(ngroups, length.out=num_plots) }
       else if (separate==0) { paletteName <- rep_len(paletteName[[1]], length.out=num_plots) ; ngroups <- rep_len(ngroups[1], length.out=num_plots) }
       else if (separate==2) { paletteName <- rep(rep_len(paletteName, length.out=num_vars), each=num_stats) 
                               ngroups <- rep(rep_len(ngroups, length.out=num_vars), each=num_stats)
       }                
       else if (separate==1) { paletteName <- rep(rep_len(paletteName, length.out=num_stats), times=num_vars) 
                               ngroups <- rep(rep_len(ngroups, length.out=num_stats), times=num_vars)
       }         
       
       
      
    }     
          

  
    #name of columns plotted   
    plotvars <- colnames(statmats[[ 1 ]][[ 1 ]])[ 2:ncol(statmats[[ 1 ]][[ 1 ]])] 


    #generate plot headers    
     if (is.null(by.pretty)) { prettynames <- "" } 
     else { prettynames <- plotvars
            
      #rename plotvars to protect against illegal names
	    #this is just to rename in plotting, does not affect output.

	    #put an A in front 
      plotvars <- paste("A", plotvars, sep="")

	    #now replace any illegal (anything but alphanumeric, underscore, or dot) characters with A. This protects from errors

      plotvars <- gsub("[^a-zA-Z0-9._]", "A", plotvars, perl=TRUE)

            for (p in 1:length(prettynames)) {
                           
              new <- paste(by.pretty, strsplit(prettynames[p], "__.__")[[1]], sep="=")
              new <- paste(new, collapse=", ")
              prettynames[p] <- new
              }
           
           #now rename plotvars, not the headers just for safety

               for (k in 1:length(statmats))  {
                   for (j in 1:length(statmats[[ k ]])) {             
                         colnames(statmats[[ k ]][[ j ]])[ 2:(length(plotvars)+1) ] <- plotvars
                       }
               }

       }#end of renaming columns
        
    
    #color breaks
    #if want to maintain the same split for all  variables and statistics      
    
    if( separate==0 ) {
      
      #use the first element of ngroups to choose number 
      #remove NAs so don't get warnings
     # styleVec <- rep(styleVec[1], length.out=num_plots)
      tmp <- lapply(sapply(statmats, "[",i=1:num_stats), "[", i=2:(length(plotvars)+1))
      
      var_noNA <- as.vector(unlist(tmp))
      var_noNA <- var_noNA[ ! is.na(var_noNA) ]
      class_div <- classInt::classIntervals(var=var_noNA, n=ngroups[1],  ...)
      class_div <- jiggleClass(x=class_div)
      breaks <- rep_len(list(class_div$brks), length.out=num_plots)

      #in case number of groups changes
      ngroups[1] <- length(breaks) -1
      
      #allow ngroups to change      
      if( usebrewer==TRUE ) {  plotclr <- RColorBrewer::brewer.pal(n=ngroups[1], name=paletteName[1])  }
      else { plotclr <- colorVec[[1]][ 1:ngroups[1] ] }
 
          
      
      plotclr <- rep_len(list(plotclr), length.out=num_plots)
      
      
    }#end of separate
    
    else if (separate==1) {
      #different for each statistic but same across variables
      
      #styleVec <- rep(rep_len(styleVec, length.out=num_stats), times=num_vars)
     
      breaks <- list()
      plotclr <- list()

      for (k in 1:num_stats) {
        
        var_noNA <- as.vector(unlist(lapply(sapply(statmats, "[", i=k), "[", i=2:(length(plotvars)+1))))
        var_noNA <- var_noNA[ ! is.na(var_noNA) ]
        class_div <- classInt::classIntervals(var=var_noNA, n=ngroups[k],  ...)
        class_div <- jiggleClass(x=class_div)
        breaks[[ k ]] <- class_div$brks
        #allow ngroups to change
        ngroups[k] <- length(breaks[[k]]) -1
                
        if ( usebrewer==TRUE ) { plotclr[[k]] <- RColorBrewer::brewer.pal(n=ngroups[k], name=paletteName[k]) }
        else {  plotclr[[k]] <- colorVec[[k]][ 1:ngroups[k] ]  }
      
      }
      
      #now repeat for each variable
      breaks <- rep(breaks, num_vars)
      plotclr <- rep(plotclr, num_vars)
      
    
        
      }
      
     else if (separate==2) {
       
       #same across each statistic but different for each variable
              
       #styleVec <- rep(rep_len(styleVec, length.out=num_vars), each=num_stats)
       
       breaks <- list()
       plotclr <- list()
       
       for (k in 1:num_vars) {
         
         var_noNA <- as.vector(unlist(lapply(statmats[[ k ]], "[", i=2:(length(plotvars)+1))))
         var_noNA <- var_noNA[ ! is.na(var_noNA) ]
         class_div <- classInt::classIntervals(var=var_noNA, n=ngroups[k*num_stats],  ...)
         class_div <- jiggleClass(x=class_div)
         breaks[[ k ]] <- class_div$brks
         ngroups[ k*num_stats ] <- length(breaks[[ k ]]) -1
    
         
         if( usebrewer==TRUE ) {  plotclr[[k]] <- RColorBrewer::brewer.pal(n=ngroups[k*num_stats], name=paletteName[k*num_stats])  }
         else { plotclr[[k]] <- colorVec[[k*num_stats]][ 1:ngroups[ k*num_stats ] ] }
        
       }
       
         breaks <- rep(breaks, each=num_stats)
         plotclr <- rep(plotclr, each=num_stats)
     } 
         

    else if (separate==3) {
       #each plot is separate
       # styleVec <- rep_len(styleVec, length.out=num_plots)
      
        breaks <- list()
        plotclr <- list()
      
        tmp <- lapply(sapply(statmats, "[",i=1:num_stats), "[", i=2:(length(plotvars)+1))

        for (k in 1:length(tmp)) {      
            var_noNA <- as.vector( unlist( tmp[[k]]))
            var_noNA <- var_noNA[ ! is.na(var_noNA) ]
            class_div <- classInt::classIntervals(var=var_noNA, n=ngroups[k], ...)
            class_div <- jiggleClass(x=class_div)
            breaks[[ k]] <- class_div$brks
            ngroups[ k ] <- length(breaks[[ k ]]) -1
            
            if ( usebrewer==TRUE ) { plotclr[[k]] <- RColorBrewer::brewer.pal(n=ngroups[k], name=paletteName[k]) }
            else {  plotclr[[k]] <- colorVec[[k]][ 1:ngroups[ k ] ]  }
                   
            
            }
        }#end of setting plot colors and breaks
    
      
    #set default titles
    if (is.null(titles)) {
  
      #default titles
      wt <- ""    
  
      if (wt.label == TRUE) {      
       if (wt.ind == FALSE) { wt <- "(unwtd.)" }
       else { wt <- "(wtd.)" }
      }
      
      titles <- expand.grid(paste(names(statmats[[1]]), "of") , paste(var.pretty, ifelse(num.col >1, "\n",""), sep=""))
      titles <- apply(titles, 1, paste, collapse=" ")
      titles <- paste(titles, "by", geo.pretty, wt, sep=" ")
  
    }     



  #  #add labels for geography
   if (map.label == TRUE) {
  
  sl1 <- list('panel.text', sp::coordinates(map.file)[,1], sp::coordinates(map.file)[,2], 
              labels=map.file@data[, map.label.names ], cex=cex.label, col=col.label, ...) 
  
  #this allows you to add extra polygons, etc. on top of   
  if( length(sp_layout.pars) >0 )  {
    
    sl2 <- sp_layout.pars
    sl2[[ length(sl2)+1 ]] <- sl1
    
  }
  
  else{ sl2 <- sl1 }
   }

  
  
  


  #unlist all the variables     
  unlisted <- lapply(sapply(statmats, "[",i=1:num_stats), "[", i=1:(length(plotvars)+1))


     #loop through output statistics 
     for (k in 1:length(unlisted)) {
       
      
       
       
           map_tmp <- map.file

           #merge each time
           #added map_tmp@data
           map_tmp@data <- merge(x=map_tmp@data, y=unlisted[[k]], 
                             by.x=map.geo.var, by.y=d.geo.var,
                             all.x=TRUE, sort=FALSE)
           map_tmp@data <- map_tmp@data[ order(map_tmp@data$sort.id), ]

      

           
        if (map.label==TRUE) {        
           
                 tmp_plot <- sp::spplot(obj=map_tmp, zcol=plotvars,names.attr=prettynames,
                                        col.regions=plotclr[[k]], at=breaks[[k]], 
                                        main=list(label=titles[k], cex=cex.title), sp.layout=sl2, ...)    
                      
              }

       #if no labels
       else {
         
          if( length(sp_layout.pars) >0 ) {

                 tmp_plot <- sp::spplot(obj=map_tmp, zcol=plotvars,names.attr=prettynames,
                                        col.regions=plotclr[[k]], at=breaks[[k]], 
                                        main=list(label=titles[k], cex=cex.title), sp.layout=sp_layout.pars, ...)
                      }
          else {
	 
                 tmp_plot <- sp::spplot(obj=map_tmp, zcol=plotvars, names.attr=prettynames ,
                                        col.regions=plotclr[[k]], at=breaks[[k]],  
                                        main=list(label=titles[k], cex=cex.title), ...)
                      }
 
    	}
     

        list_of_plots[[ titles[k] ]] <- tmp_plot
       

     

     }#finish looping over matrices

    
    #now rearrange if plotbyvar==FALSE, then plot by statistic

    if (plotbyvar==TRUE) { ord <- 1:length(unlisted) } 
    else {
      
       i <- c(1, num_stats+1)
       m <- rep(0:(num_stats-1), each=num_vars)
       ord <- i+m
       
    }
    


    list_of_plots[ ord ]



}

