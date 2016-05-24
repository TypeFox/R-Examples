


scpcp <- function(data, freqvar = "Freq", max.N = 1e6,
                gap = 0.2,
                sort.individual=TRUE,
                level.width=0.2,
                polygon = TRUE,
                base.colour = alpha("black",0.7),
                label = TRUE, lab.opt = list(rot = 0, col = 1, bg = TRUE, abbr = FALSE, abbr.var = 12, hide.sel = TRUE, var.labels = TRUE),
                sel =  NULL, 
                sel.hide = TRUE,
                sel.palette = NULL,
                col.opt = list(),
                plot = TRUE,
	            return.coords = !plot)
{

data <- as.data.frame(data)

 #if(inherits(data,"table")){
 	# convert to data.frame
 #	data <- subtable(data,1:length(dim(data)),allfactor=TRUE)
  #  freqvar <-  'Freq'
 #}



   ##################### FORMAT #####################

  # Expansion, if freqvar specified
  if(!is.null(freqvar)){
  	if(freqvar %in% names(data)){
  		if( (N <- sum(data[freqvar])) > max.N){
  			data[freqvar] <- round(data[freqvar]/N*max.N, 0) 
  		}
  		data <- untableSet(data, freqvar = freqvar)
  	}
  }
for(i in 1:ncol(data)){
	data[,i] <- as.factor(data[,i])
}


 orig.nm <- names(data)
 
  

  
  
################# PARAM INITIALIZATIONS #####################
if( "var.labels" %in% names(lab.opt) ){
	var.labels <- lab.opt$var.labels
}else{
	var.labels <- TRUE
}

if( "rot" %in% names(lab.opt) ){
	text.rotation <- lab.opt$rot
}else{
	text.rotation <- 0
}
if( "las" %in% names(lab.opt) ){
	las <- lab.opt$las
}else{
	las <- 1
}

if( "lab.cex" %in% names(lab.opt) ){
	lab.cex <- lab.opt$lab.cex
}else{
	lab.cex <- 1
}
if( "cex.axis" %in% names(lab.opt) ){
	cex.axis <- lab.opt$cex.axis
}else{
	cex.axis <- 1
}
if( "hide.sel" %in% names(lab.opt) ){
	hide.sel <- as.logical(lab.opt$hide.sel)
}else{
	hide.sel <- TRUE
}

if( "abbr.var" %in% names(lab.opt) ){
	abbr.var <- lab.opt$abbr.var
}else{
	abbr.var <- 12
}

if( "bg" %in% names(lab.opt) ){
	text.background <- lab.opt$bg
}else{
	text.background <- TRUE
}
if( "abbr" %in% names(lab.opt) ){
	text.abbreviation <- lab.opt$abbr
}else{
	text.abbreviation <- FALSE
}
if( "col" %in% names(lab.opt) ){
	text.colour <- lab.opt$col
}else{
	text.colour  <- "black"
}

if( !("alpha" %in% names(col.opt)) ){
	if(polygon){
		col.opt$alpha <- min(0.7,log(1000)/(log(1000)+log(nrow(data)/10)))
	}else{
		col.opt$alpha <- min(0.7,1000/nrow(data))
	}
}

  
  
  
   ##################### SELECTION for COLORBRUSH #####################
  

  # doodle kann entweder durch einen externen Vektor (sel), oder durch direkte ?bergabe der zu highlightenden F?lle im Funktionsaufruf initiiert werden
  # zb bei sel="Cont=='Low'&Infl=='Low'" 
  
  if(!is.null(sel)){
    doodle <- TRUE
    if(is.character(sel)){
      selection <- with(data, as.factor(eval(parse(text=sel))) ) 
    }else{ 
      selection <- as.factor(sel)
    }
    ndv <- length(levels(selection)) 
    
    #[ord<-order(selection,decreasing = TRUE),]
    #if(!is.null(numeric)){
    #	orig.data <- orig.data[ord,]
    #	num.data <- num.data[ord,]
    #}
  } else {
    doodle <- FALSE
    ndv <- 1
  }
  if(doodle){
  		data <- cbind(selection,data)
  }

# if !sel.hide then the selection variable is accepted as a new additional column of data
# otherwise it is removed later on, after itw as used for the reorderings/polygons by color

  sel.hide <- doodle & sel.hide

############## PARAM 2 #############

	# these are using only the factor variables!
  N <- nrow(data)
  m <- ncol(data)
  nm <- names(data)
  labels <- lapply(data, levels)

  
  if(text.abbreviation != FALSE){
  	labels.abbr <- lapply(labels, abbreviate, minlength=text.abbreviation)
  }else{
  	labels.abbr <- labels
  }
 


 ############## FUNDAMENTAL REORDERING ##################
  
 # reordering the dataset (including selection variable):
 base.order <- do.call( order, c(data, decreasing = FALSE))
 data <- data[base.order,,drop = FALSE] 



  ##################### COMPUTATIONS  #####################
  
  # point sequences
  
  sequences <- lapply(data, table)
  

  ##### line coordinates

  seq.list <- lapply(sequences, function(z){ # Berechnung aller Koordinaten f?r alle Linien (nicht nicht zusammenh?ngend)
    p <- c(0,cumsum(z/N))*(1-gap) # Proportionen der Kategorien in nicht-gap-Bereich
    k <- length(z) # Anzahl Kategorien
    gap.proportion <- gap/(k-1) # H?he jedes gap-Bereichs
    seqs <- list()
    for(i in 1:k){
      tmp <- seq(p[i],p[i+1], (p[i+1]-p[i])/max(1,(z[i]-1))) + (i-1)*gap.proportion 
      seqs[[i]] <- tmp[1:z[i]]
    }
    return(seqs)
  })
  

  ##### polyline ids

  id.list <- mapply( function(y, z){ # Durchnummerierung
    lapply(y, function(w){
      which(z == w)
    })
  },y = labels, z = data, SIMPLIFY=FALSE)


  ##### line coordinates

  lines.unsorted <- mapply(function(y,z){ # Kombinieren von seq.list und id.list um komplette Linien zu erzeugen
    ret <- rep(0,N) 
    for(i in seq_along(y)){                
      ret[ z[[i]] ] <-   y[[i]]
    }
    return(ret)  
  }, y = seq.list, z = id.list )
  


  ##################### additional recursive sort  ##################### 
  
  if(sort.individual & !doodle){  # Sortiert jede Sequenz nach den R?ngen der jeweils linken Variable
    lines <- lines.unsorted
    e1 <- environment()
    for(i in 2:m){
    		# the left variable is a factor (included in lines)
    		sapply(id.list[[i]], function(z){
      			rr <- rank(lines[z,i-1])
        		e1$lines[z,i]  <- lines[z,i][rr] 
        		return(invisible(TRUE))
      		})      
    }  
  } else if(sort.individual & doodle){ #Sortiert nach der jeweils linken Variable UND dem doodle
    lines <- lines.unsorted
    e1 <- environment()
    for(i in 2:m){
    	# the left variable is a factor (included in lines)
    		sapply(id.list[[i]], function(z){
        		tapply(z,e1$data[z,1],function(y){
        			e1$lines[y,i]  <- lines[y,i][rank(e1$lines[y,i-1])]
        		})
        		return(invisible(TRUE))
        	})
    }
  }  else {
    lines<-lines.unsorted
  }
    
  
  ##################### PREPARATIONS ###################### 
  
    
  # separate selection variable

  if(sel.hide){ 
  	# the selection variable is removed after computations
    selection<-data[,1]
    data<-data[2:m]
    lines<-lines[,2:m]
    lines.unsorted <- lines.unsorted[,2:m]
    seq.list <- seq.list[2:m]
    id.list <- id.list[2:m]
    labels <- labels[2:m]
    labels.abbr <- labels.abbr[2:m]
    m <- ncol(data)
    nm <- nm[2:length(nm)]
    
    
  } else if(doodle){
    selection <- data[,1] 
  } 

  

    

  ##### additional coordinates

  middles<- as.vector(rapply(seq.list, mean)) 
if(plot){
	
	  # set margins
  if(is.character(sel)){
    par(mar=c(3,0,0,0)) 
  } else {
    par(mar=c(2,0,0,0)) 
  }  
  
  ##### color palette
  if(doodle){
  	if(is.null(sel.palette)){
  		#ndv <- length(levels(selection))
  		if(ndv < 9){
  			sel.palette <- 1:ndv
  		}else{
  			sel.palette <- "rgb"
  		}
	}
  	
      sel.colours <- getcolors(ndv, sel.palette ,col.opt = col.opt) 

  } else {
    sel.colours <- getcolors(1, base.colour ,col.opt = col.opt)
  }
  
  	if( "border" %in% names(col.opt) ){
		rect.border <- rep(col.opt$border,length(sel.colours))[1:length(sel.colours)]
	}else{
		rect.border  <- sel.colours
	}
  
  
  ##################### LINE-PLOT ###################### 
  
  if(!polygon){
    	xcoords <- c(1,1+level.width)
    	for (i in 2:m){
      		xcoords <- c(xcoords, c(i,i+level.width))
    	}
    	lines.doubled <- lines[,rep(1:m,each=2)]     
    
   # coordinate computation given level.width
   
    
   ## draw colored lines hrouped by highlighting
      plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA,panel.first ={
         	for(j in 1:ndv){# ndv = length(levels(selection))
         		if(ndv > 1){ # length(levels(selection))
         			apply(lines.doubled[which(selection==levels(selection)[j]),],1,function(z){
      					lines(xcoords,z,col=sel.colours[j])
      				})
         		}else{
         			apply(lines.doubled,1,function(z){
      					lines(xcoords,z,col=sel.colours)
      				})
         		}
      		}
         } , type="l") 
  } else {
    
    ##################### POLYGON-PLOT ######################
    
    if(polygon){     
      
      
          
      ##### POLYGONS
    
     as.pol<- function(color = 1){
        cc<-as.factor(paste(data[,1])) # polygons for first axis
        M2 <- cbind(tapply(lines[,1],cc,min),tapply(lines[,1],cc,max),tapply(lines[,1],cc,max),tapply(lines[,1],cc,min)) 
        apply(M2,1,function(z){
          polygon(x= c(1,1,1+level.width,1+level.width), y = z, col = sel.colours,border=rect.border) 
          return(invisible(TRUE))
        })

        for(i in 2:m){
          cc<-as.factor(paste(cc,data[,i])) # achsenverbindenden Polygone
          M <- cbind(
          		tapply(lines[,i-1],cc,min),
          		tapply(lines[,i-1],cc,max),
          		tapply(lines[,i],cc,max),
          		tapply(lines[,i],cc,min))
          		
          apply(M,1,function(z){
            polygon(x= c(i-1+level.width,i-1+level.width,i,i), y = z, col = sel.colours,border=rect.border)  
            return(invisible(TRUE))
          })
          
          M1<-M[,c(3,4,4,3), drop = FALSE] # polygons for all other axes
          apply(M1,1,function(z){
            polygon(x= c(i,i,i+level.width,i+level.width), y = z, col = sel.colours ,border=rect.border)
            return(invisible(TRUE))
          })
        }       
      }

      doodling <- function(){
        
        ##### POLYGONE        

        for(j in 1:ndv){
          cc<-as.factor(paste(data[which(selection==levels(selection)[j]),1]))
          M2 <- cbind(tapply(lines[which(selection==levels(selection)[j]),1],cc,min),
          tapply(lines[which(selection==levels(selection)[j]),1],cc,max),
          tapply(lines[which(selection==levels(selection)[j]),1],cc,max),
          tapply(lines[which(selection==levels(selection)[j]),1],cc,min))
          apply(M2,1,function(z){
            polygon(x= c(1,1,1+level.width,1+level.width), y = z, col = sel.colours[j],border=rect.border[j])
            return(invisible(TRUE))
          })
          for(i in 2:m){
            cc<-as.factor(paste(cc,data[which(selection==levels(selection)[j]),i]))
            M <- cbind(tapply(lines[which(selection==levels(selection)[j]),i-1],cc,min),
            tapply(lines[which(selection==levels(selection)[j]),i-1],cc,max),
            tapply(lines[which(selection==levels(selection)[j]),i],cc,max),
            tapply(lines[which(selection==levels(selection)[j]),i],cc,min))
            
            apply(M,1,function(z){
              polygon(x= c(i-1+level.width,i-1+level.width,i,i), y = z, col = sel.colours[j],border=rect.border[j])  
              return(invisible(TRUE))
            })
            
            M1<-M[,c(3,4,4,3), drop = FALSE]
            apply(M1,1,function(z){
              polygon(x= c(i,i,i+level.width,i+level.width), y = z, col = sel.colours[j],border=rect.border[j])
              return(invisible(TRUE))
            })
          }
        }
      }
      
      
##### draw
      
      if(!doodle){
        plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA,panel.first = as.pol(), type="l")
      } else {
        plot(1, xlim=c(1,m+level.width),ylim=c(0,1), axes=FALSE, xlab=NA, ylab=NA,panel.first = doodling(), type="l")
      }
    }
  }  
  
  ##################### LEVELS & LABELS ##################### 
if(label){
  middlesX <- list()
  for(i in 1:ncol(lines.unsorted)){
    middlesX[[i]] <-  tapply(lines.unsorted[,i],data[,i],mean)
  }
  
  if (text.abbreviation){
    ll <- labels.abbr
  } else {
    ll <- labels
  }
 
xx <- 0.5*level.width + seq_along(middlesX)


xx <- rep(xx,sapply(middlesX,length))
yy <- unlist(middlesX)
zz <- unlist(ll)


   if(text.background){
        bgtext(xx, yy, zz, col="gray90", font=1, bg=text.colour, srt=text.rotation, cex = lab.cex)
     } else {
        text(xx, yy, zz, col=1, font=2, srt = text.rotation, cex = lab.cex)
        text(xx, yy, zz, col="gray90", font=1, srt = text.rotation, cex = lab.cex)
     } 
}   
  ##### PLOT-UNTERSCHRIFTEN
  
if(var.labels){
  if(abbr.var != FALSE){
  	nm <- sapply(nm, abbreviate, minlength=as.integer(abbr.var))
  }
  mtext(nm, side=1, line=0, at=(1:m)+0.5*level.width, font=2, las = las, cex = cex.axis) # Variablennamen unter die jeweiligen Achsen
  if(is.character(sel) & !hide.sel){ # falls vorhanden, Beschreibung des Highlightings
    mtext(paste("Highlight:",sel), side=1, line=1, at=((m+1+level.width)/2), font=1, cex=0.8, col= "red")  
  }  
}
}#end if plot

  ##################### RETURN #####################  
  
  if(return.coords){
  	colnames(lines) <- paste("y",colnames(lines),sep=".")
  	ret <- as.data.frame(cbind(lines,data))
  	return(invisible(ret))
  }else{
  	return(invisible(TRUE))
  }
}

bgtext <- function(x, y, labels, col='white', bg='black',
	k = 8, r=0.1, ... ) {

	theta= seq(pi/4, 2*pi, length.out=k)
	#xy <- xy.coords(x,y)
	xr <- r*strwidth('A')
	yr <- r*strheight('A')

	text( rep(x,each = k+1) + c(cos(theta)*xr,0), rep(y,each = k+1) + c(sin(theta)*yr,0), 
		rep(labels,each=k+1), col= c(rep(bg,k),col), ... )

	#text(x, y, labels, col=col, ... ) 
	return(invisible(TRUE))
}

