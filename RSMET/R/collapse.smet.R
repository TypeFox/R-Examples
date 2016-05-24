NULL

#'    Coerces a \code{smet-class} object to a data frame
#' 
#' 
#' 
#' @param x first \code{\link{smet-class}} object 
#' @param y second \code{\link{smet-class}} object 
#' @param headers header attributes to  check before collapsing 
#' @param date.field   field name used for date and time. Default is \code{"timestamp"}, as used for \code{SMET} format. 
#' @export
#'
#' @return a collapsed \code{\link{smet-class}} object
#' @examples 
#' 
#' smet1 <- as.smet(system.file("examples/zugspitze1.smet",package="RSMET"))
#' smet2 <- as.smet(system.file("examples/zugspitze2.smet",package="RSMET"))
#' 
#' newsmet <- collapse.smet(x=smet1,y=smet2)
#' 
#' newsmet_ <- collapse.smet(x=smet2,y=smet1)
#' 
#' 
#' x <- smet(system.file('examples/PIEM001114.smet',package="RSMET"))
#' y <- smet(system.file('examples/PIEM001114_20160308.smet',package="RSMET"))
#' 
#' 
#' xy <-  collapse.smet(x,y) 
#' 

collapse.smet <- function (x,y,headers=NULL,date.field="timestamp") {
	
	x <- as.smet(x) ## x and y must be smet!
	y <- as.smet(y)
	
	if (is.null(headers)) {
		
		headers <- names(x@header)
		
	}
	
	headers <- headers[headers %in% intersect(names(x@header),names(y@header))]
		
	header_values_x <- x@header[headers]
	header_values_y <- y@header[headers]
		
		
	#######
	
	
	
	xf <- fields(x)
	yf <- fields(y)
	
	ff <- union(xf,yf)
	
	multiplier <- union(x@header$units_multiplier,y@header$units_multiplier)
	offset <- union(x@header$units_offset,x@header$units_offset)
	
	cond3 <- all(names(multiplier)==names(offset)) & all(names(multiplier)==ff)
	
	if (cond3!=TRUE) {
		
		
		msg <- sprintf("Mismatch in fiels: %s , %s , %s",paste(names(multiplier),collapse=" "),
				paste(names(offset),collapse=" "),paste(ff,collapse=" "))
								
		stop(msg)
		
	}
	
	dy <- ff[!(ff %in% xf)] 
	dx <- ff[!(ff %in% yf)]
	
	
	x[,dy,multiplier=multiplier[dy],offset=offset[dy]] <- NA ## look at 

	y[,dx,multiplier=multiplier[dx],offset=offset[dx]] <- NA 
	
	x <- x[,ff]
	y <- y[,ff]
	
	
	if (!identical(header_values_x,header_values_y)) {
		
		s <- list()
		mx <- list()
		my <- list()
		cwarn <- list()
		for (it in headers) {
			mx[[it]] <- str_trim(paste(header_values_x[[it]],collapse=" "))
			my[[it]] <- str_trim(paste(header_values_y[[it]],collapse=" "))
			s[[it]] <- sprintf("%s: %s == %s",it,mx[[it]],my[[it]]) ###%(header_values_x)
			cwarn[[it]] <- (mx[[it]]!=my[[it]])
			
		}
		
		iw <- which(unlist(cwarn))
		if (length(iw)>=1) { 
			iw <- unique(c(iw,which(headers=="station_id")))
			m <- c("Collapsing SMET : header mismatch!!!",unlist(s[iw]))
			m <- paste(m,collapse="   ")
			
			warning(m)
		}
	}
	###################
	
	
	
	
	
	####################
	
	
	
	out <- x 
	

	xcd <- which(names(x@data)==date.field)
	ycd <- which(names(y@data)==date.field)
	
	x@data[,xcd] <- as.POSIXct(x@data[,xcd])
	y@data[,ycd] <- as.POSIXct(y@data[,ycd])
	## 
	##
	dd <- rbind(x@data,y@data)
	
	rownames(dd) <- NULL
	
	dc <- which(names(dd)==date.field)
	
	dd <- dd[order(dd[,date.field]),]
	uniq <- which(!duplicated(dd[,date.field]))
	
	nonuniq <- which(duplicated(dd[,date.field]))
	
	nonuniq_prev <- nonuniq-1
	
	cond <- array(FALSE,length(nonuniq))
	
	if (length(nonuniq)>=1) for (i in 1:length(nonuniq)) {
		###print(i)
		prev <- as.vector(dd[nonuniq_prev[i],])
		v <- as.vector(dd[nonuniq[i],])
	
		cond[i] <- all(prev==v)
		
		
	}
	
	
	
	toprint <- which((1:nrow(dd)) %in% c(uniq,nonuniq[cond==FALSE]))

	
	dd <- dd[toprint,]
	
	out@data <- dd
	
	out@file <- as.character(NA)
	#######
	
	
	return(out)
	
	
}