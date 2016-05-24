NULL
#' Pushes SMET object file into a \code{\link{meteoioini-class}} object
#' 
#' @param ... \code{\link{smet-class}} objects
#' @param smetlist (alternative) list of \code{\link{smet-class}} objects
#' @param ini a \code{\link{meteoioini-class}} object  
#' @param smetwdir  optional working directory used for SMET writing on disks. Default is the package internal directory \code{system.file("temp",package="RSMET")}.
#' @param force.smetdir logical value. In case it is \code{TRUE} the SMETs are written on \code{smetwdir} directory, otherwise thay can be written on \code{slot(ini,"Input")$METEOPAH} directory if specified.
#' @param new.ini.file character string. DEfault is \code{NA}. If it is specified, the returned \code{\link{meteoioini-class}} object is saved as the \code{new.ini.file} file. 
#' @param append logical value. Default is code{FALSE}. If it is \code{TRUE} SMETs will be appended to the pre-existing ones. If it is \code{TRUE}, \code{force.smetdir} and \code{smetwdir} arguments will be ingored.
#' 
#' @return An updated \code{\link{meteoioini-class}} object.
#' 
#' @export 
#' 
#' @examples
#' 
#' ini <- as.meteoioini("test")
#' 
#' AA <- as.smet("test")
#' BB <- as.smet("test")
#' 
#' 
#' newini <- pushSmetIntoIni(AA,BB,ini=ini)
#' newini <- pushSmetIntoIni(AA=AA,BB=BB,ini=ini)
#' newini <- pushSmetIntoIni(AA=AA,BB=BB,ini=ini)
#' 
#' 
#' 
#' 
#' ##  NOT YET IMPLEMENTED
#' 


pushSmetIntoIni <- function(...,
		smetlist=NULL,
		ini=as.meteoioini("test"),
		smetwdir=system.file("temp",package="RSMET"),
		force.smetdir=TRUE,
		new.ini.file=NA,
		append=FALSE) {
			
			
	
			if (is.null(smetlist)) {
				
				
				smetlist <- base::list(...)
				
				if (is.null(names(smetlist))) {
					
					names(smetlist) <- sapply(X=smetlist,FUN=function(x){x@header$station_id})
				}
				
				if (is.null(ini)) {
					
					
					warning("Meteoio Ini object is not espliciated! It is taken implicitly!!")
					len <- length(smetlist)
					ini <- smetlist[[len]]
					smetlist <- smetlist[-len]
					
					
				}
				
				
			}
			
			out <- ini
			
			input <- out@Input
			
			#### SMET
			cond <-  !is.null(input$METEO)
			if (cond==FALSE) {
				
				message("No METEO key value specified, automatically set as SMET")
				input$METEO <- "SMET"
			} 
			if (!identical(input$METEO,"SMET")) {
				
				message("METEO key value forced to SMET!")
				input$SMET <- "SMET"
				
				
			}	
			##### SMETDIR
			
			if (append==TRUE) smetwdir <- input$METEOPATH
			
			cond <-  !is.null(input$METEOPATH) 
			if (cond==FALSE) {
				
				message("No METEOPATH key value specified, automatically set as smetwdir argument")
				input$METEOPATH <- smetwdir
			} 
			
			if (force.smetdir==TRUE) {
				
				input$METEOPATH <- smetwdir
				
				if (append==TRUE) {
						append=FALSE
					    warning("Forcing SMET directory, append switched to FALSE!")
					}
				
			}
			if (!identical(input$METEO,"SMET")) {
				
				message("METEO key value forced to smetwdir argument!")
				input$SMET <- smetwdir
				
				
			}	
			
			if (append==TRUE) {
				
					start <- which(stringr::str_detect(names(input),"STATION"))
					start <- length(start)
					
			} else {
				
					toremove <-  which(stringr::str_detect(names(input),"STATION"))
					input <- input[-toremove]
					start <- 0				
				
			}
			
			newstation_keys <- paste("STATION",start+1:length(smetlist),sep="")
			
			if (is.null(names(smetlist))) names(smetlist) <- newstation_keys
			
			names(newstation_keys) <- names(smetlist)
			
			for (it in names(smetlist)) {
				
				file <- paste(smetwdir,it,sep="/")
				file <- paste(file,".smet",sep="")
				input[[newstation_keys[it]]] <- it
				smetlist[[it]]@file <- file			
				
			}
			
			
			lapply(X=smetlist,FUN=print,file="internal")
			
			
			out <- RSMET::as.meteoioini(out,Input=input,file=as.character(new.ini.file))
			
			
			
			
			
			
			return(out)
			
		}








