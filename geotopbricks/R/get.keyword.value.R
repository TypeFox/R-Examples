# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL
#' Importing a GEOtop Keyword and its Value into R
#' 
#' It returns the values of a keyword of "geotop.inpts" file or data frame with the suitable format.
#' 
#' @param keyword keyword name
#' @param inpts.frame data frame returned by \code{\link{declared.geotop.inpts.keywords}} or \code{NULL}. Default is \code{NULL}.
#' @param vector_sep character value for the separator chacter if Keyword Value must be returned as a vector, otherwise it is \code{NULL}. Default is \code{NULL}, but if \code{numeric} or \code{date} are \code{FALSE},  \code{vector_sep} is set \code{","} by default.
#' @param numeric logical value. If \code{TRUE} the Value has numeric type, otherwise it is a string or string vector. Default is \code{FALSE}.
#' @param date logical value. If \code{TRUE} the Value is retured as \code{\link{POSIXlt}} date, otherwise it is a string or string vector. Default is \code{FALSE}. 
#' @param format string format representing the date, see \code{\link{as.POSIXlt}}, used if \code{date} is \code{TRUE}. Default is \code{"\%d/\%m/\%Y \%H:\%M"} (which is the format used in \code{geotop.inpts} keyword \code{InitDateDDMMYYYYhhmm})
#' @param tz format string representing the time zone, see \code{\link{as.POSIXlt}}, used if \code{date} is \code{TRUE}. Default is \code{"Etc/GMT+1"} (until the previous version it was \code{"A"}) which meens UTC +1.
#' @param raster logical value. Default is \code{FALSE}. If \code{TRUE} function returns direclty the raster map as \code{\link{Raster-class}} object built with \code{\link{raster}} method. 
#' @param file_extension Extension to be added to the keyword if keyword is a file name. Default is \code{".asc"}
#' @param wpath working directory containing GEOtop files (included the inpts file). It is mandatory if \code{raster} is \code{TRUE}. See \code{\link{declared.geotop.inpts.keywords}}.
#' @param add_wpath logical value. Default is \code{FALSE}. If \code{TRUE}, the \code{wpath} string is attached to the keyword string value. It is automatically set \code{TRUE} if \code{raster} is \code{TRUE}.
#' @param use.read.raster.from.url logical value. Default is \code{TRUE}. If \code{TRUE} the RasterLayer are read with \code{\link{read.raster.from.url}}, istead of \code{\link{raster}} (otherwise). It is recomended in case the files whose paths are contained in \code{x} are remote and are 'http' addresses. In this cases the stand-alone method \code{raster(x)} does not always work and \code{use.read.raster.from.url} is necessary.  
#' @param data.frame logical value. It is an option for tabular data. If \code{TRUE} function returns direclty a data frame  or a list of  data frames as \code{\link{data.frame}} or \code{\link{zoo}} objects imported from the keyword-related files  using \code{\link{read.table}} function. In this case the argument \code{wpath} (see \code{\link{declared.geotop.inpts.keywords}}) is mandatory. Default is \code{FALSE}.
#' @param formatter string value. It is the decimal formatter contained in the file name and used in case the tabular data are referred at several points. Default is \code{"\%04d"} . It is used in case \code{data.frame} is \code{TRUE}. 
#' @param level integer values. Numbers incating all the identandification numbers of the files containing the requested data frames. Default is 1, correspondig to the decimal formatter \code{"0001"}. See examples. 
#' @param date_field string value. Default is "Date", otherwise defined by the value of \code{HeaderDateDDMMYYYYhhmmMeteo} geotop keyword. It is used only if the argument \code{data.frame} is \code{TRUE}. If it is \code{NULL} or \code{NA} the function return a list of generic \code{\link{data.frame}} object(s), otherwise \code{link{zoo}} object(s). See the arguments \code{tz} and \code{format} for Date formatting.
#' @param isNA numeric value indicating NA in geotop ascii files. Default is -9999.00
#' @param matlab.syntax logical value. Default is \code{FALSE}. If \code{TRUE} a vector is written in a string according to *.m file syntax. Warning: this synstax is not read by GEOtop. 
#' @param projfile fileneme of the GEOtop projection file. Default is \code{geotop.proj}.
#' @param start_date,end_date null objects or dates in \code{POSIXlt} format between which the variables are returned. It is enabled in case that \code{date_field} is not \code{NULL} or \code{NA} and \code{data.frame} is \code{TRUE}. Default is \code{NULL}.
#' @param zlayer.formatter decimal formatter. It is used if \code{data.frame==TRUE} and the columns refers to different soil depths. Default is \code{NULL}. 
#' @param z_unit z coordinate measurement unit. GEOtop values expressed in millimeters which are converted to centimeters by default. Default is \code{c("centimeters","millimeters")}. Otherwise can be the ratio between the unit and one meter. It is used if \code{zlayer.formatter=="z\%04d"} or similar.
#' @param geotop_z_unit z coordinate measurement unit used by GEOtop. Default is \code{millimeters}. It is used if \code{zlayer.formatter=="z\%04d"} or similar.
#' 
#' @param ContinuousRecovery integer value. Default is 0. It is used for tabular output data and is the number of times GEOtop simulation broke  during its running and was re-launched with 'Contiuous Recovery' option. 
#' @param ContinuousRecoveryFormatter character string. Default is \code{'_crec\%04d'}. It is used only for tabular output data and if \code{ContinuousRecovery} is equal or greater than 1. 
#' @param ... further arguments of \code{\link{declared.geotop.inpts.keywords}} 
#' 
#' @export 
#' 
#' @note If \code{inpts.frame} is \code{NULL}, \code{inpts.frame} will be obtained by calling the function \code{\link{declared.geotop.inpts.keywords}} with \code{...} arguments.
#' @return the keyword value 
#' @import stringr 
#' @import zoo
#' @examples
#' 
# library(stringr)  
#' library(geotopbricks)
#' 
#' #Simulation working path
#### wpath <- 'http://www.rendena100.eu/public/geotopbricks/simulations/panola13_run2xC_test3'
#' wpath <- 'http://www.rendena100.eu/public/geotopbricks/simulations/panola13_run2xC_test3'
 ###  wpath <- 'http://meteogis.fmach.it/idroclima//panola13_run2xC_test3'
#' prefix <- get.geotop.inpts.keyword.value("SoilLiqWaterPressTensorFile",wpath=wpath)
#' 
#' slope <- get.geotop.inpts.keyword.value("SlopeMapFile",raster=TRUE,wpath=wpath) 
#' bedrock_depth <- get.geotop.inpts.keyword.value("BedrockDepthMapFile",raster=TRUE,wpath=wpath) 
#' 
#' layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath)
#' names(layers) <- paste("L",1:length(layers),sep="")
#' 
#' ##### set van genuchten parameters to estimate water volume 
#' theta_sat <- get.geotop.inpts.keyword.value("ThetaSat",numeric=TRUE,wpath=wpath)
#' theta_res <- get.geotop.inpts.keyword.value("ThetaRes",numeric=TRUE,wpath=wpath)
#' alphaVG <-  get.geotop.inpts.keyword.value("AlphaVanGenuchten",
#' numeric=TRUE,wpath=wpath) # expressed in mm^-1
#' 
#' nVG <-  get.geotop.inpts.keyword.value("NVanGenuchten",numeric=TRUE,wpath=wpath) 
#' 
#' 
#' ##### end set van genuchten parameters to estimate water volume
#' 
#' 
#' ##### set meteo data
#' 
#' start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A") 
#' end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A") 
#' 
#' nmeteo <- get.geotop.inpts.keyword.value("NumberOfMeteoStations",numeric=TRUE,wpath=wpath)
#' level <- 1:nmeteo
#' 
#' # Uncomment the following lises to run the R code: 
#' 
#' ## set meteo data
#' 
#' 
#'  \dontrun{
#' # meteo <- get.geotop.inpts.keyword.value("MeteoFile",wpath=wpath,data.frame=TRUE,
#' #       level=level,start_date=start,end_date=end)
#' }
#' 
#' ##### end set meteo data
#' 
#' ## IMPORTING AN OUTPUT SOIL MOISTURE PROFILE: 
#' 
#'  wpath <- 'http://www.rendena100.eu/public/geotopbricks/simulations/Muntatschini_pnt_1_225_B2_004'
#' 
#' \dontrun{
#' #	SMC  <- get.geotop.inpts.keyword.value("SoilLiqContentProfileFile",
#' #          wpath=wpath,data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",
#' #          formatter="%04d")
#' #
#' #    SMCz  <- get.geotop.inpts.keyword.value("SoilLiqContentProfileFile",
#' #         wpath=wpath,data.frame=TRUE,date_field="Date12.DDMMYYYYhhmm.",
#' #          formatter="%04d",zlayer.formatter="z%04d")
#' }
#' 
#' 
#' 		
#' 
#' 
#' 
#' 
get.geotop.inpts.keyword.value <- function(keyword,inpts.frame=NULL,vector_sep=NULL,numeric=FALSE,format="%d/%m/%Y %H:%M",date=FALSE,tz="Etc/GMT+1",raster=FALSE,file_extension=".asc",add_wpath=FALSE,wpath=NULL,use.read.raster.from.url=TRUE,data.frame=FALSE,formatter="%04d",level=1,date_field="Date",isNA=-9999.000000,matlab.syntax=TRUE,projfile="geotop.proj",start_date=NULL,end_date=NULL,ContinuousRecovery=0,ContinuousRecoveryFormatter="_crec%04d",zlayer.formatter=NULL,z_unit=c("centimeters","millimeters"),geotop_z_unit="millimeters",...) {
#####	check.columns=FALSE
# Added by the author on Feb 6 2012	
	
	if (length(keyword)>1) {
		out <- NULL
		
		out <- lapply(X=keyword,FUN=get.geotop.inpts.keyword.value,inpts.frame=inpts.frame,vector_sep=vector_sep,numeric=numeric,format=format,date=date,tz=tz,raster=raster,file_extension=file_extension,add_wpath=add_wpath,wpath=wpath,use.read.raster.from.url=use.read.raster.from.url,data.frame=data.frame,formatter=formatter,level=level,date_field=date_field,isNA=isNA,matlab.syntax=matlab.syntax,projfile=projfile,...)
		names(out) <- keyword
		
		
		return(out)
	}
	
	
	
	if (is.null(inpts.frame)) inpts.frame <- declared.geotop.inpts.keywords(wpath=wpath,...)

	out <- inpts.frame$Value[inpts.frame$Keyword==keyword]

	if (length(out)==0) {
		print("Warning (get.geotop.inpts.keyword.value): keyword withot value:")
		print(keyword)
		return(NULL)
		
	}
	len <- str_length(out)
	
    
	
	if (len>0) {
		
		if ((str_sub(out,1,1)=='\"') |  (str_sub(out,1,1)=='\''))  out <- str_sub(out,2)
		len <- str_length(out)
		if ((str_sub(out,len,len)=='\"') |  (str_sub(out,len,len)=='\''))  out <- str_sub(out,end=len-1)
	}
	if ((numeric | date) & (is.null(vector_sep))) vector_sep <- "," 
	
	if (!is.null(vector_sep)) {
		if (numeric | matlab.syntax) {
			
			if ((str_sub(out,1,1)=='[') |  (str_sub(out,1,1)=='('))  out <- str_sub(out,2)
			len <- str_length(out)
			if ((str_sub(out,len,len)==']') |  (str_sub(out,len,len)==')'))  out <- str_sub(out,end=len-1)
		}
		
		out <- (str_split(out,vector_sep))[[1]]
		
		if (matlab.syntax) { 
			out <- str_replace_all(out,"\'","")
			out <- str_replace_all(out,"\"","")
		}
	}
	
	if (date) {
		
		out <- as.POSIXlt(out,format=format,tz=tz)
		
	} else if (numeric) {
		out <- as.numeric(out)
	} else if (raster) {
		add_wpath=TRUE
		
		if (!is.null(wpath)) out <- paste(wpath,out,sep="/")

		if (str_sub(file_extension,1,1)==".")  {
			filepath <- paste(out,file_extension,sep="") 
		} else { 	
			filepath <- paste(out,file_extension,sep=".") 
		}
		 if (use.read.raster.from.url) {
			 out <- read.raster.from.url(x=filepath)
		 } else {
		     out <- raster(x=filepath)
		}
		
		if (!is.null(wpath)) projfile <- paste(wpath,projfile,sep="/")
		cond <- file.exists(projfile)
		projection(out) <- getProjection(projfile,cond=cond)
	#	if (cond) {
	#		
	#		projection(out) <- readLines(projfile,warn=FALSE)
			
	#	}
		## ADD projection 
		
		
		
		
	} else if (data.frame) {
		
		if (file_extension==".asc" | file_extension=="asc") file_extension=".txt"

		
		
		keyword <- out
		out <- paste(wpath,out,sep="/")
		
		 if (is.null(formatter) | is.null(level) | length(level)<1) {
		 
			 formatter <- ""
		 } else {
			
			formatter <- array(formatter,length(level)) 
			for (i in 1:length(level)) {
				
				formatter[i] <- sprintf(formatter[i],level[i])
				
			} 
			
			
		 }	  
	  	
		 out <- paste(out,formatter,sep="")
		 
		 if (str_sub(file_extension,1,1)==".")  {
			 filepath <- paste(out,file_extension,sep="")
			 filecrec_extension=paste(ContinuousRecoveryFormatter,file_extension,sep="") ## Continous Recovery Option 
			 filecrecpath <- paste(out,filecrec_extension,sep="")
		 } else { 	
			 filepath <- paste(out,file_extension,sep=".") 
			 filecrec_extension=paste(ContinuousRecoveryFormatter,file_extension,sep=".") ## Continous Recovery Option 
			 filecrecpath <- paste(out,filecrec_extension,sep=".")
		 }
		 
		 ### 
		# ContinuousRecoveryMax <- 50
		 
		# ContinuousRecoveryFiles <- spritf(filecrecpath,1:ContinuousRecoveryMax)
		 
		 
		 
		 ContinuousRecoveryCond <- !is.na(ContinuousRecovery) & !is.null(ContinuousRecovery) & length(ContinuousRecovery)==1 & round(ContinuousRecovery)==ContinuousRecovery &  ContinuousRecovery>0 ## ec on 20143004 This condition was ContinuousRecovery>1
		 
		 if (ContinuousRecoveryCond) {
			 
			 
			 filecrecpath <- unlist(lapply(X=filecrecpath,FUN=function(x,nn) {sprintf(x,1:nn)},nn=ContinuousRecovery))
			 length_points <- length(filepath)
			 names_points <- filepath
			 
			exists <- file.exists(filecrecpath)
			filecrecpath <- filecrecpath[exists] 
			 
			 
			filepath <- c(filepath,filecrecpath)
			 
			 
		 }
		 out <- filepath
		 out <-  list()
		 
		 for (i in 1:length(filepath)) {
			 
			 if (is.null(date_field)) date_field <- NA 
			 # ADD POSSIBLE SSH CONNECTION!!! 
			# ec date 09-03-2013
			x <- filepath[i]
			
			### ec 
			
		
			
			
			
			####
#			if (str_sub(x,1,3)=='ssh' | str_sub(x,1,5)=='plink') {
#			
#				file <- pipe(x) # added line according to http://stackoverflow.com/posts/2226880/edit
#				open <- TRUE
#				
#			 }	else {
#				 
#				 file <- x 
#			 }
			
		##	 temp <- read.table(file,header=TRUE,sep=",")
			 
#			 if (check.columns==TRUE) {
#				 
#				 temp <- readLines(x)
#				 tempfolder <- system.file("temporary",package="geotopbricks")
#				 tempfile <- paste(tempfolder,"temp.csv")
#				# print(x)
#				# str(temp)
#		         templist <- str_split(temp,",")
#				 len <- length(templist[[1]])
#				 index <- which(unlist(lapply(X=templist,FUN=function(x,l) {length(x)==l},l=len)))
#				 
#				 
#		
#				 writeLines(text=temp[index],con=tempfile)
#				 file <- file(tempfile)
#		##		 stop("check01")
#				temp <- read.table(file,header=TRUE,sep=",")
#				 
#				 
#			 } else {
#				 file <- file(x)
#				 temp <- read.table(file,header=TRUE,sep=",")
#				 
#			 }
			
			 file <- file(x)
			 temp <- read.table(file,header=TRUE,sep=",")
			 
			 i_index <- which(names(temp)==date_field)
			 if (length(i_index>1)) {
				 
				 if (is.numeric(isNA) & length(isNA)==1) temp[,-i_index][temp[,-i_index]<=isNA] <- NA # added on 6 dec 2012
			 
			
			 }	
		
		#####
		
		
		
			 if (!is.null(date_field) & !is.na(date_field) & length(i_index)==1 & length(date_field)>0 & (!ContinuousRecoveryCond)) {
				
				 index <- temp[,i_index]
				 temp<- temp[,-i_index]
				 index <- as.POSIXlt(index,format=format,tz=tz)
				 temp <- as.zoo(temp)
				 index(temp) <- index
				 # insert sart date & date
				 if (!is.null(start_date) & !is.null(end_date)) { 
					 
					 #	print(index(temp)>=start_date & index(temp)<=end_date)
					 temp <- temp[index(temp)>=start_date & index(temp)<=end_date,]
					 
				 }
#	# alternatively			 
#	            index <- temp[,i_index] 
#	            index <- as.POSIXlt(index,format=format,tz=tz)
#				if (!is.null(start_date) & !is.null(end_date)) { 
#	
#	 					temp <- temp[index>=start_date & index<=end_date,]
#				}
#				
#	 			temp<- temp[,-i_index]
#	 			temp <- as.zoo(temp)
#	  			index(temp) <- index
#	#
#	#
				 
				 
				 
				 
			 }
			 
			
			 
			 out[[i]] <- temp
		 }
		
		 names(out) <- filepath 
		 
		## if (length(out)==1) out <- out[[1]] COMMENTED BY EC ON 20150313
		 
		 if (ContinuousRecoveryCond) {
			 
			 names_keys <- paste(keyword,formatter,sep="")
			 names_keys <- sprintf(names_keys,1:length(names_points))
			 out <- lapply(X=names_keys, FUN=function(x,list){
						 
				index <- str_detect(names(list),x)
				list <- list[index]	 
					
				out <- list[[1]]
					
				for (it in list[-1]) {
						
					out <- rbind(out,it)
				}
					
				return(out)	 
					 
			},list=out)  
			 
			for (i in 1:length(out)) {
				
				temp <- out[[i]]
				if (!is.null(date_field) & !is.na(date_field) & length(i_index)==1 & length(date_field)>0) {
					
					index <- temp[,i_index]
					temp<- temp[,-i_index]
					index <- as.POSIXlt(index,format=format,tz=tz)
					temp <- as.zoo(temp)
					index(temp) <- index
					# insert sart date & date
					if (!is.null(start_date) & !is.null(end_date)) { 
						
					#	print(index(temp)>=start_date & index(temp)<=end_date)
					temp <- temp[index(temp)>=start_date & index(temp)<=end_date,]
						
					}
					
				}
				
				
				itr <- which(index(temp)!=index(temp)[1])
				itr <- c(1,itr)
				temp <- temp[itr,]
			
				out[[i]] <- temp
				
			}  
			 
			 
			 
			 
		 }
		 
		 ##### 
		 
		 if (!is.null(date_field) & !is.na(date_field) & length(i_index)==1 & length(date_field)>0) {
			 
			 if (is.null(zlayer.formatter)) zlayer.formatter <- NA
			 
			
			 
			 
			 if (!is.na(zlayer.formatter)) {
				 
				 if (length(z_unit)>1) z_unit <- z_unit[1]
				 if (length(geotop_z_unit)>1) geotop_z_unit <- geotop_z_unit[1]
				 
				 if (z_unit=="millimeters") z_unit <- 0.001
				 if (z_unit=="centimeters") z_unit <- 0.01
				 
				 if (geotop_z_unit=="millimeters") geotop_z_unit <- 0.001
				 if (geotop_z_unit=="centimeters") geotop_z_unit <- 0.01
				 
				 
				 
				 zu <- geotop_z_unit/z_unit
				 
				 out <- lapply(X=out,FUN=function(x,zfrm,zu){
							 
							 out <- x[,str_detect(names(x),"X")]
							 
							 zval <- as.numeric(str_replace(names(out),"X",""))*zu
							 
							 names(out) <- sprintf(zfrm,zval)
							 
							 return(out)
					##		 zlayer.formatter
						 },zfrm=zlayer.formatter,z=zu)
				 
			 }
			 
		## ... to go on 	 
			 
		 }
		 
		 if (length(out)==1) out <- out[[1]] ## added by EC on 20150313
	} else 	if (add_wpath) {
		
		if (!is.null(wpath)) out <- paste(wpath,out,sep="/")
	}

	return(out)
	
}
