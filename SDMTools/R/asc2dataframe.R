#' Ascii Grid Files to Dataframe and Dataframe to Ascii Grid Files
#' 
#' \code{asc2dataframe} converts a list of Esri ascii grid formatted files to a
#' data.frame consisting of only locations with data.\cr\cr
#' \code{dataframe2asc} converts a data.frame or matrix with spatial data to
#' Esri ascii grid formatted files.
#' 
#' asc2dataframe: The ascii grid files can be read in gzip compress format. The
#' dataframe returned contains the X and Y coordinate columns followed by
#' columns of data.
#' 
#' dataframe2asc: If filenames is null, column names will be used. The
#' data.frame has to contain the Y and X coordinates and the data as columns.
#' The ascii grid files can be created as gzip compress format and would be
#' saved in the outdir.
#' 
#' @param filenames is a vector of file names
#' @param varnames is a vector of names for the output columns, and must be the
#' same length as files
#' @param tdata is the data.frame which has y, x coordinates (OR lat,lon) and
#' columns for the data to be output (MUST be in that order)
#' @param outdir is the output directory, the default is the current working
#' directory
#' @param gz boolean defining if the ascii grid files are gzip compressed
#' @return \item{asc2dataframe }{Returns a dataframe with XY coordinates and
#' the data of each ascii grid files, as columns.} \item{dataframe2asc }{
#' Returns an asc grid file for each data column within the data.frame.}
#' @author Lorena Falconi \email{lorefalconi@@gmail.com}
#' @examples
#' 
#' #Create 2 ascii files
#' y=seq(10,50,0.5)
#' x=seq(140,180,0.5)
#' cellsize=0.5
#' data1=sample(160,140)
#' data2=sample(158,140)
#' out1.asc=as.asc(matrix(data1,nc=y, nr=x), xll=min(x), yll=min(y), cellsize=cellsize)
#' out2.asc=as.asc(matrix(data2,nc=y, nr=x), xll=min(x), yll=min(y), cellsize=cellsize)
#' #write the ascii files to the work directory
#' write.asc(out1.asc, 'out1.asc')
#' write.asc(out2.asc, 'out2.asc')
#' #list the ascii files
#' ascfiles=c('out1.asc', 'out2.asc')
#' #generate a dataframe from the ascii files
#' tdata=asc2dataframe(ascfiles)
#' tdata
#' 
#' #remove the files
#' unlink('out1.asc'); unlink('out2.asc')
#' 
#' #convert the dataframe tdata to ascii grid files
#' dataframe2asc(tdata)
#' 
#' #remove the files
#' unlink('var.1.asc'); unlink('var.2.asc')
#' 
#' @export 
asc2dataframe = function(filenames,varnames=NULL,gz=FALSE) {
	#check values
	if (is.null(varnames)) {
		varnames = paste('var',1:length(filenames),sep='.') #if no variable names defined, create them
	} else {
		if (length(varnames)!=length(filenames)) stop('variable names must be the same length as the files vector')
		varnames = as.character(varnames)
	}
	out = NULL #define the output
	for (ii in 1:length(filenames)) { #cycle through each of the files
		tfile = filenames[ii]; varname = varnames[ii] #define the file and variable name
		cat('working with',tfile,'...\n')
		if (!file.exists(tfile)) { warning(paste(tfile,'does not exist and was not used')); next } #check if the file exists, if not move to the next file
		tasc = read.asc(tfile,gz=gz) #read in the ascii grid file
		if (is.null(out)) { #if out is still null, populate the row/col/x,y data
			out = as.data.frame(which(is.finite(tasc),arr.ind=T)) #get row column info for actual data
			out$y = getXYcoords(tasc)$y[out$col]#extract the longitudes
			out$x = getXYcoords(tasc)$x[out$row] #extract the latitudes
		}
		out[varname] = tasc[cbind(out$row,out$col)] #append the actual data
	}
	if (is.null(out)) { #if out is still null
		warning('no data was extracted...'); return(NA)
	} else { #if out has some data
		out = na.omit(out) #remove any missing data
		out$row = out$col = NULL #remove the row/col info
		attr(out,'filenames') = list(filenames,names=varnames) #set an atrtribute relating filenames with variable names
		#names(attr(out,'filenames')) = varnames #define the names of the attribute list
		return(out)
	}
}

#' @rdname asc2dataframe
#' @export
dataframe2asc = function(tdata,filenames=NULL,outdir=getwd(),gz=FALSE) {
	
	#check values
	if (is.null(filenames)) {
		filenames = colnames(tdata)[3:length(tdata)] #if no variable names defined, create them
	} else {
		if (length(filenames)!=length(3:length(tdata))) stop('variable names must be the same length as the files vector')
		filenames = as.character(filenames)
	}	
	for (ii in 3:(length(tdata))) { #cycle through each of the files
		lats=unique(tdata[,1]);lats=sort(lats);
		longs=unique(tdata[,2]);longs=sort(longs)
		cellsize = min(c(diff(lats),diff(longs)))  # set cell size
		nc=ceiling((max(lats)-min(lats))/cellsize)+1; nr=ceiling((max(longs)-min(longs))/cellsize)+1
		out.asc=as.asc(matrix(NA,nrow=nr,ncol=nc),xll=min(longs),yll=min(lats),cellsize=cellsize)
		out.asc = put.data(tdata[,c(2:1,ii)],out.asc)
		write.asc(out.asc,paste(outdir,'/',filenames[ii-2],sep=''),gz=gz)  #put the name and extention
	}
}
