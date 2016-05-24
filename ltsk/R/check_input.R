check_input <-
function(in_data,xcoord,ycoord,tcoord,zcoord)
{
   l.coords <- c(xcoord,ycoord,tcoord,zcoord)
   if( is.null(dim(in_data)))
   {
	## data is a vector
	l.data <- matrix(in_data,ncol=length(in_data),nrow=1)
	colnames(l.data) <- names(in_data)
   }
   else
   {
	## query is a matrix
	l.data <- in_data
   }
   
   for (l.coord in l.coords)
   {
	if (is.na(match(l.coord,colnames(l.data)))){
		l.txt <- colnames(l.data)
		l.data <- cbind(l.data,NA)
		##warning(l.coord,"not found in input data!\n")
		colnames(l.data) <- c(l.txt,l.coord)
	}
   }		
   r <- as.matrix(l.data[,l.coords])
   
   ## check for duplicates
   isdup <- duplicated(r)
   if(any(isdup)){
	 arg <- deparse(substitute(in_data))
	 warning(sum(isdup), " duplicates removed from ", arg,"\n")
	 r <- r[!isdup,]
   }
   
   r
}

check_na <- function(in_data,type_txt)
{
  if( is.null(dim(in_data)))
  {
    ## data is a vector
    l.data <- matrix(in_data,ncol=length(in_data),nrow=1)
    colnames(l.data) <- names(in_data)
  }
  else{
    l.data <- in_data
  }
  na.omit(l.data)
  ## no_na <- na.omit(l.data)
  
  #dd <- attr(no_na,"na.action")
  #if (length(dd)>0){
	#	cat(length(dd)," ",type_txt," location removed.\n")
	#}
	## no_na
}