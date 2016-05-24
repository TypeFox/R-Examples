NCEP.array2df <- function(wx.data, var.names=NULL){
	## This is a function to convert the data contained in one or many arrays to a dataframe ##
	## wx.data should be a list containing all of the arrays to be put in the dataframe ##
	## The names of the variables in the resulting dataframe can be specified in var.names ##
	## The order of the names should correspond to the order of the list of input arrays in wx.data ##
	
	## Make sure that wx.data is a list ##
	if(is.list(wx.data) == FALSE) { wx.data <- list(wx.data) }
	
	## Specify the number of arrays which will be made into a data.frame ##
	n.arrays <- length(wx.data)
	## And the number of dimensions in the arrays ##
	n.dims <- length(dimnames(wx.data[[1]]))
	
	## Make sure that the arrays contain the same number of dimensions ##
		for(i in 1:n.arrays){
			if(length(dimnames(wx.data[[1]])) > 3) { stop('Arrays contain more than three dimensions') }
			if(length(dimnames(wx.data[[1]])) != length(dimnames(wx.data[[i]]))) { stop('Multiple arrays must contain exactly the same dimensions') }
		}
		
	## Make sure that the dimensions are exactly the same in all arrays ##
		for(i in 1:n.arrays){
			test1 <- length(which(dimnames(wx.data[[1]])[[1]] != dimnames(wx.data[[i]])[[1]]))
			test2 <- length(which(dimnames(wx.data[[1]])[[2]] != dimnames(wx.data[[i]])[[2]]))
			if(n.dims == 3) { test3 <- length(which(dimnames(wx.data[[1]])[[3]] != dimnames(wx.data[[i]])[[3]])) } else { test3 <- 0 }
			if(any(c(test1, test2, test3) != 0)) { stop('Multiple arrays must contain exactly the same dimensions') }
		}

	## Determine the names of the variables other than datetime, latitude, and longitude ##
	if(is.null(var.names) == FALSE && length(var.names) != n.arrays) { warning('Number of variable names does not equal the number of input arrays.\nGenerating names automatically') }
	if(is.null(var.names) | length(var.names) != n.arrays) {
		for(i in 1:n.arrays){
			var.names <- append(var.names, paste("variable",i, sep=''))
			}
		}
	
	## Convert the data in these arrays to vectors ##
		for(i in 1:length(var.names)){
			assign(var.names[i], as.vector(wx.data[[i]]))
			}
		
	## Calculate the dimensions of the input arrays ##
	rows <- length(dimnames(wx.data[[1]])[[1]]) ## Latitudes
	columns <- length(dimnames(wx.data[[1]])[[2]]) ## Longitudes
	layers <- ifelse(n.dims == 3, length(dimnames(wx.data[[1]])[[3]]), 1) ## Datetimes

	## Create a dataframe from these arrays ##
	if(n.dims == 3) {
		output.df <- data.frame(datetime=rep(dimnames(wx.data[[1]])[[3]], each=rows*columns), 
				latitude=rep(as.numeric(dimnames(wx.data[[1]])[[1]]), columns*layers), 
				longitude=rep(rep(as.numeric(dimnames(wx.data[[1]])[[2]]), each=rows), layers))
				} else {
		output.df <- data.frame(latitude=rep(as.numeric(dimnames(wx.data[[1]])[[1]]), columns*layers), 
				longitude=rep(rep(as.numeric(dimnames(wx.data[[1]])[[2]]), each=rows), layers))
				}				

	## Specify that the datetime variable should be considered a character object not a factor ##
	if(n.dims == 3) { output.df$datetime <- as.character(output.df$datetime) }

	## Now include all of the data ##
	for(i in 1:n.arrays){
		eval(parse(text=paste("output.df$", var.names[i], " <- ", var.names[i], sep='')))
			}
			
	## Return the dataframe ##
	return(output.df)
	}