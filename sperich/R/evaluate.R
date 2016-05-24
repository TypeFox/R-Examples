evaluate <-
function(result.grid.one, result.grid.two, title.one="Histogramm of species richness", title.two="Histogramm of species richness", xmax=400, ymax=1000, directory=getwd(), filename="histogramm.png"){
	#check result grids
	numberofgrids <- 2
	if (is.null(result.grid.one)){
		numberofgrids <- numberofgrids - 1
	}
	if (is.null(result.grid.two)){
		numberofgrids <- numberofgrids - 1
	}

	#check directory
	direc <- unlist(strsplit(directory,""))
	if (direc[length(direc)] != "/"){
		direc <- c(direc, "/")
	}
	directory <- ""
	for (i in 1:length(direc)){
		directory <- paste(directory,direc[i], sep="")
	}
	#check filename
	file <- unlist(strsplit(filename,""))
	filename <- ""
	if ((file[length(file)-3]!=".")||(file[length(file)-2]!="p")||(file[length(file)-1]!="n")||(file[length(file)]!="g")){
		file <- c(file,".","p","n","g")
	}
	for (i in 1:length(file)){
			filename <- paste(filename,file[i], sep="")
	}
	
	#check number of grids
	if (numberofgrids == 0) {
		return(FALSE)
	}	

	#one grid -> one histogramm will be created
	if (numberofgrids == 1) {
		#which grids is not NULL
		if (!is.null(result.grid.one)){
			title <- title.one
			result.grid <- result.grid.one
		} else {
			title <- title.two
			result.grid <- result.grid.two
		}

		#create histogramm
		x <- result.grid[which(result.grid > 0)]

		png(paste(directory, filename, sep=""), width=500, height=500)
		cmrg <- ceiling(max(result.grid))
		cmrg <- ifelse((cmrg %% 5) == 0, cmrg, cmrg + (5- (cmrg %% 5)))
	 	steps <- ifelse(cmrg > 50, 5, 1)
		hist(x, breaks=seq(from=0, to=cmrg, by=steps), main=title, 
				xlab="number of species", xlim=c(0,xmax), ylim=c(0,ymax), las=1, 
				col=gray(.7))

		#quantils and sum over all cells
		result <- c(summary(x),round(sum(x), digits=2))
		names(result) <- c(names(summary(x)), "sum")

		result.string <- ""
		for (n in 1:7){
			result.string <- paste(result.string, names(result)[n],": ", result[n],"\n", sep="") 
		}

		#add summary-text to plot
		text(0.8*xmax, 0.8*ymax, result.string)

		dev.off()
		return(TRUE)
	}

	#two grids -> two histogramms will be created
	if (numberofgrids == 2){
		png(paste(directory, filename, sep=""), width=1000, height=500)
		#double window
		layout(rbind(c(1,2)),heights=c(1), respect=FALSE)

		#create the two plots
		while (numberofgrids > 0){
			numberofgrids <- numberofgrids - 1

			if (!is.null(result.grid.one)){
				title <- title.one
				result.grid <- result.grid.one
				result.grid.one <- NULL
			} else {
				title <- title.two
				result.grid <- result.grid.two
			}
	
			#create histogramm
			x <- result.grid[which(result.grid > 0)]
			cmrg <- ceiling(max(result.grid))
			cmrg <- ifelse((cmrg %% 5) == 0, cmrg, cmrg + (5- (cmrg %% 5)))
			steps <- ifelse(cmrg > 50, 5, 1)
			hist(x, breaks=seq(from=0, to=cmrg, by=steps), main=title, 
				xlab="number of species", xlim=c(0,xmax), ylim=c(0,ymax), las=1, 
				col=gray(.7))

			#quantils and sum over all cells
			result <- c(summary(x),round(sum(x), digits=2))
			names(result) <- c(names(summary(x)), "sum")

			result.string <- ""
			for (n in 1:7){
				result.string <- paste(result.string, names(result)[n],": ", result[n],"\n", sep="") 
			}
			#add summary-text to plot
			text(0.8*xmax,0.8*ymax,result.string)
		}
		dev.off()

		return(TRUE)
	}
}
