createImage <-
function(grid, landwatermask, image.title, directory, filename, shift, parts=10, resolution=1){
	if ((dim(grid)[1] != dim(landwatermask)[1])||(dim(grid)[2] != dim(landwatermask)[2])){
		return(FALSE)
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

	#create grid
	result <- ifelse(landwatermask == -1, -1, grid)
	dimension <- dim(result)

	#create dataframe
	richness <- vector(mode="integer", length=dimension[1]*dimension[2])
	y <- x <- vector(mode="numeric", length=dimension[1]*dimension[2])
	count <- 1
	for (m in 1:dimension[2]){
		for (n in 1:dimension[1]){
			richness[count] <- result[n,m]
			x[count] <- shift[1] + resolution*n
			y[count] <- shift[2] + resolution*m
			count <- count + 1
		}
	}
	
	maxrich <- ceiling(max(richness))
	steps <- round(maxrich/parts, digits=1)
	breaks <- seq(from=0, by=steps, length=(parts+1))

	data <- data.frame(richness=richness, x=x, y=y)
	
	#create colors
	if ((parts > 2) && (parts < 12)) {
		colors <- RColorBrewer::brewer.pal(n=parts, name="PuOr")
	} else {
		colors <- rainbow(parts)
	}

	#create levelplot
	png(filename=paste(directory,filename,sep=""), width=2000, height=2000, res=100)
	lattice::trellis.par.set(fontsize=list(text=25,points=8))
	picture <- lattice::levelplot(data$richness ~ data$x + data$y, data, at=c(-2,-1,breaks), contour=FALSE, 
			col.regions=c("lightgray","gray",colors), 
			xlab="Longitude",ylab="Latitude", main=image.title,
			colorkey=list(at=breaks, col=colors,
			labels=list(labels=breaks,at=breaks)))
	print(picture)
	dev.off()

	return(TRUE)
}
