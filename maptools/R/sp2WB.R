# 2005 (c) Virgilio Gomez Rubio, partly derived from
# earlier code by Thomas Jagger

# The code in this file exports and sp object into a S-Plus map
# format to be import by WinBUGS.

# file: file where output is written
# map: sp object (SpatialPolygons object)
# xScale/yScale: scales to be written in the output file

sp2WB <- function(map, filename, Xscale=1, Yscale=Xscale, plotorder=FALSE) {

# Write some tests here to ensure that all the objects passed
# are of the appropriate type

	f<-file(filename,"w")

# Get the total number of areas and Rings
	SRings<-slot(map, "polygons")
	nareas<-length(SRings)
	nRings<-sum(sapply(SRings, function(x) length(slot(x, "Polygons"))))
	IDs<-sapply(SRings, function(i) slot(i, "ID"))

# Plot header of the f
	cat(file=f, "map:",nareas, "\n", sep="")
	cat(file=f, "Xscale:", Xscale, "\n", sep="")
	cat(file=f, "Yscale:", Yscale, "\n", sep="")
	cat(file=f, "\n")

	if(plotorder)
		porder <- slot(map, "plotOrder")
	else
		porder<-1:nareas

# Different
	for(area in (1:nareas)[porder])
		cat(file=f, area, " area", IDs[area], "\n", sep="")

	cat(file=f, "\n")

	index<-1
# Loop to print all the individual rings
	for(area in (1:nareas)[porder]) {
		label<-paste("area", IDs[area], sep="")
		Rings<-slot(SRings[[area]], "Polygons")
		lRings<-length(Rings)

		if(plotorder)
			porderrings<-slot(Rings, "plotOrder")
		else
			porderrings<-1:lRings

		for(ring in (1:lRings)[porderrings]) {
			coords<-slot(Rings[[ring]], "coords")
			ncoords<-length(coords[,1])
		
# Should we check that there are only x/y coordinates?
                        xcrd <- formatC(coords[,1], format="f")
                        ycrd <- formatC(coords[,2], format="f")
			for(i in 1:ncoords)
                                cat(file=f, label, xcrd[i], ycrd[i],
                                        "\n")

			if(index<nRings)
				cat(file=f, "NA", "NA", "NA\n")	
			else
				cat(file=f, "\n")	

			index<-index+1
		}

	}

	cat(file=f, "END\n")
	close(f)
	invisible(NULL)
}

