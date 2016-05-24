`dlcross.other` <- function(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile)
{
  results <- list()

   # two possibilities - 3 files input OR 3 objects input
   if (!missing(pheobj)&!missing(phefile)) 
	message("You have input both a datafile and object; the object will be used, please leave missing if the data differs from the file")

   if (!missing(genobj)&!missing(genfile)) 
	message("You have input both a datafile and object; the object will be used, please leave missing if the data differs from the file")

   if (!missing(mapobj)&!missing(mapfile)) 
	message("You have input both a datafile and object; the object will be used, please leave missing if the data differs from the file")

   if (missing(pheobj) & !missing(phefile))
	pheobj <- read.table(phefile, header=TRUE, check.names=FALSE)

   if (missing(genobj) & !missing(genfile))
	genobj <- read.table(genfile, header=TRUE, check.names=FALSE)

   if (missing(mapobj) & !missing(mapfile))
	mapobj <- read.table(mapfile, header=TRUE, check.names=FALSE)

   if (missing(pheobj)&missing(phefile)) 
	stop("Must input either pheobj or phefile; required argument")

   if (missing(genobj)&missing(genfile)) 
	stop("Must input either genobj or genfile; required argument")

   if (missing(mapobj)&missing(mapfile)) 
	stop("Must input either mapobj or mapfile; required argument")

   if (!(is.data.frame(genobj)|is.matrix(genobj))) stop("genobj has the wrong format, should be a matrix or data frame")
   if (!(is.data.frame(pheobj)|is.matrix(pheobj))) stop("pheobj has the wrong format, should be a matrix or data frame")
   if (!(is.data.frame(mapobj)|is.matrix(mapobj))) stop("mapobj has the wrong format, should be a matrix or data frame")

   if ((ncol(mapobj)<2)|(ncol(mapobj)>3))
	stop("Map has the wrong number of columns; please examine the documentation for the correct input format more closely")

   if (any(is.na(match(mapobj[,1], colnames(genobj)[2:ncol(genobj)]))))
	stop("There are different markers in the map file than in the genotype file") 

  if (any(mapobj[,1]!=colnames(genobj)[2:ncol(genobj)]))
  {
  	message("Warning: order of markers in map file is different to that in the genotype file. Markers will be rearranged in map order")
	genobj[,2:ncol(genobj)] <- genobj[,match(as.character(mapobj[,1]), colnames(genobj)[2:ncol(genobj)])+1]
	colnames(genobj)[2:ncol(genobj)] <- as.character(mapobj[,1])
  }

  idname <- colnames(genobj)[1]
  if (idname!=colnames(pheobj)[1]) 
	stop("The first column of the genotype and phenotype files, which should contain the unique identifier, has different names in the two files")

  ### alter data for trick
  if (ncol(genobj) > nrow(genobj)*length(table(mapobj[,2]))) {
	## recenter data
	genobj[,2:ncol(genobj)] <- as.matrix(apply(genobj[,2:ncol(genobj)], 2, 
	function(x) x-mean(x, na.rm=TRUE)))

	## set missing values to 0
	genobj[is.na(genobj)] <- 0
  }

  # create map object
  map <- list()
  if (ncol(mapobj)==2) mapobj <- cbind(mapobj, rep(NA, nrow(mapobj)))
  chrnames <- names(table(mapobj[,2]))
  for (i in 1:length(chrnames)){
	map[[i]] <- mapobj[mapobj[,2]==chrnames[i],3]
	map[[i]] <- mapobj[order(mapobj[mapobj[,2]==chrnames[i], 3]), 3]
	names(map[[i]]) <- mapobj[order(mapobj[mapobj[,2]==chrnames[i], 3]), 1]
	if (all(is.na(map[[i]]))) names(map[[i]]) <- mapobj[mapobj[,2]==chrnames[i],1]
  }
  names(map) <- chrnames
  class(map) <- "map"

  dfm <- as.data.frame(genobj)
  dfe <- as.data.frame(pheobj)
  
  dfm[,1] <- dfm[,which(names(dfm)==idname)]

  # wouldn't you have the id variable as well in dfm?
  chr <- match(mapobj[,2], names(map))
  names(dfm)[2:ncol(dfm)] <- paste("C", chr, "M", c(1:length(unlist(map))), sep="")

  # need to merge genotypic and phenotypic data
  dfe <- cbind(ord=1:nrow(dfe), dfe)
  df.anal <- merge(dfe, dfm, by=idname, all.x=TRUE, sort=FALSE)
  df.anal[,(ncol(dfe)+1):ncol(df.anal)] <- apply(df.anal[,(ncol(dfe)+1):ncol(df.anal)], 2, as.numeric)
  df.anal <- df.anal[order(df.anal$ord), -2]

  results$dfMerged <- df.anal

  results$dfMrk <- dfm
  results$map <- map
  results$mapp <- map  
  results$nphe <- ncol(pheobj)
  results$idname <- idname
  results$loc <- FALSE

  class(results) <- "dlcross"
  attr(results, "type") <- "other"

  results
}

