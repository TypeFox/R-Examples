`dlcross.dlmap` <- function(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile, type, step, fixpos, estmap, ...)
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

   if (missing(type)) 
	stop("Must input type")

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



  # check if more than two genotype values 
  alleles <- names(table(unlist(genobj[,2:ncol(genobj)])))
  if (length(alleles)>3) stop("More than two genotypes have been input; this data is not supported by dlmap")

  message("Genotype ", alleles[1], " will be coded as 0 in the data")
  message("Genotype ", alleles[2], " will be coded as 1 in the data")

  # BEH 10/09/09
  if (length(alleles)>2)
  message("Genotype ", alleles[3], " will be coded as 2 in the data")

  genobj2 <- apply(genobj, 2, function(x) return(abbreviate(as.character(x))))
  genobj[,2:ncol(genobj)][genobj2[,2:ncol(genobj2)]==alleles[1]] <- 0
  genobj[,2:ncol(genobj)][genobj2[,2:ncol(genobj2)]==alleles[2]] <- 1
  if (length(alleles)==3)
  genobj[,2:ncol(genobj)][genobj2[,2:ncol(genobj2)]==alleles[3]] <- 2

  # create map object
  map <- list()
  if (ncol(mapobj)==2) mapobj <- cbind(mapobj, rep(NA, nrow(mapobj)))
  chrnames <- names(table(mapobj[,2]))
  for (i in 1:length(chrnames)){
	map[[i]] <- mapobj[mapobj[,2]==chrnames[i],3]
	map[[i]] <- mapobj[order(mapobj[mapobj[,2]==chrnames[i], 3]), 3]
	names(map[[i]]) <- mapobj[order(mapobj[mapobj[,2]==chrnames[i], 3]), 1]
  }
  names(map) <- chrnames
  class(map) <- "map"
  	
  inCross <- list()  
  inCross$geno <- list()
  inCross$pheno <- as.data.frame(pheobj[1:nrow(genobj),])

  if (nrow(genobj)<nrow(pheobj)) 
     message("Warning: input$genCross does not contain all the input phenotypic data, see input$dfMerged")
 
  if (length(setdiff(inCross$pheno[[idname]], pheobj[[idname]]))>0)
     message("Warning: Individuals exist with genotype data but no phenotypes and will not be considered in analysis")

  # 10/09/09 need to make sure format for genotypes matches properly; e.g. AA, AB for backcross, etc for F2
  for (i in 1:length(map))
  {
    inCross$geno[[i]] <- list()
  if (chrnames[i]=="X")	
    class(inCross$geno[[i]]) <- "X" else class(inCross$geno[[i]]) <- "A"
  inCross$geno[[i]]$data <- as.matrix(genobj[,match(names(map[[i]]), colnames(genobj))])
  #	colnames(inCross$geno[[i]]$data) <- names(map[[i]])
  # may be unnecessary if values are being inherited from genobj
  # note, will need to run through this with a dataset and see how it looks
  # these are not correct for bc/RIL/etc. Check read.cross
  #    	  inCross$geno[[i]]$data[inCross$geno[[i]]$data==2] <- "BB"
  #    	  inCross$geno[[i]]$data[inCross$geno[[i]]$data==1] <- "AB"
  #    	  inCross$geno[[i]]$data[inCross$geno[[i]]$data==0] <- "AA"
  inCross$geno[[i]]$map <- map[[i]]
  }
  names(inCross$geno) <- names(map)

  inCross$pheno[[idname]] <- genobj[,1]
  class(inCross) <- c(type, "cross")
  if (ncol(mapobj)==2) estmap <- TRUE

  if (is.na(match(idname, names(inCross$pheno))))
	stop("ID variable not found in marker data")
 
  # if pheobj is missing, take it from input files
  if (missing(pheobj)) pheobj <- inCross$pheno

  if (!missing(pheobj)) {
  if (is.na(match(idname, names(pheobj))))
	stop("ID variable not found in supplemental phenotype object")
  if (length(intersect(names(table(inCross$pheno[[paste(idname)]])), names(table(pheobj[[paste(idname)]]))))==0)
	stop("ID variables in marker and phenotype data do not coincide")
  if (is.na(match(idname, names(pheobj))))
	  stop("Error: ID variable is not in phenotype dataset")

  idIndex <- match(idname, names(pheobj))
  }
   
  exp.dat <- fill.geno(inCross, method="argmax") 

  if (estmap)
  {
     message("Map has been replaced by estimated map")
     map <- est.map(exp.dat)            
     exp.dat <- replace.map(exp.dat, map)
  }

  if (step<0) 
	stop("Invalid input: step size cannot be < 0")
  if (fixpos<0)
	stop("Invalid input: fixpos cannot be < 0")
  if ((step>0) & (fixpos>0))
	stop("Invalid input: only one of step and fixpos should be > 0")

  gp <- calcpos(exp.dat, step, fixpos)
  mapp <- lapply(gp$geno, function(x) return(attr(x$prob, "map")))
  coln <- unlist(lapply(mapp, names))
  chr <- vector()
  for (i in 1:length(mapp))
	chr <- c(chr, rep(i, length(mapp[[i]])))

  mmap <- vector()
  pmap <- vector()
  for (i in 1:length(map)) {
	mmap <- c(mmap, 1:length(map[[i]]))
	if (length(mapp[[i]])>length(map[[i]]))
	pmap <- c(pmap, 1:(length(mapp[[i]])-length(map[[i]])))
  }

   if (type=="f2")
   {
   dat1 <- lapply(gp$geno, function(x) return(x$prob[,,2]))
   dat2 <- lapply(gp$geno, function(x) return(x$prob[,,3]-x$prob[,,1]))
   dfm1 <- do.call("cbind", dat1)

   pos <- grep("loc", coln)
   mrk <- setdiff(1:length(coln), pos)
   colnames(dfm1)[pos] <- paste("C", chr[pos], "P", pmap, "D", sep="")
   colnames(dfm1)[mrk] <- paste("C", chr[mrk], "M", mmap, "D", sep="")

   dfm2 <- do.call("cbind", dat2)
   colnames(dfm2)[pos] <- paste("C", chr[pos], "P", pmap, "A", sep="")
   colnames(dfm2)[mrk] <- paste("C", chr[mrk], "M", mmap, "A", sep="")
  
   dfm <- matrix(nrow=nrow(dfm1), ncol=ncol(dfm1)+ncol(dfm2))
   dfm[,seq(1, ncol(dfm), 2)] <- dfm1
   dfm[,seq(2, ncol(dfm), 2)] <- dfm2
   dfm <- as.data.frame(dfm)
   names(dfm)[seq(1, ncol(dfm), 2)] <- colnames(dfm1)
   names(dfm)[seq(2, ncol(dfm), 2)] <- colnames(dfm2)
   } else {
   dat <- lapply(gp$geno, function(x) return(x$prob[,,2]))
   dfm <- as.data.frame(do.call("cbind", dat))
   pos <- grep("loc", coln)
   mrk <- setdiff(1:length(coln), pos)
   names(dfm)[pos] <- paste("C", chr[pos], "P", pmap, sep="")
   names(dfm)[mrk] <- paste("C", chr[mrk], "M", mmap, sep="")
   }

   dfe <- as.data.frame(pheobj)
   dfm <- cbind(inCross$pheno[[idname]], dfm)
   names(dfm)[1] <- idname

   # as it stands dfm contains both positions and markers; really just want mrk
   # also, need to be careful about doing the permutation on the F2

   dfe <- cbind(ord=1:nrow(dfe), dfe)
   df.anal <- merge(dfe, dfm, by=idname, all.x=TRUE, sort=FALSE)
   df.anal[,(ncol(dfe)+1):ncol(df.anal)] <- apply(df.anal[,(ncol(dfe)+1):ncol(df.anal)], 2, as.numeric)
   df.anal <- df.anal[order(df.anal$ord), -2]

   results$dfMerged <- df.anal

   results$dfMrk <- dfm[,c(1,(2:ncol(dfm))[grep("M", names(dfm)[2:ncol(dfm)])])]

   results$mapp <- mapp
   results$map <- pull.map(inCross) 
   results$genCross <- gp
   results$nphe <- ncol(pheobj)
   results[["idname"]] <- idname

   if (step==0 & fixpos==0) results$loc <- FALSE else results$loc <- TRUE

   class(results) <- "dlcross"
   attr(results, "type") <- type
  results
}

