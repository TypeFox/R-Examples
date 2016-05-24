`dlcross.cross` <- function(genobj, pheobj, idname, step, fixpos, estmap)
{
  results <- list()

  if (missing(genobj)) stop("Must input a cross object if format is rqtl")
  if (!inherits(genobj,"cross")) stop("Input object does not have class cross")
  inCross <- genobj
  
  type <- class(inCross)[1]

  if (is.null(inCross$pheno[[idname]])) {
	cat("missing idname in cross object; will generate automatically\n")
	inCross$pheno[[idname]] <- paste("S", 1:nind(inCross), sep="")
  }

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
  map <- pull.map(exp.dat)

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
   results$map <- map 
   results$genCross <- gp
   results$nphe <- ncol(pheobj)
   results[["idname"]] <- idname

   if (step==0 & fixpos==0) results$loc <- FALSE else results$loc <- TRUE

   class(results) <- "dlcross"
   attr(results, "type") <- type
  results
}

