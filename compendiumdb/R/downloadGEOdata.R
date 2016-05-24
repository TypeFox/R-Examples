downloadGEOdata <-
function(GSEid, destdir = getwd()){
	
  options(warn=-1)
  dir <- path.package("compendiumdb")

  scriptLoc <- paste(dir,"/scripts/R",sep="")
  scriptLoc <- gsub("^","\"",scriptLoc)
  scriptLoc <- gsub("$","\"",scriptLoc) 

  dataLoc <- destdir

  if(!length(dir(dataLoc,pattern="^BigMac$")))
    {
      dir.create(paste(dataLoc,"/BigMac",sep=""))
      dir.create(paste(dataLoc,"/BigMac/data",sep=""))
      dir.create(paste(dataLoc,"/BigMac/log",sep=""))
      
      dir.create(paste(dataLoc,"/BigMac/data/GEO",sep=""))
      dir.create(paste(dataLoc,"/BigMac/data/GEO/GPL",sep=""))
      dir.create(paste(dataLoc,"/BigMac/data/GEO/GSM",sep=""))
      dir.create(paste(dataLoc,"/BigMac/data/GEO/GSE",sep=""))

      dir.create(paste(dataLoc,"/BigMac/annotation",sep=""))
      dir.create(paste(dataLoc,"/BigMac/annotation/configuration",sep=""))

      dir.create(paste(dataLoc,"/BigMac/data/GEO/GDS_description",sep=""))
      
      x <- ""
      write(x,paste(dataLoc,"/BigMac/log/log_getAllForGSE.txt",sep=""))
      write(x,paste(dataLoc,"/BigMac/log/report_getAllForGSE.txt",sep=""))

      file.copy(paste(dir,"/extdata/configuration.txt",sep=""),paste(dataLoc,"/BigMac/annotation/configuration/configuration.txt",sep=""))
      file.copy(paste(dir,"/scripts/BigMac/COMPENDIUM",sep=""),paste(dataLoc,"/BigMac",sep=""),recursive=TRUE)
    }else if(!file.exists(paste(dataLoc,"/BigMac/COMPENDIUM/compendiumSUB.pm",sep=""))){
	stop("Please check if the BigMac directory has the right directory structure to load the data")
    }

  dataLoc <- gsub("^","\"",dataLoc)
  dataLoc <- gsub("$","\"",dataLoc)	

  plFile <- paste(dir,"/scripts/Perl/getAllForGSE.pl",sep="")
  plFile <- gsub("^","\"",plFile)
  plFile <- gsub("$","\"",plFile)

  system(paste("perl -I",paste(dataLoc,"/BigMac/COMPENDIUM",sep=""),plFile,GSEid,dataLoc,scriptLoc))

  files <- list.files(destdir)
  files <- paste(destdir,"/",files,sep="")
  i <- grep("GPL.*.soft$",files)
  j <- grep("GPL.*.annot.gz$",files)
  k <- grep("GDS.*.soft",files)
	
  ff <- files[c(i,j,k)]
  file.remove(ff)

  catch = file.remove(paste(destdir,"/BigMac/data/GEO/GSE/",GSEid,"_family.soft.gz",sep=""))
}

