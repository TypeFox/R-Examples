popgenreport <- function(cats=NULL,
                          
                            mk.counts=TRUE,   # this switch is to provide a population overview
                            mk.map=FALSE,        # this switch is for the map
                            maptype="satellite",
                            mapdotcolor ="blue",
                            mapdotsize=1,
                            mapdotalpha=0.4,
                            mapdottype=19 ,
                            mapzoom=NULL,
                          
#c("roadmap",
#"mobile",
#"satellite",
#"terrain",
#"hybrid",
#"mapmaker-roadmap",
#"mapmaker-hybrid")",
#                          
                          
                          mk.locihz=FALSE,     # this switch is to test for population heterozygosity
                          mk.hwe=FALSE,   # this switch is for population wide HWE

                          mk.fst=FALSE,        # this switch is to run FST tests on the full population
                          mk.gd.smouse=FALSE,   # this switch is to run the Smouse and Peakall genetic distances
                          mk.gd.kosman=FALSE,   # this switch is to run the Kosman and Leonard genetic distances

                          mk.pcoa=FALSE,
                          mk.spautocor=FALSE,
                          mk.allele.dist=FALSE, # this switch it to look at allele distributions by loci and pop
                          mk.null.all=FALSE,
                          mk.allel.rich=FALSE,

                          mk.differ.stats=FALSE ,     # this switch is to look at population differentiation statistics (Fst, Gst, etc)
                          mk.custom = FALSE,
                          fname="PopGenReport",
                          foldername="results",
                          path.pgr=NULL,
                          mk.Rcode=FALSE,       # make the code that was ran available as an R file
                          mk.complete=FALSE,    # create a full report)  
                          mk.pdf=TRUE)
{
  if (class(cats)!="genind") {cat("You did not provide a valid genind object! Script stopped!\n"); return;}
  
  # Check for combinations of populations and loci with low numbers of individuals and alleles  
  npops<-length(levels(cats@pop))
  nloci<-length(locNames(cats))
  
  # this splits bilby up into loci
  loci<-seploc(cats)
  
  # this further subdivides the loci into populations
  locipop<-lapply(loci,seppop)
  
  popsizes<-matrix(NA,nrow=nloci,ncol=npops)
  for (i in 1:nloci){
    for (j in 1:npops){
      popsizes[i,j]<-sum(!is.na(apply(locipop[[i]][[j]]@tab,1,sum)))
    }
  }
  
  for(i in 1:dim(popsizes)[2]){
    numlow<-length(which(popsizes[,i]<3))
    if(numlow>0) message("Population ",unname(popNames(cats))[i]," has ",numlow," locus/loci with less than 3 genotypes. This may cause errors in some analyses. We advice to combine or drop populations with low numbers of genotypes. ")
  }
  
### adjusting causes error with Adegenet 2.0 and it is controlled for when creating the object!!  
  #cut down length of loci names to  6  and make sure they are unique
  #locNames(cats) <- substr(locNames(cats),1,6)   
  #if (length(unique(locNames(cats)))!= length(locNames(cats))) 
  #{locNames(cats) <- paste(1:length(locNames(cats)),"-",substr(locNames(cats),1,4), sep="")
  #cat("Loci names were not unique and therefore adjusted.\n")
  #}
  #levels(cats@loc.fac) <- locNames(cats)  #make sure levels and factors are the same
#check if pop.names, indnames are unique!!!!
#adjust if necessary and issue a notification
# if (length(unique(indNames(cats)))!=length(indNames(cats))) 
#   {indNames(cats) <- paste(1:length(indNames(cats)),"-",substr(indNames(cats),1,8),sep="")
#   cat("Individual names were not unique and therefore adjusted.\n")
#   }
# 
# 
#   
# if (length(unique(popNames(cats)))!=length(popNames(cats))) 
#   {
#   popNames(cats) <- paste(1:length(popNames(cats)),"-",substr(popNames(cats),1,6),sep="")
#   cat("Subpopulation names were not unique and therefore adjusted.\n")
#   }
# 
#  
    #set directory where to save a file, defaults to tempdir (follow R policy)
  if (is.null(path.pgr)) 
  {
  path.pgr <- tempdir()

  }
  

  
#  setwd(path.pgr)
  
  #create a foldername folder if not existing...
  dirfiles <- list.dirs(path=path.pgr, recursive=FALSE)
  if (!(tolower (file.path(path.pgr,foldername))) %in% tolower(dirfiles)) {
    dir.create(file.path(path.pgr,foldername))
    cat("There is no ",foldername, " folder. I am trying to create it; \notherwise please create the folder manually. \n")
  }
  owd <-getwd()
  setwd(file.path(path.pgr, foldername))
  # conversion of lat longs to google map data (Mercator (dismo) wants to have long lat)
  

  coords=FALSE
  if (is.null(cats$other$latlong) & is.null(cats@other$xy)) coords=FALSE else {
  if (!is.null(cats@other$latlong)) cats@other$xy <- Mercator(cats@other$latlong[,c(2,1)])
  if (!is.null(cats@other$xy)) cats@other$latlong <- Mercator(cats@other$xy, inverse=TRUE)[,c(2,1)]  
  
  if ((nrow(cats@other$latlong) == length(indNames(cats))) & (nrow(cats@other$xy) == length(indNames(cats)) )) coords=TRUE
  


  } 
  # give cats a filename that can be seen in the snw chunks
  
  
  cats@other$filename<- fname
  cats@other$foldername<-foldername
  cats@other$path <- path.pgr
  #determine the type of map
  if ((mk.map==TRUE | mk.complete) & coords) 
  {
  cats@other$maptype=maptype
  cats@other$mapdotcolor =mapdotcolor
  cats@other$mapdotsize=mapdotsize
  cats@other$mapdotalpha=mapdotalpha
  cats@other$mapdottype=mapdottype
  cats@other$mapzoom=mapzoom
  
  }  
  
###################################
##### create a new environment to run knitr in it
 #pgr <- new.env(parent=.GlobalEnv)
 pgr <- new.env()
 assign("cats",cats,envir=pgr)
 #sumcats <- summary(cats, verbose=FALSE)
 #assign("sumcats",sumcats, envir=pgr )
###################################



  # save the data in a tempfile
 # save(cats, file=paste(foldername,"\\","tempcats.rdata",sep=""))
  
  #check path to the snw files
path <- NULL
  for(i in seq_along(.libPaths()))
{
  if (file.exists(paste(.libPaths()[i],"/PopGenReport/swchunks/header.snw",sep="")))  
  {
  path <-   paste(.libPaths()[i],"/PopGenReport/swchunks/", sep="" )
  break
  }
  
}
if (is.null(path)) {cat("Could not find snw files in the PopGenReport library folder. Please check if the package is installed correctly (e.g.  installed.packages()[\"PopGenReport\",2]). \n"); return;}
  #for testing:
  #path <- "d:\\bernd\\R\\popgenreport\\inst\\swchunks\\"
  #path<- "C:\\Aaron files\\popgenreport098\\PopGenReport_0.98\\PopGenReport\\swchunks\\"
  header.file <- readLines(paste(path,"header.snw",sep=""))
  required<- readLines(paste(path,"required.snw",sep=""))
  compl<-c(header.file,required) 
  
  cat("Compiling report...\n")
  if(mk.counts | mk.complete){
    cat("- General summary...\n")
    overview<-readLines(paste(path,"counts.snw",sep=""))
    compl<-c(compl,overview)
  }
  if (coords==FALSE) cat(" - No valid coordinates were provided. \n   Be aware you need to provide a coordinate (or NA) for each individual\n   and the coordinate heading in slot @other has to be 'latlong' or 'xy'.\n   Some of the analyses require coordinates and will be skipped!\n") 
  if ((mk.map==TRUE | mk.complete) & coords){
    cat("- Map of individuals...\n")  
    mapping<-  readLines(paste(path,"map.snw",sep=""))
    compl<-c(compl,mapping)
  }
  
  if (mk.locihz | mk.complete){
    cat("- Statistics on population heterogeneity ...\n")  
    popheterozygosity <- readLines(paste(path,"locihz.snw",sep=""))
    compl<-c(compl,popheterozygosity)
  }
  
  if (mk.allele.dist | mk.complete){
    cat("- Allelic distances ...\n")  
    numloci<-length(cats@loc.n.all)
    alleledistn <- readLines(paste(path,"allele.dist.snw",sep=""))
    compl<-c(compl,alleledistn)
}
if (mk.fst| mk.complete){
     cat("- Pairwise Fst ...\n")  
  popfst<-readLines(paste(path,"fst.snw",sep=""))
  compl<-c(compl,popfst)
}
if (mk.null.all | mk.complete){
     cat("- Checking for null alleles ...\n")  
     null.stat<-readLines(paste(path,"null.all.snw",sep=""))
  compl<-c(compl,null.stat)
}
if (mk.allel.rich | mk.complete){
     cat("- Allelic richness ...\n")  
     all.stat<-readLines(paste(path,"allel.rich.snw",sep=""))
  compl<-c(compl,all.stat)
}
if (mk.differ.stats | mk.complete){
     cat("- Pairwise differentiations ...\n")  
     diff.stat<-readLines(paste(path,"differ.stats.snw",sep=""))
  compl<-c(compl,diff.stat)
}

if (mk.hwe | mk.complete){
  cat("- Test for Hardy-Weinberg-Equilibrium ...\n") 
  cat("  !! You may get warnings when running HWE tests, if the test is based\n")
  cat(" on an entry in the chi-square table which is less than five.!! \n")
  popHWEll<-readLines(paste(path,"hwe.snw",sep=""))
  compl<-c(compl,popHWEll)
}

if ((mk.gd.kosman==TRUE | mk.complete) & coords){
  cat("- Kosman & Leonard 2005 genetic distances...\n")
  kosman<-readLines(paste(path,"gd.kosman.snw",sep=""))
  compl<-c(compl,kosman)
}

if ((mk.gd.smouse==TRUE | mk.complete) & coords){
  cat("- Smouse & Peakall 1999 genetic distances...\n")
  smouse<-readLines(paste(path,"gd.smouse.snw",sep=""))
  compl<-c(compl,smouse)
}
if ((mk.spautocor==TRUE | mk.complete) & coords){
  cat("- Spatial autocorrelation following Smouse & Peakall 1999 ...\n")
  spa<-readLines(paste(path,"spautocor.snw",sep=""))
  compl<-c(compl,spa)
}
if (mk.pcoa==TRUE | mk.complete){
  cat("- Principal coordinate analysis following Jombart et al. 2009...\n")
  pca<-readLines(paste(path,"pcoa.snw",sep=""))
  compl<-c(compl,pca)
}

if (mk.custom==TRUE){
  cat("- Run customised snw file, custom.snw ...\n")
  custom<-readLines(paste(path,"custom.snw",sep=""))
  compl<-c(compl,custom)
}



footer.file<-readLines(paste(path,"footer.snw",sep=""))  
compl<-c(compl,footer.file)

#compl <- c(header.file, required, loaddata, mapping, popheterozygosity, footer.file)


rnwfile <- paste(fname,".rnw",sep="")
texfile <-  paste(fname,".tex",sep="") 


zz <- file(file.path(path.pgr,foldername,rnwfile), "w")
writeLines(compl,zz)
close(zz) 


#setwd(paste(path.pgr,foldername, sep="/"))
cat(paste("Analysing data ...\n", sep=""))
#Sweave(paste(fname,".rnw",sep=""), output=paste(fname,".tex",sep=""), quiet=FALSE, driver=mydriver)
flush.console()
knit(input=rnwfile, output=texfile, quiet=TRUE, envir=pgr)

if (mk.pdf==TRUE)
{
cat(paste("Creating pdf from: ",rnwfile," ...\n",sep=""))
knit2pdf(texfile, texfile)
cat(paste("Finished.\nCheck ",fname,".pdf for results.\n", sep=""))
}

if (mk.Rcode) {
  cat(paste("Creating R code from: ",rnwfile,"...\n"), sep="")
  rfile <-paste(fname,".R",sep="")
  purl(input=rnwfile, output=rfile)
#  Stangle(paste(fname,".rnw",sep=""))
}
    

cat(paste("All files are available in the folder: \n",file.path(path.pgr, foldername),"\n",sep=""))

#reset working directory to previous
setwd(owd)
return(pgr$allresults)
}
