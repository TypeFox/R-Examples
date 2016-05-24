Bclim <-
function(pollen.loc,chron.loc,core.name="Core",time.grid=seq(0,14,length=100),rsfile,nchrons=10000,parallel=FALSE,save.as.you.go=TRUE) {
  
  if(!file.exists(pollen.loc)) stop(cat("Pollen data not found. Check",pollen.loc))
  if(!file.exists(chron.loc)) stop(cat("Chronology data not found. Check",chron.loc))
  
  pollen.data <- as.matrix(read.table(pollen.loc,header=TRUE))
  temp.chron <- read.table(chron.loc,nrows=5)
  
  # Some error handling
  if(ncol(pollen.data)!=28) stop("Number of pollen taxa should be 28")
  if(nrow(pollen.data)!=ncol(temp.chron)) stop("Number of columns in chronology should be the same 
                                               as the number of rows in the pollen data")
  
  ### Run Stage 1 - calculate layer posteriors
  cat("Stage 1 of 4: calculating layer posteriors...\n")
  step1filename <- paste(gsub(" ","_",core.name),"_step1.RData",sep="")
  if(!file.exists(step1filename)) {
    step1 <- BclimLayer(pollen.loc,rsfile)
    if(save.as.you.go) save(step1,file=paste(gsub(" ","_",core.name),"_step1.RData",sep=""))
  } else {
    load(step1filename) 
  }
  
  ### Run Stage 2
  cat("Stage 2 of 4: approximating as mixtures of Gaussians...\n")
  step2filename <- paste(gsub(" ","_",core.name),"_step2.RData",sep="")
  if(!file.exists(step2filename)) {
    if(parallel) {
      step2 <- BclimMixPar(step1)
    } else {
      step2 <- BclimMixSer(step1)
    }
    if(save.as.you.go) save(step2,file=paste(gsub(" ","_",core.name),"_step2.RData",sep=""))
  } else {
    load(step2filename)
  }
    
  ### Run Stage 3
  cat("Stage 3 of 4: estimating parameters...\n")
  step2$Chronsfile <- chron.loc
  step3filename <- paste(gsub(" ","_",core.name),"_step3.RData",sep="")
  if(!file.exists(step3filename)) {
    step3 <- BclimMCMC(step2,chron.loc,nchron=nchrons)
    if(save.as.you.go) save(step3,file=paste(gsub(" ","_",core.name),"_step3.RData",sep=""))
  } else {
    load(step3filename)
  }

  ### Run Stage 4
  cat("Stage 4 of 4: interpolating onto grid...\n")
  step4filename <- paste(gsub(" ","_",core.name),"_step4.RData",sep="")
  if(!file.exists(step4filename)) {
    step4 <- BclimInterp(step2,step3,time.grid)
    if(save.as.you.go) save(step4,file=paste(gsub(" ","_",core.name),"_step4.RData",sep=""))
  } else {
    load(step4filename) 
  }
  
  # Output results
  cat("Compiling results...\n")
  results <- BclimCompile(step1,step2,step3,step4,core.name=core.name)
  cat("Done!\n")
  return(results)
}
