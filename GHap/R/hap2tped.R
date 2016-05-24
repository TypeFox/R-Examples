#Function: ghap.hap2tped
#License: GPLv3 or later
#Modification date: 5 Feb 2016
#Written by: Yuri Tani Utsunomiya & Marco Milanesi
#Contact: ytutsunomiya@gmail.com, marco.milanesi.mm@gmail.com
#Description: Convert haplotype allele counts to tped

ghap.hap2tped<-function(
  infile,
  batchsize=500,
  outfile,
  verbose=TRUE
){
  
  #Check if output will ovewrite existing files
  tfam <- paste(outfile,"tfam",sep=".")
  tped <- paste(outfile,"tped",sep=".")
  tref <- paste(outfile,"tref",sep=".")
  if(file.exists(tfam) == TRUE){
    stop("The tfam file already exists!")
  }
  if(file.exists(tped) == TRUE){
    stop("The tped file already exists!")
  }
  if(file.exists(tref) == TRUE){
    stop("The tref file already exists!")
  }
  
  #Check input file
  hapsamples <- paste(infile,"hapsamples",sep=".")
  hapalleles <- paste(infile,"hapalleles",sep=".")
  hapgenotypes <- paste(infile,"hapgenotypes",sep=".")
  if(file.exists(hapsamples) == FALSE){
    stop(paste("Could not find",hapsamples,"\n"))   
  }
  if(file.exists(hapalleles) == FALSE){
    stop(paste("Could not find",hapalleles,"\n"))   
  }
  if(file.exists(hapgenotypes) == FALSE){
    stop(paste("Could not find",hapgenotypes,"\n"))
  }else{
    hapgenotypes  <- file(hapgenotypes, open = "r") 
  }
  
  #Load hapsamples file
  if(verbose == TRUE){
    cat("Creating .tfam file...")
  }
  hapsamples <- read.table(
    file = hapsamples,
    header = FALSE,
    sep = " ",
    colClasses = "character",
    stringsAsFactors = FALSE
  )
  write.table(x = cbind(hapsamples,"0 0 0 -9"), file = tfam, quote = FALSE, sep=" ", row.names = FALSE, col.names = FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  
  #Load hapalleles file
  if(verbose == TRUE){
    cat("Creating .tref file...", sep="")
  }
  hapalleles <- read.table(
    file = hapalleles,
    header = FALSE,
    sep = " ",
    colClasses = "character",
    stringsAsFactors = FALSE,
    col.names = c("BLOCK","CHR","BP1","BP2","ALLELE")
  )
  hapalleles$NAME <- apply(hapalleles[,c("BLOCK","BP1","BP2","ALLELE")],MARGIN = 1, paste, sep="", collapse="_")
  write.table(x = cbind(hapalleles$NAME,"H"), file = tref, quote = FALSE, sep=" ", row.names = FALSE, col.names = FALSE)
  if(verbose == TRUE){
    cat("Done.\n")
  }
  hapalleles$BP <- as.character(round((as.numeric(hapalleles$BP1)+as.numeric(hapalleles$BP2))/2,digits=0))
  hapalleles$CM <- "0"
  hapalleles <- hapalleles[,c("CHR","NAME","CM","BP")]
  hapalleles <- apply(hapalleles,MARGIN = 1,paste,collapse=" ")
  
  #Read hapgenotypes file
  if(verbose == TRUE){
    cat("Creating .tped file for", length(hapalleles), "haplotype alleles.\n")
  }
  i<-1
  while(length(line <- readLines(con=hapgenotypes, n = batchsize, warn = FALSE)) > 0){
    line <- gsub(pattern = "0", replacement = "N N", x = line)
    line <- gsub(pattern = "1", replacement = "N H", x = line)
    line <- gsub(pattern = "2", replacement = "H H", x = line)
    outline <- paste(hapalleles[i:((i+length(line))-1)],line,sep=" ")
    i <- i + length(line)
    
    tped.con  <- file(tped, open = "a")
    writeLines(text = outline, con = tped.con)
    if(verbose == TRUE){
      cat(i-1, "haplotype alleles converted.\r")
    }
    close(tped.con)
  }
  
  #Close connections
  close(hapgenotypes)
  
}
