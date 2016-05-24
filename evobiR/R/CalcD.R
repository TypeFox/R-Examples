CalcD <- function(alignment = "alignment.fasta", 
                  sig.test="N",                                                    # options are "N", "B", "J"
                  block.size = 1000,                                                # size of blocks to drop in jacknife
                  replicate=1000){
  # this function is used regardless of approach
  d.calc <- function(alignment){
    abba <- 0                                                                         #  set up my variables
    baba <- 0                                                                         #  set up my variables
    for(i in 1:ncol(alignment)){                                               #  run through all sites
      if(length(unique(alignment[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
        if(alignment[1, i] != alignment[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
          if(alignment[4, i] != alignment[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
            if(alignment[3, i] == alignment[1, i]) {baba <- baba + 1}   #  add to the count of baba sites
            if(alignment[2, i] == alignment[3, i]) {abba <- abba + 1}   #  add to the count of abba sites
          } 
        }
      }
    }
    d <- (abba - baba) / (abba + baba)   #what its all about
    results <- list()
    results[[1]] <- d
    results[[2]] <- abba
    results[[3]] <- baba
    return(results)
  }
  
  #### Test of empirical data
  alignment <- read.alignment(alignment, format = "fasta")                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  results <- d.calc(alignment.matrix)
  d <- results[[1]]
  if(is.nan(d)) d <- 0
  abba <- results[[2]]
  baba <- results[[3]]
  
## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
## THE STANDARD DEVIATION AND THUS A Z SCORE.
  if(sig.test=="B"){
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo)
    cat("\nperforming bootstrap")
    for(k in 1:replicate){
      if(k/(replicate/100) == round(k/(replicate/100))) cat(".")
      sim.matrix[1:4,1:foo] <- alignment.matrix[1:4, sample(1:foo, replace=T)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    sim.d[is.nan(sim.d)] <- 0
    z <- abs(d/sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    cat("\n\nResults from ", replicate, "bootstraps")
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
  }

## THIS SECTION WILL CALCULATE THE P-VAL BASED ON JACKKNIFING
## THE DATA IS REANALYZED WHILE DROPPING A PORTION OF THE DATA
## THE SIZE OF THE DROPPED PORTION IS DETERMINED BY THE BLOCK SIZE
## ARGUMENT THIS PROCEDURE ALLOWS US TO CALCULATE
## THE STANDARD DEVIATION AND THUS A Z SCORE. THIS APPROACH IS PARTICULARLY 
## IMPORTANT WHEN WE WILL BE USING DATA WHERE THE SNPs MAY BE IN LINKAGE WITH
## ONE ANOTHER
if(sig.test=="J"){
  
  #first lets test whether we are limited in the number of reps by block size
  max.rep <- ncol(alignment.matrix) - block.size
  if(block.size >= (ncol(alignment.matrix)/2)){
    stop(call. = F, paste("\nWith a block size of", block.size, 
              "and an alignment of", ncol(alignment.matrix), 
              "sites \nsome sites would never be included in the \nanalysis",
              "\n\nThe maximum block size is 1/2 the alignment length"))
  }
  if(max.rep < replicate){
    stop(call. = F, paste("\nWith a block size of", block.size, 
                          "and an alignment of", ncol(alignment.matrix), 
                          "sites", replicate, "replicates\nare not possible"))
  }
  
  
  if(max.rep >= replicate){
    drop.pos <- seq.int(from=1, to=(max.rep-1), length.out=replicate)
    replicate2 <- replicate
  }
  sim.d<-vector()
  foo <- ncol(alignment.matrix)
  sim.matrix<-matrix(,4,foo-block.size)
  cat("\nperforming jackknife")
  for(k in 1:replicate2){  
    if(k/2 == round(k/2)) cat(".")
    sim.matrix[1:4,1:(foo-block.size-1)] <-alignment.matrix[1:4, -drop.pos[k]:-(drop.pos[k]+block.size)]
    results <- d.calc(sim.matrix)
    sim.d[k] <- results[[1]]
  }
  sim.d[is.nan(sim.d)] <- 0
  z <- abs(d/sd(sim.d))
  new.pval <- 2 * (1 - pnorm(z))
  ## NOW WE MAKE THE OUTPUTS  
  cat("\nSites in alignment =", ncol(alignment.matrix))
  cat("\nNumber of sites with ABBA pattern =", abba)
  cat("\nNumber of sites with BABA pattern =", baba)
  cat("\nD raw statistic", d)
  cat("\nZ-score = ", z)
  cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
  cat("\nSD D statistic =", sd(sim.d))
  cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
}
  if(sig.test=="N"){
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic = ", d)
  }
}
CalcPopD <- function(alignment = "alignment.fasta"){
  ##  Now we have eqn. 2 from page 2240
  ##  input is an alignment the can take multiple sequences from each 
  ##  population of interest.  IMPORTANT MAKE SURE SEQUENCES ARE IN ORDER
  ##  P1, P2, P3, OUTGROUP!  Again we find the biallelic sites but now 
  ##  those biallelic sites need not be fixed and we will calculate frequencies
  ##  of SNP for each population.  The way the function is set up we do need to 
  ##  feed in an alignment where each sequence from a population has the same name:
  ##  pop1
  ##  AACCACAAGCCAGCTCAGCTACAG
  ##  pop1
  ##  TACAACAAGCGAGCTCAGCTACAG
  ##  pop1
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop3
  ##  TACCACAAGCCAGCTCAGCTACAG
  ##  OUTGROUP
  ##  TACCAGGAGCCAGCTCTTCTACCC
  Mode <- function(x) {                                                                      #  i need a little mode function which R is lacking ugh
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }  
  alignment<-read.alignment(alignment, format="fasta")                                       #  read in the alignment
  alignment.matrix<-matrix(,length(alignment$nam),nchar(alignment$seq[[1]])+1)               #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i,2:ncol(alignment.matrix)]<-unlist(strsplit(alignment$seq[[i]],""))    #  fill in the matrix
  }
  alignment.matrix[,1]<-alignment$nam                                                        #  get those names into our matrix row names dont work :(
  groups<-unique(alignment$nam)
  p1 <- p2 <- p3 <- p4 <- 0                                                                  #  lets just set up the variable names from the durand paper
  numerator <- denominator <- 0
  useful<-0                                                                                  #  plus some of my own
  segregating<-0                                                                             #  plus some of my own  
  seg.pos<-F                                                                                 #  plus some of my own  
  for(i in 2:ncol(alignment.matrix)){                                                        #  run through all sites
    seg.pos<-F                                                                               #  reset this switch  
    if(length(unique(alignment.matrix[,i]))==2){                                             #  unique(c(p1,p2,p3,o))==2 aka biallelic
      A <- Mode(alignment.matrix[alignment.matrix[, 1] == groups[4], i])                     #  lets treat the more common variant in the outgroup as "A"
      B <- unique(alignment.matrix[,i])[unique(alignment.matrix[, i]) != A]                  #  not purposely obfuscating... the other variant in variable "B"
      if(B %in% unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])){            #  makes sure that we have at least some indication of an ABBA/BABA pattern
        if(length(unique(alignment.matrix[alignment.matrix[, 1] %in% groups[1:2], i])) == 2){  #  makes sure that we've got some different resolutions in the ingroups
          useful <- useful + 1                                                                 #  lets just keep track of how many sites are even useful 
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[1], i])) == 2) {seg.pos<-T}#  next 5 lines are a lame way of counting sites that are segregating
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[2], i])) == 2) {seg.pos<-T}#  vs those that are fixed another words is population sampling
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])) == 2) {seg.pos<-T}#  really of any value within the data set that we are examining
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[4], i])) == 2) {seg.pos<-T}
          if(seg.pos == T){segregating <- segregating + 1} 
          #print(segregating)
          p1 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[1], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[1], i])  #  freq of A snp in first population
          p2 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[2], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[2], i])  #  freq of A snp in second population
          p3 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[3], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[3], i])  #  freq of A snp in third population
          p4 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[4], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[4], i])  #  freq of A snp in outgroup population
          #  Durands explanation of eqn 2 is lacking... as least to my feable mind!
          #  it appears to me that as written p hat is actually the frequency of SNP "B" so....
          #  snap...  vindicated my interpretation matches that found in the supplemental material of the 
          #  heliconius genome paper supplement... too cool
          p1 <- 1-p1  #convert these over from proportion A to proportion B
          p2 <- 1-p2  #convert these over from proportion A to proportion B
          p3 <- 1-p3  #convert these over from proportion A to proportion B
          p4 <- 1-p4  #convert these over from proportion A to proportion B
          numerator <- ((1 - p1) * p2 * p3 * (1 - p4)) - (p1 * (1 - p2) * p3 * (1 - p4)) + numerator              #  build up our numerator sum
          denominator <- ((1 - p1) * p2 * p3 * (1 - p4)) + (p1 * (1 - p2) * p3 * (1 - p4)) + denominator          #  build up our denominator sum
        }
      }
    }
  }
  d <- numerator / denominator    #what its all about
  
  user.result <- list()
  user.result$d.stat <- d
  user.result$pval <- "HELP"
  user.result$align.length <- ncol(alignment.matrix) - 1
  user.result$useful.sites <- useful
  user.result$seg.sites <- segregating
  print(paste("Sites in alignment =", ncol(alignment.matrix) - 1))
  print(paste("Number of sites with ABBA or BABA patterns =", useful))
  print(paste("Number of ABBA or BABA sites that are still segregating in at least one population =", segregating))  
  print(paste("D statistic =", d))
}