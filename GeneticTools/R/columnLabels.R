columnLabels <- function(x){
  alleles <- sort(unique(unlist(strsplit(unique(x),""))))

  # Only missing values:
  if((alleles=="X")&&(length(alleles)==1)){
    alleles[1] <- "A"
    alleles[2] <- "B"
    alleles[3] <- "X"
  }
  # Momomorph
  if((alleles!="X")&&(length(alleles)==1)){
    alleles[2] <- "B"
    alleles[3] <- "X"
  }
  if(length(alleles)==2){
    if(is.element("X",alleles)){
      alleles[2] <- "B"
      alleles[3] <- "X"
    } else {
      alleles[3] <- "X"
    }
  }
  genotypes <- c()
  hetOpt <- c(paste(alleles[1],alleles[2],sep=""),paste(alleles[2],alleles[1],sep=""))
  takeThis <- is.element(hetOpt,x)
  if(sum(takeThis)==0) takeThis <- 1
  genotypes[1] <- paste(alleles[1],alleles[1],sep="")
  genotypes[2] <- hetOpt[takeThis]
  genotypes[3] <- paste(alleles[2],alleles[2],sep="")
  genotypes[4] <- "XX"
  genotypes
} 
