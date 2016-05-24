fitting_gff <- function(positions,gff){

AZ        <- gff[,4:5] # Start-End

RUECK     <- rep(NA,2)
erg       <- apply(AZ,1,function(x){
 
  #region  <- x[1]:x[2]
  #from_to <- is.element(positions,region)
  #from_to <- which(from_to)
  from_to   <- (positions>=x[1]) & (positions<=x[2])
  from_to   <- which(from_to)
  if(length(from_to)==0){return(c(NA,NA))}
  if(length(from_to)==1){return(c(NA,NA))}
  RUECK[1] <- from_to[1]
  RUECK[2] <- from_to[length(from_to)]
  return(RUECK)

})

erg       <- t(erg)
gff[,4:5] <- erg

na.ids <- which(is.na(erg),arr.ind=TRUE)
na.ids <- unique(na.ids[,1])

if(length(na.ids)!=0){
  gff <- gff[-na.ids,]
}

return(gff)

}

############## FAST ##########################
GLOBAL.GFF <- new.env()

fitting_gff_fast <- function(positions,gff){

# sort gff for fast computation
# vec        <- gff[,4]
# names(vec) <- 1:length(vec)
# vec        <- sort(vec)
# sort.ids   <- as.numeric(names(vec))
# gff        <- gff[sort.ids,]
#-----------------------------#

AZ        <- gff[,4:5] # Start-End

MOD       <- .Call("fittingGFFC", as.matrix(AZ),positions)

gff[,4:5] <- MOD

na.ids <- which(MOD==-1,arr.ind=TRUE)
na.ids <- unique(na.ids[,1])


if(length(na.ids)!=0){
  gff                   <- gff[-na.ids,]
  # save the length of the region ! Important for SNP Data
  LENGTH                <- AZ[-na.ids,]
  SIZE                  <- (LENGTH[,2] - LENGTH[,1]) + 1
  GLOBAL.GFF$GFF        <- SIZE
  GLOBAL.GFF$start      <- LENGTH[,1]
  GLOBAL.GFF$end        <- LENGTH[,2]
  #GLOBAL.GFF$na.ids     <- na.ids
}else{
  LENGTH                <- AZ
  SIZE                  <- (LENGTH[,2] - LENGTH[,1]) + 1
  GLOBAL.GFF$GFF        <- SIZE
  GLOBAL.GFF$start      <- LENGTH[,1]
  GLOBAL.GFF$end        <- LENGTH[,2]
}

return(gff)

}








