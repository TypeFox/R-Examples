subDose <- function(file=NULL, vmmk=NULL, out=NULL, removeInsertions=TRUE, verbose=TRUE){
  # Avoid E-Notation 
    options("scipen"=100, "digits"=4)
  # Basic input check
    if(is.null(file)) stop("No *.dose file given.")
    if(is.null(vmmk)) stop("No Variant-Map-Master-Key file given.")
    if(is.null(out)) out <- paste("Allele-",file,sep="")

  # Read in the variant master key:
    master <- read.table(vmmk, stringsAsFactors=FALSE)
  # Read in the dose file
    dose <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
  # Output
    if(verbose) cat("I read the dose file and the VMMK.",date(),"\n")
  # Check if the order is correct, otherwise stop the script!
    if(sum(dose$marker==master$V1)!=nrow(dose)) stop("The order does not match between the files")
  # Now change the dose table
    dose$alleleA <- master$V3
    dose$alleleB <- master$V4
  # Remove the Insertions
    if(removeInsertions) dose <- dose[-which(master$V7=="I"),]
  # write out the new Dose file
    write.table(dose,out, col.names=TRUE, row.names=FALSE, quote=FALSE)
}