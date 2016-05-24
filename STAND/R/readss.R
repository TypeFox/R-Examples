readss <-
function(fn,L,bcol=NA,rto=5,pstat=NA,reverse=FALSE,p=0.95,gam=0.95,comma=FALSE) {
#  Read data from fn.txt or fn.csv and calculate all summary statistics
#   using IH.summary().   Output results to a csv file.
# USAGE: readss("fn",L=2,bcol=3,comma=FALSE)
# ARGUMENTS: fn in double quotes is file name without extension
#            L is specified limit value
#            comma is FALSE (for .csv file) or TRUE (for .txt file)

if(comma){ din<- read.csv( paste(fn,"csv",sep=".") ,T ) }
else{ din <- read.table( paste(fn,"txt",sep=".") ,T ) }

if( !is.na(bcol) ){ bv <- din[,bcol]
  if(!is.factor(bv)) bv <- factor(bv)
  bvc <- tapply(din[,2],bv,sum)
 # print(bvc)
  if( min(bvc) <3 ){print(bvc[ bvc <3] )
    stop("Number of Non-Detects less than 3 for Levels of BY variable Above")}
}
stat <- IH.summary(din,L,p,gam,bcol)

stat <- data.frame(stat)

    if( !is.na(pstat[1]) )  stat <- stat[pstat,]  #  select pstats

# If rev = TRUE  reverse rows and cols 
if(reverse) stat <- t(stat)

 stat <- round(stat,rto)
nout <- paste(fn,"out.csv",sep="")
write.csv(stat,nout)
#  return input data file using invisible
invisible(din)

}

