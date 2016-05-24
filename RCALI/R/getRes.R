#+++++++++++++++++++++++++++++++++++++++++++++++
# FUNCTION:
# Read a result file from CaliFloPP
# RETURN:
# a labelled data.frame
# 18/02/2010, revision 10/03/2010
# ++++++++++++++++++++++++++++++++++++++++++++++
getRes <- function(ficres) {

  # Read  nfunc and method:
  lu <- scan(ficres, nmax=8, what=character(0), quiet=TRUE)
  nfunc <- as.integer(lu[6])
  method <- lu[8]

  # Read by "read.table" is slower than "scan"
# tps <- system.time(res <- read.table(ficres,  skip=1))
#  print(tps)
  # Endd the number of columns to read
  ncollu <- 4 +nfunc
  if (method == "cubature")
    ncollu <- ncollu+ nfunc*5
  else
     ncollu <- ncollu+ nfunc*2
  
 res <- as.data.frame(matrix(data=scan(file=ficres, skip=1), ncol=ncollu, byrow=TRUE))

output.format <- 0
if (dim(res)[2] > (2+nfunc))
  output.format <- 1
  

# Labels of the infos which are similar to the both methods 
# and to the both output_format
lab <- c("poly1", "poly2")
if (nfunc >1) {
  moy <- paste("mean.flow.f",  seq(1,nfunc), sep="")
  moyaire <- paste(moy, "/area2", sep="")
}
else {
  moy <- "mean.flow"
  moyaire <- paste(moy, "/area2", sep="")
}
aire <-s  <- NULL
if (output.format > 0) {
aire <- c("area1", "area2")
switch(method,
       cubature = {
         if (nfunc >1) {
           s <- paste(c( "mean.flow.f",
                        "conf.int.lower.f", "conf.int.upper.f",
                        "abs.err.f", "n.eval.f"),
                      rep(seq(1,nfunc), rep(5, nfunc)), sep="")
         }
         else
           s <- c("mean.flow",
                  "conf.int.lower", "conf.int.upper", "abs.err", "n.eval")
       },
       grid= {
        if (nfunc >1) {
          s <- paste( c( "mean.flow.f","std.f"),
                     rep(seq(1,nfunc), rep(2,nfunc)), sep="")
        }
        else
          s <-  c("mean.flow","std")
       },
        stop(paste("Invalid name for a method:", method))

 ) # end switch
} # end output.format
dimnames(res)[[2]] <- c(lab, moyaire, aire, s)
return(res)
} # end getRes

# ++++++++++++++++++++++++++++++++++++++++++++++

