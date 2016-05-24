##
## Normalize PD values in NIfTI files.
##
niinorm <- 
function(srcdir=tempfile(), filename="data_V1",
         savedir=tempdir())
{
  normvf <- function(x) { norm(matrix(x,length(x),1),"f") }
  cat("Reading data ...\n")
  options(niftiAuditTrail = FALSE)
  filepath <- file.path(srcdir, filename)
  FEvol <- readNIfTI(filepath, reorient=FALSE)
  FE <- FEvol@.Data
  d <- dim(FE)
  cat("Data normalization ...\n")
  nnx <- apply(FE, c(1,2,3), normvf)
  nnx <- drop(nnx)
  for(i in 1:d[4])
    FE[,,,i] <- FE[,,,i]/nnx
  zi <- which(is.nan(FE))
  FE[zi] <- 0
  fo <- strsplit(filename,".nii")[[1]][1]
  fo <- paste(savedir,"/",fo,"n",sep="")
  writeNIfTI(FE, filename=fo, verbose=TRUE)
  cat("wrote",fo,"\n")
}

