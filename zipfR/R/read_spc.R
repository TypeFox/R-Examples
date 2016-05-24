read.spc <- function (file)
{
  tmp <- read.delim(auto.gzfile(file))
  vars <- colnames(tmp)
  if (!("m" %in% vars)) stop("required column 'm' missing from .spc file ", file)
  if (!("Vm" %in% vars)) stop("required column 'Vm' missing from .spc file ", file)
  variances <- ("VVm" %in% vars)

  m <- tmp$m
  Vm <- tmp$Vm
  if (variances) VVm <- tmp$VVm
  
  if (! is.integer(m)) stop("frequency classes 'm' must be integers in .spc file ", file)
  if (any(m < 1)) stop("frequency classes 'm' must be >= 1 in .spc file ", file)
  if ( (! is.numeric(Vm)) || any(Vm < 0) )
    stop("class sizes 'Vm' must be non-negative numbers in .spc file ", file)
  if (variances && (! is.numeric(Vm)) || any(Vm <= 0) )
    stop("variances 'VVm' must be positive numbers in .spc file ", file)
  if (variances)
    warning("variance of expected vocabulary size cannot be reconstructed from disk file")
    
  if (variances) spc(m=m, Vm=Vm, VVm=VVm) else spc(m=m, Vm=Vm)
}
