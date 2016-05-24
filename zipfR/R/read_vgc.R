read.vgc <- function (file)
{
  tmp <- read.delim(auto.gzfile(file))
  vars <- colnames(tmp)
  if (!("N" %in% vars)) stop("required column 'N' missing from .vgc file ", file)
  if (!("V" %in% vars)) stop("required column 'V' missing from .vgc file ", file)

  m <- 1:9
  idx <- sapply(m, function (x) if (paste("V", x, sep="") %in% vars) TRUE else FALSE)
  m.max <- sum(idx)            # = 0 if no spectrum elements are given in file
  if (m.max > 0) {
    if (any(m[!idx] <= m.max))
      stop("spectrum elements must form uninterrupted sequence 'V1', 'V2', ...")
  }

  variances <- ("VV" %in% vars)
  idx.v <- sapply(m, function (x) if (paste("VV", x, sep="") %in% vars) TRUE else FALSE)
  m.max.v <- sum(idx.v)
  if (m.max.v > 0) {
    if (any(m[!idx.v] <= m.max.v))
      stop("spectrum variances must form uninterrupted sequence 'VV1', 'VV2', ...")
  }

  if (!variances && m.max.v > 0)
    stop("column 'VV' must also be given if there are variances 'VV1', etc.")
  if (variances && (m.max.v != m.max))
    stop("variances 'VVm' must be given for exactly the same frequency classes as expectations 'Vm'")
  
  Vm.list <- if (m.max > 0) tmp[ paste("V", 1:m.max, sep="") ] else NULL
  ## single index -> data.frame interpreted as list -> always returns list, even for m.max=1
  if (variances) {
    VVm.list <- if (m.max > 0) tmp[ paste("VV", 1:m.max, sep="") ] else NULL
    vgc(N=tmp$N, V=tmp$V, VV=tmp$VV, Vm=Vm.list, VVm=VVm.list)
  }
  else {
    vgc(N=tmp$N, V=tmp$V, Vm=Vm.list)
  }
}
