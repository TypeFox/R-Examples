transform_incidences <-
function(x, from = c("PO","SO","01","-1+1"),
         to = c("PO","SO","01","-1+1"))
{
  from <- match.arg(from)
  to <- match.arg(to)

  ## check endorelation
  if ((length(dim(x)) != 2L) || (length(unique(dimnames(x))) != 1L))
    stop("Transformations only available for endorelations.")

  ## convert 'from' to canonical 'PO' representation
  SO_2_PO <- function(z) {
    tmp <- z + (z == t(z))
    tmp[is.na(tmp)] <- 0
    tmp
  }
  F01_2_PO <- function(z) {
    z[z == 0.5] <- NA
    SO_2_PO(z)
  }
  ret <- switch(from,
                PO = x,
                SO = SO_2_PO(x),
                "01" = F01_2_PO(x),
                "-1+1" = F01_2_PO((1 + x) / 2)
                )
  
  ## convert PO to target
  PO_2_SO <- function(z) {
    z[!z & (z == t(z))] <- NA
    z[(z == 1) & (z == t(z))] <- 0
    z
  }
  PO_2_01 <- function(z) {
    z[!z & (z == t(z))] <- 0.5
    z[(z == 1) & (z == t(z))] <- 0
    z
  }
  ret <- switch(to,
                PO = ret,
                SO = PO_2_SO(ret),
                "01" = PO_2_01(ret),
                "-1+1" = 2 * PO_2_01(ret) - 1
                )

  ret
}
