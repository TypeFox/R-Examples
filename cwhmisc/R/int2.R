str2dig <- function( str ) {
    as.numeric( unlist(strsplit(str,NULL)) )
}

int <- function(x) {
  as.integer(trunc(x))
}

intToASCII <- function(i) {
  ASCII[i %% 256];
}

intToBase <- function(i,base=2) {
  stopifnot(2 <= base & base <= 16)
  res <- HexDig[i %% base + 1]
  i <- i - i %% base
  while (i > 0) {
    i <- i %/% base
    res <- paste(HexDig[i %% base + 1],res,sep="")
    i <- i - i %% base
  }
  res
}

intToOct <- function(i) {
  intToBase(i,8)
}

intToHex <- function(i) {
  intToBase(i,16)
}
