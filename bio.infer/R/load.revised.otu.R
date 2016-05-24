"load.revised.otu" <-
function(bcnt.otu, otufname = "sum.otu.txt") {

  otutab <- read.delim(otufname)

  df1 <- merge(bcnt.otu, otutab, by = "TNAME")

  names0 <- names(bcnt.otu)

  bcnt2 <- df1[, c(names0[1:4], "otufin2")]
  names(bcnt2) <- c(names0[1:4], "OTU")

  return(bcnt2)
}

