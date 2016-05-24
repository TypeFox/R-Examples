# Cytoscape
edgecutoff <- function(dat, col) {
  a <- boxplot(col)
  if (length(a$out) > 1) {
    return(subset(dat, col > (mean(col) + sd(col))))
  }
  else {
    return(subset(dat, col > (median(col) + IQR(col))))
  }
}

edgecount <- function(dat, col1, col2) {
  ex <- dat[!mapply(grepl, col1, col2),]
  r <- reshape(ex, varying=list(colnames(dat)), direction='long',
               v.names='value', times=colnames(dat),timevar='colname')
  freq <- as.data.frame(table(r$value))
  freq <- subset(freq, freq$Freq != 0)
  return(freq)
}

dist <- function(dat, small, large) {
  q<-quantile(dat, c(small, large))
  q<-as.data.frame(q)
  return(subset(dat, dat > q$q[2] | dat < q$q[1]))
}

cytosub <- function(dat, col1, col2, text) {
  s <- substitute(text)
  f<-as.character(s)
  return(subset(dat, !grepl(f, col1) & !grepl(f, col2)))
}


# Direct from GWAS

syncheck <- function(dat, chrom, pos, col1, col2) {
  return(subset(dat, (col1 == chrom & col2 == pos)))
}

candpull <- function(setnum, setname, traitnum, traitname, threshold) {
  for (i in 1:setnum)
  {
    for (j in 1:traitnum)
    {
      x <- get(paste0(setname,i,sep=""))
      print(x[x[,paste0(traitname,j, sep="")]< threshold, 1:2])
    }
  }
}

traitcombos <- function(dat, ID){
  colnames <- colnames(dat, do.NULL = TRUE, prefix = "col")
  colnames <- colnames[which(colnames!=ID)]
  colnames2 <- colnames
  combos <- expand.grid(colnames, colnames2)
  combos <- combos[combos$Var1!=combos$Var2,]
  combos$subtract <- paste(combos$Var1, "-", combos$Var2)
  combos$add <- paste(combos$Var1, "+", combos$Var2)
  combos$multiply <- paste(combos$Var1, "*", combos$Var2)
  combos$divide <- paste(combos$Var1, "/", combos$Var2)  
  combos <- combos[,-c(1:2)]
  cvec <- as.vector(t(combos))
  dat[, cvec] <- NA
  dat[cvec] <- lapply(names(dat[cvec]), function(x) {
    with(dat, eval(parse(text = x)))
  })
  return(dat)
}