.spacC <- "  "  ## extra space in indexLine and charMat

n2c <- function(x, symm = FALSE){
  n2c0 <- function(y) {
    ifelse(y >= 10e10, "X", ifelse(y >= 1.0, 
      as.character(trunc(log(y, 10))), ifelse(y >= 0.9, "&",
      ifelse(y >= 0.8, "%", ifelse(y >= 0.7, "#", 
      ifelse(y >= 0.6, "*", ifelse(y >= 0.5, "=", ifelse(y >= 0.4, "+", 
      ifelse(y >= 0.3, "-", ifelse(y >= 0.2, ":", ifelse(y >= 0.1, ",", 
      ifelse(y >= 0.05, ".", " ")))))))))))) }
  if (is.matrix(x))
    y <- ifelse((row(x) < col(x)) & symm, " ", n2c0(abs(x)))
  else
    y <- n2c0(abs(x))
  attr(y,"legend") <-  paste(">=1Log","9&", "8%","7#","6*","5=","4+","3-","2:","1,","05.")
  y
}  ## n2c

indexLine <- function(n) {
  L <- c(rep(".",n),.spacC)  # extra space needed for charMat
  if (n>=5) {
    T <- seq(n%/%5)
    L[5*T] <- ";"
    if (n>=10) {
      T <- seq(n%/%10)
      L[10*T] <- T
    }
  }
  paste(L,collapse="")
} ## indexLine

n2cCompact <- function(x, symm=FALSE) {
  nB <- n2c(x, symm=symm)
  cc <- indexLine(ncol(x))
  c(cc,paste(apply(nB,1,paste,collapse=""),.spacC,if (is.null(rownames(x))) seq(nrow(x)) else abbreviate(rownames(x),minlength = 10),sep=""),cc,paste("legend: ",attr(nB,"legend"),sep=""))
} ## n2cCompact

charMat <- function(cc) { ## lines of type n2cCompact
  rows <- length(cc)-3    ## strip lines 1 and -2, -1
  colP <- nchar(cc[1])
  colS <- colP-nchar(.spacC)     ## without extra space in indexLine
  cc1  <- substring(cc,1,colP)  ## proper lines of character matrix
  cc2  <- substring(cc,colP+1)  ## column numbers
  cm1  <- rev(rev(cc1[-1])[-c(1,2)])  ## without indexLines and legend
  cm2  <- rev(rev(cc2[-1])[-c(1,2)])  ## useful line numbers
  U <- rep(colP,rows)  ## 
  x <- unlist(lapply(U,seq)) 
  y <- rep(seq(U),U)
  tx <- unlist(strsplit(cm1,split=""))
  tx[colP*seq(rows)] <- cm2
  list(x=x, y=y, tx=tx) 
} ## charMat

explainLegend <- function() {
  print("Show n>=1 as log(n)")
  print("     n<1 according to its first decimal 987654321 as &%#*=+-:,")
##  print("     n<1 according to first decimal as '&' '%' '#' '*' '=' '+' '-' ':' ','")
  print("    n>.05 as '.'")
}
