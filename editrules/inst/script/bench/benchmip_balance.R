library(editrules)
FILE ="benchmip_balance.txt"

#' Create account balance
generate_balance <- function(nvar=15){
  i <- seq_len(floor((nvar)/2))
  edits <- c( "x1 >= 0"
            , paste0("x",i,"==", "x", 2*i, "+", "x", 2*i + 1)
            )
  
  vars <- paste0("x",seq_len(nvar))
  
  if (nvar>1){
    #edits <- c(edits, paste0(head(vars,-1), ">=", tail(vars,-1)))
    edits <- c(edits, paste0("x1 >= ", tail(vars,-1)))
  }
  
  E <- editmatrix(edits)
  if (nvar %% 2 == 0){
    E[,-(nvar+1)]
  } else {
    E
  }
}

generate_data <- function(E, nerrors=0){
  vars <- getVars(E)
  n <- length(vars)
  
  x <- sapply(vars, function(v) 0)
  x1 <- x
  x1[seq_len(nerrors)] <- -1
  
  x2 <- x
  x2[1+n - seq_len(nerrors)] <- -1  

  x3 <- x
  x3[round((n-nerrors)/2) + seq_len(nerrors)] <- -1  
  as.data.frame(rbind(x1,x2,x3))
}


bench <- function(nvars = 10, nerrors=10, method="bb"){
  
  init <- !file.exists(FILE)
  txt <- file(FILE, "at")
  on.exit(close(txt))
  
  errorloc <- c("begin", "end", "middle")
  
  if (nerrors > nvars) stop("nvars cannot be less than nerrors")
  for (nvar in seq_len(nvars)){
    for (ne in seq(1, min(nerrors, nvar))){
        try({
        E <- generate_balance(nvar)
        data <- generate_data(E, ne)
        cat("\r nvar=", nvar, " ne=", ne, " method=", method)
        le <- localizeErrors(E, data, method=method)
        rpt <- cbind(method=method, nvar=nvar, nerrors=ne, errorloc=errorloc, le$status)
        
        write.table(rpt, file=txt, col.names=init, row.names=FALSE)
        init <- FALSE
        flush(txt)
      })
      gc()
    }
  }
}

## quick testing
start <- function(){
  if (file.exists(FILE)) file.remove(FILE)
  
  bench(100,10, method="mip")
  bench(50,10, method="bb")
  
  dat <- read.table(FILE, header=TRUE)
  library(ggplot2)
  qplot(data=dat, y=elapsed, x=nvar, color=method, facets=nerrors~method, shape=errorloc, geom=c("point", "line")) + geom_jitter()
  ggsave("benchmip_balance.png")
  sdat <- subset(dat, errorloc=="middle")
  qplot(data=sdat, x=nvar, y=elapsed, color=method, group=nerrors, geom="line", 
        facets=~method, ylim=c(0,150))
}
#View(dat)
# n <- 4
# (E <- generate_balance(n))
# (dat <- generate_data(E,1))
# localizeErrors(E, dat, method="mip")