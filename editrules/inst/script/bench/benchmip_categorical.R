library(editrules)
FILE ="benchmip_categorical.txt"

#' Create account balance
generate_E <- function(nvar=10){
  vars <- paste0("A", seq_len(nvar))
  n <- nvar-1
  
  # if all TRUE no errors. FALSE values introduce errors
  edits <- c( paste0(vars, " %in% c(TRUE,FALSE)")
            , paste0("if (",head(vars, -1), ")", tail(vars, -1))
            )
  edits <- c(edits, "A1 == TRUE")
  editarray(edits)
}

generate_data <- function(E, nerrors=0){
  vars <- getVars(E)
  n <- length(vars)
  
  x <- sapply(vars, function(v) TRUE)
  x1 <- x
  x1[seq_len(nerrors)] <- FALSE
  
  x2 <- x
  x2[1+n - seq_len(nerrors)] <- FALSE
  
  x3 <- x
  x3[round((n-nerrors)/2) + seq_len(nerrors)] <- FALSE
  as.data.frame(rbind(x1,x2,x3))
}


bench <- function(nvars = 10, nerrors=9, method="bb"){
  
  init <- !file.exists(FILE)
  txt <- file(FILE, "at")
  on.exit(close(txt))
  
  errorloc <- c("begin", "end", "middle")
  
  if (nerrors >= nvars) stop("nvars cannot be less than nerrors")
  for (ne in seq_len(nerrors)){
    for (nvar in seq(from=ne+1, to=nvars)){
      try({
        E <- generate_E(nvar)
        data <- generate_data(E, ne)
        cat("\r nvar=", nvar, " ne=", ne, " method=", method)
        le <- localizeErrors(E, data, method=method)
        rpt <- cbind(method=method, nvar=nvar, nerrors=ne, errorloc=errorloc, le$status)
        write.table(rpt, file=txt, col.names=init, row.names=FALSE)
        init <- FALSE
        flush(txt)
      })
      gc()
      #print(rpt)
    }
  }
}

  if (file.exists(FILE)) file.remove(FILE)

  bench(50,10, method="mip")
  bench(50,10, method="bb")

  dat <- read.table(FILE, header=TRUE)
  library(ggplot2)
  qplot(data=dat, y=elapsed, x=nvar, color=method, facets=nerrors~method, shape=errorloc) + geom_jitter()
  ggsave("benchmip_categorical.png")
# 


### quick testing
# (E <- generate_E())
# (dat <- generate_data(E, nerrors=3))
# localizeErrors(E,dat, method="mip")
# editrules:::errorLocalizer.mip(E,dat[1,])
