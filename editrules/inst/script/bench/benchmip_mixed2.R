library(editrules)
FILE ="benchmip_mixed2.txt"

generate_E <- function(nvar=10){
  s <- seq_len(nvar)
  num_vars <- paste0("x", s)
  edits <- "x1>=0"
  if (nvar > 1){
    edits <- c(edits, 
               paste0("if (", head(num_vars, -1)," >= 0)", tail(num_vars, -1), ">= 0"))
  }
  editset(edits)
}

error <- function(x){
  lapply(x, function(v){
    if (is.logical(v)) FALSE else -1
  })
}

generate_data <- function(E, nerrors=0){
  num_vars <- getVars(E,"num")
  
  x <- numeric(length(num_vars))
  names(x) <- num_vars
  n <- length(x)
    
  x1 <- x
  idx <- seq_len(nerrors)
  x1[idx] <- -1
  
  x2 <- x
  idx <- 1+n - seq_len(nerrors)
  x2[idx] <- -1
  
  x3 <- x
  idx <- (round(n-nerrors)/2) + seq_len(nerrors)
  x3[idx] <- -1
  
  as.data.frame(rbind(x1,x2,x3))
}


bench <- function(nvars = 10, nerrors=10, method="bb", maxduration=200){
  
  init <- !file.exists(FILE)
  txt <- file(FILE, "at")
  on.exit(close(txt))
  
  errorloc <- c("begin", "end", "middle")
  
  if (nerrors >= nvars) stop("nvars cannot be less than nerrors")
  for (nvar in seq_len(nvars)){
    E <- generate_E(nvar)
    for (ne in seq(1, min(nerrors, nvar))){
      try({
        data <- generate_data(E, ne)
        # only select middle
        data <- data[3,,drop=FALSE]
        errorloc_m <- errorloc[3]
        
        cat("\r nvar=", nvar, " ne=", ne, " method=", method)
        le <- localizeErrors(E, data, method=method, maxduration=maxduration)
        rpt <- cbind(method=method, nvar=nvar, nerrors=ne, errorloc=errorloc_m, le$status)
        
        write.table(rpt, file=txt, col.names=init, row.names=FALSE)
        
        init <- FALSE
        flush(txt)
      })
      gc()
    }
  }
}

if (file.exists(FILE)) file.remove(FILE)

bench(50,10, method="mip")
bench(50,10, method="bb")

# dat <- read.table(FILE, header=TRUE)
# library(ggplot2)
# qplot(data=dat, y=elapsed, x=nvar, color=method, facets=nerrors~method, shape=errorloc) + geom_jitter()
# ggsave("benchmip_categorical.png")


### quick testing
E <- generate_E(4)
E
generate_data(E, nerrors=2)
