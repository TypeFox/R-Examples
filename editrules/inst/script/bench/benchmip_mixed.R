library(editrules)
FILE ="benchmip_mixed.txt"

generate_E <- function(nvar=10){
  if (nvar < 1) stop("nvar needs to be bigger than 3")
  s <- seq_len(nvar)
  n_num <- ceiling(nvar/2)
  var_num <- head(s, n_num)
  var_cat <- tail(s, -n_num)
  n_cat <- length(var_cat)
  
  var_num <- paste0("x", var_num)
  if (n_cat) var_cat <- paste0("v", var_cat-n_num) else character()
  
  if (length(var_num) > 1){
    nsum <- paste(tail(var_num, -1), collapse="+")
    edits <- paste0("x1 == ", nsum)
    edits <- c(edits, paste0(head(var_num, -1)," >= ", tail(var_num, -1)))
  } else {
    edits <- "x1 == 0"
  }
  
  if (n_cat){
    edits <- c( paste0(tail(var_num, 1), ">= 0") 
              , edits
              , paste0(var_cat, " %in% c(TRUE,FALSE)")
              , paste0("if (!", var_cat, ") ", head(var_num,n_cat),"< 0")
              )
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
  
  x <- sapply(getVars(E), function(v){
    if (v %in% num_vars) 0 else TRUE
  }, simplify=FALSE)
  
  n <- length(x)
    
  x1 <- x
  idx <- seq_len(nerrors)
  x1[idx] <- error(x1[idx])
  
  x2 <- x
  idx <- 1+n - seq_len(nerrors)
  x2[idx] <- error(x2[idx])
  
  x3 <- x
  idx <- (round(n-nerrors)/2) + seq_len(nerrors)
  x3[idx] <- error(x3[idx])
  
  rbind(as.data.frame(x1),x2,x3)
}


bench <- function(nvars = 10, nerrors=10, method="bb", maxduration=200){
  
  init <- !file.exists(FILE)
  txt <- file(FILE, "at")
  on.exit(close(txt))
  
  errorloc <- c("begin", "end", "middle")
  
  if (nerrors >= nvars) stop("nvars cannot be less than nerrors")
  for (nvar in seq_len(nvars)){
    E <- generate_E(nvar)
    max_dur <- logical(nrow(generate_data(E, 1)))
    for (ne in seq(1, min(nerrors, nvar))){
      try({
        if (all(max_dur)) break
        data <- generate_data(E, ne)[!max_dur,,drop=FALSE]
        # only select middle
        data <- data[3,,drop=FALSE]
        errorloc_m <- errorloc[3]
        
        cat("\r nvar=", nvar, " ne=", ne, " method=", method)
        le <- localizeErrors(E, data, method=method, maxduration=maxduration)
        max_dur[!max_dur] <- le$status$maxDurationExceeded
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

dat <- read.table(FILE, header=TRUE)
library(ggplot2)
qplot(data=dat, y=elapsed, x=nvar, color=method, facets=nerrors~method, shape=errorloc) + geom_jitter()
ggsave("benchmip_categorical.png")


### quick testing
# E <- generate_E(2)
# generate_data(E, nerrors=1)
