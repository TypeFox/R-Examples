#########################################################
#### AUTHOR:     Arnost Komarek                      ####
####             (2005)                              ####
####                                                 ####
#### FILE:       files.Gspline.R                     ####
####                                                 ####
#### FUNCTIONS:  clean.Gspline                       ####
####             write.headers.Gspline               ####
#########################################################

## Functions to write headers for simulated G-splines
## and to clean files with simulated G-splines

### ======================================
### clean.Gspline
### ======================================
clean.Gspline <- function(dir, label, care.of.y=TRUE){
  FILES <- dir(dir)
  
  if (paste("mixmoment", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/mixmoment", label, ".sim", sep = ""))
  if (paste("mweight", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/mweight", label, ".sim", sep = ""))
  if (paste("mlogweight", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/mlogweight", label, ".sim", sep = ""))
  if (paste("mmean", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/mmean", label, ".sim", sep = ""))
  if (paste("gspline", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/gspline", label, ".sim", sep = ""))
  if (paste("lambda", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/lambda", label, ".sim", sep = ""))
  if ((paste("Y", label, ".sim", sep="") %in% FILES) & care.of.y) file.remove(paste(dir, "/Y", label, ".sim", sep = ""))
  if (paste("logposter", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/logposter", label, ".sim", sep = ""))    
  if (paste("r", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/r", label, ".sim", sep = ""))
}

### ======================================
### write.headers.Gspline
### ======================================
write.headers.Gspline <- function(dir, dim, nP, label, gparmi, store.a, store.y, store.r, care.of.y=TRUE){
  FILES <- dir(dir)
  
  sink(paste(dir, "/mixmoment", label, ".sim", sep = ""), append = FALSE)
  mname <- paste("Mean.", 1:dim, "   ", sep="")
  D <- diag(dim)
  rows <- row(D)[lower.tri(row(D), diag = TRUE)]
  cols <- col(D)[lower.tri(col(D), diag = TRUE)]            
  dname <- paste("D.", rows, ".", cols, sep = "")
  cat("k", mname, dname, "\n", sep = "   "); sink()

  total.length <- ifelse(dim == 1, 2*gparmi["K1"] + 1,
                                   (2*gparmi["K1"] + 1)*(2*gparmi["K2"] + 1))
  sink(paste(dir, "/mweight", label, ".sim", sep = ""), append = FALSE)
  cat(paste("w", 1:min(total.length, 9), sep = ""), sep = "             ")
  if (total.length >= 10){
    cat("             ")
    cat(paste("w", 10:total.length, sep = ""), "\n", sep = "            ")
  }     
  else{
    cat("\n")
  }     
  sink()

  if (store.a){
    sink(paste(dir, "/mlogweight", label, ".sim", sep = ""), append = FALSE)
    cat(paste("a", 1:min(total.length, 9), sep = ""), sep = "             ")
    if (total.length >= 10){
      cat("             ")
      cat(paste("a", 10:total.length, sep = ""), "\n", sep = "            ")
    }     
    else{
      cat("\n")
    }     
    sink()
  }
  else
    if (paste("mlogweight", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/mlogweight", label, ".sim", sep = ""))
   
  ind1 <- rep(1:total.length, rep(dim, total.length))
  ind2 <- rep(1:dim, total.length)
  sink(paste(dir, "/mmean", label, ".sim", sep = ""), append = FALSE)
  cat(paste("mu", ind1, ".", ind2, sep = ""), "\n", sep = "      "); sink()
 
  sink(paste(dir, "/gspline", label, ".sim", sep = ""), append = FALSE)
  gname <- paste(" gamma", 1:dim, sep = "")
  sname <- paste(" sigma", 1:dim, sep = "")   
  dname <- paste(" delta", 1:dim, sep = "")
  intcptname <- paste(" intercept", 1:dim, sep = "")   
  scname <- paste(" scale", 1:dim, sep = "")   
  cat(gname, sname, dname, intcptname, scname, "\n", sep = "  "); sink()

  sink(paste(dir, "/lambda", label, ".sim", sep = ""), append = FALSE)
  lname <- if(gparmi["equal.lambda"]) "lambda" else paste("lambda", 1:dim, sep = "")
  cat(lname, "\n", sep = "  "); sink()

  if (care.of.y){
    ind1 <- rep(1:nP, rep(dim, nP))
    ind2 <- rep(1:dim, nP)
    if (store.y){
      sink(paste(dir, "/Y", label, ".sim", sep = ""), append = FALSE)
      cat(paste("Y.", ind1, ".", ind2, sep = ""), "\n", sep = "      ")
      sink()
    }
    else{
      if (paste("Y", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/Y", label, ".sim", sep = ""))
    }
  }

  if (store.r){
    sink(paste(dir, "/r", label, ".sim", sep = ""), append = FALSE)
    cat(paste("r.", ind1, ".", ind2, sep = ""), "\n", sep = "      ")
    sink()
  }
  else{
    if (paste("r", label, ".sim", sep="") %in% FILES) file.remove(paste(dir, "/r", label, ".sim", sep = ""))
  }
    
  pname <- if (gparmi["equal.lambda"]) "penalty" else paste("penalty", 1:dim, sep = "")
  pname <- c("loglik      ", pname, "      logprw")
  sink(paste(dir, "/logposter", label, ".sim", sep = ""), append = FALSE)
  cat(pname, "\n", sep = "  "); sink()
}
