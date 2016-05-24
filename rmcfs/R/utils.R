###############################
#showme
###############################
showme <- function(x, size=10, show=c("tiles","head","tail","none")){
  show <- show[1]
  if(show == "tiles"){
    print(x[1:size,1:size])
    cat("\n\n")
    print(x[(nrow(x)-size):nrow(x), (ncol(x)-size):ncol(x)])
  }else if(show == "head"){
    print(head(x,size))
  }else if(show == "tail"){
    print(tail(x,size))
  }else if(show == "none"){
    #nothing
  }else
    stop(paste0("Parameter 'show' is incorrect: ", show))
  
  cat(paste0("class: '", class(x), "' size: ", nrow(x)," x ", ncol(x)))
}

###############################
#get.size.param
###############################
get.size.param <- function(size, default_value){
  if(is.na(size))
    size <- default_value
  if(is.na(size) | size <= 0){
    stop(paste0("'size' is NA or <= 0"))
  }
  return(size)
}

###############################
#my.seq
###############################
my.seq <- function(from, to, by, add.last.to = F){
  if(to <= by){
    s <- c(to)
  }else{
    s <- seq(from, to, by)
    if(add.last.to){
      if(s[length(s)]<to)
        s <- c(s,to)
    }
  }
  return (s)
}

###############################
#df.to.matrix
###############################
#  data(alizadeh)
#  d <- alizadeh
#  print(paste0("class: ", class(d), " size: ", nrow(d)," x ", ncol(d)))
#  d.matrix <- df.to.matrix(x=d, chunk.size=1000, verbose = T)
#  d.matrix <- df.to.matrix(x=d, chunk.size=0, verbose = T)
#  print(paste0("class: ", class(d), " size: ", nrow(d)," x ", ncol(d)))
#  d.matrix[1:10,1:10]
#  d.matrix[1:10,(ncol(d.matrix)-10):ncol(d.matrix)]

df.to.matrix <- function(x, chunk.size=50000, verbose = F)
{
  cols.x <- ncol(x)
  rows.x <- nrow(x)
  
  if(chunk.size<=1){
    x.matrix <- as.matrix(x)
  }else{
    m.steps <- my.seq(chunk.size, cols.x, chunk.size, T)
    begin <- 1
    chunk <- 1
    x.matrix.list <- list()
    for(end in m.steps){
      if(verbose)
        cat(paste0("Chunk ",chunk," of ",length(m.steps),"\n"))
      x.matrix.list[[chunk]] <- as.matrix(x[,begin:end])
      begin <- end + 1
      chunk <- chunk +1
    }
    if(verbose)
      cat(paste0("Binding chunks...\n"))
    x.matrix <- do.call("cbind", x.matrix.list)
  }
  return (x.matrix)
}

###############################
#string.empty.as.na
###############################
string.empty.as.na <- function(x) ifelse(x=="", NA, x)

###############################
#string.replace
###############################
string.replace <- function(x, sourceChars = c(" "), destinationChar = "_"){  
  for(i in 1:length(sourceChars))
    x <- gsub(paste0("\\", sourceChars[i]), destinationChar, x)  
  return(x)
}

###############################
#string.trim
###############################
string.trim <- function(str) gsub("^\\s+|\\s+$", "", str)

###############################
#string.starts.with
###############################
string.starts.with <- function(str, pattern, trim=FALSE, ignore.case=FALSE){
  pattern <- c(pattern)
  str <- c(str)
  if(trim)
    str <- string.trim(str)
  if(ignore.case){
    str <- tolower(str)
    pattern <- tolower(pattern)
  }
  ret <- rep(F, length(str))
  for(i in 1:length(pattern)){
    ret <- ret | substring(str, 1, nchar(pattern[i])) == pattern[i]
  }
  return(ret)    
}

###############################
#const.features
###############################
const.features <- function(x){
  same <- sapply(x, function(.col){
    all(is.na(.col))  || all(.col[1L] == .col)
  })
  which(same)
}

###############################
#normalize data - all columns (0,1)
###############################
normalize <- function(data) { 
  apply(data, 2, function(x) {xmin <- min(x, na.rm = T); (x-xmin)/(max(x, na.rm = T)-xmin)})
}

###############################
#scale.vector
###############################
#scaleX(c(-0.4,0,1,10),-10,5)
scale.vector<-function(x, min, max){
  minTmp <- min(x)
  maxTmp <- max(x)
  xTmp <- (x-minTmp)/(maxTmp-minTmp)
  xTmp <- (xTmp*(max-min))+min
  return(xTmp)
}

###############################
#get.files.names
###############################
get.files.names <- function(path, filter="*", ext=c('.csv','.rds'), fullNames=F, recursive=T){
  files <- NULL
  if(!File.exists(path)){
    stop(paste0("Directory: '",path,"' does not exist!"))
  }else{ 
    if(is.null(ext)){
      files <- list.files(path, pattern = NULL, full.names=fullNames, recursive=recursive, include.dirs=F)
    }else{
      for(i in 1:length(ext)){
        files <- c(files, list.files(path, pattern = paste0('\\',ext[i],'$'), full.names=fullNames, recursive=recursive, include.dirs=F ))
      }
    }
    #filter files
    filter <- gsub("([*])\\1+", "\\1", filter)  
    files <- files[files %in% files[grep(filter,files)]]    
  }
  return (files)
}

###############################
#File.exists
###############################
File.exists <- function(x) { 
  if (.Platform$OS == "windows" && grepl("[/\\]$", x)) {
    file.exists(dirname(x)) 
  } else file.exists(x) 
}

###############################
#File.ext
###############################
# path <- "//Users\\mdr.am.ins.ki\\Dropbox/DOCUM.ENTS//Money//ghfdjkhkj.hkfjdhk.EXD"
# file.extension(path)
file.extension <- function(x){
  file <- tail(unlist(strsplit(x, '[/\\]')),1)
  file <- tail(unlist(strsplit(file, '[.]')),1)
  file_ext <- tolower(file)
  return (file_ext)
}

###############################
#open.plot.file
###############################
open.plot.file <- function(filename, width, height, res = 72){
  ext <- file.extension(filename)
  if (ext == "png") {
    png(filename, width=width, height=height, res = 72)
  } else if (ext == "pdf") {
    # pdf size is set by default
    pdf(filename)
  } else if (ext == "svg") {
    # pdf size is set by default
    svg(filename)
  } else{
    png(paste0(filename,".png"), width=width, height=height, res = 72)
  }
  return(T)
}

###############################
#build.cmatrix
###############################
#x1 <- round(runif(100, 0.0, 3.0))
#x2 <- round(runif(100, 0.0, 3.0))
#build.cmatrix(x1,x2)
build.cmatrix <- function(real, predicted){
  cmatrix <- table(as.character(real), as.character(predicted))
  cmatrix <- as.data.frame.matrix(cmatrix)
  cmatrix <- as.matrix(cmatrix)
  colnMat <-  colnames(cmatrix)
  rownMat <- rownames(cmatrix)
  cmatrix <- cmatrix[,match(rownMat,colnMat)]
  colnames(cmatrix) <- rownMat
  rownames(cmatrix) <- rownMat
  cmatrix[is.na(cmatrix)] <- 0
  
  class(cmatrix) <- append("cmatrix", class(cmatrix))
  return(cmatrix)
}

###############################
#print.cmatrix
###############################
print.cmatrix <- function(cmatrix){
  if(!any(class(cmatrix) %in% "cmatrix"))
    stop("Input object is not 'cmatrix' class.")
  
  dg <- diag(cmatrix)
  TPR <- 100 * dg / rowSums(cmatrix)
  TPR <- round(TPR, digits=1)
  acc <- round(100 * calc.acc(cmatrix), digits = 1)
  wacc <- round(100 * calc.wacc(cmatrix), digits = 1)
  conf_matrix <- cmatrix / sum(cmatrix) * 100
  conf_matrix <- round(cmatrix, digits=1)
  writeLines("Confusion Matrix:\n")
  class(conf_matrix) <- "matrix"
  print(conf_matrix)
  writeLines("")
  writeLines("TPR:\n")
  print(as.data.frame(TPR))
  writeLines("")  
  writeLines(paste("Accuracy:",acc,"%"))
  writeLines(paste("wAccuracy:",wacc,"%"))
}

###############################
#calc.acc
###############################
calc.acc <- function(cmatrix){
  if(!any(class(cmatrix) %in% "cmatrix"))
    stop("Input object is not 'cmatrix' class.")
  
  dg <- diag(cmatrix)
  acc <- sum(dg) / sum(cmatrix)
  return(acc)
}

###############################
#calc.wacc
###############################
#calc.wacc(cmatrix)
calc.wacc <- function(cmatrix){
  if(!any(class(cmatrix) %in% "cmatrix"))
    stop("Input object is not 'cmatrix' class.")
  
  dg <- diag(cmatrix)
  TPR <- dg / rowSums(cmatrix)
  wacc <- mean(TPR)
  return(wacc)
}