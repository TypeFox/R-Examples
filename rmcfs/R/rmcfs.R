###############################
#mcfs
###############################
mcfs <- function(formula, data,
                         projections = 3000,
                         projectionSize = 0.05,
                         splits = 5,
                         splitRatio = 0.66,
                         balanceRatio = 1,
                         splitSetSize = 1000,
                         cutoffPermutations = 20,
                         cutoffMethod = c("permutations", "criticalAngle", "kmeans", "mean"),
                         buildID = TRUE,
                         finalRuleset = TRUE,
                         finalCV = TRUE,
                         finalCVSetSize = 1000,
                         finalCVRepetitions = 3,
                         u = 1, 
                         v = 1,
                         seed = NA,
                         threadsNumber = 2)
{ 
  if(!cutoffMethod[1] %in% c("permutations", "criticalAngle", "kmeans", "mean")){
    stop(paste0("Incorrect 'cutoffMethod' = ", cutoffMethod[1]))
  }
  
  start.time <- Sys.time()
  cat("Checking input data...\n")
  label <- "input"
  coldef <- formula.prepare(formula, data)
  target <- coldef$target
  cols <- coldef$cols

  data <- cbind(data[,cols, drop=FALSE], data[,target, drop=FALSE])

  if(length(unique(data[,target])) == 1){
    stop("Decision attribute contains only 1 value.")
  }
  
  tmp_dir <- paste(tempdir(), .Platform$file.sep, sep="")
  tmp_dir <- gsub("\\\\", .Platform$file.sep, tmp_dir)
  libdir <- gsub("\\\\", .Platform$file.sep, libdir)
  
  config_file <- paste0(tmp_dir, "mcfs.run")
  input_file_name <- paste0(label, ".adx")
  
  params <- default.params
  params$inputFileName=input_file_name
  params$inputFilesPATH=tmp_dir
  params$resFilesPATH=tmp_dir
  params$mcfs.projections = projections
  params$mcfs.projectionSize = projectionSize
  params$mcfs.splits = splits
  params$mcfs.splitRatio = splitRatio
  params$mcfs.balanceRatio = balanceRatio  
  params$mcfs.splitSetSize = splitSetSize
  params$mcfs.cutoffPermutations = cutoffPermutations
  params$mcfs.cutoffMethod = cutoffMethod[1]
  params$mcfs.buildID = buildID
  params$mcfs.finalRuleset = finalRuleset
  params$mcfs.finalCV = finalCV
  params$mcfs.finalCVSetSize = finalCVSetSize
  params$mcfs.finalCVRepetitions = finalCVRepetitions
  params$mcfs.u = u
  params$mcfs.v = v
  params$mcfs.threadsNumber = threadsNumber
  if(is.numeric(seed)){
    params$mcfs.seed = seed
  }else{
    params$mcfs.seed = ""
  }

  cat("Exporting params...\n")
  save.params.file(params, config_file)
  
  cat("Exporting input data...\n")
  write.adx(data, file=paste0(tmp_dir, input_file_name), target=target)
  
  cat("Running MCFS-ID...\n")
  cdir <- getwd()
  setwd(libdir)
  #jri.prepare()
  .jcall("dmLab/mcfs/MCFS", returnSig="V", "main",
         .jarray(config_file))
  .jcheck(silent = FALSE)
  setwd(cdir)
  
  cat("Reading results...\n")
  mcfsResult <- import.result(tmp_dir, label)
  
  #clean temp files
  #cat("Cleaning temporary files...\n")
  out.files <- get.files.names(tmp_dir, filter=label, ext=c('.csv','.txt','.adx'), fullNames=T, recursive=F)
  for(i in 1:length(out.files)){
    if(file.exists(out.files[i])){
      #cat(paste0("remove ", out.files[i]))
      file.remove(out.files[i])
    }
  }
  cat("Done.\n")
  end.time <- Sys.time()
  mcfsResult$exec_time <- end.time - start.time

  class(mcfsResult) <- "mcfs"
  return(mcfsResult)
}

###############################
#jri.prepare
###############################
jri.prepare <- function(){
  .jaddClassPath(system.file("jri",c("JRI.jar"), package="rJava"))
  
  .jcall("java/lang/System", returnSig="V", "setOut",
         .jnew("java/io/PrintStream",
               .jcast(.jnew("org/rosuda/JRI/RConsoleOutputStream",
                            .jengine(TRUE), as.integer(0)),
                      "java/io/OutputStream")))
  
  .jcall("java/lang/System", returnSig="V", "setErr",
         .jnew("java/io/PrintStream",
               .jcast(.jnew("org/rosuda/JRI/RConsoleOutputStream",
                            .jengine(TRUE), as.integer(1)),
                      "java/io/OutputStream")))
}

###############################
#formula.prepare
###############################
formula.prepare <- function(formula, data) {
  target <- as.character(as.list(formula)[[2]])
  if(!target %in% names(data)){
    stop(paste0("Target feature: '", target, "' does not exist in the data."))
  }
  
  #function model.frame crashes under LINUX
  #ndata <- model.frame(formula, data, na.action=NULL)
  #mcfs.default(ndata[,(2:ncol(ndata))], ndata[,1], names(ndata)[1], ...)
  
  #R code below works perfectly under LINUX and other OS and is much faster
  cols <- all.vars(formula)
  cols <- cols[2:length(cols)]
  if(length(cols)==1 && cols=="."){
    cols <- names(data)
    cols <- cols[cols != target]
  }
  var.exists <- cols %in% names(data)
  cols.missing <- cols[!var.exists]
  if(length(cols.missing)>0){
    cols.missing.text <- paste(cols.missing, collapse=", ")
    stop(paste0("Features: '", cols.missing.text, "' does not exist in the data."))
  }
  retList <- list(target=target, cols=cols)
  return(retList)
}

###############################
#default.params
###############################
default.params <- list(verbose = "false", debug = "false", mcfs.progressShow = "false",
                       inputFileName = "", inputFilesPATH = "", resFilesPATH = "",
                       mcfs.projections = 3000, mcfs.projectionSize = 0.05,
                       mcfs.splits = 5, mcfs.splitRatio = 0.66,
                       mcfs.balanceRatio = 1,
                       mcfs.splitSetSize = 1000,
                       mcfs.contrastAttr = "false",
                       mcfs.contrastAttrThreshold = 1,
                       mcfs.cutoffPermutations = 20,
                       mcfs.cutoffAlpha = 0.05,
                       mcfs.cutoffAngle = 0.01,
                       mcfs.cutoffMethod = "mean",
                       mcfs.buildID = "true",
                       mcfs.finalRuleset = "true",
                       mcfs.finalCV = "true",
                       mcfs.finalCVSetSize = 1000,
                       mcfs.finalCVRepetitions = 3,
                       mcfs.u = 1, mcfs.v = 1,
                       mcfs.threadsNumber = 4,
                       mcfs.seed = NA,
                       mcfs.progressInterval = 10, mcfs.classifier = "j48", 
                       j48.useGainRatio = "true", j48.maxConnectionDepth = 5,
                       adx.useComplexQuality = "true", adx.qMethod = 2, 
                       sliq.useDiversityMeasure = "true"
)

###############################
#save.params.file
###############################
save.params.file <- function(params, config_file) {
  for(i in 1:length(params)) {
    if(typeof(params[[i]]) == "logical") {
      params[[i]] = tolower(as.character(params[[i]]))
    }
  }
  
  if(length(params[["inputFiles"]])>1){
    params[["inputFiles"]] <- paste0('[',paste0('',params[["inputFiles"]],'', collapse=', '),']')
  }
  
  if(is.null(params[["testFileName"]])){
    params[["testFileName"]] <- ""
  }
  
  f <- file(config_file, "w")
  cat(paste(names(params), params, sep="="), file=f, sep="\n")
  close(f)
}

###############################
#write.adx
###############################
# data(alizadeh)
# d <- alizadeh
# d[5,1] <- Inf
# d[7,2] <- Inf
# d[5,3] <- NA
# d[7,4] <- NA
# write.adx(d, file="~/ali_test.arff", chunk.size=2000)
write.adx <- function(x, file = "", target = NA, chunk.size = 100000)
{
  if(file == ""){
    file <- stdout()
  }else if(is.character(file)) {
    file <- file(file, "wb")
    on.exit(close(file))
  }
  
  if(!inherits(file, "connection"))
    stop("Argument 'file' must be a character string or connection.")
  
  if (!is.data.frame(x))
    x <- data.frame(x)
  
  if(is.character(target)) 
    target <- match(target, names(x))
  if(is.na(target) || target < 1 || target > ncol(x)) 
    target <- ncol(x)
  
  # for speed create local variables
  # call these functions only once
  names.x <- names(x)
  cols.x <- ncol(x)
  rows.x <- nrow(x)

  verbose <- F
  if(cols.x * rows.x > chunk.size)
    verbose <- T

  if(verbose)
    cat(paste0("Determining numeric columns...\n"))
  numeric.cols <- sapply(x[1,], is.numeric)
  
  if(verbose)
    cat(paste0("Saving attributes meta info...\n"))
  ##### write attributes #####
  writeLines("attributes", file)
  writeLines("{", file)
  
  attrHeader <- matrix(nrow = length(names.x), ncol = 3, byrow = FALSE)
  attrHeader[,1] <- paste0(" '",names.x,"'")
  attrHeader[numeric.cols, 2] <- " numeric"
  attrHeader[!numeric.cols, 2] <- " nominal"
  attrHeader[target,3] <- " decision(all)"
  attrHeader[is.na(attrHeader[,3]),3] <- ""
  l <- apply(attrHeader, 1, paste, collapse="")
  cat(l, file=file, sep="\n")
  
  writeLines("}", file)
  writeLines("", file)

  # convert input data to matrix
  if(verbose)
    cat("Conversion of input data to character matrix...\n") 
  x <- df.to.matrix(x, chunk.size, verbose)

  ###### write events ######
  writeLines("events", file)
  writeLines("{", file)

  col.steps <- my.seq(chunk.size, cols.x, chunk.size, T)
  rstep <-  ceiling(chunk.size / cols.x)
  row.steps <- my.seq(rstep, rows.x, rstep, T)
  chunks <- length(col.steps) * length(row.steps)
  
  #i <- row.steps[1]
  #j <- col.steps[1]
  chunk <- 1
  row.begin <- 1
  for(i in row.steps) {
    col.begin <- 1
    for(j in col.steps) {
      if(verbose)
        cat(paste0("Saving events chunk: ", chunk, " of ", chunks,"\n"))
      # m must be matrix type
      m <- x[row.begin:i, col.begin:j, drop=FALSE]
      m <- fix.matrix(m)
      l <- apply(m, 1, paste, collapse=",")
      if(j > col.steps[1])
        cat(',', file=file, sep="")
      if(length(l)>1 | j==cols.x){
        cat(l, file=file, sep="\n")
      }else{
        cat(l, file=file, sep="")
      }
      rm(m,l)
      col.begin <- j + 1
      chunk <- chunk + 1
    }
    row.begin <- i + 1
  }
  writeLines("}", file)
  
  if(verbose){
    cat("Data has been exported.\n")
  }
}

###############################
#write.arff
###############################
# data(alizadeh)
# x <- alizadeh
# x$nominalAttr <- c(rep("val1",nrow(x)/2),rep("val2",nrow(x)))[1:nrow(x)]
# d[5,1] <- Inf
# d[7,2] <- Inf
# d[5,3] <- NA
# d[7,4] <- NA
# write.arff(x, file="~/ali_test.arff", chunk.size=2000)

write.arff <- function(x, file = "", target = NA, chunk.size = 100000)
{
  if(file == "")
    file <- stdout()
  else if(is.character(file)) {
    file <- file(file, "wb")
    on.exit(close(file))
  }
  
  if(!inherits(file, "connection"))
    stop("Argument 'file' must be a character string or connection.")
  
  #if (!is.data.frame(x) && !is.matrix(x))
  if (!is.data.frame(x))
    x <- data.frame(x)
  
  if(is.character(target)) 
    target <- match(target, names(x))
  if(is.na(target) || target < 1 || target > ncol(x)) 
    target <- ncol(x)
  
  if(target!=ncol(x)){
    #move decision column to the end
    decisionName <- names(x)[target]
    x <- cbind(x[,-target], x[,target])
    target <- ncol(x)
    names(x)[target] <- decisionName
  }  
  
  # for speed create local variables
  # call these functions only once
  names.x <- names(x)
  cols.x <- ncol(x)
  rows.x <- nrow(x)

  verbose <- F
  if(cols.x * rows.x > chunk.size)
    verbose <- T
  
  if(verbose)
    cat(paste0("Determining numeric columns...\n"))
  numeric.cols <- sapply(x[1,], is.numeric)
  
  if(verbose)
    cat(paste0("Saving attributes meta info...\n"))
  ## Write Attributes
  writeLines(paste0("@relation ",'"',names(x)[target],'"'), file)  
  writeLines("", file)
  
  attrHeader <- matrix(nrow = length(names.x), ncol = 3, byrow = FALSE)
  attrHeader[,1] <- "@attribute"
  attrHeader[,2] <- paste0(" '",names.x,"'")
  attrHeader[numeric.cols, 3] <- " real"
  attrHeader[!numeric.cols, 3] <-   sapply(x[,!numeric.cols, drop=F], function(x) paste0(" {", paste(unique(x), collapse=","),"}")) 
  attrHeader[is.na(attrHeader[,3]),3] <- ""
  l <- apply(attrHeader, 1, paste, collapse="")
  cat(l, file=file, sep="\n")
  
  # convert input data to matrix
  if(verbose)
    cat("Conversion of input data to character matrix...\n") 
  x <- df.to.matrix(x, chunk.size, verbose)
  
  ## Write Events
  writeLines("", file)
  writeLines("@data", file)
  writeLines("", file)
  
  col.steps <- my.seq(chunk.size, cols.x, chunk.size, T)
  rstep <-  ceiling(chunk.size / cols.x)
  row.steps <- my.seq(rstep, rows.x, rstep, T)
  chunks <- length(col.steps) * length(row.steps)
  
  chunk <- 1
  row.begin <- 1
  for(i in row.steps) {
    col.begin <- 1
    for(j in col.steps) {
      if(verbose)
        cat(paste0("Saving events chunk: ",chunk," of ", chunks,"\n"))
      # m must be matrix type
      m <- x[row.begin:i, col.begin:j, drop=FALSE]
      m <- fix.matrix(m)
      l <- apply(m, 1, paste, collapse=",")
      if(j > col.steps[1])
        cat(',', file=file, sep="")
      if(length(l)>1 | j==cols.x){
        cat(l, file=file, sep="\n")
      }else{
        cat(l, file=file, sep="")
      }
      rm(m,l)
      col.begin <- j + 1
      chunk <- chunk + 1
    }
    writeLines("", file)
    row.begin <- i + 1
  }

  if(verbose){
    cat("Data has been exported.\n")
  }
}

###############################
#fix.matrix
###############################
fix.matrix <- function(m){
  m[is.na(m)] = "?"
  m[is.infinite(m)] <- "?"
  if(class(m[1,1])=="character"){
    m <- string.trim(m)
    m[m=="Inf"] <- "?"
    m[m=="-Inf"] <- "?"
  }  
  return (m)
}

###############################
#fix.data
###############################
# x$mydate <- as.Date("2007-06-22")
# x$myposix <- as.POSIXct(x$mydate)
# x$diff <- x$mydate - as.Date("2010-06-22")
# showme(x)
# reshape2::melt(lapply(x, class))
# x <- fix.data(x)
# reshape2::melt(lapply(x, class))
# showme(x) 
fix.data <- function(x, type = c("all", "names", "values", "types"), 
                     source.chars=c(" ",",","/","|","#"), 
                     destination.char = "_", 
                     numeric.class=c("difftime"), 
                     nominal.class=c("factor", "logical", "Date", "POSIXct", "POSIXt"))
{
  if (!is.data.frame(x))
    x <- as.data.frame(x)
  
  if(type[1] %in% c("all","names") ){
    cat("Fixing names...\n")
    names(x) <- string.replace(string.trim(names(x)), source.chars, destination.char)
  }
  
  if(type[1] %in% c("all","values") ){
    cat("Fixing values...\n")
    x <- fix.data.values(x, source.chars, destination.char)
  }
  
  if(type[1] %in% c("all","types") ){
    cat("Fixing types...\n")
    x <- fix.data.types(x, numeric.class, nominal.class)
  }  
  return(x) 
}

###############################
#fix.data.values
###############################
# data(alizadeh)
# d <- alizadeh
# d <- alizadeh[4000:ncol(alizadeh)]
# d$art <- rep("aa|bb,cc##|dd",nrow(d))
# d$art2 <- rep("  aac bb  ",nrow(d))
# d$art2[3:13] <- ""
# fix.data.values(d)
fix.data.values <- function(x, sourceChars = c(" ",",","/","|","#"), destinationChar = "_")
{
  nominal.class <- c("character", "factor", "Date", "POSIXct", "POSIXt")
  df.class <- reshape2::melt(lapply(x, class))
  colnames(df.class) <- c("classname", "colname")
  nominalMask <- df.class$classname %in% nominal.class
  
  x.nominal <- x[,nominalMask, drop=F]
  #dplyr version
  #x.nominal <- x.nominal %>% mutate_each(funs(as.character)) %>% 
  #mutate_each(funs(string.trim)) %>% 
  #mutate_each(funs(string.empty.as.na))
  #apply version
  #x.nominal <- apply(x.nominal, c(1,2), FUN = as.character)
  #x.nominal <- apply(x.nominal, c(1,2), FUN = string.trim)
  #x.nominal <- apply(x.nominal, c(1,2), FUN = string.empty.as.na)

  x.nominal <- apply(x.nominal, c(1,2), FUN = string.replace, sourceChars=sourceChars, destinationChar=destinationChar)
  x[nominalMask] <- x.nominal
  
  return(x)
}

###############################
#fix.data.types
###############################
fix.data.types <- function(x, numeric.class = c("difftime"), 
                           nominal.class = c("factor", "logical", "Date", "POSIXct", "POSIXt"))
{
  df.class <- reshape2::melt(lapply(x, class))
  colnames(df.class) <- c("classname", "colname")
  fixCols.toNumeric <- unique(df.class$colname[df.class$classname %in% numeric.class])
  fixCols.toNominal <- unique(df.class$colname[df.class$classname %in% nominal.class])
  mask.toNumeric <- df.class$colname %in% fixCols.toNumeric
  mask.toNominal <- df.class$colname %in% fixCols.toNominal
  
  if(any(mask.toNumeric)){
      x.toCast <- x[,mask.toNumeric, drop=F]
      x.toCast <- apply(x.toCast, c(1,2), FUN = as.numeric)
      x[mask.toNumeric] <- x.toCast
  }
  if(any(mask.toNominal)){
      x.toCast <- x[,mask.toNominal, drop=F]
      x.toCast <- apply(x.toCast, c(1,2), FUN = as.character)
      x[mask.toNominal] <- x.toCast
  }
  return(x)
}

###############################
#filter.data
###############################
filter.data <- function(data, mcfs_result, size = NA){
  
  if(class(mcfs_result)!="mcfs")
    stop("Input object is not 'mcfs' class.")
  
  size <- get.size.param(size, mcfs_result$cutoff_value)
  if(is.null(size))
    return(NULL)
  
  fdata <- data[,names(data) %in% as.character(head(mcfs_result$RI,size)$attribute)]
  target.data.frame <- data.frame(data[,mcfs_result$target])
  if(ncol(target.data.frame)==1){
    names(target.data.frame) <- mcfs_result$target
    fdata <- cbind(fdata, target.data.frame)
  }
  return(fdata)
}

###############################
#read.ID (interdependencies)
###############################
read.ID <- function(fileName){
  interdeps <- NULL
  if(file.exists(fileName)){
      file.header <- readLines(fileName, n = 1, warn = FALSE)
      header.vector <- string.replace(unlist(strsplit(file.header, ",")), c("\""),"")
      if(all(c("edge_a","edge_b") %in% header.vector)){
        interdeps <- read.csv.result(fileName)
      }else{
        interdeps <- read.ID.list(fileName)
      }
  }
  return(interdeps)
}

###############################
#read.ID.list
###############################
read.ID.list <- function(fileName) {
  f <- file(fileName, "r")
  lines <- readLines(f)
  close(f)
  
  a <- strsplit(lines, "[,(]")
  
  process_row <- function(row) {
    a <- row[1]
    b <- row[2:length(row)]
    
    dim(b) <- c(2, length(b)/2)
    b <- t(b)
    return(cbind(a, b))
  }
  
  m <- lapply(a, process_row)
  d <- do.call("rbind", lapply(a, process_row))
  if(is.null(d))
    return (NULL)
  
  weights <- sapply(d[,3],function(s) return(as.numeric(substr(s, 1, nchar(s)-1))))
  d <- data.frame(edge_a=d[,1], edge_b=d[,2], weight=weights, stringsAsFactors = FALSE)
  d <- d[order(-d[,3]),]
  
  position <- 1:nrow(d)
  d <- cbind(position,d)
  
  return(d)
}

###############################
#read.cmatrix
###############################
read.cmatrix <- function(fileName){
  cmatrix <- read.table(fileName, sep=",", header = TRUE, stringsAsFactors = FALSE)
  cmatrix <- as.matrix(cmatrix[,2:ncol(cmatrix)])
  if(any(colnames(cmatrix) %in% c("other"))){
    otherIdx <- which(colnames(cmatrix) == c("other"))
    cmatrix <- cmatrix[-otherIdx,-otherIdx]
  }
  rownames(cmatrix) <- as.character(colnames(cmatrix))
  
  class(cmatrix) <- append("cmatrix", class(cmatrix))
  return (cmatrix)
}

###############################
#read.target
###############################
read.target <- function(fileName){
  matrix <- read.csv.result(fileName)
  target <- colnames(matrix)[1]
  return (target)
}

###############################
#read.RI
###############################
read.RI <- function(fileName){
  ranking <- read.table(fileName, sep=",", header = TRUE, stringsAsFactors = FALSE)
  ranking <- ranking[, names(ranking) %in% c('attribute', 'projections', 'classifiers', 'crudeRI', 'nodes', 'RI_norm')]
  # crudeRI is now nodes
  names(ranking)[names(ranking) %in% c('crudeRI')] <- 'nodes'
  ranking <- ranking[order(-ranking$RI_norm),]
  position <- 1:nrow(ranking)
  ranking <- cbind(position,ranking)
  return (ranking)
}

###############################
#read.csv.result
###############################
read.csv.result <- function(fileName){
  resultDataFrame <- NULL
  if(file.exists(fileName))
    resultDataFrame <- read.table(fileName, sep=",", header = TRUE, stringsAsFactors = FALSE)
  
  return (resultDataFrame)
}

###############################
#read.jrip
###############################
read.jrip <- function(fileName){
  resultText <- NULL
  if(file.exists(fileName))
    resultText <- readChar(fileName, file.info(fileName)$size)
  
  return(resultText)
}

###############################
#read.params
###############################
read.params <- function(fileName){
  params <- NULL
  if(file.exists(fileName)){
    paramsTXT <- readChar(fileName, file.info(fileName)$size)
    paramsTXT <- string.replace(paramsTXT, c("="), " : ")
    params <- yaml.load(paramsTXT)
  }
  return(params)
}

###############################
#import.result
###############################
import.result <- function(path, label){
  
  ri_file <- file.path(path, paste0(label, "__importances.csv"))
  if(!File.exists(ri_file))
    ri_file <- file.path(path, paste0(label, "__RI.csv"))
  
  id_file <- file.path(path, paste0(label, "_connections.csv"))
  if(!File.exists(id_file))
    id_file <- file.path(path, paste0(label, "_ID.csv"))
  
  distances_file <- file.path(path, paste0(label, "_distances.csv"))
  matrix_file <- file.path(path, paste0(label, "_cmatrix.csv"))
  cutoff_file <- file.path(path, paste0(label, "_cutoff.csv"))
  cv_file <- file.path(path, paste0(label, "_cv_accuracy.csv"))
  permutations_file <- file.path(path, paste0(label, "_permutations.csv"))
  topRanking_file <- file.path(path, paste0(label, "_topRanking.csv"))
  jrip_file <- file.path(path, paste0(label, "_jrip.txt"))
  params_file <- file.path(path, paste0(label, ".run"))
  
  ri <- read.RI(ri_file)
  id <- read.ID(id_file)
  cmatrix <- read.cmatrix(matrix_file)
  target <- read.target(matrix_file)
  dist <- read.csv.result(distances_file)
  cutoff <- read.csv.result(cutoff_file)
  cv_accuracy <- read.csv.result(cv_file)
  permutations <- read.csv.result(permutations_file)
  topRanking <- read.csv.result(topRanking_file)
  jrip <- read.jrip(jrip_file)
  params <- read.params(params_file)
  
  mcfsResult <- list()
  mcfsResult$target <- target
  mcfsResult$RI <- ri
  mcfsResult$ID <- id
  mcfsResult$distances <- dist
  mcfsResult$cmatrix <- cmatrix
  mcfsResult$cutoff <- cutoff
  if(nrow(topRanking)==0){
    mcfsResult$cutoff_value <- 0
  }else{
    mcfsResult$cutoff_value <- topRanking[nrow(topRanking),]$position
  }
  mcfsResult$cv_accuracy <- cv_accuracy
  mcfsResult$permutations <- permutations
  mcfsResult$jrip <- jrip
  mcfsResult$params <- params
  
  class(mcfsResult) <- "mcfs"
  return(mcfsResult)
}

###############################
#export.result
###############################
export.result <- function(mcfs_result, path, label = "rmcfs", save.rds = FALSE){
  
  if(class(mcfs_result) != "mcfs")
    stop("Input object is not 'mcfs' class.")
  
  dir.create(file.path(path), showWarnings=F, recursive=T)
  
  if(save.rds){
    saveRDS(mcfs_result,file=paste0(path, label, "_.rds"))
  }else{
    ri_file <- file.path(path, paste0(label, "__RI.csv"))
    id_file <- file.path(path, paste0(label, "_ID.csv"))
    distances_file <- file.path(path, paste0(label, "_distances.csv"))
    matrix_file <- file.path(path, paste0(label, "_cmatrix.csv"))
    cutoff_file <- file.path(path, paste0(label, "_cutoff.csv"))
    cv_file <- file.path(path, paste0(label, "_cv_accuracy.csv"))
    permutations_file <- file.path(path, paste0(label, "_permutations.csv"))    
    topRanking_file <- file.path(path, paste0(label, "_topRanking.csv"))
    jrip_file <- file.path(path, paste0(label, "_jrip.txt"))
    params_file <- file.path(path, paste0(label, ".run"))
    
    write.csv(mcfs_result$RI, file=ri_file, row.names = F)
    write.csv(mcfs_result$ID, file=id_file, row.names = F)
    write.csv(mcfs_result$distances, file=distances_file, row.names = F)
    
    cmatrix <- as.data.frame(mcfs_result$cmatrix)
    cmatrix <- cbind(rownames(cmatrix),cmatrix)
    colnames(cmatrix)[1] <- mcfs_result$target
    write.csv(cmatrix, file=matrix_file, row.names = F)
    
    write.csv(mcfs_result$cutoff, file=cutoff_file, row.names = F)
    #save cv_accuracy
    if(any(names(mcfs_result)=="cv_accuracy")){
      write.csv(mcfs_result$cv_accuracy, file=cv_file, row.names = F)
    }    
    #save permutations
    if(any(names(mcfs_result)=="permutations")){
      write.csv(mcfs_result$permutations, file=permutations_file, row.names = F)
    }
    
    #save top ranking
    topRanking <- head(mcfs_result$RI, mcfs_result$cutoff_value)
    position <- 1:nrow(topRanking)
    topRanking <- cbind(position, topRanking)
    write.csv(topRanking, file=topRanking_file, row.names = F)
    
    if(any(names(mcfs_result)=="jrip")){
      writeChar(mcfs_result$jrip, jrip_file)
    }
    
    if(any(names(mcfs_result)=="params")){
      save.params.file(mcfs_result$params, params_file)
    }
  }
}

###############################
#artificial.data
###############################
artificial.data <- function(rnd.features = 500, size = c(40, 20, 10), corruption = c(0,2,4), seed = NA){
  if(length(size) != 3){
    warning("Length of 'size' parameter does not equal to 3. Default values c(40, 20, 10) are used.")
    corruption <- c(40, 20, 10)
  }
  if(length(corruption) != 3){
    warning("Length of 'corruption' parameter does not equal to 3. Default values c(0, 2, 4) are used.")
    corruption <- c(0, 2, 4)
  }
  
  if(is.numeric(seed))
    set.seed(seed)
  class <- c(rep("A",40), rep("B",20), rep("C",10))
  A <- B <- C <- rep("0",length(class))
  A[class=="A"] <- "A"
  B[class=="B"] <- "B"
  C[class=="C"] <- "C"
  rnd <- runif(length(class))
  A[class=="A"][rnd[class=="A"]<=(sort(rnd[class=="A"]))[corruption[1]]] <- "0"
  B[class=="B"][rnd[class=="B"]<=(sort(rnd[class=="B"]))[corruption[2]]] <- "0"
  C[class=="C"][rnd[class=="C"]<=(sort(rnd[class=="C"]))[corruption[3]]] <- "0"
  d <- data.frame(matrix(runif(rnd.features*length(class)), ncol=rnd.features))
  d <- cbind(d,data.frame(A1=A, A2=A, B1=B, B2=B, C1=C, C2=C, class))
  return(d)
}

###############################
#build.idgraph
###############################
build.idgraph <- function(mcfs_result, size = NA, size_ID = NA, self_ID = FALSE, 
                          plot_all_nodes = FALSE, size_ID_mult = 3, size_ID_max = 100) {

  if(class(mcfs_result)!="mcfs")
    stop("Input object is not 'mcfs' class.")
  
  if(all(names(mcfs_result)!="ID")){
    stop("ID-Graph edges are not collected. Object 'mcfs_result$ID' does not exist.")
  }
  
  size <- get.size.param(size, mcfs_result$cutoff_value)
  if(is.null(size))
    return(NULL)
  
  min_ID <- get.min.ID(mcfs_result, size, size_ID, size_ID_mult, size_ID_max)

  plot_minW <- 1
  plot_maxW <- 7  
  vertexMinSize <- 3
  vertexMaxSize <- 12
  
  #add weightNorm and color columns to ranking
  ranking <- mcfs_result$RI  
  ranking$attribute <- as.character(ranking$attribute)  
  ranking$color <- scale.vector(ranking$RI_norm,0,1)
  ranking$color <- abs(ranking$color-1)
  
  #add weightNorm and color columns to interdeps
  interdeps <- mcfs_result$ID
  if(is.null(interdeps)){
    warning("ID-Graph is empty. Change input parameters and try to build it again.")
    return (NULL)
  }
  
  if(!self_ID)
    interdeps <- interdeps[interdeps$edge_a != interdeps$edge_b,]
  
  interdeps$weightNorm <- scale.vector(interdeps$weight,plot_minW,plot_maxW)
  interdeps$color <- scale.vector(interdeps$weight,0,1)
  interdeps$color <- abs(interdeps$color-1)
  
  #select interdeps to plot  
  nodes_to_keep <- head(ranking, size)
  interdeps <- interdeps[is.element(interdeps$edge_a, nodes_to_keep$attribute),]
  interdeps <- interdeps[is.element(interdeps$edge_b, nodes_to_keep$attribute),]
  interdeps <- interdeps[as.numeric(interdeps$weight) >= min_ID,]
  
  #select ranking to plot
  gNodes <- nodes_to_keep
  if(plot_all_nodes==F){
    gNodes <- unique(c(interdeps$edge_a,interdeps$edge_b))
    gNodes <- nodes_to_keep[nodes_to_keep$attribute %in% gNodes,]
  }
  cat(paste0("Selected ",nrow(gNodes)," nodes and ", nrow(interdeps)," edges.\n"))
  g <- igraph::graph.empty()
  if(nrow(interdeps) > 0) {    
    #add nodes
    for(i in 1:nrow(gNodes)){
      c <- rgb(1, gNodes$color[i], gNodes$color[i])
      g <- g + igraph::vertices(gNodes$attribute[i], shape="circle", color=c, size=10)
    }
    #add edges
    for(i in 1:nrow(interdeps)){
      c <- rgb(interdeps$color[i], interdeps$color[i], interdeps$color[i])
      g <- g + igraph::edge(interdeps$edge_a[i],interdeps$edge_b[i], weight=interdeps$weightNorm[i],
                            width=interdeps$weightNorm[i], color=c)
    }
    
    from <- NULL
    to <- NULL
    vertex.attributes <- NULL
    #set size of the nodes based on number of connections
    #V(g)[3]
    #E(g)[from(3)]
    #E(g)[to(3)]
    m <- 1:length(igraph::V(g))
    vertexSizeFrom <- sapply(m, function(x) length(igraph::E(g)[from(x)]))
    vertexSizeTo <- sapply(m, function(x) length(igraph::E(g)[to(x)]))
    vertexSize <- vertexSizeFrom + vertexSizeTo

    V(g)$size <- scale.vector(vertexSize,vertexMinSize,vertexMaxSize)
    #vertex.attributes(g)$size <- scale.vector(vertexSize,vertexMinSize,vertexMaxSize)
  }
  
  class(g) <- append("idgraph", class(g))
  return(g)  
}

###############################
#get.min.ID
###############################
get.min.ID <- function(mcfs_result, size = NA, size_ID = NA, size_ID_mult = 3, size_ID_max = 100){

  size <- get.size.param(size, mcfs_result$cutoff_value)
  if(is.null(size))
    return(NULL)
  
  top_attributes <- head(mcfs_result$RI,size)$attribute
  mask_attributes <- mcfs_result$ID$edge_a %in% top_attributes & mcfs_result$ID$edge_b %in% top_attributes
  
  if(is.na(size_ID)){
    top_ID <- head(mcfs_result$ID[mask_attributes,], min(size*size_ID_mult, size_ID_max))
  }else{
    top_ID <- head(mcfs_result$ID[mask_attributes,], size_ID)
  }
  min_ID <- min(top_ID$weight)
  
  return(min_ID)
}

###############################
#print.mcfs
###############################
print.mcfs <- function(x, ...){
  mcfs_result <- x
  if(class(mcfs_result)!="mcfs")
    stop("Input object is not 'mcfs' class.")

  writeLines(paste0("##### MCFS-ID result (s = ", mcfs_result$params$mcfs.projections,", t = ",mcfs_result$params$mcfs.splits, ", m = ",mcfs_result$params$mcfs.projectionSize, ") #####"))
  writeLines(paste0("Target feature: '",mcfs_result$target,"'"))
  writeLines("")
  writeLines(paste0("Top ",mcfs_result$cutoff_value," features:"))
  print(head(mcfs_result$RI[,c("position", "attribute", "RI_norm")], mcfs_result$cutoff_value), row.names = F)
  writeLines("")
  writeLines("#################################")
  writeLines("Cutoff values:")
  print(mcfs_result$cutoff, row.names = F)
  writeLines("")
  writeLines("#################################")
  writeLines("Confusion matrix obtained on randomly selected (st) datasets:")
  print.cmatrix(mcfs_result$cmatrix)
  writeLines("")
  if(any(names(mcfs_result) %in% c("jrip"))){
    writeLines("#################################")
    writeLines(paste0("JRIP classification rules created on top ", mcfs_result$cutoff_value," features:"))
    writeLines(mcfs_result$jrip)
  }  
  writeLines("#################################")
  writeLines(paste0("MCFS-ID execution time: ", format(mcfs_result$exec_time, digits=1),""))
}
