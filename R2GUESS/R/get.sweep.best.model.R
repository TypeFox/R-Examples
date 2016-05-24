get.sweep.best.model <- function(x,models=1){
  
  if(class(x)!="ESS") stop("x should be an ESS object")  
  path.output <- x$path.output
  file.in <- x$root.file.output
  model.history <- file.path(path.output, paste(file.in,"output_models_history.txt",sep="_"))
  
  Tmpmodel.history <- readLines(model.history)
  
  Tmpmodel.history <- Tmpmodel.history[-(1:2)]
  Tmpmodel.history <- Tmpmodel.history[-(1:x$burn.in)]
  
  res.model <- sapply(matrix(x$BestModels$modelName,ncol=1,nrow=x$top),FUN=function(x) unlist(strsplit(x," ")))
  names(res.model) <- 1:x$top
  res.model <- res.model[models]
  
  if(!is.null(x$MAP.file)){
    if(is.data.frame(x$MAP.file)){label.X <- as.character(x$MAP.file$SNPName)}else{
      NameMap.file <- file.path(x$path.input,x$MAP.file)
      SNPLabels <- read.table(NameMap.file,header=TRUE,stringsAsFactors = FALSE)
    label.X <- SNPLabels$SNPName
    }
  res.model <- lapply(res.model,FUN=function(x,nameX) as.character(which(label.X%in%x==TRUE)),nameX=label.X)
  }
  
  result <- NULL
  for (i in 1:(x$nsweep-x$burn.in)){
    tmpRow<-unlist(strsplit(Tmpmodel.history[i],"\t"))
    nameModel <- unlist(strsplit(tmpRow[5]," "))
    res <- sapply(res.model,FUN=function(x,TEMP){setequal(x,TEMP)},TEMP=nameModel)
    result <- rbind(result,res)
  }
  res <- list(result=result)
}




