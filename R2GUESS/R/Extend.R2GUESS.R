Extend.R2GUESS <- function(x,niter,time.limit=NULL){
  if(x$Finish==FALSE) stop("The previous run has not been reached the number of sweep \n you can use the fonction Resume.R2GUESS")
  
  if(!is.null(time.limit)){
    split.command <- unlist(strsplit(x$command," "))  
    ind <- which(split.command=="-timeLimit")
    split.command[ind+1] <- time.limit
    command <- paste(split.command,collapse=" ") 
  }else{
    command <- paste(x$command," -extend ",niter,sep="")}
  print(command)
  if (.Platform$OS.type == "unix") {
    system(command)
  } else if (.Platform$OS.type == "windows") {
    shell(command)
  }
  
  
  Namefeatures <- file.path(x$path.output, paste(x$root.file.output,"features.txt",sep="_"))
  features <- read.table(file=Namefeatures,header=TRUE)
  rownames(features) <- features[,1]
  if(features["run_finished","value"]==1){ BestModels <- get.best.models(x$path.output,x$path.input,x$root.file.output,label.X=x$label.X,p=x$p,x$MAP.file) 
                                           Finish <- TRUE
                                           cat("The run is ok","\n")
                                           cat("You can now analyse the results","\n")
  }
  else{
    cat("The run time reaches the specified time limit","\n")
    cat("You can use the function Resume.R2GUESS to resume the run ","\n")
    cat("You can also use the function PostProcess.R2GUESS to post-processing the current run ","\n")
    cat("You have also the possibility to extend your run to reach the number of sweep you wanted ","\n")
    BestModels <- NULL
    Finish <- FALSE
  }
  res <- list(dataY = x$dataY, dataX = x$dataX, path.input = x$path.input, 
              path.output = x$path.output, path.par=x$path.par, path.init=x$path.init, history = x$history, time = x$time, file.par =x$file.par,file.init=x$file.init,file.log=x$file.log,
              root.file.output = x$root.file.output, nsweep = features["last_sweep","value"], 
              top = x$top, BestModels = BestModels, label.X = x$label.X, 
              label.Y = x$label.Y, p = x$p, q = x$q, n=x$n, nb.chain = x$nb.chain,
              burn.in = x$burn.in,conf=x$conf,cuda=x$cuda,Egam=x$Egam,Sgam=x$Sgam,MAP.file=x$MAP.file,command=command,nsweep.extend=niter,seed=x$seed,Finish=Finish)
  class(res) <- "ESS"
  return(res)  
}
