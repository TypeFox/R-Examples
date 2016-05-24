Postprocess.R2GUESS <- function(x){
  if(x$Finish==TRUE) stop("You don't need to use this function as you run has already finished")
  
  command <- paste(x$command," -postProcess ",sep="")
  print(command)
  if (.Platform$OS.type == "unix") {
    system(command)
  } else if (.Platform$OS.type == "windows") {
    shell(command)
  }
  
  Namefeatures <- file.path(x$path.output, paste(x$root.file.output,"features.txt",sep="_"))
  features <- read.table(file=Namefeatures,header=TRUE)
  rownames(features) <- features[,1]


history.file <- file.path(
    x$path.output, paste(x$root.file.output, "output_models_history.txt", sep="_"))
model.history <- read.table(history.file, header=FALSE, sep="\t",
                            skip=x$burn.in + 1, stringsAsFactors=FALSE)
mppi.mc <- tabulate(
    as.integer(c(lapply(model.history[,5], strsplit, split=" ", fixed=TRUE),
               recursive=TRUE)),
    nbins=x$p
) / features["last_sweep","value"]
mppi.file <- file.path(
    x$path.output, paste(x$root.file.output, "output_marg_prob_incl_mc.txt", sep="_"))
write.table(mppi.mc, mppi.file, quote=FALSE, sep="\t",
            col.names=c("Predictor\tMarg_Prob_Incl"))


  BestModels <- get.best.models(x$path.output,x$path.input,x$root.file.output,label.X=x$label.X,p=x$p,x$MAP.file)
  res <- list(dataY = x$dataY, dataX = x$dataX, path.input = x$path.input, 
              path.output = x$path.output, path.par=x$path.par, path.init=x$path.init, history = x$history, time = x$time, file.par =x$file.par,file.init=x$file.init,file.log=x$file.log,
              root.file.output = x$root.file.output, nsweep = features["last_sweep","value"], 
              top = x$top, BestModels = BestModels, label.X = x$label.X, 
              label.Y = x$label.Y, p = x$p, q = x$q, n=x$n, nb.chain = x$nb.chain,
              burn.in = x$burn.in,conf=x$conf,cuda=x$cuda,Egam=x$Egam,Sgam=x$Sgam,MAP.file=x$MAP.file,command=command,nsweep.extend=x$nsweep.extend,seed=x$seed,Finish=x$Finish)
  class(res) <- "ESS"
  return(res)  
}
