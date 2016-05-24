as.ESS.object <-function(dataY,dataX,path.input,path.output,root.file.output,
label.X=NULL,label.Y=NULL,path.par,path.init=NULL,file.par,file.init=NULL,file.log=NULL,MAP.file=NULL,command=TRUE){

if(!is.character(dataY)) stop("dataY should be a character vector")
if(!is.character(dataX)) stop("dataY should be a character vector")

Namefeatures <- file.path(path.output, paste(root.file.output,"features.txt",sep="_"))
features <- read.table(file=Namefeatures,header=TRUE)
rownames(features) <- features[,1]

n <- features["n","value"]
p <- features["p","value"]
q <- features["q","value"]
nsweep <- features["nsweeps","value"]

time <- features["time","value"]
if(time==1) time <- TRUE else time <- FALSE
conf <- features["n_conf","value"]
top <- features["top","value"]
nb.chain <- features["nb.chain","value"]
history <- features["history","value"]
if(history==1) history <- TRUE else history <- FALSE
Egam <- features["Egam","value"]
Sgam <- features["Sgam","value"]
burn.in <- features["burn.in","value"]
cuda <- features["cuda","value"]
if(cuda==1) cuda <- TRUE else cuda <- FALSE
seed <- features["seed","value"]

if(features["run_finished","value"]==1){ BestModels <- get.best.models(path.output,path.input,root.file.output,label.X=label.X,p=p,MAP.file) 
Finish <- TRUE
if(is.na(BestModels$modelName[1])) cat(paste("the best model (posterior=",signif(BestModels$postProb[1],digits=2),") is the model without covariable \n",sep=" "))
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



path.input <- path.expand(path.input)

path.output <- path.expand(path.output)
path.par <- path.expand(path.par)
if(is.null(file.log)) file.log <- root.file.output

if(command==TRUE) command <- readLines(file.path(path.output, paste(root.file.output,"command-C.txt",sep="-"))) else command <- NULL


res <- list(dataY=dataY,dataX=dataX,path.input=path.input,path.output=path.output,path.par=path.par, path.init=path.init, history=history,time=time,file.par =file.par,file.init=file.init,file.log=file.log,
root.file.output=root.file.output,nsweep=nsweep,top=top,BestModels=BestModels,label.X=label.X,
label.Y=label.Y,p=p,q=q,n=n,nb.chain=nb.chain,burn.in=burn.in,conf=conf,cuda=cuda,Egam=Egam,Sgam=Sgam,MAP.file=MAP.file,command=command,seed=seed,Finish=Finish)
class(res) <- "ESS"
return(res)
}


