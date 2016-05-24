R2GUESS <-
function(dataY,dataX,path.input,path.output,path.par,path.init=NULL,file.par,file.init=NULL,file.log=NULL,nsweep,burn.in,Egam,Sgam,
root.file.output,time=TRUE,top=100,history=TRUE,label.X=NULL,label.Y=NULL,choice.Y=NULL,nb.chain,conf=NULL,cuda=TRUE,MAP.file=NULL,time.limit=NULL,seed=NULL)
{

##   for (i in c(1:length(.libPaths()))) {
##     if ( file.exists(paste(.libPaths()[i],"R2GUESS",sep="/")) )
##       which.lib <- i
##   }
  # Setup file paths
  pack.root <- system.file(package = "R2GUESS")
  ESS.directory <- file.path(pack.root, "bin", .Platform$r_arch)


if(burn.in>nsweep) stop("the number of burnin should be lower than the number of sweep (nsweep)")


### read data
if(is.data.frame(dataX)){
  X <- dataX
  p <- dim(X)[2]
  n <- dim(X)[1]
  newdataX <- "data-X-C-CODE.txt"
  name.X <- file.path(path.expand(path.input),newdataX)
  cat(n,"\n",p,"\n",file=name.X,sep="")
  write(t(X), ncolumns=p,append = TRUE,file=name.X,sep="\t")
  }else{
  info <- readLines(file.path(path.expand(path.input),dataX),n=2)
  n <- as.integer(info[1])
  p <- as.integer(info[2])
  newdataX <- dataX
  }

#  q <- length(choice.Y)
  if(is.data.frame(dataY)) {
  Y <- dataY
  q0 <- dim(Y)[2]
  label.Y <- colnames(Y)
  file.Y <- FALSE
  if(is.null(choice.Y)) choice.Y <- 1:q0
  q <- length(choice.Y)
  }else{
  info <- readLines(file.path(path.expand(path.input),dataY),n=2)
  ny <- as.integer(info[1])
  q0 <- as.integer(info[2])
  Y <- read.table(file.path(path.expand(path.input),dataY),skip=2)
  if(is.null(choice.Y)) choice.Y <- 1:q0
  q <- length(choice.Y)
  if(ny!=dim(Y)[1]) stop("number of subject has to be specified in the first row of the txt file")
  if(q0!=dim(Y)[2]) stop("dimension of the phenotype has to be specified in the second row of the txt file")
  if(n!=ny) stop("the two data set (predictors and phenotype) have not compatible regarding the number of subject")
  colnames(Y) <- label.Y
  if(q==q0){ file.Y <- TRUE
             newdataY <- dataY
  }}

if(is.null(label.Y))  label.Y <- paste("Y",1:q0,sep=".")
if(is.null(colnames(Y))) colnames(Y) <- label.Y


if(file.Y==FALSE){
if(q==q0) nameY <- "ALL" else{
if(is.numeric(choice.Y)) choice.Y <- label.Y[choice.Y]
nameY <- paste(choice.Y,collapse="-")}
newdataY <- paste("data-Y-",nameY,"-C-CODE.txt",sep="")
name.Y <- file.path(path.expand(path.input),newdataY)
### Choice of the component of Y to be analysed
if(dim(Y)[2]!=1){
Y <- Y[,choice.Y]
label.Y <- colnames(Y)
}
### write data in txt file for C++ code
cat(n,"\n",q,"\n",file=name.Y,sep="")
write(t(Y), ncolumns=q,append = TRUE,file=name.Y,sep="\t")
}


if(!is.null(conf)){
if(is.data.frame(conf)){
  Z <- conf
  k <- dim(Z)[2]
  n <- dim(Z)[1]
  newdataZ <- "data-Z-C-CODE.txt"
  name.Z <- file.path(path.expand(path.input),newdataZ)
  cat(n,"\n",k,"\n",file=name.Z,sep="")
  write(t(Z), ncolumns=k,append = TRUE,file=name.Z,sep="\t")
  }else{
  info.Z <- readLines(file.path(path.expand(path.input),conf),n=2)
  n <- as.integer(info.Z[1])
  k <- as.integer(info.Z[2])
  Z <- read.table(file.path(path.expand(path.input),conf),skip=2)
  }
Z <- as.matrix(Z)
beta <- (solve(t(Z)%*%Z)%*%t(Z))%*%as.matrix(Y)
residu <- as.matrix(Y)-Z%*%beta
newdataY <- paste("residu-",newdataY,sep="")
name.conf <- file.path(path.expand(path.input),newdataY)
cat(n,"\n",q,"\n",file=name.conf,sep="")
write(t(residu), ncolumns=q,append = TRUE,file=name.conf,sep="\t")
#conf <- k
}



if(time==TRUE) time1 <- " -time " else time1 <- NULL
top1 <- paste(" -top ",top,sep="")
burn.in1 <- paste(" -burn_in ",burn.in,sep="")
if(history==TRUE) history1 <- " -history " else history1 <- NULL

if(cuda==TRUE) cuda1 <- " -cuda " else cuda1 <- NULL

#ESS.directory <- path.expand(ESS.directory)
path.input <- path.expand(path.input)
path.output <- path.expand(path.output)
path.par <- path.expand(path.par)
if(is.null(file.log)) {file.log <- root.file.output
                       log1 <- NULL} else{log1 <- " -log "}

if(!is.null(seed)) {seed.opt <- paste(" -seed ",seed,sep="")} else{seed.opt <- NULL}

if(is.null(time.limit)) time.limit <- 2000
confounder <- 0

if(is.null(time.limit)) option.timelimit <- NULL else option.timelimit <- paste(" -timeLimit ",time.limit," ",sep="")

guess.executable <- ifelse(.Platform$OS.type == "unix", "GUESS", "GUESS.exe")
if(!is.null(path.init)){
path.init <- path.expand(path.init)
command <- file.path(ESS.directory, paste(guess.executable, " -X \"",file.path(path.input,newdataX),"\" -Y \"",file.path(path.input,newdataY), "\" -par \"",file.path(path.par,file.par),"\" -init \"",file.path(path.init,file.init),"\" -nsweep ",nsweep,
burn.in1, " -out_full \"",file.path(path.output,root.file.output),"\" ",top1," -nconf ",confounder,cuda1,time1,history1,option.timelimit," -Egam ",Egam," -Sgam ",Sgam," -n_chain ",nb.chain,seed.opt,log1," > \"", file.path(path.output,file.log),"_log\"",sep=""))
}else{
command <- file.path(ESS.directory, paste(guess.executable, " -X \"",file.path(path.input,newdataX),"\" -Y \"",file.path(path.input,newdataY), "\" -par \"",file.path(path.par,file.par),"\" -nsweep ",nsweep,
burn.in1, " -out_full \"",file.path(path.output,root.file.output),"\" ",top1," -nconf ",confounder,cuda1,time1,history1,option.timelimit," -Egam ",Egam," -Sgam ",Sgam," -n_chain ",nb.chain,seed.opt,log1," > \"", file.path(path.output,file.log),"_log\"",sep=""))}

write(command,file=file.path(path.output, paste(root.file.output,"command-C.txt",sep="-")))
print(command)
if (.Platform$OS.type == "unix") {
    system(command)
} else if (.Platform$OS.type == "windows") {
    shell(command)
}

Namefeatures <- file.path(path.output, paste(root.file.output,"features.txt",sep="_"))
features <- read.table(file=Namefeatures,header=TRUE)
rownames(features) <- features[,1]
if(features["run_finished","value"]==1){
BestModels <- get.best.models(path.output,path.input,root.file.output,label.X=label.X,p=p,MAP.file)


history.file <- file.path(
    path.output, paste(root.file.output, "output_models_history.txt", sep="_"))
model.history <- read.table(history.file, header=FALSE, sep="\t",
                            skip=burn.in + 1, stringsAsFactors=FALSE)
mppi.mc <- tabulate(
    as.integer(c(lapply(model.history[,5], strsplit, split=" ", fixed=TRUE),
               recursive=TRUE)),
    nbins=p
) / (features["last_sweep","value"] - burn.in + 1)
mppi.file <- file.path(
    path.output, paste(root.file.output, "output_marg_prob_incl_mc.txt", sep="_"))
write.table(mppi.mc, mppi.file, quote=FALSE, sep="\t",
            col.names=c("Predictor\tMarg_Prob_Incl"))


Finish <- TRUE
cat("The run is ok","\n")
cat("You can now analyse the results","\n")
}
else{
cat("(i) The run time reaches the specified time limit","\n")
cat("(ii) You can use the function Resume.R2GUESS to resume the run ","\n")
cat("(iii) You can also use the function PostProcess.R2GUESS to post-process the current run ","\n")
cat("(iv) You have also the possibility to extend your run to reach the number of sweep you wanted ","\n")
BestModels <- NULL
Finish <- FALSE
}

res <- list(dataY = newdataY, dataX = newdataX, path.input = path.input,
        path.output = path.output, path.par=path.par, path.init=path.init, history = history, time = time, file.par =file.par,file.init=file.init,file.log=file.log,
        root.file.output = root.file.output, nsweep = features["last_sweep","value"],
        top = top, BestModels = BestModels, label.X = label.X,
        label.Y = label.Y, p = p, q = q, n=n, nb.chain = nb.chain,
        burn.in = burn.in,conf=conf,cuda=cuda,Egam=Egam,Sgam=Sgam,MAP.file=MAP.file,command=command,seed=features["seed","value"],Finish=Finish)
class(res) <- "ESS"
return(res)
}
