R2GUESS.perm <-
function(dataY,dataX,path.inputx,path.inputy,path.output,path.par,path.init=NULL,file.par,file.init=NULL,file.log=NULL,nsweep,burn.in=1000,Egam=2,Sgam=2,
root.file.output,time=TRUE,top=100,history=TRUE,label.X=NULL,label.Y=NULL,nb.chain,conf=0,cuda=TRUE,MAP.file=NULL,p,q,n,time.limit=NULL,seed=NULL){


##   for (i in c(1:length(.libPaths()))) {
##     if ( file.exists(paste(.libPaths()[i],"R2GUESS",sep="/")) )
##       which.lib <- i
##   }
  # Setup file paths
  pack.root <- system.file(package = "R2GUESS")
  ESS.directory <- file.path(pack.root, "bin", .Platform$r_arch)


path.inputx <- path.expand(path.inputx)
path.inputy <- path.expand(path.inputy)


### read data


if(time==TRUE) time1 <- " -time " else time1 <- NULL
top1 <- paste(" -top ",top,sep="")
burn.in1 <- paste(" -burn_in ",burn.in,sep="")
if(history==TRUE) history1 <- " -history " else history1 <- NULL

if(cuda==TRUE) cuda1 <- " -cuda " else cuda1 <- NULL

#ESS.directory <- path.expand(ESS.directory)

path.par <- path.expand(path.par)
if(is.null(file.log)) {file.log <- root.file.output
                       log1 <- NULL} else{log1 <- " -log "}

if(!is.null(seed)) {seed.opt <- paste(" -seed ",seed,sep="")} else{seed.opt <- NULL}
if(is.null(time.limit)) time.limit <- 2000


if(is.null(time.limit)) option.timelimit <- NULL else option.timelimit <- paste(" -timeLimit ",time.limit," ",sep="")

guess.executable <- ifelse(.Platform$OS.type == "unix", "GUESS", "GUESS.exe")
if(!is.null(path.init)){
path.init <- path.expand(path.init)
command <- file.path(ESS.directory, paste(guess.executable, " -X \"",file.path(path.inputx,dataX),"\" -Y \"",file.path(path.inputy,dataY), "\" -par \"",file.path(path.par,file.par),"\" -init \"",file.path(path.init,file.init),"\" -nsweep ",nsweep,
burn.in1, " -out_full \"",file.path(path.output,root.file.output),"\" ",top1," -nconf ",conf,cuda1,time1,history1,option.timelimit," -Egam ",Egam," -Sgam ",Sgam," -n_chain ",nb.chain,seed.opt,log1," > \"", file.path(path.output,file.log),"_log\"",sep=""))
}else{
command <- file.path(ESS.directory, paste(guess.executable, " -X \"",file.path(path.inputx,dataX),"\" -Y \"",file.path(path.inputy,dataY), "\" -par \"",file.path(path.par,file.par),"\" -nsweep ",nsweep,
burn.in1, " -out_full \"",file.path(path.output,root.file.output),"\" ",top1," -nconf ",conf,cuda1,time1,history1,option.timelimit," -Egam ",Egam," -Sgam ",Sgam," -n_chain ",nb.chain,seed.opt,log1," > \"", file.path(path.output,file.log),"_log\"",sep=""))}


print(command)
if (.Platform$OS.type == "unix") {
  system(command)
} else if (.Platform$OS.type == "windows") {
  shell(command)
}

BestModels <- get.best.models(path.output,path.inputx,root.file.output,label.X=label.X,p=p,MAP.file)


res <- list(dataY = dataY, dataX = dataX, path.input = path.inputx,path.inputy=path.inputy,
        path.output = path.output, path.par=path.par, path.init=path.init, history = history, time = time,file.par =file.par,file.init=file.init,file.log=file.log,
        root.file.output = root.file.output, nsweep = nsweep,
        top = top, BestModels = BestModels, label.X = label.X,
        label.Y = label.Y, p = p, q = q, n=n, nb.chain = nb.chain,
        burn.in = burn.in,conf=NULL,cuda=cuda,Egam=Egam,Sgam=Sgam,MAP.file=MAP.file)
class(res) <- "ESS"
return(res)
}
