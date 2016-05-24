summary.ESS <- function(object,...){
x <- object
if (!inherits(x, "ESS")) stop("use only with \"ESS\" objects")

cat(" \n")
cat("Multiple Sparse Bayesian Regression model (MSBR)", "\n")
cat("fitted by GPU-based Evolutionary Stochastic Search algorithm (ESS)", "\n")

cat(" \n")
cat("Statistical Model:", "\n")
cat("  Dataset:","\n")
cat(paste("     Dataset Y:", x$dataY),"\n")
if(!is.null(x$label.Y)){
cat(paste("     Name of the phenotype:"),"\n")
cat(paste("  ",paste(x$label.Y,collapse=" ")),"\n")
cat(" \n")
}
cat(paste("     Dataset X:", x$dataX),"\n")
cat(paste("     Number of observation:", x$n),"\n")
if(x$q>1){
cat(paste("     Number of phenotypes  (q):", x$q),"\n")}else{
cat(paste("     Number of phenotype  (q):", x$q),"\n")}
cat(paste("     Number of predictors (p):", x$p),"\n")

cat(" \n")
cat(" \n")
cat("Parameter for ESS algorithm:", "\n")
if (x$cuda==TRUE) {
cat("     Use the Graphics Processing Unit (GPU) "," \n")
}
cat(paste("     The a priori average model size:", x$Egam),"\n")
cat(paste("     The 'a priori' standard deviation of the model size:", x$Sgam),"\n")
cat(paste("     Number of sweeps for MCMC:", x$nsweep),"\n")
cat(paste("     Number of sweeps to discarded:", x$burn.in),"\n")
cat(paste("     Number of chains:",x$nb.chain),"\n")
cat(" \n")
if (!is.null(x$path.init)) {
cat("MCMC start for the first run with the variable specified in:","\n")
cat(paste(x$path.init,x$file.init,sep=""),"\n")}

cat("Path for the folder containing the output:","\n")
cat(paste("  ",x$path.output),"\n")
cat("The Name of the files containing results of MCMC start by:","\n")
cat(paste("  ",x$root.file.output),"\n")
cat("Path for the folder containing the dataset:","\n")
cat(paste("  ",x$path.input),"\n")
cat("Path for the xml-formatted file containing the main parameter for MCMC: ","\n")
cat(paste("   ",x$path.par,x$file.par,sep=""),"\n")
cat(" \n")
cat(" \n")

res <- summarybest.model(x,...)

cat(paste("Summary of the ",res$N, "best models:"),"\n")

print(res$result)


}






