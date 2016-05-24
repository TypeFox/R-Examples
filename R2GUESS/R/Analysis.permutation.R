Analysis.permutation <- function(x,Npermut,namePermut,threshold=0.05,path.output,number.cutoff=200){
NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl.txt",sep="_"))
Marg <- read.table(NameMarg,header=TRUE)
MPI <- Marg[,2]
MPI.sort <- sort(MPI,decreasing=TRUE)

JF.best <- NULL
MPI.perm <- NULL

vecPermut <- 1:Npermut

Npermut <- length(vecPermut)
for (i in vecPermut){
print(i)
if(class(x)!="ESS"){
res.model.perm <- eval(parse(text=paste(namePermut,i,sep="")))
}else{
root.file.output <- paste(namePermut,i,sep="")
res.model.perm  <- as.ESS.object(dataY=x$dataY,dataX=x$dataX,path.input=x$path.input,path.output=path.output,
root.file.output=root.file.output,label.X=x$label.X,
label.Y=x$label.Y,path.par=x$path.par,path.init=x$path.init,file.par=x$file.par,file.init=x$file.init,file.log=x$file.log,MAP.file=x$MAP.file,command=FALSE)

NameMarg.perm <- file.path(res.model.perm$path.output, paste(res.model.perm$root.file.output,"output_marg_prob_incl.txt",sep="_"))
Marg.perm <- read.table(NameMarg.perm,header=TRUE)
MPI.perm <- c(MPI.perm,Marg.perm[,2])
}
JF.best <- c(JF.best,max(res.model.perm$BestModels$jeffries)) 
}

MPI.perm.sort <- sort(MPI.perm,decreasing=TRUE)

flag_emp  <- 0
cutoff_St <- NULL
FDR_emp_St <- NULL
s_pp <- MPI.sort
n_pp <- length(s_pp)
s_pp_S <- MPI.perm.sort 

for (s in 1:number.cutoff) {
  n_s_pp <- sum(s_pp >= s_pp[s])
  n_s_pp_S <- sum(s_pp_S >= s_pp[s]) / Npermut
  # print(c(s, s_pp[s], n_s_pp_S / n_s_pp))
  
  if (s == 1) {
    cutoff_St <- s_pp[s]
    FDR_emp_St <- n_s_pp_S / n_s_pp
  }
  
  if (s > 1) {
    cutoff_St <- c(cutoff_St, s_pp[s])
    FDR_emp_St <- c(FDR_emp_St, n_s_pp_S / n_s_pp)
    
    if (((n_s_pp_S / n_s_pp) > threshold) & flag_emp == 0) {
      n_s_pp <- sum(s_pp >= s_pp[(s - 1)])
      n_s_pp_S <-  sum(s_pp_S >= s_pp[(s - 1)]) / Npermut
      
      cutoff <- s_pp[s - 1]
      FDR_emp <- n_s_pp_S / n_s_pp
      
      flag_emp <- 1   # ***
    }
    
  }
  
}

int <- approx(FDR_emp_St, cutoff_St, n = 10000)
int_idx <- max(which(int$x < threshold))
FDR_emp_int <- int$x[int_idx]
cutoff_int <- int$y[int_idx]

return(list(cutoff.MPI = cutoff, cutoff_int = cutoff_int, cutoff_St = cutoff_St, FDR_emp = FDR_emp, FDR_emp_int = FDR_emp_int, FDR_emp_St = FDR_emp_St))
}
#################################################
