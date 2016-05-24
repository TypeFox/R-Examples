FDR.permutation <- function(x,path.input=NULL,Npermut,start.counter=1,path.output=NULL,threshold=0.05,nbcpu=NULL,number.cutoff=200) {

root.file.output <- paste("Permut",x$root.file.output,sep="-")

if(is.null(path.output)) path.output <- x$path.output else
{path.output <- path.expand(path.output)}

if(is.null(path.input)) path.input <- x$path.input else
{path.input<- path.expand(path.input)}

Y <- read.table(file.path(path.expand(x$path.input), x$dataY),skip=2)

NameMarg <- file.path(x$path.output, paste(x$root.file.output,"output_marg_prob_incl.txt",sep="_"))
Marg <- read.table(NameMarg,header=TRUE)
MPI <- Marg[,2]
MPI.sort <- sort(MPI,decreasing=TRUE)

if(is.null(nbcpu)){sfInit(parallel=FALSE, cpus=1)}else{
  sfInit(parallel=TRUE, cpus=nbcpu)
  sfLibrary(package=R2GUESS,pos=4)
}
sfExportAll()
result <- sfLapply(1:Npermut,wrapper,root.file.output=root.file.output,object=x,path.input=path.input,path.output=path.output,Y=Y,start.counter=start.counter)

sfStop()

MPI.perm <- unlist(lapply(result,function(x){x$MPI.perm}))
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
      
      flag_emp = 1   # ***
    }
    
  }
  
}

int <- approx(FDR_emp_St, cutoff_St, n = 10000)
int_idx <- max(which(int$x < threshold))
FDR_emp_int <- int$x[int_idx]
cutoff_int <- int$y[int_idx]

return(list(cutoff.MPI = cutoff, cutoff_int = cutoff_int, cutoff_St = cutoff_St, FDR_emp = FDR_emp, FDR_emp_int = FDR_emp_int, FDR_emp_St = FDR_emp_St))
}

