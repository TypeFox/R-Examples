wrapper <- function(x,root.file.output,object,path.input,path.output,Y=Y,start.counter) {
  number <- x+start.counter-1
  root.file.output.perm <- paste(root.file.output,number,sep="")
  if(object$n > 50) indY <- sample(1:object$n,object$n) else{
    cond <- TRUE
    while (cond){
      indY <- sample(1:object$n,object$n)
      if(sum(c(1:object$n)==indY)==0) cond <- FALSE}
  }
  print(x)
  Yperm <- Y[indY,]
  
  name.Y <- paste("Yperm-",number,"-C-CODE.txt",sep="")
  nameY <- file.path(path.expand(path.input),name.Y)
  
  cat(object$n,"\n",object$q,"\n",file=nameY,sep="")
  write(t(Yperm), ncolumns=object$q,append = TRUE,file=nameY,sep="\t")
  
  
res.model.perm <- R2GUESS.perm(dataY=name.Y,dataX=object$dataX,path.inputx=object$path.input,path.inputy=path.input,path.output=path.output,path.par=object$path.par,path.init=object$path.init,file.par=object$file.par,file.init=object$file.init,file.log=NULL,nsweep=object$nsweep,burn.in=object$burn.in,root.file.output=root.file.output.perm,time=
     FALSE,top=object$top,history=FALSE,Egam=object$Egam,Sgam=object$Sgam,label.Y=object$label.Y,nb.chain=object$nb.chain,conf=0,cuda=object$cuda,MAP.file=object$MAP.file,p=object$p,q=object$q,n=object$n,time.limit=NULL,seed=NULL)
  
  
  eval(parse(text=paste("modelY.Perm",eval(parse(text="x")),"<- res.model.perm",sep="")))
  
  
 
  NameMarg.perm <- file.path(res.model.perm$path.output, paste(res.model.perm$root.file.output,"output_marg_prob_incl.txt",sep="_"))
  Marg.perm <- read.table(NameMarg.perm,header=TRUE)

  result <- list(MPI.perm=Marg.perm[,2],JF.best=max(res.model.perm$BestModels$jeffries))
  return(result)
}
