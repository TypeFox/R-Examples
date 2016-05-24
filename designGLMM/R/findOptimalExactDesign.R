# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' @export
#' @importFrom stats optim
findOptimalExactDesign<-function(numunits,means,sigma,probs=c(1),criterion="D",link="Poisson",trace=FALSE,iter=10000,temp=10,tmax=10,silent=FALSE,tol=0.0001,maxtime=60){
  blksize<-numunits
  numblock<-1
  # Generate a starting design (at the moment, it may not be a good design)
  workingdes<-sample(1:(length(means)),size=numblock*blksize,replace=TRUE)
  progressvec<-c()
  prevopt<-0
  reduction<-1
  time.start<-Sys.time()

  while(reduction>tol){
    time.now<-Sys.time()
    if(time.now-time.start>maxtime) break
    if(criterion=="D"){
      optdesSANN <- optim(workingdes, objfnD_CRD, updateDesign_CRD, length(means),sigma,means,probs, method = "SANN",
                          control = list(maxit = iter, temp = temp, trace = trace, REPORT = 1, tmax=tmax))
      current<--optdesSANN$value

    } else {
      optdesSANN <- optim(workingdes, objfnA_CRD, updateDesign_CRD, length(means),sigma,means,probs, method = "SANN",
                          control = list(maxit = iter, temp = temp, trace = trace, REPORT = 1, tmax=tmax))
      current<-optdesSANN$value

    }

    workingdes<-optdesSANN$par
    reduction<- (current - prevopt)/current
    progressvec<-append(progressvec,current)
    prevopt<-current
  }

  # Store and sort design
  optdes<-matrix(optdesSANN$par,ncol=blksize,byrow = TRUE)
  optdes2<-t(apply(optdes,1,function(x) sort(x)))
  optdes2<-optdes2[ do.call( order, as.list(as.data.frame(optdes2)) ), ]
  # Print output if silent is FALSE
  if(silent==FALSE){
    cat("\n\nThe optimal design is:\n")
    print(optdes2)
    if(criterion=="D"){
    cat(paste("The determinant of the information matrix is: ",round(current,digits=5),"\n",sep=""))
    }
    if(criterion=="A"){
      cat(paste("The average variance of orthogonal contrasts is: ",round(current,digits=5),"\n",sep=""))
    }
    cat(paste("Progressvec is: ",progressvec,"\n",sep=""))
  }
  list("design"=optdes2,"value"=current,"iter"=progressvec)
}
