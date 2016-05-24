#'print.state.trans
#'
#'Generates output for state.trans object, see: \code{\link{StateTrans}}
#'
#'
#'@param x a state.trans object, that should be printed
#'@param ... further arguments passed to or from other methods.
#'@export


print.state.trans<-function(x, ...){

  if(length(x)==1){
    print(unclass(x))
    
  } else {
    if(class(x)[2]!="state.trans") warning("x should be a state.trans object!")

    l<-length(x)

    sum_tab<-x[[1]]

      for(i in 2:length(x)){
        sum_tab<-sum_tab+x[[i]]
      }

    mean_tab<-sum_tab/length(x)

    cat("\n\n A list of" , length(x), "state-transition tables! \n If you want to inspect a single case use [[case]] \n If you want to inspect several cases use [[from:to]] \n \n mean frequencies are:\n\n")
    print(mean_tab)

    return(mean_tab)
  }
}



