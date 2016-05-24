`print.Rasch` <-
function(x, ...){
 if(! substr(x$model,1,3)=="mix"){
 itemtable <- cbind(difficulty = x$item.par$delta.i,
                    se  = x$item.par$SE.delta.i, x$item.par$in.out)
 if(! x$model=="RM"){
   tautable <- cbind(step1 = x$item.par$tau[1,],
                     se1 = x$item.par$SE.tau[1,])
   for(stepper in 2:x$i.stat$steps){
    newnames <- colnames(tautable)
    tautable <- cbind(tautable,  x$item.par$tau[stepper,],
                                 x$item.par$SE.tau[stepper,])
    colnames(tautable) <- c(newnames,paste("step",stepper,sep=""),paste("se",stepper,sep=""))
   }
 }

 print(round(itemtable,3))
 cat("\n")
 if(x$i.stat$steps > 1) print(round(tautable,3))
 cat("\n")
 print(x$run.time)
 cat("\n")
 cat("Number of iterations: ", x$iter, "\n \n")
 } # end ! mix
  else { itemtable <- cbind(difficulty = x$LatentClass[[1]]$item.par$delta.i,
                             se = x$LatentClass[[1]]$item.par$SE.delta.i, x$LatentClass[[1]]$item.par$in.out)
          for(c in 2:length(x$LatentClass)) itemtable <- cbind(itemtable, 
                                             cbind(difficulty = x$LatentClass[[c]]$item.par$delta.i,
                                             se = x$LatentClass[[c]]$item.par$SE.delta.i, 
                                             x$LatentClass[[c]]$item.par$in.out))
        print(round(itemtable,3))
        cat("\n")
        print(x$run.time)
        cat("\n")
        cat("Number of iterations: ", x$iter, "\n \n")
        }
 invisible(x)
}

