


print.iBMA.intermediate.surv<- function(x, ...)
{
   cs<- x
   cat("\n iBMA intermediate result \n")
   cat("\nCall:\n", deparse(cs$call), "\n\n", sep = "")

    examined<- 1:(cs$nextVar - 1)
    examined<- examined[!(examined %in% cs$new.vars)]
    all.names<- colnames(cs$sortedX)[examined]
    current<- cs$currentSet
    nextSet<- c(cs$currentSet, cs$new.vars)
           
   cat("completed ", cs$currIter,  " iteration", c("", "s")[1 + (cs$currIter != 1)], "\n",sep="")
   cat("number of iterations in last call = ", cs$nIter,  "\n",sep="")
   cat("last variable examined = ", ifelse(sum(examined) > 0, max(examined), ""),  "\n",sep="")
   cat("variables currently selected =  ", ifelse(sum(current) > 0, paste(current, collapse = " ", sep = ""), ""),  "\n",sep="")
   cat("next set to examine = ", paste(nextSet, collapse = " ", sep = ""), "\n\n\n", sep = "")
 
   invisible(x)  
}



summary.iBMA.intermediate.surv<- function(object, ...)
{
   cs<- object
   cat("\n iBMA intermediate result \n")
   cat("\nCall:\n", deparse(cs$call), "\n\n", sep = "")
   
   ### if show_models show them here before statistics, (don't want user to think it is finished because 
   ### iteration statistics are scrolled off screen


    
    examined<- 1:(cs$nextVar - 1)
    examined<- examined[!(examined %in% cs$new.vars)]
    all.names<- colnames(cs$sortedX)[examined]
    current<- cs$currentSet
    nextSet<- c(cs$currentSet, cs$new.vars)
    dropped<- examined[!(examined %in% current)]
        
    status<- rep("currently selected", times = length(examined))
    status[dropped]<- "dropped"
    tble<- cbind( all.names,  status, cs$maxProbne0[examined], cs$nTimes[examined])
                   
    colnames(tble)<- c("variable", "status", "max Probne0", "n times in model")
    row.names(tble)<- examined    
        
### basic statistics
   cat("completed ", cs$currIter,  " iteration", c("", "s")[1 + (cs$currIter != 1)], "\n",sep="")
   cat("number of iterations in last call = ", cs$nIter,  "\n",sep="")
   cat("last variable examined = ", ifelse(sum(examined) > 0, max(examined), ""),  "\n",sep="")
   cat("variables currently selected =  ", ifelse(sum(current) > 0, paste(current, collapse = " ", sep = ""), ""),  "\n",sep="")
   cat("next set to examine = ", paste(nextSet, collapse = " ", sep = ""), "\n\n\n", sep = "")
   cat("statistics for variables examined so far \n ")
 
                       
 ### details
        print.default(tble, print.gap = 2, quote = FALSE, ...)         


                     
    invisible(object)  

}






print.iBMA.surv<- function(x, ...)
{
    cs<- x
    cat("\n iBMA final result \n")
    cat("\nCall:\n", deparse(cs$call), "\n\n", sep = "")
    
    examined<- 1:(cs$nextVar - 1)
    examined<- examined[!(examined %in% cs$last.added)]
    all.names<- colnames(cs$sortedX)[examined]
    current<- cs$currentSet
           
    cat("completed ", cs$currIter,  " iteration", c("", "s")[1 + (cs$currIter != 1)], "\n",sep="")
    cat("number of iterations in last call = ", cs$nIter,  "\n",sep="")
    cat("selected variables =  ", ifelse(sum(current) > 0, paste(current, collapse = " ", sep = ""), ""),  "\n",sep="")


    cat("\n\n Results of running BMA on final set of selected variables:\n\n")
 
    print(x$bma, ...)
    invisible(x)  
}





summary.iBMA.surv<- function(object, ...)
{
   cs<- object
   cat("\n iBMA.glm intermediate result \n")
   cat("\nCall:\n", deparse(cs$call), "\n\n", sep = "")
   
   ### if show_models show them here before statistics, (don't want user to think it is finished because 
   ### iteration statistics are scrolled off screen

    all.names<- colnames(cs$sortedX)
    status<- rep("dropped", times = cs$nVar)
    status[cs$selected]<- "selected"
    tble<- cbind( all.names,  status, cs$maxProbne0, cs$nTimes)
                   
    colnames(tble)<- c("variable", "status", "max Probne0", "n times in model")
    row.names(tble)<- 1:cs$nVar  
        
### basic statistics
   cat("completed ", cs$currIter,  " iteration", c("", "s")[1 + (cs$currIter != 1)], "\n",sep="")
   cat("number of iterations in last call = ", cs$nIter,  "\n",sep="")
   cat("variables selected =  ", ifelse(sum(cs$selected) > 0, paste(cs$selected, collapse = " ", sep = ""), ""),  "\n",sep="")
 
   cat("statistics for variables examined so far \n ")
 
 ### details
        print.default(tble, print.gap = 2, quote = FALSE, ...)         


    cat("\n\n Results of running BMA on final set of selected variables:\n\n")
 
    summary(object$bma, ...)
                     
    invisible(object)  

}


print.iBMA.intermediate.bicreg<- print.iBMA.intermediate.glm<- print.iBMA.intermediate.surv
print.iBMA.bicreg<- print.iBMA.glm<- print.iBMA.surv

summary.iBMA.intermediate.bicreg<- summary.iBMA.intermediate.glm<- summary.iBMA.intermediate.surv
summary.iBMA.bicreg<- summary.iBMA.glm<- summary.iBMA.surv
