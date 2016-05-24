vitalsens <- function(elements, vitalrates)
{
   ## check if elements is matrix expression?
     #  expression(matrix2(c(0,F,G,S)))
    # or  expression(0,F,G,S)
   #  grepl("matrix", elements)
  
   if(is.vector(vitalrates)){vitalrates<-as.list(vitalrates)}
   if(!is.list(vitalrates)){stop("Vital rates should be a vector or list")}
   if(class(elements)!="expression"){stop("Matrix elements should be an expression")}
   ## check length of expression
   n<-sqrt(length(elements))
   if(n %% 1 !=0)
   {
      stop(paste("Length of element expression is", length(elements),
        "- Expecting power of 2 like 4,9,16 to form a square matrix"))
   }
   ## get values for matrix elements - enclos=NULL is used to restrict eval to names in vitalrates
   vrs<-try(sapply(elements, eval, vitalrates, NULL), silent=TRUE)
   if(class(vrs)=="try-error")
   {
      # keep useful part of error message
      vrs<- sub("Error in eval\\(expr, envir, enclos\\) :", "", vrs[1])
      stop(paste("Cannot evaluate element expression using given vital rates:", vrs))
   }
   ## store results in data.frame
   res<-data.frame(estimate=unlist(vitalrates), sensitivity=0, elasticity=0)
   ## create projection matrix
   A<-matrix(vrs, nrow=n, byrow=TRUE)
   eig <- eigen.analysis(A)        
   ## finally, get derivatives  
   deriv.funcs <- sapply(elements, deriv, namevec=names(vitalrates), function.arg=TRUE)
   devs <- lapply(deriv.funcs, function (x) do.call(x, vitalrates))
   for (i in 1:length(vitalrates))
   {
      derivs <- matrix( as.numeric(lapply(devs, function (x) attr(x, "gradient")[i])), nrow=n, byrow=TRUE)
     res[i,2] <-  sum(derivs*eig$sensitivities)
     res[i,3] <- vitalrates[[i]]/eig$lambda1*sum(derivs*eig$sensitivities)   
   }
   res
}
