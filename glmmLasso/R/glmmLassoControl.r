glmmLassoControl<-function(nue=1,index=NULL,smooth=NULL,start=NULL,q_start=NULL,center=TRUE,standardize=TRUE,  
                           steps=1000,method="EM", overdispersion=FALSE,epsilon=1e-4,maxIter=200,
                           print.iter=FALSE,print.iter.final=FALSE,method.final="EM", 
                           eps.final=1e-4, Q.fac=5, complexity="hat.matrix",...)
{                       

list(nue = nue, index=index, smooth = smooth,start = start, q_start = q_start, center=center, standardize = standardize,
     steps = steps, method = method, overdispersion = overdispersion, epsilon = epsilon, 
     maxIter = maxIter,print.iter = print.iter, print.iter.final = print.iter.final, 
     method.final = method.final, eps.final = eps.final, Q.fac = Q.fac, complexity = complexity, ...)
}
