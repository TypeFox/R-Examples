summary.fem <- function(object,...){
 cat('* Model: ')
 if (object$call[[1]]=='sfem'){
   if (object$crit=='bic') cat('the chosen model is',object$model,'with K =',object$K,'and lambda =',object$l1,'( bic =',object$bic,')\n')
   else if (object$crit=='aic') cat('the chosen model is',object$model,'with K =',object$K,'and lambda =',object$l1,'( aic =',object$bic,')\n')
   else if (object$crit=='icl') cat('the chosen model is',object$model,'with K =',object$K,'and lambda =',object$l1,'( icl =',object$icl,')\n')
 }
 else{
   if (object$crit=='bic') cat('the chosen model is',object$model,'with K =',object$K,'( bic =',object$bic,')\n')
   else if (object$crit=='aic') cat('the chosen model is',object$model,'with K =',object$K,'( aic =',object$bic,')\n')
   else if (object$crit=='icl') cat('the chosen model is',object$model,'with K =',object$K,'( icl =',object$icl,')\n')
 }
 if (nrow(object$U)<10){ cat('* Loading matrix:\n'); print(object$U)}
 else if(object$call[[1]]=='sfem') cat('* Loading matrix:',sum(rowSums(abs(object$U))>1e-5),'variables are active over the',nrow(object$U),'original ones\n')
}