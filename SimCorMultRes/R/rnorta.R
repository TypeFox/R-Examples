rnorta <-
function(R=R,cor.matrix=cor.matrix,distr="normal")
{
 if(!is.numeric(R) | R<1) 
    stop("'R' must be greater than or equal to one")
 R <- as.integer(R)
 distrs <- c("normal","logistic","extreme","cauchit")
 if(!is.element(distr,distrs)) 
    stop("'distr' must be either 'normal','logistic','extreme' or 'cauchit'") 
 ans <- rsmvnorm(R=R,cor.matrix=cor.matrix)
 if(distr=="logistic")  ans <- qlogis(pnorm(ans)) 
 if(distr=="extreme")   ans <- qgumbel(pnorm(ans))  
 if(distr=="cauchit")   ans <- qcauchy(pnorm(ans))
 ans
}

