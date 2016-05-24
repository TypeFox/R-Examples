functional.beta.pair<-function (x, traits, index.family = "sorensen") 
{
    index.family <- match.arg(index.family, c("jaccard", "sorensen"))
    fbc<-x
    if (!inherits(x, "functional.betapart")) {
        fbc <- functional.betapart.core(x,traits, multi=FALSE, warning.time=FALSE, return.details=FALSE)
    } # end of computing core results
    
    switch(index.family, sorensen = {
        funct.beta.sim <- fbc$min.not.shared/(fbc$min.not.shared + fbc$shared)
        
        funct.beta.sne <- ((fbc$max.not.shared - fbc$min.not.shared)/((2 * fbc$shared) + fbc$sum.not.shared)) * (fbc$shared/(fbc$min.not.shared + fbc$shared))
        
        funct.beta.sor <- fbc$sum.not.shared/(2 * fbc$shared + fbc$sum.not.shared)
        
        functional.pairwise <- list(funct.beta.sim = as.dist(funct.beta.sim), funct.beta.sne = as.dist(funct.beta.sne), funct.beta.sor = as.dist(funct.beta.sor))
    								}, 
    
    					 jaccard = {
        funct.beta.jtu <- (2 * fbc$min.not.shared)/((2 * fbc$min.not.shared) + fbc$shared)
        
        funct.beta.jne <- ((fbc$max.not.shared - fbc$min.not.shared)/(fbc$shared + fbc$sum.not.shared)) * (fbc$shared/((2 * fbc$min.not.shared) + fbc$shared))
        
        funct.beta.jac <- fbc$sum.not.shared/(fbc$shared + fbc$sum.not.shared)
        
        functional.pairwise <- list(funct.beta.jtu = as.dist(funct.beta.jtu), funct.beta.jne = as.dist(funct.beta.jne), funct.beta.jac = as.dist(funct.beta.jac))
    								}
    								
    ) # end of switch
    
    return(functional.pairwise)
    
} # end of function 