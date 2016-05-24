getpvalClassif <-
function(x,y, nsim=500, test="normal"){
    if(is.function( test)){
        pval <- test(x,y)
    }else{
        if(test!="normal"){
            mud <- mean(x)-mean(y)
            cx <- c(x,y)
            ci <- c(rep(1,length(x)),rep(2,length(y)))
            muds <- numeric(nsim)
            for (sim in 1:nsim){
                ind <- sample( 1:length(cx),length(cx),replace=TRUE)
                muds[sim] <- mean(cx[ind[ ci[ind]==1]])-mean(cx[ind[ci[ind]==2]])
            }
            pval <- 2*(1-pt( abs(mud)/ sd(muds),df=nsim))
        }else{
            pval <- t.test(x,y)$p.value
        }
    }
    return(  pval  )
}
