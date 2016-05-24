`get.Splines` <-
function(expr.data){

y <- expr.data[,2]
tp <- expr.data[,1]

gamres <- vector("list",length(tp)-2)

        gamres[[1]] = gam(y~1) # constant fit
        gamres[[2]] = gam(y~tp) # straight line fit

     for ( i in 3:(length(tp)-2)){

        gamres[[i]] = gam(y~s(tp,(i-1)))

        }
        
        # anova tests for best gam
        aovtab <- NULL
        for(model in 1:length(gamres)) {
            aovtab <-
            rbind(aovtab,(as.matrix(anova(gamres[[1]],gamres[[model]],test="Chisq"))[2,]))
        }
        
        # select best curve
        ix <- which(aovtab[-1,5]==min(aovtab[-1,5]))[1]

        return(gamres[[ix+1]])
        
}

