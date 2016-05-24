`pw.assoc` <- 
function(formula, data, weights=NULL, freq0c=NULL)
{

## code for Cramer's V
    V <- function(tab){
        chi <- chisq.test(tab)
        mm <- min(nrow(tab)-1, ncol(tab)-1)
        out<-sqrt(chi$statistic/(sum(tab)*mm))
        names(out) <- NULL
        out
    }

### function for computing measures of
### proportional reduction of the variance Row|Column
    prv.rc <- function(tab){
        tab <- tab/sum(tab)
        rS <- rowSums(tab)
        cS <- colSums(tab)
## code for Goodman & Kruskal lambda(Row|Column)    
        V.r <- 1 - max(rS)
        EV.rgc <- 1-sum(apply(tab,2,max))
        lambda <- (V.r - EV.rgc)/V.r
## code for Goodman & Kruskal tau(Row|Column) 
        V.r <- 1 - sum(rS^2)
        a <- colSums(tab^2)
        EV.rgc <- 1 - sum(a/cS)
        tau <- (V.r - EV.rgc)/V.r
## code for Theil's Uncertainty(Row|Column)    
        V.r <- (-1)*sum(rS*log(rS))
        cS.mat <- matrix(cS, nrow=nrow(tab), ncol=ncol(tab), byrow=TRUE)
        EV.rgc <- sum(tab *log(tab/cS.mat))
        u <- (V.r + EV.rgc)/V.r
## output
        c(lambda.rc=lambda, tau.rc=tau, U.rc=u)
    }
###################################################
###################################################
    if(is.null(weights)) ww <- rep(1, nrow(data))
    else{
        ww <- data[,weights]
        data <- data[,setdiff(colnames(data), weights)]
    }
    df <- model.frame(formula=formula, data=data)
    lab <- colnames(df)
    p.x <- length(lab) - 1
    vV <- vlambda <- vtau <- vU <- numeric(p.x)
    for(i in 1:p.x){
        pos <- i+1
        form <- paste(lab[1], lab[pos], sep="+")
        form <- paste("ww", form, sep="~")
        tab <- xtabs(as.formula(form), data=df)
        if(is.null(freq0c)) freq0c <- 1/sum(tab)^2
        tab[tab==0] <- freq0c
        vV[i] <- V(tab)
        appo <- prv.rc(tab) 
        vlambda[i] <- appo[1]
        vtau[i] <- appo[2]
        vU[i] <- appo[3]
    }
    lab.am <- paste(lab[1], lab[-1], sep="." )
    names(vV) <- names(vlambda) <- names(vtau) <- names(vU) <- lab.am
    list(V=vV, lambda=vlambda, tau=vtau, U=vU)
}
