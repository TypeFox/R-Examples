##' textSummary function
##'
##' A function to print a text description of the inferred paramerers beta and eta from a call to the function lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars or lgcpPredictMultitypeSpatialPlusPars
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars 
##' @param digits see the option "digits" in ?format 
##' @param scientific see the option "scientific" in ?format
##' @param inclIntercept logical: whether to summarise the intercept term, default is FALSE.
##' @param printmode the format of the text to return, can be 'LaTeX' (the default) or 'text' for plain text.
##' @param ... other arguments passed to the function "format" 
##' @return A text summary, that can be pasted into a LaTeX document and later edited.
##' @export

textSummary <- function(obj,digits=3,scientific=-3,inclIntercept=FALSE,printmode="LaTeX",...){
    psu <- summary(obj)
    parnames <- rownames(psu)

    form <- function(x,...){
        xtxt <- format(x,digits=digits,scientific=scientific,...)
        if(length(grep("e",xtxt))>0){
            spl <- unlist(strsplit(xtxt,"e"))
            if(printmode=="LaTeX"){
                xtxt <- paste(spl[1],"$\\times10^{",as.character(as.numeric(spl[2])),"}$",sep="")
            }
            else if (printmode=="text"){
                xtxt <- paste(spl[1],"x10^",as.character(as.numeric(spl[2])),sep="")   
            }
            else{
                stop("Undefine printmode. See ?textSummary.")
            }
        }    
        return(xtxt)
    }
    
    parvals <- matrix(sapply(psu,form),nrow=nrow(psu),ncol=ncol(psu)) 
    
    sigf <- apply(psu[,2:3],1,function(x){ifelse((x[1]<1 & x[2]<1)|(x[1]>1 & x[2]>1),TRUE,FALSE)})

    nspatial <- ncol(obj$etasamp)
    nhazard <- ncol(obj$omegasamp)
    np <- nrow(psu) - nspatial - nhazard
    
    
    lfsen <- "A summary of the parameters of the spatial latent field is as follows. "
    for(i in 1:nspatial){
        if(i==1){
            lfsen <- paste(lfsen,"The parameter ",parnames[np+nhazard+i]," had median ",parvals[np+nhazard+i,1]," (95\\% CRI ",parvals[np+nhazard+i,2]," to ",parvals[np+nhazard+i,3],")",sep="")
            if(nspatial==1){
                lfsen <- paste(lfsen,".",sep="")
            }
            else{
                lfsen <- paste(lfsen,";",sep="")   
            }
        }
        else{
            lfsen <- paste(lfsen," the parameter ",parnames[np+nhazard+i]," had median ",parvals[np+nhazard+i,1]," (95\\% CRI ",parvals[np+nhazard+i,2]," to ",parvals[np+nhazard+i,3],")",sep="")
            if(i==nspatial){
                lfsen <- paste(lfsen,".",sep="")
            }
            else{
                lfsen <- paste(lfsen,";",sep="")   
            }
        }
    }
    
    hasen <- "A summary of the parameters of the baseline hazard function is as follows. "
    for(i in 1:nhazard){
        if(i==1){
            hasen <- paste(hasen,"The parameter ",parnames[np+i]," had median ",parvals[np+i,1]," (95\\% CRI ",parvals[np+i,2]," to ",parvals[np+i,3],")",sep="")
            if(nhazard==1){
                hasen <- paste(hasen,".",sep="")
            }
            else{
                hasen <- paste(hasen,";",sep="")   
            }
        }
        else{
            hasen <- paste(hasen," the parameter ",parnames[np+i]," had median ",parvals[np+i,1]," (95\\% CRI ",parvals[np+i,2]," to ",parvals[np+i,3],")",sep="")
            if(i==nhazard){
                hasen <- paste(hasen,".",sep="")
            }
            else{
                hasen <- paste(hasen,";",sep="")   
            }
        }
    }        

    tab <- parvals[1:np,]
    tsi <- sigf[1:np]
    efd <- sapply(psu[,1],function(x){ifelse(x>1,1,-1)})[1:np] 
    prn <- parnames[1:np]
    
    if(!inclIntercept){
        if(any(prn=="\\verb=(Intercept)=")){
            idx <- which(prn=="\\verb=(Intercept)=")
            tab <- tab[-idx,,drop=FALSE]
            tsi <- tsi[-idx] # significant?
            efd <- efd[-idx] # effect direction
            prn <- prn[-idx]
        }
    }
    
    sumrow <- function(infor){
        i <- infor[1]
        dir <- infor[2]
        direc <- "increase"
        if(dir==-1){
            direc <- "reduction"
        }
        return(paste("each unit increase in ",pr[i]," led to a ",direc," in relative risk with median ",tt[i,1]," (95\\% CRI ",tt[i,2]," to ",tt[i,3],")",sep=""))
    }
    
    pardesc <- c()
    pardesc2 <- c()
    
    if(all(tsi)){
        pardesc <- "All of the main effects were found to be significant: "
        tt <- tab[tsi,drop=FALSE,]
        ef <- efd[tsi]
        pr <- prn[tsi]
        ord <- order(ef)
        tt <- tt[ord,,drop=FALSE]
        ef <- ef[ord]
        pr <- pr[ord]
        apply(cbind(1:nrow(tt),ef),1,function(x){pardesc <<- paste(pardesc,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
    }
    else{
        if(any(tsi)){
            pardesc <- "The following effects were found to be significant: "
            tt <- tab[tsi,,drop=FALSE]
            ef <- efd[tsi]
            pr <- prn[tsi]
            ord <- order(ef)
            tt <- tt[ord,,drop=FALSE]
            ef <- ef[ord]
            pr <- pr[ord]
            apply(cbind(1:nrow(tt),ef),1,function(x){pardesc <<- paste(pardesc,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
        }
        
        
        
        if(all(!tsi)){
            pardesc2 <- "None of the main effects were found to be significant: "
        }
        else{
            pardesc2 <- "The remainder of the main effects were not found to be significant: "
        }
        tt <- tab[!tsi,,drop=FALSE]
        ef <- efd[!tsi]
        pr <- prn[!tsi]
        ord <- order(ef)
        tt <- tt[ord,,drop=FALSE]
        ef <- ef[ord]
        pr <- pr[ord]
        apply(cbind(1:nrow(tt),ef),1,function(x){pardesc2 <<- paste(pardesc2,sumrow(x),ifelse(x[1]==nrow(tt),". ","; "),sep="")})
    }      
    #browser()

    ret <- list(spatial=lfsen,hazard=hasen,sigpar=pardesc,notsigpar=pardesc2)
    class(ret) <- c("textSummary","list")

    return(ret)

}

##' print.textSummary function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method print textSummary
##' @param x an object inheriting class textSummary
##' @param ... additional arguments, not used here 
##' @return prints a text summary of 'x' to the console
##' @export

print.textSummary <- function(x,...){
    cat("\n")
    cat(x$spatial)
    cat("\n")
    cat("\n")
    cat(x$hazard)
    cat("\n")
    cat("\n")
    cat(x$sigpar)
    cat("\n")
    cat("\n")
    cat(x$notsigpar)
    cat("\n")
} 




##' Summarise function
##'
##' A function to completely summarise the output of an object of class mcmcspatsurv.
##'
##' @param obj an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars 
##' @param digits see the option "digits" in ?format 
##' @param scientific see the option "scientific" in ?format
##' @param inclIntercept logical: whether to summarise the intercept term, default is FALSE.
##' @param printmode the format of the text to return, can be 'LaTeX' (the default) or 'text' for plain text.
##' @param displaymode default is 'console' alternative is 'rstudio'
##' @param ... other arguments passed to the function "format" 
##' @return A text summary, that can be pasted into a LaTeX document and later edited.
##' @export

Summarise <- function(obj,digits=3,scientific=-3,inclIntercept=FALSE,printmode="LaTeX",displaymode="console",...){
    ts <- textSummary(obj,digits=digits,scientific=scientific,inclIntercept=inclIntercept,printmode=printmode,...)

    # if(displaymode=="console"){
    #     oldpar <- par()
    #     par(ask=TRUE)
    # }

    cat("We first look at a global convergence plot, a plot of the log posterior density over the iterations:\n\n")

    plot(obj$tarrec,type="s",xlab="Retained Iteration Number",main="",ylab="log posterior density.")

    cat("We next look at autocorrelations in the latent field:\n\n")

    frailtylag1(obj)

    cat("We next produce trace and autocorrelation plots of the parameters.\n\n")

    cat("Starting with beta:\n\n.")

    par(mfrow=c(1,2))

    for(i in 1:ncol(obj$betasamp)){
        plot(obj$betasamp[,i],main=colnames(obj$betasamp)[i],xlab="Retained Iteration Number",ylab=colnames(obj$betasamp)[i],type="s")
        acf(obj$betasamp[,i],main=colnames(obj$betasamp)[i])
    }

    cat("Next eta:\n\n.")

    par(mfrow=c(1,2))

    for(i in 1:ncol(obj$etasamp)){
        plot(obj$etasamp[,i],main=colnames(obj$etasamp)[i],xlab="Retained Iteration Number",ylab=colnames(obj$etasamp)[i],type="s")
        acf(obj$etasamp[,i],main=colnames(obj$etasamp)[i])
    }

    cat("Next omega:\n\n.")

    par(mfrow=c(1,2))

    for(i in 1:ncol(obj$omegasamp)){
        plot(obj$omegasamp[,i],main=colnames(obj$omegasamp)[i],xlab="Retained Iteration Number",ylab=colnames(obj$omegasamp)[i],type="s")
        acf(obj$omegasamp[,i],main=colnames(obj$omegasamp)[i])
    }

    cat("Next we check the information content in the data by plotting the priors against the posterior for each parameter. The idea is that the data should move the prior somewhat: if the prior and posterior are similar, then there is little information in the data on this particular parameter.")

    priorposterior(obj,pause=FALSE)

    cat("We next produce plots of the spatial covariate-adjusted relative risks.\n\n") 

    par(mfrow=c(1,1))   

    hupp <- hazardexceedance(c(1.25,1.5,2,3))
    hlow <- hazardexceedance(c(0.8,2/3,1/2,1/3),direction="lower")

    mceupp <- MCE(obj,hupp)
    mcelow <- MCE(obj,hlow)

    gr <- getGrid(obj)
    gr$ex1_25 <- mceupp[1,]
    gr$ex1_5 <- mceupp[2,]
    gr$ex2 <- mceupp[3,]
    gr$ex3 <- mceupp[4,]
    gr$ex0_8 <- mcelow[1,]
    gr$ex0_66 <- mcelow[2,]
    gr$ex0_5 <- mcelow[3,]
    gr$ex0_33 <- mcelow[4,]


    cat("A plot of the probability the covariate adjusted hazard exceeds 1.25:\n\n")

    expY <- exp(obj$Ysamp)

    BG <- getBackground(gr)

    if(any(expY>1.25)){
        suppressWarnings(spplot1(gr,"ex1_25",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard exceeds 1.5:\n\n")

    if(any(expY>1.5)){
        suppressWarnings(spplot1(gr,"ex1_5",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard exceeds 2:\n\n")

    if(any(expY>2)){
        suppressWarnings(spplot1(gr,"ex2",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard exceeds 3:\n\n")

    if(any(expY>3)){
        suppressWarnings(spplot1(gr,"ex3",OSMbg=BG,border=NA))
    }



    cat("A plot of the probability the covariate adjusted hazard is below 0.8:\n\n")

    if(any(expY<0.8)){
        suppressWarnings(spplot1(gr,"ex0_8",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard is below 2/3:\n\n")

    if(any(expY<0.66)){
        suppressWarnings(spplot1(gr,"ex0_66",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard is below 1/2:\n\n")

    if(any(expY<0.5)){
        suppressWarnings(spplot1(gr,"ex0_5",OSMbg=BG,border=NA))
    }

    cat("A plot of the probability the covariate adjusted hazard is below 1/3:\n\n")

    if(any(expY<0.33)){
        suppressWarnings(spplot1(gr,"ex0_33",OSMbg=BG,border=NA))
    }

    cat("A plot of the baseline hazard (left) and cumulative hazard (right):")

    par(mfrow=c(1,2))

    baselinehazard(obj)
    baselinehazard(obj,cumulative=TRUE)

    cat("A plot of the posterior spatial covariance function:")

    posteriorcov(obj)

    cat("Having established convergence (??), we summarise the parameter effects in this model.\n\n")

    print(obj)

    print(ts)

    # if(displaymode=="console"){
    #     par <- oldpar
    # }

}