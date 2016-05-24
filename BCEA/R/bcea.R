###INTRO#############################################################################################
## Define Classes & Methods
## v1.0. 4 January, 2012
## v1.1. 14 September, 2012
## v1.2. 17 September 2012
## v1.3-0 June, 2013
## v2.0-1 July, 2013
## v2.0-2 November, 2013
## v2.0-2b February, 2014 - ceac.plot and eib.plot: option comparison included for base graphics
## v2.0-2c July, 2014
## v2.1-0-pre1 AB September, 2014: documentation updated, Smoking dataset and ceef.plot function included, additional modifications
## v2.1.0-pre2 GB October, 2014: modifications to ceef.plot, CreateInputs, struct.psa
## v2.1.0 AB October, 2014: migrated from if(require()) to if(requireNamespace(,quietly=TRUE)); documentation updated
## v2.1.0 AB December, 2014: added threshold argument to ceef.plot function; documentation updated
## v2.1.1 GB+AH April/July 2015: new function for EVPPI using SPDE-INLA; modifications to the EVPPI functions; 
##        documentation updated; allows xlim & ylim in the ceplane.plot, contour and contour2 functions;
##	  it is now possible to run bcea for a scalar wtp; the old evppi function and method has been renamed 
## 	  evppi0, which means there's also a new plot.evppi0 method
## v2.2   GB October 2015: cleaned up and aligned with R's settings. EVPPI function polished up
## v2.2.1 GB+AH October 2015: adds the info-rank plot
## v2.2.2 AB January 2016: minor change to ceef.plot to align with ggplot2 v2.0.0
## (C) Gianluca Baio + contributions by Andrea Berardi, Chris Jackson, Mark Strong & Anna Heath

###Functions included################################################################################
#*> cat(paste0(ls(),"\n"))
# bcea
# bcea.default
# ceac.plot
# ceaf.plot
# ceef.plot
# ceplane.plot
# CEriskav
# CEriskav.default
# contour2
# contour.bcea
# CreateInputs
# eib.plot
# evi.plot
# evppi
# evppi.default
# ib.plot
# info.rank
# mce.plot
# mixedAn
# mixedAn.default
# multi.ce
# plot.bcea
# plot.CEriskav
# plot.evppi
# plot.mixedAn
# sim.table
# struct.psa
# summary.bcea
# summary.mixedAn
# evppi0
# evppi0.default
# plot.evppi0


###bcea##############################################################################################
bcea <- function(e,c,ref=1,interventions=NULL,Kmax=50000,wtp=NULL,plot=FALSE) UseMethod("bcea")

###bcea.default######################################################################################
## Default function
bcea.default <- function(e,c,ref=1,interventions=NULL,Kmax=50000,wtp=NULL,plot=FALSE) {
    ## Compute a Bayesian cost-effectiveness analysis of two or more interventions
    ## INPUTS:
    ## 1. Two objects (e,c). These can be directly computed in a simulation object "sim" from JAGS/BUGS, 
    ##    or derived by postprocessing of "sim" in R. The objects (e,c) have dimension (n.sim x number of 
    ##    interventions) and contain n.sim simulated values for the measures of effectiveness and costs 
    ##    for each intervention being compared. 
    ## 2. The reference intervention as a numeric value. Each intervention is a column in the matrices e 
    ##    and c so if ref=1 the first column is assumed to be associated with the reference intervention. 
    ##    Intervention 1 is assumed the default reference. All others are considered comparators.
    ## 3. A string vector "interventions" including the names of the interventions. If none is provided 
    ##    then labels each as "intervention1",...,"interventionN".
    ## 4. The value Kmax which represents the maximum value for the willingness to pay parameter. If none 
    ##    is provided, then it is assumed Kmax=50000.
    ## 5. A(n optional) vector wtp including the values of the willingness to pay grid. If not specified
    ##    then BCEA will construct a grid of 501 values from 0 to Kmax. This option is useful when 
    ##    performing intensive computations (eg for the EVPPI)
    ##
    ## OUTPUTS:
    ## Graphs & computed values for CE Plane, ICER, EIB, CEAC, EVPI 
    
    # Set the working directory to wherever the user is working, if not externally set
    if(!exists("working.dir")){working.dir <- getwd()}
    
    # Number of simulations & interventions analysed
    n.sim <- dim(e)[1]
    n.comparators <- dim(e)[2]
    
    # Define reference & comparator intervention (different labels can be given here if available!)
    if(is.null(interventions)){interventions <- paste("intervention",1:n.comparators)}
    ints <- 1:n.comparators
    
    # Define intervention i (where i can be a number in [1,...,n.comparators]) as the reference 
    # and the other(s) as comparator(s). Default is the first intervention (first column of e or c)
    comp <- ints[-ref]
    n.comparisons <- n.comparators-1
    
    # Compute Effectiveness & Cost differentials (wrt to reference intervention)
    delta.e <- e[,ref]-e[,comp]
    delta.c <- c[,ref]-c[,comp]
    
    # Compute the ICER
    if(n.comparisons==1) {
        ICER <- mean(delta.c)/mean(delta.e)
    }
    if(n.comparisons>1) {
        ICER <- colMeans(delta.c)/colMeans(delta.e) #apply(delta.c,2,mean)/apply(delta.e,2,mean)
    }
    
    
    # Compute and plot CEAC & EIB
    if(!exists("Kmax")){Kmax<-50000}
    # Lets you select the willingness to pay grid --- useful when doing EVPPI (computationally intensive)
    if (!is.null(wtp)) {
        wtp <- sort(unique(wtp))
        npoints <- length(wtp) - 1
        Kmax <- max(wtp)
        step <- NA
        k <- wtp
        K <- npoints+1
    } else {
        npoints <- 500
        step <- Kmax/npoints
        k <- seq(0,Kmax,step)
        K <- length(k)	
    }
    
    if(n.comparisons==1) {
        ib <- scale(k%*%t(delta.e),delta.c,scale=FALSE)
        ceac <- rowMeans(ib>0) #apply(ib>0,1,mean)
    }
    if(n.comparisons>1) { 
        ib <- array(rep(delta.e, K)*rep(k, each=n.sim*n.comparisons)-as.vector(delta.c),
                    dim=c(n.sim, n.comparisons, K))
        ib <- aperm(ib, c(3,1,2))
        ###          ib <- sweep(apply(delta.e,c(1,2),function(x) k%*%t(x)),c(2,3),delta.c,"-")
        ceac <- apply(ib>0,c(1,3),mean)
    }
    
    # Select the best option for each value of the willingness to pay parameter
    if(n.comparisons==1) {
        eib <- rowMeans(ib)  #apply(ib,1,mean)
        best <- rep(ref,K)
        best[which(eib<0)] <- comp
        ## Finds the k for which the optimal decision changes
        check <- c(0,diff(best))
        kstar <- k[check!=0]
    }
    if(n.comparisons>1) {
        eib <- apply(ib,3,function(x) apply(x,1,mean))
        if (is.null(dim(eib))) {
            tmp <- min(eib)
            tmp2 <- which.min(eib)	
        } else {
            tmp <- apply(eib,1,min)
            tmp2 <- apply(eib,1,which.min)
        }
        best <- ifelse(tmp>0,ref,comp[tmp2])
        # Finds the k for which the optimal decision changes
        check <- c(0,diff(best))
        kstar <- k[check!=0]
    }
    
    # Compute EVPI 
    U <- array(rep(e, K)*rep(k, each=n.sim*n.comparators) - as.vector(c),
               dim=c(n.sim, n.comparators, K))
    U <- aperm(U, c(1,3,2))
    rowMax <- function(x){do.call(pmax, as.data.frame(x))}
    Ustar <- vi <- ol <- matrix(NA,n.sim,K) 
    for (i in 1:K) {
        Ustar[,i] <- rowMax(U[,i,])
        cmd <- paste("ol[,i] <- Ustar[,i] - U[,i,",best[i],"]",sep="")
        eval(parse(text=cmd))     
        vi[,i] <- Ustar[,i] - max(apply(U[,i,],2,mean))
    }
    evi <- colMeans(ol)
    
    ## Outputs of the function
    he <- list(
        n.sim=n.sim,n.comparators=n.comparators,n.comparisons=n.comparisons,delta.e=delta.e,
        delta.c=delta.c,ICER=ICER,Kmax=Kmax,k=k,ceac=ceac,ib=ib,eib=eib,kstar=kstar,
        best=best,U=U,vi=vi,Ustar=Ustar,ol=ol,evi=evi,interventions=interventions,
        ref=ref,comp=comp,step=step,e=e,c=c
    )
    
    class(he) <- "bcea"
    if(plot)
        plot(he)
    return(he)
}



###summary.bcea#################################################################################################
## Summary of the results
summary.bcea <- function(object,wtp=25000,...) {
    
    if(max(object$k)<wtp) {
        wtp <- max(object$k)
        cat(paste("NB: k (wtp) is defined in the interval [",min(object$k)," - ",wtp,"]\n",sep=""))
    }
    if (!is.element(wtp,object$k)) {
        if (!is.na(object$step)) {# The user has selected a non-acceptable value for wtp, but has not specified wtp in the call to bcea
            stop(paste("The willingness to pay parameter is defined in the interval [0-",object$Kmax,
                       "], with increments of ",object$step,"\n",sep=""))
        } else { # The user has actually specified wtp as input in the call to bcea
            tmp <- paste(object$k,collapse=" ")
            stop(paste0("The willingness to pay parameter is defined as:\n[",tmp,"]\nPlease select a suitable value",collapse=" "))
        }
    }
    ind.table <- which(object$k==wtp)
    cols.u <- 1:object$n.comparators
    cols.ustar <- max(cols.u)+1
    cols.ib <- (cols.ustar+1):(cols.ustar+object$n.comparisons)
    cols.ol <- max(cols.ib)+1
    cols.vi <- cols.ol+1
    n.cols <- cols.vi
    
    Table <- matrix(NA,(object$n.sim+1),n.cols)
    Table[1:object$n.sim,cols.u] <- object$U[,ind.table,]
    Table[1:object$n.sim,cols.ustar] <- object$Ustar[,ind.table]
    if(length(dim(object$ib))==2){Table[1:object$n.sim,cols.ib] <- object$ib[ind.table,]}
    if(length(dim(object$ib))>2){Table[1:object$n.sim,cols.ib] <- object$ib[ind.table,,]}
    Table[1:object$n.sim,cols.ol] <- object$ol[,ind.table]
    Table[1:object$n.sim,cols.vi] <- object$vi[,ind.table]
    if(length(dim(object$ib))==2){
        Table[(object$n.sim+1),] <- c(apply(object$U[,ind.table,],2,mean),mean(object$Ustar[,ind.table]),
                                      mean(object$ib[ind.table,]),mean(object$ol[,ind.table]),mean(object$vi[,ind.table]))     
    }
    if(length(dim(object$ib))>2){
        Table[(object$n.sim+1),] <- c(apply(object$U[,ind.table,],2,mean),mean(object$Ustar[,ind.table]),
                                      apply(object$ib[ind.table,,],2,mean),mean(object$ol[,ind.table]),mean(object$vi[,ind.table]))
    }
    
    names.cols <- c(paste("U",seq(1:object$n.comparators),sep=""),"U*",
                    paste("IB",object$ref,"_",object$comp,sep=""),"OL","VI")
    colnames(Table) <- names.cols
    
    tab1 <- matrix(NA,object$n.comparators,1)
    tab1[,1] <- Table[object$n.sim+1,(paste("U",seq(1:object$n.comparators),sep=""))]
    colnames(tab1) <- "Expected utility"
    rownames(tab1) <- object$interventions
    
    tab2 <- matrix(NA,object$n.comparisons,3)
    tab2[,1] <- Table[object$n.sim+1,paste("IB",object$ref,"_",object$comp,sep="")]
    if (object$n.comparisons==1) {
        tab2[,2] <- sum(Table[1:object$n.sim,paste("IB",object$ref,"_",object$comp,sep="")]>0)/object$n.sim
        tab2[,3] <- object$ICER
    }
    if (object$n.comparisons>1) {
        for (i in 1:object$n.comparisons) {
            tab2[i,2] <- sum(Table[1:object$n.sim,paste("IB",object$ref,"_",object$comp[i],sep="")]>0)/object$n.sim
            tab2[i,3] <- object$ICER[i]
        }
    }
    colnames(tab2) <- c("EIB","CEAC","ICER")
    rownames(tab2) <- paste(object$interventions[object$ref]," vs ",object$interventions[object$comp],sep="")
    
    tab3 <- matrix(NA,1,1)
    tab3[,1] <- Table[object$n.sim+1,"VI"]
    rownames(tab3) <- "EVPI"
    colnames(tab3) <- ""
    
    ## Prints the summary table
    cat("\n")
    cat("Cost-effectiveness analysis summary \n")
    cat("\n")
    cat(paste("Reference intervention:  ",object$interventions[object$ref],"\n",sep=""))
    if(object$n.comparisons==1) {
        text.temp <- paste("Comparator intervention: ",object$interventions[object$comp],"\n",sep="")
        cat(text.temp)
    }
    
    if(object$n.comparisons>1) {
        text.temp <- paste("Comparator intervention(s): ",object$interventions[object$comp[1]],"\n",sep="")
        cat(text.temp)
        for (i in 2:object$n.comparisons) {
            cat(paste("                          : ",object$interventions[object$comp[i]],"\n",sep=""))
        }
    }
    cat("\n")
    if(length(object$kstar)==0 & !is.na(object$step)){
        cat(paste(object$interventions[object$best[1]]," dominates for all k in [",
                  min(object$k)," - ",max(object$k),"] \n",sep=""))
    }
    if(length(object$kstar)==1 & !is.na(object$step)){
        cat(paste("Optimal decision: choose ",object$interventions[object$best[object$k==object$kstar-object$step]],
                  " for k<",object$kstar," and ",object$interventions[object$best[object$k==object$kstar]],
                  " for k>=",object$kstar,"\n",sep=""))
    }
    if(length(object$kstar)>1 & !is.na(object$step)){
        cat(paste("Optimal decision: choose ",object$interventions[object$best[object$k==object$kstar[1]-object$step]],
                  " for k < ",object$kstar[1],"\n",sep=""))
        for (i in 2:length(object$kstar)) {
            cat(paste("                         ",object$interventions[object$best[object$k==object$kstar[i]-object$step]],
                      " for ",object$kstar[i-1]," <= k < ",object$kstar[i],"\n",sep=""))
        }
        cat(paste("                         ",object$interventions[object$best[object$k==object$kstar[length(object$kstar)]]],
                  " for k >= ",object$kstar[length(object$kstar)],"\n",sep=""))
    }
    cat("\n\n")
    cat(paste("Analysis for willingness to pay parameter k = ",wtp,"\n",sep=""))
    cat("\n")
    print(tab1,quote=F,digits=5,justify="center")
    cat("\n")
    print(tab2,quote=F,digits=5,justify="center")
    cat("\n")
    cat(paste("Optimal intervention (max expected utility) for k=",wtp,": ",
              object$interventions[object$best][object$k==wtp],"\n",sep=""))
    print(tab3,quote=F,digits=5,justify="center")
}


###sim.table##################################################################################################
# Produce a summary table with the results of simulations for the health economic variables of interest
sim.table <- function(he,wtp=25000) {
    
    if(wtp>he$Kmax){wtp=he$Kmax}
    
    if (!is.element(wtp,he$k)) {
        if (!is.na(he$step)) {# The user has selected a non-acceptable value for wtp, but has not specified wtp in the call to bcea
            stop(paste("The willingness to pay parameter is defined in the interval [0-",he$Kmax,
                       "], with increments of ",he$step,"\n",sep=""))
        } else { # The user has actually specified wtp as input in the call to bcea
            tmp <- paste(he$k,collapse=" ")
            stop(paste0("The willingness to pay parameter is defined as:\n[",tmp,"]\nPlease select a suitable value",collapse=" "))
        }
    }
    
    ind.table <- which(he$k==wtp)
    cols.u <- 1:he$n.comparators
    cols.ustar <- max(cols.u)+1
    cols.ib <- (cols.ustar+1):(cols.ustar+he$n.comparisons)
    cols.ol <- max(cols.ib)+1
    cols.vi <- cols.ol+1
    n.cols <- cols.vi
    
    Table <- matrix(NA,(he$n.sim+1),n.cols)
    Table[1:he$n.sim,cols.u] <- he$U[,ind.table,]
    Table[1:he$n.sim,cols.ustar] <- he$Ustar[,ind.table]
    if(length(dim(he$ib))==2){Table[1:he$n.sim,cols.ib] <- he$ib[ind.table,]}
    if(length(dim(he$ib))>2){Table[1:he$n.sim,cols.ib] <- he$ib[ind.table,,]}
    Table[1:he$n.sim,cols.ol] <- he$ol[,ind.table]
    Table[1:he$n.sim,cols.vi] <- he$vi[,ind.table]
    if(length(dim(he$ib))==2){
        Table[(he$n.sim+1),] <- c(apply(he$U[,ind.table,],2,mean),mean(he$Ustar[,ind.table]),
                                  mean(he$ib[ind.table,]),mean(he$ol[,ind.table]),mean(he$vi[,ind.table]))	
    }
    if(length(dim(he$ib))>2){
        Table[(he$n.sim+1),] <- c(apply(he$U[,ind.table,],2,mean),mean(he$Ustar[,ind.table]),
                                  apply(he$ib[ind.table,,],2,mean),mean(he$ol[,ind.table]),mean(he$vi[,ind.table]))
    }
    
    names.cols <- c(paste("U",seq(1:he$n.comparators),sep=""),"U*",paste("IB",he$ref,"_",he$comp,sep=""),"OL","VI")
    colnames(Table) <- names.cols
    rownames(Table) <- c(1:he$n.sim,"Average")
    
    ## Outputs of the function
    list(Table=Table,names.cols=names.cols,wtp=wtp,ind.table=ind.table)
}


##########################################ceplane.plot########################################################
## Plots the CE Plane
ceplane.plot <- function(he,comparison=NULL,wtp=25000,pos=c(1,1),
                         size=NULL,graph=c("base","ggplot2"),
                         xlim=NULL,ylim=NULL,...) {
    ### hidden options for ggplot2 ###
    # ICER.size =                    # changes ICER point size
    # label.pos = FALSE              # uses alternate position for wtp label (old specification)
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    alt.legend <- pos
    
    # Forces R to avoid scientific format for graphs labels
    options(scipen=10)
    
    # Additional/optional arguments
    exArgs <- list(...)
    if(!exists("xlab",where=exArgs)){xlab <- "Effectiveness differential"} else {xlab <- exArgs$xlab}
    if(!exists("ylab",where=exArgs)){ylab <- "Cost differential"} else {ylab <- exArgs$ylab}
    if(!exists("ICER.col",where=exArgs)){ICER.col <- "red"} else {ICER.col <- exArgs$ICER.col}
    if(!exists("title",where=exArgs)){title <- paste("Cost effectiveness plane \n",he$interventions[he$ref]," vs ",he$interventions[he$comp],sep="")} 
    else {title <- exArgs$title}
    
    if(base.graphics) {
        if(!is.null(size))
            message("option size will be ignored using base graphics")
        
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==0)
                temp <- paste0(temp,"bottom")
            else
                temp <- paste0(temp,"top")
            if(alt.legend[1]==0)
                temp <- paste0(temp,"left")
            else
                temp <- paste0(temp,"right")
            alt.legend <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="topright"
            else
                alt.legend="topleft"
        }
        
        # Encodes characters so that the graph can be saved as ps or pdf
        ps.options(encoding="CP1250")
        pdf.options(encoding="CP1250")
        
        if(he$n.comparisons==1) {
            m.e <- range(he$delta.e)[1]
            M.e <- range(he$delta.e)[2]
            m.c <- range(he$delta.c)[1]
            M.c <- range(he$delta.c)[2]
            step <- (M.e-m.e)/10
            m.e <- ifelse(m.e<0,m.e,-m.e)
            m.c <- ifelse(m.c<0,m.c,-m.c)
            x.pt <- .95*m.e
            y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
            xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
            yy <- xx*wtp
            xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
            yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
            xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
            if(!is.null(xlim)) {m.e <- xlim[1]; M.e <- xlim[2]}
            if(!is.null(ylim)) {m.c <- ylim[1]; M.c <- ylim[2]}
            plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c),
                 xlab=xlab,ylab=ylab,main=title,axes=F)
            polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
                    c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
                    col="grey95",border="black")
            #  polygon(c(xx,xx),c(yy,rev(yy)),col="grey95",border="black")
            axis(1); axis(2); box()
            points(he$delta.e,he$delta.c,pch=20,cex=.35,col="grey55")
            abline(h=0,col="dark grey")
            abline(v=0,col="dark grey")
            text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER,digits=6,nsmall=2),sep=""),cex=.95,pos=2,col=ICER.col)
            points(mean(he$delta.e),mean(he$delta.c),pch=20,col=ICER.col,cex=1)
            t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
            text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
        }
        if(he$n.comparisons>1 & is.null(comparison)==TRUE) { 
            if(!exists("title",where=exArgs)){title <- "Cost-effectiveness plane"} else {title <- exArgs$title}
            cl <- colors()
            color <- cl[floor(seq(262,340,length.out=he$n.comparators))]  # gray scale
            if(is.null(xlim)) {xlim <- range(he$delta.e)}
            if(is.null(ylim)) {ylim <- range(he$delta.c)}
            plot(he$delta.e[,1],he$delta.c[,1],pch=20,cex=.35,xlim=xlim,ylim=ylim,
                 xlab=xlab,ylab=ylab,main=title)
            for (i in 2:he$n.comparisons) {
                points(he$delta.e[,i],he$delta.c[,i],pch=20,cex=.35,col=color[i])
            }
            abline(h=0,col="dark grey")
            abline(v=0,col="dark grey")
            text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
            legend(alt.legend,text,col=color,cex=.7,bty="n",lty=1)
        }
        if(he$n.comparisons>1 & is.null(comparison)==FALSE & length(comparison)==1) { 
            if(!exists("title",where=exArgs)){title <- paste("Cost effectiveness plane \n",
                                                             he$interventions[he$ref]," vs ",he$interventions[he$comp[comparison]],sep="")} else {title <- exArgs$title}
            m.e <- range(he$delta.e[,comparison])[1]
            M.e <- range(he$delta.e[,comparison])[2]
            m.c <- range(he$delta.c[,comparison])[1]
            M.c <- range(he$delta.c[,comparison])[2]
            step <- (M.e-m.e)/10
            m.e <- ifelse(m.e<0,m.e,-m.e)
            m.c <- ifelse(m.c<0,m.c,-m.c)
            x.pt <- .95*m.e
            y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
            xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
            yy <- xx*wtp
            xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
            yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
            xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
            if(!is.null(xlim)) {m.e <- xlim[1]; M.e <- xlim[2]}
            if(!is.null(ylim)) {m.c <- ylim[1]; M.c <- ylim[2]}
            plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c),
                 xlab=xlab,ylab=ylab,main=title,axes=F)
            polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
                    c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
                    col="grey95",border="black")
            axis(1); axis(2); box()
            points(he$delta.e[,comparison],he$delta.c[,comparison],pch=20,cex=.35,col="grey55")
            abline(h=0,col="dark grey")
            abline(v=0,col="dark grey")
            text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER[comparison],digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
            points(mean(he$delta.e[,comparison]),mean(he$delta.c[,comparison]),pch=20,col="red",cex=1)
            t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
            text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
        }
        if(he$n.comparisons>1&is.null(comparison)==FALSE&length(comparison)!=1) {
            stopifnot(all(comparison %in% 1:he$n.comparisons))
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            return(ceplane.plot(he,wtp=wtp,pos=alt.legend,graph="base",size=size,...))
        }
    } #if(base.graphics)
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE) & requireNamespace("grid",quietly=TRUE))){
            message("Falling back to base graphics\n")
            ceplane.plot(he,comparison=comparison,wtp=wtp,pos=alt.legend,graph="base"); return(invisible(NULL))
        }
        # no visible binding note
        
        if(isTRUE(requireNamespace("ggplot2",quietly=TRUE) & requireNamespace("grid",quietly=TRUE))){
            
            delta.e <- delta.c <- lambda.e <- lambda.c <- NULL
            
            if(is.null(size)) {size <- ggplot2::rel(3.5)}
            
            label.pos <- TRUE
            opt.theme <- ggplot2::theme()
            ICER.size <- ifelse(he$n.comparisons==1,2,0)
            exArgs <- list(...)
            if(length(exArgs)>=1){
                if(exists("ICER.size",where=exArgs))
                    ICER.size <- exArgs$ICER.size
                if(exists("label.pos",where=exArgs))
                    if(is.logical(exArgs$label.pos))
                        label.pos <- exArgs$label.pos
                    for(obj in exArgs)
                        if(ggplot2::is.theme(obj))
                            opt.theme <- opt.theme + obj
            }
            
            if(he$n.comparisons==1) {
                kd <- data.frame(he$delta.e,he$delta.c)
                names(kd) <- c("delta.e","delta.c")
                # for scale_x_continuous(oob=)
                do.nothing=function(x,limits) return(x)
                # plot limits
                range.e <- range(kd$delta.e)
                range.c <- range(kd$delta.c)
                range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
                range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
                # ce plane data
                x1 <- range.e[1]-2*abs(diff(range.e))
                x2 <- range.e[2]+2*abs(diff(range.e))
                x3 <- x2
                x <- c(x1,x2,x3)
                y <- x*wtp; y[3] <- x1*wtp
                plane <- data.frame(x=x,y=y)
                
                # build a trapezoidal plane instead of a triangle if the y value is less than the minimum difference on costs
                if(y[1]>1.2*range.c[1]) {
                    plane <- rbind(plane,
                                   c(x2,2*range.c[1]), #new bottom-right vertex
                                   c(x1,2*range.c[1])) #new bottom-left vertex
                }
                
                # actual plot
                ceplane <- ggplot2::ggplot(kd, ggplot2::aes(delta.e,delta.c)) +
                    ggplot2::theme_bw() +
                    ggplot2::scale_x_continuous(limits=range.e,oob=do.nothing) + 
                    ggplot2::scale_y_continuous(limits=range.c,oob=do.nothing) +
                    ggplot2::scale_color_manual("",labels=paste0("ICER = ",format(he$ICER,digits=6,nsmall=2),"  "),values="red") +     
                    ggplot2::geom_line(data=plane[1:2,],ggplot2::aes(x=x,y=y),color="black",linetype=1) +
                    ggplot2::geom_polygon(data=plane,ggplot2::aes(x=x,y=y),fill="light gray",alpha=.3) +
                    ggplot2::geom_hline(ggplot2::aes(yintercept=0),colour="grey") + 
                    ggplot2::geom_vline(ggplot2::aes(xintercept=0),colour="grey") +
                    ggplot2::geom_point(size=1,colour="grey33") +
                    ggplot2::geom_point(ggplot2::aes(mean(delta.e),mean(delta.c),color=as.factor(1)),size=ICER.size)
                
                if(!label.pos) {
                    # moves the wtp label depending on whether the line crosses the y-axis
                    ceplane <- ceplane + ggplot2::annotate(geom="text",x=ifelse(range.c[1]/wtp>range.e[1],range.c[1]/wtp,range.e[1]),
                                                           #y=(1+.085)*range.c[1],
                                                           y=range.c[1],
                                                           label=paste0("k = ",format(wtp,digits=6)),hjust=-.15,size=size)
                }
                else{
                    m.e <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
                    m.c <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
                    x.pt <- .95*m.e
                    y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
                    ceplane <- ceplane + ggplot2::annotate(geom="text",x=x.pt,y=y.pt,
                                                           label=paste0("k = ",format(wtp,digits=6)),hjust=.15,size=size)
                }
            }
            
            if(he$n.comparisons>1&is.null(comparison)==TRUE) {
                
                # create dataframe for plotting
                kd <- with(he,data.frame(c(delta.e),c(delta.c)))
                names(kd) <- c("delta.e","delta.c")
                kd$comparison <- as.factor(sort(rep(1:he$n.comparisons,dim(he$delta.e)[1])))
                
                # dataset for ICERs
                means <- matrix(NA_real_,nrow=he$n.comparisons,ncol=2)
                for(i in 1:he$n.comparisons)
                    means[i,] <- colMeans(kd[kd$comparison==i,-3])
                means <- data.frame(means)
                means$comparison <- factor(1:he$n.comparisons)
                names(means) <- c("lambda.e","lambda.c","comparison")
                
                # labels for legend
                comparisons.label <- with(he,paste0(interventions[ref]," vs ",interventions[comp]))
                # vector of values for color, take out white, get integer values
                colors.label <- with(he,paste0("gray",round(seq(0,100,length.out=(n.comparisons+1))[-(n.comparisons+1)])))
                
                # polygon
                do.nothing=function(x,limits) return(x)
                # plot limits
                range.e <- range(kd$delta.e)
                range.c <- range(kd$delta.c)
                range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
                range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
                # ce plane data
                x1 <- range.e[1]-2*abs(diff(range.e))
                x2 <- range.e[2]+2*abs(diff(range.e))
                x3 <- x2
                x <- c(x1,x2,x3)
                y <- x*wtp; y[3] <- x1*wtp
                plane <- data.frame(x=x,y=y,comparison=factor(rep(he$n.comparisons+1,3)))
                
                # build a trapezoidal plane instead of a triangle if the y value is less than the minimum difference on costs
                if(y[1]>min(kd$delta.c)) {
                    plane <- rbind(plane,
                                   c(x2,2*min(kd$delta.c),he$n.comparisons+1), #new bottom-right vertex
                                   c(x1,2*min(kd$delta.c),he$n.comparisons+1)) #new bottom-left vertex
                }
                
                ceplane <-
                    ggplot2::ggplot(kd,ggplot2::aes(x=delta.e,y=delta.c,col=comparison)) +
                    ggplot2::theme_bw() +
                    ggplot2::scale_color_manual(labels=comparisons.label,values=colors.label,na.value="black") +
                    ggplot2::scale_x_continuous(limits=range.e,oob=do.nothing) +
                    ggplot2::scale_y_continuous(limits=range.c,oob=do.nothing) +
                    ggplot2::annotate("line",x=plane[1:2,1],y=plane[1:2,2],colour="black") +
                    ggplot2::annotate("polygon",plane$x,plane$y,fill="light grey",alpha=.3) +
                    ggplot2::geom_hline(ggplot2::aes(yintercept=0),colour="grey") + ggplot2::geom_vline(ggplot2::aes(xintercept=0),colour="grey") +
                    ggplot2::geom_point(size=1)+
                    ggplot2::geom_point(data=means,ggplot2::aes(x=lambda.e,y=lambda.c),colour="red",size=ICER.size)
                
                # wtp label
                if(!label.pos){
                    ceplane <- ceplane +
                        ggplot2::annotate(geom="text",
                                          x=ifelse(range.c[1]/wtp>range.e[1],range.c[1]/wtp,range.e[1]),
                                          y=range.c[1],
                                          label=paste0("k = ",format(wtp,digits=6),"  "),hjust=.15,size=size
                        )
                }
                else{
                    m.e <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
                    m.c <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
                    x.pt <- .95*m.e
                    y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
                    ceplane <- ceplane + ggplot2::annotate(geom="text",x=x.pt,y=y.pt,
                                                           label=paste0("k = ",format(wtp,digits=6)),hjust=.15,size=size)
                }
            }
            
            if(he$n.comparisons>1&is.null(comparison)==FALSE) {
                # adjusts bcea object for the correct number of dimensions and comparators
                he$comp <- he$comp[comparison]
                he$delta.e <- he$delta.e[,comparison]
                he$delta.c <- he$delta.c[,comparison]
                he$n.comparators=length(comparison)+1
                he$n.comparisons=length(comparison)
                he$interventions=he$interventions[sort(c(he$ref,he$comp))]
                he$ICER=he$ICER[comparison]
                he$ib=he$ib[,,comparison]
                he$eib=he$eib[,comparison]
                he$U=he$U[,,sort(c(he$ref,comparison+1))]
                he$ceac=he$ceac[,comparison]
                he$ref=rank(c(he$ref,he$comp))[1]
                he$comp=rank(c(he$ref,he$comp))[-1]
                he$mod <- TRUE #
                
                return(ceplane.plot(he,wtp=wtp,pos=alt.legend,graph="ggplot2",size=size,...))
            }
            
            if(!exists("title",where=exArgs)) {
                labs.title <- "Cost-Effectiveness Plane"
                labs.title <- with(he,paste0(labs.title,
                                             ifelse(n.comparisons==1,
                                                    paste0("\n",interventions[ref]," vs ",interventions[-ref]),
                                                    paste0(
                                                        ifelse(isTRUE(he$mod),
                                                               paste0("\n",interventions[ref]," vs ",
                                                                      paste0(interventions[comp],collapse=", ")),
                                                               "")))))
            } else {labs.title <- exArgs$title}
            
            ceplane <- ceplane + ggplot2::labs(title=labs.title,x=xlab,y=ylab)
            
            jus <- NULL
            if(isTRUE(alt.legend)) {
                alt.legend="bottom"
                ceplane <- ceplane + ggplot2::theme(legend.direction="vertical")
            }
            else{
                if(is.character(alt.legend)) {
                    choices <- c("left", "right", "bottom", "top")
                    alt.legend <- choices[pmatch(alt.legend,choices)]
                    jus="center"
                    if(is.na(alt.legend))
                        alt.legend=FALSE
                }
                if(length(alt.legend)>1)
                    jus <- alt.legend
                if(length(alt.legend)==1 & !is.character(alt.legend)) {
                    alt.legend <- c(1,1)
                    jus <- alt.legend
                }
            }
            
            ceplane <- ceplane + 
                ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),legend.background=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold")) +
                ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),legend.text.align=0) +
                ggplot2::theme(plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3)) +
                opt.theme
            
            if(he$n.comparisons==1)
                ceplane  <- ceplane + ggplot2::theme(legend.key.size=grid::unit(.1,"lines")) + opt.theme
            
            return(ceplane)
        }
    } # !base.graphics
    
}

###ib.plot####################################################################################################
## Plots the IB
ib.plot <- function(he,comparison=NULL,wtp=25000,bw=nbw,n=512,xlim=NULL,graph=c("base","ggplot2")){
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE)
    # comparison controls which comparator is used when more than 2 interventions are present
    # bw and n control the level of smoothness of the kernel density estimation
    options(scipen=10)
    
    if(!is.null(comparison))
        stopifnot(comparison<=he$n.comparison)
    
    if(base.graphics) {
        if(max(he$k)<wtp) {
            wtp <- max(he$k)
            cat(paste("NB: k (wtp) is defined in the interval [",min(he$k)," - ",wtp,"]\n",sep=""))
        }
        
        if (!is.element(wtp,he$k)) {
            if (!is.na(he$step)) {# The user has selected a non-acceptable value for wtp, but has not specified wtp in the call to bcea
                stop(paste("The willingness to pay parameter is defined in the interval [0-",he$Kmax,
                           "], with increments of ",he$step,"\n",sep=""))
            } else { # The user has actually specified wtp as input in the call to bcea
                tmp <- paste(he$k,collapse=" ")
                stop(paste0("The willingness to pay parameter is defined as:\n  [",tmp,"]\n  Please select a suitable value",collapse=" "))
            }
        }
        
        w <- which(he$k==wtp)
        if(he$n.comparisons==1) {
            nbw <- sd(he$ib[w,])/1.5
            d <- density(he$ib[w,],bw=bw,n=n)
            txt <- paste("Incremental Benefit distribution\n",he$interventions[he$ref],
                         " vs ",he$interventions[he$comp],sep="")
        }
        if(he$n.comparisons>1) {
            if(is.null(comparison)){
                comparison <- 1
            }
            nbw <- sd(he$ib[w,,comparison])/1.5
            d <- density(he$ib[w,,comparison],bw=bw,n=n)
            txt <- paste("Incremental Benefit distribution\n",he$interventions[he$ref],
                         " vs ",he$interventions[he$comp[comparison]],sep="")
        }
        if(is.null(xlim)==TRUE){
            xlim<-range(d$x)
        }
        plot(d$x,d$y,t="l",ylab="Density",xlab=expression(paste("IB(",bold(theta),")",sep="")),main=txt,axes=F,col="white",xlim=xlim)
        box()
        axis(1)
        ypt <- .95*max(d$y)
        xpt <- d$x[max(which(d$y>=ypt))]
        text(xpt,ypt,parse(text=paste("p(IB(",expression(bold(theta)),")>0,k==",format(wtp,digits=8,nsmall=2),")",sep="")),
             cex=.85,pos=4)
        xplus <- d$x[d$x>=0]
        yplus <- d$y[d$x>=0]
        polygon(c(0,xplus),c(0,yplus),density=20,border="white")
        points(d$x,d$y,t="l")
        abline(v=0,col="black")
        abline(h=0,col="black")
    }                   # ! base graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            message("falling back to base graphics\n")
            ib.plot(he,comparison=comparison,wtp=wtp,bw=bw,n=n,xlim=xlim,graph="base")
            return(invisible(NULL))
        }
        
        ### no visible binding note
        x <- y <- NA_real_
        
        if(is.null(comparison)) {
            comparison <- 1
        }
        
        if(max(he$k)<wtp) {
            wtp <- max(he$k)
            message(paste0("NB: k (wtp) is defined in the interval [",min(he$k)," - ",wtp,"]\n"))
        }
        
        if (!is.element(wtp,he$k)) {
            if (!is.na(he$step)) {# The user has selected a non-acceptable value for wtp, but has not specified wtp in the call to bcea
                stop(paste("The willingness to pay parameter is defined in the interval [0-",he$Kmax,
                           "], with increments of ",he$step,"\n",sep=""))
            } else { # The user has actually specified wtp as input in the call to bcea
                tmp <- paste(he$k,collapse=" ")
                stop(paste0("The willingness to pay parameter is defined as:\n  [",tmp,"]\n  Please select a suitable value",collapse=" "))
            }
        }
        
        w <- which(he$k==wtp)
        if(he$n.comparisons==1) {
            nbw <- sd(he$ib[w,])/1.5
            density <- density(he$ib[w,],bw=bw,n=n)
            df <- data.frame("x"=density$x,"y"=density$y)
        }
        if(he$n.comparisons>1) {
            nbw <- sd(he$ib[w,,comparison])/1.5
            density <- density(he$ib[w,,comparison],bw=bw,n=n)
            df <- data.frame("x"=density$x,"y"=density$y)
        }
        if(is.null(xlim)){
            xlim<-range(df$x)
        }
        ib <- ggplot2::ggplot(df,ggplot2::aes(x,y)) + ggplot2::theme_bw() +
            ggplot2::geom_vline(xintercept=0,colour="grey50",size=0.5) + ggplot2::geom_hline(yintercept=0,colour="grey50",size=0.5) +
            ggplot2::geom_ribbon(data=subset(df,x>0),ggplot2::aes(ymax=y),ymin=0,fill="grey50",alpha=.2) + ggplot2::geom_line() +
            ggplot2::annotate(geom="text",label=paste0("p(IB(theta)>0,k==",wtp,")"),parse=T,x=df$x[which.max(df$y)],y=max(df$y),hjust=-.5,vjust=1,size=3.5) + ggplot2::coord_cartesian(xlim=xlim)
        
        labs.title <- paste0("Incremental Benefit Distribution\n",he$interventions[he$ref]," vs ",
                             he$interventions[he$comp[comparison]],"")
        
        ib <- ib + 
            ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),text=ggplot2::element_text(size=11),
                           panel.grid=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank()) +
            ggplot2::labs(title=labs.title,x=parse(text="IB(theta)"),y="Density") +
            ggplot2::theme(plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        return(ib)
        
    } #! base.graphics
}


###eib.plot###################################################################################################
## Plots the EIB
eib.plot <- function(he,comparison=NULL,pos=c(1,0),size=NULL,plot.cri=NULL,graph=c("base","ggplot2"),...) {
    
    options(scipen=10)
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    ### evaluate arguments. possibility to include different values of confidence as "alpha"
    exArgs <- list(...)
    alpha <- 0.05
    cri.quantile <- TRUE
    if(length(exArgs)>=1){
        if(exists("cri.quantile",where=exArgs))
            cri.quantile <- exArgs$cri.quantile
        if(exists("alpha",where=exArgs)){
            alpha <- exArgs$alpha
            if(alpha<0 | alpha>1) {
                warning("Argument alpha must be between 0 and 1. Reset to default at 0.95")
                alpha <- 0.05
            }
            if(alpha>0.80 & cri.quantile) {
                warning("It is recommended adopting the normal approximation of the credible interval for high values of alpha. Please set the argument cri.quantile=FALSE to use the normal approsimation.")
            }
        }
    }
    
    ### function to calculate the credible intervals
    eib.plot.cri <- function(he,alpha,cri.quantile) {
        if(alpha<0 | alpha>1) {
            warning("Argument alpha must be between 0 and 1. Reset to default at 0.95")
            alpha <- 0.05
        }
        margin <- 1
        if(he$n.comparison>1) margin <- c(1,3)
        cri <- data.frame("low"=c(apply(he$ib,margin,function(x) ifelse(cri.quantile,
                                                                        quantile(x,(alpha)/2),
                                                                        mean(x)-qnorm((alpha)/2)*sd(x)
        ))),
        "upp"=c(apply(he$ib,margin,function(x) ifelse(cri.quantile,
                                                      quantile(x,1-(alpha)/2),
                                                      mean(x)-qnorm(1-(alpha)/2)*sd(x)
        ))),
        "comp"=as.factor(sort(rep(1:he$n.comparisons,length(he$k)))))
        return(cri)
    }
    ### if plot.cri is null, if comp=1 plot them otherwise do not (clutter!)
    if(is.null(plot.cri) & isTRUE(he$n.comparisons==1 | is.null(comparison)))
        plot.cri <- ifelse(he$n.comparisons==1,TRUE,FALSE)
    
    ### calculate credible intervals if necessary
    if(isTRUE(plot.cri))
        cri <- eib.plot.cri(he,alpha,cri.quantile)
    
    ### calculate plot vertical limits
    yl <- ifelse(rep(!isTRUE(plot.cri),2),
                 range(c(he$eib)),
                 range(c(he$eib),c(cri[,1:2])))
    
    if(base.graphics) {
        if(!is.null(size)){
            if(!is.na(size)){
                message("Option size will be ignored using base graphics.")
                size <- NULL
            }
        }
        
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==0)
                temp <- paste0(temp,"bottom")
            else
                temp <- paste0(temp,"top")
            if(alt.legend[1]==1)
                temp <- paste0(temp,"right")
            else
                temp <- paste0(temp,"left")
            alt.legend <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="topleft"
            else
                alt.legend="topright"
        }
        
        if(he$n.comparisons==1) {
            plot(NULL,xlab="Willingness to pay",ylab="EIB",ylim=yl,xlim=range(he$k),
                 main=paste0("Expected Incremental Benefit",ifelse(plot.cri,paste0("\nand ",format((1-alpha)*100,digits=4),"% credible intervals"),"")))
            ### x axis
            abline(h=0,col="grey")
            ### EIB
            lines(he$k,he$eib)
            ### CRI
            if(plot.cri){
                lines(he$k,cri$low,col="grey50",lty=2)
                lines(he$k,cri$upp,col="grey50",lty=2)
            }
            ### BEP
            if(length(he$kstar)>0 & is.null(size)) {
                abline(v=he$kstar,col="dark grey",lty="dotted")
                text(he$kstar,min(yl),paste("k* = ",he$kstar,sep=""))
            }
            if(isTRUE(he$mod))
                legend(alt.legend,paste0(he$interventions[he$ref]," vs ",he$interventions[he$comp]),cex=.7,bty="n",lty=1,lwd=1)
        }
        if(he$n.comparisons>1&is.null(comparison)) {
            lwd <- ifelse(he$n.comparisons>6,1.5,1)
            plot(NULL,xlab="Willingness to pay", ylab="EIB",ylim=yl,xlim=range(he$k),
                 main=paste0("Expected Incremental Benefit",ifelse(plot.cri,paste0("\nand ",format((1-alpha)*100,digits=4),"% credible intervals"),"")))
            abline(h=0,col="grey")
            for (j in 1:he$n.comparisons){
                lines(he$k,he$eib[,j],lty=j,lwd=ifelse(plot.cri,lwd+1,lwd))
                if(plot.cri){
                    lines(he$k,cri$low[cri$comp==j],lty=j,lwd=lwd,col="grey50")
                    lines(he$k,cri$upp[cri$comp==j],lty=j,lwd=lwd,col="grey50")
                }
            }
            if(length(he$kstar)>0 & is.null(size)) {
                abline(v=he$kstar,col="dark grey",lty="dotted")
                text(he$kstar,min(yl),paste("k* = ",he$kstar,sep=""))
            }
            text <- paste0(he$interventions[he$ref]," vs ",he$interventions[he$comp])
            legend(alt.legend,text,cex=.7,bty="n",lty=1:he$n.comparisons,lwd=ifelse(plot.cri,lwd+1,lwd))
        }
        if(he$n.comparisons>1&!is.null(comparison)) {
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            eib.plot(he,pos=alt.legend,graph="base",size=size,comparison=NULL,plot.cri=plot.cri,alpha=alpha,cri.quantile=cri.quantile,...)
        }
    } # base.graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("falling back to base graphics\n")
            eib.plot(he,pos=alt.legend,graph="base"); return(invisible(NULL))
        }
        
        ### no visible binding note
        k <- kstar <- low <- upp <- NA_real_
        
        if(is.null(size))
            size=ggplot2::rel(3.5)
        
        opt.theme <- ggplot2::theme()
        for(obj in exArgs)
            if(ggplot2::is.theme(obj))
                opt.theme <- opt.theme + obj
        
        if(he$n.comparisons==1) {
            # data frame
            data.psa <- with(he,data.frame("k"=k,"eib"=eib,"comparison"=as.factor(sort(rep(1:n.comparisons,length(he$k))))))
            if(plot.cri)
                data.psa <- cbind(data.psa,cri)
            eib <- ggplot2::ggplot(data.psa, ggplot2::aes(k,eib)) +
                ggplot2::theme_bw() +
                ggplot2::geom_hline(ggplot2::aes(yintercept=0),colour="grey")
            if(!isTRUE(he$mod)){
                eib <- eib + ggplot2::geom_line()
            }
            else{
                eib <- eib + ggplot2::geom_line(ggplot2::aes(linetype=comparison)) +
                    ggplot2::scale_linetype_manual("",values=1,labels=with(he,paste0(interventions[ref]," vs ",interventions[comp])))
            }
            
            if(!length(he$kstar)==0 & !is.na(size)) {
                # label
                label <- paste0("k* = ",format(he$kstar,digits=6))
                eib <- eib +
                    ggplot2::geom_vline(ggplot2::aes(xintercept=kstar),data=data.frame("kstar"=he$kstar),colour="grey50",linetype=2,size=.5) +
                    ggplot2::annotate("text",label=label,x=he$kstar,y=min(yl),hjust=ifelse((max(he$k)-he$kstar)/max(he$k)>1/6,-.1,1.1),size=size)
            }
            
            if(plot.cri){
                eib <- eib +
                    ggplot2::geom_line(ggplot2::aes(y=low),colour="grey50",lty=2) +
                    ggplot2::geom_line(ggplot2::aes(y=upp),colour="grey50",lty=2)
            }
        }
        if(he$n.comparisons>1&is.null(comparison)==TRUE) {
            data.psa <- with(he,data.frame("k"=c(k),"eib"=c(eib),"comparison"=as.factor(sort(rep(1:n.comparisons,length(he$k))))))
            if(plot.cri)
                data.psa <- cbind(data.psa,cri)
            # labels for legend
            comparisons.label <- with(he,paste0(interventions[ref]," vs ",interventions[comp]))
            
            # linetype is the indicator of the comparison.
            # 1 = solid, 2 = dashed, 3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash
            linetypes <- rep(c(1,2,3,4,5,6),ceiling(he$n.comparisons/6))[1:he$n.comparisons]
            
            eib <- 
                ggplot2::ggplot(data.psa,ggplot2::aes(x=k,y=eib,linetype=comparison)) + 
                ggplot2::geom_hline(yintercept=0,linetype=1,color="grey") + 
                ggplot2::theme_bw() +
                ggplot2::geom_line(lwd=ifelse(!plot.cri,0.5,0.75)) +
                ggplot2::scale_linetype_manual("",labels=comparisons.label,values=linetypes)
            
            if(!length(he$kstar)==0 & !is.na(size)) {
                # label
                label <- paste0("k* = ",format(he$kstar,digits=6))
                eib <-eib +
                    ggplot2::geom_vline(ggplot2::aes(xintercept=kstar),data=data.frame("kstar"=he$kstar),colour="grey50",linetype=2,size=.5) + 
                    ggplot2::annotate("text",label=label,x=he$kstar,y=min(yl),hjust=ifelse((max(he$k)-he$kstar)/max(he$k)>1/6,-.1,1.1),size=size,vjust=1)
            }
            
            if(plot.cri){
                eib <- eib +
                    ggplot2::geom_line(ggplot2::aes(y=low),colour="grey50",show_guide=F)+
                    ggplot2::geom_line(ggplot2::aes(y=upp),colour="grey50",show_guide=F)
            }
        }
        
        if(he$n.comparisons>1&is.null(comparison)==FALSE) {
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            return(eib.plot(he,pos=alt.legend,graph="ggplot2",size=size,comparison=NULL,plot.cri=plot.cri,alpha=alpha,cri.quantile=cri.quantile,...))
        }
        
        eib <- eib + 
            ggplot2::labs(x="Willingness to pay",y="EIB",
                          title=paste0("Expected Incremental Benefit",
                                       ifelse(plot.cri,paste0("\nand ",format((1-alpha)*100,digits=4),"% credible intervals"),"")))
        
        jus <- NULL
        if(isTRUE(alt.legend)) {
            alt.legend="bottom"
            eib <- eib + ggplot2::theme(legend.direction="vertical")
        }
        else{
            if(is.character(alt.legend)) {
                choices <- c("left", "right", "bottom", "top")
                alt.legend <- choices[pmatch(alt.legend,choices)]
                jus="center"
                if(is.na(alt.legend))
                    alt.legend=FALSE
            }
            if(length(alt.legend)>1)
                jus <- alt.legend
            if(length(alt.legend)==1 & !is.character(alt.legend)) {
                alt.legend <- c(1,0)
                jus <- alt.legend
            }
        }
        eib <- eib + 
            ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),text=ggplot2::element_text(size=11),
                           legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),
                           panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),legend.text.align=0,
                           plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3)) +
            opt.theme
        
        return(eib)
    } # !base.graphics
}


###ceac.plot##################################################################################################
## Plots the CEAC
ceac.plot <- function(he,comparison=NULL,pos=c(1,0),graph=c("base","ggplot2")) {
    options(scipen=10)
    
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics) {
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==1)
                temp <- paste0(temp,"top")
            else
                temp <- paste0(temp,"bottom")
            if(alt.legend[1]==0)
                temp <- paste0(temp,"left")
            else
                temp <- paste0(temp,"right")
            alt.legend <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="bottomright"
            else
                alt.legend="bottomleft"
        }
        
        if(he$n.comparisons==1) {
            plot(he$k,he$ceac,t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
                 ylim=c(0,1),main="Cost Effectiveness\nAcceptability Curve")
        }
        if(he$n.comparisons>1&is.null(comparison)) {
            color <- rep(1,he$n.comparisons); lwd <- 1
            if (he$n.comparisons>6) {
                cl <- colors()
                color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	# gray scale
                lwd <- 1.5
            }
            
            plot(he$k,he$ceac[,1],t="l",xlab="Willingness to pay",ylab="Probability of cost effectiveness",
                 ylim=c(0,1),main="Cost Effectiveness\nAcceptability Curve",lty=1,lwd=lwd)
            for (j in 2:he$n.comparisons) {
                points(he$k,he$ceac[,j],t="l",col=color[j],lty=j,lwd=lwd)
            }
            text <- paste(he$interventions[he$ref]," vs ",he$interventions[he$comp])
            legend(alt.legend,text,col=color,cex=.7,bty="n",lty=1:he$n.comparisons)
        }
        if(he$n.comparisons>1&!is.null(comparison)) {
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            ceac.plot(he,pos=alt.legend,graph="base")
        }
    } # base.graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("falling back to base graphics\n")
            ceac.plot(he,pos=alt.legend,graph="base"); return(invisible(NULL))
        }
        
        # no visible binding note
        k = NA_real_
        
        if(he$n.comparisons==1) {
            data.psa <- with(he,data.frame("k"=k,"ceac"=ceac))
            ceac <- ggplot2::ggplot(data.psa, ggplot2::aes(k,ceac)) + ggplot2::geom_line() 
        }
        if(he$n.comparisons>1 & is.null(comparison)==TRUE) {
            data.psa <- with(he,data.frame("k"=c(k),"ceac"=c(ceac),"comparison"=as.factor(sort(rep(1:n.comparisons,length(k))))))
            
            # labels for legend
            comparisons.label <- with(he,paste0(interventions[ref]," vs ",interventions[comp]))
            
            # linetype is the indicator
            linetypes <- rep(c(1,2,3,4,5,6),ceiling(he$n.comparisons/6))[1:he$n.comparisons]
            
            ceac <- ggplot2::ggplot(data.psa,ggplot2::aes(k,ceac,linetype=comparison)) +
                ggplot2::geom_line() +
                ggplot2::scale_linetype_manual("",labels=comparisons.label,values=linetypes)
        }
        
        if(he$n.comparisons>1&is.null(comparison)==FALSE) {
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            return(ceac.plot(he,pos=alt.legend,graph="ggplot2"))
        }
        
        ceac <- ceac + ggplot2::theme_bw() + 
            ggplot2::scale_y_continuous(limits=c(0,1)) +
            ggplot2::labs(title="Cost-Effectiveness Acceptability Curve",x="Willingness to pay",y="Probability of cost-effectiveness") 
        
        jus <- NULL
        if(isTRUE(alt.legend)) {
            alt.legend="bottom"
            ceac <- ceac + ggplot2::theme(legend.direction="vertical")
        }
        else{
            if(is.character(alt.legend)) {
                choices <- c("left", "right", "bottom", "top")
                alt.legend <- choices[pmatch(alt.legend,choices)]
                jus="center"
                if(is.na(alt.legend))
                    alt.legend=FALSE
            }
            if(length(alt.legend)>1)
                jus <- alt.legend
            if(length(alt.legend)==1 & !is.character(alt.legend)) {
                alt.legend <- c(1,0)
                jus <- alt.legend
            }
        }
        
        ceac <- ceac + 
            ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),text=ggplot2::element_text(size=11),
                           legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),
                           panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),legend.text.align=0,
                           plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        return(ceac)
    } # !base.graphics
}

###evi.plot###################################################################################################
## Plots the EVI
evi.plot <- function(he,graph=c("base","ggplot2")) {
    options(scipen=10)
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics){
        
        plot(he$k,he$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
             main="Expected Value of Information")
        if(length(he$kstar)==1) {
            points(rep(he$kstar,3),c(-10000,he$evi[he$k==he$kstar]/2,he$evi[he$k==he$kstar]),t="l",lty=2,col="dark grey")
            points(c(-10000,he$kstar/2,he$kstar),rep(he$evi[he$k==he$kstar],3),t="l",lty=2,col="dark grey")
        }
        if(length(he$kstar)>1) {
            for (i in 1:length(he$kstar)) {
                points(rep(he$kstar[i],3),c(-10000,he$evi[he$k==he$kstar[i]]/2,he$evi[he$k==he$kstar[i]]),
                       t="l",lty=2,col="dark grey")
                points(c(-10000,he$kstar[i]/2,he$kstar[i]),rep(he$evi[he$k==he$kstar[i]],3),t="l",lty=2,col="dark grey")
            }
        }
    } # base.graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("falling back to base graphics\n")
            evi.plot(he,graph="base"); return(invisible(NULL))
        }
        
        ### no visible binding note
        k <- NA_real_
        
        data.psa <- with(he,data.frame("k"=c(k),"evi"=c(evi)))
        
        evi <- ggplot2::ggplot(data.psa, ggplot2::aes(k,evi)) + ggplot2::geom_line() + ggplot2::theme_bw() +
            ggplot2::labs(title="Expected Value of Information",x="Willingness to pay",y="EVPI") +
            ggplot2::theme(plot.title=ggplot2::element_text(face="bold"))
        
        if(length(he$kstar)!=0) {
            kstars=length(he$kstar)
            evi.at.kstar <- numeric(kstars)
            for(i in 1:kstars) {
                evi.at.kstar[i] <- with(he,evi[which.min(abs(k-kstar[i]))])
            }
            
            for(i in 1:kstars) {
                evi <- evi +
                    ggplot2::annotate("segment",x=he$kstar[i],xend=he$kstar[i],y=evi.at.kstar[i],yend=-Inf,linetype=2,colour="grey50") +
                    ggplot2::annotate("segment",x=he$kstar[i],xend=-Inf,y=evi.at.kstar[i],yend=evi.at.kstar[i],linetype=2,colour="grey50")
            }
        }
        
        evi <- evi +
            ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"lines"),
                           legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),
                           legend.key=ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        return(evi)
    }
}

###plot.bcea##################################################################################################
## Plots the main health economics outcomes in just one graph
plot.bcea <- function(x,comparison=NULL,wtp=25000,pos=FALSE,graph=c("base","ggplot2"),...) {
    options(scipen=10)
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics) {
        op <- par(mfrow=c(2,2))
        ceplane.plot(x,comparison=comparison,wtp=wtp,pos=pos,graph="base",...)
        eib.plot(x,comparison=comparison,pos=pos,graph="base",...)
        ceac.plot(x,comparison=comparison,pos=pos,graph="base")
        evi.plot(x,graph="base")
        par(op)
    }
    else{
        
        if(!requireNamespace("ggplot2",quietly=TRUE) & !requireNamespace("grid",quietly=TRUE)){
            message("falling back to base graphics\n")
            plot.bcea(x,comparison=comparison,wtp=wtp,pos=pos,graph="base",...)
            return(invisible(NULL))
        }
        
        ####### multiplot ###### 
        # source: R graphics cookbook
        if(requireNamespace("ggplot2",quietly=TRUE) & requireNamespace("grid",quietly=TRUE)){
            multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
                plots <- c(list(...),plotlist)
                numPlots = length(plots)
                if(is.null(layout)) {
                    layout <- matrix(seq(1,cols*ceiling(numPlots/cols)),
                                     ncol=cols, nrow=ceiling(numPlots/cols))
                }
                if(numPlots==1) {
                    print(plots[[1]])
                } else {
                    grid::grid.newpage()
                    grid::pushViewport(grid::viewport(layout=grid::grid.layout(nrow(layout),ncol(layout))))
                    
                    for(i in 1:numPlots) {
                        matchidx <- as.data.frame(which(layout==i,arr.ind=TRUE))
                        print(plots[[i]],vp=grid::viewport(layout.pos.row=matchidx$row,
                                                           layout.pos.col=matchidx$col))
                    }
                }
            } #### multiplot end ####
            
            theme.multiplot <- 
                ggplot2::theme(text=ggplot2::element_text(size=9),legend.key.size=grid::unit(.5,"lines"),
                               legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),
                               legend.key=ggplot2::element_blank(),plot.title=ggplot2::element_text(lineheight=1,face="bold",size=11.5))
            
            exArgs <- list(...)
            for(obj in exArgs)
                if(ggplot2::is.theme(obj))
                    theme.multiplot <- theme.multiplot + obj
            
            ceplane.pos <- pos
            if(isTRUE(pos==FALSE)){
                ceplane.pos <- c(1,1.025)
            }
            ceplane <- ceplane.plot(x,wtp=wtp,pos=ceplane.pos,comparison=comparison,graph="ggplot2",...) +
                theme.multiplot
            eib <- eib.plot(x,pos=pos,comparison=comparison,graph="ggplot2",...) +
                theme.multiplot
            ceac <- ceac.plot(x,pos=pos,comparison=comparison,graph="ggplot2") +
                theme.multiplot
            evi <- evi.plot(x,graph="ggplot2") +
                theme.multiplot
            # then call multiplot
            multiplot(ceplane,ceac,eib,evi,cols=2)
        } # !base.graphics
    }
}

###contours##################################################################################################

## Contour plots for the cost-effectiveness plane
contour.bcea <- function(x,comparison=1,scale=0.5,nlevels=4,levels=NULL,pos=c(1,0),
                         xlim=NULL,ylim=NULL,graph=c("base","ggplot2"),...) {
    requireNamespace("MASS")
    options(scipen=10)
    # comparison selects which plot should be made
    # by default it is the first possible
    
    # Additional/optional arguments
    exArgs <- list(...)
    if(!exists("xlab",where=exArgs)){xlab <- "Effectiveness differential"} else {xlab <- exArgs$xlab}
    if(!exists("ylab",where=exArgs)){ylab <- "Cost differential"} else {ylab <- exArgs$ylab}
    if(!exists("title",where=exArgs)){title <- paste("Cost effectiveness plane contour plot\n",x$interventions[x$ref]," vs ",x$interventions[x$comp],sep="")} 
    else {title <- exArgs$title}
    
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics){
        
        if(is.null(comparison) | length(comparison) > 1){
            message("The first available comparison will be selected. To plot multiple comparisons together please use the ggplot2 version. Please see ?contour.bcea for additional details.")
            comparison <- 1
        }
        
        if (x$n.comparisons==1) {
            density <- MASS::kde2d(x$delta.e,x$delta.c,n=300,h=c(sd(x$delta.e)/scale,sd(x$delta.c)/scale))
            offset <- 1.0
            
            p.ne <- sum(x$delta.e>0 & x$delta.c>0)/x$n.sim
            p.nw <- sum(x$delta.e<=0 & x$delta.c>0)/x$n.sim
            p.sw <- sum(x$delta.e<=0 & x$delta.c<=0)/x$n.sim
            p.se <- sum(x$delta.e>0 & x$delta.c<=0)/x$n.sim
            
            m.c <- range(x$delta.c)[1]; M.c <- range(x$delta.c)[2]
            m.e <- range(x$delta.e)[1]; M.e <- range(x$delta.e)[2]
            
            # Changes the range so that the plot always shows the x and y axes
            ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
            ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
            ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
            ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)
            
            # If the user has specified the range of the graph, use those values
            if(!is.null(xlim)) {m.e <- xlim[1]; M.e <- xlim[2]}
            if(!is.null(ylim)) {m.c <- ylim[1]; M.c <- ylim[2]}
            
            plot(x$delta.e,x$delta.c,pch=20,cex=.3,col="dark grey",xlab=xlab,
                 ylab=ylab,main=title,
                 xlim=c(m.e,M.e),ylim=c(m.c,M.c))
            abline(h=0,col="dark grey")
            abline(v=0,col="dark grey")
            if(any(is.na(density$z))==FALSE) {
                if (is.null(levels)==FALSE){
                    # Normalise the density and use levels in the contour
                    density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
                    contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
                }
                if (is.null(levels)==TRUE) {
                    contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
                }
            }
            t1 <- paste("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3),sep="")
            text(offset*M.e,offset*M.c,parse(text=t1),cex=.8,pos=2)
            t2 <- paste("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3),sep="")
            text(offset*m.e,offset*M.c,parse(text=t2),cex=.8,pos=4)
            t3 <- paste("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3),sep="")
            text(offset*m.e,offset*m.c,parse(text=t3),cex=.8,pos=4)
            t4 <- paste("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3),sep="")
            text(offset*M.e,offset*m.c,parse(text=t4),cex=.8,pos=2)
        }
        
        if(x$n.comparisons>1) {
            if(!exists("title",where=exArgs)){title <- paste("Cost effectiveness plane contour plot \n",x$interventions[x$ref]," vs ",
                                                             x$interventions[x$comp[comparison]],sep="")} 
            else {title <- exArgs$title}
            
            density <- MASS::kde2d(x$delta.e[,comparison],x$delta.c[,comparison],n=300,
                                   h=c(sd(x$delta.e[,comparison])/scale,sd(x$delta.c[,comparison])/scale))
            offset <- 1.0
            
            p.ne <- sum(x$delta.e[,comparison]>0 & x$delta.c[,comparison]>0)/x$n.sim
            p.nw <- sum(x$delta.e[,comparison]<=0 & x$delta.c[,comparison]>0)/x$n.sim
            p.sw <- sum(x$delta.e[,comparison]<=0 & x$delta.c[,comparison]<=0)/x$n.sim
            p.se <- sum(x$delta.e[,comparison]>0 & x$delta.c[,comparison]<=0)/x$n.sim
            
            m.c <- range(x$delta.c[,comparison])[1]; M.c <- range(x$delta.c[,comparison])[2]
            m.e <- range(x$delta.e[,comparison])[1]; M.e <- range(x$delta.e[,comparison])[2]
            
            # Changes the range so that the plot always shows the x and y axes
            ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
            ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
            ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
            ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)
            
            # If the user has specified the range of the graph, use those values
            if(!is.null(xlim)) {m.e <- xlim[1]; M.e <- xlim[2]}
            if(!is.null(ylim)) {m.c <- ylim[1]; M.c <- ylim[2]}
            
            plot(x$delta.e[,comparison],x$delta.c[,comparison],pch=20,cex=.3,col="dark grey",
                 xlab=xlab,ylab=ylab,
                 main=title,xlim=c(m.e,M.e),ylim=c(m.c,M.c))	
            abline(h=0,col="dark grey")
            abline(v=0,col="dark grey")
            if(any(is.na(density$z))==FALSE) {
                contour(density$x,density$y,density$z,add=TRUE,drawlabels=TRUE)
                if (is.null(levels)==FALSE){
                    # Normalise the density and use levels in the contour
                    density$z <- (density$z-min(density$z))/(max(density$z)-min(density$z))
                    contour(density$x,density$y,density$z,add=TRUE,levels=levels,drawlabels=TRUE)
                }
                if (is.null(levels)==TRUE) {
                    contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE)
                }
            }
            t1 <- paste("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3),sep="")
            text(offset*M.e,offset*M.c,parse(text=t1),cex=.8,pos=2)
            t2 <- paste("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3),sep="")
            text(offset*m.e,offset*M.c,parse(text=t2),cex=.8,pos=4)
            t3 <- paste("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3),sep="")
            text(offset*m.e,offset*m.c,parse(text=t3),cex=.8,pos=4)
            t4 <- paste("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3),sep="")
            text(offset*M.e,offset*m.c,parse(text=t4),cex=.8,pos=2)
        }
    } # if base.graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("falling back to base graphics\n")
            contour.bcea(x,comparison=comparison,scale=scale,nlevels=nlevels,pos=alt.legend,levels=levels,graph="base",...)
            return(invisible(NULL))
        }
        
        if(!is.null(levels))
            message("option level will be ignored using ggplot2 graphics")
        
        # no visible binding note
        delta.e <- delta.c <- e <- z <- y <- hjust <- label <- NULL
        
        if(!is.null(nlevels)){
            nlevels <- round(nlevels)
            if(nlevels<0)
                nlevels <- 10
            if(nlevels==0)
                nlevels <- 1
        }
        
        if(x$n.comparisons==1) {
            kd <- with(x,data.frame("e"=delta.e,"c"=delta.c))
            
            # for scale_x_continuous(oob=)
            do.nothing=function(x,limits) return(x)
            # plot limits
            range.e <- range(kd$e)
            range.c <- range(kd$c)
            range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
            range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
            
            # labels
            p.ne <- sum(x$delta.e > 0 & x$delta.c > 0)/x$n.sim
            p.ne <- paste0("Pr(Delta[e]>0, Delta[c]>0)==",format(p.ne,digits=4,nsmall=3))
            p.nw <- sum(x$delta.e <= 0 & x$delta.c > 0)/x$n.sim
            p.nw <- paste0("Pr(Delta[e]<=0, Delta[c]>0)==",format(p.nw,digits=4,nsmall=3))
            p.sw <- sum(x$delta.e <= 0 & x$delta.c <= 0)/x$n.sim
            p.sw <- paste0("Pr(Delta[e]<=0, Delta[c]<=0)==",format(p.sw,digits=4,nsmall=3))
            p.se <- sum(x$delta.e > 0 & x$delta.c <= 0)/x$n.sim
            p.se <- paste0("Pr(Delta[e]>0, Delta[c]<=0)==",format(p.se,digits=4,nsmall=3))
            
            # labels dataframe
            labels.df <- data.frame(
                "x"=c(range.e[2],range.e[1],range.e[1],range.e[2]),
                "y"=c(rep(range.c[2],2),rep(range.c[1],2)),
                "label"=c(p.ne,p.nw,p.sw,p.se),
                "hjust"=as.factor(c(1,0,0,1)))
            
            # actual plot
            points.colour="grey"
            if(nlevels==1)
                points.colour="black"
            ceplane <- ggplot2::ggplot(kd, ggplot2::aes(e,c)) + ggplot2::geom_hline(ggplot2::aes(yintercept=0),colour="grey") + 
                ggplot2::geom_vline(ggplot2::aes(xintercept=0),colour="grey") + ggplot2::theme_bw() + 
                ggplot2::geom_point(size=1,color=points.colour) + ggplot2::scale_x_continuous(limits=range.e,oob=do.nothing) + 
                ggplot2::scale_y_continuous(limits=range.c,oob=do.nothing)
            if(!is.null(scale)&requireNamespace("MASS",quietly=TRUE)){
                density <- with(x,MASS::kde2d(delta.e,delta.c,n=300,h=c(sd(delta.e)/scale,sd(delta.c)/scale)))
                density <- data.frame(expand.grid("e"=density$x,"c"=density$y),"z"=as.vector(density$z))
                ceplane <- ceplane + ggplot2::geom_contour(ggplot2::aes(z=z),data=density,colour="black",bins=nlevels)
            }
            else{
                ceplane <- ceplane + ggplot2::stat_density2d(color="black")
            }
            
            
            ceplane <- ceplane + 
                ggplot2::geom_text(data=labels.df,ggplot2::aes(x=x,y=y,hjust=hjust,label=label),parse=TRUE,size=ggplot2::rel(3.5))
        }
        if(x$n.comparisons>1&is.null(comparison)==TRUE) {
            # creates dataframe for plotting
            kd <- with(x,data.frame(
                "delta.e"=c(delta.e),"delta.c"=c(delta.c),
                "comparison"=as.factor(sort(rep(1:n.comparisons,dim(delta.e)[1])))))
            
            # vector of values for color, take out white, get integer values
            colors.label <- paste0("gray",round(seq(0,100,length.out=(x$n.comparisons+1))[-(x$n.comparisons+1)]))
            comparisons.label <- paste0(x$interventions[x$ref]," vs ",x$interventions[x$comp])
            do.nothing=function(x,limits) return(x)
            # plot limits
            range.e <- range(kd$delta.e)
            range.c <- range(kd$delta.c)
            range.e[1] <- ifelse(range.e[1]<0,range.e[1],-range.e[1])
            range.c[1] <- ifelse(range.c[1]<0,range.c[1],-range.c[1])
            
            ceplane <-
                ggplot2::ggplot(kd,ggplot2::aes(x=delta.e,y=delta.c,col=comparison)) +
                ggplot2::geom_hline(yintercept=0,colour="grey") + ggplot2::geom_vline(xintercept=0,colour="grey") + ggplot2::theme_bw() +
                ggplot2::geom_point(size=1) +
                ggplot2::scale_color_manual(label=comparisons.label,values=colors.label,na.value="black") +
                ggplot2::scale_x_continuous(limits=range.e,oob=do.nothing) +
                ggplot2::scale_y_continuous(limits=range.c,oob=do.nothing)
            
            if(!is.null(scale)&requireNamespace("MASS",quietly=TRUE)) {
                densitydf <- data.frame()
                for(i in 1:x$n.comparison) {
                    temp <- with(x,MASS::kde2d(delta.e[,i],delta.c[,i],n=300,h=c(sd(delta.e[,i])/scale,sd(delta.c[,i])/scale)))
                    temp <- data.frame(expand.grid("e"=temp$x,"c"=temp$y),"z"=as.vector(temp$z))
                    densitydf <- rbind(densitydf,cbind(temp,rep(i,dim(temp)[[1]])))
                }
                names(densitydf) <- c("delta.e","delta.c","z","comparison")
                densitydf$comparison <- as.factor(densitydf$comparison)
                ceplane <- ceplane + ggplot2::geom_contour(ggplot2::aes(z=z,colour=comparison),data=densitydf,bins=nlevels) +
                    ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(linetype=0)))
            }
            else{
                ceplane <- ceplane + ggplot2::stat_density2d() + 
                    ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(linetype=0)))
            }
        }
        if(x$n.comparisons>1&is.null(comparison)==FALSE) {
            # adjusts bcea object for the correct number of dimensions and comparators
            x$comp <- x$comp[comparison]
            x$delta.e <- x$delta.e[,comparison]
            x$delta.c <- x$delta.c[,comparison]
            x$n.comparators=length(comparison)+1
            x$n.comparisons=length(comparison)
            x$interventions=x$interventions[sort(c(x$ref,x$comp))]
            x$ICER=x$ICER[comparison]
            x$ib=x$ib[,,comparison]
            x$eib=x$eib[,comparison]
            x$U=x$U[,,sort(c(x$ref,comparison+1))]
            x$ceac=x$ceac[,comparison]
            x$ref=rank(c(x$ref,x$comp))[1]
            x$comp=rank(c(x$ref,x$comp))[-1]
            x$mod <- TRUE #
            
            return(contour.bcea(x,scale=scale,pos=alt.legend,nlevels=nlevels,graph="ggplot2",comparison=NULL))
        }
        
        if(!exists("title",where=exArgs)) {
            labs.title <- "Cost-Effectiveness Plane"
            labs.title <- paste0(labs.title,
                                 ifelse(x$n.comparisons==1,
                                        paste0("\n",x$interventions[x$ref]," vs ",x$interventions[-x$ref]),
                                        paste0(
                                            ifelse(isTRUE(x$mod),
                                                   paste0("\n",x$interventions[x$ref]," vs ",
                                                          paste0(x$interventions[x$comp],collapse=", ")),
                                                   ""))))
        } else {labs.title <- exArgs$title}
        
        ceplane <- ceplane + ggplot2::labs(title=labs.title,x=xlab,y=ylab)
        
        jus <- NULL
        if(isTRUE(alt.legend)) {
            alt.legend="bottom"
            ceplane <- ceplane + ggplot2::theme(legend.direction="vertical")
        }
        else{
            if(is.character(alt.legend)) {
                choices <- c("left", "right", "bottom", "top")
                alt.legend <- choices[pmatch(alt.legend,choices)]
                jus="center"
                if(is.na(alt.legend))
                    alt.legend=FALSE
            }
            if(length(alt.legend)>1)
                jus <- alt.legend
            if(length(alt.legend)==1 & !is.character(alt.legend)) {
                alt.legend <- c(1,0)
                jus <- alt.legend
            }
        }
        
        ceplane <- ceplane + 
            ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),text=ggplot2::element_text(size=11),
                           legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),
                           panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),legend.text.align=0,
                           plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        return(ceplane)
    } # !base.graphics
}

#####
contour2 <- function(he,wtp=25000,xlim=NULL,ylim=NULL,comparison=NULL,graph=c("base","ggplot2"),...) {
    # Forces R to avoid scientific format for graphs labels
    options(scipen=10)
    
    # Additional/optional arguments
    exArgs <- list(...)
    if(!exists("xlab",where=exArgs)){xlab <- "Effectiveness differential"} else {xlab <- exArgs$xlab}
    if(!exists("ylab",where=exArgs)){ylab <- "Cost differential"} else {ylab <- exArgs$ylab}
    if(!exists("title",where=exArgs)){
        title <- paste("Cost effectiveness plane \n",he$interventions[he$ref]," vs ",he$interventions[he$comp],sep="")} 
    else {title <- exArgs$title
    }
    
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics) {
        # Encodes characters so that the graph can be saved as ps or pdf
        ps.options(encoding="CP1250")
        pdf.options(encoding="CP1250")
        
        # Selects the first comparison by default if not selected
        if(is.null(comparison)){
            message("The first available comparison will be selected. To plot multiple comparisons together please use the ggplot2 version. Please see ?contour2 for additional details.")
            comparison <- 1
        }
        
        if(he$n.comparisons>1) {
            if(!exists("title",where=exArgs)){title <- paste("Cost effectiveness plane contour plot \n",he$interventions[he$ref]," vs ",
                                                             he$interventions[he$comp[comparison]],sep="")} 
            else {title <- exArgs$title}
            
            
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$comp <- he$comp[comparison]
            he$ICER <- he$ICER[comparison]
        }
        
        m.e <- range(he$delta.e)[1]
        M.e <- range(he$delta.e)[2]
        m.c <- range(he$delta.c)[1]
        M.c <- range(he$delta.c)[2]
        step <- (M.e-m.e)/10
        
        m.e <- ifelse(m.e<0,m.e,-m.e)
        m.c <- ifelse(m.c<0,m.c,-m.c)
        
        x.pt <- .95*m.e
        y.pt <- ifelse(x.pt*wtp<m.c,m.c,x.pt*wtp)
        xx <- seq(100*m.c/wtp,100*M.c/wtp,step)
        yy <- xx*wtp          
        xx[1] <- ifelse(min(xx)<m.e,xx[1],2*m.e)
        yy[1] <- ifelse(min(yy)<m.c,yy[1],2*m.c)
        xx[length(xx)] <- ifelse(xx[length(xx)]<M.e,1.5*M.e,xx[length(xx)])
        
        # If the user has specified x- and/or y-limits then use those
        if(!is.null(xlim)) {m.e <- xlim[1]; M.e <- xlim[2]}
        if(!is.null(ylim)) {m.c <- ylim[1]; M.c <- ylim[2]}
        
        plot(xx,yy,col="white",xlim=c(m.e,M.e),ylim=c(m.c,M.c), 
             xlab=xlab,ylab=ylab,
             main=title,axes=F)
        polygon(c(min(xx),seq(min(xx),max(xx),step),max(xx)),
                c(min(yy),wtp*seq(min(xx),max(xx),step),min(yy)),
                col="grey95",border="black")
        #	polygon(c(xx,xx),c(yy,rev(yy)),col="grey95",border="black")
        axis(1); axis(2); box()
        points(he$delta.e,he$delta.c,pch=20,cex=.35,col="grey55")
        abline(h=0,col="dark grey")
        abline(v=0,col="dark grey")
        text(M.e,M.c,paste("\U2022"," ICER=",format(he$ICER,digits=6,nsmall=2),sep=""),cex=.95,pos=2,col="red")
        points(mean(he$delta.e),mean(he$delta.c),pch=20,col="red",cex=1)
        t1 <- paste("k==",format(wtp,digits=3,nsmall=2,scientific=F),sep="")
        text(x.pt,y.pt,parse(text=t1),cex=.8,pos=4)
        
        # And then plots the contour
        requireNamespace("MASS")
        offset <- 1.0
        nlevels <- 4
        scale <- 0.5
        
        density <- MASS::kde2d(he$delta.e,he$delta.c,n=300,h=c(sd(he$delta.e)/scale,sd(he$delta.c)/scale))
        
        m.c <- range(he$delta.c)[1]; M.c <- range(he$delta.c)[2]
        m.e <- range(he$delta.e)[1]; M.e <- range(he$delta.e)[2]
        
        # Changes the range so that the plot always shows the x and y axes
        ch1 <- ifelse(m.e>0,m.e<--m.e,m.e<-m.e)
        ch2 <- ifelse(M.e<0,M.e<--M.e,M.e<-M.e)
        ch3 <- ifelse(m.c>0,m.c<--m.c,m.c<-m.c)
        ch4 <- ifelse(M.c<0,M.c<--M.c,M.c<-M.c)
        
        par(new=TRUE)
        contour(density$x,density$y,density$z,add=TRUE,nlevels=nlevels,drawlabels=FALSE,lwd=1.5)
        return(invisible(NULL))
    } # end if base.graphics
    else{
        
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("falling back to base graphics\n")
            contour2(he,comparison=comparison,xlim=xlim,ylim=ylim,wtp=wtp,graph="base"); return(invisible(NULL))
        }
        scale=0.5
        nlevels=5
        
        requireNamespace("MASS")
        
        ### no visible binding note
        z <- e <- NA_real_
        
        if(he$n.comparisons==1) {
            density <- with(he,MASS::kde2d(delta.e,delta.c,n=300,h=c(sd(delta.e)/scale,sd(delta.c)/scale)))
            density <- data.frame(expand.grid("e"=density$x,"c"=density$y),"z"=as.vector(density$z))
            contour <- ceplane.plot(he,wtp=wtp,graph="ggplot2",...) +
                ggplot2::geom_contour(ggplot2::aes(z=z,x=e,y=c),data=density,colour="black",bins=nlevels)
        }
        if(he$n.comparisons>1&is.null(comparison)) {
            densitydf <- data.frame()
            for(i in 1:he$n.comparisons) {
                density <- with(he,MASS::kde2d(delta.e[,i],delta.c[,i],n=300,h=c(sd(delta.e[,i])/scale,sd(delta.c[,i])/scale)))
                densitydf <- rbind(densitydf,cbind(expand.grid(density$x,density$y),as.vector(density$z)))
            }
            names(densitydf) <- c("e","c","z")
            densitydf <- cbind(densitydf,"comparison"=as.factor(sort(rep(1:he$n.comparisons,dim(densitydf)[1]/he$n.comparisons))))
            contour <- ceplane.plot(he,wtp=wtp,graph="ggplot2",...) +
                ggplot2::geom_contour(data=densitydf,ggplot2::aes(x=e,y=c,z=z,colour=comparison),bins=nlevels,linetype=1)
        }
        if(he$n.comparisons>1&!is.null(comparison)) {
            # adjusts bcea object for the correct number of dimensions and comparators
            he$comp <- he$comp[comparison]
            he$delta.e <- he$delta.e[,comparison]
            he$delta.c <- he$delta.c[,comparison]
            he$n.comparators=length(comparison)+1
            he$n.comparisons=length(comparison)
            he$interventions=he$interventions[sort(c(he$ref,he$comp))]
            he$ICER=he$ICER[comparison]
            he$ib=he$ib[,,comparison]
            he$eib=he$eib[,comparison]
            he$U=he$U[,,sort(c(he$ref,comparison+1))]
            he$ceac=he$ceac[,comparison]
            he$ref=rank(c(he$ref,he$comp))[1]
            he$comp=rank(c(he$ref,he$comp))[-1]
            he$mod <- TRUE #
            
            return(contour2(he,wtp=wtp,xlim=xlim,ylim=ylim,comparison=NULL,graph="ggplot2",...))
        }
        
        contour <- contour + ggplot2::coord_cartesian(xlim=xlim,ylim=ylim)
        return(contour)
    } # end if !base.graphics
    
}



###CEriskav###################################################################################################
CEriskav <- function(he,r=NULL,comparison=1) UseMethod("CEriskav")

CEriskav.default <- function(he,r=NULL,comparison=1) {
    ### COMPARISON IS USED TO SELECT THE COMPARISON FOR WHICH THE ANALYSIS IS CARRIED OUT!!!
    # Reference: Baio G, Dawid AP (2011).
    # Default vector of risk aversion parameters
    if(is.null(r)==TRUE){
        r <- c(0.000000000001,0.0000025,.000005)
    }
    
    # Computes expected utilities & EVPI for the risk aversion cases
    K <- length(he$k)
    R <- length(r)
    Ur <- array(NA,c(dim(he$U),R))
    Urstar <- array(NA,c(dim(he$Ustar),R))
    for (i in 1:K) {
        for (l in 1:R) {
            for (j in 1:he$n.comparators) {
                Ur[,i,j,l] <- (1/r[l])*(1-exp(-r[l]*he$U[,i,j]))
            }
            Urstar[,i,l] <- apply(Ur[,i,,l],1,max)
        }
    }
    
    if (he$n.comparisons==1){
        IBr <- Ur[,,he$ref,] - Ur[,,he$comp,]
    }
    if (he$n.comparisons>1){
        IBr <- Ur[,,he$ref,] - Ur[,,he$comp[comparison],]
    }
    
    eibr <- apply(IBr,c(2,3),mean)
    vir <- array(NA,c(he$n.sim,K,R))
    for (i in 1:K) {
        for (l in 1:R) {
            vir[,i,l] <- Urstar[,i,l] - max(apply(Ur[,i,,l],2,mean))
        }
    }
    evir <- apply(vir,c(2,3),mean)
    
    ## Outputs of the function
    cr <- list(
        Ur=Ur,Urstar=Urstar,IBr=IBr,eibr=eibr,vir=vir,evir=evir,R=R,r=r,k=he$k
    )
    class(cr) <- "CEriskav"
    cr
}

###plot.CEriskav##################################################################################################
# Plots the EIB for the risk aversion case
plot.CEriskav <- function(x,pos=c(0,1),graph=c("base","ggplot2"),...) {
    options(scipen=10)
    
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    howplot <- NULL
    exArgs <- list(...)
    if(length(exArgs)>=1)
        if(exists("plot",where=exArgs))
            howplot <- exArgs$plot
    
    if(base.graphics) {
        
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==0)
                temp <- paste0(temp,"bottom")
            else
                temp <- paste0(temp,"top")
            if(alt.legend[1]==1)
                temp <- paste0(temp,"right")
            else
                temp <- paste0(temp,"left")
            alt.legend <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="topright"
            else
                alt.legend="topleft"
        }
        
        plot(x$k,x$eibr[,1],t="l",xlab="Willingness to pay", ylab=" ", main="EIB as a function of the risk aversion parameter",ylim=range(x$eibr))
        linetype=seq(1,x$R)
        for (l in 2:x$R) {
            points(x$k,x$eibr[,l],t="l",lty=linetype[l])
        }
        text <- paste("r = ",x$r,sep=""); 
        # If the first value for r is small enough, consider it close to 0 and print the label accordingly
        if (x$r[1]<1e-8) {
            text[1] <- expression(r%->%0)
        }
        legend(alt.legend,text,lty=seq(1:x$R),cex=.9,box.lty=0)
        abline(h=0,col="grey")
        
        # Plots the EVPI for the risk aversion case
        if(is.null(howplot)) {
            if(!isTRUE(Sys.getenv("RSTUDIO")==1))
                dev.new()
        }
        else{
            opt <- c("x11","ask","dev.new")
            howplot <- ifelse(is.na(pmatch(howplot,opt)),"dev.new",opt[pmatch(howplot,opt)])
            if(howplot=="x11")
                dev.new()
            if(howplot=="dev.new")
                dev.new()
            if(howplot=="ask")
                devAskNewPage(ask=TRUE)
        }
        
        plot(x$k,x$evir[,1],t="l",ylim=range(x$evir),xlab="Willingness to pay",ylab=" ",main="EVPI as a function of the risk aversion parameter")
        for (l in 2:x$R) {
            points(x$k,x$evir[,l],t="l",lty=linetype[l])
        }
        legend(alt.legend,text,lty=seq(1:x$R),cex=.9,box.lty=0)
        abline(h=0,col="grey")
        
        if(!is.null(howplot))
            if(howplot=="ask")
                devAskNewPage(ask=FALSE)
        
    } # base.graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            message("falling back to base graphics\n")
            plot.CEriskav(x,graph="base",pos=pos,...)
            return(invisible(NULL))
        }
        # no visible bindings note
        k <- r <- NA_real_
        
        linetypes <- rep(c(1,2,3,4,5,6),ceiling(x$R/6))[1:x$R]
        df <- data.frame(cbind(rep(x$k,x$R),c(x$eibr),c(x$evir)),as.factor(sort(rep(1:x$R,length(x$k)))))
        names(df) <- c("k","eibr","evir","r")
        
        # labels
        text <- paste0("r = ",x$r)
        # if the first value for r is small enough, consider it close to 0 and print the label accordingly
        if(x$r[1]<1e-8) {
            text[1] <- expression(r%->%0)
        }
        
        eibr <- ggplot2::ggplot(df,ggplot2::aes(x=k,y=eibr,linetype=r))+
            ggplot2::geom_hline(yintercept=0,linetype=1,colour="grey50")+ggplot2::geom_line()+
            ggplot2::scale_linetype_manual("",labels=text,values=linetypes)+ggplot2::theme_bw() +
            ggplot2::labs(title="EIB as a function of the risk aversion parameter",x="Willingness to pay",y="EIB") +
            ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"line"),
                           legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank())
        
        ### evir ###
        evir <- ggplot2::ggplot(df,ggplot2::aes(x=k,y=evir,linetype=r))+ggplot2::geom_hline(yintercept=0,linetype=1,colour="grey50")+
            ggplot2::geom_line()+ggplot2::scale_linetype_manual("",labels=text,values=linetypes)+ggplot2::theme_bw() +
            ggplot2::labs(title="EVPI as a function of the risk aversion parameter",x="Willingness to pay",y="EVPI") +
            ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"line"),
                           legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank())
        jus <- NULL
        if(isTRUE(alt.legend)) {
            alt.legend="bottom"
            eibr <- eibr + ggplot2::theme(legend.direction="vertical")
            evir <- evir + ggplot2::theme(legend.direction="vertical")
        }
        else{
            if(is.character(alt.legend)) {
                choices <- c("left", "right", "bottom", "top")
                alt.legend <- choices[pmatch(alt.legend,choices)]
                jus="center"
                if(is.na(alt.legend))
                    alt.legend=FALSE
            }
            if(length(alt.legend)>1)
                jus <- alt.legend
            if(length(alt.legend)==1 & !is.character(alt.legend)) {
                alt.legend <- c(0,1)
                jus <- alt.legend
            }
        }
        
        eibr <- eibr + 
            ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"),legend.text.align=0,
                           plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        
        evir <- evir + 
            ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"),
                           legend.text.align=0,plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        plot(eibr)
        
        if(is.null(howplot)) {
            if(!isTRUE(Sys.getenv("RSTUDIO")==1))
                dev.new()
        }
        else{
            opt <- c("x11","ask","dev.new")
            howplot <- ifelse(is.na(pmatch(howplot,opt)),"dev.new",opt[pmatch(howplot,opt)])
            if(howplot=="x11")
                dev.new()
            if(howplot=="dev.new")
                dev.new()
            if(howplot=="ask")
                devAskNewPage(ask=TRUE)
        }
        
        plot(evir)
        
        if(!is.null(howplot))
            if(howplot=="ask")
                devAskNewPage(ask=FALSE)
        
        return(invisible(list("eib"=eibr,"evi"=evir)))
    }
    
}


###mixedAn####################################################################################################
mixedAn <- function(he,mkt.shares=NULL,plot=FALSE) UseMethod("mixedAn")

mixedAn.default <- function(he,mkt.shares=NULL,plot=FALSE) {
    # mkt.shares is a vector of market shares for each comparators
    # if no value is provided, then assumes uniform distribution
    # dev is the device to which the graph should be printed
    # default is x11() --- on screen. Possibilities are jpeg and postscript
    # Reference: Baio G, Russo P (2009).
    
    Ubar <- OL.star <- evi.star <- NULL
    if(is.null(mkt.shares)==TRUE){
        mkt.shares <- rep(1,he$n.comparators)/he$n.comparators
    }
    temp <- array(NA,c(he$n.sim,length(he$k),he$n.comparators))
    for (j in 1:he$n.comparators) {
        temp[,,j] <- mkt.shares[j]*he$U[,,j]
    }
    Ubar <- apply(temp,c(1,2),sum)
    OL.star <- he$Ustar - Ubar
    evi.star <- apply(OL.star,2,mean)
    
    ## Outputs of the function
    ma <- list(
        Ubar=Ubar,OL.star=OL.star,evi.star=evi.star,k=he$k,Kmax=he$Kmax,step=he$step,
        ref=he$ref,comp=he$comp,mkt.shares=mkt.shares,n.comparisons=he$n.comparisons,
        interventions=he$interventions,evi=he$evi
    )
    class(ma) <- "mixedAn"
    if(plot) {
        plot.mixedAn(ma)
    }
    return(ma)
}

###summary.mixedAn############################################################################################
summary.mixedAn <- function(object,wtp=25000,...) {
    if(max(object$k)<wtp) {
        wtp <- max(object$k)
    }
    if(length(which(object$k==wtp))==0) {
        stop(paste("The willingness to pay parameter is defined in the interval [0-",object$Kmax,"], 
                     with increments of ",object$step,"\n",sep=""))
    }
    
    n.digits <- 2
    n.small <- 2
    cat("\n")
    cat(paste("Analysis of mixed strategy for willingness to pay parameter k = ",
              wtp,"\n",sep=""))
    cat("\n")
    cat(paste("Reference intervention: ",object$interventions[object$ref]," (",
              format(100*object$mkt.shares[object$ref],digits=n.digits,nsmall=n.small),"% market share)\n",sep=""))
    if(object$n.comparisons==1) {
        text.temp <- paste("Comparator intervention: ",object$interventions[object$comp]," (",
                           format(100*object$mkt.shares[object$comp],digits=n.digits,nsmall=n.small),"% market share)\n",sep="")
        cat(text.temp)
    }
    
    if(object$n.comparisons>1) {
        text.temp <- paste("Comparator intervention(s): ",object$interventions[object$comp[1]]," (",
                           format(100*object$mkt.shares[object$comp[1]],digits=n.digits,nsmall=n.small),"% market share)\n",sep="")
        cat(text.temp)
        for (i in 2:object$n.comparisons) {
            cat(paste("                          : ",object$interventions[object$comp[i]]," (", 
                      format(100*object$mkt.shares[object$comp[i]],digits=n.digits,nsmall=n.small),"% market share)\n",sep=""))
        }
    }
    cat("\n")
    cat(paste("Loss in the expected value of information = ",
              format(object$evi.star[object$k==wtp]-object$evi[object$k==wtp],digits=n.digits,nsmall=n.small),"\n",sep=""))
    cat("\n")
}

###plot.mixedAn###############################################################################################
plot.mixedAn <- function(x,y.limits=NULL,pos=c(0,1),graph=c("base","ggplot2"),...) {
    ## Plot the EVPI and the mixed strategy
    options(scipen=10)
    
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(is.null(y.limits)){
        y.limits=range(x$evi,x$evi.star)
    }
    
    if(base.graphics) {
        
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==0)
                temp <- paste0(temp,"bottom")
            else
                temp <- paste0(temp,"top")
            if(alt.legend[1]==1)
                temp <- paste0(temp,"right")
            else
                temp <- paste0(temp,"left")
            alt.legend <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="topright"
            else
                alt.legend="topleft"
        }
        
        plot(x$k,x$evi,t="l",xlab="Willingness to pay",ylab="EVPI",
             main="Expected Value of Information",ylim=y.limits)
        polygon(c(x$k,rev(x$k)),c(x$evi.star,rev(x$evi)),density=20,col="grey")
        points(x$k,x$evi.star,t="l",col="red")
        points(x$k,x$evi,t="l",col="black")
        txt <- c("Optimal strategy","Mixed strategy:",
                 paste("   ",x$interventions,"=",format(100*x$mkt.shares,digits=3,nsmall=2),"%",sep=""))
        cols <- c("black","red",rep("white",length(x$interventions)))
        legend(alt.legend,txt,col=cols,cex=.6,bty="n",lty=1)
    } # base.graphics
    else {
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            message("Falling back to base graphics\n")
            plot.mixedAn(x,y.limits=y.limits,pos=pos,graph="base")
            return(invisible(NULL))
        } 
        
        if(isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            ### no visible binding note
            k <- evi.star <- NA_real_
            
            # legend
            txt <- c("Optimal strategy",paste0("Mixed strategy:", paste0("\n   ",x$interventions,"=",format(100*x$mkt.shares,digits=3,nsmall=2),"%",collapse="")))
            colours=c("black","red")
            
            df <- data.frame("k"=x$k,"evi"=x$evi,"evi.star"=x$evi.star)
            evi <- ggplot2::ggplot(df,ggplot2::aes(x=k)) + ggplot2::theme_bw() +
                ggplot2::geom_ribbon(ggplot2::aes(x=k,ymin=evi,ymax=evi.star),colour="lightgrey",alpha=.2) +
                ggplot2::geom_line(ggplot2::aes(x=k,y=evi,colour=as.factor(1))) +
                ggplot2::geom_line(ggplot2::aes(x=k,y=evi.star,colour=as.factor(2))) +
                ggplot2::coord_cartesian(ylim=y.limits,xlim=c(0,max(df$k))) +
                ggplot2::scale_colour_manual("",labels=txt,values=colours) +
                ggplot2::labs(title="Expected Value of Information",x="Willingness to pay",y="EVPI") +
                ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"lines"),
                               legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),
                               legend.key=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"))
            jus <- NULL
            if(isTRUE(alt.legend)) {
                alt.legend="bottom"
                evi <- evi + ggplot2::theme(legend.direction="vertical")
            }
            else{
                if(is.character(alt.legend)) {
                    choices <- c("left", "right", "bottom", "top")
                    alt.legend <- choices[pmatch(alt.legend,choices)]
                    jus="center"
                    if(is.na(alt.legend))
                        alt.legend=FALSE
                }
                if(length(alt.legend)>1)
                    jus <- alt.legend
                if(length(alt.legend)==1 & !is.character(alt.legend)) {
                    alt.legend <- c(0,1)
                    jus <- alt.legend
                }
            }
            
            evi <- evi + 
                ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                               legend.background=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"),
                               legend.text.align=0,plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
            return(evi)
        }
    }
}

#####multi.ce##################################################################################################
multi.ce <- function(he){
    # Cost-effectiveness analysis for multiple comparison 
    # Identifies the probability that each comparator is the most cost-effective as well as the
    # cost-effectiveness acceptability frontier
    cl <- colors()
    # choose colors on gray scale
    color <- cl[floor(seq(262,340,length.out=he$n.comparators))]	
    
    rank <- most.ce <- array(NA,c(he$n.sim,length(he$k),he$n.comparators))
    for (t in 1:he$n.comparators) {
        for (j in 1:length(he$k)) {
            rank[,j,t] <- apply(he$U[,j,]<=he$U[,j,t],1,sum)
            most.ce[,j,t] <- rank[,j,t]==he$n.comparators
        }
    }
    m.ce <- apply(most.ce,c(2,3),mean)		# Probability most cost-effective
    ceaf <- apply(m.ce,1,max)			# Cost-effectiveness acceptability frontier
    
    # Output of the function
    list(
        m.ce=m.ce,ceaf=ceaf,n.comparators=he$n.comparators,k=he$k,interventions=he$interventions
    )
}

############mce.plot###################################
mce.plot <- function(mce,pos=c(1,0.5),graph=c("base","ggplot2")){
    alt.legend <- pos
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    
    if(base.graphics) {
        
        if(is.numeric(alt.legend)&length(alt.legend)==2){
            temp <- ""
            if(alt.legend[2]==0)
                temp <- paste0(temp,"bottom")
            else if(alt.legend[2]!=0.5)
                temp <- paste0(temp,"top")
            if(alt.legend[1]==1)
                temp <- paste0(temp,"right")
            else
                temp <- paste0(temp,"left")
            alt.legend <- temp
            if(length(grep("^((bottom|top)(left|right)|right)$",temp))==0)
                alt.legend <- FALSE
        }
        if(is.logical(alt.legend)){
            if(!alt.legend)
                alt.legend="topright"
            else
                alt.legend="right"
        }
        
        color <- rep(1,(mce$n.comparators+1)); lwd <- 1
        if (mce$n.comparators>7) {
            cl <- colors()
            color <- cl[floor(seq(262,340,length.out=mce$n.comparators))]	# gray scale
            lwd <- 1.5
        }
        
        plot(mce$k,mce$m.ce[,1],t="l",col=color[1],lwd=lwd,lty=1,xlab="Willingness to pay",
             ylab="Probability of most cost effectiveness",ylim=c(0,1),
             main="Cost-effectiveness acceptability curve \nfor multiple comparisons")
        for (i in 2:mce$n.comparators) {
            points(mce$k,mce$m.ce[,i],t="l",col=color[i],lwd=lwd,lty=i)
        }
        legend(alt.legend,mce$interventions,col=color,cex=.7,bty="n",lty=1:mce$n.comparators)
    } # base graphics
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            message("Falling back to base graphics\n")
            mce.plot(mce,pos=pos,graph="base")
            return(invisible(NULL))
        }
        
        if(isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))) {
            # no visible bindings note
            ceplane <- k <- ce <- comp <- NA_real_
            
            alt.legend <- pos
            lty <- rep(1:6,ceiling(mce$n.comparators/6))[1:mce$n.comparators]
            label <- paste0(mce$interventions)
            
            df <- cbind("k"=rep(mce$k,mce$n.comparators),"ce"=c(mce$m.ce))
            df <- data.frame(df,"comp"=as.factor(sort(rep(1:mce$n.comparators,length(mce$k)))))
            names(df) <- c("k","ce","comp")
            
            mceplot <- ggplot2::ggplot(df,ggplot2::aes(x=k,y=ce)) + ggplot2::theme_bw() +
                ggplot2::geom_line(ggplot2::aes(linetype=comp)) + 
                ggplot2::scale_linetype_manual("",labels=label,values=lty) +
                ggplot2::labs(title="Cost-effectiveness acceptability curve\nfor multiple comparisons",x="Willingness to pay",y="Probability of most cost effectiveness") +
                ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"lines"),
                               legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),
                               legend.key=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"))
            
            jus <- NULL
            if(isTRUE(alt.legend)) {
                alt.legend="bottom"
                mceplot <- mceplot + ggplot2::theme(legend.direction="vertical")
            }
            else{
                if(is.character(alt.legend)) {
                    choices <- c("left", "right", "bottom", "top")
                    alt.legend <- choices[pmatch(alt.legend,choices)]
                    jus="center"
                    if(is.na(alt.legend))
                        alt.legend=FALSE
                }
                if(length(alt.legend)>1)
                    jus <- alt.legend
                if(length(alt.legend)==1 & !is.character(alt.legend)) {
                    alt.legend <- c(1,0.5)
                    jus <- alt.legend
                }
            }
            
            mceplot <- mceplot + ggplot2::coord_cartesian(ylim=c(-0.05,1.05)) +
                ggplot2::theme(legend.position=alt.legend,legend.justification=jus,legend.title=ggplot2::element_blank(),
                               legend.background=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold"),
                               legend.text.align=0,plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
            return(mceplot)
        }
    }
}

#######################ceaf.plot##################################
ceaf.plot <- function(mce,graph=c("base","ggplot2")){
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE) 
    if(base.graphics) {
        plot(mce$k,mce$ceaf,t="l",lty=1,
             ylim=c(0,1),xlab="Willingness to pay",
             ylab="Probability of most cost effectiveness",
             main="Cost-effectiveness acceptability frontier")
    }
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("Falling back to base graphics\n")
            ceaf.plot(mce,graph="base")
            return(invisible(NULL))
        }
        
        # no visible binding note
        k  <- NA_real_
        
        df <- data.frame("k"=mce$k,"ceaf"=mce$ceaf)
        ceaf <- ggplot2::ggplot(df,ggplot2::aes(x=k,y=ceaf)) + ggplot2::theme_bw() +
            ggplot2::geom_line() + ggplot2::coord_cartesian(ylim=c(-0.05,1.05)) +
            ggplot2::theme(text=ggplot2::element_text(size=11),legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),plot.title=ggplot2::element_text(face="bold")) +
            ggplot2::labs(title="Cost-effectiveness acceptability frontier",x="Willingness to pay",y="Probability of most cost-effectiveness") +
            ggplot2::theme(plot.title = ggplot2::element_text(lineheight=1.05, face="bold",size=14.3))
        return(ceaf)
    }
}



######evppi###################################################################################################
evppi <- function(parameter,input,he,N=NA,plot=F,residuals=T,...) UseMethod("evppi") 
### AH Version 3 Dec 2015
evppi.default<-function (parameter, input, he, N = NA, plot = F, residuals = T, ...) {
  #Additional options to determine mesh parameters???
  #inputs must be as a data frame to extract column names for GAM regression  
  if (class(parameter[1]) == "character") {
    parameters <- array()
    for (i in 1:length(parameter)) {
      parameters[i] <- which(colnames(input) == parameter[i])
    }
  }
  else {
    parameters = parameter
    parameter = colnames(input)[parameters]
  }
  if (is.null(N)) {
    N <- he$n.sim
  }
  inputs <- data.frame(input)
  
  # Specifies the value of options 'robust' to NULL (it can only be active if 'method'="INLA")
  robust <- NULL
  
  # Checks if method is specified by the user - can also use sal or so univariate
  exArgs <- list(...)
  
  # Sets GAM regression as the default method for single parameter
  if (length(parameter) == 1 & !exists("method", where = exArgs)) {
    exArgs$method <- "GAM"
  }
  # Sets INLA as the default method for multiparameter, if not specified elsewhere
  if (length(parameter) > 1 & !exists("method", where = exArgs)) {
    exArgs$method <- "INLA"
  }
  if (exists("method", where = exArgs)) {
    # Sadatsafavi et al  
    if (exArgs$method == "sal") {
      method = "Sadatsafavi et al"
      n.blocks = NULL
      if (!exists("n.seps", where = exArgs)) {
        n.seps <- 1
      }
      else {
        n.seps <- exArgs$n.seps
      }
      if (length(parameters) == 1) {
        # Needs to modify a vector "param" with the values from the requested parameter
        d <- he$n.comparators
        n <- he$n.sim
        w <- parameters
        param <- inputs[, w]  # vector of simulations for the relevant parameter
        o <- order(param)
        param <- param[o]     # re-order the parameter vector
        nSegs <- matrix(1, d, d)
        nSegs[1, 2] <- n.seps
        nSegs[2, 1] <- n.seps
        res <- segPoints <- numeric()
        for (k in 1:length(he$k)) {
          nbs <- he$U[, k, ]  # this is a n.sims x n.interventions matrix 
                              # with the NB for each value of k
          nbs <- nbs[o, ]     # re-arrange according to the parameter order
          for (i in 1:(d - 1)) {
            for (j in (i + 1):d) {
              cm <- cumsum(nbs[, i] - nbs[, j])/n
              if (nSegs[i, j] == 1) {
                l <- which.min(cm)  # lower bound
                u <- which.max(cm)  # upper bound
                if (cm[u] - max(cm[1], cm[n]) > min(cm[1],cm[n]) - cm[l]) {
                  segPoint <- u
                }
                else {
                  segPoint <- l
                }
                if (segPoint > 1 && segPoint < n) {
                  segPoints <- c(segPoints, segPoint)
                }
              }
              if (nSegs[i, j] == 2) {
                distMaxMin <- 0
                distMinMax <- 0
                minL <- Inf
                maxL <- -Inf
                for (sims in 1:n) {
                  if (cm[sims] > maxL) {
                    maxLP <- sims
                    maxL <- cm[sims]
                  }
                  else {
                    if (maxL - cm[sims] > distMaxMin) {
                      distMaxMin <- maxL - cm[sims]
                      segMaxMinL <- maxLP
                      segMaxMinR <- sims
                    }
                  }
                  if (cm[sims] < minL) {
                    minLP <- sims
                    minL <- cm[sims]
                  }
                  else {
                    if (cm[sims] - minL > distMinMax) {
                      distMinMax <- cm[sims] - minL
                      segMinMaxL <- minLP
                      segMinMaxR <- sims
                    }
                  }
                }
                siMaxMin <- cm[segMaxMinL] + distMaxMin + 
                  (cm[n] - cm[segMaxMinR])
                siMinMax <- -cm[segMaxMinL] + distMinMax - 
                  (cm[n] - cm[segMinMaxR])
                if (siMaxMin > siMinMax) {
                  segPoint <- c(segMaxMinL, segMaxMinR)
                }
                else {
                  segPoint <- c(segMinMaxL, segMinMaxR)
                }
                if (segPoint[1] > 1 && segPoint[1] < 
                    n) {
                  segPoints <- c(segPoints, segPoint[1])
                }
                if (segPoint[2] > 1 && segPoint[2] < 
                    n) {
                  segPoints <- c(segPoints, segPoint[2])
                }
              }
            }
          }
          if (length(segPoints) > 0) {
            segPoints2 <- unique(c(0, segPoints[order(segPoints)], 
                                   n))
            res[k] <- 0
            for (j in 1:(length(segPoints2) - 1)) {
              res[k] <- res[k] + max(colSums(matrix(nbs[(1 + 
                                                           segPoints2[j]):segPoints2[j + 1], ], 
                                                    ncol = d)))/n
            }
            res[k] <- res[k] - max(colMeans(nbs))
          }
          else {
            res[k] <- 0
          }
        }
      }
      if (length(parameters) > 1) {
        res <- list()
        for (lp in 1:length(parameters)) {
          d <- he$n.comparators
          n <- he$n.sim
          w <- parameters[lp]
          param <- inputs[, w]
          o <- order(param)
          param <- param[o]
          nSegs <- matrix(1, d, d)
          nSegs[1, 2] <- n.seps
          nSegs[2, 1] <- n.seps
          temp <- segPoints <- numeric()
          for (k in 1:length(he$k)) {
            nbs <- he$U[, k, ]
            nbs <- nbs[o, ]
            for (i in 1:(d - 1)) {
              for (j in (i + 1):d) {
                cm <- cumsum(nbs[, i] - nbs[, j])/n
                if (nSegs[i, j] == 1) {
                  l <- which.min(cm)
                  u <- which.max(cm)
                  if (cm[u] - max(cm[1], cm[n]) > min(cm[1], 
                                                      cm[n]) - cm[l]) {
                    segPoint <- u
                  }
                  else {
                    segPoint <- l
                  }
                  if (segPoint > 1 && segPoint < n) {
                    segPoints <- c(segPoints, segPoint)
                  }
                }
                if (nSegs[i, j] == 2) {
                  distMaxMin <- 0
                  distMinMax <- 0
                  minL <- Inf
                  maxL <- -Inf
                  for (sims in 1:n) {
                    if (cm[sims] > maxL) {
                      maxLP <- sims
                      maxL <- cm[sims]
                    }
                    else {
                      if (maxL - cm[sims] > distMaxMin) {
                        distMaxMin <- maxL - cm[sims]
                        segMaxMinL <- maxLP
                        segMaxMinR <- sims
                      }
                    }
                    if (cm[sims] < minL) {
                      minLP <- sims
                      minL <- cm[sims]
                    }
                    else {
                      if (cm[sims] - minL > distMinMax) {
                        distMinMax <- cm[sims] - minL
                        segMinMaxL <- minLP
                        segMinMaxR <- sims
                      }
                    }
                  }
                  siMaxMin <- cm[segMaxMinL] + distMaxMin + 
                    (cm[n] - cm[segMaxMinR])
                  siMinMax <- -cm[segMaxMinL] + distMinMax - 
                    (cm[n] - cm[segMinMaxR])
                  if (siMaxMin > siMinMax) {
                    segPoint <- c(segMaxMinL, segMaxMinR)
                  }
                  else {
                    segPoint <- c(segMinMaxL, segMinMaxR)
                  }
                  if (segPoint[1] > 1 && segPoint[1] < 
                      n) {
                    segPoints <- c(segPoints, segPoint[1])
                  }
                  if (segPoint[2] > 1 && segPoint[2] < 
                      n) {
                    segPoints <- c(segPoints, segPoint[2])
                  }
                }
              }
            }
            if (length(segPoints) > 0) {
              segPoints2 <- unique(c(0, segPoints[order(segPoints)], 
                                     n))
              temp[k] <- 0
              for (j in 1:(length(segPoints2) - 1)) {
                temp[k] <- temp[k] + max(colSums(matrix(nbs[(1 + 
                                                               segPoints2[j]):segPoints2[j + 1], ], 
                                                        ncol = d)))/n
              }
              temp[k] <- temp[k] - max(colMeans(nbs))
            }
            else {
              temp[k] <- 0
            }
          }
          res[[lp]] <- temp
        }
        names(res) <- parameters
      }
      
      #Creates names for plot.evppi function
      if (class(parameter) == "numeric") {
        name = colnames(input)[parameter]
      }
      else {
        name = parameter
      }
      res <- list(evppi = res, index = parameters, parameters = name, 
                  k = he$k, evi = he$evi, method = method)
    }
    
    ## Strong and Oakley (univariate)
    if (exArgs$method == "so") {
      method = "Strong & Oakley (univariate)"
      n.seps = NULL
      if (!exists("n.blocks", where = exArgs)) {
        stop("Please specify the parameter 'n.blocks' to use the Strong and Oakley univariate method")
      }
      else {
        n.blocks <- exArgs$n.blocks
      }
      S <- he$n.sim                 # number of samples
      J <- S/exArgs$n.blocks        # constrains the parameters J,K to have product
                                    # equal to the number of simulations
      check <- S%%exArgs$n.blocks   # checks that Jxn.blocks = S with no remainder
      if (check > 0) {
        stop("number of simulations/number of blocks must be an integer. Please select a different value for n.blocks \n")
      }
      D <- he$n.comparators         # number of decision options
      # Checks how many parameters have been specified
      if (length(parameter) == 1) {
        sort.order <- order(inputs[, parameters])
        sort.U <- array(NA, dim(he$U))
        evpi <- res <- numeric()
        for (i in 1:length(he$k)) {
          evpi[i] <- he$evi[i]
          sort.U[, i, ] <- he$U[sort.order, i, ]
          U.array <- array(sort.U[, i, ], dim = c(J, 
                                                  exArgs$n.blocks, D))
          mean.k <- apply(U.array, c(2, 3), mean)
          partial.info <- mean(apply(mean.k, 1, max))
          res[i] <- partial.info - max(apply(he$U[, i, 
                                                  ], 2, mean))
        }
      }
      
      if (length(parameter) > 1) {
        res <- list()
        for (j in 1:length(parameter)) {
          sort.order <- order(inputs[, parameters[j]])
          sort.U <- array(NA, dim(he$U))
          evpi <- evppi.temp <- numeric()
          for (i in 1:length(he$k)) {
            evpi[i] <- he$evi[i]
            sort.U[, i, ] <- he$U[sort.order, i, ]
            U.array <- array(sort.U[, i, ], dim = c(J, 
                                                    n.blocks, D))
            mean.k <- apply(U.array, c(2, 3), mean)
            partial.info <- mean(apply(mean.k, 1, max))
            evppi.temp[i] <- partial.info - max(apply(he$U[, 
                                                           i, ], 2, mean))
          }
          res[[j]] <- evppi.temp
        }
        names(res) <- parameters
      }
      
      if (class(parameter) == "numeric") {
        name = colnames(input)[parameter]
      }
      else {
        name = parameter
      }
      res <- list(evppi = res, index = parameters, parameters = name, 
                  k = he$k, evi = he$evi, method = method)
    }
    
    if(exArgs$method == "INLA" || exArgs$method == "GAM"||exArgs$method == "gam"||exArgs$method == "g"||exArgs$method == "G"
       || exArgs$method == "GP"||exArgs$method == "gp"){
      
      ### Functions to Compute the EVPPI using GP regression-like methods #########################
      
      ## 1. Fit model using GAM
      fit.gam <- function(parameter, input, x, select, formula) {
        # Runs GAM regression
        # parameter = A vector of parameters for which the EVPPI should be calculated. 
        #             This can be given as a string (or vector of strings) of names or 
        #             a numeric vector, corresponding to the column numbers of important parameters.
        # input = A matrix containing the simulations for all the parameters monitored by the call to 
        #         JAGS or BUGS. The matrix should have column names matching the names of the parameters 
        #         and the values in the vector parameter should match at least one of those values.
        # x = a BCEA object with the PSA runs and relevant information on the comparison to be made
        # N = The number of PSA simulations used to calculate the EVPPI. The default uses all the 
        #     available samples.
        ### optional argument are 
        # formula: a string specifying the form of the non parametric regression in terms of the
        #          parameters involved. "Main" effects are constructed with the format s(par),
        #          where par is one of the parameters for the model being considered. Also, it
        #          is possible to consider "interaction" terms, to account for correlation 
        #          among parameters. This is done using the notation te(par1,par2), which 
        #          defines the tensor function between par1 and par2. The default is using the 
        #          tensor and cubic splines, which speed up the computation
        
        tic <- proc.time()
        d <- x$n.comparisons
        N <- min(x$n.sim, N, na.rm = T)
        if (N == x$n.sim) {
          select <- 1:x$n.sim
        }
        else {
          select <- sample(1:x$n.sim, size = N, replace = F)
        }
        inp <- names(inputs)[parameters]
        if (d == 1) {
          f.e <- update(formula(x$delta.e[select] ~ .), 
                        formula(paste(".~", form)))
          f.c <- update(formula(x$delta.c[select] ~ .), 
                        formula(paste(".~", form)))
          model.e <- mgcv::gam(f.e, data = data.frame(inputs)[select, 
                                                              ])
          model.c <- mgcv::gam(f.c, data = data.frame(inputs)[select, 
                                                              ])
          e.hat <- model.e$fitted
          c.hat <- model.c$fitted
        }
        else {
          e.hat <- c.hat <- matrix(NA, N, d)
          model.e <- model.c <- list()
          for (i in 1:d) {
            f.e <- update(formula(x$delta.e[select, i] ~ 
                                    .), formula(paste(".~", form)))
            f.c <- update(formula(x$delta.c[select, i] ~ 
                                    .), formula(paste(".~", form)))
            model.e[[i]] <- mgcv::gam(f.e, data = data.frame(inputs)[select, 
                                                                     ])
            model.c[[i]] <- mgcv::gam(f.c, data = data.frame(inputs)[select, 
                                                                     ])
            e.hat[, i] <- model.e[[i]]$fitted
            c.hat[, i] <- model.c[[i]]$fitted
          }
        }
        formula <- form
        fitted.costs <- matrix(c.hat, nrow = N, ncol = d)
        fitted.effects <- matrix(e.hat, nrow = N, ncol = d)
        fit.c <- model.c
        fit.e <- model.e
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to fit GAM regression (seconds)"
        
        # defines output of the function
        list(fitted.costs = fitted.costs, fitted.effects = fitted.effects, 
             formula = formula, fit.c = fit.c, fit.e = fit.e, 
             time = time)
      }
      
      
      post.density <- function(hyperparams, obj, input, parameter) {
        ## 2. Fit model using GP Regression
        ## Code by Mark Strong
        ## This function computes the log posterior density of the GP hyperparameters
        
        dinvgamma <- function(x, alpha, beta) {
          (beta^alpha)/gamma(alpha) * x^(-alpha - 1) * 
            exp(-beta/x)
        }
        if (class(parameter[1]) == "character") {
          parameters <- array()
          for (i in 1:length(parameter)) {
            parameters[i] <- which(colnames(input) == parameter[i])
          }
        }
        else {
          parameters = parameter
          parameter = colnames(input)[parameters]
        }
        input.matrix <- as.matrix(input[, parameters, drop = FALSE])
        colmin <- apply(input.matrix, 2, min)
        colmax <- apply(input.matrix, 2, max)
        colrange <- colmax - colmin
        input.matrix <- sweep(input.matrix, 2, colmin, "-")
        input.matrix <- sweep(input.matrix, 2, colrange, 
                              "/")
        N <- nrow(input.matrix)
        p <- ncol(input.matrix)
        H <- cbind(1, input.matrix)
        q <- ncol(H)
        a.sigma <- 0.001    ##  hyperparameters for IG prior for sigma^2
        b.sigma <- 0.001    ##  hyperparameters for IG prior for nu
        a.nu <- 0.001
        b.nu <- 1
        delta <- exp(hyperparams)[1:p]
        nu <- exp(hyperparams)[p + 1]
        A <- exp(-(as.matrix(dist(t(t(input.matrix)/delta), 
                                  upper = TRUE, diag = TRUE))^2))
        Astar <- A + nu * diag(N)
        T <- chol(Astar)
        y <- backsolve(t(T), obj, upper.tri = FALSE)
        x <- backsolve(t(T), H, upper.tri = FALSE)
        tHAstarinvH <- t(x) %*% (x)
        betahat <- solve(tHAstarinvH) %*% t(x) %*% y
        residSS <- y %*% y - t(y) %*% x %*% betahat - t(betahat) %*% 
          t(x) %*% y + t(betahat) %*% tHAstarinvH %*% betahat
        prior <- prod(dnorm(log(delta), 0, sqrt(100000))) * 
          dinvgamma(nu, a.nu, b.nu)
        l <- -sum(log(diag(T))) - 1/2 * log(det(tHAstarinvH)) - 
          (N - q + 2 * a.sigma)/2 * log(residSS/2 + b.sigma) + 
          log(prior)
        names(l) <- NULL
        return(l)
      }
      
      ## This function estimates the GP hyperparameters numerically using optim
      ## m is the number of PSA samples used to estimate the hyperparameters
      estimate.hyperparameters <- function(obj, input, parameter,n.sim) {
        # obj = either Delta.e or Delta.c as called by the rest of the function
        # input = the matrix (or data frame) with the PSA runs
        # parameter = a (possibly numeric) vector with the parameters of interest
        # n.sim = number of simulations to be performed to estimate the hyperparameters
        
        p <- length(parameter)
        D <- he$n.comparators
        hyperparameters <- vector("list", D)
        hyperparameters[[1]] <- NA
        for (i in 2:D) {
          initial.values <- rep(0, p + 1)
          repeat {
            log.hyperparameters <- optim(initial.values, 
                                         fn = post.density, obj = obj[1:n.sim, i], 
                                         input = input[1:n.sim, ], parameter = parameter, 
                                         method = "Nelder-Mead", control = list(fnscale = -1, 
                                                                                maxit = 10000, trace = 0))$par
            if (sum(abs(initial.values - log.hyperparameters)) < 
                0.01) {
              hyperparameters[[i]] <- exp(log.hyperparameters)
              break
            }
            initial.values <- log.hyperparameters
          }
        }
        return(hyperparameters)
      }
      
      fit.gp <- function(parameter, input, x, select, n.sim) {
        # Runs GP regression
        # parameter = A vector of parameters for which the EVPPI should be calculated. 
        #             This can be given as a string (or vector of strings) of names or 
        #             a numeric vector, corresponding to the column numbers of important parameters.
        # input = A matrix containing the simulations for all the parameters monitored by the call to 
        #         JAGS or BUGS. The matrix should have column names matching the names of the parameters 
        #         and the values in the vector parameter should match at least one of those values.
        # x = a BCEA object with the PSA runs and relevant information on the comparison to be made
        # select = The number of PSA simulations used to calculate the EVPPI. The default uses all the 
        #     available samples.
        # n.sim = number of simulations to compute the hyperparameters
        ### optional argument are 
        # formula: a string specifying the form of the non parametric regression in terms of the
        #          parameters involved. "Main" effects are constructed with the format s(par),
        #          where par is one of the parameters for the model being considered. Also, it
        #          is possible to consider "interaction" terms, to account for correlation 
        #          among parameters. This is done using the notation te(par1,par2), which 
        #          defines the tensor function between par1 and par2. The default is using the 
        #          tensor and cubic splines, which speed up the computation
        
        tic <- proc.time()
        if (!is.null(dim(x$e))) {
          de <- x$e[select, ] - x$e[select, 1]
          dc <- x$c[select, ] - x$c[select, 1]
        }
        else {
          de <- cbind(0, x$e[select])
          dc <- cbind(0, x$c[select])
        }
        p <- length(parameter)
        D <- x$n.comparators
        input.matrix <- as.matrix(input[select, parameters, 
                                        drop = FALSE])
        colmin <- apply(input.matrix, 2, min)
        colmax <- apply(input.matrix, 2, max)
        colrange <- colmax - colmin
        input.matrix <- sweep(input.matrix, 2, colmin, "-")
        input.matrix <- sweep(input.matrix, 2, colrange, 
                              "/")
        H <- cbind(1, input.matrix)
        q <- ncol(H)
        V.e <- V.c <- g.hat.c <- g.hat.e <- vector("list", 
                                                   D)
        g.hat.e[[1]] <- g.hat.c[[1]] <- rep(0, N)
        for (d in 2:D) {
          hyperparameters <- estimate.hyperparameters(obj = de, 
                                                      input = input, parameter = parameter, n.sim = n.sim)
          delta.hat <- hyperparameters[[d]][1:p]
          nu.hat <- hyperparameters[[d]][p + 1]
          A <- exp(-(as.matrix(dist(t(t(input.matrix)/delta.hat), 
                                    upper = TRUE, diag = TRUE))^2))
          Astar <- A + nu.hat * diag(N)
          Astarinv <- chol2inv(chol(Astar))
          rm(Astar)
          gc()
          AstarinvY <- Astarinv %*% de[, d]
          tHAstarinv <- t(H) %*% Astarinv
          tHAHinv <- solve(tHAstarinv %*% H)
          betahat <- tHAHinv %*% (tHAstarinv %*% de[, d])
          Hbetahat <- H %*% betahat
          resid <- de[, d] - Hbetahat
          g.hat.e[[d]] <- Hbetahat + A %*% (Astarinv %*% 
                                              resid)
          AAstarinvH <- A %*% t(tHAstarinv)
          sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*% 
                                     resid)/(N - q - 2)
          V.e[[d]] <- sigmasqhat * (nu.hat * diag(N) - 
                                      nu.hat^2 * Astarinv + (H - AAstarinvH) %*% 
                                      (tHAHinv %*% t(H - AAstarinvH)))
          rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv, 
             Hbetahat, resid, sigmasqhat)
          gc()
          hyperparameters <- estimate.hyperparameters(obj = dc, 
                                                      input = input, parameter = parameter, n.sim = n.sim)
          delta.hat <- hyperparameters[[d]][1:p]
          nu.hat <- hyperparameters[[d]][p + 1]
          A <- exp(-(as.matrix(dist(t(t(input.matrix)/delta.hat), 
                                    upper = TRUE, diag = TRUE))^2))
          Astar <- A + nu.hat * diag(N)
          Astarinv <- chol2inv(chol(Astar))
          rm(Astar)
          gc()
          AstarinvY <- Astarinv %*% dc[, d]
          tHAstarinv <- t(H) %*% Astarinv
          tHAHinv <- solve(tHAstarinv %*% H)
          betahat <- tHAHinv %*% (tHAstarinv %*% dc[, d])
          Hbetahat <- H %*% betahat
          resid <- dc[, d] - Hbetahat
          g.hat.c[[d]] <- Hbetahat + A %*% (Astarinv %*% 
                                              resid)
          AAstarinvH <- A %*% t(tHAstarinv)
          sigmasqhat <- as.numeric(t(resid) %*% Astarinv %*% 
                                     resid)/(N - q - 2)
          V.c[[d]] <- sigmasqhat * (nu.hat * diag(N) - 
                                      nu.hat^2 * Astarinv + (H - AAstarinvH) %*% 
                                      (tHAHinv %*% t(H - AAstarinvH)))
          rm(A, Astarinv, AstarinvY, tHAstarinv, tHAHinv, 
             Hbetahat, resid, sigmasqhat)
          gc()
        }
        g.hat.e <- matrix(unlist(g.hat.e), nrow = N, ncol = D)
        g.hat.c <- matrix(unlist(g.hat.c), nrow = N, ncol = D)
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to fit GP regression (seconds)"
        fitted.costs = matrix(g.hat.c[, 2:D], nrow = N, ncol = (D - 1))
        fitted.effects = matrix(g.hat.e[, 2:D], nrow = N, ncol = (D - 1))
        
        # defines output of the function
        list(fitted.costs = fitted.costs, fitted.effects = fitted.effects, 
             time = time, formula = NULL)
      }
      
      ## 3. Fit model using INLA
      make.proj <- function(inputs, parameters, select) {
        #Set up the projection from P dimesional space to 2 dimensions
        
        tic <- proc.time()
        comp1 <- prcomp(inputs[select, parameters])
        Datum <- inputs[select, parameters]
        Data <- cbind(as.matrix(Datum) %*% as.matrix(comp1$rotation[,2]), 
                      as.matrix(Datum) %*% as.matrix(comp1$rotation[,1]))
        Data.Scale <- cbind((Data[, 1] - mean(Data[, 1]))/sd(Data[,1]), (Data[, 2] - mean(Data[, 2]))/sd(Data[,2]))
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to fit find projections (seconds)"
        list(Data.Scale = Data.Scale, time = time)
      }
      
      #Plotting shows the mesh
      plot.mesh <- function(mesh, data, plot) {
        if (plot == TRUE || plot == T) {
          cat("\n")
          choice <- select.list(c("yes", "no"), title = "Would you like to save the graph?", 
                                graphics = F)
          if (choice == "yes") {
            exts <- c("jpeg", "pdf", "bmp", "png", "tiff")
            ext <- select.list(exts, title = "Please select file extension", 
                               graphics = F)
            name <- paste0(getwd(), "/mesh.", ext)
            txt <- paste0(ext, "('", name, "')")
            eval(parse(text = txt))
            plot(mesh)
            points(data, col = "blue", pch = 19, cex = 0.8)
            dev.off()
            txt <- paste0("Graph saved as: ", name)
            cat(txt)
            cat("\n")
          }
          cat("\n")
          plot(mesh)
          points(data, col = "blue", pch = 19, cex = 0.8)
        }
      }
      
      # Makes the mesh
      make.mesh <- function(data, convex.inner, convex.outer,cutoff) {
        tic <- proc.time()
        inner <- suppressMessages({
          INLA::inla.nonconvex.hull(data, convex = convex.inner)
        })
        outer <- INLA::inla.nonconvex.hull(data, convex = convex.outer)
        mesh.tmp <- INLA::inla.mesh.2d(boundary = list(inner), 
                                       max.edge = c(0.5), cutoff = cutoff/2)
        mesh <- INLA::inla.mesh.2d(loc = mesh.tmp$loc, boundary = list(INLA::inla.mesh.boundary(mesh.tmp), 
                                                                       outer), max.edge = c(1.2, 1.2 * 2), cutoff = cutoff)
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to fit determine the mesh (seconds)"
        list(mesh = mesh, pts = data, time = time)
      }
      
      # Fits INLA
      fit.inla <- function(parameters, inputs, x, select, mesh, 
                           data.scale, int.ord, convex.inner, convex.outer, 
                           cutoff, h.value) {
        
        tic <- proc.time()
        inputs.Scale <- scale(inputs[select,], apply(inputs[select,], 2, mean), 
                              apply(inputs[select, ], 2, sd))
        d <- x$n.comparators
        e.I <- fit.e <- list()
        e.scale <- array()
        fitted.effects <- matrix(NA, nrow = length(select),ncol = (d - 1))
        c.I <- fit.c <- list()
        c.scale <- array()
        fitted.costs <- matrix(NA, nrow = length(select),ncol = (d - 1))
        for (i in 1:(d - 1)) {
          g <- x$comp[i]
          
          ## Run INLA/SPDE for effects
          e.I[[i]] <- as.matrix(x$delta.e)[select, i]
          e.scale[i] <- 6/(range(e.I[[i]])[2] - range(e.I[[i]])[1])
          A <- INLA::inla.spde.make.A(mesh = mesh, loc = data.scale, silent = 2L)
          #Sets up an spde object to be used in the INLA call
          spde <- INLA::inla.spde2.matern(mesh = mesh,alpha = 2)
          #Data for the SPDE
          stk.real <- INLA::inla.stack(tag = "est", # tag
                                       data = list(y = e.I[[i]] *e.scale[i]), # response
                                       A = list(A, 1),  # two-projector matrix
                                       effects = list(  # two elements
                                         s = 1:spde$n.spde,
                                         data.frame(b0 = 1,x = cbind(data.scale,inputs.Scale))))
          
          #Using INLA to find the fitted values 
          data <- INLA::inla.stack.data(stk.real)
          ctr.pred <- INLA::inla.stack.A(stk.real)
          
          # Creates the formula
          inp <- names(stk.real$effects$data)[parameters + 4]
          form <- paste(inp, "+", sep = "", collapse = "")
          formula <- paste("y~0+(", form, "+0)+b0+f(s,model=spde)", 
                           sep = "", collapse = "")
          if (int.ord[1] > 1) {
            formula <- paste("y~0+(", form, "+0)^", int.ord[1], 
                             "+b0+f(s,model=spde)", sep = "", collapse = "")
          }
          
          # NB Matrix will throw a NOTE --- but suppressMessages will avoid this
          # NB2: Should we compute the DIC so that we could assess what model is the better??
          Result <- suppressMessages({
            INLA::inla(as.formula(formula), data = data, 
                       family = family, control.predictor = list(A = ctr.pred, 
                                                                 link = 1), control.inla = list(h = h.value), 
                       control.compute = list(config = T))
          })
          rescaled <- Result$summary.linear.predictor[1:length(select), 
                                                      "mean"]/e.scale[i]
          fitted.effects[, i] <- rescaled
          fit.e[[i]] <- Result
          
          ## Run INLA/SPDE for costs
          c.I[[i]] <- as.matrix(x$delta.c)[select, i]
          c.scale[i] <- 6/(range(c.I[[i]])[2] - range(c.I[[i]])[1])
          A <- INLA::inla.spde.make.A(mesh = mesh, loc = data.scale, silent = 2L)
          #Sets up an spde object to be used in the INLA call
          spde <- INLA::inla.spde2.matern(mesh = mesh, 
                                          alpha = 2)
          stk.real <- INLA::inla.stack(tag = "est", 
                                       data = list(y = c.I[[i]] * c.scale[i]), 
                                       A = list(A, 1), 
                                       effects = list(s = 1:spde$n.spde, 
                                                      data.frame(b0 = 1, x = cbind(data.scale, inputs.Scale))))
          
          data <- INLA::inla.stack.data(stk.real)
          ctr.pred <- INLA::inla.stack.A(stk.real)
          
          # Creates the formula
          formula <- paste("y~0+(", form, "+0)+b0+f(s,model=spde)", 
                           sep = "", collapse = "")
          if (int.ord[2] > 1) {
            formula <- paste("y~0+(", form, "+0)^", int.ord[2], 
                             "+b0+f(s,model=spde)", sep = "", collapse = "")
          }
          Result <- suppressMessages({
            INLA::inla(as.formula(formula), data = data, 
                       family = family, control.predictor = list(A = ctr.pred, 
                                                                 link = 1), control.inla = list(h = h.value), 
                       control.compute = list(config = T))
          })
          rescaled <- Result$summary.linear.predictor[1:length(select), 
                                                      "mean"]/c.scale[i]
          fitted.costs[, i] <- rescaled
          fit.c[[i]] <- Result
        }
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to fit INLA/SPDE (seconds)"
        
        list(fitted.costs = fitted.costs, fitted.effects = fitted.effects, 
             fit.c = fit.c, fit.e = fit.e, time = time, formula = formula, 
             mesh = list(mesh = mesh, pts = data.scale))
      }
      
      ## 4. Compute the EVPPI
      compute.evppi <- function(x, fit, method) {
        # x = bcea object
        # fit = the result of the model fitting routine
        # method = a string specifying the method used to fit the model
        
        EVPPI <- array()
        tic <- proc.time()
        for (i in 1:length(he$k)) {
          if (method == "GP") {
            NB.k.mid <- -(x$k[i] * fit$fitted.effects - 
                            fit$fitted.costs)
            NB.k <- cbind(NB.k.mid, rep(0, N))
            EVPPI[i] <- (mean(apply(NB.k, 1, max, na.rm = T)) - 
                           max(apply(NB.k, 2, mean, na.rm = T)))
          }
          else {
            NB.k.mid <- -(x$k[i] * fit$fitted.effects - 
                            fit$fitted.costs)
            NB.k <- cbind(NB.k.mid, rep(0, N))
            EVPPI[i] <- (mean(apply(NB.k, 1, max, na.rm = T)) - 
                           max(apply(NB.k, 2, mean, na.rm = T)))
          }
        }
        toc <- proc.time() - tic
        time <- toc[3]
        names(time) = "Time to compute the EVPPI (in seconds)"
        # Outputs of the function
        list(EVPPI = EVPPI, time = time)
      }
      
      ## 5. Creates names for the evppi.plot function
      #Creates names for plot.evppi function
      prepare.output <- function(parameter, input) {
        if (length(parameter) == 1) {
          if (class(parameter) == "numeric") {
            name = colnames(input)[parameter]
          }
          else {
            name = parameter
          }
        }
        else {
          if (class(parameter) == "numeric") {
            n.param <- length(parameter)
            end <- colnames(input)[parameter[n.param]]
            name.mid <- paste(colnames(input)[parameter[1:n.param - 
                                                          1]], ", ", sep = "", collapse = " ")
            name <- paste(name.mid, "and ", end, sep = "", 
                          collapse = " ")
          }
          else {
            n.param <- length(parameter)
            end <- parameter[n.param]
            name.mid <- paste(parameter[1:n.param - 1], 
                              ", ", sep = "", collapse = " ")
            name <- paste(name.mid, "and ", end, sep = "", 
                          collapse = " ")
          }
        }
        return(name)
      }
      ###############################################################################
      
      ## GAM regression
      if (exArgs$method == "GAM" || exArgs$method == "gam" || 
          exArgs$method == "G" || exArgs$method == "g") {
        method <- "GAM"
        mesh <- robust <- NULL
        
        # If mgcv is not installed, then asks for it
        if (!isTRUE(requireNamespace("mgcv", quietly = TRUE))) {
          stop("You need to install the package 'mgcv'. Please run in your R terminal:\n install.packages('mgcv')")
        }
        # If it is installed, then uses its namespace and all should work!
        if (isTRUE(requireNamespace("mgcv", quietly = TRUE))) {
          cat("\n")
          cat("Calculating fitted values for the GAM regression \n")
          if (any(!is.na(N)) & length(N) > 1) {
            select <- N
          }
          else {
            N <- min(he$n.sim, N, na.rm = T)
            if (N == he$n.sim) {
              select <- 1:he$n.sim
            }
            else {
              select <- sample(1:he$n.sim, size = N, replace = F)
            }
          }
          inp <- names(inputs)[parameters]
          # If the user wants to use the default case
          if (exists("formula", where = exArgs)) {
            form <- exArgs$formula
          }
          else {
            form <- paste("te(", paste(inp, ",", sep = "", 
                                       collapse = ""), "bs='cr')")
          }
          # Fit the GAM regression
          fit <- fit.gam(parameter = parameters, input = input, 
                         x = he, select = select, formula = form)
        }
      }
      
      ## GP Regression
      if (exArgs$method == "gp" || exArgs$method == "GP") {
        method <- "GP"
        mesh <- robust <- NULL
        cat("\n")
        cat("Calculating fitted values for the GP regression \n")
        # Allows the user to select the number of PSA runs used to estimate the hyperparameters
        if (!exists("n.sim", where = exArgs)) {
          n.sim = 500
        }
        else {
          n.sim = exArgs$n.sim
        }
        if (any(!is.na(N)) & length(N) > 1) {
          select <- N
        }
        else {
          N <- min(he$n.sim, N, na.rm = T)
          if (N == he$n.sim) {
            select <- 1:he$n.sim
          }
          else {
            select <- sample(1:he$n.sim, size = N, replace = F)
          }
        }
        fit <- fit.gp(parameter = parameter, input = input, 
                      x = he, select = select, n.sim = n.sim)
      }
      
      ## INLA/SPDE
      if (exArgs$method == "INLA") {
        method <- "INLA"
        # If INLA is not installed, then asks for it
        if (!isTRUE(requireNamespace("INLA", quietly = TRUE))) {
          stop("You need to install the packages 'INLA' and 'splancs'. Please run in your R terminal:\n install.packages('INLA', repos='http://www.math.ntnu.no/inla/R/stable')\n and\n install.packages('splancs')")
        }
        # If it is installed, then uses its namespace and all should work!
        if (isTRUE(requireNamespace("INLA", quietly = TRUE))) {
          if (!is.element("INLA", (.packages()))) {
            attachNamespace("INLA")
          }
          #Restricts the number of PSA samples used to calculate the EVPPI - speeds up computation time
          # Lets the user select either: 1) all the PSA samples in the bcea object (N=NA); or 2) a fixed number N
          # (and will sample N out of he$n.sim values at random); or 3) even a specified vector of values (length(N)>1) 
          if (any(!is.na(N)) & length(N) > 1) {
            select <- N
          }
          else {
            N <- min(he$n.sim, N, na.rm = T)
            if (N == he$n.sim) {
              select <- 1:he$n.sim
            }
            else {
              select <- sample(1:he$n.sim, size = N, replace = F)
            }
          }
          cat("\n")
          cat("Finding projections \n")
          projections <- make.proj(inputs = inputs, parameters = parameters, 
                                   select = select)
          Data.Scale <- projections$Data.Scale
          
          cat("Determining Mesh \n")
          #Find mesh that fits the scaled projections
          # NB: Includes this in a suppressMessage command to prevent R from printing stuff about loading several packages
          # Allows the user to specify the density of the mesh & how far the boundaries are from the data points
          if (!exists("cutoff", where = exArgs)) {
            cutoff = 0.3
          }
          else {
            cutoff = exArgs$cutoff
          }
          if (!exists("convex.inner", where = exArgs)) {
            convex.inner = -0.4
          }
          else {
            convex.inner = exArgs$convex.inner
          }
          if (!exists("convex.outer", where = exArgs)) {
            convex.outer = -0.7
          }
          else {
            convex.outer = exArgs$convex.outer
          }
          mesh <- make.mesh(data = Data.Scale, convex.inner = convex.inner, 
                            convex.outer = convex.outer, cutoff = cutoff)
          # Possibly makes the plot of the mesh
          plot.mesh(mesh = mesh$mesh, data = Data.Scale, 
                    plot = plot)
          
          #Calculating a set of fitted values for each decision option (INB - to reduce fitting)
          cat("Calculating fitted values for the GP regression using INLA/SPDE \n")
          if (exists("h.value", where = exArgs)) {
            h.value = exArgs$h.value
          }
          else {
            h.value = 0.00005
          }
          if (exists("robust", where = exArgs)) {
            if (exArgs$robust == TRUE) {
              family = "T"
              robust = TRUE
            }
            else {
              family = "gaussian"
              robust = FALSE
            }
          }
          else {
            family = "gaussian"
            robust = FALSE
          }
          if (exists("int.ord", where = exArgs)) {
            int.ord = exArgs$int.ord
          }
          else {
            int.ord = c(1, 1)
          }
          
          fit <- fit.inla(parameters = parameters, inputs = inputs, 
                          x = he, select = select, mesh = mesh$mesh, 
                          data.scale = Data.Scale, int.ord = int.ord, 
                          convex.inner = convex.inner, convex.outer = convex.outer, 
                          cutoff = cutoff, h.value = h.value)
        }
      }
      
      # Computes the EVPPI
      cat("Calculating EVPPI \n")
      comp <- compute.evppi(x = he, fit = fit, method = method)
      # Prepares the output for plot.evppi
      name <- prepare.output(parameter, input)
      # Formats computational time 
      if (method == "INLA") {
        time <- c(projections$time, mesh$time, fit$time, 
                  comp$time)
        time[5] <- sum(time)
        time <- as.matrix(time)
        rownames(time) = c("Finding projections", "Determining mesh", 
                           "Model fitting", "EVPPI computation", "Total")
      }
      else {
        time <- c(fit$time, comp$time)
        time[3] <- sum(time)
        time <- as.matrix(time)
        rownames(time) = c("Model fitting", "EVPPI computation", 
                           "Total")
      }
      colnames(time) = "Running time (seconds)"
      # Sets up the outcome
      if (residuals == TRUE || residuals == T) {
        res <- list(evppi = comp$EVPPI, index = parameters, 
                    k = he$k, evi = he$evi, parameters = name, time = time, 
                    formula = fit$formula, method = method, fitted.costs = fit$fitted.costs, 
                    fitted.effects = fit$fitted.effects, fit.e = fit$fit.e, 
                    fit.c = fit$fit.c, mesh = mesh)
      }
      
      else {
        res <- list(evppi = comp$EVPPI, index = parameters, 
                    k = he$k, evi = he$evi, parameters = name, time = time, 
                    formula = fit$formula, method = method)
      }
    }
  }
  ### 
  class(res) <- "evppi"
  return(res)
}


######plot.evppi##############################################################################################
plot.evppi<-function (x, pos = c(0, 0.8), graph = c("base", "ggplot2"), col = NULL, 
                      ...){
  # Plots the EVPI and the EVPPI for all the parameters being monitored
  # x is a "evppi" object obtained by a call to the function evppi
  # col is a vector specifying the colors to be used in the graphs
  #     if null, then all are in black
  options(scipen = 10)
  alt.legend <- pos
  base.graphics <- ifelse(isTRUE(pmatch(graph, c("base", "ggplot2")) == 
                                   2), FALSE, TRUE)
  stopifnot(isTRUE(class(x) == "evppi"))
  if (base.graphics) {
    if (is.numeric(alt.legend) & length(alt.legend) == 2) {
      temp <- ""
      if (alt.legend[2] == 0) 
        temp <- paste0(temp, "bottom")
      else if (alt.legend[2] != 0.5) 
        temp <- paste0(temp, "top")
      if (alt.legend[1] == 1) 
        temp <- paste0(temp, "right")
      else temp <- paste0(temp, "left")
      alt.legend <- temp
      if (length(grep("^((bottom|top)(left|right)|right)$", 
                      temp)) == 0) 
        alt.legend <- FALSE
    }
    if (is.logical(alt.legend)) {
      if (!alt.legend) 
        alt.legend = "topright"
      else alt.legend = "topleft"
    }
    plot(x$k, x$evi, t = "l", xlab = "Willingness to pay", 
         ylab = "", main = "Expected Value of Perfect Partial Information", 
         lwd = 2, ylim = range(range(x$evi), range(x$evppi)))
    if (is.null(col)) {
      cols <- colors()
      gr <- floor(seq(from = 261, to = 336, length.out = length(x$index)))
      col <- cols[gr]
    }
    else {
      if (length(col) != length(x$parameters)) {
        message("The vector 'col' must have the same number of elements as the number of parameters. Forced to black\n")
        col <- rep("black", length(x$parameters))
      }
    }
    if (length(x$index) == 1 | length(x$index) > 1 & (x$method == 
                                                      "INLA" || x$method == "GAM" || x$method == "GP")) {
      col = "black"
      points(x$k, x$evppi, t = "l", col = col, lty = 1)
    }
    cmd <- "EVPPI for the selected\nsubset of parameters"
    if (nchar(x$parameters[1]) <= 25) {
      cmd <- paste("EVPPI for ", x$parameters, sep = "")
    }
    if (length(x$index) > 1 & (x$method == "Strong & Oakley (univariate)" || 
                               x$method == "Sadatsafavi et al")) {
      for (i in 1:length(x$index)) {
        points(x$k, x$evppi[[i]], t = "l", col = col[i], 
               lty = i)
        text(par("usr")[2], x$evppi[[i]][length(x$k)], 
             paste("(", i, ")", sep = ""), cex = 0.7, pos = 2)
      }
      cmd <- paste("(", paste(1:length(x$index)), ") EVPPI for ", 
                   x$parameters, sep = "")
    }
    legend(alt.legend, c("EVPI", cmd), col = c("black", col), 
           cex = 0.7, bty = "n", lty = c(1, 1:length(x$parameters)), 
           lwd = c(2, rep(1, length(x$parameters))))
    return(invisible(NULL))
  }
  else {
    if (!isTRUE(requireNamespace("ggplot2", quietly = TRUE) & 
                requireNamespace("grid", quietly = TRUE))) {
      message("Falling back to base graphics\n")
      plot.evppi(x, pos = c(0, 0.8), graph = "base", col)
      return(invisible(NULL))
    }
    else {
      message("ggplot2 method not yet implemented for this function: falling back to base graphics\n")
      plot.evppi(x, pos = c(0, 0.8), graph = "base", col)
      return(invisible(NULL))
    }
  }
}




######diag.evppi################################################################################################
diag.evppi <- function(x,y,diag=c("residuals","qqplot"),int=1){
    # x = an evppi object
    # y = a bcea object
    # diag = the type of diagnostics required
    # int = the comparison to be assessed (default determined by the BCEA object)
    if (int>1 & dim(x$fitted.costs)[2]==1) {stop("There is only one comparison possible, so 'int' should be set to 1 (default)")}
    res <- ifelse(isTRUE(pmatch(diag,c("residuals","qqplot"))==2),FALSE,TRUE)
    if(res){
        par(mfrow=c(1,2))
        n <- dim(x$fitted.costs)[1]
        fitted <- x$fitted.costs[,int]
        residual <- as.matrix(y$delta.c)[1:n,int]-fitted
        plot(fitted,residual,xlab="Fitted values",
             ylab="Residuals",main="Residual plot for costs",cex=.8);abline(h=0)
        fitted <- x$fitted.effects[,int]
        residual <- as.matrix(y$delta.e)[1:n,int]-fitted
        plot(fitted,residual,xlab="Fitted values",
             ylab="Residuals",main="Residual plot for effects",cex=.8);abline(h=0)
    }else{
        par(mfrow=c(1,2))
        qqnorm(x$fitted.costs[,int],main="Normal Q-Q plot \n(costs)"); qqline(x$fitted.costs[,int])
        qqnorm(x$fitted.effects[,int],main="Normal Q-Q plot \n(effects)"); qqline(x$fitted.effects[,int])
    }
    
}





######CreateInputs##############################################################################################
CreateInputs <- function(x) {
    # Utility function --- creates inputs for the EVPPI
    # First checks whether the model is run with JAGS or BUGS
    # 1. checks whether each model has been run using JAGS or BUGS    
    if(class(x)=="rjags") {
        if(!isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
            stop("You need to install the package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
        }
        if (isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
            mdl <- x$BUGSoutput      
        }
    }
    if(class(x)=="bugs") {  	# if model is run using R2WinBUGS/R2OpenBUGS
        if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
            stop("You need to install the package 'R2OpenBUGS'. Please run in your R terminal:\n install.packages('R2OpenBUGS')")
        }
        if(isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
            mdl <- x
        }
    }
    
    ##stopifnot(requireNamespace("R2jags"))  # needs to load this library (which automatically also loads R2WinBUGS)
    ##cmd <- ifelse(class(x)=="rjags",mdl <- x$BUGSoutput, mdl <- x)
    
    # Defines the inputs matrix  
    inputs <- mdl$sims.matrix
    # If the deviance is computed, then removes it 
    if("deviance"%in%colnames(inputs)) {
        w <- which(colnames(inputs)=="deviance")
        inputs <- inputs[,-w]
    }
    pars <- colnames(data.frame(inputs))
    list(mat=data.frame(inputs),parameters=pars)
}


######Structural PSA#############################################################################
struct.psa <- function(models,effect,cost,ref=1,interventions=NULL,Kmax=50000,plot=F) {
    # Computes the weights to be associated with a set of competing models in order to
    #   perform structural PSA
    # model is a list containing the output from either R2jags or R2OpenBUGS/R2WinBUGS
    #   for all the models that need to be combined in the model average
    # effect is a list containing the measure of effectiveness computed from the 
    #   various models (one matrix with n.sim x n.ints simulations for each model)
    # cost is a list containing the measure of costs computed from the 
    #   various models (one matrix with n.sim x n.ints simulations for each model)
    
    n.models <- length(models)  # number of models to be combined
    if(n.models==1) {
        stop("NB: Needs at least two models to run structural PSA")
    }
    d <- w <- numeric()		# initialises the relevant vectors
    mdl <- list()		     	# and list
    for (i in 1:n.models) {
        # 1. checks whether each model has been run using JAGS or BUGS    
        if(class(models[[i]])=="rjags") {
            if(!isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
                stop("You need to install the package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
            }
            if (isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
                mdl[[i]] <- models[[i]]$BUGSoutput  # if model is run using R2jags/rjags        
            }
        }
        if(class(models[[i]])=="bugs") {		# if model is run using R2WinBUGS/R2OpenBUGS
            if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
                stop("You need to install the package 'R2OpenBUGS'. Please run in your R terminal:\n install.packages('R2OpenBUGS')")
            }
            if(isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
                mdl[[i]] <- models[[i]]
            }
        }
        mdl[[i]] <- models[[i]]
        # 2. saves the DIC in the vector d
        d[i] <- mdl[[i]]$DIC
    }
    dmin <- min(d)					# computes the minimum value to re-scale the DICs
    w <- exp(-.5*(d-dmin))/sum(exp(-.5*(d-dmin))) 	# Computes the model weights (cfr BMHE)
    
    # Now weights the simulations for the variables of effectiveness and costs in each model
    # using the respective weights, to produce the economic analysis for the average model
    e <- c <- matrix(NA,dim(effect[[1]])[1],dim(effect[[1]])[2])
    e <- w[1]*effect[[1]]
    c <- w[1]*cost[[1]]
    for (i in 2:n.models) {
        e <- e + w[i]*effect[[i]]
        c <- c + w[i]*cost[[i]]
    }
    
    # Now performs the economic analysis on the averaged model
    he <- bcea(e=e,c=c,ref=ref,interventions=interventions,Kmax=Kmax,plot=plot)
    
    # And finally saves the results
    list(he=he,w=w,DIC=d)
}


###cost-effectiveness efficiency frontier#############################################
# tan(theta)=e/c
# theta=atan(e/c)
# if theta_1<theta_2, take 1

ceef.plot <- function(he, comparators=NULL, pos=c(1,1), start.from.origins=TRUE, threshold=NULL, flip=FALSE, dominance=TRUE, relative=FALSE, print.summary=TRUE, graph=c("base", "ggplot2"), ...) {
    ### plots the cost-effectiveness efficiency frontier, together with the scatter plot of the simulations and optionally the dominance areas
    if(is.null(he$c) | is.null(he$e)) stop("Please use the bcea() function from BCEA version >=2.1-0 or attach the vectors e and c to the bcea object. Please see ?ceef.plot for additional details.")
    
    ### if threshold is NULL, then bound to pi/2, which is atan(Inf); else iff positive, bound to the increase angle given the slope
    if(is.null(threshold))
        threshold  <- base::pi/2
    else{
        if(threshold<=0){
            warning("The value of the cost-effectiveness threshold should be positive. The argument will be ignored.")
            threshold  <- base::pi/2
        }
        else
            threshold <- atan(threshold)
    }
    
    ### selects the comparators. No need for recursion
    if(!is.null(comparators)){
        stopifnot(all(comparators %in% 1:he$n.comparators))
        # adjusts bcea object for the correct number of dimensions and comparators
        he$comp <- he$comp[comparators]
        he$n.comparators=length(comparators)
        he$n.comparisons=length(comparators)-1
        he$interventions=he$interventions[comparators]
        he$ref=rank(c(he$ref,he$comp))[1]
        he$comp=rank(c(he$ref,he$comp))[-1]
        he$mod <- TRUE #
        ### bceanew
        he$e <- he$e[,comparators]
        he$c <- he$c[,comparators]
    }
    
    ### If the incremental analysis (relative to the reference) is required, needs to modify the BCEA object
    if(relative) {
        temp <- he
        temp$e <- temp$c <- matrix(NA,he$n.sim,he$n.comparators)
        temp$e[,he$ref] <- temp$c[,he$ref] <- rep(0,he$n.sim)
        temp$e[,-he$ref] <- -he$delta.e
        temp$c[,-he$ref] <- -he$delta.c
        he <- temp
    }
    
    stopifnot(he$n.comparators>=2)
    base.graphics <- ifelse(isTRUE(pmatch(graph,c("base","ggplot2"))==2),FALSE,TRUE)
    
    ### no visible binding note
    c.avg <- e.avg <- x <- y <- e <- e.orig <- c.orig <- NA_real_
    
    ### if the effectiveness is negative or !start.from.origins, rescale
    ec.min <- with(he,c(min(apply(e,2,mean)),
                        apply(c,2,mean)[which.min(apply(e,2,mean))],
                        which.min(apply(e,2,mean))))
    e.neg <- ec.min[1]<0
    c.neg <- any(apply(he$c,2,mean)<0)
    
    if(e.neg & !c.neg & start.from.origins){
        message("Benefits are negative, the frontier will not start from the origins")
        start.from.origins <- FALSE
    }
    if(!e.neg & c.neg & start.from.origins){
        message("Costs are negative, the frontier will not start from the origins")
        start.from.origins <- FALSE
    }
    if(e.neg & c.neg & start.from.origins){
        message("Costs and benefits are negative, the frontier will not start from the origins")
        start.from.origins <- FALSE
    }
    e.neg <- ifelse(start.from.origins,e.neg,TRUE)
    
    ### frontier calculation
    data.avg <- data.frame(
        "e.avg"=apply(he$e,2,mean)-ifelse(!e.neg,0,ec.min[1]),
        "c.avg"=apply(he$c,2,mean)-ifelse(!e.neg,0,ec.min[2]))
    data.avg <- cbind(data.avg,data.avg,as.factor(c(1:dim(data.avg)[1])))
    names(data.avg)[3:5] <- c("e.orig","c.orig","comp")
    orig.avg <- data.avg[,3:5]
    ### check for interventions with zero costs and effectiveness
    comp <- ifelse(any(apply(data.avg[,1:2],1,function(x) isTRUE(sum(x) == 0 & prod(x) == 0))),
                   which(apply(data.avg[,1:2],1,sum)==0 & apply(data.avg[,1:2],1,prod)==0),0)
    ### contains the points connecting the frontier. Always starts from the origins
    ceef.points <- data.frame(
        "x"=0,
        "y"=0,
        "comp"=comp)
    repeat{
        if(prod(dim(data.avg))==0) break
        theta <- with(data.avg,atan(c.avg/e.avg))
        theta.min <- min(theta,na.rm=TRUE)
        if(theta.min>threshold) break
        index <- which(theta==theta.min)
        if(length(index)>1)
            index=index[which.min(data.avg$e.avg[index])]
        ceef.points <- with(data.avg,rbind(ceef.points,c(e.orig[index],c.orig[index],comp[index])))
        data.avg[,1:2] <- data.avg[,3:4]-matrix(rep(as.numeric(data.avg[index,3:4]),dim(data.avg)[1]),ncol=2,byrow=TRUE)
        data.avg <- subset(subset(data.avg,c.avg*e.avg>0),c.avg+e.avg>0)
    }
    ceef.points$comp <- factor(ceef.points$comp)
    
    ceef.points$slope <- NA
    ### calculate slopes
    for(i in 2:dim(ceef.points)[1])
        ceef.points$slope[i] <- with(ceef.points,(y[i]-y[i-1])/(x[i]-x[i-1]))
    
    ### workaround for start.from.origins == FALSE: remove first row if slope is negative
    while(dim(ceef.points)[1]>1 & ceef.points$slope[2]<0){
        ceef.points <- ceef.points[-1,]
        ceef.points$slope[1] <- NA
    }
    
    ### set data.frame for points
    scatter.data <- data.frame(
        "e"=c(he$e),#-ifelse(!e.neg,0,ec.min[1]),
        "c"=c(he$c),#-ifelse(!e.neg,0,ec.min[2]),
        "comp"=as.factor(sort(rep(1:he$n.comparators,he$n.sim))))
    
    ### re-adjustment of data sets
    ceef.points[,1] <- ceef.points[,1]+ifelse(!e.neg,0,ec.min[1])
    ceef.points[,2] <- ceef.points[,2]+ifelse(!e.neg,0,ec.min[2])
    orig.avg[,1] <- orig.avg[,1]+ifelse(!e.neg,0,ec.min[1])
    orig.avg[,2] <- orig.avg[,2]+ifelse(!e.neg,0,ec.min[2])
    
    ### Summary table function
    ceef.summary <- function(he,ceef.points,orig.avg,include.ICER=FALSE,...){
        ## Tables adaptation and formatting
        no.ceef <- which(!1:he$n.comparators %in% ceef.points$comp)
        ## Interventions included
        if(ceef.points$comp[1]==0)
            ceef.points <- ceef.points[-1,]
        rownames(ceef.points) <- he$interventions[as.numeric(levels(ceef.points$comp)[ceef.points$comp])]
        
        if(!include.ICER){
            ceef.points[,5] <- atan(ceef.points[,4]^(1*ifelse(!flip,1,-1)))
            ceef.points <- ceef.points[,-3]
            colnames(ceef.points) <- c("Effectiveness","Costs","Increase slope","Increase angle")
        }
        else{
            ICERs <- numeric(dim(ceef.points)[1])
            index <- as.numeric(levels(ceef.points$comp)[ceef.points$comp])
            for(i in 1:length(ICERs)){
                if(ceef.points$comp[i]==he$ref)
                    ICERs[i] <- NA_real_
                else
                    ICERs[i] <- he$ICER[index[i]+ifelse(index[i]<he$ref,0,-1)]
            }
            ceef.points[,3] <- ICERs
            ceef.points[,5] <- atan(ceef.points[,4]^(1*ifelse(!flip,1,-1)))
            colnames(ceef.points) <- c("Effectiveness","Costs",paste0("ICER ",he$interventions[he$ref]," vs."),"Increase slope","Increase angle")
        }
        if(flip) colnames(ceef.points)[1:2] <- colnames(ceef.points[2:1])
        
        ## Interventions not included
        if(length(no.ceef)>0){
            noceef.points <- data.frame(matrix(NA_real_,ncol=4,nrow=length(no.ceef)))
            noceef.points[,1:2] <- orig.avg[no.ceef,-3]
            
            if(!include.ICER){
                noceef.points <- noceef.points[,-3]
                colnames(noceef.points) <- c("Effectiveness","Costs","Dominance type")
            }
            else{
                ICERs <- numeric(dim(noceef.points)[1])
                for(i in 1:length(ICERs)){
                    if(no.ceef[i]==he$ref)
                        ICERs[i] <- NA_real_
                    else
                        ICERs[i] <- he$ICER[no.ceef[i]+ifelse(no.ceef[i]<he$ref,0,-1)]
                }
                noceef.points[,3] <- ICERs
                colnames(noceef.points) <- c("Effectiveness","Costs",paste0("ICER ",he$interventions[he$ref]," vs."),"Dominance type")
            }
            
            how.dominated <- rep("Extended dominance",length(no.ceef))
            for(i in 1:length(no.ceef))
                for(j in 1:dim(ceef.points)[1]){
                    ### if the product of the deltas is negative it is dominated: cannot be dominant since not on the frontier 
                    if((noceef.points[i,1]-ceef.points[j,1])*(noceef.points[i,2]-ceef.points[j,2])<0){
                        how.dominated[i] <- "Absolute dominance"
                        ### alternative:
                        # how.dominated[i] <- paste0("Dominated by ",rownames(ceef.points)[j])
                        break
                    }
                }
            noceef.points[,ifelse(!include.ICER,3,4)] <- how.dominated
            rownames(noceef.points) <- he$interventions[no.ceef]
            if(flip) colnames(noceef.points)[1:2] <- colnames(noceef.points)[2:1]
        }
        
        ### Print the summary table
        cat("\nCost-effectiveness efficiency frontier summary \n\n")
        cat("Interventions on the efficiency frontier:\n")
        print(ceef.points,quote=F,digits=5,justify="center")
        cat("\n")
        if(length(no.ceef)>0){
            cat("Interventions not on the efficiency frontier:\n")
            print(noceef.points,quote=F,digits=5,justify="center")
        }
    }
    
    ### colours
    colour <- colours()[floor(seq(262,340,length.out=he$n.comparators))]  # gray scale
    
    ### plots
    ##### ***** base graphics ***** #####
    if(base.graphics){
        ### legend positioning
        if(is.numeric(pos)&length(pos)==2){
            temp <- ""
            if(pos[2]==0)
                temp <- paste0(temp,"bottom")
            else
                temp <- paste0(temp,"top")
            if(pos[1]==0)
                temp <- paste0(temp,"left")
            else
                temp <- paste0(temp,"right")
            pos <- temp
            if(length(grep("^(bottom|top)(left|right)$",temp))==0)
                pos <- FALSE
        }
        if(is.logical(pos)){
            if(!pos)
                pos="topright"
            else
                pos="topleft"
        }
        
        if(flip){
            temp <- scatter.data$e
            scatter.data$e <- scatter.data$c
            scatter.data$c <- temp
            
            temp <- ceef.points$x
            ceef.points$x <- ceef.points$y
            ceef.points$y <- temp
            
            temp <- orig.avg$e.orig
            orig.avg$e.orig <- orig.avg$c.orig
            orig.avg$c.orig <- temp
            
            rm(temp)
        }
        
        ### set up plot window
        xlab=ifelse((!flip & !relative),"Effectiveness",
                    ifelse((!flip & relative),"Differential effectiveness",
                           ifelse((flip & !relative),"Cost","Differential cost")))
        ylab=ifelse((!flip & !relative),"Cost",
                    ifelse((!flip & relative),"Differential cost",
                           ifelse((flip & !relative),"Effectiveness","Differential effectiveness")))
        plot(NULL,
             xlim=c(min(range(scatter.data$e)[1],0),max(range(scatter.data$e)[2],0)),
             ylim=c(min(range(scatter.data$c)[1],0),max(range(scatter.data$c)[2],0)),
             main="Cost-effectiveness efficiency frontier",
             xlab=xlab,ylab=ylab)
        
        if(dominance){
            ### add dominance regions
            for(i in 1:dim(ceef.points)[1]){
                rect(col="grey95",border=NA,
                     xleft=ifelse(!flip,-1,1)*2*max(abs(range(scatter.data$e))),xright=ceef.points$x[i],
                     ybottom=ceef.points$y[i],ytop=ifelse(!flip,1,-1)*2*max(abs(range(scatter.data$c))))
            }
            if(dim(ceef.points)[1]>1)
                for(i in 2:dim(ceef.points)[1]){
                    rect(col="grey85",border=NA,
                         xleft=ifelse(!flip,-1,1)*2*max(abs(range(scatter.data$e))),xright=ceef.points$x[ifelse(!flip,i-1,i)],
                         ybottom=ceef.points$y[ifelse(!flip,i,i-1)],ytop=ifelse(!flip,1,-1)*2*max(abs(range(scatter.data$c))))
                }
        }
        
        ### plot the axes
        abline(h=0,col="grey")
        abline(v=0,col="grey")
        
        ### plot the scatter
        for(i in 1:he$n.comparators)
            with(scatter.data,points(subset(scatter.data,comp==i)[,-3],type="p",pch=20,cex=.35,col=colour[i]))
        
        ### plot the frontier
        points(ceef.points[,1:2],type="l",lwd=2)
        ### add circles
        points(orig.avg[,-3],pch=21,cex=2,bg="white",col="black")
        ### add text; grey if not on the frontier
        for(i in 1:he$n.comparators){
            text(orig.avg[i,-3],labels=orig.avg[i,3],col=ifelse(i %in% ceef.points$comp,"black","grey60"),cex=.75)
        }
        ### legend text
        text <- paste(1:he$n.comparators,":",he$interventions)
        legend(pos,text,col=colour,cex=.7,bty="n",lty=1)
        
        ### needed for dominance areas overwriting the outer box
        box()
        
        if(print.summary)
            ceef.summary(he,ceef.points,orig.avg,...)
    }
    ##### ***** ggplot2 ***** #####
    else{
        if(!isTRUE(requireNamespace("ggplot2",quietly=TRUE)&requireNamespace("grid",quietly=TRUE))){
            message("Falling back to base graphics\n")
            ceef.plot(he,flip=flip,comparators=comparators,pos=pos,start.from.origins=start.from.origins,graph="base")
            return(invisible(NULL))
        }
        
        opt.theme <- ggplot2::theme()
        exArgs <- list(...)
        if(length(exArgs)>=1){
            for(obj in exArgs)
                if(ggplot2::is.theme(obj))
                    opt.theme <- opt.theme + obj
        }
        
        ceplane <- ggplot2::ggplot(ceef.points,ggplot2::aes(x=x,y=y))
        
        if(dominance){
            ### add dominance regions
            ceplane <- ceplane +
                ggplot2::geom_rect(data=ceef.points,ggplot2::aes(xmax=x,ymin=y),
                                   ymax=2*max(abs(range(scatter.data$c))),xmin=-2*max(abs(range(scatter.data$e))),
                                   alpha=.35,fill="grey75")
        }
        
        ceplane <- ceplane +
            ### draw axes
            ggplot2::geom_hline(yintercept=0,colour="grey")+ggplot2::geom_vline(xintercept=0,colour="grey")+
            ### add scatter points
            ggplot2::geom_point(data=scatter.data,ggplot2::aes(x=e,y=c,colour=comp),size=1)
        ### add frontier
        if(dim(ceef.points)[1]>1)
            ceplane <- ceplane + ggplot2::geom_path()
        ### add circles
        xlab=ifelse(!relative,"Effectiveness","Effectiveness differential")
        ylab=ifelse(!relative,"Cost","Cost differential")
        ceplane <- ceplane +
            ggplot2::geom_point(data=orig.avg,ggplot2::aes(x=e.orig,y=c.orig),size=5.5,colour="black")+
            ggplot2::geom_point(data=orig.avg,ggplot2::aes(x=e.orig,y=c.orig),size=4.5,colour="white")+
            ### set graphical parameters
            ggplot2::scale_colour_manual("",labels=paste0(1:he$n.comparators,": ",he$interventions),values=colour,na.value="black")+
            ggplot2::labs(title="Cost-effectiveness efficiency frontier",x=xlab,y=ylab)+
            ggplot2::theme_bw()
        ### add text into circles
        for(i in 1:he$n.comparators){
            ceplane <- ceplane + 
                ggplot2::geom_text(data=orig.avg[i,],ggplot2::aes(x=e.orig,y=c.orig,label=comp),size=3.5,
                                   colour=ifelse(i %in% ceef.points$comp, "black", "grey60"))
        }
        
        jus <- NULL
        if(isTRUE(pos)) {
            pos="bottom"
            ceplane <- ceplane + ggplot2::theme(legend.direction="vertical")
        }
        else{
            if(is.character(pos)) {
                choices <- c("left", "right", "bottom", "top")
                pos <- choices[pmatch(pos,choices)]
                jus="center"
                if(is.na(pos))
                    pos=FALSE
            }
            if(length(pos)>1)
                jus <- pos
            if(length(pos)==1 & !is.character(pos)) {
                pos <- c(1,1)
                jus <- pos
            }
        }
        
        ceplane <- ceplane + 
            ggplot2::theme(legend.position=pos,legend.justification=jus,legend.title=ggplot2::element_blank(),
                           legend.background=ggplot2::element_blank(),text=ggplot2::element_text(size=11),
                           legend.key.size=grid::unit(.66,"lines"),legend.margin=grid::unit(-1.25,"line"),
                           panel.grid=ggplot2::element_blank(),legend.key=ggplot2::element_blank(),legend.text.align=0,
                           plot.title = ggplot2::element_text(face="bold",lineheight=1.05,size=14.3)) +
            opt.theme
        
        if(flip) ceplane  <- ceplane + ggplot2::coord_flip()
        if(print.summary) ceef.summary(he,ceef.points,orig.avg,...)
        return(ceplane)
    }
}


# ###info.rank#########################################################################
 info.rank <- function(parameter,input,he,wtp=he$k[min(which(he$k>=he$ICER))],...) {
   # parameter = vector of parameters for which to make the plot
   # input = a matrix of PSA runs for all the parameters
   # he = a bcea object with the economic evaluation
   # wtp = a willingness-to-pay threshold (default at the break even point from he)
   # ... extra arguments including
   #     xlim = x-axis limits
   #     ca = cex axis label (default = 0.7)
   #     cn = cex names label (default = 0.7)
   #     mai = graphical parameter to determine the margins
   #     rel = if TRUE (default) then shows a plot of EVPPI/EVPI, if FALSE then only EVPPI
   #     N = number of PSA to be used to perform the evppi analysis
   
   # Function to actually create the bar plot
   make.barplot <- function(scores,chk2,xlab,xlim,ca,cn,mai,space) {
     col <- rep(c("blue","red"),length(chk2))
     par(mai=mai)
     res <- data.frame(parameter=names(chk2),info=scores,row.names=NULL)
     barplot(res[order(res$info),2],horiz=T,names.arg=res[order(res$info),1],cex.names=cn,las=1,col=col,cex.axis=ca,
             xlab=xlab,space=space,main=tit,xlim=xlim)
     par(mai=c(1.360000,1.093333,1.093333,0.560000))
     list(rank=data.frame(parameter=res[order(-res$info),1],info=res[order(-res$info),2]))
   }
   
   # Prevents BCEA::evppi from throwing messages
   quiet <- function(x) { 
     sink(tempfile()) 
     on.exit(sink()) 
     invisible(force(x)) 
   } 
   
   exArgs <- list(...)
   
   if(class(parameter[1])=="character"){
     parameters<-array()
     for(i in 1:length(parameter)){
       parameters[i]<-which(colnames(input)==parameter[i])
     }
   }else{
     parameters=parameter
     parameter=colnames(input)[parameters]
   }
   
   # needs to exclude parameters with weird behaviour (ie all 0s)
   w <- unlist(lapply(parameter,function(x) which(colnames(input)==x)))
   input <- input[,w]
   chk1 <- which(apply(input,2,var)>0)   # only takes those with var>0
   tmp <- list(); tmp <- apply(input,2,function(x) table(x)) # check those with <5 possible values (would break GAM)
   chk2 <- which(unlist(lapply(tmp,function(x) length(x)>=5))==TRUE)
   
   # Can do the analysis on a smaller number of PSA runs
   if(exists("N",where=exArgs)) {N <- exArgs$N} else {N <- he$n.sim}
   if(any(!is.na(N)) & length(N)>1) {
     select <- N
   } else {
     N <- min(he$n.sim,N,na.rm=T)
     if(N==he$n.sim) {select <-1:he$n.sim} else {select <- sample(1:he$n.sim,size=N,replace=F)} 
   }
   m <- he; m$k=wtp
   x <- list()
   for (i in 1:length(chk2)) {
     x[[i]] <- quiet(evppi(parameter=chk2[i],input=input,he=m,N=N))
   }
   scores <- unlist(lapply(x,function(x) x$evppi/x$evi[which(he$k==wtp)]))
   # Optional inputs
   if(exists("ca",where=exArgs)) {ca <- exArgs$ca} else {ca <- .7}
   if(exists("cn",where=exArgs)) {cn <- exArgs$cn} else {cn <- .7}
   xlab <- "Proportion of total EVPI"
   if(exists("rel",where=exArgs)) {
     if (exArgs$rel==FALSE) {
       scores <- unlist(lapply(x,function(x) x$evppi))
       xlab <- "Absolute value of the EVPPI"
     }
   }
   if(exists("xlim",where=exArgs)) {xlim=exArgs$xlim} else {xlim=c(0,range(scores)[2])}
   if(exists("mai",where=exArgs)) {mai=exArgs$mai} else {mai=c(1.36,1.5,1,1)}
   if(exists("tit",where=exArgs)) {tit=exArgs$tit} else {tit <- paste0("Info-rank plot for willingness to pay=",wtp)}
   if(exists("space",where=exArgs)) {space=exArgs$space} else {space=.5}
   
   # Makes the plot
   make.barplot(scores=scores,chk2=chk2,xlab=xlab,xlim=xlim,ca=ca,cn=cn,mai=mai,space)
 }

