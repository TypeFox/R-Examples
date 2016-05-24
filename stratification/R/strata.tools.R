`print.strata` <-
function(x,...)
{
    takenone <- if (is.null(x$args$takenone)) 0 else x$args$takenone
    L <- x$args$Ls + takenone
    rh <- x$args$rh 
    ph <- c(x$args$model.control$ptakenone,x$args$model.control$ph,x$args$model.control$pcertain)
    
    # Section des arguments fournis
    cat("Given arguments:\n")
    cat("x = "); print(x$call$x)
    if (!is.null(x$args$nclass)) cat("nclass = ",x$args$nclass,", ",sep="")
    if (!is.null(x$args$CV)) cat("CV = ",x$args$CV,", ",sep="")
    if (!is.null(x$args$n)) cat("n = ",x$args$n,", ",sep="")
    cat("Ls = ",x$args$Ls,sep="")
    if (!is.null(x$args$takenone)) cat(", takenone = ",x$args$takenone,sep="")
    if (!is.null(x$args$bias.penalty)) { if (x$args$takenone==1) cat(", bias.penalty = ",x$args$bias.penalty,sep="")}
    if (!is.null(x$args$takeall)) cat(", takeall = ",x$args$takeall,sep="")
    cat("\nallocation: q1 = ",x$args$alloc$q1,", q2 = ",x$args$alloc$q2,", q3 = ",x$args$alloc$q3,sep="")
    if (!is.null(x$args$model)) {
        cat("\nmodel = ",x$args$model,sep="")
        nparam <- length(x$args$model.control)
        if (nparam>0) {
            cat(": ")
            for (i in 1:nparam) {
                    cat(names(x$args$model.control)[i],"=",x$args$model.control[[i]],sep=" ")
                    if (i<nparam) cat(", ")
            }
        }
    }
    if (!is.null(x$args$algo)) {
             cat("\n")
             cat("algo = ", x$args$algo, ": ", sep = "")
             for (i in 1:length(x$args$algo.control)) {
                 if (i %in% c(5,9)) cat("\n              ")
                 cat(names(x$args$algo.control)[i]," = ",x$args$algo.control[[i]],sep="")
                 if (i<length(x$args$algo.control)) cat(", ")
             }
    }
    
    # Section du tableau de stratification
    tableau <- data.frame(x$meanh,x$varh,x$Nh,x$nh,ifelse(x$Nh==0,NA,x$nh/x$Nh))
    colnames(tableau) <- c("E(Y)","Var(Y)","Nh","nh","fh")
    if(!is.null(x$args$certain))
        tableau <- rbind(tableau,c(x$certain.info["meanc"],NA,x$certain.info["Nc"],x$certain.info["Nc"],1))
    tableau <- rbind(tableau,c(NA,NA,round(sum(tableau$Nh)),round(x$n),x$n/sum(tableau$Nh)))
    rownames(tableau)<-if(!is.null(x$args$certain)) c(paste("stratum",1:L),"","Total") else c(paste("stratum",1:L),"Total")
    if (grepl("strata.bh", deparse(x$call[[1]]))) {
        tableau<-cbind("bh"=c(x$args$bh,max(x$args$x)+1,rep(NA,nrow(tableau)-L)),
                       "|"=c(rep("|",nrow(tableau)-1),NA),
                       tableau)
    } else if (identical(as.character(x$call[[1]]),"strata.cumrootf")||identical(as.character(x$call[[1]]),"strata.geo")) {
        tableau<-cbind("|"=c(rep("|",nrow(tableau)-1),NA),
                       "bh"=c(x$bh,max(x$args$x)+1,rep(NA,nrow(tableau)-L)),
                       tableau)
    } else {
        tableau<-cbind("|"=c(rep("|",nrow(tableau)-1),NA),
                       "bh"=c(x$bh,max(x$args$x)+1,rep(NA,nrow(tableau)-L)),
                       tableau)
        if(is.numeric(x$args$initbh)) 
          tableau <- cbind("initbh" = c(x$args$initbh, max(x$args$x)+1, rep(NA,nrow(tableau)-L)), tableau)
    }
    rh.tab <- if(is.null(x$args$certain)) c(rh,NA) else c(rh,1,NA)
    if(takenone>0) rh.tab <- c(rep(NA, length(takenone)), rh.tab)
    tableau <- cbind("rh" = rh.tab,tableau)
    if (identical(x$args$model,"loglinear"))  tableau<-cbind("ph"=c(ph,NA),tableau)
    type <- c(rep("take-none", takenone), rep("take-some", x$args$Ls - x$takeall), rep("take-all", x$takeall))
    tableau <- cbind("|" = c(rep("|", nrow(tableau) - 1), NA),
                     "type" = if(is.null(x$args$certain)) c(type, NA) else c(type, "certain", NA), 
                     tableau)
    tableau[,unlist(lapply(tableau,is.numeric))]<-round(tableau[,unlist(lapply(tableau,is.numeric))],2)
#    tableau[,7:8]<-round(tableau[,ncol(tableau)-2:1],0)
    ### Correction pour affichier correctement les NA
    tableauc <- format(tableau)
    substr2last <- function(x) substr(x, nchar(x) -1 , nchar(x))
    for (i in 1:(nrow(tableauc)-1))
      tableauc[i,] <- ifelse(substr2last(tableauc[i,]) %in% c("NA", "aN"), "-", tableauc[i,])
    tableauc[dim(tableauc)[1],] <- ifelse(substr2last(tableauc[dim(tableauc)[1],]) == "NA", "", 
                                   ifelse(substr2last(tableauc[dim(tableauc)[1],]) == "aN", "-", tableauc[dim(tableauc)[1],]))
    ### Fin de la correction
    cat("\n\nStrata information:\n")
    print(tableauc,na.print="")
    cat("\nTotal sample size:",x$n,"\n")
    cat("Anticipated population mean:",x$mean,"\n")
    
    # Section sur les moments anticipés
    sortie <- if (is.null(x$args$takenone)) 1 else { if(0==x$args$takenone) 2 else 3 }
    if (sortie%in%c(1,2)) {
        cat("Anticipated CV:",ifelse(1==sortie,x$CV,x$RRMSE),"\n")
        if (2==sortie) cat("Note: CV=RRMSE (Relative Root Mean Squared Error) because takenone=0.\n")
    } else {
        est<-cbind(x$relativebias,x$propbiasMSE,x$RRMSE,x$args$CV)
        dimnames(est) <- if(length(est)==4)  list(" ",c("relative bias","prop bias MSE","RRMSE","target RRMSE"))
                         else list(" ",c("relative bias","prop bias MSE","RRMSE"))
        cat("\nAnticipated moments of the estimator:\n")
        print.default(est, print.gap = 2, quote = FALSE, right=TRUE)
    }

    if (!is.null(x$converge) && !is.na(x$converge)) { if (!x$converge) cat("\nWarning : the algorithm did not converge.\n") }
}


`plot.strata` <-
function(x,logscale=FALSE,drop=0,main=paste("Graphical Representation of the Stratified Design",xname),xlab,...)
{
    L <- x$args$L
    ncert <- if (is.null(x$args$certain)) 0 else 1
    if (!((length(drop)==1)&&isTRUE((drop%%1)==0)&&(isTRUE(drop>=0)||isTRUE(drop<x$Nh[L]))))
        stop("'drop' must be an integer between 0 and 'Nh[L]'-1 inclusively")
    data <- if(is.null(x$args$certain)) sort(x$args$x) else sort(x$args$x[-x$args$certain])     
    data <- data[1:(length(data)-drop)] # pour enlever les drop données les plus grandes
    if(logscale) data <- log(data)
    bh <- if (identical(as.character(x$call[[1]]),"strata.bh")) x$args$bh else x$bh
    bhfull <- if(logscale) c(min(data),log(bh),max(data)+1) else c(min(data),bh,max(data)+1)
    if (missing(xlab)) xlab <- if(logscale) "Stratification variable X on the log scale" else "Stratification variable X"
    xname <- paste(deparse(substitute(x), 500), collapse = "\n")
    
    if (L<10)
    {
    layout(matrix(c(2,1),2,1),heights=c(1,4))
    
    # Histogramme
    par(mar=c(5.1, 4.1, 0.1, 4.1))
    hist(data,breaks="FD",col=rgb(0.5,0.7,1),border=rgb(0.3,0.5,1),freq=FALSE,main="",xlab=xlab,xlim=c(bhfull[1],bhfull[L+1]))
    for (i in 2:L) abline(v=bhfull[i],lty=1)
    
    # Tableau en haut de l'histogramme
    par(mar=c(1, 4.1, 3.6, 4.1))
    plot(x=c(bhfull[1],bhfull[L+1]),y=rep(2,2),xlim=c(bhfull[1],bhfull[L+1]),ylim=c(1,3),
         type="l",xaxt="n",yaxt="n",xlab="",ylab="")
    abline(h=2)
    space <- bhfull[-1]-bhfull[-(L+1)]
    lcert <- (bhfull[L+1]-bhfull[1])/10 # largeur pour la strate certaine
    off <- (bhfull[L+1]-bhfull[1])*0.04 # offset par défaut pour les axes dans les graphiques en R 
    if (!is.null(x$args$certain)) space <- c(space[-L],space[L]-lcert,lcert)
    if (any(space<(bhfull[L+1]-bhfull[1])/15)) {
        space <- rep(bhfull[L+1]-bhfull[1]+off*2,L+ncert)/(L+ncert)
        space[1] <- space[1] - off
        space[length(space)] <- space[length(space)] - off
    }    
    abline(v=bhfull[1]+cumsum(space)[-length(space)],lty=1)
    pos <- bhfull[1]+c(0,cumsum(space)[-length(space)])+space/2
    pos[1] <- pos[1] - off/2    
    pos[length(pos)] <- pos[length(pos)] + off/2    
    text(x=pos,y=rep(2.5,length(space)),labels=as.character(c(x$Nh,x$certain.info["Nc"]))[1:(L+ncert)])
    text(x=pos,y=rep(1.5,length(space)),labels=as.character(c(x$nh,x$certain.info["Nc"]))[1:(L+ncert)])
    mtext(c("Nh","nh"),at=c(2.5,1.5),line=1,side=2,las=1,font=2)
    mtext(c(sum(x$Nh),sum(x$nh))+x$certain.info["Nc"],at=c(2.5,1.5),line=1,side=4,las=1,font=2)
    mtext(c(1:L,"certain")[1:(L+ncert)],at=pos,line=0,side=3,las=1,font=2)
    
    # Titre
    mtext(main,side=3,las=1,font=2,line=2,cex=1.2)

    par(mar=c(5.1, 4.1, 4.1, 2.1))
    layout(matrix(1,1,1))

    } else {
    # Histogramme
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    hist(data,breaks="FD",col=rgb(0.5,0.7,1),border=rgb(0.3,0.5,1),freq=FALSE,main=main,xlab=xlab,xlim=c(bhfull[1],bhfull[L+1]))
    for (i in 2:L) abline(v=bhfull[i],lty=1)
    }
    
}
