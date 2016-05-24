#####
##
## PCA plot functions
##
#####

plotVariable <- function(acp, axes=c(1,2), new.plot=FALSE, lab=NULL, lim.cos2.var=0, palette="rainbow", ...) {
    if (!is.null(lab)){
        col.lab <- as.colors(lab, palette=palette)
    }
    else {
        col.lab <- rep(1,dim(acp$var$coord)[1])
    }
    
    plot.PCA(acp, choix="var", axes=axes, title=paste("Variable Representation on axes", axes[1], "and", axes[2]) , new.plot=new.plot, col.var=col.lab, lim.cos2.var=lim.cos2.var, cex=0.8, ...)
    
    if (!missing(lab)) {
        legend("topright", legend = unique(lab), col=unique(col.lab), pch=20, bty="n", cex=0.7, text.col="gray50")
    }
}

qualitySample <- function(acp, axes = c(1:3)){
  qual <- data.frame(acp$ind$cos2)[, axes]
  name <- colnames(qual)
  ncomp <- length(axes)
  k <- ncomp+1

  for (i in 1:(ncomp-1)) {
    for (j in (i+1):ncomp) {
      qual[,k] <- qual[,i]+qual[,j]
      name <- c(name, paste(name[i], "-", name[j], sep=""))
      k <- k+1
    }
  }
  names(qual) <- name  
  return(qual)
}



plotSample <-  function(acp, axes=c(1,2), new.plot=FALSE, lab="quality", palette="rainbow", lim.cos2.sample=0, text=TRUE, lab.title=NULL, ellipse=FALSE, ...) {
    nind <- dim(acp$ind$coord)[1]
    col.PC <- rep("black", nind)
    qual <- qualitySample(acp, axes=axes)[,3]
    
    if (!is.null(lab)) {
        if ((length(lab) == 1)&&(lab == "quality")) {
            col.PC <- rep("grey", nind)
            col.PC[qual < 0.3] <- "red"
            col.PC[qual > 0.7] <- "green"
        }
        else{
            col.PC <- as.colors(lab, palette=palette)
            corresColLab <- unique(col.PC)
            names(corresColLab) <- unique(lab)
        }
    }
    
    col.PC[qual < lim.cos2.sample] <- "white" 
    if (text){
        text="all"
    }else{
        text=NULL
    }
    plot.PCA(acp, axes=c(axes[1],axes[2]), col.ind=col.PC, title=paste("Sample representation \n Axes", axes[1], "and", axes[2]), new.plot=new.plot, label=text, cex=0.8, ...)
    
    ## ellipses have to be plotted by hand
    if ((ellipse) && !is.null(lab) && (lab != "quality")){
        aux <- cbind.data.frame(lab, acp$ind$coord[, axes])
        coord.ell <- coord.ellipse(aux, bary=TRUE)$res
        nbre.ellipse <- nlevels(coord.ell[, 1])
        for (e in 1:nbre.ellipse) {
            data.elli <- coord.ell[coord.ell[, 1] == levels(coord.ell[, 1])[e], -1]
            lines(x = data.elli[, 1], y = data.elli[, 2], col = corresColLab[levels(coord.ell[, 1])[e]])
        }
    }
    
    
    if (!is.null(lab)) {
        if ((length(lab) == 1)&&(lab == "quality")) {
            legend("topright", paste("Quality index", c("<0.3",">=0.3 and <=0.7",">0.7")), col=c("red","grey","green"), pch=20, title="Quality", bty="n", cex=0.7, text.col="gray50")
        } else {
            legend("topright", legend=unique(lab), col=unique(col.PC), pch=20, title=lab.title, bty="n", cex=0.7, text.col="gray50")
        }
    }
}

plotBiplot <- function(acp, ...) {
  v<-acp$eig[1:2,2]
  cp1 <- round(v[1], digits = 2)
  cp2 <- round(v[2], digits = 2)
  lab.x <- paste("Dimension ",1," (",cp1,"%)",sep="")
  lab.y <- paste("Dimension ",2," (",cp2,"%)",sep="")
  nind<-dim(acp$ind$coord)[1]
 
  indcoord<-t(t(acp$ind$coord[,1:2])/sqrt(acp$eig)[1:2,1])/sqrt(nind)
  x <- indcoord
  y <- acp$var$coord[,1:2]*(sqrt(nind))
  biplot(x,y, ...)

  unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
  rangx1 <- unsigned.range(x[, 1])
  rangx2 <- unsigned.range(x[, 2])
  rangx1 <- rangx2 <- range(rangx1, rangx2)
  rangy1 <- unsigned.range(y[, 1])
  rangy2 <- unsigned.range(y[, 2])
  ratio <- max(rangy1/rangx1, rangy2/rangx2)
  points(indcoord*ratio,pch=3)
  abline(h=0,lty=2)
  abline(v=0,lty=2)
}

plotInertia <- function(acp, ncp=5, ...) {
    barplot(acp$eig[1:ncp, 2], main = "Inertia percentage of components", las=2, names.arg = paste("Dim", 1:ncp, sep = ""), ...)
}



#####
##
## PCA functions
##
#####

runPCA <- function(X, ncp=5, scale=TRUE, ind.sup=NULL, quanti.sup=NULL, quali.sup=NULL, sample.qual=TRUE, variable.qual=FALSE, sample.cont=TRUE, variable.cont=FALSE, plotSample=TRUE, plotVariable=FALSE, plotInertia = TRUE, plotBiplot=FALSE, lab.sample="quality", lab.var=NULL, palette="rainbow", lim.cos2.sample=0, lim.cos2.var=0, pdf=FALSE, pdfname= NULL, verbose=FALSE, ...) {
    
    ## lot of options here ! try with ... ?
    ##  if (class(X) == "ExpressionSet") X <- exprs(X)
    if (!is.null(pdfname)) {
        pdf <- TRUE
    }
    
    if (pdf){
        if (!is.null(pdfname))
            pdf(pdfname, width=8, height=8, pointsize=8)
        else {
            pdf(width=8, height=8, pointsize=8)
        }
        op <- par(mfrow=c(2,1))
    }
    
    acp <- PCA(X, ncp=min(ncp,dim(X)[2]), scale.unit=scale, ind.sup=ind.sup, quanti.sup=quanti.sup, quali.sup=quali.sup, graph=FALSE)
    
    if (verbose) { 
        message("---------------- Sample coordinates on axes -----------------")
        print(acp$ind$coord)
    }
    
    if ((sample.qual) & (verbose)) {
        rep.sample <- qualitySample(acp, axes=c(1:min(ncp, dim(X)[2])))
        message("---------------- Sample Quality -----------------")
        print(rep.sample)
    }
    
    if ((variable.qual) & (verbose)) {
        message("---------------- Variable Quality -----------------")
        print(acp$var$cos2)
    }
    
    if ((sample.cont) & (verbose)) {
        message("---------------- Sample Contribution (%) -----------------")
        print(acp$ind$contrib)
    }
    
    if ((variable.cont) & (verbose)) {
        message("---------------- Variable Contribution (%) -----------------")
        print(acp$var$contrib)
    }
        
    if (plotInertia) {
        if (!pdf)
            dev.new()
        plotInertia(acp, ncp=ncp, ...)
    }
    
    if (plotBiplot) {
        plotBiplot(acp, ...)
    }
    
    if (plotSample) {
        plotSample(acp, axes=c(1,2), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
        plotSample(acp, axes=c(1,3), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
        plotSample(acp, axes=c(2,3), new.plot=!pdf, lab=lab.sample, lim.cos2.sample=lim.cos2.sample, palette=palette, ...)
    }
    
    if (plotVariable) {
        plotVariable(acp, axes=c(1,2), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
        plotVariable(acp, axes=c(1,3), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
        plotVariable(acp, axes=c(2,3), new.plot=!pdf, lab=lab.var, lim.cos2.var=lim.cos2.var, palette=palette, ...)
    }
    
    if (pdf) {
        par(op)
        dev.off()
    }
    res <- acp
    return(res)
}

#####
##
## PLS Partial Least Square functions
##
#####

PLS <- function(E,F,n=1, scale=TRUE){

    n.ind <- nrow(E)
    if(scale){
        E <- scale(E)*sqrt(n.ind/(n.ind-1))
        F <- scale(F)*sqrt(n.ind/(n.ind-1))
        param.scale <- attr(E,"scaled:scale")/sqrt(n.ind/(n.ind-1))
        param.center <- attr(E,"scaled:center")    
    }
    
    if(!scale){
        E <- scale(E,scale=FALSE)
        F <- scale(F,scale=FALSE)
        param.scale <- NULL
        param.center <- attr(E,"scaled:center")            
    }
    
    norme2 <- function(x){
        if(dim(x)[2]!=1)
            stop("x is not a vector")
        
        return(sum(x*x))
    }
        
    E0 <- E
    F0 <- F
    W <- matrix(0,ncol(E),n)
    T <- matrix(0,nrow(E),n)
    P <- matrix(0,ncol(E),n)
    C <- matrix(0,n,2)
    
    for(i in 1:n){
        print(paste("composante",i))
        tEF <- t(E0) %*% F0
        
        ## w1 est vecteur propre de tEFtFE associée à la plus grande valeur propre
        w1 <- t(E0)%*%F0
        w1 <- w1/sqrt(norme2(w1))
        W[,i] <- w1
        
        ## calcul de t1
        t1 <- E0%*%w1
        T[,i] <- t1
        C[i,1] <- (cov(F0,t1)*(n.ind-1)/n.ind)^2
        ## C[i,2] <- calc.cov.theo(E0,F0)
        
        ## on effectue deux régressions
        p1 <- t(E0)%*%t1/norme2(t1)
        P[,i] <- p1
        r1 <- t(F0)%*%t1/norme2(t1)

        ## mise a jour des matrices
        E0 <- E0 - t1%*%t(p1)
        F0 <- F0 - t1%*%t(r1)
    }
    
    Wstar <- W%*%ginv(t(P)%*%W)
    return(list(W=W, Wstar=Wstar, T=T, P=P, C=C, E0=E0, F0=F0, param.center=param.center, param.scale=param.scale))
}

#####
##
## MFA Multiple Factor Analysis
##
####

MFAreport <- function(resMFA, file.txt=NULL, file.pdf=NULL){
    
    if (!is.null(file.txt)){
        cat("GROUPS\n", file=file.txt)
        for (i in 1:length(resMFA$group)){
            cat(paste("\n", names(resMFA$group)[i], "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$group[[i]], file=file.txt, append=TRUE, sep="\t")
        }
        cat("\n\n\n", file=file.txt, append=TRUE)
        cat(rep("#", 25), file=file.txt, append=TRUE)
        cat("\nPartial axes\n", file=file.txt, append=TRUE)
        for (i in 1:(length(resMFA$partial.axes)-1)){
            cat(paste("\n", names(resMFA$partial.axes)[i], "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$partial.axes[[i]], file=file.txt, append=TRUE, sep="\t")
        }
        cat("\n cor between partial factors\n", file=file.txt, append=TRUE)
        write.infile(resMFA$partial.axes[[4]], file=file.txt, append=TRUE, sep="\t")
        cat(rep("#", 25), file=file.txt, append=TRUE)
        cat("\nIndividuals\n", file=file.txt, append=TRUE)
        for (i in 1:length(resMFA$ind)){
            cat(paste("\n", names(resMFA$ind[i]), "\n"), file=file.txt, append=TRUE)
            write.infile(resMFA$ind[[i]], file=file.txt, append=TRUE, sep="\t")
        }
    }
    
    if (!is.null(file.pdf)){
        pdf(file=file.pdf, height=8.26, width=11.69)
        plot(resMFA, axes=c(1,2), choix="ind", new.plot=FALSE, title="global")
        plot(resMFA, axes=c(3,4), choix="ind", new.plot=FALSE, title="global")
        plot(resMFA, axes=c(1,2), choix="ind", new.plot=FALSE, title="global with partial ind", partial="all")
        plot(resMFA, axes=c(3,4), choix="ind", new.plot=FALSE, title="global with partial ind", partial="all")
        for (i in 1:length(resMFA$separate.analyses)){
            tt <- try(plot(resMFA$separate.analyses[[i]], choix="ind", main=paste("partial -", names(resMFA$separate.analyses)[i]), new.plot=FALSE, axes=c(1,2)))
            tt <- try(plot(resMFA$separate.analyses[[i]], choix="ind", main=paste("partial -", names(resMFA$separate.analyses)[i]), new.plot=FALSE, axes=c(3,4)))
        }
        plot(resMFA, choix="group", new.plot=FALSE, axes=c(1,2))
        plot(resMFA, choix="group", new.plot=FALSE, axes=c(3,4))
        dev.off()
    }
}

runMFA <- function(Data, group=NULL, ncp=5, name.group=NULL, type=NULL, ind.sup = NULL, num.group.sup = NULL, graph = TRUE, report.file=NULL, report.pdf=NULL){
    
    ## List
    if (!is.data.frame(Data) & is.list(Data)){
        nbGrp <- length(Data)
        if (nbGrp < 2) stop("Provide at least 2 groups")
        group <- c()
        Data2MFA <- NULL
        nbInd <- nrow(Data[[1]])
        for (i in 1:nbGrp){
            if (nrow(Data[[i]]) != nbInd) stop("All data sets must have the same number of rows")
            Data2MFA <- cbind(Data2MFA, Data[[i]])
            group <- c(group, ncol(Data[[i]]))
        }
        if (is.null(name.group) & !is.null(names(Data))) name.group <- names(Data)
    }
    ## dataframe
    if (is.data.frame(Data) | is.matrix(Data)){   
        if (is.null(group)) stop("No groups provided")
        Data2MFA <- Data
        nbGrp <- length(group)
    }
    
    if (is.null(type)) type <- rep("s", nbGrp)
    
    resMFA <- MFA(Data2MFA, group=group, type=type, ind.sup=ind.sup, ncp=ncp, name.group=name.group, num.group.sup=num.group.sup, graph=graph, axes = c(1,2))
    
    if (!is.null(report.file)|| !is.null(report.pdf)){
        MFAreport(resMFA, file.txt=report.file, file.pdf=report.pdf)
    }
    
    return(resMFA)
}

