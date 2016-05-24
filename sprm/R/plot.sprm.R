plot.sprm <-
function(x,type="yyp",alpha=.025,colors=list(bars="#0000AA",errorbars="red",background="#BBBBEE",abline="#21A0D2", scores="#0000AA",cutoffs="#00EEEE", badouts="darkred", modouts="black"),textsize=6,errorbar_width=1,data,yscale=NULL, ...)
  
  # plot.sprm automatically generates Sparse Partial Robust M regression Weights and Regression Coefficients Plots
  # inputs : x, a "sprm" class sparse PRM regression object 
  #          type, either "weights", "coefficients" or "yyp" 
  #          colors: a list containing the desired bar, errorbars, diagnal line and background colors as strings
  #          textsize: the text size in which to print the scores and loading names
  #          errorbar_width: a numeric containing the width of the errorbars on the regcoefficients  
  #          data (optional) a data frame containing new cases to predict and plot
  #          yscale, an optional scale vector for the yscale in the y vs y predicted plot (e.g. if two different regression plots have to be on the same scale)
  
  # written by Sven serneels, BASF Corp., January 2014. 
# Updated to include an y vs y predicted plot with case specific error bars, March 2014. 

{
#  require(ggplot2)
#  require(grid)
#  require(plsdof)
#  require(pcaPP)
#  source("dbshdyv.r")
  nam <- wt <- caseweights <- T2 <- y_original <- y_predicted <- llim <- ulim <- NULL
  if(!(class(x)=="sprm")){stop("The sparse PRM plot function can only be applied to sprm class objects")}
  if(!any(type==c("yyp", "weights", "coefficients", "dd"))){
    warning("Invalid plot type; type='yyp' was used.")
    type <- "yyp"
  }
  
  optcomp <- x$inputs$a 
  w <- x$w
  used.vars <- x$used.vars[[2*optcomp]]
  
  if(type=="weights"){
    plotweights <- data.frame(caseweights=w,nam=names(w))
    plotweights$nam <- factor(plotweights$nam, levels=names(w))
    print(ggplot(plotweights,aes(nam,caseweights)) + geom_bar(stat="identity",size=3,fill=colors$bars) + labs(title=paste(names(x$YMeans),"Sparse PRM Case Weights Plot")) +
            theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1.5),face="bold"),axis.text.x=element_text(angle=-90)))
  } else if (type=="dd") {
    if (!missing(data)){
      Xnames <- names(x$XMeans) #
      Xindex <- which(colnames(data)%in%Xnames) #
      if (length(Xindex)!=length(Xnames)){
        stop("Column names of data don't match variable names in the model.")
      }
      Xn <- scale(data[,Xindex], center=x$XMeans, scale=x$Xscales)
      if (length(rownames(data))==0){
        rownames(Xn) <- 1:dim(Xn)[1]
      } else{
        rownames(Xn) <- rownames(data) 
      }
      newscores <- Xn %*% x$R 
      plotscores <- rbind(as.matrix(x$scores[,1:optcomp]), as.matrix(newscores))
      Xscaled <- rbind(x$inputs$Xs, Xn)
    } else {
      plotscores <- as.matrix(x$scores[,1:optcomp])
      Xscaled <- x$inputs$Xs
    }
    plotscores <- as.matrix(x$scores[,1:optcomp])
    T2sprms <- diag((scale(plotscores,center=TRUE,scale=FALSE))%*%solve(cov(plotscores))%*%t(scale(plotscores,center=TRUE,scale=FALSE)))
    n <- length(T2sprms)
    cutoffT2 <- qchisq(1-alpha,optcomp)
    plotloadings <- (x$loadings[,1:optcomp])
    DModX <- x$inputs$Xs - plotscores%*%t(plotloadings)
    DModX <- as.matrix(DModX)
    DModX <- sqrt(diag(DModX%*%t(DModX)))
    DModX <- DModX^(2/3)
    locDModX <- median(DModX)
    scaDModX <- qn(DModX)
    cutoffDX <- qnorm(1-alpha,mean=locDModX,sd=scaDModX)
    plotdd <- data.frame(T2=T2sprms,DModX=DModX,cutoffT2=rep(cutoffT2,n),cutoffDX=rep(cutoffDX,n))
    plotty <- ggplot(plotdd,aes(T2,DModX)) + geom_text(size=textsize,color=colors$scores,label=rownames(plotdd)) + labs(title=paste("Sparse PRM Regression ", names(x$YMeans)," Distance-Distance Plot")) 
    plotty <- plotty + geom_hline(aes(yintercept=cutoffDX[2]),color=colors$cutoffs,linetype="dashed",size=2) 
    plotty <- plotty + geom_vline(aes(xintercept=cutoffT2[2]),color=colors$cutoffs,linetype="dashed",size=2) 
    plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1.5),face="bold"),
                             axis.text.x=element_text())
    print(plotty)  
  } else { 
    zerows <- which(w==0)
    if (length(zerows)==0){ 
      X <- x$inputs$Xs[,used.vars] * sqrt(w)
      y <- x$inputs$ys * sqrt(w)
    } else {
      X <- x$inputs$Xs[-zerows,used.vars] * sqrt(w[-zerows])
      y <- x$inputs$ys[-zerows] * sqrt(w[-zerows])
    }
    n <- length(y)
    p <- ncol(X)
    zb <- qnorm(1-alpha)
    Jaxx <- dbshdy(X,y,optcomp)
    J <- Jaxx$dbhdy
    bw <- Jaxx$b
    r <- y - scale(X)%*%bw[,optcomp]*sd(y) - mean(y)
    dYh <-  scale(X)%*%J + matrix(1,n,n)/n
    dr <- (diag(n) - dYh)
    df <- (sum(diag(t(dr)%*%dr)))
    sigma2 <- as.numeric(t(r)%*%r)/df
    covb <- (diag(1/(apply(X,2,sd)))%*%J%*%t(J)%*%diag(1/apply(X,2,sd)))
    covb2 <- diag(covb*sigma2)
    df <- (max(n,p)+1-sum(diag(dYh)))
    
    if(type=="coefficients"){
      b <- x$coefficients.scaled[used.vars,1]
      plotcoefs <- data.frame(coefficients=b,nam=names(b),llim=b -
                                zb*sqrt(covb2),ulim=b+zb*sqrt(covb2))
      plotcoefs$nam <- factor(plotcoefs$nam, levels=names(b)) 
      
      plotty <- ggplot(plotcoefs,aes(nam,coefficients)) + geom_bar(stat="identity",size=3,fill=colors$bars, position = "identity") + labs(title=paste(names(x$YMeans)," Sparse PRM Regression Coefficients Plot")) # added position = "identity" to get rid of the warning: Stacking not well defined when ymin != 0 
      plotty <- plotty + coord_flip() + geom_errorbar(aes(ymin=llim,ymax=ulim),width=errorbar_width,color=colors$errorbars)
      plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1),face="bold"),
                               axis.text.x=element_text(angle=-90),axis.title.x=element_blank(),axis.title.y=element_blank())
      print(plotty)
    } else {
      
      plotyyp <- intervals.sprm(x,data, optcomp, n, df, sigma2, covb, alpha)
 
      rsq <- 1 - sum(x$w*(plotyyp$y_predicted-plotyyp$y_original)^2)/(sum(x$w)*var(plotyyp$y_original))
      ynames <- names(x$YMeans)

      plotty <- ggplot(plotyyp)
      plotty <- plotty + geom_text(data=plotyyp[which(plotyyp$outs==0),], aes(x=y_original,y=y_predicted,label=nam),size=textsize,colour=colors$bars) 
      if(sum(plotyyp$outs==1)>0){
        plotty <- plotty + geom_text(data=plotyyp[which(plotyyp$outs==1),], aes(x=y_original,y=y_predicted,label=nam),size=textsize,colour=colors$modouts) 
      }
      if (sum(plotyyp$outs==2)>0){
        plotty <- plotty + geom_text(data=plotyyp[which(plotyyp$outs==2),], aes(x=y_original,y=y_predicted,label=nam),size=textsize,colour=colors$badouts) 
      }
      plotty <- plotty + geom_point(aes(x=y_original,y=y_predicted,label=nam),size=textsize,colour=NA) 
      plotty <- plotty + geom_errorbar(aes(x=y_original,y=y_predicted,ymin=llim,ymax=ulim),width=errorbar_width,colour=colors$errorbars)
      plotty <- plotty + labs(title=paste("SPRM Regression ", ynames," vs. ", ynames, " Predicted -", x$inputs$a," LV,  Rw = ",round(rsq,2)," df = ",round(df,1))) 
      plotty <- plotty + geom_abline(intercept=0,slope=1,color=colors$abline) 
      if(!is.null(yscale)){plotty <- plotty + ylim(yscale)}
      plotty <- plotty + theme(panel.background=element_rect(fill=colors$background),plot.title=element_text(size=rel(1)),
                               axis.text.x=element_text(),axis.title.x=element_blank(),axis.title.y=element_blank())
      print(plotty)
    }
  }}