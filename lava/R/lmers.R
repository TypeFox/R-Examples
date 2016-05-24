##v <- lmerplot(l1,varcomp=TRUE,colorkey=TRUE,lwd=0,col=rainbow(20))
lmerplot <- function(model,x,id,y,transform,re.form=NULL,varcomp=FALSE,colorbar=TRUE,mar=c(4,4,4,6),col,...) {
    if (varcomp) {
        Z <- lme4::getME(model,"Z")
        nn <- unlist(lapply(lme4::getME(model,"Ztlist"),nrow))
        ve <- lme4::getME(model,"sigma")^2
        vu <- varcomp(model,profile=FALSE)$varcomp
        L <- Matrix::Diagonal(sum(nn),rep(vu,nn))
        V <- Z%*%L%*%(Matrix::t(Z))
        Matrix::diag(V) <- Matrix::diag(V)+ve
        cV <- Matrix::cov2cor(V)
        if (colorbar) { opt <- par(mar=mar) }
        ##if (missing(col)) col <- c("white",rev(heat.colors(16)))
        if (missing(col)) col <- rev(gray.colors(16,0,1))
        image(seq(nrow(V)),seq(ncol(V)),as.matrix(cV),xlab="",ylab="",col=col,zlim=c(0,1),...)
        if (colorbar) {
            uu <- devcoords()
            xrange <- c(uu$fig.x2,uu$dev.x2)
            xrange <- diff(xrange)/3*c(1,-1)+xrange
            yrange <- c(uu$fig.y1,uu$fig.y2)
            colorbar(direction="vertical",x.range=xrange,y.range=yrange,clut=col,values=seq(0,1,length.out=length(col)),srt=0,position=2)
            par(opt)
        }
        return(invisible(V))
    }
    if (missing(y)) y <- model.frame(model)[,1]
    yhat <- predict(model)
    if (!is.null(re.form)) ymean <- predict(model,re.form=re.form)
    if (!missing(transform)) {
        yhat <- transform(yhat)
        if (!is.null(re.form)) ymean <- transform(ymean)
        y <- transform(y)
    }
    plot(y ~ x, col=Col(id,0.3), pch=16,...)
    if (!is.null(re.form)) points(ymean ~ x, pch="-",cex=4);
    for (i in unique(id)) {
        idx <- which(id==i)
        lines(yhat[idx]~x[idx],col=i)
    }
}

varcomp <- function(x,profile=TRUE,...) {
    cc <- cbind(lme4::fixef(x),diag(as.matrix(vcov(x)))^.5)
    cc <- cbind(cc,cc[,1]-qnorm(0.975)*cc[,2],cc[,1]+qnorm(0.975)*cc[,2],
                2*(1-pnorm(abs(cc[,1])/cc[,2])))
    pr <- NULL
    if (profile) pr <- confint(x)
    colnames(cc) <- c("Estimate","Std.Err","2.5%","97.5%","p-value")
    vc <- lme4::VarCorr(x)
    res <- structure(list(coef=lme4::fixef(x), vcov=as.matrix(vcov(x)),
                          coefmat=cc,
                          confint=pr,
                          varcomp=vc[[1]][,],
                          residual=attributes(vc)$sc^2
                          ),
                     class="estimate.lmer")
    res
}
