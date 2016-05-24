star.nominalglob <-
function (formula, data, conf.int = FALSE, symmetric = TRUE, 
    pred.coding = "reference", printpvalues = TRUE, test.rel = TRUE, refLevel = 1, 
    maxit = 100, scale = TRUE, nlines = NULL, select = NULL, dist.x = 1, 
    dist.y = 1, dist.cov = 1, dist.cat = 1, xpd = TRUE, main = "", 
    lwd.stars = 1.2, col.fill = "gray90", col.circle = "black", lwd.circle = 1, 
    lty.circle = "longdash", lty.conf = "dotted", cex.labels = 1, cex.cat = 0.8, 
    xlim = NULL, ylim = NULL) 
{
    if(!is.data.frame(data))
      stop("Argument data has to be of class data.frame")
    if(!is.logical(test.rel))
      stop("Argument test.rel has to be of class logical")
    if(!is.logical(conf.int))
      stop("Argument conf.int has to be of class logical")
    if(!is.logical(symmetric))
      stop("Argument symmetric has to be of class logical")
    if(!(pred.coding %in% c("reference","effect")))
      stop("Argument pred.coding has to be reference or effect")
    if(!is.logical(printpvalues))
      stop("Argument printpvalues has to be of class logical")
    if(!is.logical(scale))
      stop("Argument scale has to be of class logical")
    if(!is.logical(xpd) & !is.na(xpd))
      stop("Argument xpd has to be of class logical or NA")  
    if(sum(grep("\\*",as.character(formula)[3]))>0)
      stop("Interactions can not be regarded automatically yet. Please create 
           Interactions in the data frame before calling the function!") 
    if(!is.numeric(dist.x))
      stop("Argument dist.x has to be numeric")
    if(!is.numeric(dist.y))
      stop("Argument dist.y has to be numeric") 
    if(!is.numeric(dist.cov))
      stop("Argument dist.cov has to be numeric")
        if(!is.numeric(dist.cat))
      stop("Argument dist.cat has to be numeric")
  
    is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
    
    if(!is.null(nlines)){
      if(!is.wholenumber(nlines) & length(nlines)==1){
      stop("Argument nlines has to be NULL or a single number of type integer")}}
    
    
    par(xpd = xpd)
    response <- all.vars(formula)[1]
    covariates <- all.vars(formula)[-1]
    resp <- data[, colnames(data) == response]
 
    if(is.ordered(resp))
      stop("The response may not be of class ordered!")
    
    if(!is.factor(resp))
      stop("The response has to be of class factor!")   
    
    if(!is.null(select) & !is.vector(select, mode = "numeric"))
      stop("Argument select has to be NULL or a numeric vector")

    if(!is.null(select)){select<-sort(unique(floor(select)))}
    
    
    keylab <- levels(resp)
    z <- length(keylab)
    n <- nrow(data)
    p <- length(covariates)
    
    if(!(refLevel %in% 1:z))
      stop(paste("Argument refLevel does not contain a valid reference category from vector 1:",z,sep=""))

    #######################
    if(pred.coding=="reference"){
    #######################
      
    # compute model + likelihood + variances
    multobj <- vglm(formula, data = data, family = multinomial(refLevel = refLevel), 
        maxit = maxit)
    coefmat <- t(coefficients(multobj, T))
    l <- ncol(coefmat)
    vcov <- vcov(multobj)
    sdsmat <- matrix(sqrt(diag(vcov)), nrow = l, byrow = TRUE)
    if (test.rel) {
        loglik <- logLik(multobj)
    }

    logxy <- c()
    isfac <- c()
    labvec <- ""
    covlab <- c()
    #loop: create labels, compute likelihoods for lr-tests
    for (u in 1:length(covariates)) {
        actcov <- data[, colnames(data) == covariates[u]]
        # create labels
        if (is.factor(actcov) | is.ordered(actcov)) {
            nlev <- length(levels(actcov))
            isfac <- c(isfac, rep(TRUE, (nlev - 1)))
            covlab <- c(covlab, rep(covariates[u], nlev - 1))
            if (is.ordered(actcov)) {
                parlevs <- paste("(", sort(as.numeric(unique(actcov)))[-nlev], 
                  ": ", levels(actcov)[-nlev], sep = "")
                labvec <- c(labvec, parlevs)
            }
            else {
                if (length(levels(actcov)) == 2) {
                  parlevs <- paste("(", sort(as.numeric(unique(actcov)) - 
                    1)[-1], ": ", levels(actcov)[-1], sep = "")
                }
                else {
                  parlevs <- paste("(", sort(as.numeric(unique(actcov)))[-1], 
                    ": ", levels(actcov)[-1], sep = "")
                }
                labvec <- as.expression(c(labvec, parlevs), font = 3)
            }
        }
        else {
            isfac <- c(isfac, FALSE)
            covlab <- c(covlab, covariates[u])
            labvec <- c(labvec, "")
        }
        # compute likelihoods
        if (test.rel) {
            if (is.factor(actcov) & length(levels(actcov)) > 
                2) {
                data2 <- data
                levlength <- length(levels(actcov))
                if (is.ordered((actcov))) {
                  data2 <- data
                  for (pp in 1:(levlength - 1)) {
                    covcat <- actcov
                    covcat[covcat == levels(covcat)[pp]] <- levels(covcat)[levlength]
                    data2[, colnames(data2) == covariates[u]] <- covcat
                    logxy <- c(logxy, logLik(vglm(formula = formula, 
                      data = data2, family = multinomial(refLevel = refLevel), 
                      maxit = maxit)))
                  }
                }
                else {
                  data2 <- data
                  for (pp in 2:levlength) {
                    covcat <- actcov
                    covcat[covcat == levels(covcat)[pp]] <- levels(covcat)[1]
                    data2[, colnames(data2) == covariates[u]] <- covcat
                    logxy <- c(logxy, logLik(vglm(formula = formula, 
                      data = data2, family = multinomial(refLevel = refLevel), 
                      maxit = maxit)))
                  }
                }
            }
            else {
                  if (length(covariates) >= 2){
                formxy <- as.formula(paste(response, "~", paste(covariates[-u], 
                  collapse = " + ")))
                }else{
                  formxy <- as.formula(paste(response, "~", 1))
                }
                logxy <- c(logxy, logLik(vglm(formula = formxy, 
                  data = data, family = multinomial(refLevel = refLevel), 
                  maxit = maxit)))
            }
        }
    }
    # compute lr-test p-values
    if (test.rel) {
        lrstats <- 2 * (-logxy + loglik)
        lrpvals <- 1 - pchisq(lrstats, df = z - 1)
        lrexact <- lrpvals
        lrpvals[lrpvals < 5e-04] <- 0
        lrpvals[lrpvals >= 5e-04 & lrpvals < 0.001] <- 0.001
        lrpvals <- format(lrpvals, digits = 1, nsmall = 3)
    }
    # compute coefficients and standard deviances for symmetric side constraints
    if (symmetric) {
        K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
        diag(K) <- (z - 1)/z
        cosym <- K %*% coefmat
        coef3 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% cosym
        coef2 <- matrix(0, ncol = l, nrow = z)
        coef2[refLevel, ] <- coef3
        coef2[-refLevel, ] <- cosym
        coefs <- exp(coef2)
            semat <- matrix(0, ncol = l, nrow = z)
            s1 <- 1
            for (e in 1:l) {
                s2 <- s1 + z - 2
                cov2 <- vcov[s1:s2, s1:s2]
                cov1 <- K %*% cov2 %*% t(K)
                se1 <- diag(cov1)
                se2 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% 
                  cov1 %*% matrix(rep(-1, z - 1), ncol = 1)
                se3 <- rep(0, z)
                se3[refLevel] <- se2
                se3[-refLevel] <- se1
                semat[, e] <- sqrt(se3)
                s1 <- s2 + 1
            }
    }
    else {
      # coefficients and standard deviations with reference category
        coefs <- matrix(0, ncol = l, nrow = z)
        coefs[refLevel, ] <- rep(1, l)
        coefs[-refLevel, ] <- exp(coefmat)
            semat <- matrix(0, ncol = l, nrow = z)
            semat[refLevel, ] <- rep(0, l)
            semat[-refLevel, ] <- t(sdsmat)
    }
    x1 <- t(coefs)
    # compute confidence intervals
    if (printpvalues | conf.int) {
        x2 <- exp(log(x1) + qnorm(0.975) * t(semat))
        x3 <- exp(log(x1) - qnorm(0.975) * t(semat))
    }
    
    if(is.null(select)){
      lselect<-l
      select<-1:l
    }else{
      if(sum(select>l)>0){select<-select[select<=l]}
      lselect<-length(select)}
      plusint <- select==1
    
    # compute factors for scaling
    fac <- max(x1)
    facs <- apply(x1, 1, max)/fac
    # compute number of lines
   if (is.null(nlines)){
     if (floor(sqrt(lselect)) == sqrt(lselect)) {
            nlines <- sqrt(lselect)
        }
        else {
            nlines <- (floor(sqrt(lselect)) + 1)
        }
}
    l2 <- nlines^2
    ###########
    facmat <- matrix(facs, nrow = nrow(x1), ncol = ncol(x1))
    odds <- t(x1)
    # compute p-values
    if (printpvalues) {
        pvals <- 2*(1 - pnorm(abs(log(odds)/semat)))
    }
    # scaling
    if (conf.int) {
        if (scale) {
            x2 <- x2/facmat
            x3 <- x3/facmat
        }
    }
    if (scale) {
        x1 <- x1/facmat
    }
    # compute locations and limits
    plotmat<-x1[select,]
    if(lselect==1){plotmat<-t(as.matrix(plotmat))}
    loc <- stars(plotmat, nrow = nlines, scale = FALSE, len = 1, cex = 1.2, 
        lwd = 1, flip.labels = F, plot = FALSE) * fac
    loc[, 1] <- loc[, 1] * dist.x * 2
    loc[, 2] <- loc[, 2] * 1.9 * dist.y
    if (conf.int) {
        loc <- loc * 0.9 * (max(x2)/max(x1))
    }
    if (is.null(ylim)) {
        ylim <- c(min(loc[, 2]) - 1.5 * max(x1), max(loc[, 2]) + 
            2 * max(x1))
        if (conf.int) {
            ylim <- c(min(loc[, 2]) - 1.5 * max(x2), max(loc[, 
                2]) + 2 * max(x2))
        }
        
    }
    if (is.null(xlim)) {
        xlim <- c(min(loc[, 1]) - 1.6 * max(x1), max(loc[, 1]) + 
            1.6 * max(x1))
       }
    
    ###############

    # plot stars with or without confidence intervals
    if (!conf.int) {
        par(cex = 1, font = 3)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loc, 
            xlim = xlim, ylim = ylim, labels = "", main = main, add = FALSE)
        par(cex = 1, font = 1)
    }
    if (conf.int) {
    plotmat2<-x2[select,]
    if(lselect==1){plotmat2<-t(as.matrix(plotmat2))}
    plotmat3<-x3[select,]
    if(lselect==1){plotmat3<-t(as.matrix(plotmat3))}
        par(cex = 1, font = 3)

        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", ylim = ylim, xlim = xlim, 
            labels = "", locations = loc, main = main, add = FALSE)

        par(cex = 1, font = 1)
    }
    # compute angles in stars, plot category labels
    cosis <- c()
    sinis <- c()
    alts <- 0:(z - 1)
    angle <- alts * ((360/z) * pi)/180
    cosis <- cos(angle)
    sinis <- sin(angle)
    xlocs <- ylocs <- matrix(0, ncol = l, nrow = z)
    index <- 1
    xx<-x1
    if(conf.int){xx<-x2}
    disfac<-0.4
    if(!scale){disfac<-0.3}
    dis2<-max(xx * disfac * dist.cat, xx * disfac * dist.cat)
    for (co in select) {
        dir <- x1[co, ]
        if (conf.int) {
            dir <- x2[co, ]
        }
        locx <- loc[index, 1]
        locy <- loc[index, 2]
        xlocs[, co] <- rep(locx, z) + cosis * dir + cosis * dis2
        ylocs[, co] <- rep(locy, z) + sinis * dir + sinis * dis2
        index <- index + 1
    }
    # create labels with p-values
    if (printpvalues) {
        printps <- pvals
        printps[printps < 5e-04] <- 0
        printps[printps >= 5e-04 & printps < 0.001] <- 0.001
        printps <- format(printps, digits = 1, nsmall = 3)
        starlab <- paste(keylab, paste("(", printps, ")", sep = ""), 
            sep = "\n")
        starlab[is.nan(pvals)] <- keylab[refLevel]
    }
    else {
        starlab <- rep(keylab, l)
    }

    
         if (symmetric) {
      if(scale){
        radii <- 1/facs 
      }else{
       radii<-rep(1,l)
      }
    }
    else {
        radii <- x1[1:l, refLevel]
    }

    # plot cirlces
    if(!identical(select,1)){
    par(lty = lty.circle, lwd = lwd.circle, col = col.circle)
    symbols(c(loc[, 1]), c(loc[, 2]), circles = radii[select], 
        add = TRUE, inches = FALSE, fg = col.circle, bg = col.fill)    
    par(lty = 1, lwd = 1, col = 1)
    }
    ###################

    # plot stars with or without confidence intervals
    if (!conf.int) {
        par(cex = 1, font = 3)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loc, 
            xlim = xlim, ylim = ylim, labels = "", main = main, add = TRUE)
        par(cex = 1, font = 1)
    }
    if (conf.int) {
    plotmat2<-x2[select,]
    if(lselect==1){plotmat2<-t(as.matrix(plotmat2))}
    plotmat3<-x3[select,]
    if(lselect==1){plotmat3<-t(as.matrix(plotmat3))}
        par(cex = 1, font = 3)
        stars(plotmat2, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = lty.conf, ylim = ylim, xlim = xlim, 
            labels = "", locations = loc, main = main, add = TRUE)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", add = TRUE, flip.labels = FALSE, 
            locations = loc, labels = "")
    stars(plotmat3, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = lty.conf, locations = loc, add = TRUE, 
            labels = "")
        par(cex = 1, font = 1, col = 1)
    }
    
    # print labels
    text(x = xlocs[,select], y = ylocs[,select], labels = starlab[c(xlocs)!=0], cex = cex.cat, 
        font = 1)
   ######!!!!
    # create labels
    if (test.rel) {
        labvec[-1] <- paste(labvec[-1], ", ", lrpvals, ")", sep = "")
        labvec[-1][!isfac] <- paste("(", lrpvals[!isfac], ")", 
            sep = "")
    }
    else {
        labvec[-1][isfac] <- paste(labvec[-1][isfac], ")", sep = "")
    }
    if (conf.int) {
        fac <- max(x2)
    }
    plotlabs <- paste(c(labvec[1], covlab), "\n", c(labvec, rep("", (nrow(loc) - 1))), sep = "")
    plotlabs[1] <- paste("Intercept", "\n", "", sep = "")
    text(y = loc[, 2] + fac * 1.9 * (l2 + 20)/32 * dist.cov, x = loc[, 1], 
        labels = plotlabs[select], font = 1, cex = cex.labels)
    if (printpvalues) {
        rownames(pvals) <- keylab
        colnames(pvals) <- colnames(coefmat)
    }
    rownames(semat) <- keylab
    colnames(semat) <- colnames(coefmat)
    rownames(odds) <- keylab
    colnames(odds) <- colnames(coefmat)
# create result lists   
    if (test.rel) {
        p_rel <- matrix(lrexact, ncol = l - 1)
        colnames(p_rel) <- colnames(coefmat)[-1]
        if (printpvalues) {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, pvalues = pvals,
               p_rel = p_rel, xlim = xlim, ylim = ylim)
        }
        else {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, p_rel = p_rel, xlim = xlim, 
                ylim = ylim)
        }
    }
    else {
        if (!printpvalues) {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, xlim = xlim, ylim = ylim)
        }
        else {
            returns <- list(odds = odds, coefficients = log(odds), 
                se = semat, pvalues = pvals, xlim = xlim, ylim = ylim)
        }
    }
   
    returns
    ###################
    }else{
    ###################
      isfac <- c()
      dummydat <-dummydat2<- resp
      outputlab<-c()
      
      for (j in 1:p) {
        
        # dummy coding
        variab <- data[, colnames(data) == covariates[j]]
        levlen <- length(levels(variab))
        if (is.factor(variab) & levlen > 2) {
          outputlab<-c(outputlab,levels(variab))
          if (!is.ordered(variab)) {
            variab <- as.numeric(variab)
            help <- matrix(0, nrow = nrow(data), ncol = levlen - 1)
            help[variab==levlen,]<--1
            for (ooo in 1:(levlen-1)) {
              help[variab == ooo, ooo] <- 1
            }
            help2 <- matrix(0, nrow = nrow(data), ncol = levlen - 1)
            help2[variab==1,]<--1
            for (ooo in 2:(levlen)) {
              help2[variab == ooo, ooo-1] <- 1
            } 
          }
          else {
            help <- matrix(0, nrow = nrow(data), ncol = levlen - 
              1)
            for (ooo in 1:(levlen - 1)) {
              help[variab == ooo, ooo] <- 1
            }
            help2<-help
          }
          dummydat <- cbind(dummydat, help)
          dummydat2 <- cbind(dummydat2, help2)
        }
        else {
          if (is.factor(variab)){
            outputlab<-c(outputlab,levels(variab)[2])
            if(sum(levels(as.factor(as.numeric(variab))) == 
              c("1", "2")) == 2) {
              variab <- as.numeric(variab) - 1
            }}
          dummydat <- cbind(dummydat, variab)
          dummydat2 <- cbind(dummydat2, variab)
        }
      }
      
      dummydat <- as.data.frame(dummydat)
      colnames(dummydat) <- paste("V", 1:ncol(dummydat), sep = "")
      
      dummydat2 <- as.data.frame(dummydat2)
      colnames(dummydat2) <- paste("V", 1:ncol(dummydat2), sep = "")
      
      
      multi.df<-rep(0,length(covariates))
      multi<-c()
      multi2<-c()
      vb<-1
      
      for(u in 1:length(covariates)){
        check.cat<-data[, colnames(data) == covariates[u]]
        if(is.factor(check.cat)){
          if(length(levels(check.cat))>2){
            multi.df[u]<-length(levels(check.cat))
            multi<-c(multi,rep(vb,length(levels(check.cat))-1))
            multi2<-c(multi2,rep(length(levels(check.cat))-1,length(levels(check.cat))-1))
            vb<-vb+1
          }else{multi<-c(multi,0)
                multi2<-c(multi2,0)}
        }else{
          multi<-c(multi,0)
          multi2<-c(multi2,0)}
      }
      
      formdummy<-as.formula(paste("V1~",paste(paste("V",2:ncol(dummydat),sep=""),collapse="+")))
      
      # compute model + likelihood + variances
      multobj <- vglm(formdummy, data = dummydat, family = multinomial(refLevel = refLevel), 
                      maxit = maxit)
      multobj2 <- vglm(formdummy, data = dummydat2, family = multinomial(refLevel = refLevel), 
                       maxit = maxit)
      
      coefmat <- t(coefficients(multobj, T))
      l <- ncol(coefmat)
      vcov <- vcov(multobj)
      
      coef.after<-coefmat[,1]
      l.after<-l+vb-1
      vcov.after<-matrix(0,ncol=l.after*(z-1),nrow=l.after*(z-1))
      vcov.after[1:(z-1),1:(z-1)]<-vcov[1:(z-1),1:(z-1)]
      s1 <- z
      vb<-1
      for (u in 1:length(multi)){
        if(multi[u]!=0 & multi[u]==vb){
          ku<-multi2[u]
          
          T.mat<-matrix(rep(-diag(z-1),ku),nrow=z-1,byrow=FALSE)
          coef.h<-coefmat[,-1]
          coef.h<-coef.h[,multi==vb]
          
          coef.add<-T.mat%*%matrix(c(coef.h),ncol=1)
          coef.after<-cbind(coef.after,coef.h,coef.add)
          
          s2 <- s1 + length(c(coef.h)) -1
          cov2 <- vcov[s1:s2, s1:s2]
          cov1 <- T.mat %*% cov2 %*% t(T.mat)
          cov3<-cbind(cov2,matrix(0,ncol=z-1,nrow=nrow(cov2)))
          cov4<-rbind(cov3,cbind(matrix(0,nrow=z-1,ncol=ncol(cov2)),cov1))    
          vcov.after[s1:(s1+ncol(cov4)-1),s1:(s1+ncol(cov4)-1)]<-cov4
          
          s1 <- s1+ncol(cov4)
          vb<-vb+1
        }else{
          if(multi[u]==0){
            s2<-s1+z-2
            coef.after<-cbind(coef.after,coefmat[,u+1])
            vcov.after[s1:s2,s1:s2]<-vcov[(s1:s2)-(vb-1)*(z-1),(s1:s2)-(vb-1)*(z-1)]
            s1<-s2+1
          }}
      }
      
      
      l<-l.after
      vcov<-vcov.after
      coefmat<-coef.after
      
      sdsmat <- matrix(sqrt(diag(vcov)), nrow = l, byrow = TRUE)
      if (test.rel) {
        loglik <- logLik(multobj)
      }
      
      logxy <- c()
      isfac <- c()
      labvec <- ""
      covlab <- c()
      #loop: create labels, compute likelihoods for lr-tests
      for (u in 1:length(covariates)) {
        actcov <- data[, colnames(data) == covariates[u]]
        # create labels
        if (is.factor(actcov) | is.ordered(actcov)) {
          nlev <- length(levels(actcov))
          isfac <- c(isfac, rep(TRUE, (nlev - 1 + I(nlev>2))))
          if(nlev==2){
            covlab <- c(covlab, rep(covariates[u], nlev -1))}
          if (is.ordered(actcov)) {
            parlevs <- paste("(", sort(as.numeric(unique(actcov)))[-nlev], 
                             ": ", levels(actcov)[-nlev], sep = "")
            labvec <- c(labvec, parlevs)
          }
          else {
            if (length(levels(actcov)) == 2) {
              parlevs <- paste("(", sort(as.numeric(unique(actcov)) - 
                1)[-1], ": ", levels(actcov)[-1], sep = "")
            }
            else {
              covlab <- c(covlab, rep(covariates[u], nlev ))
              parlevs <- paste("(", sort(as.numeric(unique(actcov))), 
                               ": ", levels(actcov), sep = "")
            }
            labvec <- as.expression(c(labvec, parlevs), font = 3)
          }
        }
        else {
          isfac <- c(isfac, FALSE)
          covlab <- c(covlab, covariates[u])
          labvec <- c(labvec, "")
        }}
      # compute likelihoods
      if (test.rel) {
        for (u in 1:length(multi)){
          
          if (length(multi) >= 2){
            formxy<-as.formula(paste("V1~",paste(paste("V",(2:ncol(dummydat))[-u],sep=""),collapse="+")))
          }else{
            formxy <- as.formula(paste(response, "~", 1))
          }
          
          log.add<-logLik(vglm(formula = formxy, 
                               data = dummydat, family = multinomial(refLevel = refLevel), 
                               maxit = maxit))
          if(u ==length(multi) & multi[u]!=0){
            form.sym<-as.formula(paste("V1~",paste(paste("V",(2:ncol(dummydat))[-u],sep=""),collapse="+")))
            log.add<-c(log.add,logLik(vglm(formula = form.sym, 
                                           data = dummydat2, family = multinomial(refLevel = refLevel), 
                                           maxit = maxit)))
          }else{
            if(multi[u]!=0 & multi[u]!=multi[u+1]){
              form.sym<-as.formula(paste("V1~",paste(paste("V",(2:ncol(dummydat))[-u],sep=""),collapse="+")))
              log.add<-c(log.add,logLik(vglm(formula = form.sym, 
                                             data = dummydat2, family = multinomial(refLevel = refLevel), 
                                             maxit = maxit)))
            }}
          logxy <- c(logxy,log.add)
        }
      }
      
      # compute lr-test p-values
      if (test.rel) {
        lrstats <- 2 * (-logxy + loglik)
        lrpvals <- 1 - pchisq(lrstats, df = z - 1)
        lrexact <- lrpvals
        lrpvals[lrpvals < 5e-04] <- 0
        lrpvals[lrpvals >= 5e-04 & lrpvals < 0.001] <- 0.001
        lrpvals <- format(lrpvals, digits = 1, nsmall = 3)
      }
      # compute coefficients and standard deviances for symmetric side constraints
      if (symmetric) {
        K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
        diag(K) <- (z - 1)/z
        cosym <- K %*% coefmat
        coef3 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% cosym
        coef2 <- matrix(0, ncol = l, nrow = z)
        coef2[refLevel, ] <- coef3
        coef2[-refLevel, ] <- cosym
        ###
        
        ###
        coefs <- exp(coef2)
        semat <- matrix(0, ncol = l, nrow = z)
        s1 <- 1
        for (e in 1:l) {
          s2 <- s1 + z - 2
          cov2 <- vcov[s1:s2, s1:s2]
          cov1 <- K %*% cov2 %*% t(K)
          se1 <- diag(cov1)
          se2 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% 
            cov1 %*% matrix(rep(-1, z - 1), ncol = 1)
          se3 <- rep(0, z)
          se3[refLevel] <- se2
          se3[-refLevel] <- se1
          semat[, e] <- sqrt(se3)
          s1 <- s2 + 1
        }
      }
      else {
        # coefficients and standard deviations with reference category
        coefs <- matrix(0, ncol = l, nrow = z)
        coefs[refLevel, ] <- rep(1, l)
        coefs[-refLevel, ] <- exp(coefmat)
        semat <- matrix(0, ncol = l, nrow = z)
        semat[refLevel, ] <- rep(0, l)
        semat[-refLevel, ] <- t(sdsmat)
      }
      x1 <- t(coefs)
      # compute confidence intervals
      if (printpvalues | conf.int) {
        x2 <- exp(log(x1) + qnorm(0.975) * t(semat))
        x3 <- exp(log(x1) - qnorm(0.975) * t(semat))
      }
      
      if(is.null(select)){
        lselect<-l
        select<-1:l
      }else{
        if(sum(select>l)>0){select<-select[select<=l]}
        lselect<-length(select)}
      plusint <- select==1
      
      # compute factors for scaling
      fac <- max(x1)
      facs <- apply(x1, 1, max)/fac
      # compute number of lines
      if (is.null(nlines)){
        if (floor(sqrt(lselect)) == sqrt(lselect)) {
          nlines <- sqrt(lselect)
        }
        else {
          nlines <- (floor(sqrt(lselect)) + 1)
        }
      }
      l2 <- nlines^2
      ###########
      facmat <- matrix(facs, nrow = nrow(x1), ncol = ncol(x1))
      odds <- t(x1)
      # compute p-values
      if (printpvalues) {
        pvals <- 2*(1 - pnorm(abs(log(odds)/semat)))
      }
      # scaling
      if (conf.int) {
        if (scale) {
          x2 <- x2/facmat
          x3 <- x3/facmat
        }
      }
      if (scale) {
        x1 <- x1/facmat
      }
      # compute locations and limits
      plotmat<-x1[select,]
      if(lselect==1){plotmat<-t(as.matrix(plotmat))}
      loc <- stars(plotmat, nrow = nlines, scale = FALSE, len = 1, cex = 1.2, 
                   lwd = 1, flip.labels = F, plot = FALSE) * fac
      loc[, 1] <- loc[, 1] * dist.x * 2
      loc[, 2] <- loc[, 2] * 1.9 * dist.y
      if (conf.int) {
        loc <- loc * 0.9 * (max(x2)/max(x1))
      }
      if (is.null(ylim)) {
        ylim <- c(min(loc[, 2]) - 1.5 * max(x1), max(loc[, 2]) + 
          2 * max(x1))
        if (conf.int) {
          ylim <- c(min(loc[, 2]) - 1.5 * max(x2), max(loc[, 
                                                           2]) + 2 * max(x2))
        }
        
      }
      if (is.null(xlim)) {
        xlim <- c(min(loc[, 1]) - 1.6 * max(x1), max(loc[, 1]) + 
          1.6 * max(x1))
      }
      
      ###############
      
      # plot stars with or without confidence intervals
      if (!conf.int) {
        par(cex = 1, font = 3)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loc, 
              xlim = xlim, ylim = ylim, labels = "", main = main, add = FALSE)
        par(cex = 1, font = 1)
      }
      if (conf.int) {
        plotmat2<-x2[select,]
        if(lselect==1){plotmat2<-t(as.matrix(plotmat2))}
        plotmat3<-x3[select,]
        if(lselect==1){plotmat3<-t(as.matrix(plotmat3))}
        par(cex = 1, font = 3)
        
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = "solid", ylim = ylim, xlim = xlim, 
              labels = "", locations = loc, main = main, add = FALSE)
        
        par(cex = 1, font = 1)
      }
      # compute angles in stars, plot category labels
      cosis <- c()
      sinis <- c()
      alts <- 0:(z - 1)
      angle <- alts * ((360/z) * pi)/180
      cosis <- cos(angle)
      sinis <- sin(angle)
      xlocs <- ylocs <- matrix(0, ncol = l, nrow = z)
      index <- 1
      xx<-x1
      if(conf.int){xx<-x2}
      disfac<-0.4
      if(!scale){disfac<-0.3}
      dis2<-max(xx * disfac * dist.cat, xx * disfac * dist.cat)
      for (co in select) {
        dir <- x1[co, ]
        if (conf.int) {
          dir <- x2[co, ]
        }
        locx <- loc[index, 1]
        locy <- loc[index, 2]
        xlocs[, co] <- rep(locx, z) + cosis * dir + cosis * dis2
        ylocs[, co] <- rep(locy, z) + sinis * dir + sinis * dis2
        index <- index + 1
      }
      # create labels with p-values
      if (printpvalues) {
        printps <- pvals
        printps[printps < 5e-04] <- 0
        printps[printps >= 5e-04 & printps < 0.001] <- 0.001
        printps <- format(printps, digits = 1, nsmall = 3)
        starlab <- paste(keylab, paste("(", printps, ")", sep = ""), 
                         sep = "\n")
        starlab[is.nan(pvals)] <- keylab[refLevel]
      }
      else {
        starlab <- rep(keylab, l)
      }
      
      
      if (symmetric) {
        if(scale){
          radii <- 1/facs 
        }else{
          radii<-rep(1,l)
        }
      }
      else {
        radii <- x1[1:l, refLevel]
      }
      
      # plot cirlces
      if(!identical(select,1)){
        par(lty = lty.circle, lwd = lwd.circle, col = col.circle)
        symbols(c(loc[, 1]), c(loc[, 2]), circles = radii[select], 
                add = TRUE, inches = FALSE, fg = col.circle, bg = col.fill)    
        par(lty = 1, lwd = 1, col = 1)
      }
      ###################
      
      # plot stars with or without confidence intervals
      if (!conf.int) {
        par(cex = 1, font = 3)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loc, 
              xlim = xlim, ylim = ylim, labels = "", main = main, add = TRUE)
        par(cex = 1, font = 1)
      }
      if (conf.int) {
        plotmat2<-x2[select,]
        if(lselect==1){plotmat2<-t(as.matrix(plotmat2))}
        plotmat3<-x3[select,]
        if(lselect==1){plotmat3<-t(as.matrix(plotmat3))}
        par(cex = 1, font = 3)
        stars(plotmat2, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = lty.conf, ylim = ylim, xlim = xlim, 
              labels = "", locations = loc, main = main, add = TRUE)
        stars(plotmat, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = "solid", add = TRUE, flip.labels = FALSE, 
              locations = loc, labels = "")
        stars(plotmat3, scale = FALSE, draw.segments = FALSE, len = 1, 
              lwd = lwd.stars, lty = lty.conf, locations = loc, add = TRUE, 
              labels = "")
        par(cex = 1, font = 1, col = 1)
      }
      
      # print labels
      text(x = xlocs[,select], y = ylocs[,select], labels = starlab[c(xlocs)!=0], cex = cex.cat, 
           font = 1)
      ######!!!!
      # create labels
      if (test.rel) {
        labvec[-1] <- paste(labvec[-1], ", ", lrpvals, ")", sep = "")
        labvec[-1][!isfac] <- paste("(", lrpvals[!isfac], ")", 
                                    sep = "")
      }
      else {
        labvec[-1][isfac] <- paste(labvec[-1][isfac], ")", sep = "")
      }
      if (conf.int) {
        fac <- max(x2)
      }
      
      plotlabs <- paste(c(labvec[1], covlab), "\n", c(labvec, rep("", (nrow(loc) - 1))), sep = "")
      plotlabs[1] <- paste("Intercept", "\n", "", sep = "")
      text(y = loc[, 2] + fac * 1.9 * (l2 + 20)/32 * dist.cov, x = loc[, 1], 
           labels = plotlabs[select], font = 1, cex = cex.labels)
      
      covlab[isfac]<-paste(covlab[isfac],outputlab,sep="")
      covlab<-c("Intercept",covlab)
      colnames(coefmat)<-covlab
      if (printpvalues) {
        rownames(pvals) <- keylab
        colnames(pvals) <- colnames(coefmat)
      }
      rownames(semat) <- keylab
      colnames(semat) <- colnames(coefmat)
      rownames(odds) <- keylab
      colnames(odds) <- colnames(coefmat)
      # create result lists   
      if (test.rel) {
        p_rel <- matrix(lrexact, ncol = l - 1)
        colnames(p_rel) <- colnames(coefmat)[-1]
        if (printpvalues) {
          returns <- list(odds = odds, 
                          coefficients = log(odds), se = semat, pvalues = pvals,
                          p_rel = p_rel, xlim = xlim, ylim = ylim)
        }
        else {
          returns <- list(odds = odds, 
                          coefficients = log(odds), se = semat, p_rel = p_rel, xlim = xlim, 
                          ylim = ylim)
        }
      }
      else {
        if (!printpvalues) {
          returns <- list(odds = odds, 
                          coefficients = log(odds), se = semat, xlim = xlim, ylim = ylim)
        }
        else {
          returns <- list(odds = odds, coefficients = log(odds), 
                          se = semat, pvalues = pvals, xlim = xlim, ylim = ylim)
        }
      }
      
      returns
    }
    
}
