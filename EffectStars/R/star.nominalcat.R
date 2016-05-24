star.nominalcat <-
function (formula, data, xij, conf.int = FALSE, symmetric = FALSE, 
    pred.coding = "reference", printpvalues = TRUE, test.rel = TRUE, refLevel = 1, 
    maxit = 100, scale = TRUE, nlines = NULL, select = NULL, catstar = TRUE, 
    dist.x = 1, dist.y = 1, dist.cov = 1, dist.cat = 1, xpd = TRUE, main = "", 
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
    
    pcat<-length(xij)
    pglob<-p - pcat
    covcateg <- lapply(xij,function(a){return(all.vars(a)[1])})
    covglob <- covglob2 <- covariates[!(covariates %in% covcateg)]
    formglob <- paste(covglob,collapse="+")
  
    xijs<-c()
  for(iii in 1:pcat){
  xijs<-c(xijs,all.vars(xij[[iii]]))
  }
    
    ###############################
    if(pred.coding=="reference"){
      #############################
    
    form2<-as.formula(paste("~",paste(unique(c(xijs,covariates)),collapse="+"),sep=""))
    
    form2vars <- all.vars(form2)
    
    if(!(refLevel %in% 1:z))
      stop(paste("Argument refLevel does not contain a valid reference category from vector 1:",z,sep=""))
    
######################

    # compute model + likelihood + variances
    multobj <- vglm(formula, data = data, family = multinomial(parallel = 
      as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
      xij = xij, form2 = form2, trace = FALSE)
 
    coefmat <- t(coefficients(multobj, T))
    
    l <- ncol(coefmat)
    indcat<-apply(coefmat, 2,var)==0
    indcat2 <- rep(indcat, each=z-1)
    lcat <- sum(indcat)
    lglob <- l - lcat
    indout <- rep(((1:length(indcat))[indcat]-1)*(z-1),each=z-2)+(2:(z-1))
    indcat2 <-indcat2[-indout]
    vcovtotal <- vcov(multobj)
    vcov <- vcov(multobj)[!indcat2,!indcat2]
    sdsmat <- matrix(sqrt(diag(vcov)), nrow = lglob, byrow = TRUE)
    
    catspec<-coefmat[1,indcat]
catnamesresults<-names(catspec)
    if(pcat==1){
      catspecses<-sqrt(vcov(multobj)[indcat2,indcat2])
                  }else{
    catspecses <- sqrt(diag(vcov(multobj)[indcat2,indcat2]))}
    if (test.rel) {
        loglik <- logLik(multobj)
    }
    logxy <- c()
    isfac <- c()
    labvec <- ""
    covlab <- c()
    #loop: create labels, compute likelihoods for lr-tests
    for (u in 1:pglob) {
        actcov <- data[, colnames(data) == covglob[u]]
        # create labels
        if (is.factor(actcov) | is.ordered(actcov)) {
            nlev <- length(levels(actcov))
            isfac <- c(isfac, rep(TRUE, (nlev - 1)))
            covlab <- c(covlab, rep(covglob[u], nlev - 1))
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
            covlab <- c(covlab, covglob[u])
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
                    data2[, colnames(data2) == covglob[u]] <- covcat
                    mod2<-vglm(formula, data = data2, family = multinomial(parallel = 
                    as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                    xij = xij, form2 = form2, trace = FALSE)
                    logxy <- c(logxy, logLik(mod2))
                  }
                }
                else {
                  data2 <- data
                  for (pp in 2:levlength) {
                    covcat <- actcov
                    covcat[covcat == levels(covcat)[pp]] <- levels(covcat)[1]
                    data2[, colnames(data2) == covglob[u]] <- covcat
                    mod2<-vglm(formula, data = data2, family = multinomial(parallel = 
                    as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                    xij = xij, form2 = form2, trace = FALSE)
                    logxy <- c(logxy, logLik(mod2))
                  }
                }
            }
            else {
                  if (length(covglob) >= 2){
                formxy <- as.formula(paste(response, "~", paste(covariates[covariates!=covglob[u]], 
                  collapse = " + ")))
                form22<-as.formula(paste("~",paste(form2vars[form2vars!=covglob[u]],collapse="+"),sep=""))
                formglob2<-paste(covglob[-u],collapse="+")
                }else{
                  formxy <- as.formula(paste(response, "~", 1))
                }
                  mod2<-vglm(formxy, data = data, family = multinomial(parallel = 
                    as.formula(paste("FALSE ~", formglob2,"")) , refLevel = refLevel), maxit = maxit,
                    xij = xij, form2 = form22, trace = FALSE)
                    logxy <- c(logxy, logLik(mod2))
            }
        }
    }
    ######################################### LR Tests for catspec covariates
    likecat<-c()  
    if (test.rel) {
        for(qe in 1:pcat){
          actcov <- data[, colnames(data) == covcateg[qe]]
                if (length(covcateg) >= 2){
                formxy <- as.formula(paste(response, "~", paste(covariates[covariates!=covcateg[qe]], 
                  collapse = " + ")))
                form22<-as.formula(paste("~",paste(form2vars[form2vars!=covcateg[qe]],collapse="+"),sep=""))
                xij2<-xij[-qe]
                mod2<-vglm(formxy, data = data, family = multinomial(parallel = 
                    as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                    xij = xij2, form2 = form22, trace = FALSE)
                }else{
                  formxy <- as.formula(paste(response, "~", paste(covariates[covariates!=covcateg[qe]], 
                  collapse = " + ")))
                  mod2<-vglm(formxy, data = data, family = multinomial(refLevel = refLevel), 
        maxit = maxit)
                }
                  
                    likecat <- c(likecat, logLik(mod2))
        }
      }

    #########################################
    # compute lr-test p-values
    if (test.rel) {
        lrstats <- 2 * (-logxy + loglik)
        lrpvals <- 1 - pchisq(lrstats, df = z - 1)
        lrexact <- lrpvals
        lrpvals[lrpvals < 5e-04] <- 0
        lrpvals[lrpvals >= 5e-04 & lrpvals < 0.001] <- 0.001
        lrpvals <- format(lrpvals, digits = 1, nsmall = 3)
    }
   if (test.rel) {
        catstats <- 2 * (-likecat + loglik)
        catpvals <- 1 - pchisq(catstats, df = 1)
        catexact <- catpvals
        catpvals[catpvals < 5e-04] <- 0
        catpvals[catpvals >= 5e-04 & catpvals < 0.001] <- 0.001
        catpvals <- format(catpvals, digits = 1, nsmall = 3)
    }

    # compute coefficients and standard deviances for symmetric side constraints
    if (symmetric) {
        K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
        diag(K) <- (z - 1)/z
        cosym <- K %*% coefmat[,!indcat]
        coef3 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% cosym
        coef2 <- matrix(0, ncol = lglob, nrow = z)
        coef2[refLevel, ] <- coef3
        coef2[-refLevel, ] <- cosym
        coefs <- exp(coef2)
            semat <- matrix(0, ncol = lglob, nrow = z)
            s1 <- 1
            for (e in 1:lglob) {
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
        coefs <- matrix(0, ncol = lglob, nrow = z)
        coefs[refLevel, ] <- rep(1, lglob)
        coefs[-refLevel, ] <- exp(coefmat[,!indcat])
            semat <- matrix(0, ncol = lglob, nrow = z)
            semat[refLevel, ] <- rep(0, lglob)
            semat[-refLevel, ] <- t(sdsmat)
    }

    x1 <- t(coefs)
    # compute confidence intervals
    if (printpvalues | conf.int) {
        x2 <- exp(log(x1) + qnorm(0.975) * t(semat))
        x3 <- exp(log(x1) - qnorm(0.975) * t(semat))
    }
    
    if(is.null(select)){
      lselect<-lglob
      select<-1:lglob
    }else{
      if(sum(select>lglob)>0){select<-select[select<=lglob]}
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
    plotmat2<-plotmat
    if(catstar){plotmat2 <-rbind(plotmat,plotmat[1,])}
    loc <- stars(plotmat2, nrow = nlines, scale = FALSE, len = 1, cex = 1.2, 
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
    # plot stars with or without confidence intervals
    
    if(catstar){loccat<-loc[nrow(loc),]
                loc<-loc[-nrow(loc),]}
    
    
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
            lwd = lwd.stars, lty = lty.conf, ylim = ylim, xlim = xlim, 
            labels = "", locations = loc, main = main, add = FALSE)
        par(cex = 1, font = 1)
    }
    if(catstar){
      catmat<-t(as.matrix(exp(coefmat[1,indcat])))
      faccat<-1
      if(scale){maxcat<-max(catmat)
                faccat<-maxcat/max(plotmat)
                catmat<-catmat/faccat}
    par(cex = 1, font = 3)
        stars(catmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loccat, 
            xlim = xlim, ylim = ylim, labels = "", add=TRUE)
      par(lty = lty.circle, lwd = lwd.circle, col = col.circle,cex = 1, font = 1)
    symbols(c(loccat[, 1]), c(loccat[, 2]), circles = 1/faccat, 
        add = T, inches = F, fg = col.circle, bg = col.fill)
    par(lty = 1, lwd = 1, col = 1)
        par(cex = 1, font = 3)
        stars(catmat, scale = FALSE, draw.segments = FALSE, len = 1, 
            lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loccat, 
            xlim = xlim, ylim = ylim, labels = "", add=TRUE)
      par(lty = 1, lwd = 1, col = 1)
    ##
    alts <- 0:(pcat - 1)
    angle <- alts * ((360/pcat) * pi)/180
    cosis <- cos(angle)
    sinis <- sin(angle)
        dir <- c(catmat)
        xcat <- loccat[1,1]
        ycat <- loccat[1,2]
        dis <- max(dir * 0.4 * dist.cat, dir * 0.4 * dist.cat)
        xlocscat <- rep(xcat, pcat) + cosis * dir + cosis * dis
        ylocscat <- rep(ycat, pcat) + sinis * dir + sinis * dis
      catlabels<-catnamesresults
      if(test.rel){catlabels<-paste(catnamesresults, paste("(", catpvals, ")", sep = ""), 
            sep = "\n")}
text(x = xlocscat, y = ylocscat, labels = catlabels, cex = cex.cat, font = 1)
text(y = ycat + fac * 1.9 * (l2 + 20)/32, x = xcat, 
        labels = "Global coefficients", font = 1, cex = cex.labels)     
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
       radii<-rep(1,lglob)
      }
    }
    else {
        radii <- x1[1:lglob, refLevel]
    }
    # plot cirlces
    if(!identical(select,1)){
    par(lty = lty.circle, lwd = lwd.circle, col = col.circle)
    symbols(c(loc[, 1]), c(loc[, 2]), circles = radii[select], 
        add = TRUE, inches = FALSE, fg = col.circle, bg = col.fill)
    par(lty = 1, lwd = 1, col = 1)
    }
    
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
        par(cex = 1, font = 1)
    }
    
    # print labels
    text(x = xlocs[,select], y = ylocs[,select], labels = starlab[c(xlocs)!=0], cex = cex.cat, 
        font = 1)
    
    

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

    cnames<-colnames(coefmat)[!indcat]

    if (printpvalues) {
        rownames(pvals) <- keylab
        colnames(pvals) <- cnames
    }
    rownames(semat) <- keylab
    colnames(semat) <-cnames
    rownames(odds) <- keylab
    colnames(odds) <- cnames
    catspec<-matrix(catspec,nrow=1)
    catspecses<-matrix(catspecses,nrow=1)
    colnames(catspec)<-colnames(catspecses)<-catnamesresults

# create result lists   
    if (test.rel) {
        p_rel <- matrix(c(catexact,lrexact), ncol = l - 1)
        colnames(p_rel) <- c(covcateg,colnames(coefmat)[!indcat][-1])
        if (printpvalues) {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, pvalues = pvals, 
               catspec = catspec, catspecse = catspecses, p_rel = p_rel,
               xlim = xlim, ylim = ylim)
        }
        else {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, catspec = catspec, 
                catspecse = catspecses, p_rel = p_rel, xlim = xlim, 
                ylim = ylim)
        }
    }
    else {
        if (!printpvalues) {
            returns <- list(odds = odds, 
                coefficients = log(odds), se = semat, catspec = catspec, 
                catspecse = catspecses, xlim = xlim, ylim = ylim)
        }
        else {
            returns <- list(odds = odds, coefficients = log(odds), 
                se = semat, pvalues = pvals, catspec = catspec, 
                catspecse = catspecses, xlim = xlim, ylim = ylim)
        }
    }
   
    returns
    ###########################################
}else{
#############################################
  isfac <- c()
  dummydat <-dummydat2<- resp
  outputlab<-c()
  dummynames<-response
  
  for (j in 1:p) {
    
    # dummy coding
    variab <- data[, colnames(data) == covariates[j]]
    levlen <- length(levels(variab))
    if (is.factor(variab) & levlen > 2) {
      outputlab<-c(outputlab,levels(variab))
      dummynames<-c(dummynames,paste(covariates[j],1:(levlen-1),sep=""))
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
      dummynames<-c(dummynames,covariates[j])
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
  colnames(dummydat) <- dummynames
  dummydat<-cbind(dummydat,data[,(colnames(data) %in% xijs) & !(colnames(data) %in% covcateg)])
  
  dummydat2 <- as.data.frame(dummydat2)
  colnames(dummydat2) <- dummynames
  dummydat2<-cbind(dummydat2,data[,(colnames(data) %in% xijs) & !(colnames(data) %in% covcateg)])
  
  multi.df<-rep(0,length(covariates))
  multi<-c()
  multi2<-c()
  vb<-1
  
  for(u in 1:length(covglob)){
    check.cat<-data[, colnames(data) == covglob[u]]
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
  covglob <-dummynames[!(dummynames %in% covcateg) & dummynames!=response]
  formglob <- paste(covglob,collapse="+")
  
  formdummy<-as.formula(paste(paste(response,"~",sep=""),paste(dummynames[dummynames!=response],sep="",collapse="+")))
  
  form2<-as.formula(paste("~",paste(unique(c(xijs,dummynames[dummynames!=response])),collapse="+"),sep=""))
  
  form2vars <- all.vars(form2)
  
  if(!(refLevel %in% 1:z))
    stop(paste("Argument refLevel does not contain a valid reference category from vector 1:",z,sep=""))
  
  ######################
  
  # compute model + likelihood + variances
  multobj <- vglm(formdummy, data = dummydat, family = multinomial(parallel = 
    as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                  xij = xij, form2 = form2, trace = FALSE)
  
  coefmat <- t(coefficients(multobj, T))
  
  l <- ncol(coefmat)
  indcat<-apply(coefmat, 2,var)==0
  indcat2 <- rep(indcat, each=z-1)
  lcat <- sum(indcat)
  lglob <- l - lcat
  indout <- rep(((1:length(indcat))[indcat]-1)*(z-1),each=z-2)+(2:(z-1))
  indcat2 <-indcat2[-indout]
  vcovtotal <- vcov(multobj)
  vcov <- vcov(multobj)[!indcat2,!indcat2]
  
  ##########################
  catspec<-coefmat[1,indcat]
  catnamesresults<-names(catspec)
  if(pcat==1){
    catspecses<-sqrt(vcov(multobj)[indcat2,indcat2])
  }else{
    catspecses <- sqrt(diag(vcov(multobj)[indcat2,indcat2]))}
  
  coefmat.orig<-coefmat
  coefmat<-coefmat[,!indcat]
  coef.after<-coefmat[,1]
  l.after<-lglob+vb-1
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
  
  lglob<-l.after
  l<-lglob+lcat
  vcov<-vcov.after
  coefmat<-coef.after
  
  #############################################
  sdsmat <- matrix(sqrt(diag(vcov)), nrow = lglob, byrow = TRUE)
  
  if (test.rel) {
    loglik <- logLik(multobj)
  }
  logxy <- c()
  isfac <- c()
  labvec <- ""
  covlab <- c()
  #loop: create labels, compute likelihoods for lr-tests
  
  for (u in 1:pglob) {
    actcov <- data[, colnames(data) == covglob2[u]]
    # create labels
    if (is.factor(actcov) | is.ordered(actcov)) {
      nlev <- length(levels(actcov))
      isfac <- c(isfac, rep(TRUE, nlev-1+I(nlev>2)))
      covlab <- c(covlab, rep(covglob2[u], nlev -1+I(nlev>2)))
      if (is.ordered(actcov)) {
        parlevs <- paste("(", sort(as.numeric(unique(actcov))), 
                         ": ", levels(actcov), sep = "")
        labvec <- c(labvec, parlevs)
      }
      else {
        if (length(levels(actcov)) == 2) {
          parlevs <- paste("(", sort(as.numeric(unique(actcov)) - 
            1)[-1], ": ", levels(actcov)[-1], sep = "")
        }
        else {
          parlevs <- paste("(", sort(as.numeric(unique(actcov))), 
                           ": ", levels(actcov), sep = "")
        }
        labvec <- as.expression(c(labvec, parlevs), font = 3)
      }
    }
    else {
      isfac <- c(isfac, FALSE)
      covlab <- c(covlab, covglob2[u])
      labvec <- c(labvec, "")
    }}
  # compute likelihoods
  
  if (test.rel) {
    for (u in 1:length(multi)) {
      
      if (length(multi) >= 2){
        formxy<-as.formula(paste(paste(response,"~",sep=""),paste(c(covcateg, (dummynames[dummynames!=response & dummynames %in% covglob])[-u]),sep="",collapse="+")))
      }else{
        formxy <- as.formula(paste(response, "~", 1))
      }
      
      log.add<-logLik(vglm(formxy, data = dummydat, family = multinomial(parallel = 
        as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                           xij = xij, form2 = form2, trace = FALSE))
      if(u ==length(multi) & multi[u]!=0){
        form.sym<-as.formula(paste(paste(response,"~",sep=""),paste(c(covcateg, (dummynames[dummynames!=response & dummynames %in% covglob])[-u]),sep="",collapse="+")))
        log.add<-c(log.add,logLik(vglm(form.sym, data = dummydat2, family = multinomial(parallel = 
          as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                                       xij = xij, form2 = form2, trace = FALSE)))
      }else{
        if(multi[u]!=0 & multi[u]!=multi[u+1]){
          form.sym<-as.formula(paste(paste(response,"~",sep=""),paste(c(covcateg, (dummynames[dummynames!=response & dummynames %in% covglob])[-u]),sep="",collapse="+")))
          log.add<-c(log.add,logLik(vglm(form.sym, data = dummydat2, family = multinomial(parallel = 
            as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                                         xij = xij, form2 = form2, trace = FALSE)))
        }}
      logxy <- c(logxy,log.add)
    }}
  
  
  
  ######################################### LR Tests for catspec covariates
  likecat<-c()  
  if (test.rel) {
    for(qe in 1:pcat){
      actcov <- data[, colnames(data) == covcateg[qe]]
      if (length(covcateg) >= 2){
        
        formxy<-as.formula(paste(paste(response,"~",sep=""),paste(c(covcateg[-qe], (dummynames[dummynames!=response & dummynames %in% covglob])),sep="",collapse="+")))
        #form22<-as.formula(paste("~",paste(form2vars[form2vars!=covcateg[qe]],collapse="+"),sep=""))
        xij2<-xij[-qe]
        mod2<-vglm(formxy, data = dummydat, family = multinomial(parallel = 
          as.formula(paste("FALSE ~", formglob,"")) , refLevel = refLevel), maxit = maxit,
                   xij = xij2, form2 = form2, trace = FALSE)
      }else{
        formxy <- as.formula(paste(response, "~", paste(covariates[covariates!=covcateg[qe]], 
                                                        collapse = " + ")))
        mod2<-vglm(formxy, data = data, family = multinomial(refLevel = refLevel), 
                   maxit = maxit)
      }
      
      likecat <- c(likecat, logLik(mod2))
    }
  }
  
  #########################################
  # compute lr-test p-values
  if (test.rel) {
    lrstats <- 2 * (-logxy + loglik)
    lrpvals <- 1 - pchisq(lrstats, df = z - 1)
    lrexact <- lrpvals
    lrpvals[lrpvals < 5e-04] <- 0
    lrpvals[lrpvals >= 5e-04 & lrpvals < 0.001] <- 0.001
    lrpvals <- format(lrpvals, digits = 1, nsmall = 3)
  }
  if (test.rel) {
    catstats <- 2 * (-likecat + loglik)
    catpvals <- 1 - pchisq(catstats, df = 1)
    catexact <- catpvals
    catpvals[catpvals < 5e-04] <- 0
    catpvals[catpvals >= 5e-04 & catpvals < 0.001] <- 0.001
    catpvals <- format(catpvals, digits = 1, nsmall = 3)
  }
  
  # compute coefficients and standard deviances for symmetric side constraints
  if (symmetric) {
    K <- matrix((-1/z), ncol = z - 1, nrow = z - 1)
    diag(K) <- (z - 1)/z
    cosym <- K %*% coefmat
    coef3 <- matrix(rep(-1, z - 1), ncol = z - 1) %*% cosym
    coef2 <- matrix(0, ncol = lglob, nrow = z)
    coef2[refLevel, ] <- coef3
    coef2[-refLevel, ] <- cosym
    coefs <- exp(coef2)
    semat <- matrix(0, ncol = lglob, nrow = z)
    s1 <- 1
    for (e in 1:lglob) {
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
    coefs <- matrix(0, ncol = lglob, nrow = z)
    coefs[refLevel, ] <- rep(1, lglob)
    coefs[-refLevel, ] <- exp(coefmat[,!indcat])
    semat <- matrix(0, ncol = lglob, nrow = z)
    semat[refLevel, ] <- rep(0, lglob)
    semat[-refLevel, ] <- t(sdsmat)
  }
  
  x1 <- t(coefs)
  # compute confidence intervals
  if (printpvalues | conf.int) {
    x2 <- exp(log(x1) + qnorm(0.975) * t(semat))
    x3 <- exp(log(x1) - qnorm(0.975) * t(semat))
  }
  
  if(is.null(select)){
    lselect<-lglob
    select<-1:lglob
  }else{
    if(sum(select>lglob)>0){select<-select[select<=lglob]}
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
  plotmat2<-plotmat
  if(catstar){plotmat2 <-rbind(plotmat,plotmat[1,])}
  loc <- stars(plotmat2, nrow = nlines, scale = FALSE, len = 1, cex = 1.2, 
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
  # plot stars with or without confidence intervals
  
  if(catstar){loccat<-loc[nrow(loc),]
              loc<-loc[-nrow(loc),]}
  
  
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
          lwd = lwd.stars, lty = lty.conf, ylim = ylim, xlim = xlim, 
          labels = "", locations = loc, main = main, add = FALSE)
    par(cex = 1, font = 1)
  }
  if(catstar){
    catmat<-t(as.matrix(exp(coefmat.orig[1,indcat])))
    faccat<-1
    if(scale){maxcat<-max(catmat)
              faccat<-maxcat/max(plotmat)
              catmat<-catmat/faccat}
    par(cex = 1, font = 3)
    stars(catmat, scale = FALSE, draw.segments = FALSE, len = 1, 
          lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loccat, 
          xlim = xlim, ylim = ylim, labels = "", add=TRUE)
    par(lty = lty.circle, lwd = lwd.circle, col = col.circle,cex = 1, font = 1)
    symbols(c(loccat[, 1]), c(loccat[, 2]), circles = 1/faccat, 
            add = T, inches = F, fg = col.circle, bg = col.fill)
    par(lty = 1, lwd = 1, col = 1)
    par(cex = 1, font = 3)
    stars(catmat, scale = FALSE, draw.segments = FALSE, len = 1, 
          lwd = lwd.stars, lty = "solid", flip.labels = FALSE, locations = loccat, 
          xlim = xlim, ylim = ylim, labels = "", add=TRUE)
    par(lty = 1, lwd = 1, col = 1)
    ##
    alts <- 0:(pcat - 1)
    angle <- alts * ((360/pcat) * pi)/180
    cosis <- cos(angle)
    sinis <- sin(angle)
    dir <- c(catmat)
    xcat <- loccat[1,1]
    ycat <- loccat[1,2]
    dis <- max(dir * 0.4 * dist.cat, dir * 0.4 * dist.cat)
    xlocscat <- rep(xcat, pcat) + cosis * dir + cosis * dis
    ylocscat <- rep(ycat, pcat) + sinis * dir + sinis * dis
    catlabels<-catnamesresults
    if(test.rel){catlabels<-paste(catnamesresults, paste("(", catpvals, ")", sep = ""), 
                                  sep = "\n")}
    text(x = xlocscat, y = ylocscat, labels = catlabels, cex = cex.cat, font = 1)
    text(y = ycat + fac * 1.9 * (l2 + 20)/32, x = xcat, 
         labels = "Global coefficients", font = 1, cex = cex.labels)     
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
      radii<-rep(1,lglob)
    }
  }
  else {
    radii <- x1[1:lglob, refLevel]
  }
  # plot cirlces
  if(!identical(select,1)){
    par(lty = lty.circle, lwd = lwd.circle, col = col.circle)
    symbols(c(loc[, 1]), c(loc[, 2]), circles = radii[select], 
            add = TRUE, inches = FALSE, fg = col.circle, bg = col.fill)
    par(lty = 1, lwd = 1, col = 1)
  }
  
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
    par(cex = 1, font = 1)
  }
  
  # print labels
  text(x = xlocs[,select], y = ylocs[,select], labels = starlab[c(xlocs)!=0], cex = cex.cat, 
       font = 1)
  
  
  
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
  cnames<-covlab
  
  if (printpvalues) {
    rownames(pvals) <- keylab
    colnames(pvals) <- cnames
  }
  rownames(semat) <- keylab
  colnames(semat) <-cnames
  rownames(odds) <- keylab
  colnames(odds) <- cnames
  catspec<-matrix(catspec,nrow=1)
  catspecses<-matrix(catspecses,nrow=1)
  colnames(catspec)<-colnames(catspecses)<-catnamesresults
  
  # create result lists   
  if (test.rel) {
    p_rel <- matrix(c(catexact,lrexact), ncol = l - 1)
    colnames(p_rel) <- unlist(c(covcateg,cnames[-1]))
    if (printpvalues) {
      returns <- list(odds = odds, 
                      coefficients = log(odds), se = semat, pvalues = pvals, 
                      catspec = catspec, catspecse = catspecses, p_rel = p_rel,
                      xlim = xlim, ylim = ylim)
    }
    else {
      returns <- list(odds = odds, 
                      coefficients = log(odds), se = semat, catspec = catspec, 
                      catspecse = catspecses, p_rel = p_rel, xlim = xlim, 
                      ylim = ylim)
    }
  }
  else {
    if (!printpvalues) {
      returns <- list(odds = odds, 
                      coefficients = log(odds), se = semat, catspec = catspec, 
                      catspecse = catspecses, xlim = xlim, ylim = ylim)
    }
    else {
      returns <- list(odds = odds, coefficients = log(odds), 
                      se = semat, pvalues = pvals, catspec = catspec, 
                      catspecse = catspecses, xlim = xlim, ylim = ylim)
    }
  }
  
  returns
}
}