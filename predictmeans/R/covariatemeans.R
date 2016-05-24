covariatemeans <- function (model, modelterm=NULL, covariate, level=0.05, Df=NULL, trans=NULL, 
                            responsen=NULL, trillis=TRUE, plotord=NULL, mtitle=NULL, ci=TRUE, point=TRUE, jitterv=0, newwd=TRUE) {
  
  if (is.null(modelterm) || modelterm%in%c("NULL", "")) {
    modelterm <- covariate
    trillis=FALSE
  }
  
  vars <- unlist(strsplit(modelterm, "\\:"))
  ctr.matrix <- Kmatrix(model, modelterm, covariate)
  KK <- ctr.matrix$K
  pltdf <- ctr.matrix$fctnames
  response <- ctr.matrix$response
  preddf <- ctr.matrix$preddf
  mp <- mymodelparm(model)
  bhat <- mp$coef
  
  # We'll work only with the non-NA elements of bhat
  KK <- KK[, mp$estimable, drop=FALSE]   
  pltdf$yhat <- KK%*%bhat
  pltdf$ses <- sqrt(base::diag(KK %*% tcrossprod(mp$vcov, KK)))
  
  if (is.null(Df) || Df%in%c("NULL", "")) {
    if (class(model)[1] == "lme") {
      Df <- terms(model$fixDF)[modelterm]
    }else if (class(model)[1] == "lmerMod") {
      termlabel <- attr(terms(model),"term.labels")
      for (i in vars) termlabel <- termlabel[grep(i, termlabel)]
      termlabel <- paste(termlabel, collapse="-")
      model.b <- update( model, as.formula(paste(".~. -", termlabel)))
      Df <- getKR(KRmodcomp(model, model.b), "ddf")
    }else Df <- mp$df
    if (Df==0) stop("You need provide Df for this model!")
  }
  
  if (class(model)[1]=="glm" && is.null(trans)) trans <- model$family$linkinv
  if (class(model)[1] == "glmerMod" && is.null(trans)) trans <- slot(model, "resp")$family$linkinv
  
  Mean <- LL <- UL <- xvar <- factors <- bky <- NULL
  if (is.null(trans)) {
    pltdf$Mean <- pltdf$yhat
    pltdf$LL <- pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses
    pltdf$UL <- pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses
  }else{
    pltdf$Mean <- trans(pltdf$yhat)
    pltdf$LL <- trans(pltdf$yhat - qt(1 - level/2, df = Df) * pltdf$ses)
    pltdf$UL <- trans(pltdf$yhat + qt(1 - level/2, df = Df) * pltdf$ses)
  }
  pltdf$yhat <- pltdf$ses <- NULL 
  if (modelterm==covariate) pltdf$factors <- factor(1) else pltdf$factors <- factor(do.call("paste", c(pltdf[, vars, drop=FALSE], sep=":")))
  colnames(pltdf)[colnames(pltdf)==covariate] <- "xvar"
  
  # delete empty factor combinations  
  mdf <- model.frame(model)
  if (!(response %in% names(mdf))) mdf[,response] <- eval(parse(text=response), mdf)
  mdf <- cbind(mdf, preddf[, !names(preddf)%in%names(mdf), drop=FALSE])
  
  ndf <- data.frame(table(mdf[, vars, drop = FALSE]))
  if (any(ndf$Freq==0)) { 
    ndf0 <- ndf[ndf$Freq==0, , drop=FALSE] 
    ndf0$factors <- factor(do.call("paste", c(ndf0[, vars, drop=FALSE], sep=":")))
    pltdf <- pltdf[!pltdf$factors%in%ndf0$factors, ]
  }   
  
  
  if (is.null(mtitle) || mtitle%in%c("NULL", "")) mtitle <- paste("Fitted and observed relationship with", (1-level)*100, "% CI")
  if (is.null(trans)) {
    mdf$bky <- mdf[, response]
  }else{
    if (response %in% names(mdf)) {    ## Transformed y before modelling
      if (class(model)[1]%in%c("glm", "glmerMod")) {
        if (class(mdf[, response])=="factor") {
          mdf$bky <- as.numeric(mdf[, response])-1
        }else if (!is.null(dim(mdf[, response]))) {
          mdf$bky <- mdf[, response][,1]/rowSums(mdf[, response])
          # (is.null(responsen) || responsen%in%c("NULL", "")) stop("Please provide suitable name for response variable using option 'responsen'!")
          response <- "Probability"
        }else mdf$bky <- mdf[, response]
        if (isTRUE(all.equal(trans,function(x) x))) {
          if (class(model)[1]=="glm") mdf$bky <- model$family$linkfun(mdf$bky)
          if (class(model)[1] == "glmerMod") mdf$bky <- slot(model, "resp")$family$linkfun(mdf$bky)
          # f (is.null(responsen) || responsen%in%c("NULL", "")) stop("Please provide suitable name for response variable using option 'responsen'!")
          response <- "Response"
        }
      }else{
        mdf$bky <- trans(mdf[, response])
        # if (is.null(responsen) || responsen%in%c("NULL", ""))  stop("Please provide suitable name for response variable using option 'responsen'!")
        response <- paste("Transformed", response)
      }
    }else{       ## Transformed y within modelling
      response <- regmatches(response, regexec("\\(([^<]+)\\)", response))[[1]][2]
      if (!response %in% names(mdf)) {
        if (is.null(responsen) || responsen%in%c("NULL", ""))  stop("Please provide suitable name for response variable using option 'responsen'!")
        response <- responsen
      }
      mdf$bky <- mdf[, response]
    }
  }
  if (modelterm==covariate) mdf$factors <- factor(1) else mdf$factors <- do.call("paste", c(mdf[, vars, drop=FALSE], sep=":"))
  names(mdf)[names(mdf)==covariate] <- "xvar"
  if (!trillis) {
    if (newwd) dev.new()
    if (modelterm==covariate)  {
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), geom="line", ylab=paste(response, "\n"), data=pltdf, main=paste(mtitle, "\n")) +
        theme_bw()
      if (ci) plt <- plt + geom_smooth(aes(ymin = LL, ymax = UL), alpha = 0.2, data=pltdf, stat="identity")
      if (point) plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
    }else{
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), geom="line", ylab=paste(response, "\n"), data=pltdf, main=paste(mtitle, "\n"), colour=factors) +
        theme_bw()
      if (ci) plt <- plt + geom_smooth(aes(ymin = LL, ymax = UL, fill=factors), alpha = 0.2, data=pltdf, stat="identity")
      if (point) plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      plt <- plt+guides(col = guide_legend(modelterm), fill=guide_legend(modelterm))
    }
    print(plt)
  }else{
    if (length(vars)==1) {
      if (newwd) dev.new()
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), geom="line", ylab=paste(response, "\n"), data=pltdf, main=paste(mtitle, "\n"), colour=factors) +
          facet_wrap(~  factors)+
		  theme_bw()
      if (ci) plt <- plt + geom_smooth(aes(ymin = LL, ymax = UL, fill=factors), alpha = 0.2, data=pltdf, stat="identity")
      if (point) plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      plt <- plt+guides(col = guide_legend(modelterm), fill=guide_legend(modelterm))
      if (modelterm==covariate)  plt <- plt+ theme(legend.position="none")
      print(plt)
    }
    if (length(vars)==2) {
      if (newwd) dev.new()
      if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:2
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), ylab=paste(response, "\n"), main=paste(mtitle, "\n"), data=pltdf,
                   geom="line", colour=factor(eval(parse(text = fact1)))) +
        facet_grid(eval(parse(text = paste("~",fact2, sep=""))))+
        theme_bw()
      if (ci) plt <- plt + geom_smooth(aes(ymin = LL, ymax = UL, fill=factor(eval(parse(text = fact1)))), alpha = 0.2, data=pltdf, stat="identity")
      if (point) plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      plt <- plt+guides(col = guide_legend(fact1), fill=guide_legend(fact1))
      print(plt)
    }
    if (length(vars)==3) {
      if (newwd) dev.new()
      if (is.null(plotord) || plotord%in%c("NULL", "")) plotord <- 1:3
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      fact3 <- (vars[plotord])[3]
      plt <- qplot(xvar, Mean, xlab=paste("\n", covariate, sep=""), ylab=paste(response, "\n"), main=paste(mtitle, "\n"), data=pltdf,
                   geom="line", colour=factor(eval(parse(text = fact1)))) +
        facet_grid(eval(parse(text = paste(fact2, "~",fact3, sep=""))))+
        theme_bw()
      if (ci) plt <- plt + geom_smooth(aes(ymin = LL, ymax = UL, fill=factor(eval(parse(text = fact1)))), alpha = 0.2, data=pltdf, stat="identity")
      if (point) plt <- plt + geom_point(aes(x=xvar, y=bky), position = position_jitter(width = jitterv, height = jitterv), data=mdf)
      plt <- plt+guides(col = guide_legend(fact1), fill=guide_legend(fact1))
      print(plt)
    }
  }
  
  return(invisible(plt))
}
