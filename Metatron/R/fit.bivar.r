fit.bivar<-function(TP,FN,TN,FP,study, data, mods=NULL,covarying="both",...) {

    na.act <- getOption("na.action")

    if (missing(data)) 
        data <- NULL
    no.data <- is.null(data)
    if (is.null(data)) {
        data <- sys.frame(sys.parent())
    }
    else {
        if (!is.data.frame(data)) {
            data <- data.frame(data)
        }
    }

   mf <- match.call()
   
    mf.TP <- mf[[match("TP", names(mf))]]
    mf.FN <- mf[[match("FN", names(mf))]]
    mf.TN <- mf[[match("TN", names(mf))]]
    mf.FP <- mf[[match("FP", names(mf))]]
    mf.study <- mf[[match("study", names(mf))]]
    mf.mods <- mf[[match("mods",   names(mf))]]

    TP <- eval(mf.TP, data, enclos = sys.frame(sys.parent()))
    FN <- eval(mf.FN, data, enclos = sys.frame(sys.parent()))
    TN <- eval(mf.TN, data, enclos = sys.frame(sys.parent()))
    FP <- eval(mf.FP, data, enclos = sys.frame(sys.parent()))
    study <- eval(mf.study, data, enclos = sys.frame(sys.parent()))
    mods   <- eval(mf.mods,   data, enclos=sys.frame(sys.parent()))

       obsnum<-length(study)
       studynum<-length(levels(study))

    groupS<-data.frame(study, wellclassified=TP, malclassified=FN, group="disease" )

    groupH<-data.frame(study, wellclassified=TN, malclassified=FP, group="healthy"  )

    X_long<-rbind(groupS,groupH)

    bi.simple<-NULL
    bi.sen<-NULL
    bi.spe<-NULL
    bi.both<-NULL
    
   if (is.null(mods)) {suppressWarnings(bi.simple<-glmer(cbind(wellclassified,malclassified)~group-1+(group-1|study),family=binomial,data=X_long))
   }
      else {

    groupS<-data.frame(study, wellclassified=TP, malclassified=FN, group="disease", mods=mods)
    groupH<-data.frame(study, wellclassified=TN, malclassified=FP, group="healthy", mods=mods)
    X_long2<-rbind(groupS,groupH)

      X_long2$grp_d <- ifelse(X_long2$group == "disease", 1, 0)
      X_long2$grp_h <- ifelse(X_long2$group == "healthy", 1, 0)

      if (covarying=="both") suppressWarnings(bi.both<-glmer(cbind(wellclassified,malclassified)~mods:grp_d + mods:grp_h-1+(grp_d + grp_h -1|study),family=binomial,data=X_long2))
      else if (covarying=="only sensitivity") suppressWarnings(bi.sen<-glmer(cbind(wellclassified,malclassified)~grp_h + mods:grp_d-1+(grp_d + grp_h -1|study),family=binomial,data=X_long2))
      else if (covarying=="only specificity") suppressWarnings(bi.spe<-glmer(cbind(wellclassified,malclassified)~grp_d + mods:grp_h-1+(grp_d + grp_h -1|study),family=binomial,data=X_long2))
      else stop("Need to specify the influence of moderator")

   }
   
   if (!is.null(bi.simple)) {BIC <- BIC(bi.simple)
                             AIC <- AIC(bi.simple)
                        deviance <- deviance(bi.simple)
                          logLik <- logLik(bi.simple)
                           fixef <- fixef(bi.simple)
                           ranef <- ranef(bi.simple)
                           resid <- resid(bi.simple)
                         VarCorr <- VarCorr(bi.simple)
                           anova <- anova(bi.simple)
                            coef <- coef(bi.simple)
                            vcov <- vcov(bi.simple) 
                        levelmod <- NULL   
                         rancoef <- summary(bi.simple)$varcor$study                     
                      }

      else if (!is.null(bi.sen)){BIC <- BIC(bi.sen)
                             AIC <- AIC(bi.sen)
                        deviance <- deviance(bi.sen)
                          logLik <- logLik(bi.sen)
                           fixef <- fixef(bi.sen)
                           ranef <- ranef(bi.sen)
                           resid <- resid(bi.sen)
                         VarCorr <- VarCorr(bi.sen)
                           anova <- anova(bi.sen)
                            coef <- coef(bi.sen)
                            vcov <- vcov(bi.sen)                           
                        levelmod <- levels(bi.sen@frame$mods)
                         rancoef <- summary(bi.sen)$varcor$study
                    }

       else if (!is.null(bi.spe)){
                             BIC <- BIC(bi.spe)
                             AIC <- AIC(bi.spe)
                        deviance <- deviance(bi.spe)
                          logLik <- logLik(bi.spe)
                           fixef <- fixef(bi.spe)
                           ranef <- ranef(bi.spe)
                           resid <- resid(bi.spe)
                         VarCorr <- VarCorr(bi.spe)
                           anova <- anova(bi.spe)
                            coef <- coef(bi.spe)
                            vcov <- vcov(bi.spe) 
                        levelmod <- levels(bi.spe@frame$mods)                         
                         rancoef <- summary(bi.spe)$varcor$study
                    }

     else if (!is.null(bi.both)){
                             BIC <- BIC(bi.both)
                             AIC <- AIC(bi.both)
                        deviance <- deviance(bi.both)
                          logLik <- logLik(bi.both)
                           fixef <- fixef(bi.both)
                           ranef <- ranef(bi.both)
                           resid <- resid(bi.both)
                         VarCorr <- VarCorr(bi.both)
                           anova <- anova(bi.both)
                            coef <- coef(bi.both)
                            vcov <- vcov(bi.both)  
                        levelmod <- levels(bi.both@frame$mods)                        
                         rancoef <- summary(bi.both)$varcor$study
                    } 
      
      df.fixef<-length(fixef)
     
     
     output<-NULL
      output$covarying=covarying
      output$mods<-mods
      output$mf.mods<-mf.mods
      output$studynum<-studynum
      output$obsnum<-obsnum
      output$bi.simple<-bi.simple
      output$bi.sen<-bi.sen
      output$bi.spe<-bi.spe
      output$bi.both<-bi.both
      output$levelmod<-levelmod
      output$BIC<-BIC  
      output$AIC<-AIC
      output$deviance<-deviance
      output$logLik<-logLik
      output$fixef<-fixef
      output$df.fixef<-df.fixef
      output$rancoef<-rancoef
      output$ranef<-ranef
      output$resid<-resid
      output$VarCorr<-VarCorr
      output$anova<-anova
      output$coef<-coef
      output$vcov<-vcov
      output$X_long<-X_long
      
         class(output)<-c("fit.bivar","list")
          output
}

