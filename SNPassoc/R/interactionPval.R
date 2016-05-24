`interactionPval` <-
function (formula, data, quantitative = is.quantitative(formula, 
    data), model="codominant") 
  {
 
    if(!inherits(data,"setupSNP"))
     stop("data must be an object of class 'setupSNP'")

    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m0 <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m0)]

#
# aceptar respuesta sin formula
#
	if( length(grep("~",mf[[2]]))==0){
		formula<-as.formula(paste(mf[[2]],"~1",sep=""))
		formula.1<- list(formula)
		mode(formula.1)<-"call"
		mf[2]<-formula.1
	}

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    temp0 <- as.character(mt)
    adj <- paste(temp0[2], temp0[1], temp0[3])

    fam <- ifelse(quantitative,"gaussian","binomial")   

    model.type <- c("codominant", "dominant", "recessive", 
                "overdominant","log-additive")
    m <- charmatch(model, model.type, nomatch = 0)
    if (m == 0) 
          stop("model must be codominant dominant recessive overdominant or log-additive")

    modelOK<-switch(m,codominant,dominant,recessive,overdominant,additive)


    colSNPs<-attr(data,"colSNPs")
    if (is.vector(colSNPs) & length(colSNPs) > 0) 
        dataSNPs.sel <- data[, colSNPs, drop=FALSE]
    else stop("data should have an attribute called 'colSNPs'. Try again 'setupsNP' function")
    
    dataSNPs <- data.frame(lapply(dataSNPs.sel,function(x,model.sel) model.sel(x),model.sel=modelOK)) 

    SNPs.label <- names(dataSNPs)

    dimnames(data)[[1]]<-1:nrow(data)

    i<-1
    n<-ncol(dataSNPs)

    pval<-matrix(NA,nrow=n,ncol=n)


    while (i<=n)
     {
         nas <- sum(!is.na(dataSNPs[, i]))
         n.nas <- length(dataSNPs[, i])

         if (is.Monomorphic(dataSNPs[, i]))
          {
            pval[i,]<-rep(NA,n)
          }

         else if (nas/n.nas<0.80)
          {
            pval[i,]<-rep(NA,n)
          }
         
         else if (length(table(dataSNPs[, i]))==1)
          {
            pval[i,]<-rep(NA,n)
          }

         else
         {
          j<-i+1
          while (j<=n)
           {
            nas <- sum(!is.na(dataSNPs[, j]))
            n.nas <- length(dataSNPs[, j])

            if (is.Monomorphic(dataSNPs[, j]))
             {
                pval[i,j]<-NA
             }
            else if (nas/n.nas<0.80)
             {
              pval[i,j]<-NA
             }

            else if (length(table(dataSNPs[, j]))==1)
             {
              pval[i,j]<-NA
             }


            else
             {
              mod.i <- glm(as.formula(paste(adj, "+ dataSNPs[, i]*dataSNPs[, j]")),
                        data = data, family=fam)
              subset <- 1:nrow(data) %in% as.numeric(rownames(mod.i$model))
              mod.a <- glm(as.formula(paste(adj, "+ dataSNPs[, i]+dataSNPs[, j]")),
                         data = data, family=fam, subset=subset)

              mod.b1 <- glm(as.formula(paste(adj, "+ dataSNPs[, i]")), data = data,
                         family=fam,subset=subset)
              mod.b2 <- glm(as.formula(paste(adj, "+ dataSNPs[, j]")), data = data,
                         family=fam,subset=subset)

              if (quantitative)
               pval[i,j]<-anova(mod.a,mod.i,test="F")$"Pr(>F)"[2]
              else
               {
                 t1 <- anova(mod.a, mod.i, test="Chisq")
                 pval[i,j] <- t1[2, grep("^P.*Chi",names(t1))]
               }

              if(mod.b1$aic<=mod.b2$aic) 
               {
                if (quantitative)
                  pval[j,i] <- anova(mod.b1, mod.a, test="F")$"Pr(>F)"[2]
                else
                 {
                  t1 <- anova(mod.b1, mod.a, test="Chisq")
                  pval[j,i] <- t1[2, grep("^P.*Chi",names(t1))]
                 }
               }
              else  
               {
                if (quantitative)
                  pval[j,i] <- anova(mod.b2, mod.a, test="F")$"Pr(>F)"[2]
                else
                 {
                  t1 <- anova(mod.b2, mod.a, test="Chisq")
                  pval[j,i] <- t1[2, grep("^P.*Chi",names(t1))]
                 }
               }
             }
             j<-j+1
           }

          mod.0 <- glm(as.formula(paste(adj, "+ dataSNPs[, i]")), data = data,
                         family="gaussian")
          subset <- 1:nrow(data) %in% as.numeric(rownames(mod.0$model))
          mod.b <- glm(as.formula(paste(adj)), data = data,
                         family="gaussian",subset=subset)
          pval[i,i] <- anova(mod.b,mod.0,test="F")$"Pr(>F)"[2]
         }

          i<-i+1
       }

 dimnames(pval)[[2]]<-SNPs.label
 dimnames(pval)[[1]]<-SNPs.label

 class(pval)<-"SNPinteraction" 
 attr(pval,"model") <- model.type[m]
 attr(pval,"gen.info")<-attr(data,"gen.info")
 pval
}

