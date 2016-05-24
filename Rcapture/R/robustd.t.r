"robustd.t" <- function(X, dfreq=FALSE, vt, vm="M0", vh=list("Chao"), vtheta=2, neg=TRUE)
{
########################################################################################
# Validation des arguments fournis en entrée et changement de leur forme si nécessaire
########################################################################################

        valid.one(dfreq,"logical")
        valid.vt(vt)
        Xvalid<-valid.X(X=X, dfreq=dfreq, vt=vt)
            X <- Xvalid$X
            I <- length(vt)  ## nombre de periodes primaires
        vm <- valid.vm(vm,c("none","M0","Mt","Mh","Mth"),vt,typet=TRUE)      
        vh <- valid.vh(vh,c("Chao","Poisson","Darroch","Gamma"),vm)
        vtheta <- valid.vtheta(vtheta,vh)
        valid.one(neg,"logical")
    
##########################################################################################
# AJUSTEMENT DU MODÈLE
##########################################################################################

#-------------------------------------#
# Élaboration et ajustement du modèle #
#-------------------------------------#

        Y <- histfreq.t(X,dfreq=dfreq)
        fct.call <- match.call()
        histpos <- histpos.t(sum(na.rm=TRUE,vt))
        Xw <- Xomega(vt,vm,vh,vtheta,fct.call,typet=TRUE,histpos)     # deuxieme composante (celle de Beta) dans le modele loglineaire du robust design
        Xdelta <- matrix(0,dim(histpos)[1],I)
        for (i in 1:I)
        {
                if (i==1) { Xs <- histpos[,c(1:vt[i])] } else
                { Xs <- histpos[,c((sum(na.rm=TRUE,vt[1:(i-1)])+1):sum(na.rm=TRUE,vt[1:i]))] }
                Xdelta[,i] <- apply(Xs,1,max)
        }
        Zw <- Zdelta(Xdelta)     # premiere composante (celle de Alpha) dans ce meme modele
        mX. <- cbind(Zw,Xw$mat)      # on fusionne ces 2 composantes explicatives
        gammanames <- paste("gamma",1:(2*I-2),sep="")
        colnames(mX.) <- c(gammanames,Xw$coeffnames)
        dimX <- dim(mX.)[2]

        # Ajustement du modèle
        anaMrd <- suppressWarnings(glm(Y~mX.,family=poisson))

    #Vérification du bon ajustement du modèle loglinéaire
    if(!anaMrd$converged) stop("'glm' did not converge")
    if(any(is.na(anaMrd$coef))) warning("some loglinear parameter estimations cannot be evaluated")

#-----------------------------------------------------------------------------------#
# Réajustement du modèle en enlevant les gamma et les eta (si modèle Chao) négatifs #
#-----------------------------------------------------------------------------------#

        param<-anaMrd$coef
        ppositions <- 0
        if (neg)
        {
            # Vecteur d'indicatrices pour les paramètres d'intérêt négatifs
            indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0)))
            for (i in 1:I)
            {
                test.Chao <- if(is.function(vh[[i]])||is.null(vh[[i]])) FALSE else if (vh[[i]]=="Chao") TRUE else FALSE
                if(test.Chao) {
                    indic <- as.vector(c(indic,rep(0,Xw$nbcoeff[i]-vt[i]+2),ifelse(param[(2*I+sum(na.rm=TRUE,Xw$nbcoeff[1:i])-vt[i]+2):(2*I+sum(na.rm=TRUE,Xw$nbcoeff[1:i])-1)]<0,1,0)))
                } else {
                    indic <- as.vector(c(indic,rep(0,Xw$nbcoeff[i])))
                }
            }
            while(sum(na.rm=TRUE,indic)>0) # Répéter la boucle jusqu'à ce qu'aucun gamma approprié ne soit négatif
            {
                # Détermination de la position du premier gamma approprié négatif
                pos <- 1
                while(indic[pos]==0) pos <- pos + 1
                ppositions <- c(ppositions,pos)
                # Retrait de la bonne colonne de mX. et réajustement du modèle
                mX. <- mX.[,-(pos-sum(na.rm=TRUE,ppositions<pos))]
                anaMrd <- suppressWarnings(glm(Y~mX.,family=poisson))
                # Ajout de zéros dans le vecteur des paramètres loglinéaires
                positions <- sort(ppositions[-1])                
                param <- rep(0,dimX+1)
                param[-positions] <- anaMrd$coef 
                # Vecteur d'indicatrices pour les paramètres d'intérêt négatifs
                indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0)))
                for (i in 1:I)
                {
                    test.Chao <- if(is.function(vh[[i]])||is.null(vh[[i]])) FALSE else if (vh[[i]]=="Chao") TRUE else FALSE
                    if (test.Chao) {
                        indic <- as.vector(c(indic,rep(0,Xw$nbcoeff[i]-vt[i]+2),ifelse(param[(2*I+sum(na.rm=TRUE,Xw$nbcoeff[1:i])-vt[i]+2):(2*I+sum(na.rm=TRUE,Xw$nbcoeff[1:i])-1)]<0,1,0)))
                    } else {
                        indic <- as.vector(c(indic,rep(0,Xw$nbcoeff[i])))
                    }
                }
            }
        }
        positions <- sort(ppositions[-1]) 

#------------------------------------------------------------------------#
# Ajustement d'un modèle pour tester la présence d'émigration temporaire #
#------------------------------------------------------------------------#

        Inono <- length(vm[vm!="none"])
        idpemig <- NULL
        for (i in 2:(I-1)) if(vm[i]!="none") idpemig <- c(idpemig,i)
        
        if(Inono==3) {
            anaMrd2 <- NULL
            parap2 <- NULL

            mX3<-cbind(mX.,Xdelta[,idpemig])
            anaMrd3 <- suppressWarnings(glm(Y~mX3,family=poisson))
            parap3<-summary(anaMrd3)$coef[length(anaMrd3$coef),1:2]
        } else if(Inono>3) {
           mX2<-cbind(mX.,apply(Xdelta[,idpemig],1,sum))
           anaMrd2 <- suppressWarnings(glm(Y~mX2,family=poisson))
           parap2<-summary(anaMrd2)$coef[length(anaMrd2$coef),1:2]

           mX3<-cbind(mX.,Xdelta[,idpemig])
           anaMrd3 <- suppressWarnings(glm(Y~mX3,family=poisson))
           nn<-dim(summary(anaMrd3)$coef)[1]
           parap3<-summary(anaMrd3)$coef[(nn-Inono+3):nn,1:2]
        } else {
           anaMrd2 <- NULL
           parap2 <- NULL

           anaMrd3 <- NULL
           parap3 <- NULL
        }
     
#---------------------------------------#
# Formation des vecteurs des paramètres #
#---------------------------------------#

        # valeur de l intercept
        interc <- param[1]


        # creation des vecteurs de parametres alpha et beta
        Alpha <-rep(0,2*I-2)
        for (i in (1:(2*I-2)))
        {
                Alpha[i] <- param[i+1]
        }

        Beta <- rep(0,length(Xw$mat[1,]))
        for (i in (1:length(Xw$mat[1,])))
        {
                Beta[i] <- param[2*I-1+i]
        }


        # Vérification de la présence de paramètres gamma négatifs si l'option "neg"=FALSE
        if(!neg)
        {
            if (any(Alpha<0)) warning("one or more gamma parameters are negative,\n",
                                      "you can set them to zero with the 'neg' option.")
        }

#--------------------------------------------------------------#
# Matrice de variances-covariances des paramètres loglinéaires #
#--------------------------------------------------------------#

        if(length(positions)>0) {
             varcov <- matrix(0,dimX+1,dimX+1)
             varcov[-positions,-positions] <- summary(anaMrd)$cov.unscaled
        } else varcov <- summary(anaMrd)$cov.unscaled  

########################################################################################
# Estimation des paramètres démographiques
########################################################################################

#--------------------------------------------#
# calcul des probabilites de capture (pstar) #
#--------------------------------------------#

        pstar <- rep(0,I)
        dpstar<-matrix(rep(0,I*length(param)),ncol=I)
        beta <- Beta
        for (i in (1:I))
        {
                if (vm[i]=="none") # cas du no model
                {
                        pstar[i]<- exp(beta[1]+log(2^vt[i]-1))/(1+exp(beta[1]+log(2^vt[i]-1)))
                        dpstar[(2*I + if(i>1) sum(na.rm=TRUE,Xw$nbcoeff[1:(i-1)]) else 0):(2*I-1+sum(na.rm=TRUE,Xw$nbcoeff[1:i])),i] <- exp(beta[1]+log(2^vt[i]-1))/(1+exp(beta[1]+log(2^vt[i]-1)))^2
                        beta <- beta[-1]
                } else if (vm[i]=="Mt") # cas du Modele Mt
                {
                        pstar[i] <- 1-prod((1+exp(beta[c(1:vt[i])]))^-1)
                        dpstar[(2*I + if(i>1) sum(na.rm=TRUE,Xw$nbcoeff[1:(i-1)]) else 0):(2*I-1+sum(na.rm=TRUE,Xw$nbcoeff[1:i])),i] <- prod((1+exp(beta[c(1:vt[i])]))^-1)*exp(beta[c(1:vt[i])])/(1+exp(beta[c(1:vt[i])]))
                        beta <- beta[-c(1:vt[i])]
                } else # Tous les autres modèles
                {
                        Xpf <- Xclosedp(vt[i],vm[i],vh[[i]],vtheta[i])
                        pstar[i] <- sum(na.rm=TRUE,exp(Xpf$mat%*%beta[c(1:Xpf$nbcoeff)]))/(1+ sum(na.rm=TRUE,exp(Xpf$mat%*%beta[c(1:Xpf$nbcoeff)])))
                        dpstar[(2*I + if(i>1) sum(na.rm=TRUE,Xw$nbcoeff[1:(i-1)]) else 0):(2*I-1+sum(na.rm=TRUE,Xw$nbcoeff[1:i])),i] <- t(Xpf$mat)%*%exp(Xpf$mat%*%beta[c(1:Xpf$nbcoeff)])/(1+sum(na.rm=TRUE,exp(Xpf$mat%*%beta[c(1:Xpf$nbcoeff)])))^2
                        beta <- beta[-c(1:Xpf$nbcoeff)]
                }
        }
        varcovpstar <- t(dpstar)%*%varcov%*%dpstar
        pstarStderr <- sqrt(diag(varcovpstar))    

#---------------#
# calcul des Ui #
#---------------#

        uv <- rep(1,(I-1))
        duv<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        for (i in (1:(I-1)))
        {
                eAlpha <- exp(Alpha[I:(I+i-1)])
                unmpstar <- (1-pstar[I:(I-i+1)])
                uv[i] <- prod(eAlpha*unmpstar)*(1-exp(-Alpha[I+i-1]))
                duv[,i] <- prod(eAlpha*unmpstar)*c(rep(0,I+i-1),exp(-Alpha[I+i-1]),rep(0,(length(param)-I-i)))
                for (j in 1:i)
                {
                    duv[,i] <- duv[,i] + (1-exp(-Alpha[I+i-1]))*prod(eAlpha[-j]*unmpstar[-j])*((1-pstar[I-j+1])*c(rep(0,I+j-1),exp(Alpha[I+j-1]),rep(0,(length(param)-I-j)))-exp(Alpha[I+j-1])*dpstar[,I-j+1])
                }
        }
        uv <- c(1,uv)
        duv <- cbind(rep(0,length(param)),duv)
        varcovuv <- t(duv)%*%varcov%*%duv
        uvStderr <- sqrt(diag(varcovuv))    

 #---------------#
# calcul des Vi #
#---------------#

        vv <- rep(1,(I-1))
        dvv<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        for (i in (1:(I-1)))
        {
                eAlpha <- exp(Alpha[1:i])
                unmpstar <- (1-pstar[1:i])
                vv[i] <- prod(eAlpha*unmpstar)*(1-exp(-Alpha[i]))
                dvv[,i] <- prod(eAlpha*unmpstar)*c(rep(0,i),exp(-Alpha[i]),rep(0,(length(param)-i-1)))
                for (j in 1:i)
                {
                    dvv[,i] <- dvv[,i] + (1-exp(-Alpha[i]))*prod(eAlpha[-j]*unmpstar[-j])*((1-pstar[j])*c(rep(0,j),exp(Alpha[j]),rep(0,(length(param)-j-1)))-exp(Alpha[j])*dpstar[,j])
                }
        }
        vv <- c(1,vv)
        dvv <- cbind(rep(0,length(param)),dvv)
        varcovvv <- t(dvv)%*%varcov%*%dvv
        vvStderr <- sqrt(diag(varcovvv))    

#--------------------------------------------------------------#
# calcul des probabilites de survie entre chaque periode (phi) #
#--------------------------------------------------------------#

        phi <- rep(0,(I-1))
        dphi<-matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        phi[(I-1):1] <-1/(1+uv[2:I]/cumsum(uv[1:(I-1)]))   
        for ( i in 1:(I-1))
        {
            if(i==I-1)
            {
            dphi[,i] <- -(phi[i]^2)*duv[,I-i+1]
            } else {
            dphi[,i] <- -(phi[i]^2)*(duv[,I-i+1]*sum(na.rm=TRUE,uv[1:(I-i)])-uv[I-i+1]*rowSums(duv[,1:(I-i)]))/sum(na.rm=TRUE,uv[1:(I-i)])^2
            }
        }
        varcovphi <- t(dphi)%*%varcov%*%dphi
        phistderr <- sqrt(diag(varcovphi))    
  
#---------------------------------------#
# calcul des taille de populations (Ni) #
#---------------------------------------#

        Npop <- rep(0,I)
        dNpop<-matrix(rep(0,I*length(param)),ncol=I)
        
        Npop[1] <- exp(interc)/(prod(1-pstar)*prod(phi))   
        dprodpstar <-rep(0,length(param))
        for (i in 1:I)
        {
            dprodpstar<-dprodpstar-prod(1-pstar[-i])*dpstar[,i]
        }
        dprodphi <-rep(0,length(param))
        for (i in 1:(I-1))
        {
            dprodphi<-dprodphi+prod(phi[-i])*dphi[,i]
        }
        dNpop[,1] <- (prod(1-pstar)*prod(phi)*c(exp(interc),rep(0,(length(param)-1))) - exp(interc)*(prod(phi)*dprodpstar+prod(1-pstar)*dprodphi))/ (prod(1-pstar)*prod(phi))^2
        
        Npop[2:I]<-Npop[1]*cumprod((vv[2:I]/cumsum(vv[1:(I-1)])+1)*phi)
        for (i in 2:I)
        {
            if(i==2) {
                dNpop[,i] <- phi[i-1]*(vv[i]+1)*dNpop[,i-1] + Npop[i-1]*(vv[i]+1)*dphi[,i-1] + Npop[i-1]*phi[i-1]*dvv[,i]          
            } else {
                dNpop[,i] <- phi[i-1]*(vv[i]/sum(na.rm=TRUE,vv[1:(i-1)])+1)*dNpop[,i-1] + Npop[i-1]*(vv[i]/sum(na.rm=TRUE,vv[1:(i-1)])+1)*dphi[,i-1] + Npop[i-1]*phi[i-1]*(sum(na.rm=TRUE,vv[1:(i-1)])*dvv[,i]-vv[i]*rowSums(dvv[,1:(i-1)]))/sum(na.rm=TRUE,vv[1:(i-1)])^2
            }
        }
        varcovtpop <- t(dNpop)%*%varcov%*%dNpop
        NpopStderr <- sqrt(diag(varcovtpop))    
        NpopStderr <- sqrt(pmax(NpopStderr^2-Npop,0))          

#----------------------------#
# calcul des naissances (Bi) #
#----------------------------#

        B<-Npop[2:I]-Npop[1:(I-1)]*phi
        dB <- dNpop[,2:I] - t(phi*t(dNpop[,1:(I-1)])) - t(Npop[1:(I-1)]*t(dphi))
        varcovB <- t(dB)%*%varcov%*%dB
        BStderr <- sqrt(diag(varcovB))
          
#--------------------------------------------------------------------#
# Calcul du nombre total d'individus qui ont passé sur le territoire #
#--------------------------------------------------------------------#

        Ntot <- Npop[1]+sum(na.rm=TRUE,B)
        dNtot <- dNpop[,1] + rowSums(dB)    
        NtotStderr <- sqrt(max(t(dNtot)%*%varcov%*%dNtot-Ntot,0))

#######################################################################################
# Présentation des résultats
#######################################################################################

        modelfit <- matrix(c(anaMrd$deviance,anaMrd$df.residual,anaMrd$aic),nrow=1)
        dimnames(modelfit) <- list("fitted model",c("deviance","    df","      AIC"))
        
        emigfit <- cbind(c(anaMrd2$deviance,anaMrd3$deviance),c(anaMrd2$df.residual,anaMrd3$df.residual),c(anaMrd2$aic,anaMrd3$aic))
        if (Inono==3) { dimnames(emigfit) <- list(c("model with temporary emigration"),c("deviance","    df","      AIC"))
        } else if (Inono>3) { dimnames(emigfit) <- list(c("model with homogeneous temporary emigration","model with temporary emigration"),c("deviance","    df","      AIC"))}    
        
        titre.periode.emig<-rep(0,length(idpemig))
        for (i in 1:length(idpemig)){ titre.periode.emig[i]<-paste("period",idpemig[i])}
        titre.periode<-rep(0,I)
        for (i in 1:I){titre.periode[i]<-paste("period",i)}
        titre.inter.periode<-rep(0,I-1)
        for (i in 1:(I-1)){titre.inter.periode[i]<-paste("period",i,"->",i+1)}
        titre.i<-rep(0,I)
        for (i in 1:I){titre.i[i]<-paste("i =",i-1)}
        
        parap <- rbind(parap3,parap2)     
        pstar <- cbind(pstar,pstarStderr)
        phi <- cbind(phi,phistderr)
        Npop <- cbind(Npop,NpopStderr)
        B <- cbind(round(B,digits=6),round(BStderr,digits=6))
        Ntot <- cbind(Ntot,NtotStderr)
        loglinearpara <- cbind(param,sqrt(diag(varcov)))
        uv <- cbind(uv,uvStderr)
        vv <- cbind(vv,vvStderr)

        if (Inono==3) { dimnames(parap) <- list(titre.periode.emig,c("estimate","stderr")) 
        } else if (Inono>3) { dimnames(parap) <- list(c(titre.periode.emig,"homogenous"),c("estimate","stderr")) }
        dimnames(pstar) <- list(titre.periode,c("estimate","stderr"))
        dimnames(phi) <- list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Npop) <- list(titre.periode,c("estimate","stderr"))
        dimnames(B) <- list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Ntot) <- list("all periods",c("estimate","stderr"))
        dimnames(loglinearpara) <- list(c("intercept",gammanames,Xw$coeffnames),c("estimate","stderr"))
        dimnames(uv) <- list(titre.i,c("estimate","stderr"))
        dimnames(vv) <- list(titre.i,c("estimate","stderr"))

        models<-matrix(Xw$models,nrow=1)
        dimnames(models) <- list("model",titre.periode)

        # Matrice de variances-covariances des paramèters pstar, phi, Npop, B et Ntot
        dP <- cbind(dpstar,dphi,dNpop,dB,dNtot)
        covP <- t(dP)%*%varcov%*%dP
        titre.P <- c(rep(0,4*I-2),"Ntot")
        for (i in 1:I){
            titre.P[i]<-paste("p*",i)
            titre.P[2*I-1+i] <- paste("Npop",i)
        }
        for (i in 1:(I-1)){
            titre.P[I+i]<-paste("phi",i)
            titre.P[3*I-1+i] <- paste("B",i)
        }
        dimnames(covP)<-list(titre.P,titre.P)
      
      
        ans<-list(n=sum(na.rm=TRUE,Y),models=models,model.fit=modelfit,emig.fit=emigfit,emig.param=parap,
                  capture.prob=pstar,survivals=phi,N=Npop,birth=B,Ntot=Ntot,
                  loglin.param=loglinearpara,u.vector=uv,v.vector=vv,cov=covP, neg=positions)
        class(ans) <- "robustd"
        ans                
}


print.robustd <- function(x, ...){
        test <- rep(x$models[1],length(x$models))
        if(all(test==as.character(x$models)))
        {
            cat("\nClosed population model for every period:",x$models[1],"\n")
        } else {
            cat("\nClosed population models for each period:\n")
            print.default(x$models, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        }
        cat("\nModel fit:\n")
        x$model.fit[c(1,3)] <- round(x$model.fit[c(1,3)],3)
        x$model.fit[2] <- round(x$model.fit[2],0)
        print.default(x$model.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        if(!is.null(x$emig.fit)) {
             cat("\nTest for temporary emigration:\n")
             x$emig.fit[,c(1,3)] <- round(x$emig.fit[,c(1,3)],3)
             x$emig.fit[,2] <- round(x$emig.fit[,2],0)
             print.default(x$emig.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        }
        cat("\nCapture probabilities:\n")
        print.default(round(x$capture.prob,4), print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        cat("\nSurvival probabilities:\n")
        print.default(round(x$survivals,4), print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        cat("\nAbundances:\n")
        x$N <- round(x$N,1)
        print.default(x$N, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        cat("\nNumber of new arrivals:\n")
        x$birth <- round(x$birth,1)
        print.default(x$birth, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        cat("\nTotal number of units who ever inhabited the survey area:\n")
        x$Ntot <- round(x$Ntot,1)
        print.default(x$Ntot, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        cat("\nTotal number of captured units:",x$n,"\n")
        ###################################################
        ### 22 mai 2012 : On a décidé de ne plus imprimer ces notes car l'utilisateur ne comprend pas quel
        ### impact des parametres eta fixés à zéro ont sur ses résultats. Ça l'embête plus qu'autre chose.
        #if (length(x$neg[x$neg<2*dim(x$N)[1]])==1) cat("\nNote:",length(x$neg[x$neg<2*dim(x$N)[1]]),"gamma parameter has been set to zero\n")
        #if (length(x$neg[x$neg<2*dim(x$N)[1]])>1) cat("\nNote:",length(x$neg[x$neg<2*dim(x$N)[1]]),"gamma parameters have been set to zero\n")
        #if (length(x$neg[x$neg>=2*dim(x$N)[1]])==1) cat("\nNote:",length(x$neg[x$neg>=2*dim(x$N)[1]]),"eta parameter has been set to zero\n")
        #if (length(x$neg[x$neg>=2*dim(x$N)[1]])>1) cat("\nNote:",length(x$neg[x$neg>=2*dim(x$N)[1]]),"eta parameters have been set to zero\n")
        ###################################################
        cat("\n")
        invisible(x)
}
