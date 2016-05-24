"openp" <- function(X,dfreq=FALSE,m=c("up","ep"),neg=TRUE,keep=rep(TRUE,2^I-1))
{
############################################################################################################################################
# Validation des arguments fournis en entrée
############################################################################################################################################

        valid.one(dfreq,"logical")
        Xvalid <- valid.X(X=X, dfreq=dfreq)
            X <- Xvalid$X
            I <- Xvalid$t
        if (I<=2) stop("loglinear models for open populations require at least 3 capture occasions")
        m <- valid.vm(m,c("up","ep"),I)
        valid.one(neg,"logical")
        if(length(keep)!=2^I-1||!is.logical(keep)) stop("'keep' must be a logical vector of length 2^I-1")


############################################################################################################################################
# AJUSTEMENT DU MODÈLE
############################################################################################################################################

#-------------------------------------#
# Élaboration et ajustement du modèle #
#-------------------------------------#

        Ycomplete <- histfreq.t(X,dfreq=dfreq)   # construction du vecteur nu delta des frequences de captures par periodes
        Y <- Ycomplete
        Y[!keep] <- rep(NA,sum(na.rm=TRUE,!keep))

        gammanames <- paste("gamma",1:(2*I-2),sep="")
        if (m=="up") {
            betanames <- paste("beta",2:(I-1),sep="")
        } else betanames <- "beta"

        histpos <- histpos.t(I)
        Zd <- Zdelta(histpos)      # construction de la matrice Zdelta, premiere composante du modele
        colnames(Zd) <- gammanames
        Xd <- if (m=="up") matrix(histpos[,-c(1,I)],nrow=length(Y)) else matrix(rowSums(histpos),ncol=1)    # on supprime la premiere et dernier colonne de histpos pour des raisons d estimation des parametres
        colnames(Xd) <- betanames
        mX. <- cbind(Zd,Xd)      # on fusionne ces 2 composantes explicatives


        # Ajustement du modèle
        anaMpo <- suppressWarnings(glm(Y~mX.,family=poisson,na.action=na.omit))


        #Vérification du bon ajustement du modèle loglinéaire
        if(!anaMpo$converged) stop("'glm' did not converge")
        if(any(is.na(anaMpo$coef))) warning("the design matrix is not of full rank; some loglinear parameter estimations cannot be evaluated")


#-------------------------------------------------------#
# Réajustement du modèle en enlevant les gamma négatifs #
#-------------------------------------------------------#

# Particularité pour m='up' : 
# Les paramètres gamma1, gammaI-1 et gammaI ne sont pas estimables. On ne teste donc pas leur positivité.
# Ainsi, pour I>3, on vérifie seulement les paramètres gamma 2 à I-2 et I+1 à 2I-2.

        param<-anaMpo$coef
        ppositions <- 0
        test <- if(m=="up") neg && I > 2 else neg  
        if (test)
        {
            # Vecteur d'indicatrices pour les paramètres d'intérêt négatifs
            if (m=="up") {
                indic <- if(I==3) as.vector(c(rep(0,4),ifelse(param[5]<0,1,0),rep(0,I-2))) else as.vector(c(0,0,ifelse(param[3:(I-1)]<0,1,0),0,0,ifelse(param[(I+2):(2*I-1)]<0,1,0),rep(0,I-2)))
            } else indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0),0))
            while(sum(na.rm=TRUE,indic)>0) # Répéter la boucle jusqu'à ce qu'aucun gamma approprié ne soit négatif
            {
                # Détermination de la position du premier gamma approprié négatif
                pos <- 1
                while(indic[pos]==0) pos <- pos + 1
                ppositions <- c(ppositions,pos)
                # Retrait de la bonne colonne de mX. et réajustement du modèle
                mX. <- mX.[,-(pos-sum(na.rm=TRUE,ppositions<pos))]
                anaMpo <- suppressWarnings(glm(Y~mX.,family=poisson,na.action=na.omit))        
                # Ajout de zéros dans le vecteur des paramètres loglinéaires
                positions <- sort(ppositions[-1])                
                param <- c(anaMpo$coef[1:(positions[1]-1)],0)
                if(length(positions)>1)
                {
                    for ( i in 2:length(positions))
                    {
                        if(positions[i]==positions[i-1]+1) {
                            param <- c(param,0)
                        } else {
                            param <- c(param,anaMpo$coef[(positions[i-1]-i+2):(positions[i]-i)],0)
                        }
                    }
                }
                param <- c(param,anaMpo$coef[(positions[length(positions)]-length(positions)+1):length(anaMpo$coef)])
                # Vecteur d'indicatrices pour les paramètres d'intérêt négatifs
                if (m=="up") {
                    indic <- if(I==3) as.vector(c(rep(0,4),ifelse(param[5]<0,1,0),rep(0,I-2))) else as.vector(c(0,0,ifelse(param[3:(I-1)]<0,1,0),0,0,ifelse(param[(I+2):(2*I-1)]<0,1,0),rep(0,I-2)))
                } else indic <- as.vector(c(0,ifelse(param[2:(2*I-1)]<0,1,0),0))
            }
        }
        positions <- sort(ppositions[-1]) 


#------------------------------------------------------#
# Ajustement d'un modèle pour tester l'effet de trappe #
#------------------------------------------------------#

        # Effet de trappe homogène
        if(I==2)              ## inutile stop si I <=2
        {                     ## inutile stop si I <=2
            anaMpo2 <- NULL   ## inutile stop si I <=2
            parap2 <- NULL    ## inutile stop si I <=2
        } else {              ## inutile stop si I <=2
            trap <- rowSums(histpos[,-1]*histpos[,-I])
            mX2. <- cbind(mX.,trap)
            anaMpo2 <- suppressWarnings(glm(Y~mX2.,family=poisson,na.action=na.omit))
            parap2<-summary(anaMpo2)$coef[substr(rownames(summary(anaMpo2)$coef),5,8)=="trap",1:2]
            if (length(parap2)==0)
            {
                anaMpo2 <- NULL
                parap2 <- NULL
            }
        }                     ## inutile stop si I <=2 


        # Effet de trappe hétérogène
        if (m=="up")
        {
            if(I==2||I==3||I==4)
            {   
                anaMpo3 <- NULL
                parap3 <- NULL
            } else {
                trap <- histpos[,-1]*histpos[,-I]
                trapnames <- rep(0,I-1)
                for (i in 1:(I-1)){trapnames[i]<-paste("trap",i,"_",i+1,sep="")}
                colnames(trap) <- trapnames
                mX3.<- cbind(mX.,trap[,-c(1,I-1)])
                anaMpo3 <- suppressWarnings(glm(Y~mX3.,family=poisson,na.action=na.omit))
                parap3<-summary(anaMpo3)$coef[substr(rownames(summary(anaMpo3)$coef),5,8)=="trap",1:2]
                if (length(parap3)==0)
                {
                    anaMpo3 <- NULL
                    parap3 <- NULL
                }
            }
        } else {
            anaMpo3 <- NULL
            parap3 <- NULL
        }
                    


#---------------------------------------#
# Formation des vecteurs des paramètres #
#---------------------------------------#

        # valeur de l intercept
        interc <- param[1]


        # creation des vecteurs de parametres alpha et beta
        Alpha <-rep(0,2*I-2)
        for (i in (1:(2*I-2)))   Alpha[i] <- param[i+1]

        if (m=="up")
        {
            Beta <- rep(0,I-2)
            for (i in 1:(I-2))   Beta[i] <- param[2*I-1+i]
        } else Beta <- param[2*I]


        # Vérification de la présence de paramètres gamma négatifs si l'option "neg"=FALSE
        if(!neg)
        {
            if (m=="ep"&&any(Alpha<0)) 
              warning("one or more gamma parameters are negative,\n",
                      "you can set them to zero with the argument 'neg'.")

            if (m=="up"&&I>3&&any(param[c(3:(I-1),(I+2):(2*I-1))]<0)) 
              warning("one or more relevant gamma parameters are negative,\n",
                      "you can set them to zero with the argument 'neg'.")

            if (m=="up"&&I==3&&param[5]<0) 
              warning("one relevant gamma parameter is negative,\n",
                      "you can set it to zero with the argument 'neg'.")
        }


#--------------------------------------------------------------#
# Matrice de variances-covariances des paramètres loglinéaires #
#--------------------------------------------------------------#

        # Afin de déterminer la position des des paramètres non estimables
        NAindic <- as.vector(is.na(param))
        NApos <- 0
        pos <- 1
        while(pos<=length(param))
        {
            while(!NAindic[pos]&&pos<=length(param)) pos <- pos + 1
            if(pos<=length(param)) NApos <- c(NApos,pos)
            pos <- pos + 1
        }
        NApos <- NApos[-1]       
        
        
        # Pour insérer des lignes et colonnes de zéros pour les paramètres fixés à zéro et pour les paramètres non estimables
        covpos <- sort(c(positions,NApos))
        if(length(covpos)>0)
        {
            # Insertion de colonnes de zéros
            varcovc <-  if(covpos[1]==1) rep(0,dim(summary(anaMpo)$cov.unscaled)[1]) else cbind(summary(anaMpo)$cov.unscaled[,1:(covpos[1]-1)],rep(0,dim(summary(anaMpo)$cov.unscaled)[1]))
            if(length(covpos)>1)
            {
                for ( i in 2:length(covpos))
                {
                    if(covpos[i]==covpos[i-1]+1) {
                        varcovc <- cbind(varcovc,rep(0,dim(summary(anaMpo)$cov.unscaled)[1]))
                    } else {
                        varcovc <- cbind(varcovc,summary(anaMpo)$cov.unscaled[,(covpos[i-1]-i+2):(covpos[i]-i)],rep(0,dim(summary(anaMpo)$cov.unscaled)[1]))
                    }
                }
            }
            if(covpos[length(covpos)]<length(param)) varcovc <- cbind(varcovc,summary(anaMpo)$cov.unscaled[,(covpos[length(covpos)]-length(covpos)+1):dim(summary(anaMpo)$cov.unscaled)[2]])
            # Insertion de lignes de zéros
            varcov <- if(covpos[1]==1) rep(0,dim(varcovc)[2]) else rbind(varcovc[1:(covpos[1]-1),],rep(0,dim(varcovc)[2]))
            if(length(covpos)>1)
            {
                for ( i in 2:length(covpos))
                {
                    if(covpos[i]==covpos[i-1]+1) {
                        varcov <- rbind(varcov,rep(0,dim(varcovc)[2]))
                    } else {
                        varcov <- rbind(varcov,varcovc[(covpos[i-1]-i+2):(covpos[i]-i),],rep(0,dim(varcovc)[2]))
                    }
                }
            }
            if(covpos[length(covpos)]<length(param)) varcov <- rbind(varcov,varcovc[(covpos[length(covpos)]-length(covpos)+1):dim(varcovc)[1],]) 
        } else { varcov <- summary(anaMpo)$cov.unscaled }
        


############################################################################################################################################
# Estimation des paramètres démographiques
############################################################################################################################################

#--------------------------------------------#
# calcul des probabilites de capture (pstar) #
#--------------------------------------------#

        if (m=="up")
        {
            pstar <- rep(0,I)
            dpstar<-matrix(rep(0,I*length(param)),ncol=I)
            pstar[1] <- 0.5
            pstar[I] <- 0.5  
            if(I>2)
            {
                for (i in (2:(I-1)))
                {
                        pstar[i]<- exp(Beta[i-1])/(1+exp(Beta[i-1]))
                        dpstar[2*I+i-2,i] <- exp(Beta[i-1])/(1+exp(Beta[i-1]))^2
                }
            }
            varcovpstar <- t(dpstar)%*%varcov%*%dpstar
            pstarStderr <- sqrt(diag(varcovpstar))    
        } else {
            pstar<- rep(exp(Beta)/(1+exp(Beta)),I)
            dpstar <- matrix(rep(c(rep(0,2*I-1),exp(Beta)/(1+exp(Beta))^2),I),ncol=I)       
            varcovpstar <- t(dpstar)%*%varcov%*%dpstar
            pstarStderr <- sqrt(diag(varcovpstar))
        }


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
        dphi <- matrix(rep(0,(I-1)*length(param)),ncol=(I-1))
        
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
        phiStderr <- sqrt(diag(varcovphi))    


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

        # Correction si certains historiques ont été enlevés (avec l'option keep)
        corrkeep <- suppressWarnings(sum(na.rm=TRUE,Ycomplete[!keep]-predict(anaMpo,newdata=data.frame(mX.),type="response")[!keep]))
        Npop <- Npop + corrkeep


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

        if (m=="up")
        {
            #Programme pour former la matrice X avec des colonnes de zéros pour les paramètres fixés à zéro
            Xtemp<-cbind(rep(1,dim(mX.)[1]),mX.)        
            if(length(positions)>0)
            {
                # Insertion de colonnes de zéros
                X <- cbind(Xtemp[,1:(positions[1]-1)],rep(0,dim(Xtemp)[1]))
                if(length(positions)>1)
                {
                    for ( i in 2:length(positions))
                    {
                        if(positions[i]==positions[i-1]+1)
                        {
                            X <- cbind(X,rep(0,dim(Xtemp)[1]))
                        } else {
                            X <- cbind(X,Xtemp[,(positions[i-1]-i+2):(positions[i]-i)],rep(0,dim(Xtemp)[1]))
                        }
                    }
                }
                X <- cbind(X,Xtemp[,(positions[length(positions)]-length(positions)+1):dim(Xtemp)[2]])
            } else { X <- Xtemp}
    
    
            # Calcul du nombre total d'individus qui ont passé sur le territoire
            if (I==2)
            {
                Ntot <- sum(na.rm=TRUE,Ycomplete)
                dNtot <- t(X)%*%exp(X%*%param)
            } else if (I==3)
                {
                    Ntot <- sum(na.rm=TRUE,Ycomplete) + exp(interc+Alpha[1]+Alpha[3])-exp(interc+Alpha[1])
                    dNtot <- t(X)%*%exp(X%*%param) + exp(interc+Alpha[1]+Alpha[3])*c(1,1,0,1,rep(0,2)) - exp(interc+Alpha[1])*c(1,1,rep(0,4))
                } else if (I==4)
                    {
                        Ntot <- sum(na.rm=TRUE,Ycomplete) + exp(interc+Alpha[1]+Alpha[4]+Alpha[5]) + exp(interc+Alpha[1]+Alpha[4])*(exp(Alpha[2])-1)
                        dNtot <- ( t(X)%*%exp(X%*%param) + exp(interc+Alpha[1]+Alpha[4]+Alpha[5])*c(1,1,0,0,1,1,rep(0,3)) 
                        + exp(interc+Alpha[1]+Alpha[2]+Alpha[4])*c(1,1,1,0,1,rep(0,4)) - exp(interc+Alpha[1]+Alpha[4])*c(1,1,0,0,1,rep(0,4)) )
                    } else {
                        Ntot <- sum(na.rm=TRUE,Ycomplete) + exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-3)]))*(1+exp(sum(na.rm=TRUE,Alpha[2:(I-2)])-sum(na.rm=TRUE,Alpha[(I+1):(2*I-3)]))) - exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-4)]))
                        dNtot <- ( t(X)%*%exp(X%*%param) 
                        + exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-3)]))*exp(sum(na.rm=TRUE,Alpha[2:(I-2)])-sum(na.rm=TRUE,Alpha[(I+1):(2*I-3)]))*c(0,0,rep(1,I-3),0,0,rep(-1,I-3),rep(0,length(param)-2*I+2)) 
                        + exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-3)]))*(1+exp(sum(na.rm=TRUE,Alpha[2:(I-2)])-sum(na.rm=TRUE,Alpha[(I+1):(2*I-3)])))*c(1,1,rep(0,I-2),rep(1,I-2),rep(0,length(param)-2*I+2)) 
                        - exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-4)]))*c(1,1,rep(0,I-2),rep(1,I-3),rep(0,length(param)-2*I+3)) )
                        for (i in 2:(I-3))
                        {
                            Ntot <- Ntot + exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-3)]))*exp(sum(na.rm=TRUE,Alpha[2:i])-sum(na.rm=TRUE,Alpha[(2*I-1-i):(2*I-3)])) - exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-4)]))*exp(sum(na.rm=TRUE,Alpha[2:i])-sum(na.rm=TRUE,Alpha[(2*I-2-i):(2*I-4)]))
                            dNtot <- ( dNtot + exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-3)]))*exp(sum(na.rm=TRUE,Alpha[2:i])-sum(na.rm=TRUE,Alpha[(2*I-1-i):(2*I-3)]))*
                            (c(1,1,rep(0,I-2),rep(1,I-2),rep(0,length(param)-2*I+2)) + c(0,0,rep(1,i-1),rep(0,2*I-2*i-2),rep(-1,i-1),rep(0,length(param)-2*I+2)))
                            - exp(interc+Alpha[1]+sum(na.rm=TRUE,Alpha[I:(2*I-4)]))*exp(sum(na.rm=TRUE,Alpha[2:i])-sum(na.rm=TRUE,Alpha[(2*I-2-i):(2*I-4)]))*
                            (c(1,1,rep(0,I-2),rep(1,I-3),rep(0,length(param)-2*I+3)) + c(0,0,rep(1,i-1),rep(0,2*I-2*i-3),rep(-1,i-1),rep(0,length(param)-2*I+3))) )
                        }
                    }
            NtotStderr <- sqrt(max(t(dNtot)%*%varcov%*%dNtot-Ntot,0))
        } else {
            Ntot <- Npop[1]+sum(na.rm=TRUE,B)
            dNtot <- dNpop[,1] + rowSums(dB)    
            NtotStderr <- sqrt(max(t(dNtot)%*%varcov%*%dNtot-Ntot,0))
        }    


############################################################################################################################################
# Présentation des résultats
############################################################################################################################################
        
        modelfit <- matrix(c(anaMpo$deviance,anaMpo$df.residual,anaMpo$aic),nrow=1)
        dimnames(modelfit) <- list("fitted model",c("deviance","    df","      AIC"))
        
        trapfit <- cbind(c(anaMpo2$deviance,anaMpo3$deviance),c(anaMpo2$df.residual,anaMpo3$df.residual),c(anaMpo2$aic,anaMpo3$aic))
        if (!is.null(trapfit)) {
          if (dim(trapfit)[1]==1) { dimnames(trapfit) <- list(c("model with homogenous trap effect"),c("deviance","    df","      AIC"))
          } else if (dim(trapfit)[1]==2) { dimnames(trapfit) <- list(c("model with homogenous trap effect","model with trap effect"),c("deviance","    df","      AIC"))}
        }
                 
        titre.periode<-paste("period",1:I)
        titre.inter.periode<-paste("period",1:(I-1),"->",1:(I-1)+1)
        titre.i<-paste("i =",1:I-1)
            
        parap <- rbind(parap3,parap2)     
        pstar <- cbind(pstar,pstarStderr)
        phi <- cbind(phi,phiStderr)
        Npop <- cbind(Npop,NpopStderr)
        B <- cbind(B,BStderr)
        Ntot <- cbind(Ntot,NtotStderr)
        loglinearpara <- cbind(param,diag(varcov))
        uv <- cbind(uv,uvStderr)
        vv <- cbind(vv,vvStderr)
        
        if (length(parap)==2) { dimnames(parap) <- list("homogenous trap effect",c("estimate","stderr")) 
        } else if (length(parap)>2) { dimnames(parap) <- list(c(paste("period",substr(rownames(parap3),9,11)),"homogenous trap effect"),c("estimate","stderr")) }
        dimnames(pstar)<-list(titre.periode,c("estimate","stderr"))
        dimnames(phi)<-list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Npop)<-list(titre.periode,c("estimate","stderr"))
        dimnames(B)<-list(titre.inter.periode,c("estimate","stderr"))
        dimnames(Ntot)<-list("all periods",c("estimate","stderr"))
        dimnames(loglinearpara)<-list(c("intercept",gammanames,betanames),c("estimate","stderr"))
        dimnames(uv) <- list(titre.i,c("estimate","stderr"))
        dimnames(vv) <- list(titre.i,c("estimate","stderr"))

        #Paramètres non estimables
        if (m=="up")
        {
            pstar[1,] <- rep(NA,2)
            pstar[I,] <- rep(NA,2)
            phi[I-1,] <- rep(NA,2)
            Npop[1,] <- rep(NA,2)
            Npop[I,] <- rep(NA,2)
            B[1,] <- rep(NA,2)
            B[I-1,] <- rep(NA,2)
        }

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
        if (m=="up") covP <- covP[-c(1,I,2*I-1,2*I,3*I-1,3*I,4*I-2),-c(1,I,2*I-1,2*I,3*I-1,3*I,4*I-2)]
      
        
        ans<-list(n=sum(na.rm=TRUE,Ycomplete),model.fit=modelfit,trap.fit=trapfit,trap.param=parap,capture.prob=pstar,survivals=phi,N=Npop,birth=B,Ntot=Ntot,
                  glm=anaMpo,loglin.param=loglinearpara,u.vector=uv,v.vector=vv,cov=covP, neg=positions)
        class(ans) <- "openp"
        ans

}


print.openp <- function(x, ...){

        cat("\nModel fit:\n")
        x$model.fit[,c(1,3)] <- round(x$model.fit[,c(1,3)],3)
        x$model.fit[,2] <- as.integer(x$model.fit[,2])
        print.default(x$model.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
        if (!is.null(x$trap.fit)) {
          cat("\nTest for trap effect:\n")
          x$trap.fit[,c(1,3)] <- round(x$trap.fit[,c(1,3)],3)
          x$trap.fit[,2] <- round(x$trap.fit[,2],0)
          print.default(x$trap.fit, print.gap = 2, quote = FALSE, right=TRUE, na.print="--", ...)
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
        #if (length(x$neg[x$neg<2*dim(x$N)[1]])>1) cat("\nNote:",length(x$neg[x$neg<2*dim(x$N)[1]]),"gamma parameters has been set to zero\n")
        ###################################################
        cat("\n")
        invisible(x)
}

plot.openp <- function(x,main="Scatterplot of Pearson Residuals", ...){
    res <- (x$glm$model[,1]- fitted.values(x$glm))/sqrt(fitted.values(x$glm))
    if (length(x$glm$na.action)==0) 
    {
        plot(apply(histpos.t(dim(x$N)[1]),1,sum),res,xlab="Frequency of capture",ylab="Pearson residuals",main=main, ...)
    } else {
        plot(apply(histpos.t(dim(x$N)[1]),1,sum)[-x$glm$na.action],res,xlab="Frequency of capture",ylab="Pearson residuals",main=main, ...)
    }
}
