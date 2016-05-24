
cnca <- function(fespecies, fvambientales)
{
        if(missing(fvambientales))
        {
                msg<-("ERROR: this function requires two matrices as arguments")
                tkmessageBox(message=msg)
                stop(" this function requires two matrices as arguments")
        }#end if(missing(fvambientales))
        
        if(dim(fespecies)[1]!=dim(fvambientales)[1])
        {
                msg<-("ERROR: the same number of rows in both matrices is required")
                tkmessageBox(message=msg)
                stop("the same number of rows in both matrices is required")
        }#end if(dim(fespecies)[1]!=dim(fvambientales)[1]) 
        
        transforma <- function(tipo, matriz)
        {    
                if (tipo=="Subtract the global mean"){                        
                        Xstd <- as.matrix(matriz)
                        media <- mean(Xstd)
                        Xstd <- Xstd-media		
                }#end if (tipo=="Subtract the global mean")
                
                if (tipo=="Column centering"){
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 2, function(x){x-mean(as.matrix(x))})
                }#end if (tipo=="Column centering")
                
                if (tipo=="Standardize columns"){	
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 2, function(x){(x-mean(as.matrix(x)))/sqrt(var(as.matrix(x)))})
                }#end if (tipo=="Standardize columns")
                
                if (tipo=="Row centering"){		
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 1, function(x){x-mean(as.matrix(x))})
                        Xstd <- t(Xstd)
                }#end if (tipo=="Row centering")
                
                if (tipo=="Standardize rows"){
                        Xstd <- as.matrix(matriz)
                        Xstd <-apply(matriz, 1, function(x){(x-mean(as.matrix(x)))/sqrt(var(as.matrix(x)))})
                        Xstd <- t(Xstd)
                }#end if (tipo=="Standardize rows")
                
                if (tipo=="Double centering"){
                        Xstd <- as.matrix(matriz)
                        mediac <- colMeans(Xstd)
                        mediaf <- rowMeans(Xstd)
                        globalm <- mean(Xstd)
                        
                        mediafm <- array(unlist(rep(rowMeans(Xstd), dim(Xstd)[2])), dim=dim(Xstd))
                        mediacm <- t(array(unlist(rep(colMeans(Xstd), dim(Xstd)[1])), dim=c(dim(Xstd)[2],dim(Xstd)[1])))
                        Xstd <- Xstd - mediafm - mediacm +globalm
                        
                }#end if (tipo=="Double centering")
                
                if (tipo=="Raw data"){
                        Xstd <- as.matrix(matriz)
                }#end if (tipo=="Raw data")
                
                rownames(Xstd) <- rownames(matriz)
                return(Xstd)        
                
        }# end transforma 
        
        ex.cnca<-function(fespecies, fvambientales, nejes, Numbercuant, Numbermixto, tChoice)
        {
                ##We save the dimensions
                dimespec <- dim(fespecies)
                dimvambien <- dim(fvambientales)
                
                ##We create the data matrices
                especies <- array(data=unlist(fespecies), dim=c(dimespec[1],(dimespec[2])))
                vambientales <- array(data=unlist(fvambientales),dim=c(dimvambien[1],(dimvambien[2])))
                vambientalesst<-vambientales
                vambientalesst[,1:Numbercuant] <- transforma(tChoice, vambientales[,1:Numbercuant])
                
                
                ##labels
                textlugares <- rownames(fespecies)
                textespecies <- colnames(fespecies)
                textvariables <- colnames(fvambientales)
                vambientalesst <- array(data=unlist(vambientalesst),dim=c(dim(vambientalesst)))
                mixtas <- vambientalesst[,-c(1:Numbercuant)]
                
                
                ##We sum all elements of species
                sumatotal <- sum(especies)
                Fe <- especies/sumatotal
                
                fn <- array(dim=c(dimespec[1],1))
                
                ##We calculate the marginal totals and we save it in fn and fq
                Idq <- as.matrix(rep(1,times=dimespec[2]))
                fn <- Fe%*%Idq
                
                Idn <- as.matrix(rep(1,times=dimespec[1]))
                fq <- t(Fe)%*%Idn
                
                #Matrix whose diagonal is the marginals fn 
                Dn <- diag(as.vector(t(fn)))
                #Dn <- diag(rep(1/dimespec[1], times=dimespec[1]))
                
                invD <- ginv(Dn)
                
                Dq <- diag(as.vector(t(fq)))
                invDq <- ginv(Dq)
                
                P <- invD%*%Fe
                Pb <- P-(Idn%*%t(fq))
                
                covar <- t(vambientalesst)%*%Dn%*%vambientalesst
                
                descovar <- La.svd(covar)
                Pi <- vambientalesst%*%ginv(covar)%*%t(vambientalesst)%*%Dn
                Pbest <- Pi%*%Pb
                
                ##########################################################################################
                ###Singular Value Decomposition of Pbest
                ##########################################################################################
                
                descpbest <- La.svd(Pbest,nu=nejes,nv=nejes) 
                
                ##########################################################################################
                ###Singular Value Decomposition of Lp
                ##########################################################################################
                
                Fb <- Fe-(fn%*%t(fq))
                L <- t(Fb)%*%vambientalesst
                Lp <- (descovar$u%*%diag((descovar$d)^(-1/2))%*%descovar$v)%*%t(L)
                descom <- svd(Lp)
                suma2valprop <- sum((descom$d[1:length(descom$d)])^2)
                sumaRvalprop <- sum((descom$d)^2)
                inerciatot <- (descom$d[1:length(descom$d)])^2/sumaRvalprop
                
                desclp <- svd(Lp)
                V <- diag(desclp$d[1:nejes])
                Tn <- desclp$v[,1:nejes]
                A <- desclp$u[,1:nejes]
                #A <- Lp %*% Tn %*% solve(V)
                #RV<<-Pbest %*% Tn
                R<- Pbest %*% Tn %*% solve(V)
                #R <- descpbest$u
                #V<-diag(descpbest$d[1:nejes]) 
                #Tn<-t(descpbest$v)
                ###########################################################################################
                ##Principal coordinates of sites
                ###########################################################################################
                
                colugares<- Pbest %*% Tn
                
                ###########################################################################################
                ##Standard coordinates of species
                ###########################################################################################
                
                coespecies<- Tn
                
                ###########################################################################################
                ##Mixed data
                ###########################################################################################        			
                Idcuant <- diag(rep(1,times=Numbercuant))
                
                
                if(Numbermixto>0)
                {
                        NAinv <- ginv(diag(as.vector(t(mixtas)%*%Dn%*%t(t(rep(1,times=dimespec[1]))))))
                        N <- rbind(
                                cbind(Idcuant,array(rep(0,times=dim(Idcuant)[1]*dim(NAinv)[2]),dim=c(dim(Idcuant)[1],dim(NAinv)[2]))),
                                cbind(array(rep(0,times=dim(NAinv)[1]*dim(Idcuant)[2]),dim=c(dim(NAinv)[1],dim(Idcuant)[2])),NAinv)
                        )
                }else{
                        N <- Idcuant
                }
                
                
                ###########################################################################################
                ##Principal coordinates of environmental variables
                ###########################################################################################				
                
                covambien <- N%*%((descovar$u%*%(diag(descovar$d)^(1/2))%*%descovar$v))%*%A%*%V
                
                ##########################################################################
                #### Contributions and qualities of representation
                ##########################################################################
                ejes <- sapply(1:nejes,f<-function(x){paste("Axis", x)})
                
                ####Contributions of the sites
                colugarescuad <- (Pbest %*% desclp$v)^2
                fixi <- Dn %*% colugarescuad
                Calphai <- fixi[,1:nejes]
                coespeciescuad <- (desclp$v%*%diag(desclp$d))^2
                Calphak <- coespeciescuad[,1:nejes]
                
                for (i in 1:nejes)
                {
                        # fixi[,i] <- fn[i]*colugarescuad[,i]
                        Calphai[,i] <- (fixi[,i]*1000)/(diag(V)[i])^2
                        
                        ####Contributions of the species
                        Calphak[,i] <- (coespeciescuad[,i]*1000)/(diag(V)[i])^2
                }#end for (i in 1:nejes)
                
                rownames(Calphai) <- textlugares
                colnames(Calphai) <- ejes
                rownames(Calphak) <- textespecies
                colnames(Calphak) <- ejes
                
                ####qualities of representation respect to the projected space
                variabilidadps <- sum(desclp$d^2)
                calidadpst <- diag(V)^2/variabilidadps
                calidadps <- t(calidadpst)*1000
                colnames(calidadps) <- ejes
                
                ####qualities of representation respect to the original space
                descomPb<-svd(Dn^(1/2)%*%Pb)
                variabilidados <- sum((descomPb$d)^2)
                calidadost <- diag(V)[1:nejes]^2/variabilidados
                calidados <- t(calidadost)*1000
                colnames(calidados) <- ejes
                
                
                colugarescuados<-(Pb%*%descomPb$v)^2
                coespeciescuados<-(t(Pb)%*%Dn^(1/2)%*%descomPb$u)^2
                
                qalphaips <- (t(apply(colugarescuad,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                dalphaios <- (t(apply(colugarescuados,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                qalphakps <- (t(apply(coespeciescuad,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                qalphakos <- (t(apply(coespeciescuados,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                
                
                rownames(qalphaips) <- textlugares
                colnames(qalphaips) <- ejes
                
                rownames(dalphaios) <- textlugares
                colnames(dalphaios) <- ejes
                
                rownames(qalphakps) <- textespecies
                colnames(qalphakps) <- ejes
                
                rownames(qalphakos) <- textespecies
                colnames(qalphakos) <- ejes
                
                ##### retained inertia
                eigencuad <- (diag(V))^2
                inertia <- cbind(eigencuad,eigencuad)
                
                inertia[,1] <- (eigencuad/sum(eigencuad))*100
                inertia[,2] <- cumsum(inertia[,1])
                colnames(inertia) <- c("Retained inertia (%)","Cumulative retained inertia(%)")
                
                coindividuosnam <- as.data.frame(coespecies)
                rownames(coindividuosnam) <- textespecies
                colnames(coindividuosnam) <- ejes
                
                covariablesnam <- as.data.frame(covambien)
                rownames(covariablesnam) <- textvariables
                colnames(covariablesnam) <- ejes
                
                colugaresnam <- as.data.frame(colugares)
                rownames(colugaresnam) <- textlugares
                colnames(colugaresnam) <- ejes
                
                results<-list(colugares=colugares, coespecies=coespecies, covambien=covambien, colugaresnam=colugaresnam, covariablesnam=covariablesnam, coindividuosnam=coindividuosnam, ejes=ejes, Calphai=Calphai, 
                              Calphak=Calphak, calidadps=calidadps, calidados=calidados, qalphaips=qalphaips, qalphakos=qalphakos, qalphakps=qalphakps, 
                              dalphaios=dalphaios, inertia=inertia, descom=descom, inerciatot=inerciatot, valpro=descom$d)
                
                return(results)
        }
        
        ex.cca<-function(fespecies, fvambientales, nejes, Numbercuant, Numbermixto, tChoice)
        {
                dimespec <- dim(fespecies)
                dimvambien <- dim(fvambientales)
                
                ##We create the data matrices
                especies <- array(data=unlist(fespecies), dim=c(dimespec[1],(dimespec[2])))
                vambientales <- array(data=unlist(fvambientales),dim=c(dimvambien[1],(dimvambien[2])))
                vambientalesst<-vambientales
                vambientalesst[,1:Numbercuant] <- transforma(tChoice, vambientales[,1:Numbercuant])
                
                
                
                ##labels
                textlugares <- rownames(fespecies)
                textespecies <- colnames(fespecies)
                textvariables <- colnames(fvambientales)
                vambientalesst <- array(data=unlist(vambientalesst),dim=c(dim(vambientalesst)))
                mixtas <- vambientalesst[,-c(1:Numbercuant)]
                
                
                ##We sum all elements of species
                sumatotal <- sum(especies)
                Fe <- especies/sumatotal
                
                fn <- array(dim=c(dimespec[1],1))
                
                ##We calculate the marginal totals and we save it in fn and fq
                Idq <- as.matrix(rep(1,times=dimespec[2]))
                fn <- Fe%*%Idq
                
                Idn <- as.matrix(rep(1,times=dimespec[1]))
                fq <- t(Fe)%*%Idn
                
                #Matrix whose diagonal is the marginals fn 
                Dn <- diag(as.vector(t(fn)))
                # Dn <- diag(rep(1/dimespec[1], times=dimespec[1]))
                
                invD <- ginv(Dn)
                
                Dq <- diag(as.vector(t(fq)))
                invDq <- ginv(Dq)
                
                covar <- t(vambientalesst)%*%Dn%*%vambientalesst
                
                descovar <- La.svd(covar)
                
                ##center species matrix
                rsum<-rowSums(Fe)
                csum <- colSums(Fe)
                # Fe<- Fe-(rsum%*%t(csum))
                
                W <- (invDq %*% t(Fe)-t(rsum%*%t(csum))) %*% vambientalesst
                Wdesc <- Dq^(1/2) %*% W %*% (descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v)
                ##########################################################################################
                ###Singular Value Decomposition of Wdesc
                ##########################################################################################
                
                descom <- La.svd(Wdesc)
                
                
                R <- descom$u[,1:nejes]
                V <- diag(descom$d[1:nejes]) 
                Tn <- t(descom$v)[,1:nejes]
                
                
                suma2valprop <- sum((descom$d[1:length(descom$d)])^2)
                sumaRvalprop <- sum((descom$d)^2)
                inerciatot <- (descom$d[1:length(descom$d)])^2/sumaRvalprop
                
                
                ###########################################################################################
                ##Principal coordinates of sites
                ###########################################################################################
                
                colugares<- vambientalesst %*% (descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v)  %*% Tn 
                
                ###########################################################################################
                ##Standard coordinates of species
                ###########################################################################################
                
                coespecies<- invDq^(1/2) %*% R %*% V
                
                ###########################################################################################
                ##Mixed data
                ###########################################################################################                		
                Idcuant <- diag(rep(1,times=Numbercuant))
                
                
                if(Numbermixto>0)
                {
                        NAinv <- ginv(diag(as.vector(t(mixtas)%*%Dn%*%t(t(rep(1,times=dimespec[1]))))))
                        N <- rbind(
                                cbind(Idcuant,array(rep(0,times=dim(Idcuant)[1]*dim(NAinv)[2]),dim=c(dim(Idcuant)[1],dim(NAinv)[2]))),
                                cbind(array(rep(0,times=dim(NAinv)[1]*dim(Idcuant)[2]),dim=c(dim(NAinv)[1],dim(Idcuant)[2])),NAinv)
                        )
                }else{
                        N <- Idcuant
                }
                
                ###########################################################################################
                ##Standar coordinates of environmental variables
                ###########################################################################################				
                
                covambien <- N %*% ((descovar$u%*%(diag(descovar$d)^(1/2))%*%descovar$v)) %*% Tn %*% V
                
                ##########################################################################
                #### Contributions and qualities of representation
                ##########################################################################
                ejes <- sapply(1:nejes,f<-function(x){paste("Axis", x)})
                #         Pbest <- (invDq %*% t(Fe)-t(rsum%*%t(csum))) - R %*% V %*% t((descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v) %*% Tn) %*% t(vambientalesst) %*% Dn
                #         Pbest <- t(Pbest)
                
                Pbest<-invD^(1/2) %*% (Fe-(rsum%*%t(csum))) %*% invDq^(1/2)
                Rp <- svd(Pbest)#,nu=nejes,nv=nejes)
                
                
                
                colugaresos <- Rp$u%*%diag(Rp$d)
                PHIu <-Rp$u
                coespeciesos <- invDq^(1/2)%*%Rp$v%*%diag(Rp$d)
                GAHu<- invDq^(1/2)%*%Rp$v
                
                ####Contributions of the sites
                colugarescuad <-Dn %*% (vambientalesst %*% (descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v)  %*% t(descom$v)%*%diag(descom$d))^2
                Calphai <- colugarescuad[,1:nejes]
                coespeciescuad <- Dq %*% (invDq^(1/2) %*% descom$u %*% diag(descom$d))^2
                Calphak <- coespeciescuad[,1:nejes]
                
                colugarescuados <-(Dn %*% colugaresos)^2
                Calphaios <- colugarescuados[,1:nejes]
                coespeciescuados <- (Dq %*% coespeciesos)^2
                Calphakos <- coespeciescuados[,1:nejes]
                
                for (i in 1:nejes)
                {
                        Calphai[,i] <- (colugarescuad[,i]*1000)/(diag(V)[i])^2
                        ####Contributions of the species
                        Calphak[,i] <- (coespeciescuad[,i]*1000)/(diag(V)[i])^2
                        
                        Calphaios[,i] <- (colugarescuados[,i]*1000)/(Rp$d[i])^2
                        
                        ####Contributions of the species
                        Calphakos[,i] <- (coespeciescuados[,i]*1000)/(Rp$d[i])^2
                }#end for (i in 1:nejes)
                
                rownames(Calphai) <- textlugares
                colnames(Calphai) <- ejes
                rownames(Calphak) <- textespecies
                colnames(Calphak) <- ejes
                
                ####qualities of representation respect to the projected space
                variabilidadps <- sum(descom$d^2)
                calidadpst <- diag(V)^2/variabilidadps
                calidadps <- t(calidadpst)*1000
                colnames(calidadps) <- ejes
                
                ####qualities of representation respect to the original space
                variabilidados <- sum(Rp$d)
                calidadost <- Rp$d[1:nejes]/variabilidados
                calidados <- t(calidadost)*1000
                colnames(calidados) <- ejes
                
                qalphaips <- (t(apply(colugarescuad,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                dalphaios <- (t(apply(colugarescuados,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                qalphakps <- (t(apply(coespeciescuad,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                qalphakos <- (t(apply(coespeciescuados,1,function(x){x/sum(x)}))*1000)[,1:nejes]
                
                
                rownames(qalphaips) <- textlugares
                colnames(qalphaips) <- ejes
                
                rownames(dalphaios) <- textlugares
                colnames(dalphaios) <- ejes
                
                rownames(qalphakps) <- textespecies
                colnames(qalphakps) <- ejes
                
                rownames(qalphakos) <- textespecies
                colnames(qalphakos) <- ejes
                
                ##### retained inertia
                eigencuad <- (diag(V))^2
                inertia <- cbind(eigencuad,eigencuad)
                
                inertia[,1] <- (eigencuad/sum(eigencuad))*100
                inertia[,2] <- cumsum(inertia[,1])
                colnames(inertia) <- c("Retained inertia (%)","Cumulative retained inertia(%)")
                
                coindividuosnam <- as.data.frame(coespecies)
                rownames(coindividuosnam) <- textespecies
                colnames(coindividuosnam) <- ejes
                
                covariablesnam <- as.data.frame(covambien)
                rownames(covariablesnam) <- textvariables
                colnames(covariablesnam) <- ejes
                
                colugaresnam <- as.data.frame(colugares)
                rownames(colugaresnam) <- textlugares
                colnames(colugaresnam) <- ejes
                
                results<-list(colugares=colugares, coespecies=coespecies, covambien=covambien, colugaresnam=colugaresnam, covariablesnam=covariablesnam, coindividuosnam=coindividuosnam, ejes=ejes, 
                              Calphai=Calphai, Calphak=Calphak, calidadps=calidadps, calidados=calidados, qalphaips=qalphaips, qalphakos=qalphakos, qalphakps=qalphakps, 
                              dalphaios=dalphaios, inertia=inertia, 
                              descom=descom, inerciatot=inerciatot, valpro=descom$d)
                
                return(results)
        }
        
        
        ex.coia<-function(fespecies, fvambientales, nejes, tChoice, Numbercuant)
        {
                ##We save the dimensions
                dimespec <- dim(fespecies)
                dimvambien <- dim(fvambientales)
                
                ##We create the data matrices
                especies <- array(data=unlist(fespecies), dim=c(dimespec[1],(dimespec[2])))
                especies <- transforma("Column centering", especies)
                vambientales <- array(data=unlist(fvambientales),dim=c(dimvambien[1],(dimvambien[2])))
                vambientalesst <- transforma (tChoice, vambientales[,1:Numbercuant])
                
                ##labels
                textlugares <- rownames(fespecies)
                textespecies <- colnames(fespecies)
                textvariables <- colnames(fvambientales)
                
                ##We sum all elements of species
                sumatotal <- sum(especies)
                Fe <- especies/sumatotal
                fvam <-fvambientales/sum(fvambientales)
                fn <- array(dim=c(dimespec[1],1))
                
                ##We calculate the marginal totals and we save it in fn and fq
                Idq <- as.matrix(rep(1,times=dimespec[2]))
                fn <- Fe%*%Idq
                
                Idn <- as.matrix(rep(1,times=dimespec[1]))
                fq <- t(Fe)%*%Idn
                
                fp <- t(fvam)%*%Idn
                
                #Matrix whose diagonal is the marginals fn 
                #Dn <- diag(as.vector(t(fn)))
                Dn <- diag(rep(1/dimespec[1], times=dimespec[1]))
                Dn2<-diag(sqrt(diag(Dn)))
                invD <- ginv(Dn)
                
                Dq <- diag(as.vector(t(fq)))
                invDq <- ginv(Dq)
                
                Dp <- diag(as.vector(t(fp)))
                invDp <- ginv(Dp)
                
                #         ##center species matrix
                #         rsum<-rowSums(Fe)
                #         csum <- colSums(Fe)
                #         Fe<- Fe-(rsum%*%t(csum))
                #         
                intertabla <- t(fespecies) %*% Dn %*% vambientalesst
                colnames(intertabla)<- colnames(fvambientales)
                
                
                # extrae los vectores singulares por la izquierda y por la derecha para las dos 
                #matrices y para la coinercia
                
                descomx <- svd(Dn2 %*% vambientalesst, nu=nejes, nv=nejes)
                UX <- solve(Dn2) %*% (descomx$u) %*% diag(1/sqrt(diag(t(descomx$u)%*%(descomx$u))))
                VX <- (descomx$v) %*% diag(1/sqrt(diag(t(descomx$v)%*%(descomx$v))))
                descomy <- svd(Dn2 %*% especies, nu=nejes, nv=nejes)
                UY <- solve(Dn2) %*% (descomy$u) %*% diag(1/sqrt(diag(t(descomy$u)%*%(descomy$u))))
                VY <- (descomy$v) %*% diag(1/sqrt(diag(t(descomy$v)%*%(descomy$v))))
                descom <- svd(intertabla, nu=nejes, nv=nejes)
                Ucoia <- (descom$u) %*% diag(1/sqrt(diag(t(descom$u)%*%(descom$u))))
                Vcoia <- (descom$v) %*% diag(1/sqrt(diag(t(descom$v)%*%(descom$v))))
                Xaxes <- t(VX) %*% Vcoia
                Yaxes <- t(VY)%*%Ucoia
                
                suma2valprop <- sum((descom$d[1:length(descom$d)])^2)
                sumaRvalprop <- sum((descom$d)^2)
                inerciatot <- (descom$d[1:length(descom$d)])^2/sumaRvalprop
                
                # calcula las coordenadas 
                rowx <- vambientalesst %*% VX
                colx <- t(vambientalesst) %*% Dn %*% UX
                rowy <- especies %*% VY
                coly <- t(especies) %*% Dn %*% UY
                colugaresz <- vambientalesst %*% Vcoia
                covambien <- t(intertabla) %*% Ucoia
                colugaresy <- especies %*% Ucoia
                coespecies <- intertabla %*% Vcoia
                
                
                
                ##########################################################################
                #### Contributions and qualities of representation
                ##########################################################################
                ejes <- sapply(1:nejes,f<-function(x){paste("Axis", x)})
                
                coindividuosnam <- as.data.frame(coespecies)
                rownames(coindividuosnam) <- textespecies
                colnames(coindividuosnam) <- ejes
                
                covariablesnam <- as.data.frame(covambien)
                rownames(covariablesnam) <- textvariables
                colnames(covariablesnam) <- ejes
                
                colugaresnamz <- as.data.frame(colugaresz)
                rownames(colugaresnamz) <- textlugares
                colnames(colugaresnamz) <- ejes
                
                colugaresnamy <- as.data.frame(colugaresy)
                rownames(colugaresnamy) <- textlugares
                colnames(colugaresnamy) <- ejes
                
                results<-list(colugaresz=colugaresz, colugaresy=colugaresy, coespecies=coespecies, 
                              covambien=covambien, colugaresnamz=colugaresnamz, colugaresnamy=colugaresnamy,
                              covariablesnam=covariablesnam, coindividuosnam=coindividuosnam, ejes=nejes, 
                              descom=descom, inerciatot=inerciatot, valpro=descom$d, Xaxes=Xaxes, Yaxes=Yaxes)
                
                return(results)
        }
        
        calc_inerciatot<-function(fespecies, fvambientales, tChoice, tipocca, Numbercuant)
        {
                ##We save the dimensions
                dimespec <- dim(fespecies)
                dimvambien <- dim(fvambientales)
                
                ##We create the data matrices
                especies <- array(data=unlist(fespecies), dim=c(dimespec[1],(dimespec[2])))
                vambientales <- array(data=unlist(fvambientales),dim=c(dimvambien[1],(dimvambien[2])))
                vambientalesst<-vambientales
                vambientalesst[,1:Numbercuant] <- transforma(tChoice, vambientales[,1:Numbercuant])
                
                ##We sum all elements of species
                sumatotal <- sum(especies)
                Fe <- especies/sumatotal
                fn <- array(dim=c(dimespec[1],1))
                
                ##We calculate the marginal totals and we save it in fn and fq
                Idq <- as.matrix(rep(1,times=dimespec[2]))
                fn <- Fe%*%Idq
                
                Idn <- as.matrix(rep(1,times=dimespec[1]))
                fq <- t(Fe)%*%Idn
                
                #Matrix whose diagonal is the marginals fn 
                Dn <- diag(as.vector(t(fn)))
                #Dn <- diag(rep(1/dimespec[1], times=dimespec[1]))
                
                invD <- ginv(Dn)
                
                Dq <- diag(as.vector(t(fq)))
                invDq <- ginv(Dq)
                covar <- t(vambientalesst)%*%Dn%*%vambientalesst
                descovar <- La.svd(covar)
                
                ##########################################################################################
                ###Singular Value Decomposition of Lp
                ##########################################################################################
                if(tipocca=="CNCA")
                {
                        Fb <- Fe-(fn%*%t(fq))
                        L <- t(Fb)%*%vambientalesst
                        Lp <- ((descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v))%*%t(L)
                        descom <- svd(Lp)#,nu=nejes,nv=nejes) 
                }else{
                        if(tipocca=="CCA")
                        {
                                rsum<-rowSums(Fe)
                                csum <- colSums(Fe)
                                # Fe<- Fe-rsum%*%t(csum)
                                
                                W <- (invDq %*% t(Fe)-t(rsum%*%t(csum))) %*% vambientalesst
                                Wdesc <- Dq^(1/2) %*% W %*% (descovar$u%*%(solve(diag(descovar$d)^(1/2)))%*%descovar$v)
                                
                                ##########################################################################################
                                ###Singular Value Decomposition of Wdesc
                                ##########################################################################################
                                
                                descom <- svd(Wdesc)
                                
                        }else{
                                intertabla <- t(especies) %*% Dn %*% vambientalesst
                                descom <-svd(intertabla)
                        }                        
                }
                suma2valprop <- sum((descom$d[1:length(descom$d)])^2)
                sumaRvalprop <- sum((descom$d)^2)
                inerciatot <- (descom$d[1:length(descom$d)])^2/sumaRvalprop
                return(list(descom=descom, inerciatot=inerciatot))
        }
        
        
        resample_bootes<-function(fespecies)
        {
                mientorno <- new.env()
                dimespec<-dim(fespecies)
                #deshacer tabla de contingencia
                especvector<-paste("e", 1:dimespec[2], sep="")
                
                for(z in 1:dimespec[2])
                {
                        esvec<-c()
                        for(j in 1: dimespec[1])
                        {
                                esvec<-c(esvec, rep(j,fespecies[j,z])) 
                        }
                        
                        assign(especvector[z], esvec, envir=mientorno)
                }
                
                ##### resample
                especresample<-array(rep(0,dimespec[1]*dimespec[2]), dim=dimespec)
                
                for(z in 1:dimespec[2])
                {
                        columresam<-get(especvector[z], envir=mientorno)
                        indices <- sample(1:length(columresam), replace = T)
                        columresam<-columresam[indices]
                        columna<-table(columresam)
                        for (j in 1:length(columna))
                        {
                                especresample[as.numeric(rownames(columna)[j]),z]<-columna[j]
                        }
                        
                }
                colnames(especresample)<-colnames(fespecies)
                rownames(especresample)<-rownames(fespecies) 
                return(especresample)
        }# end resample_bootes<-function(fespecies)
        
        cal.ic <- function (muestra, liminf, limsup, valorobs, muestrajack, niter)
        {
                c.mean <- mean (muestra)
                se <- sd(muestra)
                sesgo <- c.mean - valorobs
                t.ic <- se * (-qt(liminf,(length(muestra)-1)))
                ic.t <- c(c.mean - t.ic, c.mean + t.ic)
                ic.p <- quantile (muestra,c(liminf, limsup), na.rm=TRUE)
                z0 <- qnorm(length(muestra[which(muestra<valorobs)])/as.numeric(niter))
                dent <- mean(muestrajack)- muestrajack
                acc <- sum(dent * dent * dent)/(6 * (sum(dent * dent))^1.5)
                alpha1 <- qnorm(liminf)
                alpha2 <- qnorm(limsup)
                zalpha1 <- pnorm(z0 + (z0 + alpha1)/(1 - acc * (z0 + alpha1)))
                zalpha2 <- pnorm(z0 + (z0 + alpha2)/(1 - acc * (z0 + alpha2)))
                ic.bca <- quantile (muestra,c(zalpha1, zalpha2), na.rm=TRUE)
                return(c(c.mean, se,sesgo,ic.t, ic.p, ic.bca))
        }#end cal.ic <- function (muestra, liminf, limsup, valorobs)
        
        
        #############################################################################
        ### Informative window
        #############################################################################
        symbolos <- c("*",".", "o","O","0","+","-","|","%","#")
        nejes<-3  
        dim1<-1
        dim2<-2
        dim3<-3 
        dim1ant<-1
        dim2ant<-2
        simChoicel <- NULL
        Namesit <- NULL
        Namespe <- NULL
        Namevar <- NULL
        Namemixto <- NULL
        Numbercuant <- NULL
        Numbermixto <- NULL
        dimespec <-NULL
        dimvambien <-NULL
        especies <-NULL
        vambientales <-NULL
        vambientalesst <-NULL
        simlugares <-NULL
        simespecies <-NULL
        simvariables <-NULL
        anteriorx <- NULL
        anteriory <- NULL
        xCoords <- NULL
        yCoords <- NULL
        xCoordsp <- NULL
        yCoordsp <- NULL
        zCoords <- NULL
        datos <- NULL
        textos <- NULL
        datosr <- NULL
        textosr <- NULL
        simbolos <- NULL
        colores <-NULL
        indexClosest <- NULL
        indexLabeled <- NULL                           
        indexLabeledaux <- NULL
        parPlotSize <- NULL
        usrCoords <- NULL
        tChoice <- "Raw data"
        img <- NULL
        imgbar <- NULL
        descom <- NULL
        inerciatot <- NULL
        Dn <- NULL
        fq <- NULL
        fn <- NULL
        Fe <- NULL
        Idn <- NULL
        Lp <- NULL
        colugares <- NULL
        coespecies <- NULL
        covambien <- NULL
        ejes <- NULL
        indicee <- NULL
        Namee <- NULL
        Cexe <- 1
        NameVale <- NULL
        NameCexe <- NULL
        Choicee <- NULL
        colore <- NULL
        simChoicee <- NULL
        indicev <- NULL
        Namev <- NULL
        NameValv <- NULL
        Cexv <- 1
        NameCexv <- NULL
        Choicev <- NULL
        colorv <- NULL
        simChoicev <- NULL
        indicei <- NULL
        Namei <- NULL
        NameVali <- NULL
        Cexi <- 1
        NameCexi <- NULL
        Choicei <- NULL
        colori <- NULL
        simChoicei <- NULL
        entry.Cexi <- NULL	
        entry.Namei <- NULL
        indicel <- NULL
        Namel <- NULL
        NameVall <- NULL
        Cexl <- 1
        NameCexl <- NULL
        Choicel <- NULL
        colorl <- NULL
        simChoicel <- NULL
        entry.Cexl <- NULL
        cblugares <- "0"
        cbVal <- NULL
        cblablugares <- "0"
        hescale <- "1.5"
        vescale <- "1.5"
        labelsVec <- NULL
        sizesVec <- NULL
        simChoice <- NULL
        textvariables <- c()
        textespecies <- c()
        textlugares <- c()
        colvariablesp <- c()
        colespeciesp <- c()
        collugaresp <- c()
        proj <- "normal"
        Choiceproj<- 0
        niter <- 1000
        alphaic <- 95
        cCalphaiVal <- NULL
        cCalphakVal <- NULL
        ccalidadpsVal <- NULL
        ccalidadosVal <- NULL
        cqalphaipsVal <- NULL
        cdalphaiosVal <- NULL
        cqalphakpsVal <- NULL
        cqalphakosVal <- NULL
        cccaVal <- NULL
        iner <- NULL
        colugaresnam <- NULL
        covariablesnam <- NULL
        coindividuosnam <- NULL
        Calphai <- NULL
        Calphak <- NULL
        calidadps <- NULL
        calidados <- NULL
        qalphaips <- NULL
        qalphakps <- NULL
        qalphakos <- NULL
        dalphaios <- NULL
        inertia <- NULL
        colugaresz <- NULL
        colugaresy <- NULL
        colugaresnamz <- NULL
        colugaresnamy <- NULL
        Xaxes <- NULL
        Yaxes <- NULL
        colugzaux <- NULL
        colugyaux <- NULL
        cexejes <- NULL
        entry.limix1 <- NULL
        entry.limix2 <- NULL
        entry.limiy1 <- NULL
        entry.limiy2 <- NULL
        datost <- NULL
        #ccoindividuosnamVal <- NULL
        #ccovariablesnamVal <- NULL
        #ccolugaresnamVal <- NULL
        ceigenVal <- NULL
        cpdfVal <- NULL
        cepsVal <- NULL
        ccaVal <- NULL
        typecoia <- "sg"
        tit_graph <- "Graph"
        Limix1 <-tclVar("0")
        Limix2 <-tclVar("0")
        Limiy1 <-tclVar("0")
        Limiy2 <-tclVar("0")
        
        
        colorescoor <- c("skyblue","red","green","blue","yellow","pink","orange", "navyblue",
                         "violet", "brown", "grey", "navyblue", "darkgreen", "papayawhip", "paleturquoise", "purple",
                         "seagreen", "azure", "coral", "springgreen", "steelblue", "plum", "orchid", 
                         "lemonchiffon", "lavender", "honeydew", "gold", "deeppink", "darksalmon", "darkmagenta")
        
        
        winfor <- tktoplevel()
        tkwm.title(winfor,"Two Tables Analyses")
        
        
        fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
        fontTextLabel <- tkfont.create(family="times",size=12)
        fontFixedWidth <- tkfont.create(family="courier",size=12)
        tkpack(tklabel(winfor,text="Two Tables Analyses",font=fontHeading, foreground = "blue"))
        
        frameprin<-tkframe(winfor, relief = "ridge", borderwidth = 2, background = "white")
        framecnca<-tkframe(frameprin, relief = "ridge", borderwidth = 2, background = "white")
        frameot<-tkframe(frameprin, relief = "ridge", borderwidth = 2, background = "white")
        
        
        ccca <- tkradiobutton(frameot)
        ccnca <- tkradiobutton(framecnca)
        ccoia <- tkradiobutton(frameot)
        rbccaValue <- tclVar("CNCA")
        tkconfigure(ccca,variable=rbccaValue,value="CCA")
        tkconfigure(ccnca,variable=rbccaValue,value="CNCA")
        tkconfigure(ccoia,variable=rbccaValue,value="COIA")
        tkpack(tklabel(framecnca, text="Canonical Non-Symmetrical Correspondence Analysis (CNCA)",
                       font=fontTextLabel, foreground = "blue"), ccnca,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(frameot, text="Canonical Correspondence Analysis (CCA)"), ccca,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(frameot, text="Coinertia Analysis (COIA)"), ccoia,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        
        tkpack(framecnca, frameot, expand = "TRUE", side="left", fill="both")
        tkpack(frameprin, frameot, expand = "TRUE", side="top", fill="both")
        
        ##We save the dimensions
        dimespec <- dim(fespecies)
        dimvambien <- dim(fvambientales)
        
        ###########################################################################################
        ##We create the vectors of the labels
        ###########################################################################################
        
        textlugares <- rownames(fespecies)
        textespecies <- colnames(fespecies)
        textvariables <- colnames(fvambientales)
        
        ###########################################################################################
        ##We create the vectors of the colors
        ###########################################################################################
        
        collugares <- rep("green",dimespec[1])
        colespecies <- rep("blue",dimespec[2])
        colvariables <- rep("red",dimvambien[2])
        
        collugaresp <- rep("grey",dimespec[1])
        colespeciesp <- rep("grey",dimespec[2])
        colvariablesp <- rep("grey",dimvambien[2])
        
        
        ###########################################################################################
        ##We create the vectors of the symbols
        ###########################################################################################
        
        simlugares <- rep("*",dimespec[1])
        simespecies <- rep("+",dimespec[2])
        simvariables <- rep(" ",dimvambien[2])
        
        ###########################################################################################
        ##We create the vectors of the character size
        ###########################################################################################
        
        cexlugares <- rep(1,dimespec[1])
        cexespecies <- rep(1,dimespec[2])        
        cexvariables <- rep(1,dimvambien[2])
        
        OnOKinf <- function()
        {
                #############################################################################
                ### Names of label window
                #############################################################################
                cccaVal <<- as.character(tclvalue(rbccaValue))
                if((cccaVal=="CCA" | cccaVal=="CNCA") & (dim(fespecies)[1]<dim(fespecies)[2]) & (dim(fvambientales)[1]<dim(fvambientales)[2]))
                {
                        msg<-("ERROR: the number of rows is less than the number of columns. The best method for this analysis is the Coinertia Analysis.")
                        tkmessageBox(message=msg)
                        tkdestroy(winfor)
                        stop("the number of rows is less than the number of columns. The best method for this analysis is the Coinertia Analysis")
                }#end  if((cccaVal=="CCA" | cccaVal=="CNCA") & (dim(fespecies)[1]<dim(fespecies)[2]) & (dim(fvambientales)[1]<dim(fvambientales)[2]))
                
                
                tkdestroy(winfor)
                wnames <- tktoplevel()
                tkwm.title(wnames,"Enter names")      
                
                OnOKnames <- function()
                {
                        Namesit <<- tclvalue(sitesname)
                        Namespe <<- tclvalue(speciesname)
                        Namevar <<- tclvalue(variablesname)
                        if(cccaVal!="COIA")
                        {
                                #Namemixto <<- tclvalue(mixto)
                                Numbercuant <<- tclvalue(ncuant)
                                Numbermixto <<- tclvalue(nmixto)
                                Numbermixto <<- as.numeric(Numbermixto)
                                Numbercuant <<- as.numeric(Numbercuant)
                                
                                if((Numbercuant+Numbermixto)!=dim(fvambientales)[2])
                                {
                                        msg<-("ERROR: Total number of variables not equal to sum of cuantitatives and categorical variables")
                                        tkmessageBox(message=msg)
                                        stop("Total number of variables not equal to sum of cuantitatives and categorical variables")
                                }#end if(dim(fespecies)[1]!=dim(fvambientales)[1]) 
                                
                                
                        }else{
                                Numbercuant<<-dimvambien[2]
                                Numbermixto<<-0
                        }
                        
                        tkdestroy(wnames)
                        
                        
                        
                        #####Window to change labels and colors#######
                        
                        tt <- tktoplevel()
                        tkwm.title(tt,"Options")
                        
                        #####Dropdown menu #############################
                        topMenuopt <- tkmenu(tt)
                        tkconfigure(tt, menu = topMenuopt)
                        fileMenutrans <- tkmenu(topMenuopt, tearoff = FALSE)
                        
                        tkadd(fileMenutrans,"command",label="Subtract the global mean",command=function() tChoice<<-"Subtract the global mean")
                        tkadd(fileMenutrans,"command",label="Column centering",command=function() tChoice<<-"Column centering")
                        tkadd(fileMenutrans,"command",label="Standardize columns",command=function() tChoice<<-"Standardize columns")
                        tkadd(fileMenutrans,"command",label="Row centering",command=function() tChoice<<-"Row centering")
                        tkadd(fileMenutrans,"command",label="Standardize rows",command=function() tChoice<<-"Standardize rows")
                        tkadd(fileMenutrans,"command",label="Double Centering",command=function() tChoice<<-"Double centering")
                        tkadd(fileMenutrans,"command",label="Raw data",command=function() tChoice<<-"Raw data")
                        tkadd(topMenuopt,"cascade",label="Transformations",menu=fileMenutrans)
                        
                        
                        #### List of transformations
                        framett<-tkframe(tt, relief = "flat", borderwidth = 2, background = "white")
                        framett1<-tkframe(framett, relief = "ridge", borderwidth = 2, background = "white")
                        framett2<-tkframe(framett, relief = "ridge", borderwidth = 2, background = "white")
                        framett3<-tkframe(framett, relief = "ridge", borderwidth = 2, background = "white")
                        framett4<-tkframe(framett, relief = "ridge", borderwidth = 2)#, background = "grey3")
                        
                        framet1<-tkframe(framett1, relief = "ridge", borderwidth = 2, background = "white")
                        frametext1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                        frameok1<-tkframe(framett1, relief = "ridge", borderwidth = 2, background = "white")
                        
                        framecol1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                        framecol11<-tkframe(framecol1, relief = "flat", borderwidth = 2, background = "white")
                        framecol12<-tkframe(framecol1, relief = "flat", borderwidth = 2, background = "white")
                        
                        framename1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                        framename11<-tkframe(framename1, relief = "flat", borderwidth = 2, background = "white")
                        framename12<-tkframe(framename1, relief = "flat", borderwidth = 2, background = "white")
                        
                        framecex1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")
                        framecex11<-tkframe(framecex1, relief = "flat", borderwidth = 2, background = "white")
                        framecex12<-tkframe(framecex1, relief = "flat", borderwidth = 2, background = "white")
                        
                        frames1<-tkframe(framett1, relief = "flat", borderwidth = 2, background = "white")      
                        frames11<-tkframe(frames1, relief = "flat", borderwidth = 2, background = "white")
                        frames12<-tkframe(frames1, relief = "flat", borderwidth = 2, background = "white")
                        
                        framet2<-tkframe(framett2, relief = "ridge", borderwidth = 2, background = "white")
                        frametext2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                        frameok2<-tkframe(framett2, relief = "ridge", borderwidth = 2, background = "white")
                        
                        framecol2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                        framecol21<-tkframe(framecol2, relief = "flat", borderwidth = 2, background = "white")
                        framecol22<-tkframe(framecol2, relief = "flat", borderwidth = 2, background = "white")
                        
                        framename2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                        framename21<-tkframe(framename2, relief = "flat", borderwidth = 2, background = "white")
                        framename22<-tkframe(framename2, relief = "flat", borderwidth = 2, background = "white")
                        
                        framecex2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                        framecex21<-tkframe(framecex2, relief = "flat", borderwidth = 2, background = "white")
                        framecex22<-tkframe(framecex2, relief = "flat", borderwidth = 2, background = "white")
                        
                        frames2<-tkframe(framett2, relief = "flat", borderwidth = 2, background = "white")
                        frames21<-tkframe(frames2, relief = "flat", borderwidth = 2, background = "white")
                        
                        framet3<-tkframe(framett3, relief = "ridge", borderwidth = 2, background = "white")
                        frametext3<-tkframe(framett3, relief = "flat", borderwidth = 2, background = "white")
                        frameok3<-tkframe(framett3, relief = "ridge", borderwidth = 2, background = "white")
                        
                        framecol3<-tkframe(framett3, relief = "flat", borderwidth = 2, background = "white")
                        framecol31<-tkframe(framecol3, relief = "flat", borderwidth = 2, background = "white")
                        framecol32<-tkframe(framecol3, relief = "flat", borderwidth = 2, background = "white")
                        
                        framename3<-tkframe(framett3, relief = "flat", borderwidth = 2, background = "white")
                        framename31<-tkframe(framename3, relief = "flat", borderwidth = 2, background = "white")
                        framename32<-tkframe(framename3, relief = "flat", borderwidth = 2, background = "white")
                        
                        framecex3<-tkframe(framett3, relief = "flat", borderwidth = 2, background = "white")
                        framecex31<-tkframe(framecex3, relief = "flat", borderwidth = 2, background = "white")
                        framecex32<-tkframe(framecex3, relief = "flat", borderwidth = 2, background = "white")
                        
                        frames3<-tkframe(framett3, relief = "flat", borderwidth = 2, background = "white")
                        frames31<-tkframe(frames3, relief = "flat", borderwidth = 2, background = "white")
                        frames32<-tkframe(frames3, relief = "flat", borderwidth = 2, background = "white")  
                        
                        framet4<-tkframe(framett4, relief = "flat", borderwidth = 2)#, background = "white")
                        frametext4<-tkframe(framett4, relief = "flat", borderwidth = 2)#, background = "white")
                        frameok4<-tkframe(framett4, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        framett4aux<-tkframe(framett4, relief = "ridge", borderwidth = 2)#, background = "white")
                        framett4auxgs<-tkframe(framett4, relief = "ridge", borderwidth = 2)#, background = "white")
                        
                        framecol4<-tkframe(framett4aux, relief = "flat", borderwidth = 2)#, background = "white")
                        framecol41<-tkframe(framecol4, relief = "flat", borderwidth = 2)#, background = "white")
                        framecol42<-tkframe(framecol4, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        framename4<-tkframe(framett4aux, relief = "flat", borderwidth = 2)#, background = "white")
                        framename41<-tkframe(framename4, relief = "flat", borderwidth = 2)#, background = "white")
                        framename42<-tkframe(framename4, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        framehvtitle<-tkframe(framett4auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                        framehv<-tkframe(framett4auxgs, relief = "flat", borderwidth = 2)#, background = "white")
                        framehvnames<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                        framehnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                        framevnames<-tkframe(framehvnames, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        framehvtext<-tkframe(framehv, relief = "flat", borderwidth = 2)#, background = "white")
                        framehtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                        framevtext<-tkframe(framehvtext, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        framecex4<-tkframe(framett4aux, relief = "flat", borderwidth = 2)#, background = "white")
                        framecex41<-tkframe(framecex4, relief = "flat", borderwidth = 2)#, background = "white")
                        framecex42<-tkframe(framecex4, relief = "flat", borderwidth = 2)#, background = "white")
                        
                        frames4<-tkframe(framett4, relief = "flat", borderwidth = 2)#, background = "white")
                        frames41<-tkframe(frames4, relief = "flat", borderwidth = 2)#, background = "white")  
                        framegraphic<-tkframe(tt, relief = "flat", borderwidth = 2, background = "white")
                        
                        
                        ##### Checkbox to show the axes or not  #######
                        
                        cb <- tkcheckbutton(framecol42)
                        cbValue <- tclVar("0")
                        tkconfigure(cb,variable=cbValue)
                        
                        ##### Checkbox to show the sites or not  #######
                        
                        cbl <- tkcheckbutton(framename42)
                        cblug <- tclVar("0")
                        tkconfigure(cbl,variable=cblug)
                        
                        ##### Checkbox to show the labels of the sites or not  #######
                        
                        cbll <- tkcheckbutton(framecex42)
                        cblablug <- tclVar("0")
                        tkconfigure(cbll,variable=cblablug)
                        
                        ##We save the dimensions
                        dimespec <- dim(fespecies)
                        
                        if(Numbermixto>0)
                        {
                                mixtas<-array(dim=c(dimespec[1],1))
                                for (i in 1:Numbermixto)
                                {
                                        factoresmix<-unique(fvambientales[,as.numeric(Numbercuant)+i])
                                        for(j in 1:length(factoresmix))
                                        {
                                                linea<-ifelse(fvambientales[,Numbercuant+i]==factoresmix[j], 1,0)
                                                linea<-as.data.frame(linea)
                                                colnames(linea)<-paste(colnames(fvambientales)[as.numeric(Numbercuant)+i],factoresmix[j],sep="-")
                                                mixtas<-cbind(mixtas,linea)
                                        }
                                }
                                
                                mixtas<-mixtas[,-1]
                                fvambientales <- cbind(fvambientales[,1:Numbercuant],mixtas)
                                
                        } 
                        
                        dimvambien <- dim(fvambientales)
                        
                        ###########################################################################################
                        ##We create the vectors of the labels
                        ###########################################################################################
                        
                        textlugares <- rownames(fespecies)
                        textespecies <- colnames(fespecies)
                        textvariables <- colnames(fvambientales)
                        
                        ###########################################################################################
                        ##We create the vectors of the colors
                        ###########################################################################################
                        
                        collugares <- rep("green",dimespec[1])
                        colespecies <- rep("blue",dimespec[2])
                        colvariables <- rep("red",dimvambien[2])
                        
                        collugaresp <- rep("grey",dimespec[1])
                        colespeciesp <- rep("grey",dimespec[2])
                        colvariablesp <- rep("grey",dimvambien[2])
                        
                        ###########################################################################################
                        ##We create the vectors of the symbols
                        ###########################################################################################
                        
                        simlugares <- rep("*",dimespec[1])
                        simespecies <- rep("+",dimespec[2])
                        simvariables <- rep(" ",dimvambien[2])
                        
                        ###########################################################################################
                        ##We create the vectors of the character size
                        ###########################################################################################
                        
                        cexlugares <- rep(1,dimespec[1])
                        cexespecies <- rep(1,dimespec[2])        
                        cexvariables <- rep(1,dimvambien[2])
                        
                        
                        ##### List of species###########################
                        scre <- tkscrollbar(framet1, repeatinterval=5, command=function(...)tkyview(tle,...))
                        tle<-tklistbox(framet1,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scre,...),background="white")
                        tkpack(tklabel(frametext1,text=Namespe),side="left",expand = "TRUE",fill="both")
                        
                        for (i in 1:(dimespec[2]))
                        {
                                tkinsert(tle,"end",textespecies[i])
                        }#end for (i in 1:(dimespec[2]))
                        
                        tkselection.set(tle,0) #  Indexing starts at zero.
                        
                        OnOKe <- function()
                        {
                                Choicee <<- textespecies[as.numeric(tkcurselection(tle))+1]
                                
                                ##### Color of the selected variable #############
                                indicee<<-as.numeric(tkcurselection(tle))+1
                                colore <<- colespecies[indicee[1]]
                                tkconfigure(canvase,bg=colore)
                                
                                ##### Text of the selected variable###############
                                
                                Namee <<- tclVar(textespecies[indicee[1]])
                                tkconfigure(entry.Namee,textvariable=Namee)
                                
                                ##### Size of the selected variable###############
                                
                                Cexe <<- tclVar(cexespecies[indicee[1]])
                                tkconfigure(entry.Cexe,textvariable=Cexe) 
                        }#end OnOKe <- function()
                        
                        OK.bute <-tkbutton(frameok1,text="    OK    ",command=OnOKe)
                        tkpack(tle,scre,expand = "TRUE", side="left", fill = "both")
                        tkpack.configure(scre,side="left")
                        tkpack(OK.bute,expand = "TRUE", side="left", fill = "both")
                        tkfocus(tt)
                        
                        #######Color#######################################
                        indicee <- as.numeric(tkcurselection(tle))+1
                        colore <<- colespecies[indicee[1]]
                        canvase <- tkcanvas(framecol11,width="57",height="20",bg=colore)
                        ChangeColore <- function()
                        {
                                
                                colore <<- tclvalue(tcl("tk_chooseColor",initialcolor=colespecies[indicee[1]],title="Choose a color"))
                                
                                if (nchar(colore)>0)
                                {
                                        tkconfigure(canvase,bg=colore)
                                        colespecies[indicee]<<-colore
                                }#end if (nchar(colore)>0)
                        }#end ChangeColore <- function()
                        
                        ChangeColor.buttone<- tkbutton(framecol12,text="Change Color",command=ChangeColore,width=4)
                        tkpack(canvase,expand = "TRUE", side="left", fill = "both")
                        tkpack(ChangeColor.buttone,expand = "TRUE", side="left", fill = "both")
                        
                        ######Labels  ###################################
                        Namee <- textespecies[indicee[1]]
                        entry.Namee <-tkentry(framename11,width="10",textvariable=Namee, bg="white")
                        
                        OnOKle <- function()
                        {
                                NameVale <<- tclvalue(Namee)
                                textespecies[indicee[1]] <-NameVale
                                
                                #####Values of listbox###############################
                                for (i in 1:(dimespec[2]))
                                {
                                        tkdelete(tle,0)
                                }#end for (i in 1:(dimespec[2]))
                                
                                for (i in (1:(dimespec[2])))
                                {
                                        tkinsert(tle,"end",textespecies[i])
                                }#end for (i in (1:(dimespec[2])))
                        }#end OnOKle <- function()
                        
                        OK.butle <-tkbutton(framename12,text="Change label",command=OnOKle,width=4)
                        tkbind(entry.Namee, "<Return>",OnOKle)
                        tkpack(entry.Namee,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butle,expand = "TRUE", side="left", fill = "both")
                        
                        ###### Sizes  ###################################
                        Cexe <- cexespecies[indicee[1]]
                        entry.Cexe <-tkentry(framecex11,width="10",textvariable=Cexe, bg="white")
                        
                        OnOKce <- function()
                        {
                                NameCexe <<- tclvalue(Cexe)
                                cexespecies[indicee] <<-NameCexe
                        }#end OnOKce <- function()
                        
                        OK.butce <-tkbutton(framecex12,text="Change size",command=OnOKce,width=4)
                        tkbind(entry.Cexe, "<Return>",OnOKce)
                        tkpack(entry.Cexe,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butce,expand = "TRUE", side="left", fill = "both")
                        
                        ######Symbols###################################
                        
                        comboBoxe <- tkwidget(frames11,"ComboBox",editable=FALSE,values=symbolos,width=7, background="white")
                        
                        chang.syme <- function()
                        {
                                simChoicee <<- symbolos[as.numeric(tclvalue(tcl(comboBoxe,"getvalue")))+1]
                                simespecies[indicee] <<-simChoicee
                        }#end chang.syme <- function()
                        
                        Change.symbole <-tkbutton(frames12,text="Change symbol",command=chang.syme,width=4)
                        tkpack(comboBoxe,side="left",expand="TRUE", fill="both")
                        tkpack(Change.symbole,side="left",expand="TRUE", fill="both")
                        
                        ##### List of environmental v. ###########################
                        
                        scrv <- tkscrollbar(framet2, repeatinterval=5, command=function(...)tkyview(tlv,...))
                        tlv<-tklistbox(framet2,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scrv,...),background="white")
                        tkpack(tklabel(frametext2,text=Namevar),side="left",expand = "TRUE",fill="both")
                        
                        for (i in (1:(dimvambien[2])))
                        {
                                tkinsert(tlv,"end",textvariables[i])
                        }#end for (i in (1:(dimvambien[2])))
                        
                        tkselection.set(tlv,0) #  Indexing starts at zero.
                        
                        OnOKv <- function()
                        {
                                Choicev <<- textvariables[as.numeric(tkcurselection(tlv))+1]
                                
                                ##### Color of the selected variable  #############
                                indicev<<-as.numeric(tkcurselection(tlv))+1
                                colorv <<- colvariables[indicev[1]]
                                tkconfigure(canvasv,bg=colorv)
                                
                                ##### Text of the selected variable  ###############
                                
                                Namev <<- tclVar(textvariables[indicev[1]])
                                tkconfigure(entry.Namev,textvariable=Namev)
                                
                                ##### Size of the selected variable  ###############
                                
                                Cexv <<- tclVar(cexvariables[indicev[1]])
                                tkconfigure(entry.Cexv,textvariable=Cexv)
                        }#end OnOKv <- function()
                        
                        OK.butv <-tkbutton(frameok2,text="    OK    ",command=OnOKv)
                        
                        tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")
                        
                        Graphics <- function()
                        {
                                
                                tkdestroy(tt)
                                hescale <<- tclvalue(entryvalueh)
                                vescale <<- tclvalue(entryvaluev)
                                
                                barvp<-tktoplevel()
                                tkwm.title(barvp,"Eigenvalues")
                                
                                
                                plotbar<-function()
                                {
                                        iner <<- calc_inerciatot(fespecies, fvambientales, tChoice, cccaVal,Numbercuant)
                                        inerciatot <<- iner[[2]]
                                        barplot(iner[[1]]$d, col="blue", xlab="", ylab="", names.arg=round(inerciatot, digits=2))
                                }#end plotbar<-function()
                                
                                imgbar <- tkrplot(barvp,fun=plotbar,hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                                msginertia <- "Proportion of inertia explained by each axis:"
                                for (i in 1:length(iner[[1]]$d))
                                {
                                        msginertia <- paste(msginertia, "\n",i, "\t", round(inerciatot[i]*100, digits=2), "%")
                                }#end for (i in 1:length(iner[[1]]$d))
                                tk2tip(imgbar, msginertia)
                                Onaxis <- function()
                                {
                                        nejes <- tclvalue(numaxis)
                                        nejes <- as.numeric(nejes)
                                        if (nejes > length(iner[[1]]$d))
                                        {
                                                msg <- paste("The maximum number of dimensions is ",length(iner[[1]]$d))
                                                tkmessageBox(message=msg)
                                                
                                        }else{
                                                tkdestroy(barvp)
                                                nejes <- as.integer(nejes)
                                                
                                                if (cccaVal=="CNCA")
                                                {
                                                        cnca_res <- ex.cnca(fespecies, fvambientales, nejes, Numbercuant, Numbermixto, tChoice)
                                                        colugares <<- cnca_res$colugares
                                                        coespecies <<- cnca_res$coespecies
                                                        covambien <<- cnca_res$covambien
                                                        colugaresnam <<-cnca_res$colugaresnam
                                                        covariablesnam <<- cnca_res$covariablesnam
                                                        coindividuosnam <<-cnca_res$coindividuosnam
                                                        ejes <<- cnca_res$ejes
                                                        Calphai <<- cnca_res$Calphai 
                                                        Calphak <<- cnca_res$Calphak
                                                        calidadps <<- cnca_res$calidadps
                                                        calidados <<- cnca_res$calidados
                                                        qalphaips <<- cnca_res$qalphaips
                                                        qalphakos <<- cnca_res$qalphakos
                                                        qalphakps <<- cnca_res$qalphakps 
                                                        dalphaios <<- cnca_res$dalphaios 
                                                        inertia <<-cnca_res$inertia
                                                        descom <<- cnca_res$descom
                                                        inerciatot <<- cnca_res$inerciatot
                                                        
                                                }else{
                                                        if(cccaVal=="CCA")
                                                        {
                                                                cca_res <- ex.cca(fespecies, fvambientales, nejes, Numbercuant, Numbermixto, tChoice)
                                                                colugares <<- cca_res$colugares
                                                                coespecies <<- cca_res$coespecies
                                                                covambien <<- cca_res$covambien
                                                                colugaresnam <<-cca_res$colugaresnam
                                                                covariablesnam <<- cca_res$covariablesnam
                                                                coindividuosnam <<-cca_res$coindividuosnam
                                                                ejes <<- cca_res$ejes
                                                                Calphai <<- cca_res$Calphai 
                                                                Calphak <<- cca_res$Calphak
                                                                calidadps <<- cca_res$calidadps
                                                                calidados <<- cca_res$calidados
                                                                qalphaips <<- cca_res$qalphaips
                                                                qalphakos <<- cca_res$qalphakos
                                                                qalphakps <<- cca_res$qalphakps 
                                                                dalphaios <<- cca_res$dalphaios 
                                                                inertia <<-cca_res$inertia
                                                                descom <<- cca_res$descom
                                                                inerciatot <<- cca_res$inerciatot
                                                                
                                                        }else{
                                                                coia_res <- ex.coia(fespecies, fvambientales, nejes, tChoice,Numbercuant)
                                                                colugaresz <<- coia_res$colugaresz
                                                                colugaresy <<- coia_res$colugaresy
                                                                coespecies <<- coia_res$coespecies
                                                                covambien <<- coia_res$covambien
                                                                colugaresnamz <<-coia_res$colugaresnamz
                                                                colugaresnamy <<-coia_res$colugaresnamy
                                                                covariablesnam <<- coia_res$covariablesnam
                                                                coindividuosnam <<-coia_res$coindividuosnam
                                                                ejes <<- coia_res$ejes
                                                                inerciatot <<- coia_res$inerciatot
                                                                Xaxes <<- coia_res$Xaxes
                                                                Yaxes <<- coia_res$Yaxes
                                                        }
                                                        
                                                }
                                                
                                                
                                                cat("File saved in:    ",file="Results.txt")
                                                cat(getwd(),file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                cat("\n",file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                cat("\n",file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                if(cccaVal!="COIA")
                                                {
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("Proportion of inertia explained by each axis (projected space):\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(calidadps, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("Proportion of inertia explained by each axis (original space):\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(calidados, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("CONTRIBUTIONS:\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat(paste("Contributions of", Namesit,"to each factor:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(Calphai, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste("Contributions of", Namespe,"to each factor:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(Calphak, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("Relative contributions of factors to elements:\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namesit,"respect to projected space:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(qalphaips, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namesit,"respect to original space:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(dalphaios, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namespe,"respect to projected space:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(qalphakps, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namespe,"respect to original space:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(qalphakos, digits=3),file="temp.txt", sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namespe,"coordinates:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(coindividuosnam, digits=3),file="temp.txt",sep="\t",dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namevar,"coordinates:\n"),file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(covariablesnam, digits=3), file="temp.txt", sep="\t", dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                        
                                                        cat("\n",file="temp.txt")
                                                        file.append("Results.txt","temp.txt")
                                                        cat(paste(Namesit,"coordinates:\n"),file="temp.txt")        
                                                        file.append("Results.txt","temp.txt")
                                                        write.table(round(colugaresnam, digits=3), file="temp.txt", sep="\t", dec=",")
                                                        file.append("Results.txt","temp.txt")
                                                }
                                                
                                                
                                                cat("\n",file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                cat("Eigen values: \n",file="temp.txt")
                                                file.append("Results.txt","temp.txt")
                                                write.table(round(iner[[1]]$d, digits=3), file="temp.txt", sep="\t", dec=",")
                                                file.append("Results.txt","temp.txt")
                                                
                                                file.show("Results.txt")
                                                file.remove("temp.txt")
                                                
                                                ###########################################################################################
                                                ##Rescale
                                                ###########################################################################################
                                                
                                                if(cccaVal!="COIA")
                                                {
                                                        sumalugares <- sum(colugares^2)
                                                        sumaespecies <- sum(coespecies^2)
                                                        sumavariables <- sum(covambien^2)
                                                        
                                                        #slg <- sumalugares/(dim(colugares)[1])
                                                        sesp <- sumaespecies/(dim(coespecies)[1])
                                                        svar <- sumavariables/(dim(covambien)[1])
                                                        
                                                        scfev <- ((svar/sesp)^(1/2))^(1/2)
                                                        #scflv <- ((slg/sesp)^(1/2))^(1/2)
                                                        
                                                        coespecies <- coespecies*scfev
                                                        covambien <- covambien/scfev
                                                        #colugares <- (colugares/scflv)*scfev
                                                        datos <- rbind(coespecies,covambien,colugares)
                                                        limix <- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                        limiy <- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                        textos <- datos
                                                        datosr <- rbind(coespecies,covambien)
                                                        textosr <- datosr
                                                        simbolos <<- c(simespecies, simvariables, simlugares)
                                                        
                                                }else{
                                                        
                                                        if(typecoia=="sg") 
                                                        {
                                                                datos <<-rbind(colugaresy, colugaresz)
                                                                textos <<-colugaresy
                                                                for (d in 1:ejes)
                                                                {
                                                                        dimension <- cbind(colugaresy[,d], colugaresz[,d])
                                                                        textos[,d]<<-rowMeans(dimension)
                                                                }        
                                                        }
                                                        
                                                        if(typecoia=="eg")
                                                        {
                                                                datos<<-coespecies
                                                                textos<<-datos
                                                        }                
                                                        
                                                        if(typecoia=="vg")
                                                        {
                                                                datos<<-covambien
                                                                textos<<-datos
                                                        }
                                                        
                                                        limix <- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                        limiy <- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                        
                                                }
                                                
                                                
                                                centro <- c(0,0)
                                                
                                                ################ Show axes or not
                                                cbVal <<- as.character(tclvalue(cbValue))
                                                if(cccaVal!="COIA")
                                                {
                                                        cblugares <<- as.character(tclvalue(cblug))
                                                        cblablugares <<- as.character(tclvalue(cblablug))
                                                }
                                                
                                                
                                                
                                                wgr <- tktoplevel()
                                                tkwm.title(wgr,"Graph")
                                                parPlotSize <- c()
                                                usrCoords <- c()
                                                
                                                
                                                normalLine <- function(variableele, colorele, coordenadas, colorcoor) 
                                                { 
                                                        A <- variableele
                                                        B <- centro
                                                        
                                                        slopeAB <- (B[[2]] - A[[2]])/(B[[1]] - A[[1]]) 
                                                        slopeNorm <- -1/slopeAB 
                                                        a <- A[[2]] - slopeAB * A[[1]]
                                                        
                                                        b<-c()
                                                        xintersect<-c()
                                                        yintersect<-c()
                                                        
                                                        for(i in 1:dim(coordenadas)[1])
                                                        {
                                                                b[i] <- coordenadas[i,2] - slopeNorm * coordenadas[i,1] 
                                                                
                                                                xintersect[i] <- (b[i] - a)/(slopeAB - slopeNorm) 
                                                                yintersect[i] <- b[i] + slopeNorm * xintersect[i] 
                                                                
                                                        }#end for (i in 1:dim(coordenadas)[1])
                                                        
                                                        abline(a =a,b=slopeAB, col=colorele, lwd=3)
                                                        
                                                        for(i in 1:dim(coordenadas)[1])
                                                                segments(xintersect[i], yintersect[i], coordenadas[i,1], coordenadas[i,2],lty=2, col=colorcoor[i]) 
                                                } #end normalLine
                                                
                                                plotFunction <- function(screen=TRUE)
                                                {  
                                                        tclvalue(Limix1) <- limix[1]
                                                        tclvalue(Limix2) <- limix[2]
                                                        tclvalue(Limiy1) <- limiy[1]
                                                        tclvalue(Limiy2) <- limiy[2]
                                                        
                                                        if(cccaVal!="COIA")
                                                        {
                                                                ################ Show sites or not
                                                                if ((cblugares=="1")|(cblablugares=="1")){
                                                                        xCoords <<- textos[,dim1]
                                                                        yCoords <<- textos[,dim2]
                                                                        if(nejes>2)
                                                                                zCoords <<- textos[,dim3]
                                                                        
                                                                        if (cblablugares=="1"){
                                                                                labelsVec <<- c(textespecies, textvariables,textlugares)
                                                                                sizesVec <<- c(cexespecies, cexvariables,cexlugares)
                                                                        }else{
                                                                                labelsVec <<- c(textespecies, textvariables)
                                                                                sizesVec <<- c(cexespecies, cexvariables)
                                                                        }#end if (cblablugares=="1")
                                                                        
                                                                        colores <<- c(colespecies,colvariables,collugares)
                                                                }else{
                                                                        xCoords <<- textosr[,dim1]
                                                                        yCoords <<- textosr[,dim2]
                                                                        if(nejes>2)
                                                                                zCoords <<- textosr[,dim3]
                                                                        labelsVec <<- c(textespecies, textvariables)
                                                                        sizesVec <<- c(cexespecies, cexvariables)
                                                                        colores <<- c(colespecies,colvariables)
                                                                }#end if (cblugares=="1")
                                                                
                                                                if (is.null(indexLabeled))
                                                                        indexLabeled <<- c(1:length(xCoords))
                                                                indexLabeledaux <- c()
                                                                labeledPoints <- list()
                                                                params <- par(bg="white")
                                                                if ((cblugares=="1")|(cblablugares=="1")){
                                                                        xCoords<<-textos[,dim1]
                                                                        yCoords<<-textos[,dim2]
                                                                        plot(datos[,c(dim1,dim2)], main= tit_graph, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limix * 1.05, ylim = limiy * 1.05)
                                                                }else{
                                                                        xCoords<<-textosr[,dim1]
                                                                        yCoords<<-textosr[,dim2]
                                                                        plot(datosr[,c(dim1,dim2)],main= tit_graph, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limix * 1.05,  ylim = limiy * 1.05)
                                                                }#end if (cblugares=="1")
                                                                
                                                                if (proj=="normal")
                                                                {
                                                                        points(coespecies[,dim1],coespecies[,dim2],pch=simespecies, col=colespecies)
                                                                        arrows(centro[1],centro[2],covambien[1:Numbercuant,dim1],covambien[1:Numbercuant,dim2],col=colvariables[1:Numbercuant],lty=1, length=0.08)
                                                                         
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        if ((cblugares=="1")|(cblablugares=="1")){
                                                                                points(colugares[,dim1],colugares[,dim2],pch=simlugares, col=collugares)
                                                                        }#end if (cblugares=="1")
                                                                        
                                                                        if (length(indexLabeled)>0)
                                                                        {
                                                                                for (i in 1:length(indexLabeled))
                                                                                {
                                                                                        indexClosest <- indexLabeled[i]
                                                                                        text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= colores[indexClosest], cex= as.numeric(sizesVec[indexClosest]))
                                                                                }#end for
                                                                                
                                                                                if(Numbermixto>0)
                                                                                {
                                                                                        rangosup<-length(textespecies)+Numbercuant+1
                                                                                        rangoinf<-length(textespecies)+length(textvariables)
                                                                                        for(i in rangoinf:rangosup)
                                                                                        {
                                                                                                boxed.labels(xCoords[i],yCoords[i],labels=labelsVec[i], col= colores[i], cex= as.numeric(sizesVec)[i],border = colores[i])
                                                                                        }#end for
                                                                                }#end if
                                                                        }#end if
                                                                        
                                                                }else{
                                                                        colvariablesp[Choiceproj] <-colvariables[Choiceproj]
                                                                        if(proj=="e")
                                                                        {       
                                                                                indexLabeledp<-c(1:(length(textespecies)+length(textvariables)))
                                                                                coloresp<-c(colespecies, colvariablesp)
                                                                                points(coespecies[,dim1],coespecies[,dim2],pch=simespecies, col=colespecies)
                                                                                points(covambien[1:Numbercuant,dim1],covambien[1:Numbercuant,dim2],col=colvariablesp[1:Numbercuant])
                                                                                if(cbVal=="1")
                                                                                {
                                                                                        abline(h=centro[2],v=centro[1],lty="dotted")
                                                                                }
                                                                                if (length(indexLabeledp)>0)
                                                                                {
                                                                                        for (i in 1:length(indexLabeledp))
                                                                                        {
                                                                                                indexClosest <- indexLabeledp[i]
                                                                                                text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= coloresp[indexClosest], cex= as.numeric(sizesVec)[indexClosest])
                                                                                        }#end for
                                                                                        
                                                                                        if(Numbermixto>0)
                                                                                        {
                                                                                                rangosup<-length(textespecies)+Numbercuant+1
                                                                                                rangoinf<-length(textespecies)+length(textvariables)
                                                                                                for(i in rangoinf:rangosup)
                                                                                                {
                                                                                                        boxed.labels(xCoords[i],yCoords[i],labels=labelsVec[i], col= coloresp[i], cex= as.numeric(sizesVec)[i],border = coloresp[i])
                                                                                                }#end for
                                                                                        }#end if
                                                                                }#end if
                                                                                
                                                                                normalLine(covambien[Choiceproj,c(dim1,dim2)], colvariables[Choiceproj], coespecies[,c(dim1,dim2)], colespecies)
                                                                                
                                                                        }else{
                                                                                if(cblugares=="1" | cblablugares=="1")
                                                                                {
                                                                                        indexLabeledp<-c((1+length(textespecies)):(length(textespecies)+length(textvariables)+length(textlugares)))
                                                                                        coloresp<-c(colespecies, colvariablesp, collugares)
                                                                                        points(colugares[,dim1],colugares[,dim2],pch=simlugares, col=collugares)
                                                                                        points(covambien[1:Numbercuant,dim1],covambien[1:Numbercuant,dim2],col=colvariablesp[1:Numbercuant])
                                                                                        if(cbVal=="1")
                                                                                        {
                                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                                        }
                                                                                        if (length(indexLabeledp)>0)
                                                                                        {
                                                                                                for (i in (1:length(indexLabeledp)))
                                                                                                {
                                                                                                        indexClosest <- indexLabeledp[i]
                                                                                                        text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= coloresp[indexClosest], cex= as.numeric(sizesVec)[indexClosest])
                                                                                                }#end for
                                                                                                
                                                                                                if(Numbermixto>0)
                                                                                                {
                                                                                                        rangosup<-length(textespecies)+Numbercuant+1
                                                                                                        rangoinf<-length(textespecies)+length(textvariables)
                                                                                                        for(i in rangoinf:rangosup)
                                                                                                        {
                                                                                                                boxed.labels(xCoords[i],yCoords[i],labels=labelsVec[i], col= coloresp[i], cex= as.numeric(sizesVec)[i],border = coloresp[i])
                                                                                                        }#end for
                                                                                                }#end if
                                                                                        }#end if
                                                                                        
                                                                                        normalLine(covambien[Choiceproj,c(dim1,dim2)], colvariables[Choiceproj], colugares[,c(dim1,dim2)], collugares)
                                                                                        
                                                                                }else{
                                                                                        proj<<-"normal" 
                                                                                        
                                                                                        points(coespecies[,dim1],coespecies[,dim2],pch=simespecies, col=colespecies)
                                                                                        arrows(centro[1],centro[2],covambien[1:Numbercuant,dim1],covambien[1:Numbercuant,dim2],col=colvariables[1:Numbercuant],lty=1, length=0.08)
                                                                                        
                                                                                        if (cbVal=="1"){
                                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                                        }#end if (cbVal=="1")
                                                                                        if ((cblugares=="1")|(cblablugares=="1")){
                                                                                                points(colugares[,dim1],colugares[,dim2],pch=simlugares, col=collugares)
                                                                                        }#end if (cblugares=="1")
                                                                                        
                                                                                        if (length(indexLabeled)>0)
                                                                                        {
                                                                                                for (i in 1:length(indexLabeled))
                                                                                                {
                                                                                                        indexClosest <- indexLabeled[i]
                                                                                                        text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= colores[indexClosest], cex= as.numeric(sizesVec[indexClosest]))
                                                                                                }#end for
                                                                                                
                                                                                                if(Numbermixto>0)
                                                                                                {
                                                                                                        rangosup<-length(textespecies)
                                                                                                        rangoinf<-Numbercuant+1
                                                                                                        for(i in rangoinf:rangosup)
                                                                                                        {
                                                                                                                boxed.labels(xCoords[i],yCoords[i],labels=labelsVec[i], col= colores[i], cex= as.numeric(sizesVec)[i],border = colores[i])
                                                                                                        }#end for
                                                                                                }#end if
                                                                                        }#end if
                                                                                }#end if
                                                                        }#end if (proj=="e")
                                                                        
                                                                }#end if (proj=="normal")
                                                                
                                                        }else{
                                                                
                                                                xCoords <<- textos[,dim1]
                                                                yCoords <<- textos[,dim2]
                                                                if(nejes>2)
                                                                        zCoords <<- textos[,dim3]
                                                                if (is.null(indexLabeled))
                                                                        indexLabeled <<- c(1:length(xCoords))
                                                                indexLabeledaux <- c()
                                                                labeledPoints <- list()
                                                                params <- par(bg="white")
                                                                plot(datos[,c(dim1,dim2)], main= tit_graph, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limix * 1.05,  ylim = limiy * 1.05)
                                                                if (cbVal=="1"){
                                                                        abline(h=centro[2],v=centro[1],lty="dotted")
                                                                }#end if (cbVal=="1")
                                                                
                                                                if((typecoia=="eg"))
                                                                {
                                                                        labelsVec <<- textespecies
                                                                        sizesVec <<- cexespecies
                                                                        colores <<- colespecies
                                                                        arrows(centro[1],centro[2],datos[,dim1],datos[,dim2],col=colores,lty=1, length=0.08) 
                                                                }
                                                                
                                                                if((typecoia=="vg"))
                                                                {
                                                                        labelsVec <<- textvariables
                                                                        sizesVec <<- cexvariables
                                                                        colores <<- colvariables
                                                                        arrows(centro[1],centro[2],datos[,dim1],datos[,dim2],col=colores,lty=1, length=0.08) 
                                                                }
                                                                if(typecoia=="sg")
                                                                {
                                                                        labelsVec <<- textlugares
                                                                        sizesVec <<- cexlugares
                                                                        colores <<- collugares
                                                                        colugyaux <<- datos[1:length(textlugares),]
                                                                        colugzaux <<- datos[-(1:length(textlugares)),]
                                                                        arrows(colugzaux[,dim1],colugzaux[,dim2],colugyaux[,dim1],colugyaux[,dim2],col=colores,lty=1, length=0.08)
                                                                        
                                                                }
                                                                
                                                                if(typecoia=="ef")
                                                                {
                                                                        labelsVec <<- paste("Axis", 1:nejes)
                                                                        cexejes<<-rep(1,length(ejes))
                                                                        sizesVec <<- cexejes
                                                                        xCoords<<-datos[,dim1]
                                                                        yCoords<<-datos[,dim2]
                                                                        colores <<- rep("black", times=nejes)
                                                                        plot(0,0, type="n", xlim=c(-1,1), ylim=c(-1,1), asp=1/1, main=paste(Namespe, " axes"), xlab="", ylab="")
                                                                        draw.circle(0,0,1, border="black", lty=1,lwd=1)
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        arrows(centro[1],centro[2],datos[,dim1],datos[,dim2],lty=1, length=0.08)
                                                                }
                                                                
                                                                if(typecoia=="vf")
                                                                {
                                                                        labelsVec <<- paste("Axis", 1:nejes)
                                                                        cexejes<<-rep(1,length(ejes))
                                                                        sizesVec <<- cexejes
                                                                        xCoords<<-datos[,dim1]
                                                                        yCoords<<-datos[,dim2]
                                                                        colores <<- rep("black", times=nejes)
                                                                        plot(0,0, type="n", xlim=c(-1,1), ylim=c(-1,1), asp=1/1, main=paste(Namespe, " axes"), xlab="", ylab="")
                                                                        draw.circle(0,0,1, border="black", lty=1,lwd=1)
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        arrows(centro[1],centro[2],datos[,dim1],datos[,dim2],lty=1, length=0.08)
                                                                }
                                                                
                                                                if (length(indexLabeled)>0)
                                                                {
                                                                        for (i in 1:length(indexLabeled))
                                                                        {
                                                                                indexClosest <- indexLabeled[i]
                                                                                if(typecoia=="sg")
                                                                                {
                                                                                        boxed.labels(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= colores[indexClosest], cex= as.numeric(sizesVec[indexClosest]),border = colores[indexClosest])
                                                                                }else{
                                                                                        if(typecoia%in% c("ef", "vf"))
                                                                                        {                                              
                                                                                                text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= colores[indexClosest], cex= as.numeric(sizesVec[indexClosest]), pos=1)                
                                                                                        }else{
                                                                                                text(xCoords[indexClosest],yCoords[indexClosest],labels=labelsVec[indexClosest], col= colores[indexClosest], cex= as.numeric(sizesVec[indexClosest]))                
                                                                                        }
                                                                                }
                                                                        }#end for          
                                                                }#end if
                                                                
                                                        }      
                                                        
                                                        parPlotSize <<- par("plt")
                                                        usrCoords   <<- par("usr")
                                                        par(params)   
                                                        
                                                }#end plotFunction <- function(screen=TRUE)
                                                
                                                img <- tkrplot(wgr,fun=plotFunction,hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                                                framedim1<-tkframe(wgr, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                framecomb<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                framescal<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                framescala<-tkframe(framescal, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                framescalb<-tkframe(framescal, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                framerefr<-tkframe(framedim1, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                                
                                                comboBoxdim1 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=15, text= dim1, background="white")
                                                comboBoxdim2 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=15, text= dim2, background="white")
                                                if (nejes>2)
                                                        comboBoxdim3 <- tkwidget(framecomb,"ComboBox",editable=FALSE,values=rep(1:nejes),width=15, text= dim3, background="white")
                                                
                                                chang.symdim1 <- function()
                                                {
                                                        dim1 <<-as.numeric(tclvalue(tcl(comboBoxdim1,"getvalue")))+1
                                                        dim2 <<-as.numeric(tclvalue(tcl(comboBoxdim2,"getvalue")))+1
                                                        if (nejes>2)
                                                                dim3 <<-as.numeric(tclvalue(tcl(comboBoxdim3,"getvalue")))+1
                                                        if((dim1==dim1ant)&(dim2==dim2ant))
                                                        {
                                                                limix[1] <<- as.numeric(tclvalue(Limix1))
                                                                limix[2] <<- as.numeric(tclvalue(Limix2))
                                                                limiy[1] <<- as.numeric(tclvalue(Limiy1))
                                                                limiy[2] <<- as.numeric(tclvalue(Limiy2))
                                                        }else{
                                                                dim1ant<<-dim1
                                                                dim2ant<<-dim2
                                                                if (((cblugares=="1")|(cblablugares=="1"))|cccaVal=="COIA")
                                                                {
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                }else{
                                                                        limix <<- round(c(min(datosr[,dim1],0), max(datosr[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datosr[,dim2],0), max(datosr[,dim2],0)), digits=2)
                                                                }
                                                                
                                                        }
                                                        tkrreplot(img)
                                                }#end chang.symdim1 <- function()
                                                
                                                Change.symboldim1 <-tkbutton(framerefr,text="Refresh",command=chang.symdim1, bg= "lightblue", width=15, foreground = "navyblue")
                                                if(nejes>2){
                                                        tkpack(tklabel(framecomb, text="Select X, Y and Z axes numbers:"), expand="FALSE", side= "top", fill ="both")
                                                        tkpack(comboBoxdim1, comboBoxdim2, comboBoxdim3, side="top", fill="x")
                                                             }
                                                Limix1 <- tclVar(limix[1])
                                                entry.limix1 <<-tkentry(framescalb,width=10,textvariable=Limix1, bg="white")
                                                tkconfigure(entry.limix1,textvariable=Limix1)          	
                                                tkbind(entry.limix1, "<Return>",chang.symdim1)
                                                
                                                Limix2 <- tclVar(limix[2])
                                                entry.limix2 <<-tkentry(framescalb,width=10,textvariable=Limix2, bg="white")
                                                tkconfigure(entry.limix2,textvariable=Limix2)                  
                                                tkbind(entry.limix2, "<Return>",chang.symdim1)
                                                
                                                Limiy1 <- tclVar(limiy[1])
                                                entry.limiy1 <<-tkentry(framescalb,width=10,textvariable=Limiy1, bg="white")
                                                tkconfigure(entry.limiy1,textvariable=Limiy1)                  
                                                tkbind(entry.limiy1, "<Return>",chang.symdim1)
                                                
                                                Limiy2 <- tclVar(limiy[2])
                                                entry.limiy2 <<-tkentry(framescalb,width=10,textvariable=Limiy2, bg="white")
                                                tkconfigure(entry.limiy2,textvariable=Limiy2) 
                                                tkbind(entry.limiy2, "<Return>",chang.symdim1)
                                                tkpack(tklabel(framescala, text="-X"),
                                                       tklabel(framescala, text="+X"),
                                                       tklabel(framescala, text="-Y"),
                                                       tklabel(framescala, text="+Y"),
                                                       expand = "TRUE",side="top", fill="both")
                                                tkpack(entry.limix1,entry.limix2,entry.limiy1,entry.limiy2,side="top", expand="TRUE", fill="both")
                                                tkpack(tklabel(framescal, text="Zoom:"), expand = "TRUE",side="top", fill="both")
                                                tkpack(framescala, framescalb, side="left", expand ="TRUE", fill="both")
                                                tkpack(Change.symboldim1, side="top", expand="TRUE", fill="both")
                                                tkpack(framecomb, framescal, framerefr, side="top", expand="TRUE", fill="both")
                                                tkpack(framedim1, side="left", expand="FALSE", fill="x")
                                                tkpack(img, side="top", expand="TRUE", fill="both")
                                                
                                                
                                                labelClosestPointd <- function(xClick,yClick,imgXcoords,imgYcoords)
                                                {
                                                        squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                        indexClosest <<- which.min(squared.Distance)
                                                        mm<-tktoplevel() 
                                                        tkwm.title(mm, labelsVec[indexClosest])
                                                        framemm1<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                        if(indexClosest<=length(labelsVec))
                                                        {
                                                                framemm2<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                                framemm3<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")
                                                        }
                                                        framemm4<-tkframe(mm, relief = "groove", borderwidth = 2, background = "white")                                                                
                                                        
                                                        
                                                        colori <- colores[indexClosest]
                                                        canvasi <- tkcanvas(framemm1,width="70",height="20",bg=colori)
                                                        
                                                        ChangeColori <- function()
                                                        {
                                                                colori <<- tclvalue(tcl("tk_chooseColor",initialcolor=colores[indexClosest],title="Choose a color"))
                                                                if (nchar(colori)>0)
                                                                {
                                                                        tkconfigure(canvasi,bg=colori)
                                                                        colores[indexClosest]<<-colori
                                                                        if(cccaVal!="COIA")
                                                                        {
                                                                                if (cblugares=="1"){
                                                                                        colespecies<<-colores[1:length(colespecies)]
                                                                                        colvariables<<-colores[(length(colespecies)+1):(length(colespecies)+length(colvariables))]
                                                                                        collugares<<-colores[(length(colespecies)+length(colvariables)+1):(length(colespecies)+length(colvariables)+length(collugares))]
                                                                                }else{
                                                                                        colespecies<<-colores[1:length(colespecies)]
                                                                                        colvariables<<-colores[(length(colespecies)+1):(length(colespecies)+length(colvariables))]
                                                                                }#end if (cblugares=="1")        
                                                                        }else{
                                                                                if(typecoia=="sg")
                                                                                {
                                                                                        collugares<<-colores                                                                                        
                                                                                }else{
                                                                                        if((typecoia=="eg")|(typecoia=="ef"))
                                                                                        {
                                                                                                colespecies<<-colores
                                                                                        }else{
                                                                                                if((typecoia=="vg")|(typecoia=="vf"))
                                                                                                {
                                                                                                        colvariables<<-colores
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                        
                                                                }#end if (nchar(colori)>0)
                                                                
                                                                tkrreplot(img)
                                                                tkdestroy(mm)
                                                        }#end ChangeColori <- function()
                                                        
                                                        ChangeColor.buttoni<- tkbutton(framemm1,text="Change Color",command=ChangeColori,width=12)
                                                        tkpack(canvasi,ChangeColor.buttoni,expand = "TRUE", side="left", fill = "both")
                                                        if(indexClosest<=length(labelsVec))
                                                        {
                                                                Namei <<- labelsVec[indexClosest]
                                                                tclvalue(Namei) <<- labelsVec[indexClosest]
                                                                entry.Namei <<-tkentry(framemm2,width="11",textvariable=Namei, bg="white")
                                                                NameVali <<- Namei 
                                                                
                                                                OnOKli <- function()
                                                                {
                                                                        NameVali <<- tclvalue(Namei)
                                                                        labelsVec[indexClosest]<<-NameVali
                                                                        
                                                                        if(cccaVal!="COIA")
                                                                        {
                                                                                if (cblablugares=="1"){
                                                                                        textespecies <<-labelsVec[1:length(textespecies)]
                                                                                        textvariables <<-labelsVec[(length(textespecies)+1):(length(textespecies)+length(textvariables))]
                                                                                        textlugares <<-labelsVec[(length(textespecies)+length(textvariables)+1):(length(textespecies)+length(textvariables)+length(collugares))]
                                                                                }else{
                                                                                        textespecies <<-labelsVec[1:length(textespecies)]
                                                                                        textvariables <<-labelsVec[(length(textespecies)+1):(length(textespecies)+length(textvariables))]
                                                                                }#end if (cblugares=="1")
                                                                        }else{
                                                                                if(typecoia=="sg")
                                                                                {
                                                                                        textlugares<<-labelsVec                                                                                        
                                                                                }else{
                                                                                        if((typecoia=="eg")|(typecoia=="ef"))
                                                                                        {
                                                                                                textespecies<<-labelsVec
                                                                                        }else{
                                                                                                if((typecoia=="vg")|(typecoia=="vf"))
                                                                                                {
                                                                                                        textvariables<<-labelsVec
                                                                                                }
                                                                                        }
                                                                                }
                                                                                
                                                                        }
                                                                        
                                                                        tkrreplot(img)
                                                                        tkdestroy(mm)
                                                                }#end OnOKli <- function()
                                                                
                                                                OK.butli <-tkbutton(framemm2,text="Change label",command=OnOKli,width=12)
                                                                tkbind(entry.Namei, "<Return>",OnOKli)
                                                                tkpack(entry.Namei,OK.butli,expand = "TRUE", side="left", fill = "both")
                                                                
                                                                Cexi <<- sizesVec[indexClosest]
                                                                
                                                                tclvalue(Cexi) <<- sizesVec[indexClosest]
                                                                entry.Cexi <<-tkentry(framemm3,width="11",textvariable=Cexi, bg="white")
                                                                NameCexi <<- Cexi 
                                                                
                                                                OnOKci <- function()
                                                                {
                                                                        NameCexi <<- tclvalue(Cexi)
                                                                        sizesVec[indexClosest]<<-NameCexi
                                                                        if(cccaVal!="COIA")
                                                                        {
                                                                                if (cblablugares=="1"){
                                                                                        cexespecies <<-sizesVec[1:length(cexespecies)]
                                                                                        cexvariables <<-sizesVec[(length(cexespecies)+1):(length(cexespecies)+length(cexvariables))]
                                                                                        cexlugares <<-sizesVec[(length(cexespecies)+length(cexvariables)+1):(length(cexespecies)+length(cexvariables)+length(cexlugares))]
                                                                                }else{
                                                                                        cexespecies <<-sizesVec[1:length(cexespecies)]
                                                                                        cexvariables <<-sizesVec[(length(cexespecies)+1):(length(cexespecies)+length(cexvariables))]
                                                                                }#end if (cblugares=="1")
                                                                                
                                                                        }else{
                                                                                if(typecoia=="sg")
                                                                                {
                                                                                        cexlugares<<-sizesVec                                                                                        
                                                                                }else{
                                                                                        if((typecoia=="eg")|(typecoia=="ef"))
                                                                                        {
                                                                                                cexespecies<<-sizesVec
                                                                                        }else{
                                                                                                if((typecoia=="vg")|(typecoia=="vf"))
                                                                                                {
                                                                                                        cexvariables<<-sizesVec
                                                                                                }
                                                                                        }
                                                                                }
                                                                        }
                                                                        
                                                                        tkrreplot(img)        
                                                                        tkdestroy(mm)
                                                                }#end OnOKci <- function()
                                                                
                                                                OK.butci <-tkbutton(framemm3,text="Change size",command=OnOKci,width=12)
                                                                tkbind(entry.Cexi, "<Return>",OnOKci)
                                                                tkpack(entry.Cexi,OK.butci,expand = "TRUE", side="left", fill = "both")
                                                                
                                                        }        
                                                        
                                                        if(cccaVal!="COIA")
                                                                comboBox <- tkwidget(framemm4,"ComboBox",editable=FALSE,values=symbolos,width=8, text= simbolos[indexClosest], background="white")
                                                        
                                                        chang.sym <- function()
                                                        {
                                                                if(cccaVal!="COIA")
                                                                {
                                                                        simChoice <<-symbolos[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
                                                                        simbolos[indexClosest]<<-simChoice
                                                                        
                                                                        if (cblugares=="1"){
                                                                                simespecies <<-simbolos[1:length(simespecies)]
                                                                                simvariables <<-simbolos[(length(simespecies)+1):(length(simespecies)+length(simvariables))]
                                                                                simlugares <<-simbolos[(length(simespecies)+length(simvariables)+1):(length(simespecies)+length(simvariables)+length(collugares))]
                                                                        }else{
                                                                                simespecies <<-simbolos[1:length(simespecies)]
                                                                                simvariables <<-simbolos[(length(simespecies)+1):(length(simespecies)+length(simvariables))]
                                                                        }#end if (cblugares=="1")
                                                                        
                                                                        
                                                                }
                                                                
                                                                tkrreplot(img)
                                                                tkdestroy(mm)
                                                        }#end chang.sym <- function()
                                                        if(cccaVal!="COIA")
                                                        {
                                                                if(indexClosest %in% c(length(simespecies)+1):(length(simespecies)+length(simvariables)))
                                                                {}else{
                                                                        Change.symbol <-tkbutton(framemm4,text="Change symbol",command=chang.sym,width=12)
                                                                        tkpack(comboBox,Change.symbol,side="left",expand="TRUE", fill="both")        
                                                                }       
                                                        }
                                                        
                                                        
                                                        
                                                        if(proj=="normal")
                                                        {
                                                                if(indexClosest<= length(labelsVec))
                                                                {
                                                                        tkpack(framemm1,framemm2,framemm3,framemm4, expand = "TRUE", side="top", fill = "both")
                                                                }else{
                                                                        tkpack(framemm1,framemm4, expand = "TRUE", side="top", fill = "both")        
                                                                }
                                                        }else{
                                                                tkpack(framemm1,framemm2, expand = "TRUE", side="top", fill = "both")
                                                        }
                                                        if(cccaVal=="COIA")
                                                        {
                                                                tkpack(framemm1,framemm2,framemm3, expand = "TRUE", side="top", fill = "both")
                                                        }
                                                }#end labelClosestPointd <- function(xClick,yClick,imgXcoords,imgYcoords)
                                                
                                                
                                                #############################################################################
                                                ### Functions to save the graph
                                                #############################################################################
                                                SaveFileJPG <- function() {
                                                        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Jpeg files} {.jpg .jpeg}} {{All files} *}"))
                                                        if (nchar(FileName)) {
                                                                nn <- nchar(FileName)
                                                                if (nn < 5 || substr(FileName, nn - 3, nn) != ".jpg") 
                                                                        FileName <- paste(FileName, ".jpg", sep = "")
                                                                jpeg(FileName, width = 8, height = 8, units = "in", res = 96, quality = 100)
                                                                plotFunction(screen = FALSE)
                                                                dev.off()
                                                        }#end if (nchar(FileName))
                                                }#end SaveFileJPG <- function()
                                                
                                                SaveFilePDF <- function() {
                                                        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{PDF files} {.pdf}} {{All files} *}"))
                                                        if (nchar(FileName)) {
                                                                nn <- nchar(FileName)
                                                                if (nn < 5 || substr(FileName, nn - 3, nn) != ".pdf") 
                                                                        FileName <- paste(FileName, ".pdf", sep = "")
                                                                pdf(FileName, width = 7, height = 7, useDingbats=FALSE)
                                                                plotFunction(screen = FALSE)
                                                                dev.off()
                                                        }#end if (nchar(FileName)) 
                                                }#end SaveFilePDF <- function()
                                                
                                                SaveFileeps <- function() {
                                                        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Eps files} {.eps}} {{All files} *}"))
                                                        if (nchar(FileName)) {
                                                                nn <- nchar(FileName)
                                                                if (nn < 5 || substr(FileName, nn - 3, nn) != ".eps") 
                                                                        FileName <- paste(FileName, ".eps", sep = "")
                                                                postscript(FileName, width = 8, height = 8)
                                                                plotFunction(screen = FALSE)
                                                                dev.off()
                                                        }#end if (nchar(FileName))
                                                }#end SaveFilePng <- function()
                                                
                                                                                               
                                                SaveFilePng <- function() {
                                                        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Png files} {.png}} {{All files} *}"))
                                                        if (nchar(FileName)) {
                                                                nn <- nchar(FileName)
                                                                if (nn < 5 || substr(FileName, nn - 3, nn) != ".png") 
                                                                        FileName <- paste(FileName, ".png", sep = "")
                                                                png(FileName, width = 8, height = 8, units = "in", res = 96)
                                                                plotFunction(screen = FALSE)
                                                                dev.off()
                                                        }#end if (nchar(FileName))
                                                }#end SaveFilePng <- function() 
                                                
                                                g3d<-function()
                                                {
                                                        if (nejes>2)
                                                        { 
                                                                bg3d("white")
                                                                aspect3d("iso")
                                                                lims <- par3d("bbox")
                                                                
                                                                if (cbVal=="1"){
                                                                        axes3d()
                                                                }#end if (cbVal=="1")
                                                                
                                                                
                                                                if(cccaVal!="COIA")
                                                                {
                                                                        if (cblugares=="1"){
                                                                                zCoords<<-datos[,dim3]
                                                                                
                                                                                if(cblablugares=="1"){
                                                                                        labelsVec3d <-c(textespecies, textvariables, textlugares)
                                                                                        sizesVec3d <-c(cexespecies, cexvariables, cexlugares)
                                                                                }else{
                                                                                        labelsVec3d <-c(textespecies, textvariables, rep("", times =length(textlugares)))
                                                                                        sizesVec3d <-c(cexespecies, cexvariables, rep("", times =length(cexlugares)))
                                                                                }#end if(cblablugares=="1")
                                                                                
                                                                        }else{
                                                                                zCoords<<-datosr[,dim3]
                                                                                labelsVec3d <-c(textespecies, textvariables)
                                                                                sizesVec3d <-c(cexespecies, cexvariables)
                                                                        }#end if (cblugares=="1")
                                                                        points3d(xCoords,yCoords,zCoords, color=colores)
                                                                        texts3d(xCoords, yCoords, zCoords,labelsVec3d,color=colores, cex = as.numeric(sizesVec3d))
                                                                        
                                                                        for (i in 1:(dim(covambien)[1]))
                                                                        {
                                                                                linea <-rbind(covambien[i,c(dim1, dim2, dim3)],c(0,0,0))
                                                                                segments3d(linea[,1],linea[,2], linea[,3],color=colvariables[i])
                                                                        }#end for (i in 1:(dim(covambien)[1]))
                                                                        
                                                                }else{
                                                                        xCoords<<-datos[,dim1]
                                                                        yCoords<<-datos[,dim2]
                                                                        zCoords<<-datos[,dim3]
                                                                        
                                                                        labelsVec3d <-labelsVec
                                                                        sizesVec3d <-sizesVec
                                                                        if(typecoia=="sg")
                                                                        {
                                                                                points3d(xCoords,yCoords,zCoords, color=colores)
                                                                                colugyaux <<- datos[1:length(textlugares),]
                                                                                colugzaux <<- datos[-(1:length(textlugares)),]
                                                                                
                                                                                for (i in 1:(dim(colugyaux)[1]))
                                                                                {
                                                                                        linea <-rbind(colugyaux[i,c(dim1, dim2, dim3)],colugzaux[i,c(dim1, dim2, dim3)])
                                                                                        segments3d(linea[,1],linea[,2], linea[,3],color=colores[i])
                                                                                }#end for (i in 1:(dim(colugyaux)[1]))
                                                                                texts3d(textos[,1], textos[,2], textos[,3],labelsVec3d,color=colores, cex = as.numeric(sizesVec3d))        
                                                                        }
                                                                        if(typecoia=="eg" | typecoia=="vg")
                                                                        {
                                                                                for (i in 1:(dim(datos)[1]))
                                                                                {
                                                                                        linea <-rbind(datos[i,c(dim1, dim2, dim3)],c(0,0,0))
                                                                                        segments3d(linea[,1],linea[,2], linea[,3],color=colores[i])
                                                                                }#end for (i in 1:(dim(datos)[1]))
                                                                                texts3d(textos[,1], textos[,2], textos[,3],labelsVec3d,color=colores, cex = as.numeric(sizesVec3d))        
                                                                                
                                                                        }
                                                                        if(typecoia=="ef" | typecoia=="vf")
                                                                        {
                                                                                rgl.bg(sphere=TRUE, color=c("white","green"), lit=FALSE, back="lines" )
                                                                                for (i in 1:(dim(datos)[1]))
                                                                                {
                                                                                        linea <-rbind(datos[i,c(dim1, dim2, dim3)],c(0,0,0))
                                                                                        segments3d(linea[,1],linea[,2], linea[,3],color=colores[i])
                                                                                }#end for (i in 1:(dim(datos)[1]))
                                                                                texts3d(textos[,1], textos[,2], textos[,3],labelsVec3d,color=colores, cex = as.numeric(sizesVec3d))        
                                                                                
                                                                        }
                                                                        
                                                                }
                                                                rgl.bringtotop()
                                                        }else{
                                                                msg <- "You have selected less than 3 dimensions. 3D-graph not available"
                                                                tkmessageBox(message=msg)
                                                        }#end if (nejes>2)
                                                }#end g3d<-function()
                                                
                                                bootcnca<-function()
                                                {
                                                        wboot<-tktoplevel()
                                                        tkwm.title(wboot,"Bootstrap")
                                                        #### Frames
                                                        
                                                        framewi<-tkframe(wboot, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi1<-tkframe(framewi, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi2<-tkframe(framewi, relief = "ridge", borderwidth = 2, background = "white")
                                                        
                                                        framewi21<-tkframe(framewi2, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi21i<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi21a<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi21f<-tkframe(framewi21, relief = "ridge", borderwidth = 2, background = "white")
                                                        
                                                        framewi21ft<-tkframe(framewi21f, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi21fb<-tkframe(framewi21f, relief = "ridge", borderwidth = 2, background = "white")
                                                        
                                                        framewi21fl<-tkframe(framewi21fb, relief = "ridge", borderwidth = 2, background = "white")
                                                        framewi21fr<-tkframe(framewi21fb, relief = "ridge", borderwidth = 2, background = "white")
                                                        
                                                        framewi22<-tkframe(framewi2, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi22c<-tkframe(framewi22, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi22l<-tkframe(framewi22, relief = "flat", borderwidth = 2, background = "white")
                                                        
                                                        framewi221c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi222c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi223c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi224c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi225c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi226c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi227c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi228c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi229c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2210c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2211c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2212c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2213c<-tkframe(framewi22c, relief = "flat", borderwidth = 2, background = "white")
                                                        
                                                        framewi221l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi222l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi223l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi224l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi225l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi226l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")	
                                                        framewi227l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi228l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi229l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2210l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2211l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2212l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewi2213l<-tkframe(framewi22l, relief = "flat", borderwidth = 2, background = "white")
                                                        framewigr<-tkframe(wboot, relief = "flat", borderwidth = 2, background = "white")
                                                        
                                                        fontHeading <- tkfont.create(family="times",size=24,weight="bold",slant="italic")
                                                        fontFixedWidth <- tkfont.create(family="courier",size=12)
                                                        tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                        tkpack(tklabel(framewi1,text="BOOTSTRAP",font=fontHeading, foreground = "blue"), expand = "TRUE", side="left", fill = "both")
                                                        tkpack(tklabel(framewi1, text="    "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                        
                                                        
                                                        ######Iterations   ###################################
                                                        Niter <- tclVar(niter)
                                                        entry.Niter <-tkentry(framewi21i,width=10,textvariable=Niter, bg="white")
                                                        tkconfigure(entry.Niter,textvariable=Niter)  		
                                                        
                                                        tkpack(tklabel(framewi21i, text="Number of resamples"),entry.Niter, expand = "TRUE", side="left", fill = "both")
                                                        
                                                        ######alpha confidence intervals ###################################
                                                        Nalpha <- tclVar(alphaic)
                                                        entry.Nalpha <-tkentry(framewi21a,width=10,textvariable=Nalpha, bg="white")
                                                        tkconfigure(entry.Nalpha,textvariable=Nalpha)
                                                        
                                                        tkpack(tklabel(framewi21a, text="Confidence Level     "),entry.Nalpha, expand = "TRUE", side="left", fill = "both")
                                                        
                                                        
                                                        tkpack(framewi21fl, framewi21fr, expand = "TRUE",side="left", fill="both")
                                                        tkpack(framewi21ft, framewi21fb, expand = "TRUE",side="top", fill="both")
                                                        tkpack(framewi21i, framewi21a, framewi21f, expand = "TRUE",side="top", fill="both")
                                                        
                                                        
                                                        ###### Parameters to estimate   ###################################
                                                        tkpack(tklabel(framewi221l, text="Calculate confidence intervals for:"), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                        tkpack(tklabel(framewi221c, text=" "), expand = "TRUE", side="left",expand="TRUE", fill = "both")
                                                        
                                                        ##### Checkbox Contributions to namesit  #######
                                                        
                                                        cCalphai <- tkcheckbutton(framewi222c)
                                                        cCalphaiValue <- tclVar("0")
                                                        tkconfigure(cCalphai,variable=cCalphaiValue)
                                                        
                                                        ##### Checkbox Contributions to namespe  #######
                                                        
                                                        cCalphak <- tkcheckbutton(framewi223c)
                                                        cCalphakValue <- tclVar("0")
                                                        tkconfigure(cCalphak,variable=cCalphakValue)
                                                        
                                                        ##### Checkbox Proportion of inertia explained by each axis (projected space)  #######
                                                        
                                                        ccalidadps <- tkcheckbutton(framewi224c)
                                                        ccalidadpsValue <- tclVar("0")
                                                        tkconfigure(ccalidadps,variable=ccalidadpsValue)
                                                        
                                                        ##### Checkbox Proportion of inertia explained by each axis (original space)  #######
                                                        
                                                        ccalidados <- tkcheckbutton(framewi225c)
                                                        ccalidadosValue <- tclVar("0")
                                                        tkconfigure(ccalidados,variable=ccalidadosValue)
                                                        
                                                        ##### Checkbox Qualities of representation relative to Namesit respect to projected space  #######
                                                        
                                                        cqalphaips <- tkcheckbutton(framewi226c)
                                                        cqalphaipsValue <- tclVar("0")
                                                        tkconfigure(cqalphaips,variable=cqalphaipsValue)
                                                        
                                                        ##### Checkbox Qualities of representation relative to Namesit respect to original space #######
                                                        
                                                        cdalphaios <- tkcheckbutton(framewi227c)
                                                        cdalphaiosValue <- tclVar("0")
                                                        tkconfigure(cdalphaios,variable=cdalphaiosValue)
                                                        
                                                        ##### Checkbox Qualities of representation relative to Namespe respect to projected space #######
                                                        
                                                        cqalphakps <- tkcheckbutton(framewi228c)
                                                        cqalphakpsValue <- tclVar("0")
                                                        tkconfigure(cqalphakps,variable=cqalphakpsValue)
                                                        
                                                        ##### Checkbox Qualities of representation relative to Namespe respect to original space #######
                                                        
                                                        cqalphakos <- tkcheckbutton(framewi229c)
                                                        cqalphakosValue <- tclVar("0")
                                                        tkconfigure(cqalphakos,variable=cqalphakosValue)
                                                        
                                                        tkpack(tklabel(framewi21ft,text="Save files as:"),
                                                               expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        cpdf <- tkradiobutton(framewi21fl)
                                                        cpdfb <- tkradiobutton(framewi21fl)
                                                        rbpdfValue <- tclVar("Color pdf")
                                                        tkconfigure(cpdf,variable=rbpdfValue,value="Color pdf")
                                                        tkconfigure(cpdfb,variable=rbpdfValue,value="Black and white pdf")
                                                        tkpack(tklabel(framewi21fl,text="Color pdf"),cpdf,
                                                               expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        tkpack(tklabel(framewi21fl,text="Black and white pdf"),cpdfb,
                                                               expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        
                                                        
                                                        ceps <- tkradiobutton(framewi21fr)
                                                        cepsb <- tkradiobutton(framewi21fr)
                                                        rbepsValue <- tclVar("Color eps")
                                                        tkconfigure(ceps,variable=rbepsValue,value="Color eps")
                                                        tkconfigure(cepsb,variable=rbepsValue,value="Black and white eps")
                                                        tkpack(tklabel(framewi21fr,text="Color eps"),ceps,
                                                               expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        tkpack(tklabel(framewi21fr,text="Black and white eps"),cepsb,
                                                               expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        
                                                       
                                                        ##### Checkbox Eigenvalues #######
                                                        
                                                        ceigen <- tkcheckbutton(framewi2213c)
                                                        ceigenValue <- tclVar("0")
                                                        tkconfigure(ceigen,variable=ceigenValue)
                                                        
                                                        tkpack( tklabel(framewi222l, text=paste("-Contributions to factor of", Namesit), anchor="nw"), 
                                                                tklabel(framewi223l, text=paste("-Contributions to factor of", Namespe), anchor="nw"), 
                                                                tklabel(framewi224l, text="-Proportion of inertia explained (projected space)", anchor="nw"),
                                                                tklabel(framewi225l, text="-Proportion of inertia explained (original space)", anchor="nw"),
                                                                tklabel(framewi226l, text=paste("-Relative contribution of each factor to ", Namesit, "respect to projected space"), anchor="nw"),
                                                                tklabel(framewi227l, text=paste("-Relative contribution of each factor to ", Namesit, "respect to original space"), anchor="nw"),
                                                                tklabel(framewi228l, text=paste("-Relative contribution of each factor to ", Namespe, "respect to projected space"), anchor="nw"),
                                                                tklabel(framewi229l, text=paste("-Relative contribution of each factor to ", Namespe, "respect to original space"), anchor="nw"),
                                                                tklabel(framewi2213l, text="-Eigenvalues", anchor="nw"),
                                                                expand = "FALSE", side="top",expand="TRUE", fill = "both")
                                                        
                                                        tkpack(cCalphai,cCalphak,ccalidadps,ccalidados,cqalphaips,cdalphaios,cqalphakps,cqalphakos,
                                                               ceigen,
                                                               expand = "TRUE",side="top", fill="both")
                                                        tkpack(framewi221l,framewi222l,framewi223l,framewi224l,framewi225l,framewi226l,framewi227l,framewi228l,framewi229l,
                                                               framewi2213l, expand = "TRUE",side="top", fill="both")
                                                        tkpack(framewi221c,framewi222c,framewi223c,framewi224c,framewi225c,framewi226c,framewi227c,framewi228c,framewi229c,
                                                               framewi2213c, expand = "TRUE",side="top", fill="both")
                                                        tkpack(framewi22c, framewi22l, expand = "TRUE",side="left", fill="both")
                                                        
                                                        OnOKboot<-function()
                                                        {
                                                                tkdestroy(wboot) 
                                                                niter <<- tclvalue(Niter)
                                                                alphaic <<- tclvalue(Nalpha)
                                                                cCalphaiVal <<- as.character(tclvalue(cCalphaiValue))
                                                                cCalphakVal <<- as.character(tclvalue(cCalphakValue))
                                                                ccalidadpsVal <<- as.character(tclvalue(ccalidadpsValue))
                                                                ccalidadosVal <<- as.character(tclvalue(ccalidadosValue))
                                                                cqalphaipsVal <<- as.character(tclvalue(cqalphaipsValue))
                                                                cdalphaiosVal <<- as.character(tclvalue(cdalphaiosValue))
                                                                cqalphakpsVal <<- as.character(tclvalue(cqalphakpsValue))
                                                                cqalphakosVal <<- as.character(tclvalue(cqalphakosValue))
                                                                ceigenVal <<- as.character(tclvalue(ceigenValue))
                                                                
                                                                cpdfVal <<- as.character(tclvalue(rbpdfValue))
                                                                cepsVal <<- as.character(tclvalue(rbepsValue))
                                                                
                                                                
                                                                
                                                                # jakknife submuestra
                                                                
                                                                jakespecies <- vector("list",(sum(fespecies)-1))
                                                                z<-1
                                                                for(i in 1:dim(fespecies)[1])
                                                                {
                                                                        for (j in 1:dim(fespecies)[2])
                                                                        {
                                                                                if(fespecies[i,j]>0)
                                                                                {
                                                                                        for (k in 1:(fespecies[i,j]))
                                                                                        {
                                                                                                especaux<-as.matrix(fespecies)
                                                                                                especaux[i,j]<-fespecies[i,j]-1
                                                                                                jakespecies[[z]]<-especaux
                                                                                                z<-z+1    
                                                                                        }
                                                                                }
                                                                        }
                                                                }
                                                                
                                                                ## muestra bootstrap
                                                                especiesrep<-rep(list(fespecies),niter)
                                                                especresample<-lapply(especiesrep, resample_bootes)
                                                                if(cccaVal=="CNCA")
                                                                {
                                                                        bootresult<-mapply(ex.cnca, especresample, fvambientales=rep(list(fvambientales),niter), nejes=nejes, Numbercuant=Numbercuant, Numbermixto=Numbermixto, tChoice=tChoice)        
                                                                }else{
                                                                        bootresult<-mapply(ex.cca, especresample, fvambientales=rep(list(fvambientales),niter), nejes=nejes, Numbercuant=Numbercuant, Numbermixto=Numbermixto, tChoice=tChoice)
                                                                }
                                                                eigenalm<-bootresult[19,]
                                                                calidadpsalm<-bootresult[10,]
                                                                calidadosalm<-bootresult[11,]
                                                                Calphakalm<-bootresult[9,]
                                                                qalphakpsalm<-bootresult[14,]
                                                                qalphakosalm<-bootresult[13,]
                                                                Calphaialm<-bootresult[8,]
                                                                qalphaipsalm<-bootresult[12,]
                                                                dalphaiosalm<-bootresult[15,]
                                                                
                                                                if(cccaVal=="CNCA")
                                                                {
                                                                        maresult<-mapply(ex.cnca, jakespecies, fvambientales=rep(list(fvambientales),sum(fespecies)), nejes=nejes, Numbercuant=Numbercuant, Numbermixto=Numbermixto, tChoice=tChoice)
                                                                }else{
                                                                        maresult<-mapply(ex.cca, jakespecies, fvambientales=rep(list(fvambientales),sum(fespecies)), nejes=nejes, Numbercuant=Numbercuant, Numbermixto=Numbermixto, tChoice=tChoice)
                                                                }
                                                                descomjackr<-maresult[19,]
                                                                calidadpsjackr<-maresult[10,]
                                                                calidadosjackr<-maresult[11,]
                                                                Calphakjackr<-maresult[9,]
                                                                qalphakpsjackr<-maresult[14,]
                                                                qalphakosjackr<-maresult[13,]
                                                                Calphaijackr<-maresult[8,]
                                                                qalphaipsjackr<-maresult[12,]
                                                                dalphaiosjackr<-maresult[15,]
                                                                
                                                                
                                                                
                                                                
                                                                ####crear vectores con coordenadas de variables
                                                                coorvarrot<-bootresult[3,]
                                                                coorvarrot<-array(c(unlist(covambien),unlist(coorvarrot)), dim=c(dim(fvambientales)[2], nejes, as.numeric(niter)+1))
                                                                
                                                                out.var<-procGPA(coorvarrot, reflect=TRUE, distances=FALSE, pcaoutput=FALSE)
                                                                plot(out.var$rotated[,dim1,], out.var$rotated[,dim2,], type="n", main=paste("Bootstrap Coordinates (", Namevar, ")"), xlab=paste("Dimension", dim1), ylab=paste("Dimension", dim2), asp=1/1)
                                                                
                                                                for(i in 1:dim(covariablesnam)[1])
                                                                {
                                                                        points(out.var$rotated[i,dim1,],out.var$rotated[i,dim2,], col=colorescoor[i], pch=20)
                                                                        
                                                                }
                                                                
                                                                abline(v=0)
                                                                abline(h=0)
                                                                
                                                                
                                                                for(i in 1:dim(covariablesnam)[1])
                                                                {
                                                                        hpts <- chull(t(out.var$rotated[i,c(dim1,dim2),]))
                                                                        hpts <- c(hpts, hpts[1])
                                                                        lines(t(out.var$rotated[i,c(dim1,dim2),hpts]), col=colorescoor[i])
                                                                }
                                                                text(out.var$rotated[,dim1,1], out.var$rotated[,dim2,1], labels=textvariables)
                                                                
                                                                
                                                                
                                                                ####crear vectores con coordenadas de especies
                                                                coorsperot<-bootresult[2,]
                                                                coorsperot<-array(c(unlist(coespecies), unlist(coorsperot)), dim=c(dim(fespecies)[2], nejes, as.numeric(niter)+1))
                                                                out.spe<-procGPA(coorsperot, reflect=TRUE, distances=FALSE, pcaoutput=FALSE)
                                                                
                                                                plot(out.spe$rotated[,dim1,], out.spe$rotated[,dim2,], type="n", main=paste("Bootstrap Coordinates (", Namespe, ")"), xlab=paste("Dimension", dim1), ylab=paste("Dimension", dim2), asp=1/1)
                                                                
                                                                for(i in 1:dim(coindividuosnam)[1])
                                                                {
                                                                        points(out.spe$rotated[i,dim1,],out.spe$rotated[i,dim2,], col=colorescoor[i], pch=20)
                                                                        
                                                                }
                                                                
                                                                abline(v=0)
                                                                abline(h=0)
                                                                
                                                                
                                                                for(i in 1:dim(coindividuosnam)[1])
                                                                {
                                                                        hpts <- chull(t(out.spe$rotated[i,c(dim1,dim2),]))
                                                                        hpts <- c(hpts, hpts[1])
                                                                        lines(t(out.spe$rotated[i,c(dim1,dim2),hpts]), col=colorescoor[i])
                                                                }
                                                                text(out.spe$rotated[,dim1,1], out.spe$rotated[,dim2,1], labels=textespecies)
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                ####crear vectores con coordenadas de lugares
                                                                coorlugrot<-bootresult[1,]
                                                                coorlugrot<-array(c(unlist(colugares),unlist(coorlugrot)), dim=c(dim(fvambientales)[1], nejes, as.numeric(niter)+1))
                                                                
                                                                out.sit<-procGPA(coorlugrot, reflect=TRUE, distances=FALSE, pcaoutput=FALSE)
                                                                plot(out.sit$rotated[,dim1,], out.sit$rotated[,dim2,], type="n", main=paste("Bootstrap Coordinates (", Namesit, ")"), xlab=paste("Dimension", dim1), ylab=paste("Dimension", dim2), asp=1/1)
                                                                
                                                                for(i in 1:dim(colugaresnam)[1])
                                                                {
                                                                        points(out.sit$rotated[i,dim1,],out.sit$rotated[i,dim2,], col=colorescoor[i], pch=20)
                                                                        
                                                                }
                                                                
                                                                abline(v=0)
                                                                abline(h=0)
                                                                
                                                                
                                                                for(i in 1:dim(colugaresnam)[1])
                                                                {
                                                                        hpts <- chull(t(out.sit$rotated[i,c(dim1,dim2),]))
                                                                        hpts <- c(hpts, hpts[1])
                                                                        lines(t(out.sit$rotated[i,c(dim1,dim2),hpts]), col=colorescoor[i])
                                                                }
                                                                text(out.sit$rotated[,dim1,1], out.sit$rotated[,dim2,1], labels=textlugares)
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                #########################################################################
                                                                ### Guardar resultados
                                                                #########################################################################
                                                                
                                                                cat("File saved in:    ",file="Resultsbootstrap.txt")
                                                                cat(getwd(),file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")	
                                                                cat("\n",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                cat("\n",file="temp.txt")					
                                                                file.append("Resultsbootstrap.txt","temp.txt")		
                                                                
                                                                alphaic <<- tclvalue(Nalpha)
                                                                alphaic <<- as.numeric(alphaic)
                                                                liminf <- (1 - alphaic*0.01) / 2
                                                                limsup <- 1 - (1 - alphaic*0.01) / 2
                                                                
                                                                
                                                                titulo <-c("Obs. Value","Mean","SE","Bias","IC t-boot inf","IC t-boot sup","IC perc inf","IC perc sup","IC BCa inf","IC BCa sup")
                                                                
                                                                ### Eigenvalues
                                                                
                                                                if (ceigenVal=="1")
                                                                {
                                                                        
                                                                        calc.ceigen <-c()
                                                                        ceigen.mean <-c()
                                                                        se.ceigen <-c()
                                                                        sesgo.ceigen <-c()
                                                                        ic.t.ceigeninf <-c()
                                                                        ic.t.ceigensup <-c()
                                                                        ic.p.ceigeninf <-c()
                                                                        ic.p.ceigensup <-c()
                                                                        ic.bca.ceigeninf <-c()
                                                                        ic.bca.ceigensup <-c()
                                                                        
                                                                        
                                                                        for (i in 1:length(eigenalm[[1]]))
                                                                        { 
                                                                                calc.ceigen <-cal.ic(sapply(eigenalm, function(x) x[i]), liminf, limsup, descom$d[i], sapply(descomjackr, function(x) x[i]), niter)
                                                                                ceigen.mean <- c(ceigen.mean, calc.ceigen[1])
                                                                                se.ceigen <- c(se.ceigen,calc.ceigen[2])
                                                                                sesgo.ceigen <- c(sesgo.ceigen,calc.ceigen[3])
                                                                                ic.t.ceigeninf <- c(ic.t.ceigeninf,calc.ceigen[4])
                                                                                ic.t.ceigensup <- c(ic.t.ceigensup,calc.ceigen[5])
                                                                                ic.p.ceigeninf <- c(ic.p.ceigeninf,calc.ceigen[6])
                                                                                ic.p.ceigensup <- c(ic.p.ceigensup,calc.ceigen[7])
                                                                                ic.bca.ceigeninf <- c(ic.bca.ceigeninf,calc.ceigen[8])
                                                                                ic.bca.ceigensup <- c(ic.bca.ceigensup,calc.ceigen[9])
                                                                                
                                                                                pdf(paste("Histogram of eigenvalue", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(eigenalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                                                
                                                                                if(cpdfVal=="Color pdf")
                                                                                {
                                                                                        abline(v=ceigen.mean[i], lwd=2, col="blue")
                                                                                        abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=ceigen.mean[i], lwd=2)
                                                                                        abline(v=descom$d[i], lty =2, lwd=2)
                                                                                }        
                                                                                qqnorm(sapply(eigenalm, function(x) x[i]))
                                                                                dev.off()
                                                                                
                                                                                
                                                                                postscript(paste("Histogram of eigenvalue", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(eigenalm, function(x) x[i]), main="Histogram", xlab=paste("Eigenvalue", i))
                                                                                
                                                                                if(cepsVal=="Color eps")
                                                                                {
                                                                                        abline(v=ceigen.mean[i], lwd=2, col="blue")
                                                                                        abline(v=descom$d[i], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=ceigen.mean[i], lwd=2)
                                                                                        abline(v=descom$d[i], lty =2, lwd=2)
                                                                                }        
                                                                                qqnorm(sapply(eigenalm, function(x) x[i]))
                                                                                dev.off()                                                                               
                                                                        }#end for (i in 1:length(eigenalm))
                                                                        
                                                                        
                                                                        calc.ceigen <-array(cbind(descom$d, ceigen.mean, se.ceigen, sesgo.ceigen, ic.t.ceigeninf, ic.t.ceigensup, ic.p.ceigeninf, ic.p.ceigensup, ic.bca.ceigeninf, ic.bca.ceigensup),
                                                                                            dim=c(length(descom$d),10))
                                                                        calc.ceigen <- as.data.frame(calc.ceigen)
                                                                        colnames(calc.ceigen) <- titulo
                                                                        
                                                                        nombreseig<-c()
                                                                        for (i in 1: length(descom$d))
                                                                        {
                                                                                nombreseig <-c(nombreseig, paste("Eigenvalue",i, sep=""))
                                                                        }#end for (i in 1: length(descom$d))
                                                                        
                                                                        rownames(calc.ceigen) <- nombreseig
                                                                        
                                                                        cat("\n",file="temp.txt")					
                                                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                                                        cat("Eigenvalues: \n",file="temp.txt")					
                                                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                                                        write.table(round(calc.ceigen, digits=3), file="temp.txt", sep="\t", dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end if (ceigenVal=="1")
                                                                
                                                                ### calidadps
                                                                
                                                                if (ccalidadpsVal=="1")
                                                                {
                                                                        
                                                                        calc.ccalps <-c()
                                                                        ccalps.mean <-c()
                                                                        se.ccalps <-c()
                                                                        sesgo.ccalps <-c()
                                                                        ic.t.ccalpsinf <-c()
                                                                        ic.t.ccalpssup <-c()
                                                                        ic.p.ccalpsinf <-c()
                                                                        ic.p.ccalpssup <-c()
                                                                        ic.bca.ccalpsinf <-c()
                                                                        ic.bca.ccalpssup <-c()
                                                                        
                                                                        for (i in 1:length(calidadpsalm[[1]]))
                                                                        { 
                                                                                
                                                                                calc.ccalps <-cal.ic(sapply(calidadpsalm, function(x) x[i]), liminf, limsup, calidadps[i], sapply(calidadpsjackr, function(x) x[i]), niter)
                                                                                ccalps.mean <- c(ccalps.mean, calc.ccalps[1])
                                                                                se.ccalps <- c(se.ccalps,calc.ccalps[2])
                                                                                sesgo.ccalps <- c(sesgo.ccalps,calc.ccalps[3])
                                                                                ic.t.ccalpsinf <- c(ic.t.ccalpsinf,calc.ccalps[4])
                                                                                ic.t.ccalpssup <- c(ic.t.ccalpssup,calc.ccalps[5])
                                                                                ic.p.ccalpsinf <- c(ic.p.ccalpsinf,calc.ccalps[6])
                                                                                ic.p.ccalpssup <- c(ic.p.ccalpssup,calc.ccalps[7])
                                                                                ic.bca.ccalpsinf <- c(ic.bca.ccalpsinf,calc.ccalps[8])
                                                                                ic.bca.ccalpssup <- c(ic.bca.ccalpssup,calc.ccalps[9])
                                                                                
                                                                                pdf(paste("Histogram of inertia by axis (PS)", i, ".pdf", sep = ""), height = 500, width = 500, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(calidadpsalm, function(x) x[i]),main="Histogram", xlab=paste("Axis", i))
                                                                                
                                                                                if(cpdfVal=="Color pdf")
                                                                                {
                                                                                        abline(v=ccalps.mean[i], lwd=2, col="blue")
                                                                                        abline(v=calidadps[i], lty =2, lwd=2, col="red")
                                                                                }else{
                                                                                        abline(v=ccalps.mean[i], lwd=2)
                                                                                        abline(v=calidadps[i], lty =2, lwd=2)     
                                                                                }
                                                                                
                                                                                qqnorm(sapply(calidadpsalm, function(x) x[i]))
                                                                                dev.off()
                                                                                
                                                                                
                                                                                postscript(paste("Histogram of inertia by axis (PS)", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(calidadpsalm, function(x) x[i]),main="Histogram", xlab=paste("Axis", i))
                                                                                
                                                                                if(cepsVal=="Color eps")
                                                                                {
                                                                                        abline(v=ccalps.mean[i], lwd=2, col="blue")
                                                                                        abline(v=calidadps[i], lty =2, lwd=2, col="red")
                                                                                        
                                                                                }else{
                                                                                        abline(v=ccalps.mean[i], lwd=2)
                                                                                        abline(v=calidadps[i], lty =2, lwd=2)        
                                                                                }
                                                                                
                                                                                qqnorm(sapply(calidadpsalm, function(x) x[i]))
                                                                                dev.off()
                                                                        }#end for (i in 1:length(calidadpsalm))
                                                                        
                                                                        
                                                                        calc.ccalps <-array(cbind(t(calidadps), ccalps.mean, se.ccalps, sesgo.ccalps, ic.t.ccalpsinf, ic.t.ccalpssup, ic.p.ccalpsinf, ic.p.ccalpssup, ic.bca.ccalpsinf, ic.bca.ccalpssup),
                                                                                            dim=c(length(calidadps),10))
                                                                        calc.ccalps <- as.data.frame(calc.ccalps)
                                                                        colnames(calc.ccalps) <- titulo
                                                                        
                                                                        nombrescalps<-c()
                                                                        for (i in 1: length(calidadps))
                                                                        {
                                                                                nombrescalps <-c(nombrescalps, paste("Axis",i, sep=""))
                                                                        }#end for (i in 1: length(calidadps))
                                                                        
                                                                        rownames(calc.ccalps) <- nombrescalps
                                                                        
                                                                        cat("\n",file="temp.txt")                  		
                                                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                                                        cat("Proportion of inertia explained by each axis (PS): \n",file="temp.txt")					
                                                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                                                        write.table(round(calc.ccalps, digits=3), file="temp.txt", sep="\t", dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end if (ccalpsVal=="1")
                                                                
                                                                
                                                                
                                                                ### calidados
                                                                
                                                                if (ccalidadosVal=="1")
                                                                {
                                                                        
                                                                        calc.ccalos <-c()
                                                                        ccalos.mean <-c()
                                                                        se.ccalos <-c()
                                                                        sesgo.ccalos <-c()
                                                                        ic.t.ccalosinf <-c()
                                                                        ic.t.ccalossup <-c()
                                                                        ic.p.ccalosinf <-c()
                                                                        ic.p.ccalossup <-c()
                                                                        ic.bca.ccalosinf <-c()
                                                                        ic.bca.ccalossup <-c()
                                                                        
                                                                        for (i in 1:length(calidadosalm[[1]]))
                                                                        { 
                                                                                calc.ccalos <-cal.ic(sapply(calidadosalm, function(x) x[i]), liminf, limsup, calidados[i], sapply(calidadosjackr, function(x) x[i]), niter)
                                                                                ccalos.mean <- c(ccalos.mean, calc.ccalos[1])
                                                                                se.ccalos <- c(se.ccalos,calc.ccalos[2])
                                                                                sesgo.ccalos <- c(sesgo.ccalos,calc.ccalos[3])
                                                                                ic.t.ccalosinf <- c(ic.t.ccalosinf,calc.ccalos[4])
                                                                                ic.t.ccalossup <- c(ic.t.ccalossup,calc.ccalos[5])
                                                                                ic.p.ccalosinf <- c(ic.p.ccalosinf,calc.ccalos[6])
                                                                                ic.p.ccalossup <- c(ic.p.ccalossup,calc.ccalos[7])
                                                                                ic.bca.ccalosinf <- c(ic.bca.ccalosinf,calc.ccalos[8])
                                                                                ic.bca.ccalossup <- c(ic.bca.ccalossup,calc.ccalos[9])
                                                                                
                                                                                pdf(paste("Histogram of inertia by axis (OS)", i, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(calidadosalm, function(x) x[i]),main="Histogram", xlab=paste("Axis", i))
                                                                                
                                                                                if(cpdfVal=="Color pdf")
                                                                                {
                                                                                        abline(v=ccalos.mean[i], lwd=2, col="blue")
                                                                                        abline(v=calidados[i], lty =2, lwd=2, col="red")        
                                                                                }else{
                                                                                        abline(v=ccalos.mean[i], lwd=2)
                                                                                        abline(v=calidados[i], lty =2, lwd=2)
                                                                                }
                                                                                qqnorm(sapply(calidadosalm, function(x) x[i]))
                                                                                dev.off()
                                                                                
                                                                                
                                                                                postscript(paste("Histogram of inertia by axis (OS)", i, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                par(mfrow=c(1,2))
                                                                                hist(sapply(calidadosalm, function(x) x[i]),main="Histogram", xlab=paste("Axis", i))
                                                                                
                                                                                if(cepsVal=="Color eps")
                                                                                {
                                                                                        abline(v=ccalos.mean[i], lwd=2, col="blue")
                                                                                        abline(v=calidados[i], lty =2, lwd=2, col="red")        
                                                                                }else{
                                                                                        abline(v=ccalos.mean[i], lwd=2)
                                                                                        abline(v=calidados[i], lty =2, lwd=2)
                                                                                }
                                                                                qqnorm(sapply(calidadosalm, function(x) x[i]))
                                                                                dev.off()
                                                                        }#end for (i in 1:length(calidadosalm))
                                                                        
                                                                        
                                                                        calc.ccalos <-array(cbind(t(calidados), ccalos.mean, se.ccalos, sesgo.ccalos, ic.t.ccalosinf, ic.t.ccalossup, ic.p.ccalosinf, ic.p.ccalossup, ic.bca.ccalosinf, ic.bca.ccalossup),
                                                                                            dim=c(length(calidados),10))
                                                                        calc.ccalos <- as.data.frame(calc.ccalos)
                                                                        colnames(calc.ccalos) <- titulo
                                                                        
                                                                        nombrescalos<-c()
                                                                        for (i in 1: length(calidados))
                                                                        {
                                                                                nombrescalos <-c(nombrescalos, paste("Axis",i, sep=""))
                                                                        }#end for (i in 1: length(calidados))
                                                                        
                                                                        rownames(calc.ccalos) <- nombrescalos
                                                                        
                                                                        cat("\n",file="temp.txt")    			
                                                                        file.append("Resultsbootstrap.txt","temp.txt")	
                                                                        cat("Proportion of inertia explained by each axis (OS): \n",file="temp.txt")					
                                                                        file.append("Resultsbootstrap.txt","temp.txt")					
                                                                        write.table(round(calc.ccalos, digits=2), file="temp.txt", sep="\t", dec=",")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                }#end if (ccalosVal=="1")
                                                                
                                                                
                                                                ### Calphak
                                                                
                                                                if (cCalphakVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions to", Namespe, ":",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(Calphak)[1])
                                                                        {
                                                                                calc.ccalk <-c()
                                                                                ccalk.mean <-c()
                                                                                se.ccalk <-c()
                                                                                sesgo.ccalk <-c()
                                                                                ic.t.ccalkinf <-c()
                                                                                ic.t.ccalksup <-c()
                                                                                ic.p.ccalkinf <-c()
                                                                                ic.p.ccalksup <-c()
                                                                                ic.bca.ccalkinf <-c()
                                                                                ic.bca.ccalksup <-c()
                                                                                
                                                                                for (j in  1:dim(Calphak)[2])
                                                                                {
                                                                                        calphakaux <- sapply(Calphakalm, function(x) x[i,j])
                                                                                        calc.ccalk <-cal.ic(calphakaux, liminf, limsup, Calphak[i,j], sapply(Calphakjackr, function(x) x[i,j]), niter)
                                                                                        ccalk.mean <- c(ccalk.mean, calc.ccalk[1])
                                                                                        se.ccalk <- c(se.ccalk, calc.ccalk[2])
                                                                                        sesgo.ccalk <- c(sesgo.ccalk,calc.ccalk[3])
                                                                                        ic.t.ccalkinf <- c(ic.t.ccalkinf,calc.ccalk[4])
                                                                                        ic.t.ccalksup <- c(ic.t.ccalksup,calc.ccalk[5])
                                                                                        ic.p.ccalkinf <- c(ic.p.ccalkinf,calc.ccalk[6])
                                                                                        ic.p.ccalksup <- c(ic.p.ccalksup,calc.ccalk[7])
                                                                                        ic.bca.ccalkinf <- c(ic.bca.ccalkinf,calc.ccalk[8])
                                                                                        ic.bca.ccalksup <- c(ic.bca.ccalksup,calc.ccalk[9])
                                                                                        
                                                                                        pdf(paste("Histogram of contributions to ", textespecies[i]," of axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(calphakaux, main="Histogram", xlab=paste("Contribution to ", textespecies[i],"\n of axis ", j))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=ccalk.mean[j], lwd=2, col="blue")
                                                                                                abline(v=Calphak[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=ccalk.mean[j], lwd=2)
                                                                                                abline(v=Calphak[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(calphakaux)
                                                                                        dev.off()
                                                                                        
                                                                                        
                                                                                        postscript(paste("Histogram of contributions to ", textespecies[i]," of axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(calphakaux, main="Histogram", xlab=paste("Contribution to ", textespecies[i],"\n of axis ", j))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=ccalk.mean[j], lwd=2, col="blue")
                                                                                                abline(v=Calphak[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=ccalk.mean[j], lwd=2)
                                                                                                abline(v=Calphak[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(calphakaux)
                                                                                        dev.off()
                                                                                }#end for (j in  1:dimCalphak)[2])
                                                                                
                                                                                calc.ccalk <-array(cbind(unlist(Calphak[i,]), ccalk.mean, se.ccalk, sesgo.ccalk, ic.t.ccalkinf, ic.t.ccalksup, ic.p.ccalkinf, ic.p.ccalksup, ic.bca.ccalkinf, ic.bca.ccalksup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.ccalk <- as.data.frame(calc.ccalk)
                                                                                colnames(calc.ccalk) <- titulo
                                                                                rownames(calc.ccalk) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textespecies[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.ccalk, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(Calphak)[1])
                                                                }#end if (cCalphakVal=="1")
                                                                
                                                                
                                                                ### qalphakps
                                                                
                                                                if (cqalphakpsVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions of ", Namespe, " (PS):",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(qalphakps)[1])
                                                                        {
                                                                                calc.cqkps <-c()
                                                                                cqkps.mean <-c()
                                                                                se.cqkps <-c()
                                                                                sesgo.cqkps <-c()
                                                                                ic.t.cqkpsinf <-c()
                                                                                ic.t.cqkpssup <-c()
                                                                                ic.p.cqkpsinf <-c()
                                                                                ic.p.cqkpssup <-c()
                                                                                ic.bca.cqkpsinf <-c()
                                                                                ic.bca.cqkpssup <-c()
                                                                                
                                                                                for (j in  1:dim(qalphakps)[2])
                                                                                {
                                                                                        qalphakpsaux <- sapply(qalphakpsalm, function(x) x[i,j])
                                                                                        calc.cqkps <-cal.ic(qalphakpsaux, liminf, limsup, qalphakps[i,j], sapply(qalphakpsjackr, function(x) x[i,j]), niter)
                                                                                        cqkps.mean <- c(cqkps.mean, calc.cqkps[1])
                                                                                        se.cqkps <- c(se.cqkps, calc.cqkps[2])
                                                                                        sesgo.cqkps <- c(sesgo.cqkps,calc.cqkps[3])
                                                                                        ic.t.cqkpsinf <- c(ic.t.cqkpsinf,calc.cqkps[4])
                                                                                        ic.t.cqkpssup <- c(ic.t.cqkpssup,calc.cqkps[5])
                                                                                        ic.p.cqkpsinf <- c(ic.p.cqkpsinf,calc.cqkps[6])
                                                                                        ic.p.cqkpssup <- c(ic.p.cqkpssup,calc.cqkps[7])
                                                                                        ic.bca.cqkpsinf <- c(ic.bca.cqkpsinf,calc.cqkps[8])
                                                                                        ic.bca.cqkpssup <- c(ic.bca.cqkpssup,calc.cqkps[9])
                                                                                        
                                                                                        pdf(paste("Histogram of contribution of ", textespecies[i]," (PS)", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphakpsaux,main="Histogram", xlab=paste("Contribution of ", textespecies[i]," (PS)"))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=cqkps.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphakps[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqkps.mean[j], lwd=2)
                                                                                                abline(v=qalphakps[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphakpsaux)
                                                                                        dev.off()
                                                                                        
                                                                                        
                                                                                        postscript(paste("Histogram of contribution of ", textespecies[i]," (PS)", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphakpsaux,main="Histogram", xlab=paste("Contribution of ", textespecies[i]," (PS)"))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=cqkps.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphakps[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqkps.mean[j], lwd=2)
                                                                                                abline(v=qalphakps[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphakpsaux)
                                                                                        dev.off()
                                                                                }#end for (j in  1:dim (qalphakps)[2])
                                                                                
                                                                                calc.cqkps <-array(cbind(unlist(qalphakps[i,]), cqkps.mean, se.cqkps, sesgo.cqkps, ic.t.cqkpsinf, ic.t.cqkpssup, ic.p.cqkpsinf, ic.p.cqkpssup, ic.bca.cqkpsinf, ic.bca.cqkpssup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.cqkps <- as.data.frame(calc.cqkps)
                                                                                colnames(calc.cqkps) <- titulo
                                                                                rownames(calc.cqkps) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textespecies[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.cqkps, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(qalphakps)[1])
                                                                }#end if (cqalphakpsVal=="1")
                                                                
                                                                
                                                                ### qalphakos
                                                                
                                                                if (cqalphakosVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions of ", Namespe, " (OS):",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(qalphakos)[1])
                                                                        {
                                                                                calc.cqkos <-c()
                                                                                cqkos.mean <-c()
                                                                                se.cqkos <-c()
                                                                                sesgo.cqkos <-c()
                                                                                ic.t.cqkosinf <-c()
                                                                                ic.t.cqkossup <-c()
                                                                                ic.p.cqkosinf <-c()
                                                                                ic.p.cqkossup <-c()
                                                                                ic.bca.cqkosinf <-c()
                                                                                ic.bca.cqkossup <-c()
                                                                                
                                                                                for (j in  1:dim(qalphakos)[2])
                                                                                {
                                                                                        qalphakosaux <- sapply(qalphakosalm, function(x) x[i,j])
                                                                                        calc.cqkos <-cal.ic(qalphakosaux, liminf, limsup, qalphakos[i,j], sapply(qalphakosjackr, function(x) x[i,j]), niter)
                                                                                        cqkos.mean <- c(cqkos.mean, calc.cqkos[1])
                                                                                        se.cqkos <- c(se.cqkos, calc.cqkos[2])
                                                                                        sesgo.cqkos <- c(sesgo.cqkos,calc.cqkos[3])
                                                                                        ic.t.cqkosinf <- c(ic.t.cqkosinf,calc.cqkos[4])
                                                                                        ic.t.cqkossup <- c(ic.t.cqkossup,calc.cqkos[5])
                                                                                        ic.p.cqkosinf <- c(ic.p.cqkosinf,calc.cqkos[6])
                                                                                        ic.p.cqkossup <- c(ic.p.cqkossup,calc.cqkos[7])
                                                                                        ic.bca.cqkosinf <- c(ic.bca.cqkosinf,calc.cqkos[8])
                                                                                        ic.bca.cqkossup <- c(ic.bca.cqkossup,calc.cqkos[9])
                                                                                        
                                                                                        pdf(paste("Histogram of contribution of ", textespecies[i]," (OS)", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphakosaux,main="Histogram", xlab=paste("Contribution of ", textespecies[i]," (OS)"))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=cqkos.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphakos[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqkos.mean[j], lwd=2)
                                                                                                abline(v=qalphakos[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphakosaux)
                                                                                        dev.off()
                                                                                        
                                                                                        
                                                                                        postscript(paste("Contribution of ", textespecies[i]," (OS)", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphakosaux,main="Histogram", xlab=paste("Contribution of ", textespecies[i]," (OS)"))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=cqkos.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphakos[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqkos.mean[j], lwd=2)
                                                                                                abline(v=qalphakos[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphakosaux)
                                                                                        dev.off()
                                                                                }#end for (j in  1:dim (qalphakos)[2])
                                                                                
                                                                                calc.cqkos <-array(cbind(unlist(qalphakos[i,]), cqkos.mean, se.cqkos, sesgo.cqkos, ic.t.cqkosinf, ic.t.cqkossup, ic.p.cqkosinf, ic.p.cqkossup, ic.bca.cqkosinf, ic.bca.cqkossup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.cqkos <- as.data.frame(calc.cqkos)
                                                                                colnames(calc.cqkos) <- titulo
                                                                                rownames(calc.cqkos) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textespecies[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.cqkos, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(qalphakos)[1])
                                                                }#end if (cqalphakosVal=="1")
                                                                
                                                                
                                                                
                                                                ### Calphai
                                                                
                                                                if (cCalphaiVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions to", Namesit, ":",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(Calphai)[1])
                                                                        {
                                                                                calc.ccali <-c()
                                                                                ccali.mean <-c()
                                                                                se.ccali <-c()
                                                                                sesgo.ccali <-c()
                                                                                ic.t.ccaliinf <-c()
                                                                                ic.t.ccalisup <-c()
                                                                                ic.p.ccaliinf <-c()
                                                                                ic.p.ccalisup <-c()
                                                                                ic.bca.ccaliinf <-c()
                                                                                ic.bca.ccalisup <-c()
                                                                                
                                                                                for (j in  1:dim(Calphai)[2])
                                                                                {
                                                                                        calphaiaux <- sapply(Calphaialm, function(x) x[i,j])
                                                                                        calc.ccali <-cal.ic(calphaiaux, liminf, limsup, Calphai[i,j], sapply(Calphaijackr, function(x) x[i,j]), niter)
                                                                                        ccali.mean <- c(ccali.mean, calc.ccali[1])
                                                                                        se.ccali <- c(se.ccali, calc.ccali[2])
                                                                                        sesgo.ccali <- c(sesgo.ccali,calc.ccali[3])
                                                                                        ic.t.ccaliinf <- c(ic.t.ccaliinf,calc.ccali[4])
                                                                                        ic.t.ccalisup <- c(ic.t.ccalisup,calc.ccali[5])
                                                                                        ic.p.ccaliinf <- c(ic.p.ccaliinf,calc.ccali[6])
                                                                                        ic.p.ccalisup <- c(ic.p.ccalisup,calc.ccali[7])
                                                                                        ic.bca.ccaliinf <- c(ic.bca.ccaliinf,calc.ccali[6])
                                                                                        ic.bca.ccalisup <- c(ic.bca.ccalisup,calc.ccali[7])
                                                                                        
                                                                                        pdf(paste("Histogram of contributions to ", textlugares[i]," of axis ", j, ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(calphaiaux,main="Histogram", xlab=paste("Contribution to ", textlugares[i],"\n of axis ", j))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=ccali.mean[j], lwd=2, col="blue")
                                                                                                abline(v=Calphai[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=ccali.mean[j], lwd=2)
                                                                                                abline(v=Calphai[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(calphaiaux)
                                                                                        dev.off()
                                                                                        
                                                                                        
                                                                                        postscript(paste("Histogram of contributions to ", textlugares[i]," of axis ", j, ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(calphaiaux,main="Histogram", xlab=paste("Contribution to ", textlugares[i],"\n of axis ", j))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=ccali.mean[j], lwd=2, col="blue")
                                                                                                abline(v=Calphai[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=ccali.mean[j], lwd=2)
                                                                                                abline(v=Calphai[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(calphaiaux)
                                                                                        dev.off()
                                                                                }#end for (j in  1:dimCalphai)[2])
                                                                                
                                                                                calc.ccali <-array(cbind(unlist(Calphai[i,]), ccali.mean, se.ccali, sesgo.ccali, ic.t.ccaliinf, ic.t.ccalisup, ic.p.ccaliinf, ic.p.ccalisup, ic.bca.ccaliinf, ic.bca.ccalisup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.ccali <- as.data.frame(calc.ccali)
                                                                                colnames(calc.ccali) <- titulo
                                                                                rownames(calc.ccali) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textlugares[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.ccali, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(Calphai)[1])
                                                                }#end if (cCalphaiVal=="1")
                                                                
                                                                
                                                                ### qalphaips
                                                                
                                                                if (cqalphaipsVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions of ", Namesit, " (PS):",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(qalphaips)[1])
                                                                        {
                                                                                calc.cqips <-c()
                                                                                cqips.mean <-c()
                                                                                se.cqips <-c()
                                                                                sesgo.cqips <-c()
                                                                                ic.t.cqipsinf <-c()
                                                                                ic.t.cqipssup <-c()
                                                                                ic.p.cqipsinf <-c()
                                                                                ic.p.cqipssup <-c()
                                                                                ic.bca.cqipsinf <-c()
                                                                                ic.bca.cqipssup <-c()
                                                                                
                                                                                for (j in  1:dim(qalphaips)[2])
                                                                                {
                                                                                        qalphaipsaux <- sapply(qalphaipsalm, function(x) x[i,j])
                                                                                        calc.cqips <-cal.ic(qalphaipsaux, liminf, limsup, qalphaips[i,j], sapply(qalphaipsjackr, function(x) x[i,j]), niter)
                                                                                        cqips.mean <- c(cqips.mean, calc.cqips[1])
                                                                                        se.cqips <- c(se.cqips, calc.cqips[2])
                                                                                        sesgo.cqips <- c(sesgo.cqips,calc.cqips[3])
                                                                                        ic.t.cqipsinf <- c(ic.t.cqipsinf,calc.cqips[4])
                                                                                        ic.t.cqipssup <- c(ic.t.cqipssup,calc.cqips[5])
                                                                                        ic.p.cqipsinf <- c(ic.p.cqipsinf,calc.cqips[6])
                                                                                        ic.p.cqipssup <- c(ic.p.cqipssup,calc.cqips[7])
                                                                                        ic.bca.cqipsinf <- c(ic.bca.cqipsinf,calc.cqips[8])
                                                                                        ic.bca.cqipssup <- c(ic.bca.cqipssup,calc.cqips[9])
                                                                                        
                                                                                        pdf(paste("Histogram of contribution of ", textlugares[i]," (PS)", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphaipsaux,main="Histogram", xlab=paste("Contribution of ", textlugares[i]," (PS)"))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=cqips.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphaips[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqips.mean[j], lwd=2)
                                                                                                abline(v=qalphaips[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphaipsaux)
                                                                                        dev.off()
                                                                                        
                                                                                        
                                                                                        postscript(paste("Histogram of contribution of ", textlugares[i]," (PS)", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphaipsaux,main="Histogram", xlab=paste("Contribution of ", textlugares[i]," (PS)"))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=cqips.mean[j], lwd=2, col="blue")
                                                                                                abline(v=qalphaips[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqips.mean[j], lwd=2)
                                                                                                abline(v=qalphaips[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphaipsaux)
                                                                                        dev.off()
                                                                                        
                                                                                }#end for (j in  1:dim (qalphaips)[2])
                                                                                
                                                                                
                                                                                calc.cqips <-array(cbind(unlist(qalphaips[i,]), cqips.mean, se.cqips, sesgo.cqips, ic.t.cqipsinf, ic.t.cqipssup, ic.p.cqipsinf, ic.p.cqipssup, ic.bca.cqipsinf, ic.bca.cqipssup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.cqips <- as.data.frame(calc.cqips)
                                                                                colnames(calc.cqips) <- titulo
                                                                                rownames(calc.cqips) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textlugares[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.cqips, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(qalphaips)[1])
                                                                }#end if (cqalphaipsVal=="1")
                                                                
                                                                
                                                                ### dalphaios
                                                                
                                                                if (cdalphaiosVal=="1")
                                                                {
                                                                        cat("\n",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        cat("Contributions of ", Namesit, " (OS):",file="temp.txt")
                                                                        file.append("Resultsbootstrap.txt","temp.txt")
                                                                        
                                                                        for(i in 1: dim(dalphaios)[1])
                                                                        {
                                                                                calc.cqios <-c()
                                                                                cqios.mean <-c()
                                                                                se.cqios <-c()
                                                                                sesgo.cqios <-c()
                                                                                ic.t.cqiosinf <-c()
                                                                                ic.t.cqiossup <-c()
                                                                                ic.p.cqiosinf <-c()
                                                                                ic.p.cqiossup <-c()
                                                                                ic.bca.cqiosinf <-c()
                                                                                ic.bca.cqiossup <-c()
                                                                                
                                                                                for (j in  1:dim(dalphaios)[2])
                                                                                {
                                                                                        qalphaiosaux <- sapply(dalphaiosalm, function(x) x[i,j])
                                                                                        calc.cqios <-cal.ic(qalphaiosaux, liminf, limsup, dalphaios[i,j], sapply(dalphaiosjackr, function(x) x[i,j]), niter)
                                                                                        cqios.mean <- c(cqios.mean, calc.cqios[1])
                                                                                        se.cqios <- c(se.cqios, calc.cqios[2])
                                                                                        sesgo.cqios <- c(sesgo.cqios,calc.cqios[3])
                                                                                        ic.t.cqiosinf <- c(ic.t.cqiosinf,calc.cqios[4])
                                                                                        ic.t.cqiossup <- c(ic.t.cqiossup,calc.cqios[5])
                                                                                        ic.p.cqiosinf <- c(ic.p.cqiosinf,calc.cqios[6])
                                                                                        ic.p.cqiossup <- c(ic.p.cqiossup,calc.cqios[7])
                                                                                        ic.bca.cqiosinf <- c(ic.bca.cqiosinf,calc.cqios[8])
                                                                                        ic.bca.cqiossup <- c(ic.bca.cqiossup,calc.cqios[9])
                                                                                        
                                                                                        pdf(paste("Histogram of contribution of ", textlugares[i]," (OS)", ".pdf", sep = ""), height = 7, width = 7, useDingbats=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphaiosaux,main="Histogram", xlab=paste("Contribution of ", textlugares[i]," (OS)"))
                                                                                        
                                                                                        if(cpdfVal=="Color pdf")
                                                                                        {
                                                                                                abline(v=cqios.mean[j], lwd=2, col="blue")
                                                                                                abline(v=dalphaios[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqios.mean[j], lwd=2)
                                                                                                abline(v=dalphaios[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphaiosaux)
                                                                                        dev.off()
                                                                                        
                                                                                        postscript(paste("Histogram of contribution of ", textlugares[i]," (OS)", ".eps", sep = ""), height = 600, width = 1200, horizontal=FALSE)
                                                                                        par(mfrow=c(1,2))
                                                                                        hist(qalphaiosaux,main="Histogram", xlab=paste("Contribution of ", textlugares[i]," (OS)"))
                                                                                        
                                                                                        if(cepsVal=="Color eps")
                                                                                        {
                                                                                                abline(v=cqios.mean[j], lwd=2, col="blue")
                                                                                                abline(v=dalphaios[i,j], lty =2, lwd=2, col="red")        
                                                                                        }else{
                                                                                                abline(v=cqios.mean[j], lwd=2)
                                                                                                abline(v=dalphaios[i,j], lty =2, lwd=2)        
                                                                                        }
                                                                                        qqnorm(qalphaiosaux)
                                                                                        dev.off()
                                                                                }#end for (j in  1:dim (dalphaios)[2])
                                                                                
                                                                                calc.cqios <-array(cbind(unlist(dalphaios[i,]), cqios.mean, se.cqios, sesgo.cqios, ic.t.cqiosinf, ic.t.cqiossup, ic.p.cqiosinf, ic.p.cqiossup, ic.bca.cqiosinf, ic.bca.cqiossup),
                                                                                                   dim=c(nejes,10))
                                                                                calc.cqios <- as.data.frame(calc.cqios)
                                                                                colnames(calc.cqios) <- titulo
                                                                                rownames(calc.cqios) <- ejes
                                                                                
                                                                                cat("\n",file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                cat(textlugares[i],"\n", file="temp.txt")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                                write.table(round(calc.cqios, digits=2),file="temp.txt", sep="\t",dec=",")
                                                                                file.append("Resultsbootstrap.txt","temp.txt")
                                                                        }#end for(i in 1: dim(dalphaios)[1])
                                                                }#end if (cdalphaiosVal=="1")
                                                                
                                                                
                                                                
                                                                
                                                                
                                                                file.show("Resultsbootstrap.txt")
                                                                file.remove("temp.txt")
                                                        }
                                                        
                                                        OK.butboot <-tkbutton(framewigr,text="   OK   ",command=OnOKboot, bg= "lightblue", width=20, foreground = "navyblue")
                                                        
                                                        tkpack(OK.butboot, expand="TRUE", side= "left", fill ="both")
                                                        tkpack(framewi21,framewi22, expand = "TRUE",side="left", fill="x")
                                                        tkpack(framewi1,framewi2, expand = "TRUE",side="top", fill="both")
                                                        tkpack(framewi,framewigr, expand = "TRUE",side="top", fill="y")
                                                        
                                                        
                                                        
                                                }#end bootcnca
                                                
                                                gcoia <-function(typeg)
                                                {
                                                        if(typeg=="all")
                                                        {
                                                                winall<-tktoplevel()
                                                                tkwm.title(winall,"All Graphs")
                                                                
                                                                plotall <- function()
                                                                {
                                                                        #nomejes<-c()
                                                                        #for (i in 1: nejes)
                                                                        #{
                                                                        #        nomejes <<-c(nomejes, paste("Axis",i, sep=""))
                                                                        #}#end for (i in 1:nejes)
                                                                        
                                                                        layout(matrix(c(1, 1, 2, 1, 1, 3, 4, 5, 6), 3, 3))
                                                                        par(mar=c(4,4,2,2))
                                                                        ##coinercia
                                                                        datosaux<-rbind(colugaresz, colugaresy)
                                                                        limixa <- round(c(min(datosaux[,dim1],0), max(datosaux[,dim1],0)), digits=2)
                                                                        limiya <- round(c(min(datosaux[,dim2],0), max(datosaux[,dim2],0)), digits=2)
                                                                        
                                                                        plot(datosaux[,c(dim1,dim2)], main=Namesit, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limixa * 1.01, ylim = limiya * 1.01)
                                                                        arrows(colugaresz[,dim1],colugaresz[,dim2],colugaresy[,dim1],colugaresy[,dim2],col=collugares,lty=1, length=0.05)
                                                                        
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        
                                                                        datost <<-colugaresy
                                                                        for (d in 1:ejes)
                                                                        {
                                                                                dimension <- cbind(colugaresy[,d], colugaresz[,d])
                                                                                datost[,d]<<-rowMeans(dimension)
                                                                        }  
                                                                        for (i in 1: length(textlugares))
                                                                        {
                                                                                boxed.labels(datost[i,dim1],datost[i,dim2],labels=textlugares[i], col= collugares[i], cex= as.numeric(cexlugares[i]), border = collugares[i])
                                                                        }
                                                                        
                                                                        
                                                                        ##especies
                                                                        datosaux<-coespecies
                                                                        limixa <- round(c(min(datosaux[,dim1],0), max(datosaux[,dim1],0)), digits=2)
                                                                        limiya <- round(c(min(datosaux[,dim2],0), max(datosaux[,dim2],0)), digits=2)
                                                                        
                                                                        plot(datosaux[,c(dim1,dim2)], main= Namespe, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limixa * 1.01, ylim = limiya * 1.01)
                                                                        arrows(centro[1],centro[2],coespecies[,dim1],coespecies[,dim2],col=colespecies,lty=1, length=0.08)
                                                                        
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        text(datosaux[,dim1],datosaux[,dim2],labels=textespecies, col= colespecies, cex= as.numeric(cexespecies))
                                                                        
                                                                        
                                                                        ##variables
                                                                        datosaux<-covambien
                                                                        limixa <- round(c(min(datosaux[,dim1],0), max(datosaux[,dim1],0)), digits=2)
                                                                        limiya <- round(c(min(datosaux[,dim2],0), max(datosaux[,dim2],0)), digits=2)
                                                                        plot(datosaux[,c(dim1,dim2)], main= Namevar, type="n",xlab=paste(round(inerciatot[dim1]*100, digits=2),"%"),ylab=paste(round(inerciatot[dim2]*100,digits=2),"%"), asp=1/1, xlim = limixa * 1.01, ylim = limiya * 1.01)
                                                                        arrows(centro[1],centro[2],covambien[,dim1],covambien[,dim2],col=colvariables,lty=1, length=0.08)
                                                                        
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        text(datosaux[,dim1],datosaux[,dim2],labels=textvariables, col= colvariables, cex= as.numeric(cexvariables))
                                                                        
                                                                        
                                                                        ##ejes especies
                                                                        plot(0,0, type="n", xlim=c(-1,1), ylim=c(-1,1), asp=1/1, main=paste(Namespe, " axes"), xlab="", ylab="")
                                                                        draw.circle(0,0,1, border="black", lty=1,lwd=1)
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        arrows(centro[1],centro[2],Xaxes[,dim1],Xaxes[,dim2],lty=1, length=0.08)
                                                                        text(Xaxes[,dim1],Xaxes[,dim2],labels=paste("Axis", 1:nejes))
                                                                        
                                                                        
                                                                        
                                                                        ##ejes variables
                                                                        plot(0,0, type="n", xlim=c(-1,1), ylim=c(-1,1), asp=1/1, main=paste(Namevar, " axes"), xlab="", ylab="")
                                                                        draw.circle(0,0,1, border="black", lty=1,lwd=1)
                                                                        if (cbVal=="1"){
                                                                                abline(h=centro[2],v=centro[1],lty="dotted")
                                                                        }#end if (cbVal=="1")
                                                                        arrows(centro[1],centro[2],Yaxes[,dim1],Yaxes[,dim2],lty=1, length=0.08)
                                                                        text(Yaxes[,dim1],Yaxes[,dim2],labels=paste("Axis", 1:nejes))
                                                                        
                                                                        #layout(matrix(c(1),1, 1))
                                                                        
                                                                }
                                                                imgall <- tkrplot(winall,fun=plotall,hscale=as.numeric(hescale),vscale=as.numeric(vescale))
                                                                tkpack(imgall, expand="TRUE", fill="both")  
                                                                tkfocus(winall)
                                                                
                                                        }else{
                                                                
                                                                typecoia<<-typeg
                                                                if(typecoia=="sg") 
                                                                {
                                                                        datos <<-rbind(colugaresy, colugaresz)
                                                                        textos <<-colugaresy
                                                                        for (d in 1:ejes)
                                                                        {
                                                                                dimension <- cbind(colugaresy[,d], colugaresz[,d])
                                                                                textos[,d]<<-rowMeans(dimension)
                                                                        }  
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                        tclvalue(Limix1) <<- limix[1]
                                                                        tclvalue(Limix2) <<- limix[2]
                                                                        tclvalue(Limiy1) <<- limiy[1]
                                                                        tclvalue(Limiy2) <<- limiy[2]
                                                                }
                                                                
                                                                if(typecoia=="eg")
                                                                {
                                                                        datos<<-coespecies
                                                                        textos<<-datos
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                        tclvalue(Limix1) <<- limix[1]
                                                                        tclvalue(Limix2) <<- limix[2]
                                                                        tclvalue(Limiy1) <<- limiy[1]
                                                                        tclvalue(Limiy2) <<- limiy[2]
                                                                }                
                                                                
                                                                if(typecoia=="vg")
                                                                {
                                                                        datos<<-covambien
                                                                        textos<<-datos
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                        tclvalue(Limix1) <<- limix[1]
                                                                        tclvalue(Limix2) <<- limix[2]
                                                                        tclvalue(Limiy1) <<- limiy[1]
                                                                        tclvalue(Limiy2) <<- limiy[2]
                                                                }
                                                                
                                                                if(typecoia=="ef")
                                                                {
                                                                        datos <<-Xaxes
                                                                        textos <<-datos
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                        tclvalue(Limix1) <<- limix[1]
                                                                        tclvalue(Limix2) <<- limix[2]
                                                                        tclvalue(Limiy1) <<- limiy[1]
                                                                        tclvalue(Limiy2) <<- limiy[2]
                                                                }
                                                                
                                                                if(typecoia=="vf")
                                                                {
                                                                        datos <<-Yaxes
                                                                        textos <<-datos
                                                                        limix <<- round(c(min(datos[,dim1],0), max(datos[,dim1],0)), digits=2)
                                                                        limiy <<- round(c(min(datos[,dim2],0), max(datos[,dim2],0)), digits=2)
                                                                        tclvalue(Limix1) <<- limix[1]
                                                                        tclvalue(Limix2) <<- limix[2]
                                                                        tclvalue(Limiy1) <<- limiy[1]
                                                                        tclvalue(Limiy2) <<- limiy[2]
                                                                }
                                                                
                                                                tkrreplot(img)
                                                                
                                                        }
                                                }#end gcoia
                                                
                                                projections <- function(project)
                                                {					
                                                        proj<<-project
                                                        if(proj=="normal")
                                                        {
                                                                tkrreplot(img)
                                                                return()
                                                        }#end if(proj=="normal")
                                                        if(proj=="e" | (proj=="s" & cblugares=="1")) 
                                                        {
                                                                wproj <- tktoplevel()
                                                                tkwm.title(wproj, paste("Select ",Namevar))
                                                                
                                                                scrproj <- tkscrollbar(wproj, repeatinterval=5, command=function(...)tkyview(tlproj,...))
                                                                tlproj<-tklistbox(wproj,height=6,width=42,yscrollcommand=function(...)tkset(scrproj,...),background="white")
                                                                
                                                                for (i in 1:Numbercuant)
                                                                {
                                                                        tkinsert(tlproj,"end",textvariables[i])
                                                                }#end for (i in 1:Numbercuant)
                                                                
                                                                tkselection.set(tlproj,0) #  Indexing starts at zero.
                                                                
                                                                OnOKproj <- function()
                                                                {
                                                                        Choiceproj <<- as.numeric(tkcurselection(tlproj))+1
                                                                        tkdestroy(wproj)
                                                                        tkrreplot(img)
                                                                }#end OnOKproj <- function()
                                                                
                                                                OK.butp <- tkbutton(wproj, text = "   OK   ", command = OnOKproj)
                                                                tkpack(tlproj,scrproj,expand = "TRUE", side="left", fill = "both")
                                                                tkpack.configure(scrproj,side="left")
                                                                tkpack(OK.butp,expand = "TRUE", side="left", fill = "both")
                                                                tkfocus(wproj)
                                                                tkwait.window(wproj)
                                                        }#end if       
                                                }#end projections <-function()
                                                
                                                showaxes <- function()
                                                {
                                                        if(cbVal=="1")
                                                        {
                                                                cbVal<<-"0"
                                                                tkrreplot(img)
                                                        }else{
                                                                cbVal<<-"1"
                                                                tkrreplot(img)
                                                        }#end if(cbVal=="1")        
                                                }#end showaxes
                                                
                                                showsites <- function()
                                                {
                                                        if(cblugares=="1")
                                                        {
                                                                cblugares<<-"0"
                                                                cblablugares<<-"0"
                                                                tkrreplot(img)
                                                        }else{
                                                                cblugares<<-"1"
                                                                tkrreplot(img)
                                                        }#end if(cblugares=="1")        
                                                }#end showsites
                                                
                                                showlabsites <- function()
                                                {
                                                        if(cblablugares=="1")
                                                        {
                                                                cblablugares<<-"0"
                                                                tkrreplot(img)
                                                        }else{
                                                                cblablugares<<-"1"
                                                                cblugares<<-"1"
                                                                indexLabeled <<- c(indexLabeled,
                                                                                   (length(textespecies)+length(textvariables)+1):(length(textespecies)+length(textvariables)+length(textlugares)))
                                                                tkrreplot(img)
                                                        }#end if(cblablugares=="1")        
                                                }#end showlabsites
                                                
                                                changetit <- function()
                                                {
                                                        ctwin<-tktoplevel()
                                                        tkwm.title(ctwin,"Change title")
                                                        OnOKchantit <- function()
                                                        {
                                                                tit_graph <<- tclvalue(tit_gr)
                                                                tkrreplot(img)
                                                                tkdestroy(ctwin)
                                                                
                                                        }
                                                        OK.butchantit<-tkbutton(ctwin,text=" Change ", command=OnOKchantit,  bg= "lightblue", width=20, foreground = "navyblue")
                                                        tkbind(OK.butchantit, "<Return>",OnOKchantit)
                                                        
                                                        tit_gr<-tclVar(tit_graph)
                                                        entry.tit <-tkentry(ctwin, width="50",textvariable=tit_gr, bg="white")
                                                        tkbind(entry.tit, "<Return>",OnOKchantit)
                                                        
                                                        
                                                        tkpack(tklabel(ctwin,text="New title:    "),entry.tit, expand = "TRUE", side="left", fill = "both")
                                                        tkpack(OK.butchantit)
                                                        
                                                        tkfocus(ctwin)
                                                        
                                                }#end changetit
                                                
                                                topMenugr <- tkmenu(wgr)
                                                tkconfigure(wgr, menu = topMenugr)
                                                menuFile <- tkmenu(topMenugr, tearoff = FALSE)
                                                menuSaveAs <- tkmenu(topMenugr, tearoff = FALSE)
                                                menu3d <- tkmenu(topMenugr, tearoff = FALSE)
                                                menuproj <- tkmenu(topMenugr, tearoff = FALSE)          				
                                                menuboot <- tkmenu(topMenugr, tearoff = FALSE)
                                                menuopt <-tkmenu(topMenugr, tearoff = FALSE)
                                                menucoia <-tkmenu(topMenugr, tearoff = FALSE)
                                                
                                                tkadd(menuFile, "command", label = "Copy image", command = function() {tkrreplot(img)})
                                                tkadd(menuFile, "cascade", label = "Save image", menu = menuSaveAs)
                                                tkadd(menuSaveAs, "command", label = "Pdf file", command = function() {SaveFilePDF()})
                                                tkadd(menuSaveAs, "command", label = "Eps file", command = function() {SaveFileeps()})
                                                #tkadd(menuSaveAs, "command", label = "Bmp file", command = function() {SaveFileBmp()})
                                                tkadd(menuSaveAs, "command", label = "Png file", command = function() {SaveFilePng()})
                                                tkadd(menuSaveAs, "command", label = "Jpg/Jpeg file", command = function() {SaveFileJPG()})
                                                tkadd(menuFile, "separator")
                                                tkadd(menuFile, "command", label = "Exit", command = function() {tkdestroy(wgr)})
                                                tkadd(menu3d, "command", label = "3D", command = function() {g3d()})
                                                tkadd(menuboot, "command", label = "Bootstrap", command = function() {bootcnca()})
                                                tkadd(menuopt, "command", label = "Change title", command = function() {changetit()})
                                                tkadd(menuopt, "command", label = "Show/Hide axes", command = function() {showaxes()})
                                                if (cccaVal!="COIA")
                                                {
                                                        tkadd(menuopt, "command", label = paste("Show/Hide ", Namesit), command = function() {showsites()})
                                                        tkadd(menuopt, "command", label = paste("Show/Hide labels for", Namesit), command = function() {showlabsites()})       
                                                }
                                                
                                                tkadd(topMenugr, "cascade", label ="File", menu = menuFile)
                                                tkadd(menuFile, "separator")
                                                tkadd(topMenugr, "cascade", label = "3D", menu = menu3d)
                                                if (cccaVal!="COIA")
                                                { 
                                                        tkadd(menuproj, "command", label = Namespe, command = function() {projections(project="e")})
                                                        tkadd(menuproj, "command", label = Namesit, command = function() {projections(project="s")})        					
                                                        tkadd(menuproj, "command", label = "Back to original graph", command = function() {projections(project="normal")})						
                                                        tkadd(topMenugr, "cascade", label = "Projections", menu = menuproj)
                                                        tkadd(topMenugr, "cascade", label = "Bootstrap", menu = menuboot)
                                                        
                                                }else{
                                                        tkadd(topMenugr, "cascade", label = "COIA", menu = menucoia) 
                                                        tkadd(menucoia, "command", label = "Show all graphs", command = function() {gcoia(typeg="all")})
                                                        tkadd(menucoia, "command", label = paste(Namesit, " graph"), command = function() {gcoia(typeg="sg")})
                                                        tkadd(menucoia, "command", label = paste(Namespe, " graph"), command = function() {gcoia(typeg="eg")})
                                                        tkadd(menucoia, "command", label = paste(Namevar, " graph"), command = function() {gcoia(typeg="vg")})
                                                        tkadd(menucoia, "command", label = paste(Namespe, " axes"), command = function() {gcoia(typeg="ef")})
                                                        tkadd(menucoia, "command", label = paste(Namevar, " axes"), command = function() {gcoia(typeg="vf")})
                                                        
                                                }
                                                tkadd(topMenugr, "cascade", label = "Options", menu = menuopt)
                                                
                                                
                                                
                                                OnLeftClick.up <- function(x,y)
                                                {
                                                        if(cccaVal!="COIA" | typecoia %in% c("all", "sg", "eg", "vg"))
                                                        {
                                                                
                                                                msg <- ("-To change the label press Yes.\n-To remove it press No.\n-If you do not want to do anything press Cancel.")
                                                                mbval<- tkmessageBox(title="Change of label", message=msg,type="yesnocancel",icon="question")
                                                                if (tclvalue(mbval)=="yes"){  
                                                                        indexLabeled <<- c(indexLabeled,indexClosest)
                                                                }#end if (tclvalue(mbval)=="yes")
                                                                
                                                                if(tclvalue(mbval)=="no"){
                                                                        
                                                                        indexLabeledaux<<-c()
                                                                        for (i in (1:length(indexLabeled)))
                                                                        {       
                                                                                if (indexLabeled[i]!=indexClosest)
                                                                                        indexLabeledaux <<- c(indexLabeledaux,indexLabeled[i])
                                                                        }#end for (i in (1:length(indexLabeled)))
                                                                        
                                                                        indexLabeled<<-indexLabeledaux 
                                                                }#end if(tclvalue(mbval)=="no")
                                                                
                                                                if(tclvalue(mbval)=="cancel"){
                                                                        if ((cblugares=="1")|(cblablugares=="1")|(cccaVal=="COIA")){
                                                                                
                                                                                textos[indexClosest,dim1] <<- anteriorx
                                                                                textos[indexClosest,dim2] <<- anteriory
                                                                        }else{
                                                                                textosr[indexClosest,dim1] <<- anteriorx
                                                                                textosr[indexClosest,dim2] <<- anteriory
                                                                        }#end if (cblugares=="1")
                                                                }#end if(tclvalue(mbval)=="cancel")
                                                                
                                                                tkrreplot(img)
                                                        }
                                                        
                                                }#end OnLeftClick.up <- function(x,y)
                                                
                                                
                                                OnLeftClick.move <- function(x,y)
                                                {
                                                        if(cccaVal!="COIA" | typecoia %in% c("all", "sg", "eg", "vg"))
                                                        {
                                                                xClick <- x
                                                                yClick <- y
                                                                width  = as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                                height = as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                                
                                                                xMin = parPlotSize[1] * width
                                                                xMax = parPlotSize[2] * width
                                                                yMin = parPlotSize[3] * height
                                                                yMax = parPlotSize[4] * height
                                                                
                                                                rangeX = usrCoords[2] - usrCoords[1]
                                                                rangeY = usrCoords[4] - usrCoords[3]
                                                                
                                                                imgXcoords = (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                                imgYcoords = (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                                
                                                                xClick <- as.numeric(xClick)+0.5
                                                                yClick <- as.numeric(yClick)+0.5
                                                                yClick <- height - yClick
                                                                
                                                                xPlotCoord = usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                                yPlotCoord = usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                                
                                                                if ((cblugares=="1")|(cblablugares=="1")|(cccaVal=="COIA")){
                                                                        textos[indexClosest,dim1]<<-xPlotCoord
                                                                        textos[indexClosest,dim2]<<-yPlotCoord
                                                                }else{
                                                                        textosr[indexClosest,dim1]<<-xPlotCoord
                                                                        textosr[indexClosest,dim2]<<-yPlotCoord
                                                                }#end if (cblugares=="1")
                                                                
                                                                tkrreplot(img) 
                                                        }
                                                }#end OnLeftClick.move <- function(x,y)
                                                
                                                
                                                OnLeftClick.down <- function(x,y)
                                                {
                                                        if(cccaVal!="COIA" | typecoia %in% c("all", "sg", "eg", "vg"))
                                                        {
                                                                anteriorx <- NULL
                                                                anteriory <- NULL
                                                                
                                                                xClick <- x
                                                                yClick <- y
                                                                width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                                height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                                
                                                                xMin <- parPlotSize[1] * width
                                                                xMax <- parPlotSize[2] * width
                                                                yMin <- parPlotSize[3] * height
                                                                yMax <- parPlotSize[4] * height
                                                                
                                                                rangeX <- usrCoords[2] - usrCoords[1]
                                                                rangeY <- usrCoords[4] - usrCoords[3]
                                                                
                                                                imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                                imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                                
                                                                xClick <- as.numeric(xClick)+0.5
                                                                yClick <- as.numeric(yClick)+0.5
                                                                yClick <- height - yClick
                                                                
                                                                xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                                yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                                
                                                                squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                                                indexClosest <<- which.min(squared.Distance) 
                                                                
                                                                if ((cblugares=="1")|(cblablugares=="1")|(cccaVal=="COIA")){
                                                                        anteriorx <<- textos[indexClosest,dim1]
                                                                        anteriory <<- textos[indexClosest,dim2]
                                                                }else{
                                                                        anteriorx <<- textosr[indexClosest,dim1]
                                                                        anteriory <<- textosr[indexClosest,dim2]
                                                                }#end if (cblugares=="1")
                                                        }
                                                }#end OnLeftClick.down <- function(x,y)
                                                
                                                OnRightClick <- function(x,y)
                                                {	
                                                        if(cccaVal!="COIA" | typecoia %in% c("all", "sg", "eg", "vg"))
                                                        {
                                                                xClick <- x
                                                                yClick <- y
                                                                width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
                                                                height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
                                                                
                                                                xMin <- parPlotSize[1] * width
                                                                xMax <- parPlotSize[2] * width
                                                                yMin <- parPlotSize[3] * height
                                                                yMax <- parPlotSize[4] * height
                                                                
                                                                rangeX <- usrCoords[2] - usrCoords[1]	
                                                                rangeY <- usrCoords[4] - usrCoords[3]
                                                                
                                                                imgXcoords <- (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
                                                                imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
                                                                
                                                                xClick <- as.numeric(xClick)+0.5
                                                                yClick <- as.numeric(yClick)+0.5
                                                                yClick <- height - yClick
                                                                
                                                                xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
                                                                yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
                                                                
                                                                labelClosestPointd(xClick,yClick,imgXcoords,imgYcoords)
                                                        }
                                                }#end OnRightClick <- function(x,y)
                                                
                                                tkbind(img, "<B1-Motion>",OnLeftClick.move)
                                                tkbind(img, "<ButtonPress-1>",OnLeftClick.down)
                                                tkbind(img, "<ButtonRelease-1>",OnLeftClick.up)
                                                tkconfigure(img,cursor="pencil")
                                                tkbind(img, "<Button-3>",OnRightClick)
                                                tkconfigure(img,cursor="pencil")
                                        }#end if (nejes > length(descom$d))
                                }#end Onaxis<-function()
                                
                                numaxis <- tclVar( 1 )
                                enumaxis <-tkentry(barvp,width="50",textvariable=numaxis)
                                but.axis <-tkbutton(barvp,text="Choose",command=Onaxis, bg= "lightblue", width=10, foreground = "navyblue")
                                tkpack(imgbar, expand="TRUE", fill="both")  
                                tkpack(tklabel(barvp,text="Select the number of axes:"),
                                       enumaxis,but.axis,expand="FALSE", side= "left", fill ="both")
                                tkfocus(barvp)
                                
                                
                        }#end Graphics <- function()
                        graphic.button <-tkbutton(framegraphic,text="Graph",command=Graphics, bg= "lightblue", width=20, foreground = "navyblue")
                        
                        tkpack(tlv,scrv,expand = "TRUE", side="left", fill = "both")
                        tkpack.configure(scrv,side="left")
                        tkpack(OK.butv,expand = "TRUE", side="left", fill = "both")
                        tkfocus(tt)
                        
                        #######Color#######################################
                        indicev<-as.numeric(tkcurselection(tlv))+1
                        colorv <- colvariables[indicev[1]]
                        canvasv <- tkcanvas(framecol21,width="57",height="20",bg=colorv)
                        
                        ChangeColorv <- function()
                        {
                                colorv <<- tclvalue(tcl("tk_chooseColor",initialcolor=colvariables[indicev[1]],title="Choose a color"))
                                
                                if (nchar(colorv)>0)
                                {
                                        tkconfigure(canvasv,bg=colorv)
                                        colvariables[indicev]<<-colorv
                                }#end if (nchar(colorv)>0)
                        }#end ChangeColorv <- function()
                        
                        ChangeColor.buttonv<- tkbutton(framecol22,text="Change Color",command=ChangeColorv,width=4)
                        tkpack(canvasv,expand = "TRUE", side="left", fill = "both")
                        tkpack(ChangeColor.buttonv,expand = "TRUE", side="left", fill = "both")
                        
                        ######Labels  ###################################
                        Namev <- textvariables[indicev[1]]
                        entry.Namev <-tkentry(framename21,width="10",textvariable=Namev, bg="white")
                        
                        OnOKlv <- function()
                        {
                                NameValv <<- tclvalue(Namev)
                                textvariables[indicev[1]] <<-NameValv
                                
                                #####Values of listbox###############################
                                for (i in 1:dimvambien[2])
                                {
                                        tkdelete(tlv,0)
                                }#end for (i in 1:dimvambien[2])
                                
                                for (i in (1:(dimvambien[2])))
                                {
                                        tkinsert(tlv,"end",textvariables[i])
                                }#end for (i in (1:(dimvambien[2])))
                        }#end OnOKlv <- function()
                        
                        OK.butlv <-tkbutton(framename22,text="Change label",command=OnOKlv,width=4)
                        tkbind(entry.Namev, "<Return>",OnOKlv)
                        tkpack(entry.Namev,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butlv,expand = "TRUE", side="left", fill = "both")
                        
                        ###### Sizes  ###################################
                        Cexv <- cexvariables[indicev[1]]
                        entry.Cexv <-tkentry(framecex21,width="10",textvariable=Cexv, bg="white")
                        
                        OnOKcv <- function()
                        {
                                NameCexv <<- tclvalue(Cexv)
                                cexvariables[indicev] <<-NameCexv
                        }#end OnOKcv <- function()
                        
                        OK.butcv <-tkbutton(framecex22,text="Change size",command=OnOKcv,width=4)
                        tkbind(entry.Cexv, "<Return>",OnOKlv)
                        tkpack(entry.Cexv,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butcv,expand = "TRUE", side="left", fill = "both")
                        
                        ######Symbols###################################
                        
                        tkpack(tklabel(frames21,text="      ",width=27),expand = "TRUE", side="left", fill = "both")
                        
                        ##### List of sites ###########################
                        
                        scrl <- tkscrollbar(framet3, repeatinterval=5, command=function(...)tkyview(tll,...))
                        tll<-tklistbox(framet3,height=6,width=42,selectmode="multiple",yscrollcommand=function(...)tkset(scrl,...),background="white")
                        tkpack(tklabel(frametext3,text=Namesit),side="left",expand = "TRUE",fill="both")
                        
                        for (i in (1:(dimespec[1])))
                        {
                                tkinsert(tll,"end",textlugares[i])
                        }#end for (i in (1:(dimespec[1])))
                        
                        tkselection.set(tll,0) #  Indexing starts at zero.
                        
                        OnOKl <- function()
                        {
                                Choicel <<- textespecies[as.numeric(tkcurselection(tll))+1]
                                
                                ##### Color of the selected variable #############
                                indicel<<-as.numeric(tkcurselection(tll))+1
                                colorl <- collugares[indicel[1]]
                                tkconfigure(canvasl,bg=colorl)
                                
                                ##### Text of the selected variable  ###############
                                
                                Namel <<- tclVar(textlugares[indicel[1]])
                                tkconfigure(entry.Namel,textvariable=Namel)
                                tkconfigure(comboBoxl,text=simlugares[indicel[1]])
                                
                                ##### Size of the selected variable  ###############
                                
                                Cexl <<- tclVar(cexlugares[indicel[1]])
                                tkconfigure(entry.Cexl,textvariable=Cexl)
                                tkconfigure(comboBoxl,text=simlugares[indicel[1]])
                        }#end OnOKl <- function()
                        
                        OK.butl <-tkbutton(frameok3,text="    OK    ",command=OnOKl)
                        tkpack(tll,scrl,expand = "TRUE", side="left", fill = "both")
                        tkpack.configure(scrl,side="left")
                        
                        tkpack(OK.butl,expand = "TRUE", side="left", fill = "both")
                        tkfocus(tt)
                        
                        #######Color#######################################
                        indicel<-as.numeric(tkcurselection(tll))+1
                        colorl <- collugares[indicel[1]]
                        canvasl <- tkcanvas(framecol31,width="57",height="20",bg=colorl)
                        
                        ChangeColorl <- function()
                        {
                                colorl <<- tclvalue(tcl("tk_chooseColor",initialcolor=collugares[indicel[1]],title="Choose a color"))
                                
                                if (nchar(colorl)>0)
                                {
                                        tkconfigure(canvasl,bg=colorl)
                                        collugares[indicel]<<-colorl
                                }#end if (nchar(colorl)>0)
                        }#end ChangeColorl <- function()
                        
                        ChangeColor.buttonl<- tkbutton(framecol32,text="Change Color",command=ChangeColorl,width=4)
                        tkpack(canvasl,expand = "TRUE", side="left", fill = "both")
                        tkpack(ChangeColor.buttonl,expand = "TRUE", side="left", fill = "both")
                        
                        ######Labels  ###################################
                        Namel <- textlugares[indicel[1]]
                        entry.Namel <-tkentry(framename31,width="10",textvariable=Namel, bg="white")
                        
                        OnOKll <- function()
                        {
                                NameVall <<- tclvalue(Namel)
                                textlugares[indicel[1]] <-NameVall
                                
                                #####Values of listbox###############################
                                for (i in 1:dimespec[1])
                                {
                                        tkdelete(tll,0)
                                }#end for (i in 1:dimespec[1])
                                
                                for (i in (1:(dimespec[1])))
                                {
                                        tkinsert(tll,"end",textlugares[i])
                                }#end for (i in (1:(dimespec[1])))
                        }#end OnOKll <- function()
                        
                        OK.butll <-tkbutton(framename32,text="Change label",command=OnOKll,width=4)
                        tkbind(entry.Namel, "<Return>",OnOKll)
                        tkpack(entry.Namel,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butll,expand = "TRUE", side="left", fill = "both")
                        
                        ######Labels  ###################################
                        Cexl <- cexlugares[indicel[1]]
                        entry.Cexl <-tkentry(framecex31,width="10",textvariable=Cexl, bg="white")
                        
                        OnOKcl <- function()
                        {
                                NameCexl <<- tclvalue(Cexl)
                                cexlugares[indicel] <<-NameCexl
                        }#end OnOKcl <- function()
                        
                        OK.butcl <-tkbutton(framecex32,text="Change size",command=OnOKcl,width=4)
                        tkbind(entry.Cexl, "<Return>",OnOKll)
                        tkpack(entry.Cexl,expand = "TRUE", side="left", fill = "both")
                        tkpack(OK.butcl,expand = "TRUE", side="left", fill = "both")
                        
                        ######Symbols###################################
                        
                        comboBoxl <- tkwidget(frames31,"ComboBox",editable=FALSE,values=symbolos,width=7, background="white")
                        
                        chang.syml <- function()
                        {
                                simChoicel <<-symbolos[as.numeric(tclvalue(tcl(comboBoxl,"getvalue")))+1]
                                simlugares[indicel] <<-simChoicel
                        }#end chang.syml <- function()
                        
                        Change.symboll <-tkbutton(frames32,text="Change symbol",command=chang.syml,width=4)
                        tkpack(comboBoxl,side="left",expand="TRUE", fill="both")
                        tkpack(Change.symboll,side="left",expand="TRUE", fill="both")
                        
                        tkpack(tklabel(frametext4,text=""),side="right",expand = "TRUE",fill="both")
                        tkpack(tklabel(framet4,text=""),side="right",expand = "TRUE",fill="both")
                        tkpack(tklabel(frameok4,text=""),side="right",expand = "TRUE",fill="both")
                        tkpack(tklabel(frames4,text=""),side="right",expand = "TRUE",fill="both")
                        
                        tkpack(tklabel(framehvtitle,text="Graph size"),side="right",expand = "TRUE",fill="both")
                        tkpack(tklabel(framehnames,text="Horizontal"),side="right",expand = "TRUE",fill="both")
                        tkpack(tklabel(framevnames,text="Vertical"),side="right",expand = "TRUE",fill="both")
                        
                        #####  Textbox to change the size of the graph window #####
                        entryvalueh <- tclVar(hescale)
                        entryh <-tkentry(framehtext,width="10",textvariable=entryvalueh, bg="white")
                        tkbind(entryh, "<Return>",Graphics)
                        
                        entryvaluev <- tclVar(vescale)
                        entryv <-tkentry(framevtext,width="10",textvariable=entryvaluev, bg="white")
                        tkbind(entryv, "<Return>",Graphics)
                        
                        tkpack(entryh, entryv, expand = "TRUE",side="top", fill="both")
                        
                        tkpack(tklabel(framecol41,text="Show axes"),expand = "TRUE", side="left", fill = "both")
                        tkpack(cb,expand = "TRUE", side="left",expand="TRUE", fill = "both")
                        
                        tkpack(tklabel(framename41,text=paste("Show", Namesit)),expand = "TRUE", side="left", fill = "both")
                        tkpack(cbl,expand = "TRUE", side="left", fill = "both")
                        
                        tkpack(tklabel(framecex41,text=paste("Show labels for", Namesit)),expand = "TRUE", side="left", fill = "both")
                        tkpack(cbll,expand = "TRUE", side="left", fill = "both")
                        tkpack(tklabel(frames41,text="      ",width=27),expand = "TRUE", side="left", fill = "both")
                        
                        tkpack(framecol41,framecol42, side="left", expand = "TRUE", fill="both")
                        if(cccaVal!="COIA")
                        {
                                tkpack(framename41,framename42, side="left", expand = "TRUE", fill="both")
                                tkpack(framecex41,framecex42, side="left", expand = "TRUE", fill="both")
                        }
                        #tkpack(frames41, side="left", expand = "TRUE", fill="both")
                        
                        tkpack(framecol11,framecol12, side="left", expand = "TRUE", fill="both")
                        tkpack(framename11,framename12, side="left", expand = "TRUE", fill="both")
                        tkpack(framecex11,framecex12, side="left", expand = "TRUE", fill="both")
                        tkpack(frames11,frames12, side="left", expand = "TRUE", fill="both")
                        tkpack(framecol21,framecol22, side="left", expand = "TRUE", fill="both")
                        tkpack(framename21,framename22, side="left", expand = "TRUE", fill="both")
                        tkpack(framecex21,framecex22, side="left", expand = "TRUE", fill="both")
                        tkpack(frames21, side="left", expand = "TRUE", fill="both")
                        tkpack(framecol31,framecol32, side="left", expand = "TRUE", fill="both")
                        tkpack(framename31,framename32, side="left", expand = "TRUE", fill="both")
                        tkpack(framecex31,framecex32, side="left", expand = "TRUE", fill="both")
                        tkpack(frames31,frames32, side="left", expand = "TRUE", fill="both")
                        
                        tkpack(frametext1,framet1,frameok1,framecol1,framename1,framecex1,frames1,expand = "TRUE", fill="both")
                        tkpack(frametext2,framet2,frameok2,framecol2,framename2,framecex2,frames2,expand = "TRUE", fill="both")
                        tkpack(frametext3,framet3,frameok3,framecol3,framename3,framecex3,frames3,expand = "TRUE", fill="both")
                        tkpack(framecol4,framename4,framecex4,expand = "FALSE", fill="both")
                        tkpack(framehnames, framevnames, expand = "FALSE", fill="both")
                        tkpack(framehtext, framevtext, expand = "FALSE", fill="both")
                        tkpack(framehvnames, framehvtext,expand = "FALSE",side="left", fill="both")
                        tkpack(framehvtitle, framehv, expand = "FALSE", fill="both")
                        tkpack(frametext4,framet4,frameok4,frames4,framett4aux,framett4auxgs, expand = "FALSE", fill="both")
                        
                        tkpack(graphic.button, expand="TRUE", side= "left", fill ="both")
                        
                        tkpack(framett1,framett2,framett3, framett4,expand = "TRUE",side="left", fill="both")
                        tkpack(framett,framegraphic,expand = "TRUE",side="top", fill="y")
                }#end OnOKnames<-function()
                
                framenames<-tkframe(wnames, relief = "flat", borderwidth = 2, background = "white")
                framee<-tkframe(framenames, relief = "flat", borderwidth = 2, background = "white")
                framel<-tkframe(framenames, relief = "flat", borderwidth = 2, background = "white")
                frameoknames<-tkframe(wnames, relief = "flat", borderwidth = 2, background = "white")
                
                OK.butnames<-tkbutton(frameoknames,text="   OK   ", command=OnOKnames,  bg= "lightblue", width=20, foreground = "navyblue")
                tkbind(OK.butnames, "<Return>",OnOKnames)
                
                sitesname<-tclVar("Sites")
                speciesname<-tclVar("Species")
                variablesname<-tclVar("Environmental v")
                # mixto<-tclVar("NO")
                ncuant<-tclVar(dim(fvambientales)[2])
                nmixto<-tclVar(0)
                
                entry.Namesit <-tkentry(framee,width="50",textvariable=sitesname, bg="white")
                tkbind(entry.Namesit, "<Return>",OnOKnames)
                
                entry.Namespec <-tkentry(framee,width="50",textvariable=speciesname, bg="white")
                tkbind(entry.Namespec, "<Return>",OnOKnames)
                
                entry.Namev <-tkentry(framee,width="50",textvariable=variablesname, bg="white")
                tkbind(entry.Namev, "<Return>",OnOKnames)
                
                if (cccaVal!="COIA")
                {
                        
                        #                         entry.Namemix <-tkentry(framee,width="50",textvariable=mixto, bg="white")
                        #                         tkbind(entry.Namemix, "<Return>",OnOKnames)
                        #                         
                        entry.Namencuant <-tkentry(framee,width="50",textvariable=ncuant, bg="white")
                        tkbind(entry.Namencuant, "<Return>",OnOKnames)
                        
                        entry.Namenmix <-tkentry(framee,width="50",textvariable=nmixto, bg="white")
                        tkbind(entry.Namenmix, "<Return>",OnOKnames)
                        tkpack(tklabel(framel,text="Names for sites:"),
                               tklabel(framel,text="Names for species:"),
                               tklabel(framel,text="Names for environmental v:"),
                               #   tklabel(framel,text="Mixed environmental variables?:"),
                               tklabel(framel,text="Number of cuantitative variables:"),
                               tklabel(framel,text="Number of categorical variables:"), expand = "TRUE", side="top", fill = "both")
                        tkpack(OK.butnames)
                        tkpack(entry.Namesit, entry.Namespec, entry.Namev, #entry.Namemix, 
                               entry.Namencuant, entry.Namenmix,
                               expand = "TRUE",side="top", fill="both")
                        
                }else{
                        tkpack(tklabel(framel,text="Names for sites:"),
                               tklabel(framel,text="Names for species:"),
                               tklabel(framel,text="Names for environmental v:"),
                               expand = "TRUE", side="top", fill = "both")
                        tkpack(OK.butnames)
                        tkpack(entry.Namesit, entry.Namespec, entry.Namev,
                               expand = "TRUE",side="top", fill="both")
                        
                }
                tkpack(framel, framee, expand = "TRUE",side="left", fill="y")
                tkpack(framenames, frameoknames, expand = "TRUE",side="top", fill="y")  
                tkfocus(wnames)
        }#end OnOKinf <- function()
        
        OK.butinf <-tkbutton(winfor,text="   OK   ",command=OnOKinf, bg= "lightblue", width=20, foreground = "navyblue")
        tkbind(OK.butinf, "<Return>",OnOKinf)
        
        
        tkpack(OK.butinf)
        tkfocus(winfor)
}#end function
