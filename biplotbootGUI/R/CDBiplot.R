#indexr<-sample(1:1599,500)
#indexw<-sample(1:4898,500)
#vinos<-rbind(winequality.red[indexr,], winequality.white[indexw,])

CDBiplot<-function(data, clase)
{
        cbiplotint<-function(data, P, Q, tol, iter, times, clase, showgr)
        {
                
                Fmax <- 0 
                itermax <- 0
                for(t in 1:times)
                {
                        I<-dim(data)[1]
                        J<-dim(data)[2]
                        
                        data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normlized data, with variance divided by (I-1)
                        X <- as.matrix(data.sd*sqrt(I/(I-1)))  # matrix of normlized data (with var divided by I)
                        
                        #matriz de pertenencia a cluster
                        # if (randallo==FALSE)
                        #        {
                        #               U0 <- U 
                        #      }else{
                        U0<-array(rep(0, I*P), dim=c(I,P))
                        
                        for(i in 1:I)
                        {
                                p<-sample(1:P,1)
                                U0[i,p]<-1
                        } 
                        
                        
                        sumaU <- colSums(U0)
                        while (sum(sumaU == 0) > 0) {
                                ind.max <- which.max(sumaU)
                                su.max <- max(sumaU)
                                ind.min <- which.min(sumaU)
                                su.min <- min(sumaU)
                                ind.nzU <- which(U0[, ind.max] == 1)
                                ind.sel <- ind.nzU[1:floor(su.max)/2]
                                U0[ind.sel, ind.min] <- 1
                                U0[ind.sel, ind.max] <- 0
                                sumaU <- colSums(U0)
                        }  # end while
                        
                        ## matriz de centroides
                        Xc0 <- ginv((t(U0)%*%U0))%*%t(U0)%*%X
                        ##matriz que contiene centroides en lugar de individuos
                        Z0 <- U0%*%Xc0
                        if(I>=dim(Z0)[1])
                        {
                                desW0<-svd(Z0)
                                ## coordenadas hj biplot para las variables
                                B0<-desW0$v[,1:Q]%*%diag(desW0$d[1:Q])
                                L0<-desW0$d[1:Q]
                        }else{
                                W0p <- Z0%*%t(Z0)
                                B0 <- t(Z0)%*%eigen(W0p)$vectors[, 1:Q]
                                L0 <- sqrt(eigen(W0p)$values[1:Q])
                        }
                        
                        A0 <-X%*%B0%*%solve(diag(L0))
                        ### coordenadas de centroides en el espacio reducido
                        Ac0 <- Xc0%*%B0%*%solve(diag(L0))
                        
                        F0 <-  sum(diag(t(U0%*%Ac0%*%solve(diag(L0))%*%t(B0))%*%(U0%*%Ac0%*%solve(diag(L0))%*%t(B0))))
                        dev0 <- F0/sum(diag(t(A0)%*%A0))
                        
                        Fk <- F0
                        A <- A0
                        Ac <- Ac0
                        B <- B0
                        
                        
                        for(k in 1:iter)
                        {
                                print(paste("N of repetition", t, "iteration", k, sep=" "))
                                U <- matrix(0, I, P)
                                # update U
                                for (i in 1:I)
                                {
                                        dist <- rep(0, P)
                                        for (p in 1:P) 
                                        {          
                                                dist[p] <- sum((A[i, ]-Ac[p, ])^2)
                                        } # end for 
                                        min.dist <- which.min(dist)
                                        U[i, min.dist] <- 1
                                }  # end for
                                sumaU <- colSums(U)
                                while (sum(sumaU == 0) > 0) {
                                        ind.max <- which.max(sumaU)
                                        su.max <- max(sumaU)
                                        ind.min <- which.min(sumaU)
                                        su.min <- min(sumaU)
                                        ind.nzU <- which(U[, ind.max] == 1)
                                        ind.sel <- ind.nzU[1:floor(su.max)/2]
                                        U[ind.sel, ind.min] <- 1
                                        U[ind.sel, ind.max] <- 0
                                        sumaU <- colSums(U)
                                }  # end while
                                
                                
                                Xc <- ginv((t(U)%*%U))%*%t(U)%*%X
                                Z <- U%*%Xc
                                
                                
                                desW<-svd(Z)
                                ## coordenadas hj biplot para las variables
                                if(I>=dim(Z)[1])
                                {
                                        B <- desW$v[,1:Q]%*%diag(desW$d[1:Q])
                                        L <- desW$d[1:Q]
                                }else{
                                        W0p <- Z%*%t(Z)
                                        B <- t(Z)%*%eigen(W0p)$vectors[, 1:Q]
                                        L <- sqrt(eigen(W0p)$values[1:Q])
                                }
                                
                                A <-X%*%B%*%solve(diag(L))
                                
                                ### coordenadas de centroides en el espacio reducido
                                Ac <- Xc%*%B%*%solve(diag(L))
                                
                                Fk1<-sum(diag(t(U%*%Ac%*%solve(diag(L))%*%t(B))%*%(U%*%Ac%*%solve(diag(L))%*%t(B))))
                                dev <- Fk1/sum(diag(t(A)%*%A))
                                
                                if(abs(Fk1-Fk)<tol & (Fk1>Fmax))
                                {
                                        if(k>itermax)
                                                itermax<-k
                                        break   
                                }else{
                                        Fk<-Fk1
                                }
                        }
                } 
                
                varA <- var(A)
                vp <- diag(varA)/sum(diag(var(X)))
                orden <- order(vp, decreasing = TRUE)
                varexp <- data.frame(1:Q, vp[orden])
                varexp[,2]<-round(varexp[,2]*100, 2)
                colnames(varexp) <- c("Axis", "Expl. Var (%)")
                Aorden <- A[,orden]
                Acorden <- Ac[,orden]
                Borden <-B[,orden]
                Lorden <- diag(L[orden])
                
                Umax <- U
                Amax <- Aorden
                Acmax <- Acorden
                Fmax <- Fk1
                devmax <- dev
                varexpmax <- varexp
                Lmax <-Lorden
                
                tabpseu<-as.data.frame(Umax%*%matrix(1:ncol(Umax)))
                colnames(tabpseu)<-"CBiplot class."
                if(exists("clase"))
                {
                        clase<-data.frame(clase)
                        colnames(clase)<-"Real Class"
                        tabaux<-as.data.frame(c(clase, tabpseu))
                        tabpseu<-table(tabaux)       
                }
                resultados<-list(U=Umax, B=Borden,A=Amax,Ac=Acmax,Fk=Fmax,expvar=varexpmax,Lmax=Lmax, classif=tabpseu, itermax=itermax)
                print(resultados)
                cat("File saved in:    ",file="Results.txt")
                cat(getwd(),file="temp.txt")                        		
                file.append("Results.txt","temp.txt")	
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")		
                cat("Row clasification:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                write.table(Umax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Eigenvalues:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Lmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Explained variability by each component:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(varexpmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                cat("Row coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Amax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Centroid coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Acmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Variable coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Borden,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                        		
                file.append("Results.txt","temp.txt")	
                cat("Pseudoconfusion matrix:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(tabpseu,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum number of iterations:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(itermax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum F:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(Fmax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                file.show("Results.txt")
                file.remove("temp.txt")       
                return(resultados)
        }
        
        
        
        
        dbiplotint<-function(data, Q, tol, iter, times, showgr)
        {
                
                
                calc_biplot<-function(pos)
                {
                        if(length(pos)>0)
                        {
                                
                                W <- X[,pos,drop=FALSE]
                                if(I>=dim(W)[1])
                                {
                                        desW<-svd(W)
                                        L<-desW$d[1]
                                        ## coordenadas hj biplot para las variables
                                        b<-desW$v[,1]*desW$d[1]      
                                }else{
                                        Wp <- W%*%t(W)
                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                        L <- sqrt(eigen(Wp)$values[1])
                                }
                                
                                return(list(b, L))
                        }
                        
                }
                
                Fmax <- 0  
                itermax <- 0
                
                for(t in 1:times)
                {
                        I<-dim(data)[1]
                        J<-dim(data)[2]
                        
                        data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normlized data, with variance divided by (I-1)
                        X <- as.matrix(data.sd*sqrt(I/(I-1)))  # matrix of normlized data (with var divided by I)
                        
                        ##matriz de pertenencia a eje
                        V0<-array(rep(0, J*Q), dim=c(J,Q))
                        
                        for(j in 1:J)
                        {
                                q<-sample(1:Q,1)
                                V0[j,q]<-1
                        }
                        
                        
                        sumaV <- colSums(V0)
                        while (sum(sumaV == 0) > 0) {
                                ind.max <- which.max(sumaV)
                                su.max <- max(sumaV)
                                ind.min <- which.min(sumaV)
                                su.min <- min(sumaV)
                                ind.nzV <- which(V0[, ind.max] == 1)
                                ind.sel <- ind.nzV[1:floor(su.max)/2]
                                V0[ind.sel, ind.min] <- 1
                                V0[ind.sel, ind.max] <- 0
                                sumaV <- colSums(V0)
                        }  # end while
                        ## espacio donde proyectar
                        B0 <-V0
                        L0<-rep(0,Q)
                        for(i in 1:Q)
                        {
                                ## para cada eje se extrae la submatriz de variables que van a cargar en el
                                W0<-X[,which(V0[,i]>0), drop=FALSE]
                                if(dim(W0)[2]>0)
                                {
                                        
                                        ## coordenadas hj biplot para las variables
                                        
                                        if(I>=dim(W0)[1])
                                        {
                                                desW0 <- svd(W0)
                                                b0 <- desW0$v[,1]*desW0$d[1]
                                                L0[i] <- desW0$d[1]
                                                
                                        }else{
                                                W0p <- W0%*%t(W0)
                                                b0 <- t(W0)%*%eigen(W0p)$vectors[,1]
                                                L0[i] <- sqrt(eigen(W0p)$values[1])
                                        }
                                        B0[which(B0[,i]>0),i]<-b0    
                                        
                                }
                        }
                        
                        A0 <-X%*%B0%*%solve(diag(L0))
                        
                        F0 <-  sum(diag(t(A0%*%solve(diag(L0))%*%t(B0))%*%(A0%*%solve(diag(L0))%*%t(B0))))
                        dev0 <- F0/sum(diag(t(A0)%*%A0))
                        
                        
                        V <- V0
                        Fk <- F0
                        A <- A0
                        B <- B0
                        
                        
                        for(k in 1:iter)
                        {
                                print(paste("N of repetition", t, "iteration", k, sep=" "))
                                
                                B<-array(rep(0, J*Q), dim=c(J,Q))
                                
                                
                                V[1,1]<-1
                                V[1,-1]<-0
                                
                                sumaV <- colSums(V)
                                while (sum(sumaV == 0) > 0) {
                                        ind.max <- which.max(sumaV)
                                        su.max <- max(sumaV)
                                        ind.min <- which.min(sumaV)
                                        su.min <- min(sumaV)
                                        ind.nzV <- which(V[, ind.max] == 1)
                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                        V[ind.sel, ind.min] <- 1
                                        V[ind.sel, ind.max] <- 0
                                        sumaV <- colSums(V)
                                }  # end while
                                
                                Baux<-V
                                Laux<-rep(0,Q)
                                for(i in 1:Q)
                                {
                                        W<-X[,which(V[,i]>0),drop=FALSE]
                                        if(dim(W)[2]>0)
                                        {
                                                if(I>=dim(W)[1])
                                                {
                                                        desW <- svd(W)
                                                        b <- desW$v[,1]*desW$d[1]
                                                        Laux[i] <- desW$d[1]
                                                        
                                                }else{
                                                        Wp <- W%*%t(W)
                                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                                        Laux[i] <- sqrt(eigen(Wp)$values[1])
                                                }
                                                Baux[,i]<-0
                                                Baux[which(V[,i]>0),i]<-b    
                                                
                                        }
                                }
                                
                                for(j in 1:J)
                                {
                                        f<-c()
                                        for(q in 1:Q)
                                        {
                                                V[j,-q]<-0
                                                V[j,q]<-1
                                                
                                                sumaV <- colSums(V)
                                                while (sum(sumaV == 0) > 0) {
                                                        ind.max <- which.max(sumaV)
                                                        su.max <- max(sumaV)
                                                        ind.min <- which.min(sumaV)
                                                        su.min <- min(sumaV)
                                                        ind.nzV <- which(V[, ind.max] == 1)
                                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                                        V[ind.sel, ind.min] <- 1
                                                        V[ind.sel, ind.max] <- 0
                                                        sumaV <- colSums(V)
                                                }  # end while
                                                
                                                if(q!=1)
                                                {
                                                        indices<-c(q-1,q)
                                                        calc_pos<-apply(V[,indices], 2, posic<-function(v){which(v>0)})
                                                        calc_b<-vector("list", 2)
                                                        calc_b<-lapply(calc_pos, calc_biplot)
                                                        for (i in 1:2)
                                                        {
                                                                if (!is.null(calc_b[[i]][[1]]))
                                                                {
                                                                        ind<-indices[i]
                                                                        Baux[,ind]<-0
                                                                        Laux[ind]<-calc_b[[i]][[2]]
                                                                        Baux[which(V[,ind]>0),ind]<-calc_b[[i]][[1]]
                                                                }
                                                        }
                                                        
                                                }
                                                
                                                A <-X%*%Baux%*%solve(diag(Laux))
                                                f <- c(f,sum(diag(t(A%*%solve(diag(Laux))%*%t(Baux))%*%(A%*%solve(diag(Laux))%*%t(Baux)))))           
                                        }
                                        V[j,]<-0
                                        if(length(which(f==max(f)))==1)
                                        {        
                                                V[j,which(f==max(f))]<-1
                                        }else{
                                                V[j,which(f==max(f))[1]]<-1
                                        }
                                        
                                }
                                
                                sumaV <- colSums(V)
                                while (sum(sumaV == 0) > 0) {
                                        ind.max <- which.max(sumaV)
                                        su.max <- max(sumaV)
                                        ind.min <- which.min(sumaV)
                                        su.min <- min(sumaV)
                                        ind.nzV <- which(V[, ind.max] == 1)
                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                        V[ind.sel, ind.min] <- 1
                                        V[ind.sel, ind.max] <- 0
                                        sumaV <- colSums(V)
                                }  # end while
                                B<-V
                                L<-rep(0,Q)
                                for(i in 1:Q)
                                {
                                        W<-X[,which(V[,i]>0), drop=FALSE]
                                        if(dim(W)[2]>0)
                                        {
                                                if(I>=dim(W)[1])
                                                {
                                                        desW <- svd(W)
                                                        b <- desW$v[,1]*desW$d[1]
                                                        L[i] <- desW$d[1]
                                                        
                                                }else{
                                                        Wp <- W%*%t(W)
                                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                                        L[i] <- sqrt(eigen(Wp)$values[1])
                                                }
                                                B[which(B[,i]>0),i]<-b    
                                                
                                        }
                                }
                                
                                A <-X%*%B%*%solve(diag(L))
                                
                                
                                
                                Fk1<-sum(diag(t(A%*%solve(diag(L))%*%t(B))%*%(A%*%solve(diag(L))%*%t(B))))
                                dev <- Fk1/sum(diag(t(A)%*%A))
                                
                                if(abs(Fk1-Fk)<tol & (Fk1>Fmax))
                                {
                                        if(k>itermax)
                                                itermax<-k
                                        break   
                                }else{
                                        Fk<-Fk1
                                }
                        }
                } 
                
                varA <- var(A)
                vp <- diag(varA)/sum(diag(var(X)))
                orden <- order(vp, decreasing = TRUE)
                varexp <- data.frame(1:Q, vp[orden])
                varexp[,2]<-round(varexp[,2]*100, 2)
                colnames(varexp) <- c("Axis", "Expl. Var (%)")
                Aorden <- A[,orden]
                Vorden <- V[,orden]
                Borden <-B[,orden]
                Lorden <- diag(L[orden])
                
                Vmax <- Vorden
                Amax <- Aorden
                Fmax <- Fk1
                devmax <- dev
                varexpmax <- varexp
                Lmax <-Lorden
                
                resultados<-list(V=Vmax, B=Borden,A=Amax,Fk=Fmax,expvar=varexpmax,Lmax=Lmax, cor=cor(Amax), itermax=itermax)
                print(resultados)
                cat("File saved in:    ",file="Results.txt")
                cat(getwd(),file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")		
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Variable clasification:\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                write.table(Vmax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Eigenvalues:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Lmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Explained variability by each component:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(varexpmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                cat("Row coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Amax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Variable coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Borden,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Correlations matrix:\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                write.table(round(cor(Amax),digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum number of iterations:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(itermax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum F:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(Fmax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                file.show("Results.txt")
                file.remove("temp.txt")
                return(resultados)
        }
        
        
        
        cdbiplotint<-function(data, P, Q, tol, iter, times, clase, showgr)
        {
                calc_biplot<-function(pos)
                {
                        if(length(pos)>0)
                        {
                                
                                W <- X[,pos,drop=FALSE]
                                if(I>=dim(W)[1])
                                {
                                        desW<-svd(W)
                                        L<-desW$d[1]
                                        ## coordenadas hj biplot para las variables
                                        b<-desW$v[,1]*desW$d[1]      
                                }else{
                                        Wp <- W%*%t(W)
                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                        L <- sqrt(eigen(Wp)$values[1])
                                }
                                
                                return(list(b, L))
                        }
                        
                }
                Fmax <- 0  
                itermax <- 0
                for(t in 1:times)
                {
                        I<-dim(data)[1]
                        J<-dim(data)[2]
                        
                        data.sd <- data.frame(scale(data, center = TRUE, scale = TRUE))  # normlized data, with variance divided by (I-1)
                        X <- as.matrix(data.sd*sqrt(I/(I-1)))  # matrix of normlized data (with var divided by I)
                        
                        ##matriz de pertenencia a cluster
                        # if (randallo==FALSE)
                        #        {
                        #               U0 <- U 
                        #      }else{
                        U0<-array(rep(0, I*P), dim=c(I,P))
                        
                        for(i in 1:I)
                        {
                                p<-sample(1:P,1)
                                U0[i,p]<-1
                        } 
                        #     }   
                        
                        sumaU <- colSums(U0)
                        while (sum(sumaU == 0) > 0) {
                                ind.max <- which.max(sumaU)
                                su.max <- max(sumaU)
                                ind.min <- which.min(sumaU)
                                su.min <- min(sumaU)
                                ind.nzU <- which(U0[, ind.max] == 1)
                                ind.sel <- ind.nzU[1:floor(su.max)/2]
                                U0[ind.sel, ind.min] <- 1
                                U0[ind.sel, ind.max] <- 0
                                sumaU <- colSums(U0)
                        }  # end while
                        
                        ##matriz de pertenencia a eje
                        V0<-array(rep(0, J*Q), dim=c(J,Q))
                        
                        for(j in 1:J)
                        {
                                q<-sample(1:Q,1)
                                V0[j,q]<-1
                        }
                        
                        ## matriz de centroides
                        Xc0 <- ginv((t(U0)%*%U0))%*%t(U0)%*%X
                        ##matriz que contiene centroides en lugar de individuos
                        Z0 <- U0%*%Xc0
                        
                        sumaV <- colSums(V0)
                        while (sum(sumaV == 0) > 0) {
                                ind.max <- which.max(sumaV)
                                su.max <- max(sumaV)
                                ind.min <- which.min(sumaV)
                                su.min <- min(sumaV)
                                ind.nzV <- which(V0[, ind.max] == 1)
                                ind.sel <- ind.nzV[1:floor(su.max)/2]
                                V0[ind.sel, ind.min] <- 1
                                V0[ind.sel, ind.max] <- 0
                                sumaV <- colSums(V0)
                        }  # end while
                        ## espacio donde proyectar
                        B0 <-V0
                        L0<-rep(0,Q)
                        for(i in 1:Q)
                        {
                                ## para cada eje se extrae la submatriz de variables que van a cargar en el
                                W0<-Z0[,which(V0[,i]>0), drop=FALSE]
                                if(dim(W0)[2]>0)
                                {
                                        if(I>=dim(W0)[1])
                                        {
                                                desW0<-svd(W0)
                                                ## coordenadas hj biplot para las variables
                                                b0<-desW0$v[,1]*desW0$d[1]
                                                L0[i]<-desW0$d[1]
                                        }else{
                                                W0p <- W0%*%t(W0)
                                                b0 <- t(W0)%*%eigen(W0p)$vectors[,1]
                                                L0[i] <- sqrt(eigen(W0p)$values[1])    
                                        }        
                                        B0[which(B0[,i]>0),i]<-b0                
                                }
                        }
                        
                        A0 <-X%*%B0%*%solve(diag(L0))
                        ### coordenadas de centroides en el espacio reducido
                        Ac0 <- Xc0%*%B0%*%solve(diag(L0))
                        
                        F0 <-  sum(diag(t(U0%*%Ac0)%*%(U0%*%Ac0)))
                        dev0 <- F0/sum(diag(t(A0)%*%A0))
                        
                        
                        V <- V0
                        Fk <- F0
                        A <- A0
                        Ac <- Ac0
                        B <- B0
                        
                        
                        for(k in 1:iter)
                        {
                                print(paste("N of repetition", t, "iteration", k, sep=" "))
                                U <- matrix(0, I, P)
                                # update U
                                for (i in 1:I)
                                {
                                        dist <- rep(0, P)
                                        for (p in 1:P) 
                                        {          
                                                dist[p] <- sum((A[i, ]-Ac[p, ])^2)
                                        } # end for 
                                        min.dist <- which.min(dist)
                                        U[i, min.dist] <- 1
                                }  # end for
                                sumaU <- colSums(U)
                                while (sum(sumaU == 0) > 0) {
                                        ind.max <- which.max(sumaU)
                                        su.max <- max(sumaU)
                                        ind.min <- which.min(sumaU)
                                        su.min <- min(sumaU)
                                        ind.nzU <- which(U[, ind.max] == 1)
                                        ind.sel <- ind.nzU[1:floor(su.max)/2]
                                        U[ind.sel, ind.min] <- 1
                                        U[ind.sel, ind.max] <- 0
                                        sumaU <- colSums(U)
                                }  # end while
                                
                                
                                Xc <- ginv((t(U)%*%U))%*%t(U)%*%X
                                Z <- U%*%Xc
                                
                                B<-array(rep(0, J*Q), dim=c(J,Q))
                                
                                
                                V[1,1]<-1
                                V[1,-1]<-0
                                
                                sumaV <- colSums(V)
                                while (sum(sumaV == 0) > 0) {
                                        ind.max <- which.max(sumaV)
                                        su.max <- max(sumaV)
                                        ind.min <- which.min(sumaV)
                                        su.min <- min(sumaV)
                                        ind.nzV <- which(V[, ind.max] == 1)
                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                        V[ind.sel, ind.min] <- 1
                                        V[ind.sel, ind.max] <- 0
                                        sumaV <- colSums(V)
                                }  # end while
                                
                                Baux<-V
                                Laux<-rep(0,Q)
                                for(i in 1:Q)
                                {
                                        W<-Z[,which(V[,i]>0),drop=FALSE]
                                        if(dim(W)[2]>0)
                                        {
                                                if(I>=dim(W)[1])
                                                {
                                                        desW<-svd(W)
                                                        ## coordenadas hj biplot para las variables
                                                        b<-desW$v[,1]*desW$d[1]
                                                        Laux[i]<-desW$d[1]
                                                }else{
                                                        Wp <- W%*%t(W)
                                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                                        Laux[i] <- sqrt(eigen(Wp)$values[1])  
                                                }
                                                
                                                Baux[,i]<-0
                                                Baux[which(V[,i]>0),i]<-b
                                        }
                                }
                                
                                for(j in 1:J)
                                {
                                        f<-c()
                                        for(q in 1:Q)
                                        {
                                                V[j,-q]<-0
                                                V[j,q]<-1
                                                
                                                sumaV <- colSums(V)
                                                while (sum(sumaV == 0) > 0) {
                                                        ind.max <- which.max(sumaV)
                                                        su.max <- max(sumaV)
                                                        ind.min <- which.min(sumaV)
                                                        su.min <- min(sumaV)
                                                        ind.nzV <- which(V[, ind.max] == 1)
                                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                                        V[ind.sel, ind.min] <- 1
                                                        V[ind.sel, ind.max] <- 0
                                                        sumaV <- colSums(V)
                                                }  # end while
                                                
                                                if(q!=1)
                                                {
                                                        indices<-c(q-1,q)
                                                        calc_pos<-apply(V[,indices], 2, posic<-function(v){which(v>0)})
                                                        calc_b<-vector("list", 2)
                                                        calc_b<-lapply(calc_pos, calc_biplot)
                                                        for (i in 1:2)
                                                        {
                                                                if (!is.null(calc_b[[i]][[1]]))
                                                                {
                                                                        ind<-indices[i]
                                                                        Baux[,ind]<-0
                                                                        Laux[ind]<-calc_b[[i]][[2]]
                                                                        Baux[which(V[,ind]>0),ind]<-calc_b[[i]][[1]]
                                                                }
                                                        }
                                                        
                                                }
                                                
                                                A <-X%*%Baux%*%solve(diag(Laux))
                                                Ac <-Xc%*%Baux%*%solve(diag(Laux))
                                                f <- c(f,sum(diag(t(U%*%Ac)%*%(U%*%Ac))))           
                                        }
                                        V[j,]<-0
                                        if(length(which(f==max(f)))==1)
                                        {        
                                                V[j,which(f==max(f))]<-1
                                        }else{
                                                V[j,which(f==max(f))[1]]<-1
                                        }
                                        
                                }
                                
                                sumaV <- colSums(V)
                                while (sum(sumaV == 0) > 0) {
                                        ind.max <- which.max(sumaV)
                                        su.max <- max(sumaV)
                                        ind.min <- which.min(sumaV)
                                        su.min <- min(sumaV)
                                        ind.nzV <- which(V[, ind.max] == 1)
                                        ind.sel <- ind.nzV[1:floor(su.max)/2]
                                        V[ind.sel, ind.min] <- 1
                                        V[ind.sel, ind.max] <- 0
                                        sumaV <- colSums(V)
                                }  # end while
                                B<-V
                                L<-rep(0,Q)
                                for(i in 1:Q)
                                {
                                        W<-Z[,which(V[,i]>0), drop=FALSE]
                                        if(dim(W)[2]>0)
                                        {
                                                if(I>=dim(W)[1])
                                                {
                                                        desW<-svd(W)
                                                        ## coordenadas hj biplot para las variables
                                                        b<-desW$v[,1]*desW$d[1]
                                                        L[i]<-desW$d[1]
                                                }else{
                                                        Wp <- W%*%t(W)
                                                        b <- t(W)%*%eigen(Wp)$vectors[,1]
                                                        L[i] <- sqrt(eigen(Wp)$values[1])  
                                                }
                                                
                                                B[which(B[,i]>0),i]<-b
                                                
                                                
                                        }
                                }
                                
                                A <-X%*%B%*%solve(diag(L))
                                
                                ### coordenadas de centroides en el espacio reducido
                                Ac <- Xc%*%B%*%solve(diag(L))
                                
                                
                                Fk1<-sum(diag(t(U%*%Ac)%*%(U%*%Ac)))
                                dev <- Fk1/sum(diag(t(A)%*%A))
                                
                                if(abs(Fk1-Fk)<tol & (Fk1>Fmax))
                                {
                                        if(k>itermax)
                                                itermax<-k
                                        break   
                                }else{
                                        Fk<-Fk1
                                }
                        }
                } 
                
                varA <- var(A)
                vp <- diag(varA)/sum(diag(var(X)))
                orden <- order(vp, decreasing = TRUE)
                varexp <- data.frame(1:Q, vp[orden])
                varexp[,2]<-round(varexp[,2]*100, 2)
                colnames(varexp) <- c("Axis", "Expl. Var (%)")
                Aorden <- A[,orden]
                Acorden <- Ac[,orden]
                Vorden <- V[,orden]
                Borden <-B[,orden]
                Lorden <- diag(L[orden])
                
                Umax <- U
                Vmax <- Vorden
                Amax <- Aorden
                Acmax <- Acorden
                Fmax <- Fk1
                devmax <- dev
                varexpmax <- varexp
                Lmax <-Lorden
                
                tabpseu<-as.data.frame(Umax%*%matrix(1:ncol(Umax)))
                colnames(tabpseu)<-"CDBiplot class."
                if(exists("clase"))
                {
                        clase<-data.frame(clase)
                        colnames(clase)<-"Real Class"
                        tabaux<-as.data.frame(c(clase, tabpseu))
                        tabpseu<-table(tabaux)       
                }
                resultados<-list(U=Umax,V=Vmax, B=Borden,A=Amax,Ac=Acmax,Fk=Fmax,expvar=varexpmax,Lmax=Lmax, cor=cor(Amax), classif=tabpseu, itermax=itermax)
                print(resultados)
                cat("File saved in:    ",file="Results.txt")
                cat(getwd(),file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")		
                cat("Row clasification:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                write.table(Umax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Variable clasification:\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                write.table(Vmax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Eigenvalues:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Lmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Explained variability by each component:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(varexpmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")	
                cat("Row coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Amax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                cat("Centroid coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Acmax,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Variable coordinates:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(round(Borden,digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                        		
                file.append("Results.txt","temp.txt")	
                cat("Pseudoconfusion matrix:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(tabpseu,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                			
                file.append("Results.txt","temp.txt")	
                cat("Correlations matrix:\n",file="temp.txt")        				
                file.append("Results.txt","temp.txt")	
                write.table(round(cor(Amax),digits=3),file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum number of iterations:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(itermax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                cat("\n",file="temp.txt")                                	
                file.append("Results.txt","temp.txt")	
                cat("Maximum F:\n",file="temp.txt")					
                file.append("Results.txt","temp.txt")					
                write.table(Fmax,file="temp.txt", sep="\t",dec=",")
                file.append("Results.txt","temp.txt")
                
                file.show("Results.txt")
                file.remove("temp.txt")
                
                return(resultados)
        }
        
        cVal <- NULL
        cvVal <- NULL
        indicador <- NULL
        indicadorant <- NULL
        cbVal <- NULL
        resultados <- NULL
        tit_graph <- NULL
        Umax <- NULL
        Amax <- NULL
        Acmax <- NULL
        Borden <- NULL
        varexpmax <- NULL
        datos <- NULL
        etiquetas <- NULL
        textos <- NULL
        indexClosest <- NULL
        indexLabeledaux <- NULL
        anteriorx <- NULL
        anteriory <- NULL
        parPlotSize <- NULL
        usrCoords <- NULL
        xCoords <- NULL
        yCoords <- NULL
        
        wtipo <-tktoplevel()
        tkwm.title(wtipo,"Clustering and/or Disjoint Biplot")
        
        fontHeading <- tkfont.create(family="times",size=20,weight="bold",slant="italic")
        fontFixedWidth <- tkfont.create(family="courier",size=12)
        
        frametb1<-tkframe(wtipo, relief = "ridge", borderwidth = 2, background = "white")
        frametb2<-tkframe(wtipo, relief = "ridge", borderwidth = 2)#, background = "white")
        
        cbiplot <- tkradiobutton(frametb1)
        dbiplot <- tkradiobutton(frametb1)
        cdbiplot <- tkradiobutton(frametb1)
        rbValue <- tclVar("CDBiplot")
        tkconfigure(cbiplot,variable=rbValue,value="CBiplot")
        tkconfigure(dbiplot,variable=rbValue,value="DBiplot")
        tkconfigure(cdbiplot,variable=rbValue,value="CDBiplot")
        tkpack(tklabel(frametb1, text="Clustering Biplot (CBiplot)"), cbiplot,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(frametb1, text="Disjoint Biplot (DBiplot)"), dbiplot,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        tkpack(tklabel(frametb1, text="Clustering Disjoint Biplot (CDBiplot)"), cdbiplot,
               expand = "FALSE", side="top",expand="TRUE", fill = "both")
        
        
        OnOKtipo <- function()
        {
                tkdestroy(wtipo)
                P<-0
                Q<-0
                iter<-500
                times<-1000
                tol<-0.00001
                showgr <-"Y"
                cVal <<- as.character(tclvalue(rbValue))
                cvVal<<-"1"
                indicador<<-"0"
                indicadorant<<-"1"
                cbVal<<-"1"
                centro<-c(0,0)
                cchVal<-"1"
                fill<-"0"
                clabVal<-"1"
                indexlabeled<-0
                colores<-NULL
                
                winfor <- tktoplevel()
                
                if(cVal=="CBiplot")
                {
                        tkwm.title(winfor,"Clustering Biplot")
                }else{
                        if(cVal=="DBiplot")
                        {
                                tkwm.title(winfor,"Disjoint Biplot")
                        }else{
                                tkwm.title(winfor,"Clustering Disjoint Biplot")
                        }
                }
                
                framet1<-tkframe(winfor, relief = "ridge", borderwidth = 2, background = "white")
                framet2<-tkframe(winfor, relief = "ridge", borderwidth = 2, background = "white")
                framet11<-tkframe(framet1, relief = "ridge", borderwidth = 2, background = "white")
                framet12<-tkframe(framet1, relief = "ridge", borderwidth = 2, background = "white")
                
                
                OnOK <- function()
                {
                        tkdestroy(winfor)
                        if(cVal=="CBiplot")
                        {
                                P <<- as.numeric(tclvalue(pnames))
                                Q <<- as.numeric(tclvalue(qnames))
                        }else{
                                if(cVal=="DBiplot")
                                {
                                        Q <<- as.numeric(tclvalue(qnames))
                                }else{
                                        P <<- as.numeric(tclvalue(pnames))
                                        Q <<- as.numeric(tclvalue(qnames))
                                }
                        }
                        
                        
                        iter <<- as.numeric(tclvalue(iternames))
                        times <<- as.numeric(tclvalue(timesnames))
                        tol <<- as.numeric(tclvalue(tolnames))
                        
                        if(cVal=="CBiplot")
                        {
                                resultados <<- cbiplotint(data, P, Q, tol, iter, times, clase, showgr)
                                tit_graph <<- "CBiplot"
                        }else{
                                if(cVal=="DBiplot")
                                {
                                        resultados <<- dbiplotint(data, Q, tol, iter, times, showgr)
                                        tit_graph <<- "DBiplot"
                                }else{
                                        resultados <<- cdbiplotint(data, P, Q, tol, iter, times, clase, showgr)
                                        tit_graph <<- "CDBiplot"
                                }
                        }
                        
                        Umax <<-resultados$U
                        Amax <<-resultados$A
                        Acmax <<-resultados$Ac
                        Borden <<-resultados$B
                        varexpmax <<-resultados$expvar
                        
                        if (showgr=="Y")
                        {
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
                                
                                
                                showaxes <- function()
                                {
                                        if(cbVal=="1")
                                        {
                                                cbVal<<-"0"
                                                indicador<<-"0"
                                                indicadorant<-"0"
                                                tkrreplot(img)
                                        }else{
                                                cbVal<<-"1"
                                                indicador<<-"0"
                                                indicadorant<-"0"
                                                tkrreplot(img)
                                        }#end if(cbVal=="1")        
                                }#end showaxes
                                
                                
                                showvar <- function()
                                {
                                        if(cvVal=="1")
                                        {
                                                cvVal<<-"0"
                                                datos<<-rbind(Amax)
                                                colores<<-c()
                                                etiquetas<<-c()
                                                textos<<-c()
                                                
                                                if(clabVal=="1")
                                                {
                                                        textos<<-Amax
                                                        etiquetas<<-rownames(data)
                                                        colores<<-etiquetas
                                                        
                                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                                        {
                                                                for(i in 1: dim(Umax)[2])
                                                                {
                                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                                }
                                                        }else{
                                                                colores<<-rep(3, dim(Amax)[1])
                                                        }
                                                }
                                                
                                                
                                        }else{
                                                cvVal<<-"1" 
                                                
                                                datos<<-rbind(Amax, Borden)
                                                colores<<-c()
                                                etiquetas<<-c()
                                                textos<<-c()
                                                
                                                if(clabVal=="1")
                                                {
                                                        textos<<-rbind(Amax, Borden)
                                                        etiquetas<<-c(rownames(data), colnames(data))
                                                        colores<<-rownames(data)
                                                        
                                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                                        {
                                                                for(i in 1: dim(Umax)[2])
                                                                {
                                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                                }
                                                        }else{
                                                                colores<<-rep(3, dim(Amax)[1])
                                                        }
                                                        colores<<-c(colores,rep(1, dim(Borden)[1]))
                                                }
                                                
                                                
                                                if(clabVal=="0")
                                                {
                                                        textos<<-Borden
                                                        etiquetas<<-colnames(data)
                                                        colores<<-rep(1, dim(Borden)[1])
                                                }
                                                
                                        }#end if(cvVal=="1")   
                                        if(length(colores)>0){
                                                indexlabeled <<-1:length(colores)
                                        }else{
                                                indexlabeled<<-NULL
                                        }
                                        
                                        
                                        tkrreplot(img)
                                        
                                }#end showvar
                                
                                showlab <- function()
                                {
                                        if(clabVal=="0")
                                        {
                                                clabVal<<-"1"
                                                colores<<-c()
                                                etiquetas<<-c()
                                                textos<<-c()
                                                
                                                if(cvVal=="1")
                                                {
                                                        textos<<-rbind(Amax, Borden)
                                                        etiquetas<<-c(rownames(data), colnames(data))
                                                        colores<<-rownames(data)
                                                        
                                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                                        {
                                                                for(i in 1: dim(Umax)[2])
                                                                {
                                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                                }
                                                        }else{
                                                                colores<<-rep(3, dim(Amax)[1])
                                                        }
                                                        colores<<-c(colores,rep(1, dim(Borden)[1]))
                                                }
                                                
                                                if(cvVal=="0")
                                                {
                                                        textos<<-Amax
                                                        etiquetas<<-rownames(data)
                                                        colores<<-etiquetas
                                                        
                                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                                        {
                                                                for(i in 1: dim(Umax)[2])
                                                                {
                                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                                }
                                                        }else{
                                                                colores<<-rep(3, dim(Amax)[1])
                                                        }
                                                }
                                                
                                        }else{
                                                clabVal<<-"0"
                                                if(cvVal=="1")
                                                {
                                                        textos<<-Borden
                                                        etiquetas<<-colnames(data)
                                                        colores<<-rep(1, dim(Borden)[1])
                                                }else{
                                                        textos<<-c()
                                                        etiquetas<<-c()
                                                        colores<<-c()   
                                                }
                                        }#end if(clabVal=="1")  
                                        if(length(colores)>0){
                                                indexlabeled <<-1:length(colores)
                                        }else{
                                                indexlabeled<<-NULL
                                        } 
                                        
                                        tkrreplot(img)
                                        
                                }#end showlab
                                
                                
                                convexhull <- function()
                                {
                                        wch <-tktoplevel()
                                        tkwm.title(wch,"Convex hull")
                                        
                                        framech<-tkframe(wch, relief = "ridge", borderwidth = 2, background = "white")
                                        framech1<-tkframe(framech, relief = "ridge", borderwidth = 2, background = "white")
                                        framech2<-tkframe(framech, relief = "ridge", borderwidth = 2, background = "white")
                                        framech3<-tkframe(wch, relief = "ridge", borderwidth = 2, background = "white")
                                        framech21<-tkframe(framech2, relief = "ridge", borderwidth = 2, background = "white")
                                        framech22<-tkframe(framech2, relief = "ridge", borderwidth = 2, background = "white")
                                        
                                        chyes <- tkradiobutton(framech21)
                                        chno <- tkradiobutton(framech21)
                                        chfill <- tkradiobutton(framech22)
                                        chempt <- tkradiobutton(framech22)
                                        rbchValue <- tclVar("Yes")
                                        rbfeValue <- tclVar("Empty")
                                        tkconfigure(chyes,variable=rbchValue,value="Yes")
                                        tkconfigure(chno,variable=rbchValue,value="No")
                                        
                                        tkconfigure(chfill,variable=rbfeValue,value="Filled")
                                        tkconfigure(chempt,variable=rbfeValue,value="Empty")
                                        
                                        tkpack(tklabel(framech21, text="Yes"), chyes,
                                               expand = "FALSE", side="left",expand="TRUE", fill = "both")
                                        tkpack(tklabel(framech21, text="No"), chno,
                                               expand = "FALSE", side="left",expand="TRUE", fill = "both")
                                        
                                        tkpack(tklabel(framech22, text="Filled"), chfill,
                                               expand = "FALSE", side="left",expand="TRUE", fill = "both")
                                        tkpack(tklabel(framech22, text="Empty"), chempt,
                                               expand = "FALSE", side="left",expand="TRUE", fill = "both")
                                        
                                        tkpack(tklabel(framech1, text="CONVEX HULL"), side="top",expand="TRUE", fill = "both")
                                        
                                        onokch<-function()
                                        {
                                                tkdestroy(wch)
                                                cchVal <<- as.character(tclvalue(rbchValue))
                                                fill <<- as.character(tclvalue(rbfeValue))
                                                tkrreplot(img)
                                        }
                                        OK.butch<-tkbutton(framech3,text="   OK   ", command=onokch,  bg= "lightblue", width=20, foreground = "navyblue")
                                        tkbind(OK.butch, "<Return>",onokch)
                                        
                                        tkpack(framech21, framech22, side="top",expand="TRUE", fill = "both")
                                        tkpack(framech1, framech2, side="left",expand="TRUE", fill = "both")
                                        tkpack(OK.butch, side="top",expand="TRUE", fill = "both")
                                        tkpack(framech, framech3, side="top",expand="TRUE", fill = "both")
                                        
                                        
                                }#end convexhull
                                
                                OnLeftClick.up <- function(x,y)
                                {
                                        msg <- ("-To change the label press Yes.\n-To remove it press No.\n-If you do not want to do anything press Cancel.")
                                        mbval<- tkmessageBox(title="Change of label", message=msg,type="yesnocancel",icon="question")
                                        if (tclvalue(mbval)=="yes"){  
                                                indexlabeled <<- c(indexlabeled,indexClosest)
                                        }#end if (tclvalue(mbval)=="yes")
                                        
                                        if(tclvalue(mbval)=="no"){
                                                indexLabeledaux<<-c()
                                                for (i in (1:length(indexlabeled)))
                                                {
                                                        if (indexlabeled[i]!=indexClosest)
                                                                indexLabeledaux <<- c(indexLabeledaux,indexlabeled[i])
                                                }#end for (i in (1:length(indexlabeled)))
                                                indexlabeled<<-indexLabeledaux 
                                        }#end if(tclvalue(mbval)=="no")
                                        
                                        if(tclvalue(mbval)=="cancel"){
                                                textos[indexClosest,dim1] <<- anteriorx
                                                textos[indexClosest,dim2] <<- anteriory
                                        }#end if(tclvalue(mbval)=="cancel")
                                        tkrreplot(img)
                                }#end OnLeftClick.up <- function(x,y)
                                
                                OnLeftClick.move <- function(x,y)
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
                                        
                                        
                                        textos[indexClosest,dim1]<<-xPlotCoord
                                        textos[indexClosest,dim2]<<-yPlotCoord
                                        tkrreplot(img) 
                                }#end OnLeftClick.move <- function(x,y)
                                
                                OnLeftClick.down <- function(x,y)
                                {
                                        anteriorx <- NULL
                                        anteriory <- NULL
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
                                        
                                        squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
                                        indexClosest <<- which.min(squared.Distance) 
                                        
                                        anteriorx <<- textos[indexClosest,dim1]
                                        anteriory <<- textos[indexClosest,dim2]          
                                }#end OnLeftClick.down <- function(x,y)
                                
                                
                                dim1 <- 1
                                dim2 <- 2
                                dim3 <- 3
                                
                                
                                datos<<-rbind(Amax, Borden)
                                colores<<-c()
                                etiquetas<<-c()
                                textos<<-c()
                                
                                if(clabVal=="1" & cvVal=="1")
                                {
                                        textos<<-rbind(Amax, Borden)
                                        etiquetas<<-c(rownames(data), colnames(data))
                                        colores<<-rownames(data)
                                        
                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                        {
                                                for(i in 1: dim(Umax)[2])
                                                {
                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                }
                                        }else{
                                                colores<<-rep(3, dim(Amax)[1])
                                        }
                                        colores<<-c(colores,rep(1, dim(Borden)[1]))
                                }
                                
                                if(clabVal=="1" & cvVal=="0")
                                {
                                        textos<<-Amax
                                        etiquetas<<-rownames(data)
                                        colores<<-etiquetas
                                        
                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                        {
                                                for(i in 1: dim(Umax)[2])
                                                {
                                                        colores[which(Umax[,i]==1)]<<-i+1       
                                                }
                                        }else{
                                                colores<<-rep(3, dim(Amax)[1])
                                        }
                                }
                                
                                if(clabVal=="0" & cvVal=="1")
                                {
                                        textos<<-Borden
                                        etiquetas<<-colnames(data)
                                        colores<<-rep(1, dim(Borden)[1])
                                }
                                
                                if(length(colores)>0){
                                        indexlabeled <<-1:length(colores)
                                }else{
                                        indexlabeled<<-NULL
                                }
                                
                                
                                plotFunction<-function(screen=TRUE)
                                {   
                                        xCoords<<-textos[,dim1]
                                        yCoords<<-textos[,dim2]
                                        
                                        plot(datos[,c(dim1, dim2)], type="n", main=tit_graph, xlab=paste("Dim ",dim1, "(",varexpmax[dim1,2],"%)"), ylab=paste("Dim ",dim2, "(",varexpmax[dim2,2],"%)"))
                                        
                                        
                                        if(cbVal=="1")
                                        {
                                                abline(h=centro[2],v=centro[1],lty="dotted")        
                                        }
                                        
                                        if(cVal=="CBiplot" | cVal=="CDBiplot")
                                        {
                                                
                                                for (p in 1:dim(Umax)[2])
                                                {
                                                        if(cchVal=="Yes")
                                                        {
                                                                clusteri <-Amax[which(Umax[,p]==1),c(dim1,dim2), drop=FALSE]
                                                                clusteri <- t(t(clusteri))
                                                                hpts <- chull(clusteri)
                                                                hpts <- c(hpts, hpts[1]) 
                                                                polygon(clusteri[hpts,],border=p+1)
                                                                if(fill=="Filled")
                                                                {
                                                                        polygon(clusteri[hpts,],col=p+1, border=p+1)
                                                                }else{
                                                                        points(Amax[which(Umax[,p]==1),dim1], Amax[which(Umax[,p]==1),dim2], col=p+1, pch=8)
                                                                        polygon(clusteri[hpts,],border=p+1)
                                                                }
                                                        }else{
                                                                points(Amax[which(Umax[,p]==1),dim1], Amax[which(Umax[,p]==1),dim2], col=p+1, pch=8)
                                                        }
                                                        
                                                        points(Acmax[p,dim1], Acmax[p,dim2], col=p+1, pch=18)
                                                }
                                                
                                                
                                                
                                        }else{
                                                points(Amax[,dim1], Amax[,dim2],col=3, pch=8)
                                        }
                                        
                                        
                                        
                                        if(cvVal=="1")
                                        {
                                                suppressWarnings( arrows(centro[1],centro[2],Borden[,dim1],Borden[,dim2],lty=1, length=0.08))
                                        }
                                        
                                        if (length(indexlabeled)>0)
                                                for (i in (1:length(indexlabeled)))
                                                {
                                                        indexClosest <- indexlabeled[i]
                                                        text(textos[indexClosest,dim1],textos[indexClosest,dim2], labels=etiquetas[indexClosest], col=colores[indexClosest])
                                                }#end for (i in (1:length(indexLabeled)))
                                        
                                        parPlotSize <<- par("plt")
                                        usrCoords   <<- par("usr")
                                        
                                }
                                
                                
                                wgr <- tktoplevel()
                                tkwm.title(wgr,tit_graph)
                                
                                
                                g3d<-function()
                                {
                                        if (Q>2)
                                        { 
                                                bg3d("white")
                                                aspect3d("iso")
                                                lims <- par3d("bbox")
                                                
                                                if(cbVal=="1")
                                                {
                                                        axes3d()
                                                }
                                                
                                                if(cVal!="DBiplot")
                                                {
                                                        
                                                        for (p in 1:dim(Umax)[2])
                                                        {
                                                                points3d(Amax[which(Umax[,p]==1),dim1], Amax[which(Umax[,p]==1),dim2],Amax[which(Umax[,p]==1),dim3], col=p+1)
                                                        }
                                                        points3d(Acmax[p,dim1], Acmax[p,dim2], Acmax[p,dim3], col=p+1, pch=18)
                                                        
                                                }else{
                                                        points3d(Amax[,dim1], Amax[,dim2],Amax[,dim3], col=3)
                                                }
                                                
                                                if(cvVal=="1")
                                                {
                                                        for (i in 1:(dim(Borden)[1]))
                                                        {
                                                                linea<-rbind(Borden[i,c(dim1, dim2, dim3)],c(0,0,0))	
                                                                segments3d(linea[,1],linea[,2], linea[,3],color=1)
                                                        }#end for (i in 1:(dim(Borden)[1]))
                                                        
                                                }
                                                
                                                if(length(colores)>0)
                                                        suppressWarnings(texts3d(textos[indexlabeled,dim1], textos[indexlabeled,dim2], textos[indexlabeled,dim3],etiquetas[indexlabeled],color=colores))
                                                
                                                rgl.bringtotop()
                                        }else{
                                                msg <- "You have selected less than 3 dimensions. 3D-graph not available"
                                                tkmessageBox(message=msg)
                                        }#end if (Q>2)
                                }#end g3d<-function()
                                
                                
                                topMenugr <- tkmenu(wgr)
                                tkconfigure(wgr, menu = topMenugr)
                                menuFile <- tkmenu(topMenugr, tearoff = FALSE)
                                menuSaveAs <- tkmenu(topMenugr, tearoff = FALSE)
                                menu3d <- tkmenu(topMenugr, tearoff = FALSE)
                                menuopt <-tkmenu(topMenugr, tearoff = FALSE)
                                menuclus <-tkmenu(topMenugr, tearoff = FALSE)
                                
                                tkadd(menuFile, "command", label = "Copy image", command = function() {tkrreplot(img)})
                                tkadd(menuFile, "cascade", label = "Save image", menu = menuSaveAs)
                                tkadd(menuSaveAs, "command", label = "Pdf file", command = function() {SaveFilePDF()})
                                tkadd(menuSaveAs, "command", label = "Eps file", command = function() {SaveFileeps()})
                                tkadd(menuSaveAs, "command", label = "Png file", command = function() {SaveFilePng()})
                                tkadd(menuSaveAs, "command", label = "Jpg/Jpeg file", command = function() {SaveFileJPG()})
                                tkadd(menuFile, "separator")
                                tkadd(menuFile, "command", label = "Exit", command = function() {tkdestroy(wgr)})
                                tkadd(menuclus, "command", label = "Convex-hull", command = function() {convexhull()})
                                
                                tkadd(topMenugr, "cascade", label ="File", menu = menuFile)
                                tkadd(menuFile, "separator")
                                tkadd(topMenugr, "cascade", label = "3D", menu = menu3d)
                                tkadd(topMenugr, "cascade", label = "Options", menu = menuopt)
                                if(cVal!="DBiplot")
                                        tkadd(topMenugr, "cascade", label = "Clusters", menu = menuclus)
                                
                                tkadd(menu3d, "command", label = "3D", command = function() {g3d()})
                                tkadd(menuopt, "command", label = "Change title", command = function() {changetit()})
                                tkadd(menuopt, "command", label = "Show/Hide axes", command = function() {showaxes()})
                                tkadd(menuopt, "command", label = "Show/Hide variables", command = function() {showvar()})
                                tkadd(menuopt, "command", label = "Show/Hide row labels", command = function() {showlab()})
                                
                                img <- tkrplot(wgr,fun=plotFunction,hscale=1.5,vscale=1.5)
                                framedim1<-tkframe(wgr, relief = "ridge", borderwidth = 2, background = "whitesmoke", height=40)
                                
                                comboBoxdim1 <- tkwidget(framedim1,"ComboBox",editable=FALSE,values=rep(1:Q),width=15, text= dim1)
                                comboBoxdim2 <- tkwidget(framedim1,"ComboBox",editable=FALSE,values=rep(1:Q),width=15, text= dim2)
                                
                                chang.symdim1 <- function()
                                {
                                        dim1 <<-as.numeric(tclvalue(tcl(comboBoxdim1,"getvalue")))+1
                                        dim2 <<-as.numeric(tclvalue(tcl(comboBoxdim2,"getvalue")))+1
                                        if (Q>2)
                                                dim3 <<-as.numeric(tclvalue(tcl(comboBoxdim3,"getvalue")))+1
                                        tkrreplot(img)
                                }#end chang.symdim1 <- function()
                                
                                Change.symboldim1 <-tkbutton(framedim1,text="Choose",command=chang.symdim1, bg= "lightblue", width=15, foreground = "navyblue")
                                
                                if (Q>2){
                                        comboBoxdim3 <- tkwidget(framedim1,"ComboBox",editable=FALSE,values=rep(1:Q),width=15, text= dim3)
                                        tkpack(tklabel(framedim1, text="Select X, Y and Z axes numbers:"), expand="FALSE", side= "top", fill ="both")
                                        tkpack(comboBoxdim1, comboBoxdim2, comboBoxdim3, Change.symboldim1, side="top", fill="x")
                                        tkpack(framedim1, side="left", expand="FALSE", fill="x")
                                        
                                }
                                tkpack(img, side="top", expand="TRUE", fill="both")
                                
                                tkbind(img, "<B1-Motion>",OnLeftClick.move)
                                tkbind(img, "<ButtonPress-1>",OnLeftClick.down)
                                tkbind(img, "<ButtonRelease-1>",OnLeftClick.up)
                                tkconfigure(img,cursor="pencil")
                                
                        }#end showgr
                        
                        
                }#end OnOK <- function()
                
                OK.but<-tkbutton(framet2,text="   OK   ", command=OnOK,  bg= "lightblue", width=20, foreground = "navyblue")
                tkbind(OK.but, "<Return>",OnOK)
                
                pnames<-tclVar(P)
                qnames<-tclVar(Q)
                tolnames<-tclVar(tol)
                iternames<-tclVar(iter)
                timesnames<-tclVar(times)
                
                if(cVal!="DBiplot")
                {
                        
                        entry.pnames <-tkentry(framet12,width="50",textvariable=pnames, bg="white")
                        tkbind(entry.pnames, "<Return>",OnOK)        
                }
                
                entry.qnames <-tkentry(framet12,width="50",textvariable=qnames, bg="white")
                tkbind(entry.qnames, "<Return>",OnOK)
                
                entry.tolnames <-tkentry(framet12,width="50",textvariable=tolnames, bg="white")
                tkbind(entry.tolnames, "<Return>",OnOK)
                
                entry.iternames <-tkentry(framet12,width="50",textvariable=iternames, bg="white")
                tkbind(entry.iternames, "<Return>",OnOK)
                
                entry.timesnames <-tkentry(framet12,width="50",textvariable=timesnames, bg="white")
                tkbind(entry.timesnames, "<Return>",OnOK)
                
                if(cVal!="DBiplot")
                {
                        tkpack(tklabel(framet11,text="Number of clusters:"),
                               tklabel(framet11,text="Number of components:"),
                               tklabel(framet11,text="Tolerance:"),
                               tklabel(framet11,text="Iterations:"),
                               tklabel(framet11,text="Repetitions of the algorithm:"),
                               expand = "TRUE", side="top", fill = "both")
                        tkpack(OK.but)
                        tkpack(entry.pnames , entry.qnames, entry.tolnames, entry.iternames,
                               entry.timesnames, 
                               expand = "TRUE",side="top", fill="both")
                }else{
                        tkpack( tklabel(framet11,text="Number of components:"),
                                tklabel(framet11,text="Tolerance:"),
                                tklabel(framet11,text="Iterations:"),
                                tklabel(framet11,text="Repetitions of the algorithm:"),
                                expand = "TRUE", side="top", fill = "both")
                        tkpack(OK.but)
                        tkpack(entry.qnames, entry.tolnames, entry.iternames,
                               entry.timesnames, 
                               expand = "TRUE",side="top", fill="both")      
                }
                
                tkpack(framet11, framet12, expand = "TRUE",side="left", fill="y")
                tkpack(framet1, framet2, expand = "TRUE",side="top", fill="y")  
                tkfocus(winfor)     
                
        }#end OnOKtipo
        
        OK.butipo <-tkbutton(frametb2,text="   OK   ",command=OnOKtipo, bg= "lightblue", width=20, foreground = "navyblue")
        tkbind(OK.butipo, "<Return>",OnOKtipo)
        tkpack(OK.butipo)
        
        tkpack(tklabel(wtipo, text="  CLUSTERING AND/OR DISJOINT BIPLOT  ",font=fontHeading, foreground = "blue"),frametb1, frametb2, expand = "TRUE", side="top", fill="both")
}
