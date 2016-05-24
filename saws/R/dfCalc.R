`dfCalc` <-
function(Psi,u,omega,H,test,vm){
        K<-dim(u)[[1]]
        p<-dim(u)[[2]]
        ## July 9, 2010: fix bug when p=1
        ## Psi is treated as vector unless specifically make it an array
        ## if p>1, was already an array
        if (p==1){
            Psi<-array(Psi,c(K,p,p))
        }
        if (is.vector(test)){
            r<-1
            test<-matrix(test,1,length(test))
        } else r<-dim(test)[[1]]
        G<-Psimat<-matrix(0,p*K,p*K)
        for (i in 1:K){
       index<- ((i-1)*p + 1): (i*p)
       omegaivm<-omega[i,,] %*% vm
            G[index,index]<- diag(p) 
            G[index,]<- G[index,] - matrix(rep(omegaivm,K),p,p*K)
            Psimat[index,index]<-Psi[i,,]
        }
        df<-rep(NA,r)
        for (j in 1:r){
            M<-matrix(0,p*K,p*K)
            for (i in 1:K){
           index<- ((i-1)*p + 1): (i*p)
                M[index,index]<- diag(H[i,]) %*% vm %*% matrix(test[j,],p,1) %*% matrix(test[j,],1,p) %*% vm %*% diag(H[i,]) 
            }
            Psi.B1<-Psimat %*% t(G) %*% M %*% G
            df[j]<-sum(diag( Psi.B1 ))^2 / sum(as.vector(Psi.B1)*as.vector(t(Psi.B1)))
        }
        df
    }

