`PsiTildeCalc` <-
function(u,H,test,omega,vm,vminv,p,K){
        PsiHatSum<-array(0,c(p,p))
        r<-dim(test)[[1]] 
        w<-matrix(NA,K,r)
        PsiTilde<-array(NA,c(r,K,p,p))
        HU<-H*u
        for (i in 1:K){
            PsiHatSum<-PsiHatSum + matrix(HU[i,],ncol=1) %*% matrix(HU[i,],nrow=1)
           w[i,]<- diag( test %*% ( solve( vminv - omega[i,,]) -vm ) %*% t(test) )
        }
        wstandardized<-t( t(w)/apply(w,2,sum))
        for (i in 1:K){
            for (j in 1:r){
                PsiTilde[j,i,,]<-w[i,j]*PsiHatSum
            }
        }
        PsiTilde
    }

