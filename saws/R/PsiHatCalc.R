`PsiHatCalc` <-
function(u,H){
        K<-dim(H)[[1]]
        p<-dim(H)[[2]]
        out<-array(NA,c(K,p,p))
        HU<-H*u
        for (i in 1:K){
            out[i,,]<-matrix(HU[i,],ncol=1) %*% matrix(HU[i,],nrow=1)
        }
        out
    }

