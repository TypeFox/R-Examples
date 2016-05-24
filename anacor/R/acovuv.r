acovuv <- function(gs, prop, nr, nc, N, ndim)
{
# computes VC-matrices for row and column scores

gdl<-gs$dl
acovd<-gdl%*%(prop*t(gdl))/N
acovu<-array(0,c(ndim,ndim,nr))
gdx<-gs$dx

for (i in 1:nr) {
        gdxi <- rbind(gdx[i,,])
        
        #gdxi<-gdx[i,,]
	gdxs<-drop(gdxi%*%prop)
     
        acovu[,,i]<-(gdxi%*%(prop*t(gdxi))-outer(gdxs,gdxs))/N
	}

acovv<-array(0,c(ndim,ndim,nc))
gdz<-gs$dz

for (i in 1:nc) {
	gdzi<- rbind(gdz[i,,])
        
        #gdzi<-gdz[i,,]
        gdzs<-drop(gdzi%*%prop)
	acovv[,,i]<-(gdzi%*%(prop*t(gdzi))-outer(gdzs,gdzs))/N
	}
return(list(acovd=acovd,acovu=acovu,acovv=acovv))
}
