# gsvdScal: generalized singular value decomposition with singular vector
# length scaled to generalized singular value (in four ways)

gsvdScal<-function(gs,scaling)
{
	gd<-gs$gd
  gu<-gs$gu
  gv<-gs$gv
	dl<-gs$dl
  dx<-gs$dx
  dz<-gs$dz
	dxs<-dx
  dzs<-dz
  gus<-gu
  gvs<-gv
	ind<-gs$ind
  neval<-length(ind)

  #--------- Benzecri ------------
  if (scaling[1]=="Benzecri") {
		for (i in length(gd)) {
			gus[,i]<-gu[,i]*gd[i]

			}
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-outer(gu[,j],dl[i,])+gd[j]*dx[,i,]
			}
	}
  if (scaling[2]=="Benzecri") {
		for (i in length(gd)) gvs[,i]<-gv[,i]*gd[i]
		for (i in 1:neval) {
			j<-ind[i]
			dzs[,i,]<-outer(gv[,j],dl[i,])+gd[j]*dz[,i,]
			}
  }
  #------ end Benzecri ----------

  #--------- Goodman ------------
  if (scaling[1] == "Goodman") {
		for (i in length(gd)) gus[,i]<-gu[,i]*sqrt(gd[i])
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-1/(2*sqrt(gd[j]))*outer(gu[,j],dl[i,])+sqrt(gd[j])*dx[,i,]
			}
  }
	if (scaling[2] == "Goodman") {
		for (i in length(gd)) gvs[,i]<-gv[,i]*sqrt(gd[i])
		for (i in 1:neval) {
			j<-ind[i]
			dzs[,i,]<-1/(2*sqrt(gd[j]))*outer(gv[,j],dl[i,])+sqrt(gd[j])*dz[,i,]
    }
	}
  #---------- end Goodman ---------

  #------------ centroids ----------
  if (scaling[2] == "centroid") {
		for (i in length(gd)) gvs[,i]<-gv[,i]*gd[i]
		for (i in 1:neval) {
			j<-ind[i]
			dzs[,i,]<-outer(gv[,j],dl[i,])+gd[j]*dz[,i,]
			}
  }
	if (scaling[1] == "centroid") {
		for (i in length(gd)) gus[,i] <- gu[,i]*gd[i]
		for (i in 1:neval) {
			j<-ind[i]
			dxs[,i,]<-outer(gu[,j],dl[i,])+gd[j]*dx[,i,]
		}
	}
	#------------ end centroid -------------
	
	#--------------- standard --------------
	if (scaling[1] == "standard") {
    gus <- gu
    dxs <- dx
  }
  if (scaling[2] == "standard") {
    gvs <- gv
    dzs <- dz
  }
  #-------------- end standard -----------
	
return(list(gd=gd,gu=gus,gv=gvs,dl=dl,dx=dxs,dz=dzs,ind=ind))
}