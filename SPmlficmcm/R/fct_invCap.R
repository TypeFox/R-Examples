fct_invCap <-
function(dat,Ni,vrzc,vrze,genm,genc,theta.start){
            #dat : donnee
            #vrzc variable outcome; vrze : variable explicative, genm : gnotype de la mere genc : genotype de lenfant
            n<-dim(dat)[1];n1<-dim(dat[dat[vrzc]==1,])[1];n0<-n-n1;ni=c(n0,n1);di=Ni-ni;N=sum(Ni)
            r1=(di[2]+n1)/N;r0=1-r1
            # nv
            ind<-rep(1,dim(dat)[1]);dat1<-data.frame(dat,ind)
            Nijmc<-fpol1(dat1,c(vrzc,vrze,genm,genc),"ind","Cjm")
            vep<-unique(Nijmc[,genm])
            dgm<-Prgm_HW1(as.numeric(unique(dat[,genm])),theta.start)
            # calcul de Pgm|y
            da1<-dat[dat[,vrzc]==1,];da0<-dat[dat[,vrzc]==0,];m0<-dim(da0)[1];m1<-dim(da1)[1]
            if(length(vep)==3){
                               nm<-c(sum(Nijmc[Nijmc[,genm]==0,][,"Cjm"]),sum(Nijmc[Nijmc[,genm]==1,][,"Cjm"]),sum(Nijmc[Nijmc[,genm]==2,][,"Cjm"])) 
                                # construction de l estimateur de Nm
                               pgi1<-di[2]*(dim(da1[da1[,genm]==0,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==0,])[1]/m0)
                               pgi2<-di[2]*(dim(da1[da1[,genm]==1,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==1,])[1]/m0)
                               pgi3<-di[2]*(dim(da1[da1[,genm]==2,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==2,])[1]/m0)
                               pg<-c(pgi1,pgi2,pgi3)
                               # les diff grand du deno
                               Nem_sp<-nm+(N-n)
                               # sup
                               Nem20<-(nm+pg)*r0
                               Nem21<-(nm+pg)*r1
             
                               #borne sup
                               brs=rbind(1-(di[1]*dgm)/(Nem_sp*r0),1-(di[2]*dgm)/(Nem_sp*r1))
                               #born inf
                               brf=rbind(1-(di[1]*dgm)/(nm*r0),1-(di[2]*dgm)/(nm*r1))
                               # tt danns la borne
                               ma.u=rbind(c(runif(1,brf[1,1],brs[1,1]),runif(1,brf[1,2],brs[1,2]),runif(1,brf[1,3],brs[1,3])),
                               c(runif(1,brf[2,1],brs[2,1]),runif(1,brf[2,2],brs[2,2]),runif(1,brf[2,3],brs[2,3])))
                               # estimateur empirique (avec correction)
                               mat.emp=rbind(1-(di[1]*dgm)/(Nem20),1-(di[2]*dgm)/(Nem21))   
                                 }
            if(length(vep)==2 & sum(vep)==1){
                                               nm<-c(sum(Nijmc[Nijmc[,genm]==0,][,"Cjm"]),sum(Nijmc[Nijmc[,genm]==1,][,"Cjm"]))  
                                # construction de lestimateur de Nm^*
                               pgi1<-di[2]*(dim(da1[da1[,genm]==0,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==0,])[1]/m0)
                               pgi2<-di[2]*(dim(da1[da1[,genm]==1,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==1,])[1]/m0)
                               pg<-c(pgi1,pgi2)
                               # les diff grand du deno
                               Nem_sp<-nm+(N-n)
                               # sup
                               Nem20<-(nm+pg)*r0
                               Nem21<-(nm+pg)*r1
             
                               #borne sup
                               brs=rbind(1-(di[1]*dgm)/(Nem_sp*r0),1-(di[2]*dgm)/(Nem_sp*r1))
                               #born inf
                               brf=rbind(1-(di[1]*dgm)/(nm*r0),1-(di[2]*dgm)/(nm*r1))
                               # tt danns la borne
                               ma.u=rbind(c(runif(1,brf[1,1],brs[1,1]),runif(1,brf[1,2],brs[1,2])),
                               c(runif(1,brf[2,1],brs[2,1]),runif(1,brf[2,2],brs[2,2])))
                               # estimateur empirique (avec correction)
                               mat.emp=rbind(1-(di[1]*dgm)/(Nem20),1-(di[2]*dgm)/(Nem21)) 
                                               }
            if(length(vep)==2 & sum(vep)==2){
                                               nm<-c(sum(Nijmc[Nijmc[,genm]==0,][,"Cjm"]),sum(Nijmc[Nijmc[,genm]==2,][,"Cjm"]))  
                                # construction de lestimateur de Nm^*
                               pgi1<-di[2]*(dim(da1[da1[,genm]==0,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==0,])[1]/m0)
                               pgi2<-di[2]*(dim(da1[da1[,genm]==2,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==2,])[1]/m0)
                               pg<-c(pgi1,pgi2)
                               # les diff grand du deno
                               Nem_sp<-nm+(N-n)
                               # sup
                               Nem20<-(nm+pg)*r0
                               Nem21<-(nm+pg)*r1
             
                               #borne sup
                               brs=rbind(1-(di[1]*dgm)/(Nem_sp*r0),1-(di[2]*dgm)/(Nem_sp*r1))
                               #born inf
                               brf=rbind(1-(di[1]*dgm)/(nm*r0),1-(di[2]*dgm)/(nm*r1))
                               # tt danns la borne
                               ma.u=rbind(c(runif(1,brf[1,1],brs[1,1]),runif(1,brf[1,2],brs[1,2])),
                               c(runif(1,brf[2,1],brs[2,1]),runif(1,brf[2,2],brs[2,2])))
                               # estimateur empirique (avec correction)
                               mat.emp=rbind(1-(di[1]*dgm)/(Nem20),1-(di[2]*dgm)/(Nem21)) 
                                               }
            if(length(vep)==2 & sum(vep)==3){
                                               nm<-c(sum(Nijmc[Nijmc[,genm]==1,][,"Cjm"]),sum(Nijmc[Nijmc[,genm]==2,][,"Cjm"])) 
                                # construction de lestimateur de Nm
                               pgi1<-di[2]*(dim(da1[da1[,genm]==1,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==1,])[1]/m0)
                               pgi2<-di[2]*(dim(da1[da1[,genm]==2,])[1]/m1)+di[1]*(dim(da0[da0[,genm]==2,])[1]/m0)
                               pg<-c(pgi1,pgi2)
                               # les diff grand du deno
                               Nem_sp<-nm+(N-n)
                               # sup
                               Nem20<-(nm+pg)*r0
                               Nem21<-(nm+pg)*r1
             
                               #borne sup
                               brs=rbind(1-(di[1]*dgm)/(Nem_sp*r0),1-(di[2]*dgm)/(Nem_sp*r1))
                               #born inf
                               brf=rbind(1-(di[1]*dgm)/(nm*r0),1-(di[2]*dgm)/(nm*r1))
                               # tt danns la borne
                               ma.u=rbind(c(runif(1,brf[1,1],brs[1,1]),runif(1,brf[1,2],brs[1,2])),
                               c(runif(1,brf[2,1],brs[2,1]),runif(1,brf[2,2],brs[2,2])))
                               # estimateur empirique (avec correction)
                               mat.emp=rbind(1-(di[1]*dgm)/(Nem20),1-(di[2]*dgm)/(Nem21)) 
                                               }
            return(list(brf=brf,brs=brs,ma.u=ma.u,mat.emp=mat.emp))
            }
