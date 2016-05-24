plotvolu2d<-function(vd,theta=NULL,phi=NULL,typer="flat")
{
# typer "dens"/"flat"

if (is.null(phi)) phi<-30

if (vd$type2=="slice"){

if (vd$type=="radius"){

 if (typer=="flat"){
     if (is.null(theta)) theta<-50
     persp(vd$x,vd$y,vd$z,
     xlab="level",ylab="",zlab="radius",ticktype="detailed",
     phi=phi,theta=theta)
 }
 else{
     levnum<-length(vd$x)
     ynumold<-length(vd$y)
     maksi<-max(vd$z)
     gnum<-100
     step<-maksi/(gnum-1)
     xnew<-seq(0,maksi,step)
     znew<-matrix(0,gnum,ynumold)
     for (i in 1:levnum){
        for (j in 1:ynumold){
           highness<-round(gnum*vd$z[i,j]/maksi)
           znew[1:highness,j]<-vd$x[i]  #level[i]
        }
     }
     if (is.null(theta)) theta<-40
     persp(xnew,vd$y,znew,
     xlab="radius",ylab="",zlab="level",
     ticktype="detailed",
     phi=phi,theta=theta)
 }
}

if (vd$type=="proba"){
if (is.null(theta)) theta<--130
if (vd$norma) xlab<-"normalized volume" else xlab<-"volume"
persp(vd$x,vd$y,vd$z,
xlab=xlab,ylab="",zlab="radius",ticktype="detailed",
phi=phi,theta=theta)
}

}   #type2=="slice"

else{ #type2=="boundary"

if (is.null(theta)) theta<-50
persp(vd$x,vd$y,vd$z,
xlab="",ylab="",zlab="level",ticktype="detailed",
phi=phi,theta=theta)

}


}




