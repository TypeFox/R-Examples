lifetable<-function(age=NULL,numbers=NULL,r=NULL, type=1){

lt<-as.data.frame(cbind(age,numbers))
tx<-length(age)
dx<-function(nx,tx){
       for(t in 1:tx){
           if(t<tx) lt$dx[t]<<-nx[t]-nx[t+1]
           if(t==tx) lt$dx[t]<<-NA
         }
    }
qx<-function(dx,nx){
   lt$qx<<-ifelse(is.nan(dx/nx),NA,dx/nx)
   lt$px<<-1-lt$qx
  }
Lx<-function(nx,tx){
   for(t in 1:tx) {
    if(t<tx) lt$Lx[t]<<-(nx[t]+nx[t+1])/2
    if(t==tx) lt$Lx[t]<<-NA
   }
 }

Tex<-function(Lx,nx,tx){
for(t in 1:tx){
   if(t<tx){
      lt$Tx[t]<<-sum(Lx[t:as.numeric(tx-1)],na.rm=T)
      lt$ex[t]<<-lt$Tx[t]/nx[t]
    }
   if(t==tx) {lt$Tx[t]<<-NA;lt$ex[t]<<-NA} 
  }
}
if(type==1){
	names(lt)[2]<-c("dx")
	n0<-sum(lt$dx,na.rm=T)
	lt$nx[1]<-n0
	for(t in 2:tx){lt$nx[t]<-lt$nx[t-1]-lt$dx[t-1]}
      lt$lx<-lt$nx/n0
	qx(lt$dx,lt$nx)
	Lx(lt$nx,tx)
	Tex(lt$Lx,lt$nx,tx)
  }
if(type==2){
  	names(lt)[2]<-c("nx")
      dx(lt$nx,tx)
      qx(lt$dx,lt$nx)
	lt$lx<-lt$nx/lt$nx[1]
      Lx(lt$nx,tx)
	Tex(lt$Lx,lt$nx,tx)
  }
if(type==3){
	names(lt)[2]<-c("dx")
      lt$dx<-lt$dx*exp(r*lt$age)
	n0<-sum(lt$dx,na.rm=T)
	lt$nx[1]<-n0
	for(t in 2:tx){lt$nx[t]<-lt$nx[t-1]-lt$dx[t-1]}
      lt$lx<-lt$nx/n0
	qx(lt$dx,lt$nx)
	Lx(lt$nx,tx)
	Tex(lt$Lx,lt$nx,tx)
  }
if(type==4){
      names(lt)[2]<-c("nx")
      lt$nx<-lt$nx*exp(r*lt$age)
      dx(lt$nx,tx)
      if(any(lt$dx<0)) stop("some dx values are negative - have to smooth new nx*e(rx) values outside of program")
      qx(lt$dx,lt$nx)
	lt$lx<-lt$nx/lt$nx[1]
      Lx(lt$nx,tx)
	Tex(lt$Lx,lt$nx,tx)
   }
 return(lt)
}


