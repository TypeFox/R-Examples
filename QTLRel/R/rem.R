
# matrix associated with random term z*b where b~norm(0,s)
rem.1<- function(b,z=1){
   b<- as.matrix(b)[,1]
   nn<- length(b)
   if(length(z)==1){
      if(z!=1) stop("Check formula for errors.")
      z<- rep(1,nn)
   }else{
      z<- as.matrix(z)[,1]
      if(!is.numeric(z) || length(z)!=nn)
         stop("Check formula for errors.")
   }
   mtr<- matrix(0,nrow=nn,ncol=nn)
   buv<- unique(b)
   if(length(buv)<2)
      stop("No variation in the group variable?")
   for(i in 1:length(buv)){
      idx<- b==buv[i]
      mtr[idx,idx]<- z[idx]%o%z[idx]
   }

   mtr
}

# matrix associated with random term z*b where b~norm(0,s) and a/b
rem.2<- function(b,a,z=1){
   a<- as.matrix(a)[,1]
   b<- as.matrix(b)[,1]
   nn<- length(b)
   if(length(z)==1){
      if(z!=1) stop("Check formula for errors.")
      z<- rep(1,nn)
   }else{
      z<- as.matrix(z)[,1]
      if(!is.numeric(z) || length(z)!=nn)
         stop("Check formula for errors.")
   }
   mtr<- matrix(0,nrow=nn,ncol=nn)
   auv<- unique(a)
   if(length(auv)<2)
      stop("No variation in the group variable?")
   for(i in 1:length(auv)){
      idx<- a==auv[i]
      buv<- unique(b[idx])
      for(j in 1:length(buv)){
         idxTmp<- idx & (b==buv[j])
         mtr[idxTmp,idxTmp]<- z[idxTmp]%o%z[idxTmp]
      }
   }

   mtr
}

rem<- function(formula,data){
   mlst<- list()
   formula<- deparse(formula)
   if(length(grep("~",formula,))<1)
      stop("Check formula for errors.")
   tms<- strsplit(formula,split="~",fixed=TRUE)[[1]]
      tms<- tms[-1]
      tms<- strsplit(tms,split="+",fixed=TRUE)[[1]]
   ii<- 0
   for(i in 1:length(tms)){
      tt<- tms[i]
         tt<- strsplit(tt,"|",fixed=TRUE)[[1]]
      if(length(tt)<1 || length(tt)>2){
         stop("Check formula for errors.")
      }else if(length(tt)==1){
         zz<- 1
      }else if(length(tt)==2){
         zz<- formula(paste("~",tt[1]))
            zz<- attr(terms(zz),"term.labels")
            zz<- data[,zz]
            if(length(zz)<1) zz<- 1
         tt<- tt[2]
      }
      tt<- strsplit(tt,split="/")[[1]]
         tt<- formula(paste("~",paste(tt,collapse="+")))
            tt<- attr(terms(tt),"term.labels")
      for(j in 1:length(tt)){
         ii<- ii+1
         if(j==1){
            mlst[[ii]]<- rem.1(data[,tt[j]],zz)
         }else mlst[[ii]]<- rem.2(data[,tt[j]],data[,tt[j-1]],zz)
         names(mlst)[ii]<- tt[j]
      }
   }
   mlst
}

