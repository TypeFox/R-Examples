
M.empirical<-function(Linf=NULL,Winf=NULL,Kl=NULL,Kw=NULL,T=NULL,tmax=NULL,tm=NULL,GSI=NULL,
      Wdry=NULL,Wwet=NULL,Bl=NULL, method=c(1,2,3,4,5,6,7,8,9)){
   if(any(method==1) & any(is.null(Linf),is.null(Kl),is.null(T)))
           stop("Method 1 requires Linf, Kl, and T")
    if(any(method==2) & any(is.null(Winf),is.null(Kw),is.null(T)))
           stop("Method 2 requires Winf, Kw, and T")
    if(any(method==3) & is.null(tmax))
           stop("Method 3 requires tmax")
    if(any(method==4) & any(is.null(tmax),is.null(Kl)))
           stop("Method 4 requires Kl and tmax")
    if(any(method==5) & any(is.null(tm),is.null(Kl)))
           stop("Method 5 requires Kl and tm")
    if(any(method==6) & is.null(GSI))
           stop("Method 6 requires GSI")
    if(any(method==7) & is.null(Wdry))
           stop("Method 7 requires Wdry")
    if(any(method==8) & is.null(Wwet))
           stop("Method 8 requires Wwet")
   if(any(method==9) & any(is.null(Linf),is.null(Kl),is.null(Bl)))
           stop("Method 9 requires Linf, Kl, and Bl")


    n<-length(method)
    if(any(method==3)) n<-n+1
    out<-matrix(NA,n,1L)
    dimnames(out)<-list(rep(NA,n),c("M"))
    cnt<-0
   if(any(method==1)){
      cnt<-cnt+1
      out[cnt,1]<-round(10^(-0.0066-0.279*log10(Linf)+0.6543*log10(Kl)+0.4634*log10(T)),3)
      dimnames(out)[[1]][cnt]<-list("Pauly (1980) - Length Equation")
      if (T<4 || T>30) warning ("Temperature value seems wrong -- <4 or >30")
     }
   if(any(method==2)){ 
       cnt<-cnt+1
       out[cnt,1]<-round(10^(-0.2107-0.0824*log10(Winf)+0.6757*log10(Kw)+0.4627*log10(T)),3)
       dimnames(out)[[1]][cnt]<-list("Pauly (1980) - Weight Equation")
       if (T<4 || T>30) warning ("Temperature value seems wrong -- <4 or >30")
    }
   if(any(method==3)){
        if (tmax<0.5 || tmax>300)
		stop ("Error: maximum age value(s) < 0.5 or > 300.")
        cnt<-cnt+1
        out[cnt,1]<-round(4.22/(tmax^0.982),3)
        dimnames(out)[[1]][cnt]<-list("Hoenig (1983) - Joint Equation")
        cnt<-cnt+1
        out[cnt,1]<-round(exp(1.46 - 1.01*log(tmax)),3)
        dimnames(out)[[1]][cnt]<-list("Hoenig (1983) - Fish Equation")
    }
   if(any(method==4)){
        cnt<-cnt+1 
        out[cnt,1]<-round((3*Kl)/(exp(Kl*(0.38*tmax))-1),3)
        dimnames(out)[[1]][cnt]<-list("Alverson and Carney (1975)")
    }
   if(any(method==5)){
        cnt<-cnt+1 
        out[cnt,1]<-round((3*Kl)/(exp(Kl*tm)-1),3)
        dimnames(out)[[1]][cnt]<-list("Roff (1984)")
     }
   if(any(method==6)){
        cnt<-cnt+1
        out[cnt,1]<-round(0.03+1.68*GSI,3)
        dimnames(out)[[1]][cnt]<-list("Gunderson and Dygert (1988)")
     }
   if(any(method==7)){
        cnt<-cnt+1
        if (Wdry<0.00001 || Wdry>2000) warning ("Dry weight may be outside of range used to derive equation.")
        out[cnt,1]<-round(1.92*(Wdry^-0.25),3)
        dimnames(out)[[1]][cnt]<-list("Peterson and Wroblewski (1984)")
     }
   if(any(method==8)){
      cnt<-cnt+1
      out[cnt,1]<-round(3.0*(Wwet^-0.288),3)
      dimnames(out)[[1]][cnt]<-list("Lorenzen (1996)")
     }
     if(any(method==9)){
      cnt<-cnt+1
      out[cnt,1]<-round(exp(0.55-1.61*log(Bl)+1.44*log(Linf)+log(Kl)),3)
      dimnames(out)[[1]][cnt]<-list("Gislason et al. (2010)")
     }

    return(out)
}

   








