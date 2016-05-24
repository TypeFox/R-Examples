## group-sequential repeated confidence bounds
cb.r.gsd <-function(GSD, GSDo, level=NULL) {
    if(!is.null(level) && (level != GSD$al)) {
        GSD$al <- level
        GSD$b  <- compBounds(t=GSD$t/GSD$t[GSD$K], t2 = GSD$t*GSD$Imax, iuse = GSD$SF, asf = NULL, 
                             alpha = GSD$al, phi=ifelse(is.null(GSD$phi),0,GSD$phi), 
                             ztrun = 8)
    }
    (GSDo$z - GSD$b[GSDo$T]) / sqrt(GSD$t[GSDo$T]*GSD$Imax)
}


cb.r.ad <- function(pT,iD,sT,sTo,level=NULL){
################################################################################
# sT$SF = 6, b_(j,u) according to stage wise ordering (exact second stage p-value; no flexible monitoring in secondary trial)
# sT$SF = 7, b_{j,u} according the ordering of the repeated confidence bound (not implemented yet)
################################################################################
          pT$K<-length(pT$t);
          sT$K<-length(sT$t);
          if(ifelse(!is.null(level),level!=pT$al,FALSE)){
                 pT$al <- level;
                 pT$b  <- compBounds(t=pT$t/pT$t[pT$K], t2 = pT$t*pT$Imax, iuse = pT$SF, asf = NULL, 
                                   alpha = pT$al, phi =ifelse(is.null(pT$phi),0,pT$phi), 
                                   ztrun = 8)#$upper.bounds;
                  sT$ce0 <- rCER(0,pT,iD);
                  if(sT$K>1){
                    sT$b  <-compBounds(t=sT$t/sT$t[sT$K], t2 = sT$t*sT$Imax, iuse = sT$SF, asf = NULL, 
                                   alpha = sT$ce0, phi =ifelse(is.null(sT$phi),0,sT$phi), 
                                   ztrun = 8)#$upper.bounds;
                   } else { sT$b <- qnorm(1-rCER(h=0,pT,iD)) }
          };
          
          if(sTo$z>=sT$b[sTo$T]){ hl<-0; hu<-(sTo$z-min(qnorm(1-pT$al),sT$b[sTo$T]))/sqrt(sT$t[sTo$T]*sT$Imax)
             }else{ hu<-0; hl<-(sTo$z-sT$b[sTo$T])/sqrt(sT$t[sTo$T]*sT$Imax) };
          if(sT$K==1) {
             res=uniroot(function(h){ sTo$z-h*sqrt(sT$t[sTo$T]*sT$Imax)-qnorm(1-rCER(h,pT,iD))
                                },c(hl,hu)
             )$root
          } else {
            if(ifelse(!is.null(sT$SF),sT$SF==6,FALSE)){
               res=uniroot(function(h){P.so.gsd(h,GSD=sT)-rCER(h,pT,iD)},c(hl,hu))$root
            } else  {
                 res=uniroot(function(h){
                       sTo$z-h*sqrt(sT$t[sTo$T]*sT$Imax)-
                         compBounds(t=sT$t/sT$t[sT$K], t2 = sT$t*sT$Imax, iuse = sT$SF, asf = NULL, 
                                   alpha = rCER(h,pT,iD), phi =ifelse(is.null(sT$phi),0,sT$phi), 
                                   ztrun = 8)[sTo$T]#$upper.bounds[sTo$T]                         
                       },c(hl,hu)
                  )$root
               }
            }
return(res)
}

## group-sequential stage-wise ordering confidence bounds
cb.so.gsd <- function(GSD, GSDo, level=NULL) {
    GSD$K <- length(GSD$t)
    
    if(!is.null(level) && (GSD$al != level)) {
        GSD$al <- level
        GSD$alab <- comp.alab(GSD=GSD)
    } else {
        if(is.null(GSD$alab)) GSD$alab <- comp.alab(GSD=GSD) 
    }
    
    if(GSDo$T == 1) res <- ((GSDo$z - qnorm(1-GSD$al)) / sqrt(GSD$t[GSDo$T]*GSD$Imax))
    else {
        if(GSDo$z == GSD$b[GSDo$T]) res <- GSD$alab[GSDo$T]
        else {
            if(GSDo$z > GSD$b[GSDo$T]) hl <- GSD$alab[GSDo$T]
            else {
                if(is.null(GSD$als)) GSD$als <- comp.als(GSD=GSD)
                hl <- (GSDo$z - qnorm(1-GSD$al+GSD$als[GSD$K-1])) / sqrt(GSD$t[GSDo$T]*GSD$Imax)
            }

            res <- uniroot(f=function(h) { seqmon(a=GSD$a[1:GSDo$T],
                               b=pbounds(h=h,pT=list(t=GSD$t[1:GSDo$T],b=c(GSD$b[1:(GSDo$T-1)],GSDo$z),Imax=GSD$Imax)),
                               t=GSD$t[1:GSDo$T]*GSD$Imax,int=500*array(c(1),GSDo$T))[2*GSDo$T] - GSD$al},
                           interval=c(hl, GSD$alab[GSDo$T-1])
                          )$root
                         
        }
    }   
    res
}


cb.so.ad <- function(pT, iD, sT, sTo, level=NULL) {
    error <- 0;

                                        #if(!is.loaded(symbol.C("AGSDest"))) {
                                        #lib.file <- file.path(paste("AGSDest", .Platform$dynlib.ext, sep=""))
                                        #dyn.load(paste("../library/AGSDest/libs/",lib.file, sep=""))
                                        #cat(" -Loaded ", lib.file, "\n")
                                        #}
    
    if(is.null(level)) {
        level=pT$al
        theta_1=NULL;
    } else theta_1 <- comp.alab(GSD=list(a=pT$a,b=pT$b,t=pT$t,al=level,Imax=pT$Imax))

    if(iD$T > length(pT$a)) {
        print("iD$T > number stages pT")
        error=1
    }
    
    if(sTo$T > length(sT$a)) {
        print("sTo$T > number stages sT")
        error=1
    }

    k_1 <- length(pT$a)

    if(is.null(pT$a)) {
        print("pT$a is missing")
        error=1
    } else a_1 <- pT$a
    
    if(is.null(pT$b)) {
        print("pT$b is missing");
        error=1
    } else b_1 <- pT$b
    
    if(is.null(pT$t)) {
        print("pT$t is missing")
        error=1
    } else t_1 <- pT$t*pT$Imax

    if(is.null(pT$alab)) {
        print("pT$alab is missing")
        error=1
    } else {
        if(is.null(theta_1)) theta_1 <- pT$alab
    }

    k_2 <- length(sT$a);
    
    if(is.null(sT$a)) {
        print("sT$a is missing")
        error=1
    } else a_2 <- sT$a
    
    if(is.null(sT$b)) {
        print("sT$b is missing")
        error=1
    } else b_2 <- sT$b
    
    if(is.null(sT$t)) {
        print("sT$t is missing")
        error=1
    } else t_2 <- sT$t*sT$Imax

    if(is.null(iD$T)) {
        print("iD$T is missing")
        error=1
    } else T_1 <- iD$T
    
    if(is.null(iD$z)) {
        print("iD$z is missing")
        error=1
    } else z_1 <- iD$z

    if(is.null(sTo$T)) {
        print("sTo$T is missing")
        error=1
    } else T_2 <- sTo$T

    if(is.null(sTo$z)) {
        print("sTo$z is missing")
        error=1
    } else zT_2 <- sTo$z

    e<-c(0,0);
    out <- 0

    if(error == 0) {
        out <- .C(mainf,
                  k = as.integer(k_1),
                  a = as.numeric(a_1),
                  b = as.numeric(b_1),
                  t = as.numeric(t_1),
                  alab = as.numeric(theta_1),
                  al = as.numeric(level),
                  T = as.integer(T_1),
                  z = as.numeric(z_1),
                  k2 = as.integer(k_2),
                  a2 = as.numeric(a_2),
                  b2 = as.numeric(b_2),
                  t2 = as.numeric(t_2),
                  T2 = as.integer(T_2),
                  zT = as.numeric(zT_2),
                  erg = as.numeric(e), PACKAGE="AGSDest")
        out$erg[2]
    } else {
        print("Values are missing.")
        return(0);
    }
}

