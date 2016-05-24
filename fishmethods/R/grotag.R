grotag<-function(L1=NULL,L2=NULL,T1=NULL,T2=NULL,alpha=NULL,beta=NULL,
               design=list(nu=0,m=0,p=0,sea=0),
               stvalue=list(sigma=0.9,nu=0.4,m=-1,p=0.01,u=0.4,w=0.4),
               upper=list(sigma=5,nu=1,m=2,p=1,u=1,w=1),
			lower=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0),
              gestimate=TRUE,st.ga=NULL,st.gb=NULL,st.galow=NULL,st.gaup=NULL,
		    st.gblow=NULL,st.gbup=NULL,
              control=list(maxit=10000)){
     if(is.null(L1)) stop ("L1 vector does not exist")
     if(is.null(L2)) stop ("L2 vector does not exist")
     if(is.null(T1)) stop ("T1 vector does not exist")
     if(is.null(T2)) stop ("T2 vector does not exist")
     if(is.null(alpha)) stop ("alpha does not exist")
     if(is.null(beta)) stop ("beta does not exist")
     if(alpha >= beta) stop("Error: parameter alpha must be smaller than beta")
     if(length(L1)!=length(L2)) stop ("L1 and L2 vectors of different lengths")
     if(length(T1)!=length(T2)) stop ("T1 and T2 vectors of different lengths")
     if(length(L1)!=length(T1)) stop ("L1 and T1 vectors of different lengths")
     if(length(L2)!=length(T2)) stop ("L2 and T2 vectors of different lengths")
     if(gestimate==FALSE & (is.null(st.ga)|is.null(st.gb))) stop("Enter values for st.ga and st.gb")
     if(gestimate==FALSE & (is.null(st.galow)|is.null(st.gaup))) stop("Enter values for st.galow and st.gaup")
     if(gestimate==FALSE & (is.null(st.gblow)|is.null(st.gbup))) stop("Enter values for st.gblow and st.gbup")

    # Test if design values are either zero or 1
       if(any(!design$nu %in% c(0,1),!design$m %in% c(0,1),
            !design$p %in% c(0,1),!design$sea %in% c(0,1))) stop("Only 0 or 1 allowed in design list")
        
    #Compare default design to entered design, change 
          defdesign=list(nu=0,m=0,p=0,u=0,w=0)
          if(any(names(design)=="m")){ 
            if(defdesign$m!=design$m) defdesign$m<-design$m
           }
            if(any(names(design)=="p")){ 
             if(defdesign$p!=design$p) defdesign$p<-design$p
            }
           if(any(names(design)=="nu")){
            if(defdesign$nu!=design$nu) defdesign$nu<-design$nu
          }
           if(any(names(design)=="sea")){
            if(design$sea==1){
                     defdesign$u<-1
                     defdesign$w<-1
            }
          }
          defdesign$sigma<-1
            
    #Compare stvalues, change 
		stdef=list(sigma=0.5,nu=0,m=0,p=0.01,u=0,w=0)
         if(any(names(stvalue)=="sigma")){ 
 		  if(stdef$sigma!=stvalue$sigma) stdef$sigma<-stvalue$sigma
          }
         if(any(names(stvalue)=="nu")){ 
   		  if(stdef$nu!=stvalue$nu) stdef$nu<-stvalue$nu
          }
         if(any(names(stvalue)=="m")){ 
 		  if(stdef$m!=stvalue$m) stdef$m<-stvalue$m
         }
         if(any(names(stvalue)=="p")){ 
		  if(stdef$p!=stvalue$p) stdef$p<-stvalue$p
         }
         if(any(names(stvalue)=="u")){ 
 		 if(stdef$u!=stvalue$u) stdef$u<-stvalue$u
         }
         if(any(names(stvalue)=="w")){ 
  		  if(stdef$w!=stvalue$w) stdef$w<-stvalue$w
         }
          
   #Compare maxvalues, change 
 	  maxdef=list(sigma=5,nu=5,m=2,p=1,u=1,w=1)
        if(any(names(upper)=="sigma")){
	      if(maxdef$sigma!=upper$sigma) maxdef$sigma<-upper$sigma
        }
       if(any(names(upper)=="nu")){
	  	if(maxdef$nu!=upper$nu) maxdef$nu<-upper$nu
       }
       if(any(names(upper)=="m")){
	    if(maxdef$m!=upper$m) maxdef$m<-upper$m
       }
       if(any(names(upper)=="u")){
	    if(maxdef$u!=upper$u) maxdef$u<-upper$u
       }
       if(any(names(upper)=="w")){
 	    if(maxdef$w!=upper$w) maxdef$w<-upper$w
       }
       if(any(names(upper)=="p")){
 	    if(maxdef$p!=upper$p) maxdef$p<-upper$p
       }
       maxs<-c(defdesign$sigma*maxdef$sigma,defdesign$nu*maxdef$nu,
           defdesign$m*maxdef$m,defdesign$p*maxdef$p,defdesign$u*maxdef$u,
           defdesign$w*maxdef$w)

    #Compare minvalues, change 
	  mindef=list(sigma=0,nu=0,m=-2,p=0,u=0,w=0)
       if(any(names(lower)=="sigma")){
         if(mindef$sigma!=lower$sigma) mindef$sigma<-lower$sigma
       }
      if(any(names(lower)=="nu")){
         if(mindef$nu!=lower$nu) mindef$nu<-lower$nu
       }
       if(any(names(lower)=="m")){
         if(mindef$m!=lower$m) mindef$m<-lower$m
       }
       if(any(names(lower)=="u")){
         if(mindef$u!=lower$u) mindef$u<-lower$u
       }
       if(any(names(lower)=="w")){
         if(mindef$w!=lower$w) mindef$w<-lower$w
       }
       if(any(names(lower)=="p")){
         if(mindef$p!=lower$p) mindef$p<-lower$p
       }
#Check upper and inital values     
   if(defdesign$sigma==1){
      if(maxdef$sigma<stdef$sigma) stop("upper sigma is < starting sigma")
      if(maxdef$sigma==stdef$sigma) warning("upper sigma and starting sigma are equal")
     }
    if(defdesign$nu==1){
      if(maxdef$nu<stdef$nu) stop("upper nu is < starting nu")
      if(maxdef$nu==stdef$nu) warning("upper nu and starting nu are equal")
     }
 	if(defdesign$u==1){
      if(maxdef$u<stdef$u) stop("upper u is < starting u")
      if(maxdef$u==stdef$u) warning("upper u and starting u are equal")
     }
	if(defdesign$w==1){
      if(maxdef$w<stdef$w) stop("upper w is < starting w")
      if(maxdef$w==stdef$w) warning("upper w and starting w are equal")
     }
	if(defdesign$m==1){
      if(maxdef$m<stdef$m) stop("upper m is < starting m")
      if(maxdef$m==stdef$m) warning("upper m and starting m are equal")
     }
    if(defdesign$p==1){
      if(maxdef$p<stdef$p) stop("upper p is < starting p")
      if(maxdef$p==stdef$p) warning("upper p and starting p are equal")
     }

 #Check lower and inital values     
   if(defdesign$sigma==1){
      if(mindef$sigma>stdef$sigma) stop("lower sigma is > starting sigma")
      if(mindef$sigma==stdef$sigma) warning("lower sigma and starting sigma are equal")
     }
    if(defdesign$nu==1){
      if(mindef$nu>stdef$nu) stop("lower nu is >  starting nu")
      if(mindef$nu==stdef$nu) warning("lower nu and starting nu are equal")
     }
 	if(defdesign$u==1){
      if(mindef$u>stdef$u) stop("lower u is > starting u")
      if(mindef$u==stdef$u) warning("lower u and starting u are equal")
     }
	if(defdesign$w==1){
      if(mindef$w>stdef$w) stop("lower w is > starting w")
      if(mindef$w==stdef$w) warning("lower w and starting w are equal")
     }
	if(defdesign$m==1){
      if(mindef$m>stdef$m) stop("lower m is > starting m")
      if(mindef$m==stdef$m) warning("lower m and starting m are equal")
     }
     if(defdesign$p==1){
      if(mindef$p>stdef$p) stop("lower p is > starting p")
      if(mindef$p==stdef$p) warning("lower p and starting p are equal")
     }

           
#get data 
     x<-as.data.frame(cbind(L1,T1,L2,T2)) 
     x<-x[!is.na(x$L1) & !is.na(x$T1) & !is.na(x$L2) & !is.na(x$T2),]  
   
#calculate deltaL and deltaT
     x$delta.L<-x$L2-x$L1
     x$dt<-x$T2-x$T1
     
#calcalate ga and gb
    if(gestimate==TRUE){
	  ga.init <-  with(subset(x, L1 > (0.8 * alpha) & L1 < (1.2 * alpha) & dt > median(dt)), mean(delta.L / dt))
       gb.init <-  with(subset(x, L1 > (0.8 * beta)  & L1 < (1.2 * beta)  & dt > median(dt)), mean(delta.L / dt))
     }
    if(gestimate==FALSE){
      ga.init<-st.ga
      gb.init<-st.gb
    }

#set up design
   ests<-c(1,1,defdesign$sigma*1,defdesign$nu*1,defdesign$m*1,defdesign$p*1,defdesign$u*1,defdesign$w*1)
   allparms<-as.data.frame(cbind(seq(1,length(ests),1),ests))
   names(allparms)<-c("pars","est")
   if(gestimate==TRUE){
     lower1<-c(0.5 * ga.init,0.5 * gb.init,
                 mindef$sigma,mindef$nu,mindef$m,mindef$p,mindef$u,
                 mindef$w)
     upper1<-c( 1.5 * ga.init,1.5 * gb.init,maxs[1],maxs[2],
             maxs[3],maxs[4],maxs[5],maxs[6])
   }
   if(gestimate==FALSE){
     lower1<-c(st.galow,st.gblow,
                 mindef$sigma,mindef$nu,mindef$m,mindef$p,mindef$u,
                 mindef$w)
     upper1<-c(st.gaup,st.gbup,maxs[1],maxs[2],
             maxs[3],maxs[4],maxs[5],maxs[6])
   }




   getlow<-cbind(allparms,lower1);getlow<-getlow[getlow$est==1,3]
   gethigh<-cbind(allparms,upper1);gethigh<-gethigh[gethigh$est==1,3]
   trall<-allparms[allparms$est>0,]	

   parms<-c(ga.init,gb.init,stdef$sigma,
         ifelse(defdesign$nu==0,NA,stdef$nu),
         ifelse(defdesign$m==0,NA,stdef$m),
         ifelse(defdesign$p==0,NA,stdef$p),
         ifelse(defdesign$u==0,NA,stdef$u),
         ifelse(defdesign$w==0,NA,stdef$w))
   parms<-parms[!is.na(parms)]

# Likelihood function
   f1<-function(p){
        ee<-cbind(trall,p)
        ff<-merge(allparms,ee,by.x="pars",by.y="pars",all.x=TRUE,all.y=TRUE) 
        ff$p<-ifelse(is.na(ff$p),0,ff$p)
        R <-range(x$delta.L)[2] - range(x$delta.L)[1] 
        phi1<-ff$p[7]*(sin(2*pi*(x$T1-ff$p[8])))/(2*pi)
        phi2<-ff$p[7]*(sin(2*pi*(x$T2-ff$p[8])))/(2*pi)
 	   mu<-((beta * ff$p[1] - alpha * ff$p[2])/(ff$p[1] - ff$p[2]) - x$L1) * (1 - ( 1 + (ff$p[1] - ff$p[2]) / (alpha - beta)) ^ (x$dt+(phi2-phi1)))
        std.dev<-mu
        std.dev<-replace(std.dev, which(std.dev < 0), ff$p[4])
        sdev<-sqrt((ff$p[4] * std.dev)^2+ ff$p[3]^2)
       -1*sum(log((1-ff$p[6]) * dnorm(x = x$delta.L, 
     	    mean = ((beta * ff$p[1] - alpha * ff$p[2])/(ff$p[1] - ff$p[2]) - x$L1) * (1 - ( 1 + (ff$p[1] - ff$p[2]) / (alpha - beta)) ^ (x$dt+(phi2-phi1)))+ff$p[5],
      		sd = sdev) + ff$p[6] / R)) 
     }
    mod1<-try(optim(par=parms,fn=f1,lower=getlow,upper=gethigh,
          method="L-BFGS-B",hessian=TRUE, control=control),silent=TRUE)

  if(class(mod1)!="try-error"){
     AIC<-round(2*mod1$value+2*length(mod1$par),1)
     SE<-sqrt(diag(solve(mod1$hessian)))
      result = as.data.frame(matrix(NA, ncol=1, nrow = 12))
      dimnames(result)[[1]] = c("Parameters","Mean growth rates   ga","Mean growth rates   gb", 
       "Seasonal variation  u","                    w","Growth variability  nu", "Measurement error   s", "                    m", "Outliers            p"," ","-Log likelihood","AIC")

      #Get estimated parms 
      ee<-cbind(trall,mod1$par)
      names(ee)<-c("pars","est","p")
      kk<-merge(allparms,ee,by.x="pars",by.y="pars",all.x=TRUE,all.y=TRUE) 
      kk$p<-ifelse(is.na(kk$p),0,kk$p)
      result = cbind(result[,-1], Estimate = c(NA,round(kk$p[1],2), round(kk$p[2],2),round(kk$p[7],3),round(kk$p[8],3),round(kk$p[4],3),round(kk$p[3],3),round(kk$p[5],3),round(kk$p[6],3),NA,round(mod1$value,1),AIC))
      phi1<-kk$p[7]*(sin(2*pi*(x$T1-kk$p[8])))/(2*pi)
      phi2<-kk$p[7]*(sin(2*pi*(x$T2-kk$p[8])))/(2*pi)
      pred<-((beta * kk$p[1] - alpha * kk$p[2])/(kk$p[1] - kk$p[2]) - x$L1) * (1 - ( 1 + (kk$p[1] - kk$p[2]) / (alpha - beta)) ^ (x$dt+(phi2-phi1)))
      resid<-x$delta.L-pred
      #Get estimated SE 
      ee<-cbind(trall,SE)
      names(ee)<-c("pars","est","p")
      kk<-merge(allparms,ee,by.x="pars",by.y="pars",all.x=TRUE,all.y=TRUE) 
      kk$p<-ifelse(is.na(kk$p),0,kk$p)
      result = cbind(result, SE = c(NA,round(kk$p[1],2), round(kk$p[2],2),round(kk$p[7],3),round(kk$p[8],3),round(kk$p[4],3),round(kk$p[3],3),round(kk$p[5],3),round(kk$p[6],4),NA,NA,NA))
	Linf<-(beta*mod1$par[1]-alpha*mod1$par[2])/(mod1$par[1]-mod1$par[2])*-1
      K<--1*log(1+(mod1$par[1]-mod1$par[2])/(alpha-beta))
      result1 = as.data.frame(matrix(NA, ncol=1, nrow = 3))
      dimnames(result1)[[1]] = c("Parameters","Linf","K ")
      names(result1)<-"Estimate"
      result1$Estimate[2]<-Linf
      result1$Estimate[3]<-K
     dd<-list(result,result1,cov2cor(solve(mod1$hessian)),pred,resid)
     names(dd)<-c("table","VBparms","correlation","predicted","residuals")
     return(dd)
  }

  if(class(mod1)=="try-error"){
    return("Fit Failed.")
   }
}#end function




