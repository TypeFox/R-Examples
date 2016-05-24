 loc_est<-function(xtab, ytab, x, h, method="u")
 {
   stopifnot(length(xtab)==length(ytab), method%in%c("u","m"))
      
   fitt<-rep(0,length(x))

   for (i in 1:length(x))
   {
   yloc<-ytab[(rep(x[i],length(xtab))-h<=xtab) & (xtab<=rep(x[i],length(xtab))+h)]
   xloc<-xtab[(rep(x[i],length(xtab))-h<=xtab) & (xtab<=rep(x[i],length(xtab))+h)]

   if(length(xloc)!=0)
   { l_l_end<- x[i]-h
     r_l_end<- x[i]+h

     opt_coef<-c((r_l_end-l_l_end),0.5*((r_l_end-x[i])^2-(l_l_end-x[i])^2))
     desM<- cbind(1,xloc-x[i]) # design matrix

     # formulation of linear programming
     obj<-opt_coef
     mat<-desM
     rhs<-yloc
     dir <- c(rep(">=",length(xloc)))

     if(method=="u")
     {bounds <- list(lower = list(ind = 1:2, val = c(-Inf,-Inf)),upper=list(ind = 1:2, val = c(Inf,Inf)))
     }
     else
     {bounds <- list(lower = list(ind = 1:2, val = c(-Inf,0)),upper=list(ind = 1:2, val = c(Inf,Inf)))
     }     
     Sol<-Rglpk_solve_LP(obj, mat, dir, rhs, bounds,types=NULL,max=FALSE)
     OPT<-Sol$sol

     fitt[i]<-OPT[1]
     }
     else
     {fitt[i]<-NA
     }
   }
   return(fitt)
 }