 quad_spline_est<-function(xtab, ytab, x, kn=ceiling((length(xtab))^(1/4)),method="u", all.dea=FALSE)
 {

  #----------------------------------------------------------------------
  #% inputs:
  #     (xtab,ytab) = sample of (nx2) obs
  #       x = vector of evaluation points
  #       kn = number of internal knots
  #       method = "u" (unconstrained estimator),  "m" (montonicity constraint),
  #       "mc" (monotonicity and concavity constraints) 
  #       all.dea = FALSE  # if method="mc", then it keeps all the points 
  #       defining the DEA envelope as knots 
  #----------------------------------------------------------------------

    if(method=="m")
     {cv=0} else{
      if(method=="mc")
       {cv=1}
       else{
       cv<- -1}
      } 
   
    stopifnot(method%in%c("u","m","mc"),!(all.dea & cv%in%c(-1,0)),length(xtab)==length(ytab))

     spl_ord <- 3 # quadratic spline    
     
     # support = [l_end,r_end]
     l_end<-min(x)
     r_end<-max(x)

    if(method=="u")
    {dea_fdh_index<- 1:(length(xtab))
    }else{
     dea_fdh_index<-which(dea(xtab,ytab,RTS=cv,ORIENTATION="out")$eff==1)
    } 
    
     # change the number of knots to select all the DEA enveloppe points
     kn<-ifelse(all.dea,length(dea_fdh_index),kn)
     
     lgrid<-(r_end-l_end)/(kn+1)
  
     if(!all.dea)
     {knots <- c(l_end,quantile(xtab[dea_fdh_index],prob=1:kn/(kn+1),type=7),r_end)
     }
     else
     {knots <- c(l_end,sort(xtab[dea_fdh_index]),r_end)
     }
  
     rind<-which(diff(knots)<0.001)+1
  
     if(length(rind>0))
     {knots<-knots[-rind]}
  
     kn<-length(knots)-2 # corrected number of internal knots
     knots <- c(-2*lgrid+l_end, -lgrid+l_end, knots, r_end+lgrid, r_end+2*lgrid)

     x_cal <- seq(l_end,r_end, length.out=5001)
     der <- rep(0,length(x_cal))
     bb <- spline.des(knots, ord =spl_ord, derivs=der, x=x_cal, outer.ok = TRUE)

     opt_coef<-apply(bb$design,2,sum)/(length(x_cal)-1)

     # design matrix
     desM<- spline.des(knots, ord =spl_ord, x=xtab, outer.ok = TRUE)$design


     # monotonicity constraint
     nknots<-knots[c(spl_ord:(length(knots)-(spl_ord-1)))]
     der_n<-rep(1,length(nknots))
     monoC<- spline.des(knots, ord =spl_ord, derivs=der_n, x=nknots, outer.ok = TRUE)$design


     if(method=="u")
     {mat<-desM
     rhs<-ytab
     dir <- rep(">=",length(xtab))
     }

     if (cv==0)
     {
     # formulation of linear programming
     mat<-rbind(desM,monoC)
     rhs<-c(ytab,rep(0,length(nknots)))
     dir <- c(rep(">=",length(xtab)),rep(">=",length(nknots)))
     }

     if (cv==1)
     {
      # concavity constraint
      mknots<-NULL
      dfd<-min(diff(nknots))/4
      convC<-NULL

      for ( ni in 1:(length(nknots)-1))
      {
      mknots<-c(mknots,mean(nknots[ni:(ni+1)]))
      tk<-c(mean(nknots[ni:(ni+1)])-dfd,mean(nknots[ni:(ni+1)])+dfd)
      der_tk<-rep(1,length(tk))
      Bm<-spline.des(knots, ord =spl_ord, derivs=der_tk, x=tk, outer.ok = TRUE)$design
      convC<-rbind(convC,apply(Bm,2,diff)/(2*dfd))
      }

      # formulation of linear programming
      mat<-rbind(desM,monoC,convC)
      rhs<-c(ytab,rep(0,length(nknots)),rep(0,length(mknots)))
      dir <- c(rep(">=",length(xtab)),rep(">=",length(nknots)),rep("<=",length(mknots)))
    }
          
      bounds <- list(lower = list(ind = 1:length(opt_coef), val = rep(-Inf,length(opt_coef))),upper=list(ind = 1:length(opt_coef), val = rep(Inf,length(opt_coef))))
      Sol<-Rglpk_solve_LP(opt_coef, mat, dir, rhs, bounds,types=NULL,max=FALSE)
      
      OPT<-Sol$sol

      xgrid <- x
      der <- rep(0,length(xgrid))
      fitt<-spline.des(knots, ord =spl_ord, derivs=der, x=xgrid, outer.ok = TRUE)$design %*% OPT
        
    return(fitt)
 }