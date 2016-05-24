Spmlficmcm <-
function(fl,N,gmname,gcname,DatfE,typ,start,p=NULL){ 
                  # ============================================================ 
                  # Etape 1 estimation des valeurs initiales 
                  # Valeurs initiales des parametre du modele 
                  # baterie des tests
                  yname = all.vars(fl)[1]
                  if(missing(N))
                  {
                  if (is.null(p))
                  {
                  	print(p)
                  stop("Missing prevalence or N=c(N0,N1)")
                    }
                  else
                  {
                  	if (p > 0.5) stop ("Disease prevalence needs to be <= 0.5")
                  	if (p < 0) stop ("Negative disease prevalence")
                  	clb<- model.frame(fl, data = DatfE)
                    # extraction de la variable reponse
                      outcb<-model.extract(clb,"response")
                    # nombre de cas
                  	n1 = sum(outcb)
                  	N1 = 5*n1 
                  	N0 = round(N1 * (1-p)/p)
                  	N<-c(N0,N1)
                    #print (N)
                    }
                    }
                  #else
                  #{
                    # modele
                    clb1<- model.frame(fl, data = DatfE)
                    # extraction de la variable reponse
                    outcb1<-model.extract(clb1,"response")
                    # conditions
                    n1<-sum(outcb1)
                    n0<-length(outcb1)-n1
                    if(N[1]<n0 | N[2]<n1){print(N)
                     stop("The cases number or the controls in the study population is lower than cases or controls of study sample")
                    }else{
                    genom<-DatfE[gmname];genom1<-genom[is.na(genom)!=TRUE]
                    genoen<-DatfE[gcname];genoen1<-genoen[is.na(genoen)!=TRUE]
                    nagm<-genom[is.na(genom)==TRUE]
                    dat<-DatfE[is.na(DatfE[gcname])!=TRUE,]
                    gt1<-dat[gmname];
                    gt2<-dat[gcname];
                    teg1<-ifelse(gt1==0 & gt2==2,1,0)
                    teg2<-ifelse(gt1==2 & gt2==0,1,0)
                    tegk<-teg1+teg2
                    
                    if(length(nagm)>0){print(gmname)
                                      stop("missing mother genotype value")} 
                    if(min(genom1)<0){print(gmname)
                                 stop("mother's genotype is negative")}
                    if(max(genom1)>2){print(gmname)
                                 stop("mother's genotype is greater than 2")}
                    if(min(genoen1)<0){print(gmname)
                                 stop("child's genotype is negative")}
                    if(max(genoen1)>2){print(gmname)
                                 stop("child's genotype is greater than 2")} 
                    if(max(tegk)>0){print(gmname)
                               stop("mother and child genotypes are not compatible")
                               }else{
                               if(typ==1){
                               vIn<-Est.Inpar(fl,N,gmname,gcname,DatfE,1)
                               }else{
                               vIn<-Est.Inpar(fl,N,gmname,gcname,DatfE,2)
                                }   
                    if(missing(start)){parms<-vIn$parms
                                             }else{parms<-start}
                    p = length(parms)
                    beta.start=parms[1:(p-1)];
                    theta.start=parms[p]
                  # Valeurs initiales du systeme d quation non linaire
                  ma.u<-vIn$ma.u
                  vecma.u=c(ma.u[1,],ma.u[2,]) 
                  #=============================================================
                  # Etape 2 resolution du systeme d equation  
                  RSeq<-Nlsysteq(fl,DatfE,N,gmname,gcname,yname,beta.start,theta.start)
                  SS<-nleqslv(vecma.u,RSeq) 
                  vecma.u<-SS$x
                  #=============================================================
                  # Etape 3 ecriture de la vraisemblance profile
                  if(typ==1){
                    ftlh<-ft_likhoodCas1(fl,DatfE,N,gmname,gcname,yname,vecma.u)
                            }else{
                    ftlh<-ft_likhoodCasM(fl,DatfE,N,gmname,gcname,yname,vecma.u)
                                  }        
                  #=============================================================
                  # Etape 4 calcul des estimateurs  
                  # calcul du gradient
                   if(typ==1){
				   fctgrad<-ft_gradientCas1(fl,DatfE,N,gmname,gcname,yname,vecma.u)
				   }
				   else{
				   fctgrad<-ft_gradientCasM(fl,DatfE,N,gmname,gcname,yname,vecma.u)
				   }	
				   Grad<-fctgrad(parms)
                   #calcul de la hessian
                   delta <- 1e-5;
                   hess <- matrix(0, p, p);
                   for(gg in 1:p){
                   delta_ggamma <- parms;
                   delta_ggamma[gg]<- delta_ggamma[gg] + delta;
                   hess[, gg]<-(fctgrad(delta_ggamma) - Grad)/delta;
                   }
                  
                   # estimateur des parametres
                   Parms.est=parms-solve(hess)%*%Grad
                   # calcul de la variance
    
				   Grad1<-fctgrad(Parms.est)
				   # Code de debuggage
				   # print(Grad1)
                   Hes1 <- matrix(0, p, p);
                   for(gg in 1:p){
                   delta_ggamma <- Parms.est;
                   delta_ggamma[gg] <- delta_ggamma[gg] + delta;
                   Hes1[, gg] <- (fctgrad(delta_ggamma) - Grad1)/delta;
                   }
                  matv<-(-1)*solve(Hes1)
                  var.par<-sqrt(diag(matv))
                  #=============================================================
                  # Etape 5 preparation des resultats
                  mats<-cbind(Parms.est,var.par)
                  nma<-c("Intercept",attr(terms.formula(fl),"term.labels"),"theta")
                  nac<-c("Estimate","Std.Error")
                  colnames(mats)<-nac
                  rownames(mats)<-nma
                  loglik<-ftlh(mats[,"Estimate"])
                  rr<-list(N=N,Uim=vecma.u,MatR=mats,Matv=matv,Lhft=ftlh,Value_loglikh=loglik)
                  return(rr)
                  }
                  }
                    #}
                  }
