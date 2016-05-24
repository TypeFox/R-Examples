Est.Inpar <-
function(fl,N,gnma,gnch,tab1,typ,p=NULL){
                            # fl formule of the equation
                            # tab1 : database
                            # N =(N0,N1) le nombres individus eligibles N0 temoins et N1 cas
                            # gnma : genotype de la mere 
                            # gnch : genotype de enfant
                            # typ : indique si nous avons les donnees manquantes ou non. 1 indique pas des donnees manquantes
                            # p : la prevalence de la maladie 
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
			                  	clb<- model.frame(fl, data = tab1)
			                    # extraction de la variable reponse
			                      outcb<-model.extract(clb,"response")
			                    # nombre de cas
			                  	n1 = sum(outcb)
			                  	N1 = 5*n1 
			                  	N0 = round(N1 * (1-p)/p)
			                  	N = c(N0,N1)
			                    }
			                    }
			                N0<-N[1];N1<-N[2];
                            varz0<-all.vars(fl)[-1];vrze<-varz0[-which(varz0%in%c(gnma,gnch))]
                            outc=as.character(attr(terms.formula(fl),"variables")[[2]]) 
                            fit<-glm(fl,data=tab1,family=binomial)
                            if(typ==1){tab=tab1}else{tab=tab1[is.na(tab1[gnch])!=TRUE,]}
                            n1<-dim(tab[tab[outc]==1,])[1];
                            n0<-dim(tab[tab[outc]==0,])[1]  
                            a=coef(fit)[1]+log(N1/n1)-log(N0/n0)
                            beta.start=c(a,coef(fit)[-1]);
                            # calcul de la valeur initiale de la distribution du genotype
                            d=c(N0/n0,N1/n1)
                            # ecriture de l equation (3)
                            fgt<-fgp_mf1(tab,d,gmname=gnma,gcname=gnch,outc=outc)     
                            ss <- nleqslv(0.1,fgt)
                            theta.start=ss$x;
                            # solution des parametres finaux 
                            parms=c(beta.start,theta.start)
                            # calcul des valeurs intiales du systeme non lineair  
                            lst_Matin2<-fct_invCap(tab,N,outc,vrze,gnma,gnch,theta.start)
                            Mat.sup=lst_Matin2$brs
                            return(list(parms=parms,ma.u=Mat.sup))
                            }
