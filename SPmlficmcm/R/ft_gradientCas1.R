ft_gradientCas1 <-
function(fl,data1,N,gmname,gcname,yname,vecma.u,HW=TRUE)
                   {
                    # Arguments specifiques  la fonction
                    # d: vecteur des proportion de cas et de temoins echantillonnes 
                    # gmname: nom de la variable contenant le genotype de la mere
                    # gcname: nom de la variable contenant le genotype de l enfant
                    # varz: nom des vriable de l environement c("","")
                    # ma.u est la matrice des parametre u  
                    
                    # creation de deux table une avec tous les donnees complet de c et l autre avec les seules c observables
                      lstdat<-fctcd(data1,gcname,yname)
                      datMod<-lstdat$datdmcp
                      datNo<-lstdat$datnmv
                      
                    # 1-Generer le "frame" du modele
                      cl<- model.frame(fl, data = datMod)
                      vx <- model.matrix(fl,data=datMod)
                    # extarction de la variable reponse
                      outc<-model.extract(cl,"response")
                    # Extraction du vecteur de genotypes de la mere
                      gm <- vx[,gmname]
                    # Extraction du vecteur de genotypes de l enfant
                      gc <- vx[,gcname]
                    
                    # noms des labels 
                      varz0<-all.vars(fl)[-1];varz<-varz0[-which(varz0%in%c(gmname,gcname))]
                    # les valeurs possibles du geotype de la mere
                      gm1<-gm[is.na(gm)!=TRUE]
                      frq<-unique(gm1);np1<-length(frq)
                      noutc=as.character(terms.formula(fl)[[2]])
                    # selon les cas il nous donne la fonction indique a utilise pour le genotype de la mere
                      #indfg<-IndF3(gm)
                      pp<-length(vecma.u)/2;ppd<-2*pp;uu<-pp+1
                      ma.u<-rbind(vecma.u[1:pp],vecma.u[uu:ppd])
                      
                    # 2-Construction du systeme non lineaire =======================================
                    # 2.1-creation de la table A
                        matd<-cbind(outc,vx)
                        np<-dim(vx)[2]
                        d<-vector()
                        d[1]<-N[1]-dim(datNo[datNo[noutc]==0,])[1]
                        d[2]<-N[2]-dim(datNo[datNo[noutc]==1,])[1]
                       
                        #construction de Cjm 
                        mat.cjm<-fpol1(datMod,c(varz,gmname),"vdcop","Cjm")
                          
                    ## la fonction de gradient de la log-vraisemblence 
                    gradient_prof<-function(parms){
                                    beta.start<-parms[1:np];n1<-np+1;theta.start<-parms[n1:length(parms)]
                                    
                                    # nv
                                    eta = vx%*%beta.start
                                    Pijmc<-((1/(exp(-eta)+1))^outc)*((1-1/(exp(-eta)+1))^(1-outc))
                                    dPijmc.eta = (-1)^(1-outc)*S1(eta)
                                    dlogPijmc.eta = outc - S0(-eta)
                                                                                                                                           
                                    # construction de la distribution conditionnelle du genotype de l enfant sachant la mere et celle de la mere selon que nous sommes sous HW ou non
                                    # nv
                                    # prob conditionel de l'enf sachant la mere Pc/m matd
                                    Pgcm<-Prgcm_HW1(matd[,c(gmname,gcname)],theta.start)
                                    # Derivee de la prob conditionel de l'enf sachant la mere                                    
                                    dPgcm.theta<-dPrgcm_HW1.theta(matd[,c(gmname,gcname)],theta.start)
                                    dlogPgcm.theta<-dlogPrgcm_HW1.theta(matd[,c(gmname,gcname)],theta.start)
                                                                                                            
                                    # calcul de la fonction hijmc 
                                    Hijmc<-Pijmc*Pgcm
                                    # calcul des derivees de la fonction hijmc                                    
                                    dHijmc.eta<-dPijmc.eta * Pgcm
                                    dHijmc.theta<-Pijmc * dPgcm.theta   
                                                                     
                                    nva<-vx[,varz]
                                    
                                    # calcul du genotype 
                                    Pgm<-Prgm_HW1(matd[,gmname],theta.start)
                                    # calcul de sa derivee
                                    dPgm.theta<-dPrgm_HW1.theta(matd[,gmname],theta.start)
                                                                        
                                    vdcop<-datMod["vdcop"]
                                    
                                    # data.frame A
                                    # On a besoin des variables de vx pour le gradient
                                    #nam<-c("outc",varz,gmname,gcname,"vdcop","Pijmc","Pgcm","Pgm","Hijmc","dPijmc.eta", "dlogPijmc.eta","dPgcm.theta","dlogPgcm.theta","dPgm.theta","dHijmc.eta","dHijmc.theta",colnames(vx))
                                    #matA.comp0<-data.frame(outc,nva,gm,gc,vdcop,Pijmc,Pgcm,Pgm,Hijmc,dPijmc.eta,dlogPijmc.eta,dPgcm.theta, dlogPgcm.theta,dPgm.theta,dHijmc.eta,dHijmc.theta,vx)
                                    nam<-c("outc","vdcop","Pijmc","Pgcm","Pgm","Hijmc","dPijmc.eta", "dlogPijmc.eta","dPgcm.theta","dlogPgcm.theta","dPgm.theta","dHijmc.eta","dHijmc.theta",colnames(vx))
                                    matA.comp0<-data.frame(outc,vdcop,Pijmc,Pgcm,Pgm,Hijmc,dPijmc.eta,dlogPijmc.eta,dPgcm.theta, dlogPgcm.theta,dPgm.theta,dHijmc.eta,dHijmc.theta,vx)
                                    names(matA.comp0)<-nam;
                                    
#                                    mat.dq1<-fpol1(matA.comp0,c("outc",varz,gmname,gcname,"Pgm","dPgm.theta","Pgcm","dlogPijmc.eta", "dlogPgcm.theta", colnames(vx)),"vdcop","nijmc")
                                    mat.dq1<-fpol1(matA.comp0,c("outc",varz,gmname,gcname,"Pgm"),"vdcop","nijmc", c(setdiff(colnames(vx),c(varz,gmname,gcname)),"dPgm.theta","Pgcm","dlogPijmc.eta", "dlogPgcm.theta"))
                                    # compte les modalite i,j,m,c
                                    # premier terme 
                                    dq1.beta<-apply(mat.dq1[,colnames(vx)]*mat.dq1[,"dlogPijmc.eta"]*mat.dq1[,"nijmc"],2,sum)
                                    dq1.theta<-sum(mat.dq1[,"dlogPgcm.theta"]*mat.dq1[,"nijmc"])
                                    dq1 = c(dq1.beta,dq1.theta)
                                                                                
                                    # 2.2-calcule de q2
                                    mat.Hijmc1<-fpol1(matA.comp0,c("outc",varz,gmname,gcname,"Pgm"), "vdcop","nijmc",c(setdiff(colnames(vx),c(varz,gmname,gcname)),"Hijmc","dHijmc.eta","dHijmc.theta","dPgm.theta"))
                                    mat.Hijmc<-mat.Hijmc1[mat.Hijmc1[,"Hijmc"]!=0,]
                                    # Somme sur c de Hijmc
                                    mat.Hijm<-fpol1(mat.Hijmc,c("outc",varz,gmname,"Pgm"),"Hijmc","Hijm","dPgm.theta")
                                    dHijmc.names = paste("dHijmc.",colnames(vx),sep="") 
                                    # Multiplication des termes de derivee a sommer par la matrix de design    
                                    mat.Hijmc[,dHijmc.names] = mat.Hijmc[,colnames(vx)]*mat.Hijmc$dHijmc.eta                                    
                                    # Somme sur c de dHijmc.beta
                                    dHijm.names = paste("dHijm.",colnames(vx),sep="") 
                                    mat.dHijm.beta<-fpolm(mat.Hijmc,c("outc",varz,gmname,"Pgm"),dHijmc.names,dHijm.names)
                                    # Somme sur c de dHijmc.theta
                                    mat.dHijm.theta<-fpol1(mat.Hijmc,c("outc",varz,gmname,"Pgm"),"dHijmc.theta","dHijm.theta")
                                    # nv ** 6
                                    matHijm.uim<-function(ma.u){
                                                rr<-as.matrix(mat.Hijm[,c("outc",gmname)])
                                                Uim<-(ifelse(rr[,1]==0 & rr[,2]==0,ma.u[1,1],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==1,ma.u[1,2],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==2,ma.u[1,3],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==0,ma.u[2,1],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==1,ma.u[2,2],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==2,ma.u[2,3],0))
                                                Fijm<-mat.Hijm$Hijm*Uim
                                                # Attention! Dans le data.frame dFijm.beta, les variables gardent les noms dHijm.names
                                                dFijm.beta<-mat.dHijm.beta[,dHijm.names]*Uim
                                                # On change les noms
                                                dFijm.names = paste("dFijm.",colnames(vx),sep="")                                               
                                                dFijm.theta<-mat.dHijm.theta$dHijm.theta*Uim                                                
                                                tmp = data.frame(mat.Hijm,mat.dHijm.beta[,dHijm.names],mat.dHijm.theta$dHijm.theta, Fijm,dFijm.beta,dFijm.theta)
                                                names(tmp) = c(names(mat.Hijm),dHijm.names,"dHijm.theta","Fijm",dFijm.names,"dFijm.theta")
                                                return(tmp)
                                                }
                                    tab.fijm<-matHijm.uim(ma.u)
                                    
                                    # nv 
                                    # Ici on garde les dHijm.names, dHijm.theta et les variables x pour besoins futurs                                 
                                    dFijm.names = paste("dFijm.",colnames(vx),sep="")                                    
                                    dFjm.names = paste("dFjm.",colnames(vx),sep="")                                    
                                    mat.fjm<-fpolm(tab.fijm,c(varz,gmname,"Pgm"),c("Fijm",dFijm.names,"dFijm.theta"), c("Fjm",dFjm.names,"dFjm.theta"),"dPgm.theta")
                                      
                                    tab.cjm<-merge(mat.cjm,mat.fjm,by=c(varz,gmname))
                                    qjm<-tab.cjm$Cjm/tab.cjm$Fjm
                                    # termes de derivee a sommer sur j pour obtenir dNm.beta et dNm.theta
                                    djm.beta<-(-tab.cjm$Cjm/tab.cjm$Fjm^2)*tab.cjm[,dFjm.names]                                    
                                    djm.theta<-(-tab.cjm$Cjm/tab.cjm$Fjm^2)*tab.cjm$dFjm.theta
                               
                                    tab.qjm<-data.frame(tab.cjm,qjm,djm.beta,djm.theta)
                                    djm.names = paste("djm.",colnames(vx),sep="")                                    
                                    names(tab.qjm) = c(names(tab.cjm),"qjm",djm.names,"djm.theta")
                                   
                                    # nv
                                    mat.Nm<-fpol1(tab.qjm,c(gmname,"Pgm"),"qjm","Nm")                                    
                                    mat.dNm.theta<-fpol1(tab.qjm,c(gmname,"Pgm"),"djm.theta","dNm.theta","dPgm.theta")
                                    dNm.names = paste("dNm.",colnames(vx),sep="")                                    
                                    mat.dNm.beta<-fpolm(tab.qjm,c(gmname,"Pgm"),djm.names,dNm.names)
                                                                        
                                    dlogNm.beta = mat.dNm.beta[,dNm.names]/mat.Nm$Nm
									mat.dlogNm.beta = data.frame(mat.dNm.beta,dlogNm.beta)
									dlogNm.names = paste("dlogNm.",colnames(vx),sep="")									
									names(mat.dlogNm.beta) = c(names(mat.dNm.beta),dlogNm.names)
									tab.dlogNm.beta = merge(tab.qjm,mat.dlogNm.beta,by=gmname)
									sumdlogNm.beta = apply(tab.dlogNm.beta[,dlogNm.names]*tab.dlogNm.beta$Cjm,2,sum)																		
                                    dlogFjm.beta = apply((tab.qjm[,dFjm.names]/tab.qjm$Fjm)*tab.qjm$Cjm,2,sum)
                                    
                                    dlogNm.theta = mat.dNm.theta$dNm.theta/mat.Nm$Nm                                    
									mat.dlogNm.theta = data.frame(mat.dNm.theta,dlogNm.theta)
									tab.dlogNm.theta = merge(tab.qjm,mat.dlogNm.theta,by=c(gmname))
									sumdlogNm.theta = sum(tab.dlogNm.theta$dlogNm.theta*tab.dlogNm.theta$Cjm)
																										                                    
                                    dlogFjm.theta = sum((tab.qjm$dFjm.theta/tab.qjm$Fjm)*tab.qjm$Cjm)
                                    
                                    dlogPgm.theta = sum((tab.cjm[,"dPgm.theta"]/tab.cjm[,"Pgm"])*tab.cjm[,"Cjm"])

									dq2.beta = -dlogFjm.beta - sumdlogNm.beta 
																		
									dq2.theta = dlogPgm.theta - dlogFjm.theta - sumdlogNm.theta
																		
									dq2 = c(dq2.beta,dq2.theta)
																
                                    # calcul de q3
                                    tab.Zijm<-merge(tab.fijm,tab.qjm,by=c(varz,gmname))
                                    Zijm<-tab.Zijm$Hijm*tab.Zijm$qjm
                                    tab.Zijm$Zijm<-Zijm;
                                    tab.Zijm$Fijm<-NULL
                                    tab.Zijm$Pgm.y<-NULL
                                    tab.Zijm$qjm<-NULL                                                                       
                                    
                                    dg312j.beta = (tab.Zijm$Fjm*tab.Zijm[,dHijm.names] - tab.Zijm[,dFjm.names]*tab.Zijm$Hijm)*tab.Zijm$Cjm/tab.Zijm$Fjm^2
									# On renome les variables du data.frame dg312j.beta
                                    names(dg312j.beta) = paste("dg312j.",colnames(vx),sep="")                                                                        
                                    dg312j.theta = (tab.Zijm$Fjm*tab.Zijm$dHijm.theta - tab.Zijm$dFjm.theta*tab.Zijm$Hijm)*tab.Zijm$Cjm/tab.Zijm$Fjm^2
                                    
                                    tab.dg312 = data.frame(tab.Zijm,dg312j.beta,dg312j.theta)
                                    names(tab.dg312) = c(names(tab.Zijm),names(dg312j.beta),"dg312j.theta")
                                    
                                    # nv 
                        			mat.Zim<-fpol1(tab.Zijm,c("outc",gmname),"Zijm","Zim")
                        			dg312.names = paste("dg312.",colnames(vx),sep="")
                        			# Attention Ici, pour une raison inexpliquee, il faut que c(gmname,"outc") soient dans cet ordre pour un resultat correct                                                                        
                        			mat.dg312<-fpolm(tab.dg312,c(gmname,"outc"),c(names(dg312j.beta),"dg312j.theta"),c(dg312.names,"dg312.theta"))
                        			
			
                                    Rm<-as.numeric(mat.Nm[,2])/as.numeric(mat.Nm$Nm)
                                    mat.Rm<-data.frame(mat.Nm,Rm)
                        			tab.wim<-merge(mat.Rm,mat.Zim,by=c(gmname))
                        			wim<-(tab.wim$Rm)*(tab.wim$Zim);tab.wim$wim<-wim;
                        			tab.wim$Pgm<-NULL;tab.wim$Nm<-NULL;
                        			tab.wim$Pgm.x<-NULL

                                    # Termes 311
                                    dg311.beta = -mat.dNm.beta$Pgm*mat.dNm.beta[,dNm.names]/mat.Nm$Nm^2
                                    names(dg311.beta) = paste("dg311.",colnames(vx),sep="")                                     
                                    dg311.theta = (mat.dNm.theta$dPgm.theta*mat.Nm$Nm - mat.dNm.theta$Pgm*mat.dNm.theta$dNm.theta)/mat.Nm$Nm^2
                                    # Ici, prendre mat.dNm.beta ou mat.dNm.theta ne devrait rien changer. Voir s'il faut fusionner ces 2 mat                                    
                                    mat.dg311 = data.frame(mat.dNm.beta,dg311.beta,dg311.theta)
                                    names(mat.dg311) = c(names(mat.dNm.beta),names(dg311.beta),"dg311.theta") 
                                    
                                    tmp = merge(tab.wim,mat.dg312,by=c("outc",gmname))                                     
                                    mat.dg31 = merge(tmp,mat.dg311,by=gmname)
                                    # Attention En sommant les variables dg311.beta et dg312.names, le resultat garde les noms dg311.beta                                  
                                    dg31.beta = mat.dg31[,names(dg311.beta)]*mat.dg31$Zim + mat.dg31$Rm*mat.dg31[,dg312.names]
                                    dg31.theta = mat.dg31$dg311.theta*mat.dg31$Zim + mat.dg31$Rm*mat.dg31$dg312.theta                                    

                                    dq3.beta = d[1]*apply(dg31.beta[tab.wim$outc==0,names(dg311.beta)],2,sum)/sum(tab.wim[tab.wim$outc==0,"wim"]) + d[2]*apply(dg31.beta[tab.wim$outc==1,names(dg311.beta)],2,sum)/sum(tab.wim[tab.wim$outc==1,"wim"])
                                    dq3.theta = d[1]*sum(dg31.theta[tab.wim$outc==0])/sum(tab.wim[tab.wim$outc==0,"wim"]) + d[2]*sum(dg31.theta[tab.wim$outc==1])/sum(tab.wim[tab.wim$outc==1,"wim"])
                                    dq3 = c(dq3.beta,dq3.theta)
                                                                                                                                                                                                                                                                                   
                        			return(dq1+dq2+dq3)    
                                    }                      
                   return(gradient_prof)
                    }
