ft_likhoodCas1 <-
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
                        mat.geno<-cbind(gm,gc);d<-vector()
                        d[1]<-N[1]-dim(datNo[datNo[noutc]==0,])[1]
                        d[2]<-N[2]-dim(datNo[datNo[noutc]==1,])[1]
                       
                        #construction de Cjm 
                        mat.cjm<-fpol1(datMod,c(varz,gmname),"vdcop","Cjm")
                          
                    ## calcule de hijmc
                    ## la fonction de vraisemblence 
                    liklihood_prof<-function(parms){
                                    beta.start<-parms[1:np];n1<-np+1;theta.start<-parms[n1:length(parms)]
                                    
                                    # nv
                                    eta = vx%*%beta.start
                                    Pijmc<-((1/(exp(-eta)+1))^outc)*((1-1/(exp(-eta)+1))^(1-outc))
                               
                                    # construction de la distribution conditionnelle du genotype de l enfant sachant la mere et celle de la mere selon que nous sommes sous HW ou non
                                    # nv
                                    # prob conditionel de l enf sachant la mere Pc/m matd
                                    Pgcm<-Prgcm_HW1(matd[,c(gmname,gcname)],theta.start)
                                    
                                    # calcul de la fonction hijmc 
                                    Hijmc<-Pijmc*Pgcm
                                    nva<-vx[,varz]
                                    
                                    # calcul du genotype 
                                    Pgm<-Prgm_HW1(matd[,gmname],theta.start)
                                    vdcop<-datMod["vdcop"]
                                    
                                    # data.frame A
                                    nam<-c("outc",varz,gmname,gcname,"vdcop","Pijmc","Pgcm","Pgm","Hijmc")
                                    matA.comp0<-data.frame(outc,nva,gm,gc,vdcop,Pijmc,Pgcm,Pgm,Hijmc)
                                    names(matA.comp0)<-nam;
                                    
                                    mat.Hijmc1<-fpol1(matA.comp0,c("outc",varz,gmname,gcname,"Pgm","Hijmc"),"vdcop","nijmc")
                                    mat.Hijmc<-mat.Hijmc1[mat.Hijmc1[,"Hijmc"]!=0,]
                                    
                                    # compte les modalite i,j,m,c
                                    # premier terme 
                                    q1<-sum(log(mat.Hijmc[,"Hijmc"])*(mat.Hijmc[,"nijmc"]))  
                                            
                                    # 2.2-calcule de q2
                                    mat.Hijm<-fpol1(mat.Hijmc,c("outc",varz,gmname,"Pgm"),"Hijmc","Hijm")
                                    Hijm<-mat.Hijm$Hijm
                                    # nv ** 6
                                     matHijm.uim<-function(ma.u){
                                                rr<-as.matrix(mat.Hijm[,c("outc",gmname)])
                                                Uim<-(ifelse(rr[,1]==0 & rr[,2]==0,ma.u[1,1],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==1,ma.u[1,2],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==2,ma.u[1,3],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==0,ma.u[2,1],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==1,ma.u[2,2],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==2,ma.u[2,3],0))
                                                Fijm<-Hijm*Uim;return(data.frame(mat.Hijm,Fijm))
                                                }
                                    tab.fijm<-matHijm.uim(ma.u);Fijm<-tab.fijm$Fijm
                                    
                                    # nv 
                                    mat.fjm<-fpol1(tab.fijm,c(varz,gmname,"Pgm"),"Fijm","Fjm")
                                      
                                    tab.cjm<-merge(mat.cjm,mat.fjm,by=c(varz,gmname));
                                    qjm<-tab.cjm$Cjm/tab.cjm$Fjm
                                    tab.qjm<-data.frame(tab.cjm,qjm)
                                   
                                    # nv
                                    mat.Nm<-fpol1(tab.qjm,c(gmname,"Pgm"),"qjm","Nm")
                                    
                                    Rm<-as.numeric(mat.Nm[,2])/as.numeric(mat.Nm[,3])
                                    mat.Rm<-data.frame(mat.Nm,Rm)
                                    tab.rm<-merge(tab.qjm,mat.Rm,by=c(gmname))
                                    q2<-sum(tab.rm["Cjm"]*log(tab.rm["Rm"]/tab.rm["Fjm"]))
                                    
                                    # calcul de q3
                                    tab.Zijm<-merge(tab.fijm,tab.qjm,by=c(varz,gmname))
                                    Zijm<-tab.Zijm$Hijm*tab.Zijm$qjm
                                    tab.Zijm$Zijm<-Zijm;tab.Zijm$Hijm<-NULL;tab.Zijm$Fijm<-NULL
                                    tab.Zijm$Pgm.y<-NULL;tab.Zijm$Fjm<-NULL;tab.Zijm$Cjm<-NULL
                                    tab.Zijm$qjm<-NULL
                                    
                                    # construction hijm (somme sur c)
                                    # nv 
                        			mat.Zim<-fpol1(tab.Zijm,c("outc",gmname),"Zijm","Zim")
                        			
                        			tab.wim<-merge(mat.Rm,mat.Zim,by=c(gmname))
                        			wim<-(tab.wim$Rm)*(tab.wim$Zim);tab.wim$wim<-wim;
                        			tab.wim$Pgm<-NULL;tab.wim$Nm<-NULL;tab.wim$Rm<-NULL;
                        			tab.wim$Pgm.x<-NULL;tab.wim$Zim<-NULL
                        			
                        			q3<-d[1]*log(sum(tab.wim[tab.wim$outc==0,]["wim"]))+d[2]*log(sum(tab.wim[tab.wim$outc==1,]["wim"]))
                        			return(q1+q2+q3)    
                                    }                      
                   return(liklihood_prof)
                    }
