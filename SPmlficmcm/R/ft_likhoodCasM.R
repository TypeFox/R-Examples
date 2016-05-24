ft_likhoodCasM <-
function(fl,data1,N,gmname,gcname,yname,vecma.u,HW=TRUE)
                   {
                    # Arguments specifiques  la fonction
                    # d: vecteur des proportion de cas et de temoins echantillonnes 
                    # gmname: nom de la variable contenant le genotype de la mere
                    # gcname: nom de la variable contenant le genotype de enfant
                    # varz: nom des vriable de environement c("","")
                    # ma.u est la matrice des parametre u  
                    
                    # creation de deux table une avec tous les donnees complet de c et autre avec les seules c observables
                      lstdat<-fctcd(data1,gcname,yname)
                      datMod<-lstdat$datdmcp # data des donnees complets (tout ijm * 3)
                      datNo<-lstdat$datnmv # data complet c'est-a-dire pour chaque patient nous avons ijmc
                      datMv<-lstdat$datdm  # data complete c'est a dire pour chaque coupe ijm on donne toutes les combi de c
                      datrmv<-lstdat$datmv # data de donnee manquantes 
                      
                    # Nom des variables
                      varz00<-all.vars(fl)
                      # # noms des labels 
                      varz0<-all.vars(fl)[-1];varz<-varz0[-which(varz0%in%c(gmname,gcname))]
                      # re<-rep(1,dim(datrmv)[1])
                      # datrmv1<-data.frame(datrmv,re)
                      # datrmv2<-fpol1(datrmv1,c(yname,varz,gmname),"re","nijm")
                      # vec<-c(0,1,2);#datrmvc<-datrmv2
                      # gcamv<-rep(vec,dim(datrmv2)[1])
                      # datrmv3<-datrmv2[rep(1:nrow(datrmv2),rep(3,nrow(datrmv2))),]
                      # datrmv3[gcname]<-gcamv
                      
                    # A partie des donnees complet  
                    # A.1-Generer la data.frame du modele
                      cl<- model.frame(fl, data = datMod)
                      vx <- model.matrix(fl,data=datMod)
                    # extarction de la variable reponse
                      outc<-model.extract(cl,"response")
                    # Extraction du vecteur de genotypes de la mere
                      gm <- vx[,gmname]
                    # Extraction du vecteur de genotypes de l enfant
                      gc <- vx[,gcname]
                      vdcop<-datMod[,"vdcop"]
                    
                    # B partie des donnees imcomplet  
                    # B.1-Generer le "frame" du modele
#                      clb<- model.frame(fl, data = datrmv3)
#                      vxb <- model.matrix(fl,data=datrmv3)
                      clb<- model.frame(fl, data = datMv)
                      vxb <- model.matrix(fl,data=datMv)
                    # extarction de la variable reponse
                      outcb<-model.extract(clb,"response")
                    # Extraction du vecteur de genotypes de la mere
                      gmb <- vxb[,gmname]
                    # Extraction du vecteur de genotypes de lenfant
                      gcb <- vxb[,gcname]
                      vdcopb<-datMv[,"vdcop"]
                      idm = datMv["id"]                      
                    
                    # C partite commune 
                    # les valeurs possibles du geotype de la mere
                      gm1<-gm[is.na(gm)!=TRUE]
                      frq<-unique(gm1);np1<-length(frq)
                      noutc=as.character(terms.formula(fl)[[2]])
                    # selon les cas il nous donne la fonction indique a utilise pour le genotype de la mere
                      #indfg<-IndF3(gm)
                      ppd<-length(vecma.u);pp<-ppd/2;uu<-pp+1
                      ma.u<-rbind(vecma.u[1:pp],vecma.u[uu:ppd])

					# Donnees de tous les sujets                      
                      vxt = rbind(vx,vxb)
                    # Extraction du vecteur de genotypes de la mere
                      gmt <- vxt[,gmname]
                    # Extraction du vecteur de genotypes de lenfant
                      gct <- vxt[,gcname]
                    # Vecteur de reponse
                      outct = c(outc,outcb)
                      
                    # 2-Construction du systeme non lineaire =======================================
                    # A-2.1-creation de la table A
                        matd<-cbind(outc,vx)
                        np<-dim(vx)[2]
                    
                    # B-2.1-creation de la table B
                        matdb<-cbind(outcb,vxb)
                        
                        d<-vector()
                        d[1]<-N[1]-sum(data1[,noutc]==0)
                        d[2]<-N[2]-sum(data1[,noutc]==1)
                                               
                   # A- construction des indices
                        or<-rep(1,dim(data1)[1])
                        data_ori<-data.frame(data1,or)
                        
                        #construction de Cjm 
                         mat.cjm<-fpol1(data_ori,c(varz,gmname),"or","Cjm")  

                    ## calcule de hijmc
                    ## la fonction de vraisemblence 
                    liklihood_prof<-function(parms){
                                    beta.start<-parms[1:np];n1<-np+1;theta.start<-parms[n1:length(parms)]
                                    
                                    # A-nv 
                                    Pijmc<-((1/(exp((-1)*vx%*%beta.start)+1))^outc)*((1-1/(exp((-1)*(vx%*%beta.start))+1))^(1-outc))
                                    
                                    # B-nv 
                                    Pijmcb<-((1/(exp((-1)*vxb%*%beta.start)+1))^outcb)*((1-1/(exp((-1)*(vxb%*%beta.start))+1))^(1-outcb))
                                    
                                    # construction de la distribution conditionnelle du genotype de l enfant sachant la mere et celle de la mere selon que nous sommes sous HW ou non
                                    # A-nv
                                    # prob conditionel de l enf sachant la mere Pc/m matd
                                    Pgcm<-Prgcm_HW1(matd[,c(gmname,gcname)],theta.start)
                                    
                                    # B-nv
                                    # prob conditionel de l enf sachant la mere Pc/m matd
                                    Pgcmb<-Prgcm_HW1(matdb[,c(gmname,gcname)],theta.start)
                                    
                                    # A-calcul de la fonction hijmc 
                                    Hijmc<-Pijmc*Pgcm
                                    nva<-vx[,varz]
                                    
                                    # B-calcul de la fonction hijmc 
                                    Hijmcb<-Pijmcb*Pgcmb
                                    nvab<-vxb[,varz]
                                    
                                    # A-calcul du genotype 
                                    Pgm<-Prgm_HW1(matd[,gmname],theta.start)
                                    
                                    # B-calcul du genotype 
                                    Pgmb<-Prgm_HW1(matdb[,gmname],theta.start)
                                    
                                    # data.frame A
                                    nam<-c("outc",varz,gmname,gcname,"vdcop","Pijmc","Pgcm","Pgm","Hijmc")
                                    mat.Hijmc<-data.frame(outc,nva,gm,gc,vdcop,Pijmc,Pgcm,Pgm,Hijmc)
                                    names(mat.Hijmc)<-nam;
                                    matHijmc.nijmc<-fpol1(mat.Hijmc,c("outc",varz,gmname,gcname,"Pgm","Hijmc"),"vdcop","nijmc") 
                                    matHijmc.Nijmc<-matHijmc.nijmc[matHijmc.nijmc[,"Hijmc"]!=0,]
                                    
                                    #B-data.frame 
                                    namb<-c("idm","outc",varz,gmname,gcname,"Pijmc","Pgcm","Pgm","Hijmc")
                                    mat.Hijmcb<-data.frame(idm,outcb,nvab,gmb,gcb,Pijmcb,Pgcmb,Pgmb,Hijmcb)
                                    names(mat.Hijmcb)<-namb;
#                                    mat.hijmcb1<-mat.Hijmcb[mat.Hijmcb[,"Hijmc"]!=0,]
                                    mat.hijmcb1<-mat.Hijmcb[mat.Hijmcb[,"Hijmc"]!=0&vdcopb==1,]
#                                    mat.hijmb<-fpol1(mat.hijmcb1,c("outc",varz,gmname,"Pgm"),"Hijmc","Hijm")
                                    mat.hijmb<-fpol1(mat.hijmcb1,c("idm","outc",varz,gmname,"Pgm"),"Hijmc","Hijm")
#                                    mat.hijmd<-merge(datrmv2,mat.hijmb,by=c("outc",varz,gmname))
                                    
                                    #C-data.frame
                                    vdcopt = c(vdcop,vdcopb)
                                    nvat<-vxt[,varz]
                                    Pijmct = c(Pijmc,Pijmcb)                                    
                                    Pgcmt = c(Pgcm,Pgcmb)                                    
                                    Pgmt = c(Pgm,Pgmb)
                                    Hijmct = c(Hijmc,Hijmcb)                                    
                                    mat.Hijmct<-data.frame(outct,nvat,gmt,gct,vdcopt,Pijmct,Pgcmt,Pgmt,Hijmct)
                                    names(mat.Hijmct)<-nam;
                                    # Ici on ne se sert pas de nijmc. Tout ce que ceci fait, c'est enlever
                                    # les observations redondantes de mat.Hijmct
                                    matHijmct.nijmc<-fpol1(mat.Hijmct,c("outc",varz,gmname,gcname,"Pgm","Hijmc"),"vdcop","nijmc")
                                    matHijmct.Nijmc<-matHijmct.nijmc[matHijmct.nijmc[,"Hijmc"]!=0,]
                                                                        
                                    # A-compte les modalite i,j,m,c
                                    # premier terme 
                                    q1<-sum(log(matHijmc.Nijmc[,"Hijmc"])*(matHijmc.Nijmc[,"nijmc"]))  
                                    
                                    # B-compte les modalite i,j,m,c
                                    # premier terme 
#                                    q1b<-sum(log(mat.hijmd[,"Hijm"])*(mat.hijmd[,"nijm"]))
									q1b<-sum(log(mat.hijmb[,"Hijm"]))  
                                            
                                    # 2.2-calcule de q2          
                                    mat.Hijm<-fpol1(matHijmct.Nijmc,c("outc",varz,gmname,"Pgm"),"Hijmc","Hijm")
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
                        			return(q1+q1b+q2+q3)
                        			#return(q1+q2+q3)    
                                    }                      
                   return(liklihood_prof)
                    }
