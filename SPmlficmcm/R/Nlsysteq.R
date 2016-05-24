Nlsysteq <-
function(fl,data1,N,gmname,gcname,yname,beta.start=NULL,theta.start=NULL,HW=TRUE)
                   {
                    # Arguments specifiques la fonction
                    # d: vecteur des proportions des cas et des temoins echantillonnes (verifie Moliere)
                    # gmname: nom de la variable contenant le genotype de la mere
                    # gcname: nom de la variable contenant le genotype de l enfant
                    # varz: nom des vriables de l environement c("","")
                    # beta.start: valeurs initiales du vecteur de coefficients beta
                    # theta.start: valeurs initiales du vecteur de coefficients theta
                    # bateerie des testes 
                    genom<-data1[gmname];genom1<-genom[is.na(genom)!=TRUE]
                    genoen<-data1[gcname];genoen1<-genoen[is.na(genoen)!=TRUE]
                    nagm<-genom[is.na(genom)==TRUE]
                    dat<-data1[is.na(data1[gcname])!=TRUE,]
                    gt1<-dat[gmname];
                    gt2<-dat[gcname];
                    teg1<-ifelse(gt1==0 & gt2==2,1,0)
                    teg2<-ifelse(gt1==2 & gt2==0,1,0)
                    tegk<-teg1+teg2
                    # conditions
                    if(length(nagm)>0){
                                 print(gmname)
                                 stop("missing mother genotype value")
                                 } 
                    if(min(genom1)<0){
                                 print(gmname)
                                 stop("mother's genotype is negative")
                                 }
                    if(max(genom1)>2){
                                 print(gmname)
                                 stop("mother's genotype is greater than 2")
                                 }
                    if(min(genoen1)<0){
                                 print(gcname)
                                 stop("child's genotype is negative")
                                 }
                    if(max(genoen1)>2){
                                 print(gcname)
                                 stop("genoen1's genotype is greater than 2")
                                 } 
                    if(max(tegk)>0){
                                 print(gcname)
                                 stop("mother and child genotypes are not compatible")
                                 }else{
                                 
                    # creation de deux table une avec tous les donnees complet de c et l autre avec les seules c observables
                      lstdat<-fctcd(data1,gcname,yname)
                      datMod<-lstdat$datdmcp
                      datNo<-lstdat$datnmv
                    # 1-Generer le "frame" du modele
                      cl<- model.frame(fl,data=datMod)
                      vx <- model.matrix(fl,data=datMod)
                    # extarction de la variable reponse
                      outc<-model.extract(cl,"response")
                    # Extraction du vecteur de genotypes de la mere
                      gm <- vx[,gmname]
                    # Extraction du vecteur de genotypes de l enfant
                      gc <- vx[,gcname]
                    
                    # noms des labels 
                      varz0<-all.vars(fl)[-1];varz<-varz0[-which(varz0%in%c(gmname,gcname))]
                      #varz0<-attr(terms.formula(fl),"term.labels");varz<-varz0[-which(varz0%in%c(gmname,gcname))]
                      var.outc<-all.vars(fl)[1]
                    # les valeurs possibles du geotype de la mere
                      gm1<-gm[is.na(gm)!=TRUE]
                      frq<-unique(gm1);np1<-length(frq)
                                  
                    # Construction du systeme non lineaire =======================================
                      # 2 Construction de la premiere table ou nous avons la variable reponse et la matrice de design du modele
                      
                      # Matrice de designe
                        matd<-cbind(outc,vx)
                      # calcul de proba 
                        Pijmc<-((1/(exp((-1)*vx%*%beta.start)+1))^outc)*((1-1/(exp((-1)*(vx%*%beta.start))+1))^(1-outc))
                      
                      # construction de la distribution conditionnelle du genotype de l enfant sachant la mere et celle de la mere selon que nous sommes sous HW ou non
                                     
                        # prob conditionel de l enf sachant la mere
                        Pgcm<-Prgcm_HW1(matd[,c(gmname,gcname)],theta.start)
                        Hijmc<-Pijmc*Pgcm
                        nva<-vx[,varz]
                        # calcul du distribution du genotype de la mere 
                        Pgm<-Prgm_HW1(matd[,gmname],theta.start)
                        
                        vdcop<-datMod["vdcop"]
                        # construction de la table Hijmc
                        nam<-c("outc",varz,"gm","gc","vdcop","Pgm","Hijmc")
                        matA.comp0<-data.frame(outc,nva,gm,gc,vdcop,Pgm,Hijmc)
                        names(matA.comp0)<-nam;
                        
                        mat.Hijmc1<-fpol1(matA.comp0,c("outc",varz,"gm","gc","Pgm","Hijmc"),"vdcop","nijmc")
                        mat.Hijmc<-mat.Hijmc1[mat.Hijmc1[,"Hijmc"]!=0,]
                        d<-vector()
                        d[1]<-N[1]-dim(datNo[datNo[var.outc]==0,])[1]
                        d[2]<-N[2]-dim(datNo[datNo[var.outc]==1,])[1]
                                 
                    # 2.2-Creation de la tabe qui contient la fonction hijm
                    
                        # construction de tab.cjm
                        mat.cjm<-fpol1(mat.Hijmc,c(varz,"gm"),"nijmc","Cjm")
                           
                        # tab.hijm
                        mat.Hijm<-fpol1(mat.Hijmc,c("outc",varz,"gm","Pgm"),"Hijmc","Hijm")
                       
                 # 2.4 preparation de la fonction du systeme non lineaire
                       # la fonction qui a chaque vecteur de 6 parametres renvoie la valeur de chaque equation
                      qim<-function(vecMa.u){
                                            if(np1==3){ma.u<-rbind(vecMa.u[1:3],vecMa.u[4:6])
                                                      }else{ma.u<-rbind(vecMa.u[1:2],vecMa.u[3:4])}
                                           
                                           # construction des tables
                                            rr<-as.matrix(mat.Hijm[,c("outc","gm")])
                                                Uim<-(ifelse(rr[,1]==0 & rr[,2]==0,ma.u[1,1],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==1,ma.u[1,2],0)
                                                +ifelse(rr[,1]==0 & rr[,2]==2,ma.u[1,3],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==0,ma.u[2,1],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==1,ma.u[2,2],0)
                                                +ifelse(rr[,1]==1 & rr[,2]==2,ma.u[2,3],0))
                                                
                                                Fijm<-as.numeric(mat.Hijm[,"Hijm"])*Uim;
                                                mat.fijm<-data.frame(mat.Hijm,Fijm)
                                                
                                                mat.fjm1<-fpol1(mat.fijm,c(varz,"gm","Pgm"),"Fijm","Fjm")
                                                mat.fjm<-merge(mat.cjm,mat.fjm1,by=c(varz,"gm"))
                                                
                                                qjm<-as.numeric(mat.fjm[,"Cjm"])/as.numeric(mat.fjm[,"Fjm"])
                                                mat.qjm<-data.frame(mat.fjm,qjm)
                                                mat.Nm1<-fpol1(mat.qjm,c("gm","Pgm"),"qjm","Nm")
                                                Rm<-as.numeric(mat.Nm1[,"Pgm"])/as.numeric(mat.Nm1[,"Nm"])
                                                mat.Nm<-data.frame(mat.Nm1,Rm)
                                                
                                                mat.Zijm1<-merge(mat.Hijm,mat.qjm,by=c(varz,"gm","Pgm"))
                                                Zijm<-as.numeric(mat.Zijm1[,"qjm"])*as.numeric(mat.Zijm1[,"Hijm"])
                                                mat.Zijm<-data.frame(mat.Zijm1,Zijm)
                                                mat.Zim<-fpol1(mat.Zijm,c("outc","gm","Pgm"),"Zijm","Zim")
                                                
                                                mat.wim1<-merge(mat.Zim,mat.Nm,by=c("gm","Pgm"))
                                                wim<-(as.numeric(mat.wim1[,"Pgm"])*as.numeric(mat.wim1[,"Zim"]))/as.numeric(mat.wim1[,"Nm"])
                                                mat.wim<-data.frame(mat.wim1,wim)
                                                mat.wi<-fpol1(mat.wim,c("outc"),"wim","wi")
                                          # mat finale       
                                          D1<-as.numeric(mat.Nm[,"Nm"])*as.numeric(mat.wi[mat.wi["outc"]==1,]["wi"])
                                          D0<-as.numeric(mat.Nm[,"Nm"])*as.numeric(mat.wi[mat.wi["outc"]==0,]["wi"])
                                          gma<-as.numeric(mat.Nm[,"gm"])
                                          vgm<-as.numeric(mat.Nm[,"Pgm"])
                                          MatR1<-data.frame(gma,vgm,D0,D1)
                                         
                                        # ecriture des equation
                                        if(np1==3){ 
                                                   ## le systeme non lineaire 6 equation
                                                   
                                                   Eq1<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D0"])*(1-ma.u[1,1])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                   Eq2<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D0"])*(1-ma.u[1,2])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                   Eq3<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D0"])*(1-ma.u[1,3])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                   Eq4<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D1"])*(1-ma.u[2,1])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                   Eq5<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D1"])*(1-ma.u[2,2])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                   Eq6<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D1"])*(1-ma.u[2,3])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                   
                                                   Seqnl<-list(Eq1=Eq1,Eq2=Eq2,Eq3=Eq3,Eq4=Eq4,Eq5=Eq5,Eq6=Eq6)
                                                   }
                                        if(np1==2 & sum(frq)==1){
                                                                   ## le systeme non lineaire 4 equation
                                                                   Eq11<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D0"])*(1-ma.u[1,1])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                                   Eq12<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D0"])*(1-ma.u[1,2])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                                   Eq13<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D1"])*(1-ma.u[2,1])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                                   Eq14<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D1"])*(1-ma.u[2,2])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                                   Seqnl<-list(Eq1=Eq11,Eq2=Eq12,Eq3=Eq13,Eq4=Eq14)
                                                                   }
                                        if(np1==2 & sum(frq)==3){
                                                                   ## le systeme non lineaire
                                                                  
                                                                   Eq12<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D0"])*(1-ma.u[1,2])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                                   Eq13<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D0"])*(1-ma.u[1,2])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                                   Eq15<-as.numeric(MatR1[MatR1[,"gma"]==1,]["D1"])*(1-ma.u[2,1])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==1,]["vgm"])
                                                                   Eq16<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D1"])*(1-ma.u[2,2])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                                   
                                                                   Seqnl<-list(Eq1=Eq12,Eq2=Eq13,Eq3=Eq15,Eq4=Eq16)
                                                                   }
                                        if(np1==2 & sum(frq)==2){
                                                                 ## le systeme non lineaire
                                                                 
                                                                 Eq11<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D0"])*(1-ma.u[1,1])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                                 Eq13<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D0"])*(1-ma.u[1,2])-d[1]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                                 Eq14<-as.numeric(MatR1[MatR1[,"gma"]==0,]["D1"])*(1-ma.u[2,1])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==0,]["vgm"])
                                                                 Eq16<-as.numeric(MatR1[MatR1[,"gma"]==2,]["D1"])*(1-ma.u[2,2])-d[2]*as.numeric(MatR1[MatR1[,"gma"]==2,]["vgm"])
                                                                 Seqnl<-list(Eq1=Eq11,Eq2=Eq13,Eq3=Eq14,Eq4=Eq16)
                                                                 }
                                        return(Seqnl)
                                        }
                                         
                    ## 2.5-le systeme non lineair
                            if(np1==3){
                                       Model1<-function(ma.u){RR<-qim(ma.u);return(c(f1=RR$Eq1,f2=RR$Eq2,f3=RR$Eq3,f4=RR$Eq4,f5=RR$Eq5,f6=RR$Eq6))}
                                       }else{
                                             Model1<-function(ma.u){RR<-qim(ma.u);return(c(f1=RR$Eq1,f2=RR$Eq2,f3=RR$Eq3,f4=RR$Eq4))}
                                             }                    
                            return(Model1)
                    }
                    }
