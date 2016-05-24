
 
  
 
############################################################################
##              1      betafam-tests
############################################################################ 
##   source("betafam.R")  
##   example.ped<-read.table("example.ped",head=1,stringsAsFactors=F) 
##   library(glmnet)
##   test<-betafam(ped=example.ped,trace=TRUE)
##   test$elasticnet.LC.pvalue 
############################################################################ 



'betafam'<-function(ped,group.threshold=-1,fix.group.index=NULL, fix.weight=NULL,
                     mute.SMM=TRUE,trait=c("binary","qtl"),LC.test=c("LC.true","LC","sig.LC","LC.mreg","LC.lasso","LC.elasticnet"),
                     sig.LC.cutoff=0.1,true.beta=NULL,ped2multifam=FALSE,useParInRegression=FALSE,trace=FALSE){ # group.threshold=0.5
   
    

   start.time<-Sys.time()
   ################################################# 
   #(1) find nuclear families ********************** 
   #################################################
   head<-6
   M<-ncol(ped)-head 
   fam.ID<-unique(ped[,"FID"])
 

   parents.name<-NULL 
   for (i in 1:nrow(ped)) parents.name[i]<-paste(ped[i,c("FA","MO")],collapse="+")
   parents.name.u<-unique(parents.name[which(parents.name!="NA+NA")])
     
   
   
   fam.info<-NULL  #1.1 all nuclear fam
   for (i in 1:length(parents.name.u)){
        p<-parents.name.u[i] 
        fid<-ped[which(parents.name==p),"FID"]
        if (length(unique(fid))!=1) stop("multiple FID number for same parents-offspring")
        which.i<-which(parents.name==p)
        nchild<-length(which.i)
        p.split<-unlist(strsplit(p,"\\+"))
        child<-ped[which(parents.name==p),"IID"]
        par.vec<-c(as.numeric(p.split[1] %in% ped[,"IID"]),as.numeric(p.split[2] %in% ped[,"IID"]))
        if (par.vec[1]==1 & par.vec[2]==1 & nchild>=1) complete<-1 else complete<-0
        line<-c(complete,fid[1],p,par.vec,nchild,p.split,paste(child,collapse="|") )
        fam.info<-rbind(fam.info,line) 
   }
   colnames(fam.info)<-c("complete","FID","parents","nfa","nmo","nchi","FA","MO","Children")
   

   fam.complete.info<-fam.info[which(fam.info[,"complete"]==1),] #1.2 all complete nuclear fam
   dim(fam.complete.info)
   if (!ped2multifam){ #1.3 find a largest nuclear fam in each pedigree
          fid.u<-unique(fam.complete.info[,"FID"])
          select.fam<-NULL
          for (i in 1:length(fid.u)){
              which.i<-which(fam.complete.info[,"FID"]==fid.u[i]) 
              which.i.max<-which.max(as.numeric(fam.complete.info[which.i,"nchi"]))
              select.fam<-c(select.fam,which.i[which.i.max])
           }
    fam.complete.info<-fam.complete.info[select.fam,]
   }
   nfam<-nrow(fam.complete.info)
   if (nfam==0) stop("Number of complete fam is zero, so stop here!")
   maxNsib<-max(as.numeric(fam.complete.info[,"nchi"]))
   if (trace) print(paste("We got ",nfam, "  complete families",sep=""))
  
   ################################################# 
   #(2) calculate MAF and sqrt(n*f*(1-f)) as a weight
   #################################################

    
   founder.list<- which(is.na(ped[,"FA"]) & is.na(ped[,"MO"]) ) 
   ngeno<-NULL; code.freq<-NULL; weight<-NULL
   for (k in 1:M){
     ngeno[k]<-length(which(!is.na(ped[founder.list, head+k])))
     code.freq[k]<-sum(as.numeric(ped[founder.list, head+k]),na.rm=TRUE)/(2*ngeno[k])
     weight[k]<-sqrt(ngeno[k]*code.freq[k]*(1-code.freq[k]))
   }
    delete.snp<-which(code.freq<=0)
    if (length(delete.snp)>=1) {
        M<-M-length(delete.snp)
        ngeno<-ngeno[-delete.snp] 
        code.freq<-code.freq[-delete.snp]  
        weight<-weight[-delete.snp]  
        if (!is.null(true.beta)) true.beta<-true.beta[-delete.snp] 
        if (!is.null(fix.weight)) fix.weight<-fix.weight[-delete.snp] 
        if (!is.null(fix.group.index)) fix.group.index<-fix.group.index[-delete.snp]  
        if (trace) print(paste(c("We delete nopolymorphism snps: ",delete.snp),collapse="   "))
    }
   
   if (!is.null(fix.weight)) { 
       if (length(weight)==1) weight<-rep(fix.weight,M) 
       if (length(weight)>1 & length(weight)!=M) stop("fix.weigth do not match in dimension") else
       weight<-fix.weight
       if (trace) print(paste(c("Get allele weight =",weight[1],"..."),collapse=" "))
   }
 

   ################################################# 
   #(3) collapsing
   #################################################
   group.index<-NULL
   if (is.null(fix.group.index)){
      rare.snps<-which(code.freq<=group.threshold)
      Nrare<-length(rare.snps)
      if (Nrare>1) { 
          group.index[rare.snps]<- -9
          group.index[setdiff(1:M,rare.snps)]<-1:(M-Nrare) 
         } else {  
          group.index<-1:M 
        }
       if (trace) print(paste(c("In collapsing step, ",Nrare,"snps are rare (<=",group.threshold,"). Then group.index=",group.index),collapse=" "))
    } else {
    group.index<-fix.group.index 
    }
   group.ID<-sort(unique(group.index))
   Ngroup<-length(group.ID) 
 
   ################################################# 
   #(4) calculate Z[i,j,k]=T*(X_ijk-E(X_ijk)) for each marker
   #################################################

 
   Z<-array(NA,dim=c(nfam,maxNsib,Ngroup)) #all matrix is for 2 siblings.
   EX<-array(NA,dim=c(nfam,maxNsib,Ngroup)) 
   Den<-array(NA,dim=c(nfam,maxNsib,Ngroup)) 
   T<-array(NA,dim=c(nfam,maxNsib)) 
   nsib<-NULL  #i=1;j=1;i2=3;g=1;g2=1;   i=4

   for (i in 1:nfam){  # ith fam, jth individual
       
         f<-which(ped[,"IID"]==fam.complete.info[i,"FA"])
         m<-which(ped[,"IID"]==fam.complete.info[i,"MO"])
         child<-unlist(strsplit(fam.complete.info[i,"Children"],"\\|"))
         c<-which(ped[,"IID"] %in% child)
         for (j in 1:length(child)){
                i2<-c[j]
                T[i,j]<-as.numeric(ped[i2,"PHENO"]) 
                for (g in 1:Ngroup){
                           g.snp<-which(group.index==group.ID[g])
                           if (length(g.snp)==1){  
                              k<-g.snp 
                              X_ijk<-as.numeric(ped[i2,head+k])
                              moment<-call.moment(father=ped[f,head+k],mother=ped[m,head+k])
                              Z[i,j,g]<-T[i,j]*(X_ijk-moment$mean)
                              EX[i,j,g]<- moment$mean  
                              Den[i,j,g]<-T[i,j]^2*moment$var  
                             } else { 
                              Z_ijg<-NULL;X_ijg<-NULL
                              EX[i,j,g]<-0  
                              Den[i,j,g]<-0
                              for (g2 in 1:length(g.snp)){
                                k<-g.snp[g2] #kth snp in 1:M  = g2 th snp in this group
                                X_ijg[g2]<-as.numeric(ped[i2,head+k])
                                moment<-call.moment(father=ped[f,head+k],mother=ped[m,head+k])
                                Z_ijg[g2]<-weight[k]*(X_ijg[g2]-moment$mean)
                                EX[i,j,g]<-EX[i,j,g]+weight[k]* moment$mean   
                                Den[i,j,g]<-Den[i,j,g]+weight[k]^2* T[i,j]^2*moment$var 
                               }#g2
                              Z[i,j,g]<-T[i,j]*sum(Z_ijg,na.rm=T) 
                              if (trace) print(c(i,j,g,T[i,j],Z[i,j,g])) 
                            } 
                } # g 
        } #j
    }# i
 

    Zk.vec<-NULL;Zk.var<-NULL
    for (g in 1:Ngroup) {Zk.vec[g]<-sum(Z[,,g],na.rm=TRUE)
                         Zk.var[g]<-sum(Den[,,g],na.rm=TRUE)
                        }
    Zk.stat<-Zk.vec/sqrt(Zk.var) 
    single.P<-2*(1-pnorm(abs(Zk.stat)))
    minP<-min(single.P,na.rm=TRUE)
  
    if (trace) print("got Z and EX")  
   
   #################################################
   #(5) covariance matrix: sigma term 
   # D=diag( Var(Us) ), and  VA = sqrt(D)[diag(VE)^(-1/2) VE diag(VE)^(-1/2)] sqrt(D)
   #################################################   
 
   sigma<-array(NA,dim=c(Ngroup,Ngroup))  
   for (k1 in 1:Ngroup){  # k1=1; k2=2; i=3
      for (k2 in 1:Ngroup){ 
          sigma[k1,k2]<-0
          for (i in 1:nfam){ 
             sigma_i<-sum(Z[i,,k1],na.rm=T)*sum(Z[i,,k2],na.rm=T)
             sigma[k1,k2]<- sigma[k1,k2]+sigma_i
            }
      }
   }
   positive.elem<-function(x) {x[which(x<0 |x==Inf)]<-0 ; return(x)}
   D1<-diag(   positive.elem(1/sqrt(diag(sigma)))    )
   sigma_laird<-diag(sqrt(Zk.var))%*%( D1 %*%sigma%*%  D1 )%*%diag(sqrt(Zk.var))
  
   if (trace) print(paste(c("dim(sigma)=",dim(sigma_laird)),collapse=" ")) 


   #################################################
   #(6)  S_MM  
   #################################################   
   
   inv.sigma<-NA
   SMM.stat<-NA
   SMM.pvalue<-NA
   why.na<-"mute.SMM" 
   if ( !mute.SMM ){
    inv.sigma<-try(solve(sigma_laird)) 
    if (is.matrix(inv.sigma)){
        SMM.stat<-t(Zk.vec)%*% inv.sigma %*%(Zk.vec) #which is same as below.  
        SMM.pvalue<-pchisq(as.numeric(SMM.stat),df=Ngroup,lower.tail = FALSE)
        why.na<-"OK" } else { 
        if (trace) print("* * * * The error message is due to SMM test's sigular inverse VE")
        why.na<-"Singular.Sigma"
       }   
     if (trace) print(paste( "Get SMM.pvalue=:", SMM.pvalue ,sep=" ")) 
   }
     

   #################################################
   #(7.1) Get parents matrix, Tij and Xij for belowing regression.
   #################################################   
   com.par.list<-as.vector(c(fam.complete.info[,c("FA","MO")]))
   EX.par<-array(NA,dim=c(length(com.par.list),head+Ngroup)) 
   colnames(EX.par)<-c(colnames(ped)[1:head],paste("group",1:Ngroup,sep=""))
   EX.par[,1:head]<-as.matrix(ped[com.par.list,1:head])

   for (g in 1:Ngroup){
       g.snp<-which(group.index==group.ID[g])
       if (length(g.snp)==1){  
              k<-g.snp 
              EX.par[,head+g]<- ped[com.par.list,head+k] 
            } else {  
              EX.par[,head+g]<-0  
              for (g2 in 1:length(g.snp)){
                  k<-g.snp[g2] #kth snp in 1:M  = g2 th snp in this group
                  EX.par[,head+g]<-as.numeric(EX.par[,head+g])+weight[k]*as.numeric(ped[com.par.list,head+k])
                 }
       }
   EX.par[,head+g]<-as.numeric(EX.par[,head+g]) 
    }# g 
   if (trace) print(paste(c("parents matrix for lc has dimension:", dim(EX.par)),collapse=" "))
 

   #################################################
   #(7.2) linear com (LC) method by single regression 
   #################################################
  
   if ("LC" %in% LC.test){ #use the single regression for each marker
    LC.beta.hat<-NULL
    LC.prior.p<-NULL
    for (k in 1:Ngroup){  # k=1;  i=1;  j=1
      Xdata<-NULL; Ydata<-NULL;c<-0
      for (i in 1:nfam){ 
         for (j in 1:maxNsib){ 
            if (!is.na(T[i,j]) & !is.na(EX[i,j,k])){
                 c<-c+1
                 Xdata[c]<- EX[i,j,k] 
                 Ydata[c]<- T[i,j]  
                }
          }
       }#i
      if (useParInRegression) { #do not use the parents traits in regression
       
        for (i in 1:nrow(EX.par)){
         if ((trait=="binary" & EX.par[i,"PHENO"]!=0) |(trait=="qtl" & !is.na(EX.par[i,"PHENO"])) ){#*****************do not use parents
                 c<-c+1 
                 print(c)
                 Xdata[c]<-as.numeric( EX.par[i,head+k]  )
                 Ydata[c]<-as.numeric( EX.par[i,"PHENO"] )
          }  
        } #i
      } #if useParInRegression
     
      lm0<- lm(formula=Ydata~Xdata ) 
      test<-coefficients(summary(lm0))
 
      ii<-which(rownames(test)=="Xdata")
      if (length(ii)==1){
          LC.prior.p[k]<-as.numeric(test[ii,"Pr(>|t|)"])
          LC.beta.hat[k]<-as.numeric(test[ii,"Estimate"]) } else {
          LC.prior.p[k]<-99
          LC.beta.hat[k]<-NA 
         }
     } #end k for snps
     LC.beta.hat[is.na(LC.beta.hat)]<-0 #change NA, otherwise, zlc.stat=na. 
     LC.stat<-t(LC.beta.hat)%*% Zk.vec/ sqrt(t(LC.beta.hat)%*%sigma_laird %*%LC.beta.hat)  #2.32
     LC.pvalue<- 2*pnorm(abs(as.numeric( LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
     if (trace) print(paste("Get the pvalue on lc single marker regression-test: ",LC.pvalue,sep=" "))      
    
     #[##########test: single marker regression, use beta on sig.snps##
     if ("sig.LC" %in% LC.test){
          sig.LC.beta.hat<-LC.beta.hat
          sig.LC.beta.hat[which(LC.prior.p>sig.LC.cutoff)]<-0
          sig.LC.beta.hat[is.na(sig.LC.beta.hat)]<-0 #change NA, otherwise, zlc.stat=na. 
          sig.LC.stat<-t(sig.LC.beta.hat)%*% Zk.vec/ sqrt(t(sig.LC.beta.hat)%*%sigma_laird %*%sig.LC.beta.hat)  #2.32
          sig.LC.pvalue<- 2*pnorm(abs(as.numeric( sig.LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
          if (trace) print(paste("Get the pvalue on sig.LC single marker: ",sig.LC.pvalue,sep=" "))      
         } else {
          sig.LC.beta.hat<-NA 
          sig.LC.stat<-NA 
          sig.LC.pvalue<-NA 
      }
     #--------------]

  } else {
     LC.beta.hat<-NA 
     LC.stat<-NA 
     LC.pvalue<-NA 
  }
   
   #################################################
   #(7.3) LC method by  multiple regression without penalty 
   #################################################
   aic.lambda<-NA
   if ("LC.true" %in% LC.test |"LC.mreg" %in% LC.test | "LC.lasso" %in% LC.test | "LC.elasticnet" %in% LC.test){  

      ########(7.3.1 prepare the regression input)########
      Xdata<-array(NA,dim=c(nfam*(2+maxNsib),Ngroup))  
      Ydata<-NULL;c<-0
      for (i in 1:nfam){ 
         for (j in 1:maxNsib){ 
            if (!is.na(T[i,j]) & sum(as.numeric(is.na(EX[i,j,])))==0 ){
                 c<-c+1
                 Xdata[c,]<- EX[i,j,] 
                 Ydata[c]<- T[i,j]  
                }
          }
       }#i
     if (useParInRegression) { #do not use the parents traits in regression   
      for (i in 1:nrow(EX.par)){
         if ((trait=="binary" & EX.par[i,"PHENO"]==0) |(trait=="qtl" & is.na(EX.par[i,"PHENO"])) ){
                 c<-c+1 
                 Xdata[c,]<- EX.par[i,-(1:head)]  
                 Ydata[c]<- EX.par[i,"PHENO"]  
          }  
       } #i
      } #if useparent
      Xdata<-Xdata[1:c,] 
      if (trace) print("Get the input data for multiple regression.")

      #######(test 1)#######
      if ("LC.true" %in% LC.test & length(true.beta)==Ngroup){ #for debug,let true.beta=1:Ngroup
         true.beta.hat<-true.beta  
         if (t(true.beta.hat)%*% Zk.vec==0)  true.LC.stat<-0 else
             true.LC.stat<-t(true.beta.hat)%*% Zk.vec/ sqrt(t(true.beta.hat)%*%sigma_laird %*%true.beta.hat)  #2.32
         true.LC.pvalue<- 2*pnorm(abs(as.numeric( true.LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
         if (trace) print(paste("Get the pvalue on true-test: ",true.LC.pvalue,sep=" "))
         } else {
         true.beta.hat<-NA 
         true.LC.stat<-NA 
         true.LC.pvalue<-NA 
       }



      #######(test 2)#######
      if ("LC.mreg" %in% LC.test){ 
         X.input<-cbind(1,Xdata)
         dim(X.input)
         colnames(X.input)<-c("Intercept",paste("snp",1:Ngroup,sep="_")) 
         lm1<- glm.fit(x=X.input, y=Ydata, family = gaussian(), intercept = TRUE)
         mreg.beta.hat<-as.numeric(lm1$coeff[-1]) 
         
         mreg.beta.hat[is.na(mreg.beta.hat)]<-0 #change NA, otherwise, zlc.stat=na.
      
         mreg.LC.stat<-t(mreg.beta.hat)%*% Zk.vec/ sqrt(t(mreg.beta.hat)%*%sigma_laird %*%mreg.beta.hat)  #2.32
         mreg.LC.pvalue<- 2*pnorm(abs(as.numeric( mreg.LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
         if (trace) print(paste("Get the pvalue on mreg-test: ",mreg.LC.pvalue,sep=" "))
         } else {
         mreg.beta.hat<-NA 
         mreg.LC.stat<-NA 
         mreg.LC.pvalue<-NA 
       }

         
      #######(test 3)####### 
      if ("LC.lasso" %in% LC.test) { #install.packages("glmnet") 
          lm2<-glmnet(x=Xdata,y=Ydata,alpha=1)   #(1-a)beta^2+a*|beta|
          df<-lm2$df #without intercept
          RSS<- lm2$nulldev - lm2$dev * lm2$nulldev   
          AIC<-2*RSS+2*df
          AIC.i<-which.min(AIC)  
          aic.lambda<-lm2$lambda[AIC.i]
          lasso.beta.hat<-as.numeric(lm2$beta[,AIC.i])
          lasso.beta.hat[is.na(lasso.beta.hat)]<-0 #change NA, otherwise, zlc.stat=na.
      
          lasso.LC.stat<-t(lasso.beta.hat)%*% Zk.vec/ sqrt(t(lasso.beta.hat)%*%sigma_laird %*%lasso.beta.hat)  #2.32
          lasso.LC.pvalue<- 2*pnorm(abs(as.numeric( lasso.LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
          if (trace) print(paste("Get the pvalue on lasso-test: ",lasso.LC.pvalue,sep=" "))
          } else {
          lasso.beta.hat<-NA 
          lasso.LC.stat<-NA 
          lasso.LC.pvalue<-NA 
       }


      #######(test 4)####### 
      if ("LC.elasticnet" %in% LC.test) { #install.packages("MASS);library(MASS)
          lm3<-glmnet(x=Xdata,y=Ydata,alpha=0.05)   #(1-a)/2 beta^2+a*|beta|
          df<-lm2$df #without intercept
          RSS<- lm2$nulldev - lm2$dev * lm2$nulldev   
          AIC<-2*RSS+2*df
          AIC.i<-which.min(AIC)  
          aic.lambda<-lm2$lambda[AIC.i]
          elasticnet.beta.hat<-as.numeric(lm3$beta[,AIC.i])
          elasticnet.beta.hat[is.na(elasticnet.beta.hat)]<-0 #change NA, otherwise, zlc.stat=na.
      
          elasticnet.LC.stat<-t(elasticnet.beta.hat)%*% Zk.vec/ sqrt(t(elasticnet.beta.hat)%*%sigma_laird %*%elasticnet.beta.hat)  #2.32
          elasticnet.LC.pvalue<- 2*pnorm(abs(as.numeric( elasticnet.LC.stat )),mean=0,sd=1,lower.tail = FALSE)  #fbat p definition, one-side test
          if (trace) print(paste("Get the pvalue on elasticnet-test: ",elasticnet.LC.pvalue,sep=" "))
          } else {
          elasticnet.beta.hat<-NA 
          elasticnet.LC.stat<-NA 
          elasticnet.LC.pvalue<-NA 
       } 
 
   } #end all LC multiple regression tests.
 
   #################################################
   #(8) output
   #################################################   
   time.use<-Sys.time()-start.time
   if (trace) print(time.use)   
  
   ans<-list(single.P=single.P,minP=minP,Z=Zk.vec,Z.stat=Zk.stat,Zk.var=Zk.var,allele.weight=weight,group.index=group.index,Ngroup=Ngroup,
             sigma=sigma_laird,inv.sigma=inv.sigma,SMM.stat=as.numeric(SMM.stat),SMM.pvalue=SMM.pvalue,why.SMM.na=why.na,
             LC.beta=LC.beta.hat,LC.stat=as.numeric(LC.stat),LC.pvalue=LC.pvalue,
             sig.LC.beta=sig.LC.beta.hat,sig.LC.stat=as.numeric(sig.LC.stat),sig.LC.pvalue=sig.LC.pvalue,
          
             true.LC.beta=true.beta.hat,true.LC.stat=as.numeric(true.LC.stat),true.LC.pvalue=true.LC.pvalue,
             mreg.LC.beta=mreg.beta.hat,mreg.LC.stat=as.numeric(mreg.LC.stat),mreg.LC.pvalue=mreg.LC.pvalue,
             lasso.LC.beta=lasso.beta.hat,lasso.LC.stat=as.numeric(lasso.LC.stat),lasso.LC.pvalue=lasso.LC.pvalue,
             elasticnet.LC.beta=elasticnet.beta.hat,elasticnet.LC.stat=as.numeric(elasticnet.LC.stat),elasticnet.LC.pvalue=elasticnet.LC.pvalue,
             runtime=time.use,fam.info=fam.info
           )
   return(ans)
} #end  main

 

 

############################################################################
##            subfunction                       #########
############################################################################ 
  
call.moment<-function(father,mother){#father=2;mother=2; 
  
   if (is.na(father) | is.na(mother)){
     expect<-NA
     variance<-NA } else {
     father<-as.numeric(father)
     mother<-as.numeric(mother)
     if (father==0 & mother==0 ) expect<-1*0
     if (father==1 & mother==0 ) expect<-0.5*0+0.5*1 
     if (father==2 & mother==0 ) expect<-1*1
     if (father==0 & mother==1 ) expect<-0.5*0+0.5*1 
     if (father==1 & mother==1 ) expect<-0.25*0+0.5*1+0.25*2 
     if (father==2 & mother==1 ) expect<-0.5*1+0.5*2
     if (father==0 & mother==2 ) expect<-1*1
     if (father==1 & mother==2 ) expect<-0.5*1+0.5*2
     if (father==2 & mother==2 ) expect<-1*2


     if (father==0 & mother==0 ) square<-1*0^2
     if (father==1 & mother==0 ) square<-0.5*0^2+0.5*1^2 
     if (father==2 & mother==0 ) square<-1*1^2
     if (father==0 & mother==1 ) square<-0.5*0^2+0.5*1^2 
     if (father==1 & mother==1 ) square<-0.25*0^2+0.5*1^2+0.25*2^2 
     if (father==2 & mother==1 ) square<-0.5*1^2+0.5*2^2
     if (father==0 & mother==2 ) square<-1*1^2
     if (father==1 & mother==2 ) square<-0.5*1^2+0.5*2^2
     if (father==2 & mother==2 ) square<-1*2^2

     variance<-square-expect^2

   }
   ans<-list(mean=expect,var=variance)
   return(ans) 
}
 
 
############################################################################################
#      end 
############################################################################################
 
 

