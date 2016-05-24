cov.sel.high.sim.res<-function(object){

N<-object[[1]]$N      
Setting<-object[[1]]$Setting  
rep<-object[[1]]$rep 
Models<-object[[1]]$Models
type<-object[[1]]$type 
varnames<-object[[1]]$varnames
resmat<-matrix(NA,ncol=50,nrow=rep)
covarcol<-1:(length(varnames)-2)  
ycol<-length(varnames)-1
Tcol<-length(varnames)

##True ATE   
if(Models=="Binary"){
  if(Setting==1){beta<-0.1235}
  if(Setting==2){beta<-0.0826}
}else{
  beta<-2
}
if(Setting==1){ 
  uc1<-varnames[c(1,2,7)]
  uc2<-varnames[c(1,2,8)]
  XTS<-varnames[c(1,2,3,4,7)]
  Q1S<-Q0S<-QS<-varnames[c(1,2,7)]
  X1S<-X0S<-XS<-varnames[c(1,2,5,6,8)]
  Z1S<-Z0S<-ZS<-varnames[c(1,2,8)]
  XdS<-varnames[1:8]
}

if(Setting==2){
  uc1<-varnames[c(1,2,4,7)]
  uc2<-varnames[c(1,2,4,8)]
  XTS<-varnames[c(1,2,3,4,7)]
  Q1S<-Q0S<-QS<-varnames[c(1,2,4,7)]
  X1S<-X0S<-XS<-varnames[c(1,2,5,6,8)]
  Z1S<-Z0S<-ZS<-varnames[c(1,2,8)]
  XdS<-varnames[1:8]
  collider<-"x9"
}
  


for(i in 1:rep){
     XT<-object[[i]]$X.T
     Q<-object[[i]]$Q
     XY<-object[[i]]$X.Y
     Z<-object[[i]]$Z
     XD<-object[[i]]$X.TY
     cards<-numeric(5)
     for(j in 1:5){
                  cards[j]<-object[[i]]$cardinalities[[j]]
                  }
     betahatest<-betahatse<-numeric(6)
     for(j in 1:6){
                  betahatest[j]<-object[[i]]$est[[j]]
                  betahatse[j]<-object[[i]]$se[[j]]
                  }
     
    #Common cause criterion (Waernbaum, de Luna and Richardson)
    ##Algorithm 1
    ##Subset X.T
    if(length(XT)==0){
          XTuc<-SinXT<-XTeqS<-0
          }else{
                XTuc1<-length(which(match(uc1,XT)!="NA"))
                XTuc2<-length(which(match(uc2,XT)!="NA"))
                if(Setting==1){
                        XTuc<-ifelse(XTuc1==length(uc1) || XTuc2==length(uc2),1,0)
                        }else{ 
                              XTuc0<-ifelse(XTuc1==length(uc1) || XTuc2==length(uc2),1,0)
                              XTuc<-ifelse((XTuc0==1) && (length(which(match(collider,XT)!="NA")))==0,1,0)
                              }
                SinXT<-ifelse(length(which(match(XTS,XT)!="NA"))==length(XTS),1,0)
                XTeqS<-ifelse(SinXT==1 && length(XT)==length(XTS),1,0)
                }
    ###Subset Q
    if(length(Q)==0){
          Quc<-SinQ<-QeqS<-0
          }else{
                Quc1<-length(which(match(uc1,Q)!="NA"))
                Quc2<-length(which(match(uc2,Q)!="NA"))
                if(Setting==1){
                        Quc<-ifelse(Quc1==length(uc1) || Quc2==length(uc2),1,0)
                        }else{ 
                              Quc0<-ifelse(Quc1==length(uc1) || Quc2==length(uc2),1,0)
                              Quc<-ifelse((Quc0==1) && (length(which(match(collider,Q)!="NA")))==0,1,0)
                              }
              SinQ<-ifelse(length(which(match(QS,Q)!="NA"))==length(QS),1,0)
              QeqS<-ifelse(SinQ==1 && length(Q)==length(QS),1,0)
              }


    ##Algortithm 2
    ##Subset X.Y
    if(length(XY)==0){
          Xuc<-SinX<-XeqS<-0
          }else{
                Xuc1<-length(which(match(uc1,XY)!="NA"))
                Xuc2<-length(which(match(uc2,XY)!="NA"))
                if(Setting==1){
                        Xuc<-ifelse(Xuc1==length(uc1) || Xuc2==length(uc2),1,0)
                        }else{
                              Xuc0<-ifelse(Xuc1==length(uc1) || Xuc2==length(uc2),1,0)
                              Xuc<-ifelse((Xuc0==1) && (length(which(match(collider,XY)!="NA")))==0,1,0)
                              }
                SinX<-ifelse(length(which(match(XS,XY)!="NA"))==length(XS),1,0)
                XeqS<-ifelse(SinX==1 && length(XY)==length(XS),1,0)
                }

    ##Subset Z
    if(length(Z)==0){
          Zuc<-SinZ<-ZeqS<-0
          }else{
                Zuc1<-length(which(match(uc1,Z)!="NA"))
                Zuc2<-length(which(match(uc2,Z)!="NA"))
                if(Setting==1){
                        Zuc<-ifelse(Zuc1==length(uc1) || Zuc2==length(uc2),1,0)
                        }else{ 
                              Zuc0<-ifelse(Zuc1==length(uc1) || Zuc2==length(uc2),1,0)
                              Zuc<-ifelse((Zuc0==1) && (length(which(match(collider,Z)!="NA")))==0,1,0)
                              }
                SinZ<-ifelse(length(which(match(ZS,Z)!="NA"))==length(ZS),1,0)
                ZeqS<-ifelse(SinZ==1 && length(Z)==length(ZS),1,0)
              }
                  
    #Disjunctive cause criterion (VanderWeele and Shpitser)
    ##Subset X.D
    if(length(XD)==0){
          Xduc<-SinXd<-XdeqS<-0
          }else{             
                Xduc1<-length(which(match(uc1,XD)!="NA"))
                Xduc2<-length(which(match(uc2,XD)!="NA"))
                if(Setting==1){
                        Xduc<-ifelse(Xduc1==length(uc1) || Xduc2==length(uc2),1,0)
                        }else{ 
                              Xduc0<-ifelse(Xduc1==length(uc1) || Xduc2==length(uc2),1,0)
                              Xduc<-ifelse((Xduc0==1) && (length(which(match(collider,XD)!="NA")))==0,1,0)
                              }
                        SinXd<-ifelse(length(which(match(XdS,XD)!="NA"))==length(XdS),1,0)
                        XdeqS<-ifelse(SinXd==1 && length(XD)==length(XdS),1,0)
                }
 
                  
    ##Confidence intervals for ATE estimates                
    ciL<-betahatest-1.96*betahatse
    ciU<-betahatest+1.96*betahatse
    ## Coverage indicator
    cimat<-matrix(c(ciL,ciU),ncol=2)

    cifunc<-function(cimat){cicov<-ifelse(cimat[1]<beta && cimat[2]>beta,1,0)}
    betahat_cicov<-apply(cimat,1,cifunc) 
    
    #Matrix with results for every replication
    resmat[i,]<-c(XTuc,Quc,Xuc,Zuc,Xduc,
                  SinXT,SinQ,SinX,SinZ,SinXd,
                  XTeqS,QeqS,XeqS,ZeqS,XdeqS,
                  cards, betahatest, betahatse, betahat_cicov, ciL, ciU)
      
    }
#Summarizes the simulation results            

##Proportion of replications where the three different subset condtions are fullfilled
ss_uc<-matrix(colMeans(resmat)[1:15],nrow=5)
##Median cardinalities
med_cards<-apply(resmat[,16:20],2,median) 
##ATE bias
betahat_bias<-colMeans(resmat[,21:26])-beta
##ATE SD
betahat_sd<-apply(resmat[,21:26],2,sd) 
##ATE MSE
betahat_mse<-betahat_bias^2+apply(resmat[,21:26],2,var)
##Mean CI coverage
betahat_meancoverage<-colMeans(resmat)[33:38]
##Mean lower CI
mean_ciL<-colMeans(resmat[,39:44])
##Mean upper CI
mean_ciU<-colMeans(resmat[,45:50])
    
##List containing simulation results                              
summary_resmat<-list(Subset_selection=ss_uc,Median_cardinality=med_cards,
                         Betahat_bias=betahat_bias,Betahat_sd=betahat_sd,
                         Betahat_MSE=betahat_mse,Betahat_CI_coverage=betahat_meancoverage,
                         Betahat_mean_lower_CI=mean_ciL,Betahat_mean_upper_CI=mean_ciU)
##Table for producing LaTeX table with simulation results                              
xtab1<-data.frame(matrix(c(" ", N, rep(" ",4), " ", type, rep(" ", 4)), ncol=2))
xtab2<-data.frame(matrix(c("$X$", "$\\hat{X}_{\\rightarrow T}$", "$\\hat{Q}_{\\rightarrow T}$" ,"$ \\hat{X}_{\\rightarrow Y}$" , "$\\hat{Z}_{\\rightarrow Y}$",
                               "$\\hat{X}_{\\rightarrow T, Y}$" ), ncol=1))
if(Setting==1){
xtab3<-data.frame(rbind(c(1,1,0),summary_resmat$Subset_selection)*100)
}else{
xtab3<-data.frame(rbind(c(0,1,0),summary_resmat$Subset_selection)*100)
  
}
xtab4<-data.frame(matrix(c(100,summary_resmat$Median_cardinality,summary_resmat$Betahat_bias,
                               summary_resmat$Betahat_sd,summary_resmat$Betahat_MSE,
                               summary_resmat$Betahat_CI_coverage*100, summary_resmat$Betahat_mean_lower_CI,
                               summary_resmat$Betahat_mean_upper_CI),ncol=7))
xtab<-cbind(xtab1,xtab2, xtab3,xtab4)
digits<-c(rep(0,4),rep(1,3),0,rep(3,3),1,rep(3,2))
align<-c(rep("l",4),rep("r",10))
colnames(xtab)<-c("n","Method","$\\hat{S}$" ," $Y_t \\perp\\!\\!\\!\\perp T \\mid \\hat{S}$",
                  "$S \\subseteq \\hat{S}$","$S=\\hat{S}$",  "\\#", "Bias", "SD" ,"MSE", "CP" ,"CIL", "CIU")
xtab_res<-xtable(x=xtab,digits=digits,align=align)
                                  
return(list(resmat=resmat,summary_resmat=summary_resmat,xtable=xtab_res))
}
