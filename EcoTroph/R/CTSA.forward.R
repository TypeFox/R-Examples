#-------------------------------------------------------------------------------
#                 CTSA : Catch Trophic Spectrum Analysis
#                        "VPA made in EcoTroph"
#-------------------------------------------------------------------------------

# selectivity function ---------------------------------------------------------
Sel.=function(x,asymptote,TL50,slope){return(asymptote/(1+exp(-slope*(x-TL50))))}
# Kinetic function -------------------------------------------------------------
K.=function(tl,temp){return(20.19*(tl**(-3.26))*exp(.041*temp))}# Gascuel et al. 2008 # all marin groups
K.acc=function(tl,temp){return(2.31*(tl**(-1.72))*exp(.053*temp))}# Gascuel et al. 2008 # fishes
# TopDown functions ------------------------------------------------------------
a13.eq.ctsa=function(compteur,kin.virg,biom,biom.virg,fish.m,TopD,FormD,range.TLpred){
  kin=kin.virg[compteur] * (1 + TopD[compteur] * (sum(biom[(compteur + range.TLpred[1]):(compteur + range.TLpred[2])])^FormD[compteur] - sum(biom.virg[(compteur+ range.TLpred[1]):(compteur + range.TLpred[2])])^FormD[compteur])/(sum(biom.virg[(compteur + range.TLpred[1]):(compteur + range.TLpred[2])])^FormD[compteur])) +  fish.m[compteur]
  return(kin)  
}
# for the highest TLs (no predators).
regPB.ctsa=function(compteur,kin.,TL,range.highTL){
  x. <- TL[(range.highTL[1]):(range.highTL[2])]
  y <- log(kin.[(range.highTL[1]):(range.highTL[2])])
  reg <- coef(lm(y ~ x.))
  reg. <- exp(reg[1] + reg[2] * TL[compteur])
  return(reg.)
}

#-------------------------------------------------------------------------------
# 1. CTSA forward starting from lowest trophic levels --------------------------
#-------------------------------------------------------------------------------

CTSA.forward=function(catch,Prod.TL1=NULL,Biomass.TL1=NULL,TE=NULL,Kin=NULL,Kin_acc=NULL,temp=NULL,
                      Selec=NULL,TL50=NULL,asymptote=NULL,slope=NULL,TopD=NULL,FormD=NULL){
  
  n.fleet=length(catch)
  # for(j in 1:n.fleet){catch.=apply(catch[[j]],1,sum)
  #                      if(j==1){catch.tot=catch.}else{catch.tot=catch.tot+catch.}
  #  }
  catch.tot=Reduce(rowSums,catch)
  #catch.tot[1]=.01
  TL=as.numeric(names(catch.tot));n.TL=length(TL);names(TL)=1:n.TL
  delta.tl=TL[2:n.TL]-TL[1:(n.TL-1)]
  
  # Calculation of TL range for the two TD equations
  tll=names(TL[TL>=2.8 & TL<=3.3])
  range.TLpred=as.numeric(c(tll[1],tll[length(tll)]))-2
  high.tl=abs(TL-5.6)
  lim.high.TL=as.numeric(names(high.tl[high.tl==min(high.tl)[1]]))
  tlll=names(TL[TL>=(TL[lim.high.TL]-0.5) & TL<=(TL[lim.high.TL])])
  range.highTL=as.numeric(c(tlll[1],tlll[length(tlll)]))
  
  if (length(TE)==1){TE=rep(TE,n.TL)}
  N_loss=-log(TE/100)
  if (is.null(Selec)){Selec=sapply(TL,Sel.,asymptote=asymptote,TL50=TL50,slope=slope)}
  if (is.null(TopD)){TopD <- rep(.2,n.TL)}else{if(length(TopD)==1){TopD=rep(TopD,n.TL)}}
  if (is.null(FormD)){FormD <- rep(.5,n.TL)}else{if(length(FormD)==1){FormD=rep(FormD,n.TL)}}
  if (is.null(Kin)){Kin=sapply(TL,K.,temp);names(Kin)=TL;Kin.virg=Kin}else{names(Kin)=TL;Kin.virg=Kin}
  if (is.null(Kin_acc)){Kin_acc=sapply(TL,K.acc,temp);names(Kin_acc)=TL;Kin_acc.virg=Kin_acc}else{names(Kin_acc)=TL;Kin_acc.virg=Kin_acc}
    
  # ------------------------------------------------------------------------------  
  # --- I. case of unexploited primary production --------------------------------
  # ------------------------------------------------------------------------------
  if(catch.tot[1]==0){
    P=rep(NA,n.TL);names(P)=TL
    
    # biomass and flow at virgin state (no fishing)
    P.virg=P
    if(!is.null(Biomass.TL1)){P.virg[1]=Biomass.TL1*Kin.virg[1]}else{P.virg[1]=Prod.TL1}
    # hypothesis the first trophic level is unexploited
    P.virg[2]=P.virg[1]*exp(-N_loss[1]*delta.tl[1])*delta.tl[2]
    for(t in 3:n.TL){P.virg[t]=P.virg[t-1]*exp(-N_loss[t]*delta.tl[t-1])}
    B.virg=P.virg/Kin.virg
    
    F_loss=rep(NA,n.TL)
    
    SC=3
    # iterations
    i=0
    while(!SC==0){
      i=i+1
      temp_Kin=Kin
      if(!is.null(Biomass.TL1)){
        #virgin state
        # flow
        P[1]=Biomass.TL1*Kin[1]
        P.virg[1]=P[1]
        P.virg[2]=P.virg[1]*exp(-N_loss[1]*delta.tl[1])*delta.tl[2]
        for(t in 3:n.TL){P.virg[t]=P.virg[t-1]*exp(-N_loss[t]*delta.tl[t-1])}
        # biomass
        B.virg=P.virg/Kin.virg
        # initialisation current flow 
        P[1]=Biomass.TL1*Kin[1]
      }else{P[1]=Prod.TL1}
      
      # equivalent to Pope approximation : Gascuel Guénette Pauly 2011 equation (7)
      
        for(t in 2:(n.TL-1)){P[t]=P[t-1]*exp(-N_loss[t-1]*delta.tl[t-1])*(delta.tl[t]/delta.tl[t-1])-catch.tot[t-1]*(delta.tl[t])*exp(-N_loss[t-1]*(delta.tl[t-1])/2)}
        # delta.tl(n) n'est pas défini TL(n+1)-TL(n), on considère que la dernière classe de est aussi de dtl=0.1
        t=n.TL
        P[t]=P[t-1]*exp(-N_loss[t-1]*delta.tl[t-1])-catch.tot[t-1]*(delta.tl[t-1])*exp(-N_loss[t-1]*(delta.tl[t-1])/2)
      
      # biomass, Gascuel Guénette Pauly 2011 equation (3)
      B=P/Kin
      
      # Fishing mortality
      Fish_mort=catch.tot/B #
      
      # Kinetic, topDown equation : Gascuel Guénette Pauly 2011 equation (5) or A13
      Kin[1]=Kin.virg[1]*(1+TopD[1]*(sum(B[TL[TL>=2&TL<=2.3]])^FormD[1]-sum(B.virg[TL[TL>=2&TL<=2.3]])^FormD[1])/(sum(B.virg[TL[TL>=2&TL<=2.3]])^FormD[1]))+Fish_mort[1]
      Kin[2:lim.high.TL]=sapply(2:lim.high.TL,a13.eq.ctsa,Kin.virg,B,B.virg,Fish_mort,TopD,FormD,range.TLpred)
      Kin[(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB.ctsa,Kin,TL,range.highTL)
      
      SC=round((sum(temp_Kin)-sum(Kin))*1E3)
      # print(i)
      # print(SC)
    }# end of iterations
    
    # ------------------------------------------------------------------------------  
    # --- II. case of exploited primary production ---------------------------------
    # ------------------------------------------------------------------------------
  }else{
    P=rep(NA,n.TL);names(P)=TL
    
    # biomass and flow at virgin state 
    P.virg=P
    # Flow.virg[1]=Flow.current[1]*exp(F_loss[1], with F_loss=catch.tot/FLow.current
    if(!is.null(Biomass.TL1)){P.virg[1]=Biomass.TL1*Kin.virg[1]*exp(catch.tot[1]/(Biomass.TL1*Kin.virg[1]))}else{P.virg[1]=Prod.TL1*exp(catch.tot[1]/Prod.TL1)}
    # flow virg computation
    P.virg[2]=P.virg[1]*exp(-N_loss[1]*delta.tl[1])*delta.tl[2]
    for(t in 3:n.TL){P.virg[t]=P.virg[t-1]*exp(-N_loss[t]*delta.tl[t-1])}
    B.virg=P.virg/Kin.virg
    
    F_loss=rep(NA,n.TL)
    
    SC=3
    # iterations
    i=0
    
    while(!SC==0){
      i=i+1
      temp_Kin=Kin
      if(!is.null(Biomass.TL1)){P[1]=Biomass.TL1*Kin[1]}else{P[1]=Prod.TL1}
      
      # equivalent to Pope approximation : Gascuel Guénette Pauly 2011 equation (7)
      
        for(t in 2:(n.TL-1)){P[t]=P[t-1]*exp(-N_loss[t-1]*delta.tl[t-1])*(delta.tl[t]/delta.tl[t-1])-catch.tot[t-1]*(delta.tl[t])*exp(-N_loss[t-1]*(delta.tl[t-1])/2)}
        t=n.TL
        P[t]=P[t-1]*exp(-N_loss[t-1]*delta.tl[t-1])-catch.tot[t-1]*(delta.tl[t-1])*exp(-N_loss[t]*(delta.tl[t-1])/2)
      
      # biomass, Gascuel Guénette Pauly 2011 equation (3)
      B=P/Kin
      
      # Fishing mortality
      Fish_mort=catch.tot/B #
      
      # Kinetic, topDown equation : Gascuel Guénette Pauly 2011 equation (5) or A13
      Kin[1]=Kin.virg[1]*(1+TopD[1]*(sum(B[TL[TL>=2&TL<=3.3]])^FormD[1]-sum(B.virg[TL[TL>=2&TL<=3.3]])^FormD[1])/(sum(B.virg[TL[TL>=2&TL<=3.3]])^FormD[1]))+Fish_mort[1]
      Kin[2:lim.high.TL]=sapply(2:lim.high.TL,a13.eq.ctsa,Kin.virg,B,B.virg,Fish_mort,TopD,FormD,range.TLpred)
      Kin[(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB.ctsa,Kin,TL,range.highTL)
      
      SC=round((sum(temp_Kin)-sum(Kin))*1E3)
      # print(i)
      # print(SC)
    }# end of iterations  
  }  
  # ------------------------------------------------------------------------------  
  # --- III. compute other variables/parameters ----------------------------------
  # ------------------------------------------------------------------------------
  
  # ET_Main Table
  Y_tot=catch.tot
  # B_acc
  B_acc=B*Selec
  # Fish_mort
  Fish_mort=catch.tot/B
  # Fish_mort_acc
  Fish_mort_acc=catch.tot/B_acc
  # F_loss & F_loss_acc
  F_loss[1]=log(P[1]/P[2]*delta.tl[2])-N_loss[1]
  for(t in 2:(n.TL-1)){F_loss[t]=log(P[t]/P[t+1])/delta.tl[t]-N_loss[t]} # Gascuel Guénette Pauly 2011 
  F_loss_acc=F_loss/Selec
  # Kin_acc
  Kin_acc[1]=Kin_acc.virg[1]*(1+TopD[1]*(sum(B[TL[TL>=2&TL<=3.3]])^FormD[1]-sum(B.virg[TL[TL>=2&TL<=3.3]])^FormD[1])/(sum(B.virg[TL[TL>=2&TL<=3.3]])^FormD[1]))+Fish_mort_acc[1]
  Kin_acc[2:lim.high.TL]=sapply(2:lim.high.TL,a13.eq.ctsa,Kin_acc.virg,B,B.virg,Fish_mort_acc,TopD,FormD,range.TLpred)
  Kin_acc[(lim.high.TL+1):n.TL]=sapply((lim.high.TL+1):n.TL,regPB.ctsa,Kin_acc,TL,range.highTL)
  # P_acc
  P_acc=B_acc*Kin_acc
  
  # Natural accessible loss rate
  N_loss_acc=rep(NA,n.TL);names(N_loss_acc)=TL # éq (A14)
  if(!F_loss_acc[1]%in%c(0,NA)){# biomasse exploitée au niveau 1 
    N_loss_acc[1]=log(P_acc[1]/P_acc[2]*delta.tl[2])-F_loss_acc[1]
  }
  for(t in 2:(n.TL-1)){N_loss_acc[t]=log(P_acc[t]/P_acc[t+1])/delta.tl[t]-F_loss_acc[t]} # Gascuel Guénette Pauly 2011 
  
  ET_Main=data.frame(cbind(B,B_acc,P,P_acc,Kin,Kin_acc,N_loss,N_loss_acc,F_loss,F_loss_acc,Fish_mort,Fish_mort_acc,Y_tot,Selec))
  x=list(ET_Main=ET_Main,Y=catch)
  class(x)='ET.CTSA.forward'
  return(x)
}

#-------------------------------------------------------------------------------
# 2. CTSA backward starting from highest trophic levels ------------------------
#-------------------------------------------------------------------------------