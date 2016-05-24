#------------------------------------------------------------------------------------------
# bring the lange data
#
library(gamlss.nl)
# bring the lange data
data(la) 
plot(PET60~bflow,data=la)

#source("C:/GAMLSS/New Functions/logNO.R")

#----------------------------------------------------------------------------------------
# fiiting the log Normal 

modLOGNO<- nlgamlss(y=PET60, mu.fo=~log(bflow)+log(1-(1-exp(p1))*exp(-p2/bflow)), 
               sigma.formula=~1, mu.start = c(-.99, 110), sigma.start= -2, 
               family=LOGNO, data=la)
if(abs(deviance(modLOGNO)-2293.9) > 0.1) stop("error in nl gamlss log-NO")
#-----------------------------------
# fitting the Normal 
#with original parameterization
modNO<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
                  mu.start = c(.6, 90), sigma.start= 3, data=la) 
if(abs(deviance(modNO)-2278.7) > 0.1) stop("error in nl gamlss NO")
# with different parameterizaion
modNO1<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), 
               sigma.formula=~1, mu.start = c(-.9, 90), sigma.start= 3, data=la) 
if(abs(deviance(modNO1)-2278.68) > 0.1) stop("error in nl gamlss NO")
# as an example of using function
funnl<- function(p)  bflow*(1-p[1]*exp(-p[2]/bflow))
modNO2<- nlgamlss(y=PET60, mu.fo= funnl, sigma.formula=~1, 
                  mu.start = c(.6, 90), sigma.start= 3, data=la)
if(abs(deviance(modNO2)-2278.68) > 0.1) stop("error in nl gamlss NO") 
#-----------------------------------
# fitting the Gumbel
modGU<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
                    mu.start = c(.6, 110), sigma.start= 3, family=GU, data=la) 
if(abs(deviance(modGU)-2382.9) > 0.1) stop("error in nl gamlss GU")
modGU1<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), 
         sigma.formula=~1, mu.start = c(-.6, 90), sigma.start= 3, family=GU, data=la) 
if(abs(deviance(modGU1)-2382.9) > 0.1) stop("error in nl gamlss GU")
#-----------------------------------
# fitting the reverse Gumber
modRG<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
                 mu.start = c(.6, 110), sigma.start= 3, family=RG, data=la) 
if(abs(deviance(modRG)-2249.82) > 0.1) stop("error in nl gamlss RG")
#-----------------------------------
# fitting Gamma
modGA<- nlgamlss(y=PET60,  mu.fo= ~log(bflow)+log(1-(1-exp(p1))*exp(-p2/bflow)), 
       sigma.formula=~1, mu.start = c(-.99, 90),, sigma.start= -.5, family=GA, data=la) 
if(abs(deviance(modGA)-2299.87) > 0.1) stop("error in nl gamlss GA")
#-----------------------------------
# fitting Inverse Gaussian
modIG<- nlgamlss(y=PET60,  mu.fo= ~log(bflow)+log(1-(1-exp(p1))*exp(-p2/bflow)), 
       sigma.formula=~1, mu.start = c(-.99, 90),, sigma.start= -.5, family=IG, data=la) 
if(abs(deviance(modIG)-2408.33) > 0.1) stop("error in nl gamlss IG")
#-----------------------------------
# getting the AIC
AIC(modLOGNO,modNO, modGU, modRG, modGA, modIG, k=0)
#----------------------------------------------------------------------------------------
# three  parameters distributions
modTF<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, nu.fo=~1, 
             mu.start = c(.6, 110), sigma.start= 3, nu.start=2.5 ,family=TF,  data=la) 
if(abs(deviance(modTF)-2273.5) > 0.1) stop("error in nl gamlss TF")
#-----------------------------------
modPE<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
nu.fo=~1, mu.start = c(.6, 110), sigma.start= 3, nu.start=0.6 ,family=PE, data=la) 
if(abs(deviance(modPE)-2275.89) > 0.1) stop("error in nl gamlss PE")
#-----------------------------------
modBCCG<- nlgamlss(y=PET60, mu.fo=~bflow*(1-(1-exp(p1))*exp(-p2/bflow)),sigma.formula=~1, 
nu.fo=~1, mu.start = c(-.9, 90), sigma.start= -2.3, nu.start=0, family=BCCG, data=la) 
if(abs(deviance(modBCCG)- 2293.74) > 0.1) stop("error in nl gamlss BCCG")
AIC(modTF,modPE, modBCCG,k=0)
#----------------------------------------------------------------------------------------
#  four parameters
# SEP
modSEP<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
nu.fo=~1, mu.start = c(.6, 110), sigma.start= 3, nu.start=1, tau.start=0.6, 
family=SEP, data=la) 
if(abs(deviance(modSEP)-2273.75) > 0.1) stop("error in nl gamlss SEP")
#------------------------------------
# the BCT 
modBCT<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), sigma.formula=~1,
  nu.fo=~1, mu.start=c(-.9, 90), sigma.start= -2.3, nu.start=0, tau.start=log(2.5),  
  family=BCT, data=la, control=NL.control(hessian=FALSE))  
if((deviance(modBCT)-2293.74) > 0.1) stop("error in nl gamlss BCT")
#------------------------------------
# BCPE
modBCPE<- nlgamlss(y=PET60, mu.fo=~bflow*(1-(1-exp(p1))*exp(-p2/bflow)),sigma.formula=~1, 
           mu.start = c(-.9, 90), sigma.start= -2.3, nu.start=0, tau.start=log(2.5),  
           family=BCPE, data=la)
if((deviance(modBCPE)- 2292.81) > 0.1) stop("error in nl gamlss BCPE")
#------------------------------------
# ST3
#source("C:/GAMLSS/Distributions_26_11_04/Continuous/ST3/ST3.R")
#modST3<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
#nu.fo=~1, mu.start = c(.6, 110), sigma.start= 3, nu.start=1, tau.start=0.6, 
#family=ST3, data=la) 
#if(abs(deviance(modST3)-2266.67) > 0.1) stop("error in nl gamlss ST3")
#------------------------------------
# JSU
modJSU<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), sigma.formula=~1, 
nu.fo=~1, mu.start = c(.6, 110), sigma.start= 3, nu.start=1, tau.start=0.6, 
family=JSU, data=la) 
if(abs(deviance(modJSU)-2253.45) > 0.1) stop("error in nl gamlss JSU")
AIC(modSEP,modBCPE,modJSU,  k=0) # modST3,modBCT,
#------------------------------------
#AIC(modSEP,modBCT,modBCPE,modST3,modJSU,modLOGNO,modNO, modGU, modRG, modGA, 
#    modIG,modTF,modPE, modBCCG,  k=3)
#----------------------------------------------------------------------------------------
# modelling the scale parameter
#------------------------------------
pp<-poly(la$bflow,2)
#------------------------------------
# two parameters 
#----------------------------------------------------------------------------------------
# RG
modRGSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), 
          sigma.formula=~p3+p4*pp[,1]+p5*pp[,2], 
          mu.start = c(-.9, 90), sigma.start= c(3,-1,-2), family=RG, data=la) 
if(abs(deviance(modRGSc) - 2244.87) > 0.1) stop("error in nl gamlss RG")
#------------------------------------
# LOGNO
modLOGNOSc<- nlgamlss(y=PET60, mu.fo=~log(bflow)+log(1-(1-exp(p1))*exp(-p2/bflow)), 
                sigma.formula=~p3+p4*pp[,1]+p5*pp[,2],  mu.start = c(-.99, 90), 
                sigma.start= c(3,-1,-2), 
               family=LOGNO, data=la)
if(abs(deviance(modLOGNOSc) - 2221.56) > 0.1) stop("error in nl gamlss LOGNO")
# NO
modNOSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), 
               sigma.formula=~p3+p4*pp[,1]+p5*pp[,2], mu.start = c(-.9, 90), 
                sigma.start= c(3,-1,-2),  data=la) 
if(abs(deviance(modNOSc)-2269.99) > 0.1) stop("error in nl gamlss NO")
AIC(modRGSc,modLOGNOSc,modNOSc, k=0 )
#----------------------------------------------------------------------------------------
# three parametrs
# TF
modTFSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), 
               sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2], nu.fo=~1, 
               mu.start = c(.6, 110), sigma.start= c(3,-1,-2), 
               nu.start=2.5 ,family=TF, data=la) 
if(abs(deviance(modTFSc)- 2267.42) > 0.1) stop("error in nl gamlss TF")
#----------------------------------
# PE
modPESc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), 
               sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2], nu.fo=~1, 
               mu.start = c(.6, 110), sigma.start= c(3,-1,-2), 
               nu.start=log(2) ,family=PE, data=la) 
if(abs(deviance(modPESc)- 2269.11) > 0.1) stop("error in nl gamlss TF")
#----------------------------------
# BCCG
modBCCGSc<- nlgamlss(y=PET60, mu.fo=~bflow*(1-(1-exp(p1))*exp(-p2/bflow)),
            sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2], nu.fo=~1, 
            mu.start = c(-.9, 90), sigma.start=  c(3,-1,-2), nu.start=0,   
           family=BCCG, data=la)
if((deviance(modBCCGSc)- 2221.56) > 0.1) stop("error in nl gamlss BCCGSc")
AIC(modTFSc,modPESc,modBCCGSc,  k=0)
#---------------------------------------------------------------------------------------
# four parameters 
#----------------------------------
# SEP
#modSEPSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), 
#    sigma.formula=~p3+p4*pp[,1]+p5*pp[,2],
#    nu.fo=~1, mu.start = c(.6, 110), 
#   sigma.start= c(3,-1,-2), nu.start=1, tau.start=0.6, 
#    family=SEP, data=la) 
#if(abs(deviance(modSEPSc)-2232.781) > 0.1) stop("error in nl gamlss SEP")
#----------------------------------
# BCT
modBCTSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-(1-exp(p1))*exp(-p2/bflow)), 
    sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2],
   nu.fo=~1, mu.start=c(-.9, 90), sigma.start= c(3,-1,-2), nu.start=0, 
   tau.start=log(2.5),   family=BCT, data=la, control=NL.control(hessian=FALSE)) 
if((deviance(modBCTSc)-2221.6) > 0.1) stop("error in nl gamlss BCT")
#----------------------------------
# BCPE
modBCPESc<- nlgamlss(y=PET60, mu.fo=~bflow*(1-(1-exp(p1))*exp(-p2/bflow)),
           sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2],
           mu.start = c(-.9, 90), sigma.start= c(3,-1,-2), nu.start=0, 
           tau.start=log(2.5), family=BCPE, data=la)
if((deviance(modBCPESc)- 2212.9) > 0.1) stop("error in nl gamlss BCPE")
#----------------------------------
#  JSU
#modJSUSc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), 
#           sigma.formula=~p3+p4*pp[,1]+p5*pp[,2], nu.fo=~1, 
#           mu.start = c(.6, 110), sigma.start= c(3,-1,-2), 
#           nu.start=1, tau.start=0.6, ,family=JSU, data=la) 
#if(abs(deviance(modJSUSc)-2246.1) > 0.1) stop("error in nl gamlss JSU")
#----------------------------------
# ST3 
#modST3Sc<- nlgamlss(y=PET60, mu.fo= ~bflow*(1-p1*exp(-p2/bflow)), 
#          sigma.formula=~p3+p4*pp[,1]+p5*pp[,2], nu.fo=~1, mu.start = c(.6, 110), 
#          sigma.start= c(3,-1,-2), nu.start=1, tau.start=0.6, ,family=ST3, data=la) 
#if(abs(deviance(modST3Sc)-2245.7) > 0.1) stop("error in nl gamlss ST")
#----------------------------------
#AIC(modSEPSc,modBCPESc,modJSUSc,  k=0) # modST3Sc, modBCTSc

#AIC(modSEPSc,modBCTSc,modBCPESc,modST3Sc,modJSUSc,modTFSc,modPESc,modBCCGSc,modRGSc,
#    modLOGNOSc,modNOSc,  k=3)

#AIC(modSEPSc,modBCTSc,modBCPESc,modST3Sc,modJSUSc,modTFSc,modPESc,modBCCGSc,modRGSc,
#    modLOGNOSc,modNOSc,modLOGNO,modNO, modRG, modGA,modTF,modPE, modBCCG, modSEP,
#    modBCT,modBCPE,modST3,modJSU,
#      k=3) #log(251) 

#----------------------------------------------------------------------------------------
# final model 

modFinal<- nlgamlss(y=PET60, mu.fo=~bflow*(1-(1-exp(p1))*exp(-p2/bflow)),
           sigma.formula=~ p3+p4*pp[,1]+p5*pp[,2],
           mu.start = c(-.9, 90), sigma.start= c(3,-1,-2), nu.start=0, 
           tau.start=log(2.5), nu.fix=TRUE, family=BCPE, data=la)
if((deviance(modBCPESc)- 2212.9) > 0.1) stop("error in nl gamlss BCPE")
#----------------------------------
mylogy<-log(la$PET60)
modFinal1<- nlgamlss(y=mylogy, mu.fo=~log(bflow)+log(1-p1*exp(-p2/bflow)), 
      sigma.formula=~p3+p4*pp[,1]+p5*pp[,2], mu.start = c(.6, 40), 
      sigma.start= c(3,-1,-2),nu.start=log(2.5), family=PE, data=la)
#----------------------------------
