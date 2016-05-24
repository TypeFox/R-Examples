#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.RothC<- structure(
    function #Effects of moisture on decomposition rates according to the RothC model
        ###Calculates the effects of moisture (precipitation and pan evaporation) on decomposition rates according to the RothC model.
        ##references<< Coleman, K., and D. S. Jenkinson (1999), RothC-26.3 A model for the turnover of carbon in soil: 
        ##model description and windows user guide (modified 2008), 47 pp, IACR Rothamsted, Harpenden.
    (P,           ##<< A vector with monthly precipitation (mm).
     E,           ##<< A vector with same length with open pan evaporation or evapotranspiration (mm).
     S.Thick=23,  ##<< Soil thickness in cm. Default for Rothamsted is 23 cm.
     pClay=23.4,  ##<< Percent clay.
     pE=0.75,     ##<< Evaporation coefficient. If open pan evaporation is used pE=0.75. If Potential evaporation is used, pE=1.0.
     bare=FALSE   ##<< Logical. Under bare soil conditions, bare=TRUE. Dafault is set under vegetated soil.
     )
    {  
     B=ifelse(bare == FALSE,1,1.8)
     Max.TSMD=-(20+1.3*pClay-0.01*(pClay^2))*(S.Thick/23)*(1/B)
     M=P-E*pE

     Acc.TSMD=NULL
       for(i in 2:length(M)){
          Acc.TSMD[1]=ifelse(M[1] > 0, 0, M[1])
          if(Acc.TSMD[i-1]+M[i] < 0){
             Acc.TSMD[i]=Acc.TSMD[i-1]+M[i]
          }
            else(Acc.TSMD[i]=0)
   
         if(Acc.TSMD[i]<=Max.TSMD) {
            Acc.TSMD[i]=Max.TSMD
         }
       }
 
     b=ifelse(Acc.TSMD > 0.444*Max.TSMD,1,(0.2+0.8*((Max.TSMD-Acc.TSMD)/(Max.TSMD-0.444*Max.TSMD))))
     return(data.frame(Acc.TSMD,b))
     ### A data.frame with accumulated top soil moisture deficit 
     ### (Acc.TSMD) and the rate modifying factor b. 
     }
     ,
    ex=function(){
       P=c(74,59,62,51,52,57,34,55,58,56,75,71) #Monthly Precipitation (mm)
       E=c(8,10,27,49,83,99,103,91,69,34,16,8)  #Monthly open pan evaporation (mm)
       
       Rothamsted=fW.RothC(P,E)
       data.frame(month.name,P,E,0.75*E,P-0.75*E,Rothamsted)  
       # This reproduces Table 1 in the RothC documentation (Coleman and Jenkinson 1999)
       
    }        
)
