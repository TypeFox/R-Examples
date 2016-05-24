

    #########################
    #### function Hwyhat ####
    #########################

 ## Descripcion: 
 ##    funció intermitja que es cridada si els metodes són "ocv",..., tots menys 
 ##    el que defineix el rank l'usuari (eff.rank).
 ##    Implementada per desminuir el temps d'execució (al repetir-se en molts 
 ##    casos fins a n vegades). Evita doncs, tots els càlculs innecessaris de la 
 ##    funció principal dblm.
 ##         call HwProject:
 ##              - hat matrix Hw
 ##              - eff.rank 
 ##              - Fwplus
 ##         estimators:
 ##              - fitted values yhat
 ##              - resStand.err
 ##              - ordinary cross-validation
 ##              - generalized cross-validation
 ##              - residual degree freedom
 ##
 
 
 Hwyhat<-function(G,n,Dsqw,weights,rk,epsilon,y,y0,cvyes=FALSE,ori_weights,...){
 
    # call the HwProjector to know the Hat matrix, Fwplus and the eff.rank 
    #(that normally is the same that rk)
    HwFwplus<- HwProject(G,Dsqw,rk=rk,epsilon=epsilon,cvyes=cvyes,n) 
    Hw<- HwFwplus$Hw
    Fwplus<- HwFwplus$Fwplus
    eff.rank<-  HwFwplus$eff.rank
    used_rel.gvar<-HwFwplus$used_rel.gvar
    
    # fitted values,considering the weights,the hat matix, and the centered y (y0).
    yhat <- sum(weights*y) + Hw %*% y0
    
    # resid. stand. desviacion. use the original weights (not the percentual weights) 
    # to achive the same results that lm (see that the weights are only multiplied in the numerator). 
    if(n-1-eff.rank!=0) resStand.err<-sqrt(sum(ori_weights*(yhat-y)^2)/(n-1-eff.rank))
    else resStand.err<-0
    
    # calculing the ordinary cross-validation estimator
    ocv<-  sum(weights*((y-yhat)/(1-diag(Hw)))^2)
      
    # calculing the generalized cross-validation estimator
    gcv <-sum(weights*(yhat-y)^2)/((1-mean(diag(Hw)))^2)  
    
    # Return a list with the fitted values, the effective rank, the inverse of 
    # Fw, the hat matrix, the ocv and gcv estimators, the residual standard 
    # error, and the residual degree of freedom
    ans<-list(yhat=yhat,eff.rank=eff.rank,Fwplus=Fwplus,Hw=Hw,ocv=ocv,gcv=gcv,
      resStand.err=resStand.err,rdf=(n-1-eff.rank),used_rel.gvar=used_rel.gvar)
    
    return (ans)
 }
