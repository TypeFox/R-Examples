# 'coucou

MSE2.0<-function(X_appr=X_appr,B=B,Z=Z,X_test=X_test,Y_appr=Y_appr,Y_test=Y_test,A=A,type="lar",moyenne=T,compl_lars=T,expl_lars=T,pred_lars=T){
  #mod?le complet####
  res_lars_compl=lars(X_appr,Y_appr,type=type)
  mod_complBIC=meilleur_lars(res_lars_compl,X=X_appr,Y=Y_appr)
  BIC_compl=BIC(mod_complBIC)
  
  #mod?le explicatif####
  quiI1=which(colSums(Z)==0)
  X_appr_expl=X_appr[,quiI1]
  
  res_lars_expl=lars(x=X_appr_expl,y=Y_appr,type=type)
  mod_explBIC=meilleur_lars(res_lars_expl,X=X_appr_expl,Y=Y_appr)
  BIC_expl=BIC(mod_explBIC)
  #mod?le pr?dictif
  Ytilde_larsBIC=as.matrix(mod_explBIC$residuals)
  Xtilde_larsBIC=as.matrix(extract_Xtilde(X=X_appr,B=B))
  res_lars_pred=lars(Xtilde_larsBIC,Ytilde_larsBIC,intercept=FALSE,type=type)
  mod_inj_lars=meilleur_lars(res_lars_pred,X=Xtilde_larsBIC,Y=Ytilde_larsBIC,intercept=FALSE)
  
  
  compl_lars=mean((Y_test-predict(mod_complBIC,newdata=data.frame(X_test)))^2)
  expl_lars=mean((Y_test-predict(mod_explBIC,newdata=data.frame(X_test)))^2)
  #    mod_pred=mon_predict_pred2(mod_expl=mod_explBIC,mod_inj=mod_inj_lars,newdata=data.frame(X_test),B=B)
  #    pred_lars=mean((Y_test-mod_pred$hatY)^2)
  #    theta=mod_pred$hatA
  #    BIC_pred=calcul_BIC_theta(theta=as.matrix(c(apply(Y_test-mod_pred$hatY,2,sd),mod_pred$hatA)),X=X_appr,Y=Y_appr,k=(1+length(mod_pred$hatA)))
  pred_lars=NULL
  BIC_pred=NULL
  vraimod=mean((Y_test-cbind(rep(1,times=nrow(X_test)),X_test)%*%A)^2)
  BIC_vrai=calcul_BIC_theta(theta=as.matrix(c(apply(Y_appr-cbind(rep(1,times=nrow(X_appr)),X_appr)%*%A,2,sd),A)),X=X_appr,Y=Y_appr,k=(1+length(A[A!=0])))
  
  moyenne=mean((Y_test-mean(Y_appr))^2)
  #bilan####  
  return(list(compl.lars=compl_lars,expl.lars=expl_lars,pred.lars=pred_lars,vraimod=vraimod,moyenne=moyenne,BIC_expl=BIC_expl,BIC_compl=BIC_compl,BIC_pred=BIC_pred,BIC_vrai=BIC_vrai))
}