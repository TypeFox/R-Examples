get.coef <-
function(event,trt,marker, study.design, rho, link, ...){  
  if(link == "risks_provided"){
    return(NULL)
  }else{
  myglm<- glm(event~trt+ marker+trt*marker,  family=binomial(link = link))
  #myglm <- glm.fit(cbind(1, trt, marker, trt*marker), event, family=binomial(link = link))
  mycoef <- myglm$coefficients

  if(substr(study.design, 1, 4) == "nest") mycoef[1] <- mycoef[1] - logit(mean(event)) + logit(rho[3]) #adjust the intercept
  else if( substr(study.design, 1,5) == "strat" ) {
   
   #need to adjust the interecept (mycoef[1]) and a1 (mycoef[2])
    p0 <- rho[3]/(rho[2]+rho[3])
    p1 <- rho[5]/(rho[4]+rho[5])
    mycoef[1] <- mycoef[1] - logit(mean(event[trt==0])) + logit(p0)
           
    mycoef[2] <- mycoef[2] + logit(mean(event[trt==0])) - logit(mean(event[trt==1]))  - logit(p0)+ logit(p1)
                                                                                                         
  }
 # browser()
  coefficients <- summary(myglm)$coefficients
  coefficients[,1] <- mycoef
  
  }
  return(coefficients)
}
