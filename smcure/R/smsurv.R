smsurv <-
function(Time,Status,X,beta,w,model){    
    death_point <- sort(unique(subset(Time, Status==1)))
	if(model=='ph') coxexp <- exp((beta)%*%t(X[,-1]))  
    lambda <- numeric()
    event <- numeric()
      for(i in 1: length(death_point)){
       event[i] <- sum(Status*as.numeric(Time==death_point[i]))
                 if(model=='ph')  temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
       		if(model=='aft')  temp <- sum(as.numeric(Time>=death_point[i])*w)
                  temp1 <- event[i]
       lambda[i] <- temp1/temp
        }
    HHazard <- numeric()
    for(i in 1:length(Time)){
        HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
        if(Time[i]>max(death_point))HHazard[i] <- Inf
        if(Time[i]<min(death_point))HHazard[i] <- 0
        }
   survival <- exp(-HHazard)
   list(survival=survival)
}

