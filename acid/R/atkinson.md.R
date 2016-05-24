atkinson.md <-
  function(n,epsilon=1,dist1,dist2,theta,p0,p1,p2,dist.para.table,zero.approx){
    theta1<-theta[1:dist.para.table[which(dist.para.table$dist==dist1),3]]
    theta2<-theta[(length(theta1)+1):(length(theta1)+dist.para.table[which(dist.para.table$dist==dist2),3])]
    rdist1 <- get(paste("r",dist1,sep=""), mode = "function", envir = parent.frame())
    rdist2 <- get(paste("r",dist2,sep=""), mode = "function", envir = parent.frame())
    if(round(p0+p1+p2,digits=5)!=1) print(paste("Warning: The probabilities don't add up to unity! sum(probs)=",p0+p1+p2))
    y.sim<-rep(NA,n)
    s.sim<-rep(NA,n)
    for(i in 1:n){
      s.sim[i]<-si<-sample(0:2,size=1,prob=c(p0,p1,p2))
      if(si==1){
        if(length(theta1)==1){      y.sim[i]<-rdist1(1,theta1[1]) 
        } else if(length(theta1)==2){ y.sim[i]<-rdist1(1,theta1[1],theta1[2]) 
        } else if(length(theta1)==3){ y.sim[i]<-rdist1(1,theta1[1],theta1[2],theta1[3]) 
        } else if(length(theta1)==4){ y.sim[i]<-rdist1(1,theta1[1],theta1[2],theta1[3],theta1[4])                               
        } else print("The number of parameters cannot exceed 4.")
      }else if(si==2){ 
        if(length(theta2)==1){        y.sim[i]<-rdist2(1,theta2[1]) 
        } else if(length(theta2)==2){ y.sim[i]<-rdist2(1,theta2[1],theta2[2]) 
        } else if(length(theta2)==3){ y.sim[i]<-rdist2(1,theta2[1],theta2[2],theta2[3]) 
        } else if(length(theta2)==4){ y.sim[i]<-rdist2(1,theta2[1],theta2[2],theta2[3],theta2[4])                               
        } else print("The number of parameters cannot exceed 4.")
      }else y.sim[i]<-0
    }
    y.sim2<-y.sim
    zero.repl<-y.sim2<=0
    y.sim2[zero.repl] <-zero.approx
    A <- atkinson(y.sim2,epsilon=epsilon)
    list(AIM=A,epsilon=epsilon,y=y.sim,y2=y.sim2,zero.replace=zero.repl,stat=s.sim)
  }