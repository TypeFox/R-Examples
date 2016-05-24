getModel <-
function(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,...){
# s2Treat is the name of the stage 2 treatment variable; interact is the name of the history variables which interact with s2Treat
# ... corresponds to other LM arguements

s2_data<-subset(completeData,completeData[,s2Indicator]==1)
lm.stage2<-lm(s2Formula,data=s2_data,...)
s2_var<-names(as.data.frame(model.matrix(s2Formula,data=s2_data,...)))[-1]
interaction<-paste(interact,s2Treat,sep=":")
interaction2<-paste(s2Treat,interact,sep=":")
reverse_match<-match(interaction2,s2_var)
if(sum(!is.na(reverse_match))>0){
match_num<-which(!is.na(reverse_match))
interaction[match_num]<-interaction2[match_num]
}
s2_var_noint<-s2_var[-match(s2Treat,s2_var)]
s2_var_noint<-s2_var_noint[-match(interaction,s2_var_noint)]
peudo_formula<- as.formula(paste(" ~ ", paste(s2_var_noint, collapse= "+")))
MM <- model.matrix(peudo_formula,data=completeData)
a<-match(names(as.data.frame(MM)),names(completeData))
if (sum(is.na(a))>0) {
MM <- MM[match(rownames(completeData),rownames(MM)),] 
rownames(MM) <- rownames(completeData)
extra<-MM[,is.na(a)]
completeData<-cbind(completeData,extra)
}

#Calculate the degree of nonregularity
frame<-model.frame(s2Formula,data=s2_data,...)
Y<-model.response(frame)
X2<-model.matrix(s2Formula,data=s2_data,...)
Yprime<-coef(lm.stage2)[s2Treat]+as.matrix(frame[,interact])%*%coef(lm.stage2)[interaction]
n2 <- dim(frame)[1]
Sigma2 <- n2*solve(t(X2)%*%X2)
Z2 <- diag(as.vector(Y - X2%*%coef(lm.stage2)))%*%X2%*%Sigma2/sqrt(n2-dim(X2)[2])
Cov2 <- t(Z2)%*%Z2
nonreg_term<-c(s2Treat,interaction)
nonreg_col<-match(nonreg_term,names(as.data.frame(X2)))
Sigma_stage2 <- Cov2[nonreg_col,nonreg_col]
h <- cbind(1,as.matrix(s2_data[match(rownames(X2),rownames(s2_data)),interact]))
sigma2 <- diag(h%*%Sigma_stage2%*%t(h))/n2
TS <- abs(Yprime)/sqrt(sigma2)
alpha <- 0.001            # pretest level (\nu, in paper)
cutoff <- qnorm(1 - alpha/2)
nonregularity <- (TS <= cutoff)
p <- mean(nonregularity)   # estimated value of the degree of non-regularity

s1_var<-names(model.frame(s1Formula,data=completeData,...))
peudo_data<-completeData

# Construct peusdo outcome.
peudo_data[,s1_var[1]]<-peudo_data[,s1_var[1]]+coef(lm.stage2)[1]+as.matrix(peudo_data[,s2_var_noint])%*%coef(lm.stage2)[s2_var_noint]+abs(coef(lm.stage2)[s2Treat]+as.matrix(peudo_data[,interact])%*%coef(lm.stage2)[interaction])

lm.stage1<-lm(s1Formula,data=peudo_data,...)
return(list(s2Model=lm.stage2,s1Model=lm.stage1,pHat=p))
}

