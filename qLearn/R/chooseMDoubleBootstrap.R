chooseMDoubleBootstrap <-
function(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,alpha=0.05,boot1Num=500,boot2Num=500,...){
s2_data<-subset(completeData,completeData[,s2Indicator]==1)
model<-getModel(s2Formula,s1Formula,completeData,s2Treat,interact,s2Indicator,...)
p<-model[[3]]
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
s1_var<-names(model.frame(s1Formula,data=completeData,...))
n1 <- dim(completeData)[1]
xi <- 0.025
m <- ceiling(n1^(1-p*(xi/(1+xi))))
while (xi <= 0.15){    
    # print(xi)  # for debugging
    coverage <- 0    
    boot_est2_hard <- matrix(0,nrow=boot2Num, ncol=1)    
    for (b1 in 1:boot1Num){        
        print(b1) # for debugging
        index_b1 <- sample(1:n1,n1,replace=TRUE)
        s1_boot1 <- completeData[index_b1,]
        s2_boot1 <- subset(s1_boot1,s1_boot1[,s2Indicator]==1)
lm.stage2_boot1 <- lm(s2Formula,data=s2_boot1,...)
s1_boot1[,s1_var[1]]<-s1_boot1[,s1_var[1]]+coef(lm.stage2_boot1)[1]+as.matrix(s1_boot1[,s2_var_noint])%*%coef(lm.stage2_boot1)[s2_var_noint]+abs(coef(lm.stage2_boot1)[s2Treat]+as.matrix(s1_boot1[,interact])%*%coef(lm.stage2_boot1)[interaction])

# The next line is an ad hoc fix to the problem of NAs (sometimes the above lm model gives NA coefficients)
# ds_S1_boot1$L2_Pseudo <- ifelse(is.na(ds_S1_boot1$L2_Pseudo),-ds_S1_boot1$L2_End_QCTOT,ds_S1_boot1$L2_Pseudo)

lm.stage1_boot1 <- lm(s1Formula,data=s1_boot1,...)
boot_est1_hard <- coef(lm.stage1_boot1)[2]
            
        for (b2 in 1:boot2Num){            
            #print(b2)  # for debugging
            #print(m)
            index_b2 <- sample(1:n1,m,replace=TRUE)  #here
            s1_boot2 <- s1_boot1[index_b2,]
            s2_boot2 <- subset(s1_boot2,s1_boot2[,s2Indicator]==1)
            
            #print(summary(ds_S2_boot2))  # for debugging           
            lm.stage2_boot2 <- lm(s2Formula,data=s2_boot2,...)
s1_boot2[,s1_var[1]]<-s1_boot2[,s1_var[1]]+coef(lm.stage2_boot2)[1]+as.matrix(s1_boot2[,s2_var_noint])%*%coef(lm.stage2_boot2)[s2_var_noint]+abs(coef(lm.stage2_boot2)[s2Treat]+as.matrix(s1_boot2[,interact])%*%coef(lm.stage2_boot2)[interaction])

# The next line is an ad hoc fix to the problem of NAs (sometimes the above lm model gives NA coefficients)
# ds_S1_boot1$L2_Pseudo <- ifelse(is.na(ds_S1_boot1$L2_Pseudo),-ds_S1_boot1$L2_End_QCTOT,ds_S1_boot1$L2_Pseudo)


#print(summary(ds_S1_boot2))  # for debugging

lm.stage1_boot2 <- lm(s1Formula,data=s1_boot2,...)
      boot_est2_hard[b2,] <- coef(lm.stage1_boot2)[2]        
        }
    
        qq_hard <- quantile(boot_est2_hard, probs=c(0.05,0.95), na.rm=TRUE)
        true_param<-coef(model[[2]])[2]
        coverage <- coverage + ((2*boot_est1_hard - qq_hard[2] <= true_param) & (2*boot_est1_hard - qq_hard[1] >= true_param))   
    }
  
    coverage <- coverage/boot1Num
    # print(coverage) # here
    cat("coverage = ", coverage, "\n")
    
    if (coverage >= ((1-alpha)-qnorm(1-alpha/2)*sqrt(alpha*(1-alpha)/boot1Num))){ print('here1'); break }
    else { 
      print('here2')
        xi <- xi + 0.025
        print(xi)
        m <- ceiling(n1^(1-p*(xi/(1+xi))))
        print(m)
    } 
}
return(m)
}

