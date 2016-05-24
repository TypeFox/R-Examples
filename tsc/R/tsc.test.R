tsc.test <-
function(x,y,method="DBEL",t_m=2,mc=3000) {
  #if ((method==" ") || (!exists(method)) method="DBEL"
  #Make sure method is equal to "TAS" or "TAS&Shapiro-wilk" or "DBEL"
  if ( ( method !="TAS" ) & (  method !="TAS&SW" ) & (  method !="DBEL" )) {
    print("tsc.test(): method must be equal to TAS or TAS&SW or DBEL") 
    return(NaN);
  }
  
    # Make sure x and y are numeric.
    if ( ( ! is.numeric(x) ) || ( ! is.numeric(y) ) ) {
      print("tsc.test(): input must be numeric!")
      return(NaN);
    }
    
    # Make sure x and y are vectors.
    if ( ( ! is.vector(x) ) || ( ! is.vector(y) ) ) {
      print("tsc.test(): input must be vectors!")
      return(NaN);
    }
  
    if (method=="TAS") {
      cat("\n The LRT method\n")
      cat("\n...Working on test_stat \n")
      test<-plrt(x,y)
      cat("\n...Working on p_value \n") 
      final <- list(test_stat=test$test_stat, p_value=test$p_value) 
      print(final)
    }
    if (method=="TAS&SW") {
      cat("\nThe combined test based on LRT and the Shapiro-Wilk test\n")
      cat("\n...Working on test_stat \n") 
      test<-plrt(x,y)
      cat("\n...Working on p_value \n") 
      p1<-shapiro.test(x)$p.value
      p2<-shapiro.test(y)$p.value
      p3<-min(p.adjust(c(test$p_value,p1,p2),method="bonferroni"))
      final <- list(test_stat=test$test_stat, p_value=p3)
      print(final)
    } 
  
    
    if (method=="DBEL"){
    cat("\n The DBEL ratio test\n")
    cutval<-NULL
    rm(list="cutval")
    table<-cutval
    delta<-0.1
    length1=length(x)
    length2=length(y)
    indx<-which(table[,1]==length1 & table[,2]==length2)
    #if the sample sizes exist in the table, we can directly get the p-value based on the table.
    if (length(indx)==1){
    cat("\n...Working on test_stat \n")
    cat("\nThe sample size of x and the sample size of y exist in the tabulated value, so the interpolation method based on regression method  is used to obtain the p-value. \n")
    alpha<-seq(0.01,0.1,0.01)
    log<-log(alpha/(1-alpha))
    cv<-as.numeric(as.character((c(t(table[indx,3:12])))))
    new=.C("CWrapper1", n1=length1,n2=length2,y1=x,y2=y,test_stat=as.double(1))
    model<-as.matrix(as.numeric(glm(log~cv)$coefficients ))
    cat("\n...Working on p_value \n")
    if (new$test_stat+100<min(cv))p=0.00001
    else p<-exp(model[1]+model[2]*new$test_stat)/(1+exp(model[1]+model[2]*new$test_stat))
    final <- list(test_stat=new$test_stat, p_value=p)
    print(final)
    }
    #if the sample size doesn't exist in the table, using MC or table or hybrid method or bayesian method to get the p_value.  
    if (length(indx)==0){
    #if the method isn't specified to obtain p-value, then hybrid is default to use.
    #if ((t_m==" ") || (!exists(method))(!exists(t_m))) t_m=2
    if (t_m == 1) {
      cat("\n...Working on test_stat \n")
      cat("\n...Working on p-value \n")
      cat("\nMonte Carlo method is used to obtain p-value\n")
      test<-.C("CWrapper", n1=length1,n2=length2,y1=x,y2=y,mc=mc,test_stat=as.double(1),p_value=as.double(1))
      final <- list(test_stat=test$test_stat, p_value=test$p_value)
      print(final)
    }
   
    #Interpulation and MC combined.
    if (t_m == 2 ) {
      cat("\n...Working on test_stat \n")
      cat("\n...Working on p-value \n")
      cat("\nThe interpolation method based on regression technique and tabulated critical values is used to obtain the p-value\n")
      
    #if the sample sizes don't exist in the table, interpolation will be used to get the p-value
    #if sample size less than 425, there will be 80 obs included into the regression model
     if (length1<=225 & length2<=225) {
     indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1],table[findInterval(length1,table[,1])+25,1],table[findInterval(length1,table[,1])+80,1])
     indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+2,2])
     
     n2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
     m2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
     combine<-merge(n2,m2)
     for(i in 1:2809) if ((length1+length2==combine[i,1]+combine[i,2])&((combine[i,2]<indx2[1])||(combine[i,2]>indx2[4]))&((combine[i,1]<indx1[1])||(combine[i,1]>indx1[4]))&(combine[i,1]>3 & combine[i,2]>3) & (combine[i,1]<=425 & combine[i,2]<=425) ) { length3<-combine[i,1]
                                                                                                                                                   length4<-combine[i,2]      
     } 
     #include 40 more observations to make the regression line accurate. include rule: n1+n2=m1+m2
     indx4<-c(table[findInterval(length3,table[,1])-60,1],table[findInterval(length3,table[,1]),1],table[findInterval(length3,table[,1])+40,1],table[findInterval(length3,table[,1])+80,1])
     indx5<-c(table[findInterval(length4,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2])+1,2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2])+2,2])  
     cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12],table[table$V1 == indx1[3] & table$V2 == indx2[3],3:12],table[table$V1 == indx1[4] & table$V2 == indx2[4],3:12],
                                     table[table$V1 == indx4[1] & table$V2 == indx5[1],3:12],table[table$V1 == indx4[2] & table$V2 == indx5[2],3:12],table[table$V1 == indx4[3] & table$V2 == indx5[3],3:12],table[table$V1 == indx4[4] & table$V2 == indx5[4],3:12]))))
     #sample size used in the model.
     #1/sqrt(n1)
     n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10),rep.int(1/sqrt(indx1[3]),10),rep.int(1/sqrt(indx1[4]),10),
            rep.int(1/sqrt(indx4[1]),10),rep.int(1/sqrt(indx4[2]),10),rep.int(1/sqrt(indx4[3]),10),rep.int(1/sqrt(indx4[4]),10))
     #1/sqrt(n2)
     n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10),rep.int(1/sqrt(indx2[3]),10),rep.int(1/sqrt(indx2[4]),10),
            rep.int(1/sqrt(indx5[1]),10),rep.int(1/sqrt(indx5[2]),10),rep.int(1/sqrt(indx5[3]),10),rep.int(1/sqrt(indx5[4]),10))
     ncombine1<-n11+n22
     #sqrt(n1)
     n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10),rep.int(sqrt(indx1[3]),10),rep.int(sqrt(indx1[4]),10),
            rep.int(sqrt(indx4[1]),10),rep.int(sqrt(indx4[2]),10),rep.int(sqrt(indx4[3]),10),rep.int(sqrt(indx4[4]),10))
    #sqrt(n2)
     n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10),rep.int(sqrt(indx2[3]),10),rep.int(sqrt(indx2[4]),10),
            rep.int(sqrt(indx5[1]),10),rep.int(sqrt(indx5[2]),10),rep.int(sqrt(indx5[3]),10),rep.int(sqrt(indx5[4]),10))
     ncombine2<-n33+n44
     alpha<-seq(0.01,0.1,0.01)
     log<-rep(log(alpha/(1-alpha)),8)
     #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
     model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log)$coefficients ))
     new<-.C("CWrapper1", n1=length1,n2=length2,y1=x,y2=y,test_stat=as.double(1))$test_stat
     #1/sqrt(length1)+1/sqrt(length2)
     n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
     #sqrt(length1)+sqrt(length2)
     n_sqrt<-sqrt(length1)+sqrt(length2)
     m<-exp((new-model[1]-model[2]*n_sqrt1-model[3]*n_sqrt)/model[4])
     if (new+100<min(cv)) {p=0.000001} else {p<-m/(1+m)} 
     final <- list(test_stat=new, p_value=p)
     print(final)
     }    
     #if sample size is less or equal to 450 and greater than 425, there will be 60 obs included into the regression model
     if (((length1<300 & length1> 225) & (length2<300 & length2 > 225))||((length1<225) & (length2<300 & length2 > 225))||((length2<225) & (length1<300 & length1 > 225))) {
      indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1],table[findInterval(length1,table[,1])+40,1])
      indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+1,2])
      n2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
      m2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
   
      cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12],table[table$V1 == indx1[3] & table$V2 == indx2[3],3:12]))))
      #sample size used in the model.
      #1/sqrt(n1)
      n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10),rep.int(1/sqrt(indx1[3]),10))
      #1/sqrt(n2)
      n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10),rep.int(1/sqrt(indx2[3]),10))
      ncombine1<-n11+n22
      #sqrt(n1)
      n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10),rep.int(sqrt(indx1[3]),10))
      #sqrt(n2)
      n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10),rep.int(sqrt(indx2[3]),10))
      ncombine2<-n33+n44
      alpha<-seq(0.01,0.1,0.01)
      log<-rep(log(alpha/(1-alpha)),3)
      #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
      model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log-1)$coefficients ))
      new<-.C("CWrapper1", n1=length1,n2=length2,y1=x,y2=y,test_stat=as.double(1))$test_stat
      #1/sqrt(length1)+1/sqrt(length2)
      n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
      #sqrt(length1)+sqrt(length2)
      n_sqrt<-sqrt(length1)+sqrt(length2)
      m<-exp((new-model[1]*n_sqrt1-model[2]*n_sqrt)/model[3])
      if (new+100<min(cv)) {p=0.000001} else {p<-m/(1+m)}
      final <- list(test_stat=new, p_value=p)
      print(final)
      }
      #if sample size is greater than 450, there will be 20 obs included into the regression model
      if ((length1>=300) || (length2>=300)) {
      indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1])
      indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2])
      cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12]))))
      #sample size used in the model.
      #1/sqrt(n1)
      n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10))
      #1/sqrt(n2)
      n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10))
      ncombine1<-n11+n22
      #sqrt(n1)
      n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10))
      #sqrt(n2)
      n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10))
      ncombine2<-n33+n44
      alpha<-seq(0.01,0.1,0.01)
      log<-rep(log(alpha/(1-alpha)),2)
      #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
      #summary(glm(cv~ncombine1+ncombine2+log-1))
      model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log-1)$coefficients ))
     # new<-teststatistics(x,y)
      new<-.C("CWrapper1", n1=length1,n2=length2,y1=x,y2=y,test_stat=as.double(1))$test_stat
      #1/sqrt(length1)+1/sqrt(length2)
      n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
      #sqrt(length1)+sqrt(length2)
      n_sqrt<-sqrt(length1)+sqrt(length2)
      m<-exp((new-model[1]*n_sqrt1-model[2]*n_sqrt)/model[3])
      if (new+100<min(cv)) {p=0.000001} else {p<-m/(1+m)}
      final <- list(test_stat=new, p_value=p)
      print(final)
    }
    }
    #method interpulation and MC combination
    if (t_m==3) {
      cat("\n...Working on test_stat \n")
      cat("\n...Working on p-value \n")
      cat("\nA Bayesian method is used to obtain the p-value\n")
      #step I :getting u0, sigma0
      #if sample size less than 275, there will be 80 obs included into the regression model
        if (length1<=225 & length2<=225) {  
          indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1],table[findInterval(length1,table[,1])+40,1],table[findInterval(length1,table[,1])+80,1])
          indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+2,2])
          n2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
          m2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
          combine<-merge(n2,m2)
          for(i in 1:2809) if ((length1+length2==combine[i,1]+combine[i,2])&((combine[i,2]<indx2[1])||(combine[i,2]>indx2[4]))&((combine[i,1]<indx1[1])||(combine[i,1]>indx1[4]))&(combine[i,1]>3 & combine[i,2]>3) & (combine[i,1]<=425 & combine[i,2]<=425) ) { length3<-combine[i,1]
                                                                                                                                                                                                                                                                  length4<-combine[i,2]      
          } 
          
          #include 40 more observations to make the regression line accurate. include rule: n1+n2=m1+m2
          indx4<-c(table[findInterval(length3,table[,1])-60,1],table[findInterval(length3,table[,1]),1],table[findInterval(length3,table[,1])+40,1],table[findInterval(length3,table[,1])+80,1])
          indx5<-c(table[findInterval(length4,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2])+1,2],table[findInterval(length4,table[which(table[,1]==indx1[1]),2])+2,2])  
          cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12],table[table$V1 == indx1[3] & table$V2 == indx2[3],3:12],table[table$V1 == indx1[4] & table$V2 == indx2[4],3:12],
                                          table[table$V1 == indx4[1] & table$V2 == indx5[1],3:12],table[table$V1 == indx4[2] & table$V2 == indx5[2],3:12],table[table$V1 == indx4[3] & table$V2 == indx5[3],3:12],table[table$V1 == indx4[4] & table$V2 == indx5[4],3:12]))))
          #sample size used in the model.
          #1/sqrt(n1)
          n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10),rep.int(1/sqrt(indx1[3]),10),rep.int(1/sqrt(indx1[4]),10),
                 rep.int(1/sqrt(indx4[1]),10),rep.int(1/sqrt(indx4[2]),10),rep.int(1/sqrt(indx4[3]),10),rep.int(1/sqrt(indx4[4]),10))
          #1/sqrt(n2)
          n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10),rep.int(1/sqrt(indx2[3]),10),rep.int(1/sqrt(indx2[4]),10),
                 rep.int(1/sqrt(indx5[1]),10),rep.int(1/sqrt(indx5[2]),10),rep.int(1/sqrt(indx5[3]),10),rep.int(1/sqrt(indx5[4]),10))
          ncombine1<-n11+n22
          #sqrt(n1)
          n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10),rep.int(sqrt(indx1[3]),10),rep.int(sqrt(indx1[4]),10),
                 rep.int(sqrt(indx4[1]),10),rep.int(sqrt(indx4[2]),10),rep.int(sqrt(indx4[3]),10),rep.int(sqrt(indx4[4]),10))
          #sqrt(n2)
          n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10),rep.int(sqrt(indx2[3]),10),rep.int(sqrt(indx2[4]),10),
                 rep.int(sqrt(indx5[1]),10),rep.int(sqrt(indx5[2]),10),rep.int(sqrt(indx5[3]),10),rep.int(sqrt(indx5[4]),10))
          ncombine2<-n33+n44
          alpha<-seq(0.01,0.1,0.01)
          log<-rep(log(alpha/(1-alpha)),8)
          #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
          model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log)$coefficients ))
          #1/sqrt(length1)+1/sqrt(length2)
          n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
          #sqrt(length1)+sqrt(length2)
          n_sqrt<-sqrt(length1)+sqrt(length2)
          alpha1<-0.05
          log1<-log(alpha1/(1-alpha1))
          u0<-model[1]+model[2]*n_sqrt1+model[3]*n_sqrt+model[4]*log1
       #   sigma0<-var(glm(cv~ncombine1+ncombine2+log)$residuals)
          sigma0<-var(glm(cv~ncombine1+ncombine2+log)$residuals)*sqrt(length(n11))    #check whether it should be time sqrt(n11) or just n11          
        } 
        #if sample size is less or equal to 300 and greater than 275 
        if (((length1<300 & length1> 225) & (length2<300 & length2 > 225))||((length1<225) & (length2<300 & length2 > 225))||((length2<225) & (length1<300 & length1 > 225))){
          indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1],table[findInterval(length1,table[,1])+40,1])
          indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2])+1,2])
         # n2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
         # m2<-c(2:30,35,40,45,50,55,60,70,80,90,100,120,150,170,200,225,250,275,300,325,350,375,400,425,450)
          cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12],table[table$V1 == indx1[3] & table$V2 == indx2[3],3:12]))))
          #sample size used in the model.
          #1/sqrt(n1)
          n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10),rep.int(1/sqrt(indx1[3]),10))
          #1/sqrt(n2)
          n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10),rep.int(1/sqrt(indx2[3]),10))
          ncombine1<-n11+n22
          #sqrt(n1)
          n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10),rep.int(sqrt(indx1[3]),10))
          #sqrt(n2)
          n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10),rep.int(sqrt(indx2[3]),10))
          ncombine2<-n33+n44
          alpha<-seq(0.01,0.1,0.01)
          log<-rep(log(alpha/(1-alpha)),3)
          #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
          #summary(glm(cv~ncombine1+ncombine2+log-1))
          model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log-1)$coefficients ))
        #  eqn = summary(lm(cv~ncombine1+ncombine2+log))
        #  coeff = as.numeric(eqn$coefficients[,1])
        #  std.e = eqn$sigma
          #1/sqrt(length1)+1/sqrt(length2)
          n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
          #sqrt(length1)+sqrt(length2) 
          n_sqrt<-sqrt(length1)+sqrt(length2)
          alpha1<-0.05
          log1<-log(alpha1/(1-alpha1))
      #    vl = coeff[1]*n_sqrt1+ coeff[2]*n_sqrt + coeff[3]*log1 + std.e
       #   sigma<-std.e^2
       #   sigma0<-var(glm(cv~ncombine1+ncombine2+log)$residuals)
          u0<-model[1]*n_sqrt1+model[2]*n_sqrt+model[3]*log1
          log1<-log(alpha1/(1-alpha1))
          sigma0<-var(glm(cv~ncombine1+ncombine2+log-1)$residuals)*sqrt(30)
       #   u0<-model[1]+model[2]*n_sqrt1+model[3]*n_sqrt+model[4]*log1+sqrt(sigma0)
       #   sigma1<-sigma0*sqrt(30)
          }
         if ((length1>=300) || (length2>=300)) {
          indx1<-c(table[findInterval(length1,table[,1])-60,1],table[findInterval(length1,table[,1]),1])
          indx2<-c(table[findInterval(length2,table[which(table[,1]==indx1[1]),2])-1,2],table[findInterval(length2,table[which(table[,1]==indx1[1]),2]),2])
          cv<-as.numeric(as.character(t(c(table[table$V1 == indx1[1] & table$V2 == indx2[1],3:12],table[table$V1 == indx1[2] & table$V2 == indx2[2],3:12]))))
          #sample size used in the model.
          #1/sqrt(n1)
          n11<-c(rep.int(1/sqrt(indx1[1]),10),rep.int(1/sqrt(indx1[2]),10))
          #1/sqrt(n2)
          n22<-c(rep.int(1/sqrt(indx2[1]),10),rep.int(1/sqrt(indx2[2]),10))
          ncombine1<-n11+n22
          #sqrt(n1)
          n33<-c(rep.int(sqrt(indx1[1]),10),rep.int(sqrt(indx1[2]),10))
          #sqrt(n2)
          n44<-c(rep.int(sqrt(indx2[1]),10),rep.int(sqrt(indx2[2]),10))
          ncombine2<-n33+n44
          alpha<-seq(0.01,0.1,0.01)
          log<-rep(log(alpha/(1-alpha)),2)
          #model with term 1/sqrt(n1), 1/sqrt(n2),sqrt(n1),sqrt(n2),log(alpha/(1-alpha)) and the coefficient for sqrt(n1) and sqrt(n2) are the same,the coefficient for 1/sqrt(n1) and 1/sqrt(n2) are the same 
          #summary(glm(cv~ncombine1+ncombine2+log-1))
          model<-as.matrix(as.numeric(glm(cv~ncombine1+ncombine2+log-1)$coefficients ))
          #1/sqrt(length1)+1/sqrt(length2)
          n_sqrt1<-1/sqrt(length1)+1/sqrt(length2)
          #sqrt(length1)+sqrt(length2)
          n_sqrt<-sqrt(length1)+sqrt(length2)
          alpha1<-0.05
          log1<-log(alpha1/(1-alpha1))
          sigma0<-var(glm(cv~ncombine1+ncombine2+log-1)$residuals)*sqrt(20)  
          u0<-model[1]*n_sqrt1+model[2]*n_sqrt+model[3]*log1
        #  sigma1<-sigma0*sqrt(20)
          } #step II get the q-hat
          #get 200 test-statistics under the null hypothesis
          o.nsims<-200
          nsims<-200
          teststats1<-sort(simulatedata(x,y,nsims))
      #    teststats1_<-simulatedata(x,y,nsims)
          qhat<-Qhat(teststats1,0.05,nsims,u0,sigma0)
      #    qhat<-Qhat(teststats1,0.05,nsims,u0,sigma1)

          #step III get the variance estimator 
          variance_result<-variance(x,y,o.nsims,nsims,qhat)
#           variance_result_<-varianceEstimator(x,y, nsims, qhat, o.nsims)
          teststats2<- variance_result[[1]] 
#           teststats2_<-variance_result_[[1]]
          variance_estimate<-variance_result[[2]]
#          variance_estimate_<-variance_result_[[2]]
          #step IV 
          while((nsims<=35000) & (variance_estimate>sigma0)){
          # Step 4 compute new qhat 
          nsims = length(teststats1)+length(teststats2)
 #         nsims_ = length(teststats1_)+length(teststats2_)
          teststats1 = sort(c(teststats1, teststats2))
#          teststats1_ = c(teststats1_, teststats2_)    
          new.qhat = Qhat(teststats1, 0.05, nsims, u0,sigma0)
#          new.qhat_= Qhat(teststats1_, 0.05, nsims_, u0,sigma0)
          var = variance(x,y, o.nsims,nsims, new.qhat)
#          var_=varianceEstimator(x,y, nsims, new.qhat_,o.nsims)
          teststats2 = var[[1]]
#          teststats2_ = var_[[1]]
          variance_estimate=var[[2]]
#          variance_estimate_=var_[[2]]
          }
          if((nsims>35000) & (variance_estimate>sigma0)){
          cat("Exceded simulation limits without reaching acceptable variance \n")
          qhatC = NA
          }else{  
          # compute final qhat
          nsims = length(teststats1)+length(teststats2)
          teststats1 = sort(c(teststats1, teststats2))
          qhatC = Qhat(teststats1, 0.05, nsims, u0,sigma0)
          }
          output2<-.C("CWrapper1",n1=as.integer(length1),n2=as.integer(length2),y1=as.double(x),y2=as.double(y),test_stat=as.double(1))$test_stat;
          f <-function(alpha_r){
              Jmin = floor(max(c(2,((1-alpha_r)*nsims-sqrt(nsims)*log(nsims)))))
              Jmax = ceiling(min(c(nsims,((1-alpha_r)*nsims+sqrt(nsims)*log(nsims)))))
              Jint = as.numeric(Jmin:Jmax)
              num = exp((-nsims/(2*alpha_r*(1-alpha_r)))*((1-alpha_r - Jint/nsims)^2))*(sqrt(sigma0/(2*pi))*exp((-((teststats1[Jint-1]-u0)^2)/(2*sigma0))-exp(-((teststats1[Jint]-u0)^2)/(2*sigma0)))+u0*(pnorm(teststats1[Jint], u0, sqrt(sigma0))-pnorm(teststats1[Jint-1], u0, sqrt(sigma0))))
              den = exp((-nsims/(2*alpha_r*(1-alpha_r)))*((1-alpha_r - Jint/nsims)^2))*((pnorm(teststats1[Jint], u0, sqrt(sigma0))-pnorm(teststats1[Jint-1], u0, sqrt(sigma0))))
              value.cal = sum(num)/sum(den)
              result=value.cal-output2
              #list=c(result,result^2)
              return(result)
              #return(list)
             }
         f1 <-function(alpha_r) f(alpha_r)^2 
          if (output2>qhatC) {
          Imi<-0.001
          Ima<-0.04
          if ((f(Imi)!="NaN") & (f(Ima)!="NaN")) {
          if (f(Imi)*f(Ima)<=0) p<-uniroot(f,lower=Imi,upper=Ima,tol=0.00001)$root
          if (f(Imi)*f(Ima)>0) p<-optimize(f1,c(Imi,Ima),tol=0.00001)$minimum} else p=0.001
          
          }  
          if (output2<=qhatC) {
          Imi<-0.1
          Ima<-0.5
          if ((f(Imi)!="NaN") & (f(Ima)!="NaN")) {

          if (f(Imi)*f(Ima)<=0) p<-uniroot(f,lower=Imi,upper=Ima,tol=0.00001)$root
          if (f(Imi)*f(Ima)>0) p<-optimize(f1,c(Imi,Ima),tol=0.00001)$minimum} else p=0.5

          }      
          final <- list(test_stat=output2, p_value=p)
          print(final)
          }
}
}
}
