BootsModelTest<-function(x,k,m,B,Alpha,Beta,tol=1.0e-7, maxit=50,trace=FALSE)
{
    if(B>=10) cat("It will take a while to do the bootstrap test.\nPlease be patient.\n")
    Delta<-Deltan.ModelTest(x,k,m,Alpha,Beta)
#    cat("Delta=",Delta,"\n")
    if(is.vector(x)){
        n<-length(x)
        x<-matrix(x,n,1)
    }
    if(!is.matrix(x)) x<-as.matrix(x)
    if(is.matrix(x)) n<-nrow(x)
    P<-NULL; Q<-NULL
    P<-1/((m-k)*exp(Alpha+x%*%Beta)+(n-m+k))
    Q<-exp(Alpha+x%*%Beta)*P
    y<-matrix(0,n,ncol(x)); Freq<-0   
    for(b in 1:B){
        if(k>=1 && k<n){
            temp<-sample(1:n, k,replace = T,prob = P)
            y[1:k,]<-as.matrix(x[temp,])
            if(m<n){
                temp<-sample(1:n, (n-m),replace = T,prob = P)
                y[(m+1):n,]<-as.matrix(x[temp,])
            }
            temp<-sample(1:n,(m-k),replace = T,prob = Q)
            y[(k+1):m,]<-as.matrix(x[temp,])
        }
        else{
            temp<-sample(1:n, n,replace = T, prob = P)
            y<-as.matrix(x[temp,])
        }
        if(m<n){
            Temp<-SemiparChangePoint(y, alternative = "epidemic",tol = 1.0e-7, maxit = 50,trace = F)
            m.temp<-Temp$m.hat; 
#            cat("Iteration ",b,"k-hat=",Temp$k.hat,"m-hat=",Temp$m.hat,Temp$alpha.hat,Temp$beta.hat)
        }
        else{
            Temp<-SemiparChangePoint(y, alternative = "one.change",tol = 1.0e-7, maxit = 50,trace = F)
            m.temp<-n
#            cat("Iteration ",b,"k-hat=",Temp$k.hat,Temp$alpha.hat,Temp$beta.hat)
        }
        k.temp<-Temp$k.hat; 
        delta<-Deltan.ModelTest(y,k.temp,m.temp,Temp$alpha.hat,Temp$beta.hat)
        Freq<-Freq+(delta>Delta)
#        cat("delta=",delta,"P-value=",Freq/b,"\n")
    }
    Pvalue<-Freq/B
    list(Delta = Delta, Pvalue = Pvalue)
}
