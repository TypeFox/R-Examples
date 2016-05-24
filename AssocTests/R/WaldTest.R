WaldTest <- function(x.mat, y, num.test=3, C.mat)
{

    x.mat<-as.matrix(x.mat)
    num.subj<-dim(x.mat)[1]/num.test

    num.para<-dim(x.mat)[2]

    ## to get the estimate

    model.full<-glm(y~-1+., data=data.frame(x.mat), family=binomial)
    
    model.sum<-summary(model.full)

    parms<-model.sum$coefficients[,1]
    
    parms<-matrix(parms, ncol=1, nrow=num.para)
    bread.mat<-model.sum$cov.scaled

    res.vec<-y-model.full$fitted.values    

    score.mat<-scale(t(x.mat), center=F, scale=1/res.vec)

    ## accumulate seveal m of them

    list.mat<-matrix(0, ncol=num.test, nrow=num.subj)

    list.1<-num.test *(c(1:num.subj)-1)+1
    for (i in 1:num.test)
    {
       list.mat[,i]<-list.1 + (i-1)
    }

    acc.mat<-matrix(0, ncol=num.subj, nrow=num.para)

    for (i in 1:num.test)
    {
      acc.mat<-acc.mat + score.mat[, list.mat[,i]]
    }


    meat.mat<-acc.mat %*% t(acc.mat)

    sand.cov<-bread.mat %*% meat.mat %*% bread.mat

    para.select<-C.mat %*% parms
    cov.select<-C.mat %*% sand.cov %*% t(C.mat)

    beta.cov <- cov.select
    beta.sd  <- sqrt(diag(beta.cov))
    T <- para.select/beta.sd
    cor.mat <- beta.cov/(cbind(beta.sd)%*%rbind(beta.sd))

    list(T, cor.mat)
}
