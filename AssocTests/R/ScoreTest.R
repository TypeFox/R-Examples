## id.null gives the column id under the null model
ScoreTest <- function(x.mat, y, num.test=3, C.mat, id.null)
{
    x.mat<-as.matrix(x.mat)

    num.subj<-dim(x.mat)[1]/num.test

    num.para<-dim(x.mat)[2]
 
    list.1<-num.test *(c(1:num.subj)-1)+1
 
    ## to get the estimate from the null model

    x.mat.null<-data.frame(x.mat[list.1, id.null])

    y0 <- y[list.1]

    model.null<-glm(y0~-1+., data=x.mat.null, family=binomial)
    
    model.sum<-summary(model.null)

  
    ## to get the bread.mat under the null

    p.vec<-rep(model.null$fitted.values, each=num.test)

    D.mat<-p.vec *(1-p.vec)
    
    temp.0<-scale(t(x.mat), center=F, scale=1/D.mat)

    temp.1<-temp.0 %*% x.mat

    ### inverse of information matrix
    bread.mat<-solve(temp.1)


    ## get the meat under the null
    res.vec<-y-p.vec    
    
    score.mat<-scale(t(x.mat), center=F, scale=1/res.vec)

    score.vec<-apply(score.mat, 1, sum)

    score.vec<-matrix(score.vec, ncol=1, nrow=num.para)
    
    ## see Generalized estiation equations book page 87

    list.mat<-matrix(0, ncol=num.test, nrow=num.subj)
    
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

    score.select<-C.mat %*% score.vec

    bread.select<-C.mat %*% bread.mat %*% t(C.mat)
    sand.select<-C.mat %*% sand.cov %*% t(C.mat)

    temp.2<-solve(bread.select)

    cov.select<-temp.2 %*% sand.select %*% t(temp.2)

    #list(score.select=score.select, cov.select=cov.select)

    beta.cov <- cov.select
    beta.sd  <- sqrt(diag(beta.cov))
    T <- score.select/beta.sd
    cor.mat <- beta.cov/(cbind(beta.sd)%*%rbind(beta.sd))

    list(T, cor.mat)
}
