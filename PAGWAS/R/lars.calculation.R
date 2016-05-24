lars.calculation <-
function(y,x,z,s2){
    y=y-mean(y)
    u.star=c()
    z=z^(1/2)
    x=x%*%diag(1/z)
    a=lars(y=y,x=x,use.Gram=TRUE,normalize=TRUE)
    u=predict.lars(a,s=2*s2,type="coef",mode="lambda")
    u.star=u$coefficients*(1/z)
    rm(z,x,a,u)
    u.star
}
