Andrews.Chen1 <-
function(x,p,h)
{
x <- as.matrix(x)
n <- nrow(x)
b <- LSE(x,p)$coef
alphau <- alpha.u0(b,p,n)
if(p == 1 & alphau == 1) newb <- rbind(1,0)
else newb <- rbind(alphau,estmf0(x,p,alphau))

tem1 <- newb[1]
tem <- newb; tem[1,1] <- b[1,1] 
{
if(p == 1 | alphau == 1) newb2 <- newb
else
    {
    for( i in 1:10)
    {
    alphau <- alpha.u0(tem,p,n)
    newb2 <- rbind(alphau,estmf0(x,p,alphau))
    tem2 <- newb2[1]
    {
    if (abs(tem1-tem2) < 0.01) break
    else
    tem1 <- tem2; tem <- newb2; tem[1,1] <- b[1]}
    }
    }

}
b<-arlevel(newb2,p)
e <- RESID(x,b)
f <- {}
if(h > 0)
f <- AR.Fore(x,b,h)
return(list(coef=b,ecmcoef=newb2,resid=e,forecast=f))
}
