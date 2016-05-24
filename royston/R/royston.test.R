royston.test <-
function(a){
p=dim(a)[2]
n=dim(a)[1]
z <- matrix(nrow=p, ncol = 1)
z <- as.data.frame(z)
w <- matrix(nrow=p, ncol = 1)
w <- as.data.frame(w)

if(n<=3) 
{
stop("n must be greater than 3")
}

else if (n >=4 || n <= 11) 
{
  x = n
    g = -2.273 + 0.459*x
    m = 0.5440 - 0.39978*x + 0.025054*x^2 - 0.0006714*x^3
    s = exp(1.3822 - 0.77857*x + 0.062767*x^2 - 0.0020322*x^3)

    for(i in 1:p){

a2=a[,i]
{

if (kurtosis(a2) > 3) 
{
w=sf.test(a2)$statistic
}

else
{

w=shapiro.test(a2)$statistic

}

}

z[i,1]=(-log(g - (log(1 - w))) - m)/s

            }
}
if(n>2000) 
{
stop("n must be less than 2000")
}

else if (n >=12 || n <= 2000) 

{
   x = log(n)
    g = 0
    m = -1.5861 - 0.31082*x - 0.083751*x^2 + 0.0038915*x^3
    s = exp(-0.4803 -0.082676*x + 0.0030302*x^2)

       for(i in 1:p){

a2=a[,i]

{

if (kurtosis(a2) > 3) #Shapiro-Francia test is better for leptokurtic samples
{
w=sf.test(a2)$statistic
}

else #Shapiro-Wilk test is better for platykurtic samples
{

w=shapiro.test(a2)$statistic

}

}
 z[i,1] =  ((log(1 - w)) + g - m)/s
}
}

else {
            stop("n is not in the proper range")
        }

u = 0.715
v = 0.21364 + 0.015124*(log(n))^2 - 0.0018034*(log(n))^3
l = 5
C = cor(a) #correlation matrix
NC = (C^l)*(1 - (u*(1 - C)^u)/v)#transformed correlation matrix
T = sum(sum(NC)) - p# %total
mC = T/(p^2 - p)# %average correlation
edf = p/(1 + (p - 1)*mC)# %equivalent degrees of freedom
Res <- matrix(nrow=p, ncol = 1)
Res <- as.data.frame(Res)
 for(i in 1:p){
Res=(qnorm((pnorm( - z[,]))/2))^2
}
Sa=cov(a)
D2 <- mahalanobis(a, colMeans(a), Sa)
qqplot(D2,qchisq(ppoints(n), df=p) , main="Mahalanobis Distance vs Chi-Square",
xlab="Mahalanobis Distance", ylab="Chi-Square")
abline(0, 1, col = 'black')
RH = (edf*(sum(Res)))/p
pv=pchisq(RH, edf,lower.tail = FALSE)

{
if (pv<0.05){
dist.check=("Data analyzed have a non-normal distribution.")
}
else
{
dist.check=("Data analyzed have a normal distribution.")
}
}
res = structure(list(H=RH,p.value=pv,distribution=dist.check))
res
}
