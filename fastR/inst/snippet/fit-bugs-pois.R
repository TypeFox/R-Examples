o <- c(2,10,16,11,5,3,3)
o.collapsed <- c(2+10,16,11,5,3+3)
n <- sum(o)
m <- sum(o * 0:6 ) / n     # mean count = MLE for lambda (full data)
p <- dpois(0:6,m)  
p.collapsed <- c(p[1] + p[2], p[3:5], 1-sum(p[1:5]))   # collapsed probs
e.collapsed <- p.collapsed * n
print(cbind(o.collapsed,p.collapsed,e.collapsed))
lrt  <- 2 * sum(o.collapsed * log(o.collapsed/e.collapsed)); lrt
pearson <- sum( (o.collapsed-e.collapsed)^2/e.collapsed ); pearson
1-pchisq(lrt, df=3)
1-pchisq(pearson,df=3)
