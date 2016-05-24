ufun <-
function(xdat,grp,nperm=50,delta=0.1,lim=c(-5,5), 
       normalize=FALSE)
{

#require(rsmooth)
# discretizations
 zz = seq(lim[1], lim[2],by = delta)
 zk = zz[-1] - delta/2
 nk = length(zk)

if (normalize){
  phi0 = dnorm(zk)
  u1 = dnorm(zk)*(zk^2 -1)
  u1 = u1/sqrt(sum(u1^2))
  return(list(zk=zk,phi0=phi0,u1=u1, delta=delta))
}

# .....................  if using permutation based result
Y= NULL
for (i in 1:nperm){
  grp.star = sample(grp)
  stat = tstatistics2(xdat, grp.star)$tstat 
  # grouped data  
  y = table(cut(stat,zz))
  Y = cbind(Y,y)
}

# .............. direct SVD analysis
Ym = apply(Y,1,mean)
Yc = Y -Ym   # mean centered
svdy = svd(Yc)

phi0 = Ym/sum(Ym)/delta  
#d1 = rsmooth(zk,svdy$u[,1], k=1, lambda=.001)$y
d1 = lowess(zk,svdy$u[,1])$y
sign1=d1[nk/2]/abs(d1[nk/2])

return(list(zk=zk, phi0=phi0, u1=-d1/sign1, delta=delta, Y=Y))
}
