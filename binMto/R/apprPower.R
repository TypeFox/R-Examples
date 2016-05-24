"apprPower" <-
function(n, pH1, alpha=0.05, alternative="greater", method="Add4")
{


nc <- n[1]
nx <- n[-1]

pc <- pH1[1]
px <- pH1[-1]

k=length(px)

if(k!=length(nx)) {stop("p and n differ in length")}


if(alternative=="less" & all(px-pc > 0)) 
 {cat("Arguments pH1 and alternative mis-specified","\n", "pi-p0 should be less than 0 for alternative='less'", "\n")}

if(alternative=="greater" & all(px-pc < 0))
 {cat("Arguments pH1 and alternative mis-specified","\n", "pi-p0 should be greater than 0 for alternative='greater'", "\n")}



# # because with I want to show that px-pc < 0 for the diff to be detected 
# #    if alternative="less"-> upper bound
# # because with I want to show that px-pc > 0 for the diff to be detected 
# #    if alternative="greater"-> lower bound


# # define the add-4 p and n-tilde
if(method=="Add4")
 {
  pxI <- (px*nx+1)/(nx+2)
  pcI <- (pc*nc+1)/(nc+2)
  nxI <- nx+2
  ncI <- nc+2
 }
else{
 if(method=="Add2")
  {
   pxI <- (px*nx+0.5)/(nx+1)
   pcI <- (pc*nc+0.5)/(nc+1)
   nxI <- nx+1
   ncI <- nc+1
  }
else{
 if(method=="Wald")
  {
   pxI <- px
   pcI <- pc

   nxI <- nx
   ncI <- nc
  }
}}

expectH1 <- (pxI - pcI) / sqrt( (pxI*(1-pxI)/nxI)+(pcI*(1-pcI)/ncI) )

# # # function for the correlation matrix

 corrmat <- function(n0, nx, p0, px)
 {
 k=length(nx)
 # here, k=the number of treatment groups
 biv <- numeric(length=k)
 # correct for zeros
 if(p0==0){p0 <- (0.5)/n0}
 if(p0==1){p0 <- (n0-0.5)/n0}
 
  for(i in 1:k)
   {
   biv[i] <- 1/sqrt(1 + n0/nx[i] * ( ( px[i]*(1-px[i]) ) / ( p0*(1-p0) ) )  ) 
   }

 corr <- diag(x=1, nrow=k)
 # the columns:
 for(e in 1:k)
 {
  for(a in 1:k)
   {   
   if( e !=a ) {corr[e,a] <- biv[e]*biv[a]}
   }
 }

 return(corr)
 }
 # end of corrmat

corrmatH1<-corrmat(n0=nc, nx=nx, p0=pc, px=px)


if(alternative=="greater")
{
H0quant<-qmvnorm(1-alpha, corr=corrmatH1, tail="lower")$quantile
vec <- (H0quant-expectH1)
power <- 1-pmvnorm(upper=vec, corr=corrmatH1)[1]
}


if(alternative=="less")
{
H0quant<-qmvnorm(1-alpha, corr=corrmatH1, tail="upper")$quantile
vec <- (H0quant-expectH1)
power <- 1-pmvnorm(lower=vec, corr=corrmatH1)[1]
}


if(alternative=="two.sided")
{

H0quant<-qmvnorm(1-alpha, corr=corrmatH1, tail="both")$quantile
veclower <- (-H0quant-abs(expectH1))
vecupper <- H0quant-abs(expectH1)

power <- 1-pmvnorm(lower=veclower, upper=vecupper, corr=corrmatH1)[1]
}

return(power)
}

