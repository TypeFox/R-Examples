"precision.grubbs.ncb.od" <-
function (x,beta.bar.x=beta.bar(x)) 
{
# Grubbs Precision estimates for NonConstant Bias using Original Data
# Jaech, p. 184

  sqrt(diag(var(x))-beta.bar.x^2*process.sd(x)^2)

}
