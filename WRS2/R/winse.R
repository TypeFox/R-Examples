winse<-function(x, tr = .2){
#
# Estimate the standard error of the Winsorized mean
#
x=elimna(x)
n=length(x)
h=n-2*floor(tr*n)
top=(n-1)*sqrt(winvar(x,tr=tr))
bot=(h-1)*sqrt(n)
se=top/bot
se
}