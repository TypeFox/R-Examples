winvarN <-
function(x,tr=.2){
#
# rescale the winsorized variance so that it equals one for the standard
# normal distribution
#
x=elimna(x)
cterm=NULL
if(tr==0)cterm=1
if(tr==0.1)cterm=0.6786546
if(tr==0.2)cterm=0.4120867
if(is.null(cterm))cterm=area(dnormvar,qnorm(tr),qnorm(1-tr))+2*(qnorm(tr)^2)*tr
bot=winvar(x,tr=tr)/cterm
bot
}
