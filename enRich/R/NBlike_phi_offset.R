NBlike_phi_offset <-
function(para, data1, ppr1, ppr0, mus,mub, k=k)
{
#Log-likelihood based on mixture of NB(mu, size), maximised using optim() function
   N=length(data1)
   phis=para[1]
   phib=para[2]
   temp1=sum(ppr1*ifelse(data1<k, 0, dnbinom(data1-k, phis, ,mus, log=TRUE)))
   temp2=sum(ppr0*dnbinom(data1, phib, ,mub, log=TRUE))
   logl=temp1+temp2
   return(-logl)
}
