s2 <-
function(gamm,nstrata , nn1 , nn0 , n1a , n0a , grpa , repp , xx , ofs , yy)
{
	l <- lagrange(gamm,nstrata , nn1 , nn0 , n1a , n0a , grpa , repp , xx , ofs , yy)
  sum(l^2)
}
