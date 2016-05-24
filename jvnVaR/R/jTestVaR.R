jTestVaR <-
function(Ret, VaR, p,test_significant,type){
v <- jSmallerVaRCount(Ret, VaR)
n <- length(Ret)
k <- length(v)
switch(type,
         p_value = jPValue (n,k,p,test_significant),
         pof = jPofTest (n,k,p,test_significant),
 tuff = jTuffTest (n,v[1],p,test_significant),
 mixkup = jMixKupTest (n,v,p,test_significant), 
  )
}
