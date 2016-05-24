CTMC.simulate <-
function(rate.fn, jumpLim.fn,T.time, init.state)
{
  ## note using rate.fn as a func rather htan matrix allows
  ## to not worry about indexing issues.  all integer states are fine.
 t.sim = T.time
 idx = 1
 v.sz = 16
 v = matrix(0, ncol=2, nrow=v.sz)
 cur = init.state
 while(t.sim > 0)
   {
     if(idx > v.sz)
       {
         tmp = matrix(0, ncol=2, nrow=v.sz*2)
         tmp[1:v.sz,] = v
         v = tmp
         v.sz = v.sz * 2
       }
     v[idx, 1] = T.time - t.sim
     v[idx, 2] = cur
     idx = idx + 1     
     curRate <- rate.fn(cur,cur);
     if (curRate==0) {
       t.sim <- 0
     }
     else {
       wait.time = rexp(1, rate=-curRate)
       t.sim = t.sim - wait.time
       jumpLims <- jumpLim.fn(cur);
       cur.rates = rate.fn(cur, jumpLims[1]:jumpLims[2])
       cur.rates[cur-jumpLims[1]+1] = 0
       cur = which(rmultinom(1, 1, cur.rates) == 1) + jumpLims[1]-1;
     }
   }
 v = v[1:(idx-1),]
 if (is.vector(v)) {
   v <- matrix(v,nrow=1,byrow=TRUE) #not sure of better way
 }
 return( new("CTMC", states=v[,2],times=v[,1], T=T.time)  )
}

