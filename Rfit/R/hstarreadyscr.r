#hstarreadyscr = function(ehat,asc=wilasc,ascpr=wilascpr){
# jk 06-19-2011: removed defaults to asc, ascpr 
hstarreadyscr = function(ehat,asc,ascpr){
   n<-length(ehat)
   orde = sort(ehat)
   const = asc[n] - asc[1]
   n = length(ehat)
   ic = (1:n)/(n+1)
   vphi = rep(0,n)
   cn = 0
   for(i in 1:n){
            vphi[i] = ascpr[i]
            cn = cn + ascpr[i]/const
   }
   vad = pairup(orde)
   vwt = pairup(vphi)/const
   vad2 = abs(vad[,1]-vad[,2])
   vwt2 = vwt[,1]+vwt[,2]
   iord = order(vad2)
   absdifford = vad2[iord]
   wtsord = vwt2[iord]
   list(absdifford=absdifford,wtsord=wtsord,cn=cn)
}


