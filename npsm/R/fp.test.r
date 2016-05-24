fp.test = function(x,y,delta0=0,alternative='two.sided'){
    y = y - delta0
    ds = matrix(outer(y,x,"-"),ncol=1)
    m = length(x)
    n = length(y)
    mu = m*n/2
    ind = rep(0,m*n)
    ind[ds  > 0] = 1
    ind[ds  == 0] = .5
    ts = sum(ind) - mu
    
#    pl = rep(0,m)
#    ql = rep(0,n)
#    for(i in 1:m){
#       tmp = x[i]
#       pl[i] = place(tmp,y)
#    }
   pl <- apply(as.matrix(x),1,place,y) # placement of x in y
#    for(j in 1:n){
#       tmp = y[j]
#       ql[j] = place(tmp,x)
#    }
   ql <- apply(as.matrix(y),1,place,x) # placement of y in x
    v1 = (m-1)*var(pl)
    v2 = (n-1)*var(ql)
   sig = sqrt(v1+v2+(mean(pl)*mean(ql)))
   zp = ts/sig
   std = zp
   zp = (ts - .5)/sig
   zn = (ts  + .5)/sig
   if(alternative=='greater'){
        pval = 1 - pnorm(zp)
        zs = zp
   }
    if(alternative=='less'){
        pval = pnorm(zn)
        zs = zn
   }
    if(alternative=='two.sided'){
        if(ts >= 0){
            pval = 2*(1 - pnorm(abs(zp)))
            zs = zp
        } else {
            pval = 2*pnorm(zn)
            zs = zn
        }
   }
res<-list(statistic = std,p.value=pval,numerator=ts,denominator=sig)
class(res)<-'rank.test'
res
}

place = function(x,y){
#
#    Placement of x in the vector y
#
     ic = 0
     n = length(y)
     ac = 0
     ys = sort(y)
     i = 1
     while(ac==0){
        tmp = ys[i]
        if(x > tmp){
              i = i + 1
              ic = ic + 1
        } else { 
              ac=1 
        }
        if(i > n){ac=1}
     }
     ic
}
