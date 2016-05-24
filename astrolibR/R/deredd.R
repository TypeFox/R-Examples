deredd = function(eby,by,m1,c1,ub) {

 rm1 = -0.33
 rc1 = 0.19
 rub = 1.53 
 eby0 = eby >0

 by0 = by - eby0
 m0 = m1 - rm1*eby0
 c0 = c1 - rc1*eby0
 ub0 = ub - rub*eby0
 
 return(list(by0=by0,m0=m0,c0=c0,ub0=ub0))
}
