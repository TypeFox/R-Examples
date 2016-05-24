q.Johansen <-
function(p,nComT,type=c("z","c","lc","l","ql") )
{
x=nComT
t=type
#
e3=c(x^2,x,1,0,0,x^0.5)
e1=c(x^2,x,1,1,0,x^0.5)
e2=c(x^2,x,1,0,1,x^0.5)
#
v3=c(x^2,x,1,0,0)
v1=c(x^2,x,1,1,0)
v2=c(x^2,x,1,0,1)
#
P=t(matrix(c(2,2,2,2,2,
-1,2.01,1.05,4.05,2.85,
0.07,0,-1.55,0.5,-5.1,
0.07,0.06,-0.5,-0.23,-0.1,
0,0.05,-0.23,-0.07,-0.06,
0,0,0,0,1.35,
3,3,3,3,3,
-0.33,3.6,1.8,5.7,4,
-0.55,0.75,0,3.2,0.8,
0,-0.4,-2.8,-1.3,-5.8,
0,-0.3,-1.1,-0.5,-2.66),5,11))
#
if (x>=3)          { e=e3
   v=v3} else if(x==2) { e=e2
   v=v2} else if(x==1) { e=e1
   v=v1} else { stop("nComT is the number of of common trends, so it has to be an positive integer.")
   }
if (t=="z")           { ap=P[,1]
   } else if (t=="c") { ap=P[,2]
   } else if (t=="lc"){ ap=P[,3]
   } else if (t=="l") { ap=P[,4]
   } else if (t=="ql"){ ap=P[,5]
   } else { stop("type must be either z,c,lc,l or ql.")
   }
#
E=t(ap[1:6])%*%e
V=t(ap[7:11])%*%v
#
qgamma((1-p),(E^2)/V,E/V)
}