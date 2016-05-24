"decay.spec" <-
function(qm,p,lambda)
  {
	ld<-seq(1:p)^-lambda;			# regular lag decay (note ^-lambda) 
	if(qm==12)
                 { j<-ceiling(p/3)^-lambda;   	# last quarter (rounded up) eg. l2=13=>xx2=5
                   b<-0
                   if(p > 1)
                     { b<-(log(1)-log(j))/(1-p) }
                   a<-exp(-b);
                   ld<-a*exp(b*seq(1:p));	# Tao Zha's lag decay to match 13th lag 
                 }
   (ts(ld))
      }

