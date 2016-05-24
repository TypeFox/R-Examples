Burgman.fun <-
function(dd){
				names(dd)<-c("yrs", "sights")
				N<-seq(1,sum(dd$sights))
				CT<-length(dd$yrs-dd$yrs[1])
				ss<-subset(dd, dd$sights>0)
				s1<-c(dd$yrs[length(dd$yrs)], dd$sights[length(dd$sights)])
				ss<-rbind(ss, s1)
				r<-max(diff(ss$yrs))
				Sum<-0
				for(j in 1:length(N)){
				  TT<-seq(1, j)
				  KK<-seq(1, j+1)
				  for(k in 1:length(KK)){
				      if(k<=(CT/r)){
				      fallfac = 1
				      for(indexn in 1:length(TT)){
					        n = indexn - 1
					        fallfac = fallfac*(CT - (r*k) - n)
					      }
				      dummy = 0
				      for(indexi in 1:length(KK)){
				        i = indexi-1
				        dummy = dummy + (((-1)^i)*choose(j,i)*((j-i)^length(N)))
      					}
 			     Sterl = (1/factorial(j))*dummy
		   		 calc = ((-1)^(k+1))*choose(j+1,k)*fallfac*Sterl
      			}
      			else{calc = 0}
      			Sum = Sum + calc
      			}
			}
			p = (CT^(-length(N)))*Sum
			return(p)
			}
