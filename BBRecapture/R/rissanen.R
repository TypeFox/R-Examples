#########################################################
rissanen=function(n) {
#########################################################
  
	n=ceiling(n)
	res=0
	a=log(n,2)
        aa=max(a)
	while((2^aa)>1) {
		res=res+a
                a=log(a,2)*(log(a,2)>0)
		aa=log(aa,2)
}
      
        ind=which(is.na(res))
        
        if(any(ind)){
        temp=c()

res.temp=res
n.single=n[ind]
      
for(j in 1:length(ind)){

       n.new=n.single[j]
	n.new=ceiling(n.new)
	res=0
	a=log(n.new,2)
        aa=max(a)

	while((2^aa)>1) {
		res=res+a
                a=log(a,2)*(log(a,2)>0)
		aa=log(aa,2)
              }

temp=c(temp,res)

     }

res=res.temp
res[ind]=temp

     }

		return(2^(-res))
}
