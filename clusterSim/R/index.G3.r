index.G3<-function(d,cl)
{
	d<-data.matrix(d)
     	DU<-0
	r<-0
	v_max<-array(1,max(cl))
	v_min<-array(1,max(cl))
	for (i in 1:max(cl))
	{
		
		n<-sum(cl==i)
		if (n>1)
		{
			t<-d[cl==i,cl==i]
			DU=DU+sum(t)/2
			v_max[i]=max(t)
			if (sum(t==0)==n)		# nie ma zer poza przekatna
				v_min[i]<-min(t[t!=0])
			else
				v_min[i]<-0
			r<-r+n*(n-1)/2
		}
	}
	Dmin=min(v_min)
	Dmax=max(v_max)
	if(Dmin==Dmax)
		result<-NA
	else
		result<-(DU-r*Dmin)/(Dmax*r-Dmin*r)
	result		
}
