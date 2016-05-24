index.S<-function(d,cl,singleObject=0)
{
	d<-as.matrix(d)
	Si<-0
	for(k in 1:max(cl))
	{
		if ((sum(cl==k))<=1)
			Sil<-singleObject
		else
		{
			Sil<-0
			for(i in 1:length(cl))
			{
				if(cl[i]==k)
				{
					ai<-sum(d[i,cl==k])/(sum(cl==k)-1)
					dips<-NULL
					for(j in 1:max(cl))
						if (cl[i]!=j)
							if(sum(cl==j)!=1)
								dips<-cbind(dips,c( (sum(d[i,cl==j])) /(sum(cl==j)) ))
							else
								dips<-cbind(dips,c( (sum(d[i,cl==j]))))
					bi<-min(dips)
					Sil<-Sil+(bi-ai)/max(c(ai,bi))
				}
			}
#			Sil<-Sil/sum(cl==k)
		}
		Si<-Si+Sil
	}
#	Si/max(cl)
	Si/length(cl)
}
