freqs<-function(pp)
{
	m<-union(pp$marks,NULL)
	h<-NULL
	for(i in m) h<-c(h, sum(pp$marks==i))
	names(h)<-m
	h
}
####################################################################################
clean.up.data<-function(pp,dbh=10,atleast=10)
{
	h<-freqs(pp)
	m<-union(pp$marks,NULL)
	i<-pp$dbh>dbh & pp$marks%in%m[h>atleast]
	p<-pp[i]
        marks<-p$marks
        p$marks<-NULL
		p$markformat<-"none"
	p$dbh<-pp$dbh[i]
        p$marks<-factor(marks)
	p
}

####################################################################################
shake<-function(pp, a=0.001)
{
	pp$x<-pp$x+runif(pp$n, -a,a)
	pp$y<-pp$y+runif(pp$n, -a,a)
	pp
}
####################################################################################
minusID<-function(pp, minusR, dbh, atleast=0)#, lifeform, status)
{
	id<-( pp$x<(pp$window$x[2]-minusR) ) & (pp$x>(pp$window$x[1]+minusR)) & (pp$y<(pp$window$y[2]-minusR)) & (pp$y>(pp$window$y[1]+minusR)) 
	if(!is.null(pp$dbh) & !missing(dbh))
		id<-id & pp$dbh>dbh
	if(atleast>0)
	{
		fqs<-freqs(pp)
		marks<-names(fqs)
		for(m in 1:length(marks))	
		{
			bad<-which(pp$marks==marks[m])
			if(fqs[m]<atleast)
				id[bad]<- FALSE
		}
	}
#	if(!missing(lifeform))
#	{
#		sp<-levels(pp$marks)
#		lfs<-bcispecies.lifeform(lifeform=lifeform)
#		bad<-sp[!(sp%in%lfs)]
#		for(s in bad)
#			id[pp$marks==s]<-FALSE
#	}
#	if(!missing(status))
#	{
#		id[!(pp$status%in%status)]<-FALSE
#	}
	if(sum(is.na(id))>0) warning("Vector contains NA's.")
	id
}
