####### guess starting values (beta-version) 
etas.starting = function(cat.orig,
m0=2.5,
p.start=1,
a.start=1.5,
gamma.start=0.5,
q.start=2,
longlat.to.km=TRUE,
sectoday=TRUE
)
{
cat=cat.orig[cat.orig$magn1>m0,]
if(sectoday)cat$time=cat$time/86400
t=cat$time

a=cat[order(t),]
t=a$time
		if (longlat.to.km){
		
		radius=6371.3
		

		ycat.km      =   radius*cat$lat*pi/180
		xcat.km      =   radius*cat$long*pi/180

		cat$long=xcat.km
		cat$lat=ycat.km
		
		  }

n=nrow(a)
dt=diff(t)
ds=diff(ycat.km)^2+diff(xcat.km)^2
c.start=quantile(dt,0.25)
d.start=quantile(ds,0.05)
mu.start=n*0.5/diff(range(t))
# approximate evaluation of integrals
tmax=max(t)
it=log((c.start+tmax-t)/c.start)
em=exp((a.start-gamma.start)*(a$magn1-m0))
is=d.start^(1-q.start)*exp(gamma.start*(a$magn1-m0))*pi/(q.start-1)
k0.start=n*0.5/sum(em*it*is)
return(list(
mu.start=mu.start,
k0.start=k0.start,
c.start=c.start,
p.start=p.start,
a.start=a.start,
gamma.start=gamma.start,
d.start=d.start,
q.start=q.start,
longlat.to.km=longlat.to.km,
sectoday=sectoday
))
}


