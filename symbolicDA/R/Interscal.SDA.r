.data2dist<-function(x){
#  dolne<-x[,seq(2,ncol(dane),2)]
#  gorne<-x[,seq(3,ncol(dane),2)]
  dolne<-x[,,1]
  gorne<-x[,,2]
  zmienne<-ncol(dolne)
  obiekty<-nrow(dolne)
  dist_min<-array(,c(obiekty,obiekty))
  dist_max<-array(,c(obiekty,obiekty))
  for (i in 1:obiekty)
  {
    for (j in 1:obiekty)
    {
      for (k in 1:zmienne)
      {
        dist_max[i,j]<-dist_max[j,i]<-0.5*sqrt( sum( ( (gorne[i,k]-dolne[i,k])+(gorne[j,k]-dolne[j,k])+2*abs( (gorne[i,k]-dolne[i,k])/2-(gorne[j,k]-dolne[j,k])/2 ) )^2 ) )
        dist_min[i,j]<-dist_min[j,i]<-0.25*sqrt( sum( ( (gorne[i,k]-dolne[i,k])+(gorne[j,k]-dolne[j,k])+2*abs( (gorne[i,k]-dolne[i,k])/2-(gorne[j,k]-dolne[j,k])/2 )- !((gorne[i,k]-dolne[i,k])+(gorne[j,k]-dolne[j,k])+2*abs( (gorne[i,k]-dolne[i,k])/2-(gorne[j,k]-dolne[j,k])/2 )) )^2 ) )
      }
    }
  }
  list(dist_max=dist_max,dist_min=dist_min)
}

interscal.SDA<-function(x,d=2,calculateDist=FALSE){
if(.is.symbolic(x)){
x<-SO2Simple(x)
}
if(calculateDist){
  minMax<-.data2dist(x)
  dist_min<-minMax$dist_min
  dist_max<-minMax$dist_max
}
else{
  dist_min<-x[,,1]
  dist_max<-x[,,2]
}
z<-ncol(dist_min)
delta<-array(,c(z*2,z*2))
for (i in 1:z)
{
	for (j in 1:z)
	{
		if(i==j)
		{
			delta[2*i,2*j]<-delta[2*i-1,2*j-1]<-0
			delta[2*i,2*j-1]<-delta[2*i-1,2*j]<-dist_max[i,j]
		}
		else
		{
			delta[2*i,2*j]<-dist_max[i,j]
			delta[2*i-1,2*j-1]<-dist_min[i,j]
			delta[2*i,2*j-1]<-delta[2*i-1,2*j]<-(dist_max[i,j]+dist_min[i,j])/2			
		}
	}
}
b<-array(,c(z*2,z*2))
k<-z*2
suma_d<-sum(delta^2)
for (i in 1:k)
{
	for (j in 1:k)
	{
	
	b[i,j]<-(-0.5)*(delta[i,j]^2-(1/k)*sum(delta[,j]^2)-(1/k)*sum(delta[i,]^2)+(1/(k^2))*suma_d)
		
	}
}

wartosci<-as.matrix(eigen(b)$values)
wektory<-as.matrix(eigen(b)$vectors)
k_2m<-array(,c(z*2,d))
for (r in 1:k)
  {
	for (i in 1:d)
	{
		k_2m[r,i]<-sqrt(abs(wartosci[r,]))*wektory[i,r]
	}
}
k_dolne<-array(,c(z,d))
k_gorne<-array(,c(z,d))
for (i in 1:z)
{
	for (j in 1:d)
	{
		k<-2*i-1
		t<-2*i
			k_dolne[i,j]<-min(k_2m[k,j],k_2m[t,j])
			k_gorne[i,j]<-max(k_2m[k,j],k_2m[t,j])
	}
}
xprim<-array(,c(z,d,2))
for (i in 1:z)
{
	for (j in 1:d){
    xprim[i,j,1]<-k_dolne[i,j]
    xprim[i,j,2]<-k_gorne[i,j]
	}
}


minimum<-min(k_dolne)
maksimum<-max(k_gorne)
mini<-abs(minimum)
maks<-abs(maksimum)
x_srodki<-array(,c(z,d))
r_bokow<-array(,c(z,d))
	for (i in 1:z)
	{
		for (s in 1:d)
		{
		x_srodki[i,s]<-abs(k_gorne[i,s]-k_dolne[i,s])
		r_bokow[i,s]<-(k_dolne[i,s]+(x_srodki[i,s]/2))
		}
	}
macierz_d_dolne<-array(,c(z,z))
macierz_d_gorne<-array(,c(z,z))

	for (i in 1:z)
	{
		for (j in 1:z)
		{
			for (s in 1:d)
			{
				macierz_d_dolne[i,j]<-macierz_d_dolne[j,i]<-sqrt(sum(max(0,(abs(x_srodki[i,s]-x_srodki[j,s])+r_bokow[i,s]+r_bokow[j,s])^2)))
				macierz_d_gorne[i,j]<-macierz_d_gorne[j,i]<-sqrt(sum(abs(x_srodki[i,s]-x_srodki[j,s])+r_bokow[i,s]+r_bokow[j,s])^2)
			}
		}
	}
	
STRESSSym<-0
	for (i in 1:z)
	{
		for (j in 1:z)
		{
		if(i<j) STRESSSym<-STRESSSym+(1/z^2*(dist_max[i,j]-macierz_d_gorne[i,j])^2)+(1/z^2*(dist_min[i,j]-macierz_d_dolne[i,j])^2)
		}
	}
		mian1<-0  
    mian2<-0
			for (k in 1:(z-1))
    {
      for (l in (k+1):z)
      {
        mian1<-mian1+(1/z^2*macierz_d_dolne[k,l]^2)
        mian2<-mian2+(1/z^2*macierz_d_gorne[k,l]^2)
      }
      }
      STRESSSym<-STRESSSym/(mian1+mian2)
      res<-list(xprim=xprim,stress.sym<-STRESSSym)
      res
  
}
