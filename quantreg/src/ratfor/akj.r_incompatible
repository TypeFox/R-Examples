subroutine akj(x,z,p,iker,dens,psi,score,nx,nz,h,alpha,kappa,xlam)
# univariate kernel density-score estimator
# the algorithm is basically from Silverman as adapted for Portnoy and Koenker
# Annals paper on adaptive L-estimation in regression.

# x--pts used for centers of kernel assumed to be sorted!!!!
# z--pts at which density is calculated
# p--probability associated with x's
# dens--f(z), the densities at z
# psi--f'(z)/f(z) the score at z
# score--(log(f(z)))'', the J fn at z
# nx--number of pts in x
# nz--number of pts in z
# iker--kernel
#          0=gaussian
#	   1=cauchy
# h--initial window size (overall)--choose zero for default
# kappa--constant determining initial (default) window width
# xlam--Silverman's lambda, window adjustment for each x

double precision dens(nz),score(nz),psi(nz),h,kappa
double precision z(nz),x(nx),xlam(nx),p(nx),qrange,pi
double precision con1,sum,sqsum,xsd,a,fifth,hinv,half
double precision xn,xker,dxker,ddxker,fact,xponen,alpha,glog,zero,one,two
parameter( zero = 0.d0)
parameter( one = 1.d0)
parameter( two = 2.d0)
parameter( four = 4.d0)
parameter( half = 0.5d0)
parameter( fifth = 0.2d0)
parameter( pi = 3.141593d0)
xn=nx

# call srtad(x,1,nx) # port sort routine now done in S interface.
if(iker==0) con1=one/sqrt(2.0*pi)
else if(iker==1) con1=one/pi

# if no h is provided, calculate a default
if(h<=0.) {
	sum=0.
	sqsum=0.
	do i=1,nx
		{ sqsum=sqsum+x(i)*x(i)*p(i)
		sum=sum+x(i)*p(i)
		}
	xsd=dsqrt(sqsum-sum*sum)
	sum=zero
	for(i=1;i<nx;i=i+1) {
		sum=sum+p(i)
		if(sum<.25) next
		else { qrange=x(i);break }
		}
	sum=one
	for(i=nx;i>0;i=i-1) {
		sum=sum-p(i)
		if(sum>.75) next
		else { qrange=x(i)-qrange;break }
		}
	a=min(xsd,qrange/1.34)
	h=kappa*a/(xn**fifth) # see Silverman p 48
	}
hinv=one/h
# Stage one:  compute pilot estimate of density
	do j=1,nx {
		xker=0.
		if(iker==0) {
			do i=1,nx {
				xponen=(x(j)-x(i))*hinv
				xponen=half*xponen**2
				xker=xker+p(i)*exp(-xponen)*hinv
				}
			}
		else if(iker==1) {
			do i=1,nx {
				xponen=(x(j)-x(i))*hinv
				xker=xker+p(i)*hinv/(1+xponen**2)
				}
			}
		xlam(j)=con1*xker
		}
# Stage two:  Automatic window widths (Silverman p101)
	glog=zero
	do i=1,nx
		{ glog=glog+p(i)*log(xlam(i))}
	g=exp(glog)
	ginv=one/g
	do i=1,nx {
            xlam(i)=hinv/((xlam(i)*ginv)**(-alpha))
        }
# notice xlam no longer its own self at this pt! xlam is 1/(h*lambda(i))
# substitution of * for / thus achieved speeds things up

# Stage two:  new density-score estimates
	do j=1,nz {
		xker=zero
		dxker=zero
		ddxker=zero
		if(iker==0) {  # gaussian kernel
                    do i=1,nx {
                        xponen=(z(j)-x(i))*xlam(i)
                        fact=exp(-half*xponen*xponen)*xlam(i)
                        xker=xker+p(i)*fact
                        dxker=dxker-p(i)*fact*xponen*xlam(i)
                        ddxker=ddxker-
                            p(i)*fact*(one - xponen**2)*xlam(i)**2
                    }
                }
		else if(iker==1) {  # cauchy kernel
                    do i=1,nx {
                        xponen=(z(j)-x(i))*xlam(i)
                        fact=xlam(i)/(one+xponen**2)
                        xker=xker+p(i)*fact
                        dxker=dxker-p(i)*two*xponen*fact**2
                        ddxker=ddxker-
                            p(i)*two*(fact**2)*(xlam(i)-
                                                four*(xponen**2)*fact)
                    }
                }
		dens(j)=con1*xker
		psi(j)=-(dxker/xker)
		score(j)=(dxker/xker)**2-ddxker/xker
		}
return
end
