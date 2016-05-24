fma.em.alg <-
function(y,numobs,r,k,p,x.z,q.z,x.w,q.w,phi,it,H,w,Beta,sigma,eps,psi,lik)
                     
                    
{ 
                  
likelihood<-NULL 
hh<-0 
ratio<-1000
lik<--10000

chsi<-array(0,c(k,r,r)) 
roy<-array(0,c(k,r,numobs)) 
ph.y<-matrix(0,k,numobs) 
py.h<-matrix(0,k,numobs)

while ((hh < it) & (ratio > eps )) {
 hh<-hh+1 
    
#### E step determino p(z|y,h=i)

Ezz.hy<-array(0,c(k,r,r,numobs))

sigma.tot<-array(0,c(k,p,p)) 
cambia<-rep(FALSE,k)
if ((k==1) & (r==1)) {sigma.tot[1,,]<-H%*%sigma[1]%*%t(H)+psi
                       if (sum(diag(sigma.tot[1,,])<0)) cambia<-TRUE
                       if (det(sigma.tot[1,,])<0) cambia<-TRUE 
                        } else for (i in 1:k) {sigma.tot[i,,]<-H%*%sigma[i,,]%*%t(H)+psi
                                                if (sum(diag(sigma.tot[i,,])<0)) cambia[i]<-TRUE
                                                if (det(sigma.tot[i,,])<0) cambia[i]<-TRUE}
                                                
if ((k==1) & (r==1)) { 
                chsi<-1/(t(H)%*%ginv(psi)%*%H+solve(sigma[1]))
                roy[1,,]<-chsi%*%(t(H)%*%ginv(psi)%*%t(y)+solve(sigma[1])%*%Beta[1,,]%*%t(x.z))
                Ezz.hy[1,,,]<-matrix(chsi,ncol=numobs)+roy[1,,]^2
                if (cambia) sigma.tot[1,,]<-var(y)
                py.h[1,]<-w[1]*dmvnorm(y-t(H%*%t(x.z%*%Beta[1,,])),sigma=sigma.tot[1,,])  
 } else for (i in 1:k) {
                chsi[i,,]<-ginv(t(H)%*%ginv(psi)%*%H+ginv(sigma[i,,]))
                roy[i,,]<-chsi[i,,]%*%(t(H)%*%ginv(psi)%*%t(y)+ginv(sigma[i,,])%*%Beta[i,,]%*%t(x.z))
                if (r>1) {
                            temp<-(t(roy[i,,]))%o%(roy[i,,])
                            temp<-apply(temp,c(2,3),diag)
                            temp<-aperm(temp,c(2,3,1))                              
                            temp2<-array(chsi[i,,],c(r,r,numobs))                                            
                            Ezz.hy[i,,,]<-temp+temp2 } else Ezz.hy[i,r,r,]<-roy[i,,]^2+rep(chsi[i,,],numobs)
                            
                if (cambia[i]) sigma.tot[i,,]<-var(y)
                py.h[i,]<-(dmvnorm(y-t(H%*%Beta[i,,]%*%t(x.z)),sigma=sigma.tot[i,,]))
                ifelse(is.na(py.h[i,]),mean(py.h[i,],na.rm=TRUE),py.h[i,])
                if (sum(is.na(py.h[i,]))==numobs) py.h[i,]<-0
}
    
Ez.hy<-roy


#### E step determino E(z|y)        
Ez.y<-matrix(0,r,numobs) 
Ezz.y<-array(0,c(r,r,numobs))

for (i in 1:k) {
                ph.y[i,]<-w[i,]*py.h[i,]/diag(t(w)%*%py.h)
                ph.y[i,]<-ifelse(is.na(ph.y[i,]),mean(ph.y[i,],na.rm=TRUE),ph.y[i,]) 
                Ez.y<-Ez.y+t(matrix(ph.y[i,],numobs,r))*Ez.hy[i,,]
                Ezz.y<-Ezz.y+aperm(array(ph.y[i,],c(numobs,r,r)),c(2,3,1))*Ezz.hy[i,,,]}
           
Ez.y<-ifelse(is.na(Ez.y),rowMeans(Ez.y,na.rm=TRUE),Ez.y)      
Ezz.y<-ifelse(is.na(Ezz.y),rowMeans(Ezz.y,na.rm=TRUE),Ezz.y)    
                        
## MSTEP

EEzz.y<-rowMeans(Ezz.y, na.rm=TRUE, dims=2) 
H<-(t((y))%*%t(Ez.y)%*%ginv(EEzz.y))/numobs


psi<-(t(y)%*%y-t(y)%*%t(Ez.y)%*%t(H))/numobs
psi<-diag(diag(psi))


####stima dei parametri della mistura
ph.y<-ifelse(ph.y==0,0.0000000000001,ph.y) 

w<- matrix(t(rowMeans(ph.y,na.rm=TRUE)))
w<-w/sum(w)
w<-matrix(w,nrow=k,ncol=numobs)


if (q.w > 1) {
if (k>1) for (i in 2:k) {

                a<-phi[i,]
                for (g in 1:20) {A<-pi.greco.hess(a,i,ph.y,x.w,phi)
                                b<-pi.greco.grad(a,i,ph.y,x.w,phi)
                                a<-a-solve(A)%*%t(b)
                                a<-c(a)}
            phi[i,]<-a
            }

w<-(exp(phi%*%t(x.w))/colSums(exp(phi%*%t(x.w))))
}

for (i in 1:k) 
{Beta[i,,]<-(t(matrix(ph.y[i,],numobs,r))*Ez.hy[i,,])%*%x.z%*%ginv(t(x.z)%*%x.z)/rowMeans(t(matrix(ph.y[i,],numobs,r)))
                Beta[i,,]<-ifelse(is.na(Beta[i,,]),0,Beta[i,,])}


temp1<-matrix(0,r,r) 
temp2<-matrix(0,r,r) 
temp3<-matrix(0,r,1)


for (i in 1:k) {
                sigma[i,,]<-apply(aperm(array(ph.y[i,],c(numobs,r,r)),c(2,3,1))*(array(Ezz.hy[i,,,],c(r,r,numobs))-array(Beta[i,,]%*%t(x.z)%*%x.z%*%matrix(t(Beta[i,,]),nrow=q.z)/numobs,c(r,r,numobs))),1,rowMeans)/mean(ph.y[i,])
                temp1<-temp1+mean(w[i,])*sigma[i,,]
                dep<-Beta[i,,]%*%t(x.z)
                dep<-(dep%*%t(dep))/numobs
                temp2<-temp2+mean(w[i,])*(dep)
                temp3<-temp3+matrix(mean(w[i,])*rowMeans(Beta[i,,]%*%t(x.z)))
}


### identifiability correction
var.z<-temp1+temp2-temp3%*%t(temp3) 
A<-(chol(var.z)) 
for (i in 1:k) {sigma[i,,]<-t(ginv(A))%*%sigma[i,,]%*%ginv(A)
                Beta[i,,]<-t(ginv(A))%*%Beta[i,,]
}
                
H<-H%*%t(A)


             
temp<- sum(log(colSums(w*py.h))) 
if (is.infinite(temp)) temp<-2*lik
likelihood<-c(likelihood,temp) 
ratio<-abs((temp-lik)/lik) 
if ((temp < lik) & (hh > 20)) ratio<-eps 
lik<-temp
                                    }
out<-list(H=H,w=w,Beta=Beta,psi=psi,likelihood=likelihood,sigma=sigma,ph.y=ph.y,py.h=py.h,phi=phi) 
return(out) }
