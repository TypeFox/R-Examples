"rUL_fc" <-
function( E=Z-XB(X,b) ) { 

## sample U, column mean of U, and L from full conditionals 
## needs pp_mu, pm_mu

  U<-UL$U ; L<-UL$L

  if(R>0) {

  mean_u<<-rnorm(R, 
  ( (1/var_u)*apply(UL$U,2,sum)+pp_mu*pm_mu ) / ( n/var_u+pp_mu ),
  1/sqrt(n/var_u + pp_mu ) 
                 )


  for(i in sample(1:n,size=n) ) {
    w<-E[i,-i]
    vr<-  solve( t(L)%*%t(U[-i,])%*%U[-i,]%*%L + diag(1/var_u,nrow=dim(U)[2]) )
    mn<-  vr%*%( t(L)%*%t(U[-i,])%*%w  + mean_u/var_u)
    U[i,]<- rmvnorm(mn,vr)        }

  vr<-solve(  ( (t(U)%*%U)^2 - t(U^2)%*%(U^2)) /2 + diag(pp_l,nrow=dim(U)[2])  )
  E[upper.tri(E,diag=TRUE)]<-0
  mn<-vr%*% diag(t(U)%*%E%*%U)
  L<-diag(c(rmvnorm(mn,vr)),nrow=dim(U)[2])
 
           }

  list(U=U,L=L)
             }

