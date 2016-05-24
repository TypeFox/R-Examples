PGFrailty_SE.h <-
function(res1,nrand,q,qcum,dord=1,varfixed=FALSE) {
x<-res1[1][[1]]
z<-res1[2][[1]]
y<-res1[3][[1]]
del<-res1[4][[1]]
Mi<-res1[5][[1]]
idx2<-res1[6][[1]]
t2<-res1[7][[1]]
di<-res1[8][[1]]
beta_h<-res1[9][[1]]
v_h<-res1[10][[1]]
beta_h1<-res1[11][[1]]
v_h1<-res1[12][[1]]
alpha_h<-res1[13][[1]]
alpha_h1<-res1[14][[1]]
dft<-res1[15][[1]]
Hinv<-res1[16][[1]]
clam0<-res1[17][[1]]
H<-res1[18][[1]]
mat<-res1[19][[1]]
H0<-res1[22][[1]]

################################################
######## SE for frailty parameter ###############
################################################
    n<-nrow(x)
    p<-ncol(x)
    u_h1<-exp(v_h1)
    mat11<-t(x)%*%mat%*%x
    mat12<-t(x)%*%mat%*%z
    mat13<-t(z)%*%mat%*%z
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h1[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    U <- iD%*%diag(u_h1[,1])
    mmat<-mat11-mat12%*%solve(mat13+U)%*%t(mat12)
    hminv<-solve(mmat)
    done<-matrix(1,idx2,1)
    muh<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(muh)
    Wi<-diag(expeta[,1])
    cla0<-di/(t(Mi)%*%expeta)
    Ai<-diag(cla0[,1])
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    Bi<-diag(clam0[,1])
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)
    Dinv0<-solve(t(z)%*%mat%*%z+U)
    se_lam<-matrix(0,nrand,1)
  if (varfixed==FALSE) {
    for (i in 1:nrand){
    ial<-1/alpha_h1[i]
    index1<-qcum[i]+1
    index2<-qcum[i+1]
    vv_h1<-matrix(0,q[i],1)
    uu_h1<-matrix(0,q[i],1)
    vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
    uu_h1[1:q[i],1]<-u_h1[index1:index2,1]
    c_vh <- ial^2*(uu_h1-1)
    dv<-solve(t(z)%*%mat%*%z+U)%*%c_vh
    dexpeta<-expeta*(z%*%dv)
    dcla0<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta)
    dWi<-diag(dexpeta[,1])
    dAi<-diag(dcla0[,1])
    temp4<-Mi%*%dAi%*%done
    dBi<-diag(temp4[,1])
    dvec<-2*(cla0*dcla0)
    temp4<-dvec/di
    dAs<-diag(temp4[,1])
    dmat<-(dWi%*%Bi)+(Wi%*%dBi)-(dWi%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi)
    uad1<-uu_h1*dv
    temp4<-(-ial^2*uu_h1)+(ial*uad1) 
    dia_d<-diag(temp4[,1])
    Hd <- rbind(cbind(t(x)%*%dmat%*%x,t(x)%*%dmat%*%z),cbind(t(z)%*%dmat%*%x,t(z)%*%dmat%*%z+ dia_d))
    term<-(t(z)%*%mat%*%z+U)
    invt<-solve(term)
    dterm<-t(z)%*%dmat%*%z + dia_d
    ddv<- - invt%*%(dterm)%*%invt%*%c_vh+ invt%*%( -2*ial^3*(uu_h1-1)+ ial^2*uad1 )
    ddcla0<- -dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
    ddAi<-diag(ddcla0[,1])
    ddclam0<-Mi%*%ddAi%*%done
    ddBi<-diag(ddclam0[,1])
    temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
    ddWi<-diag(temp4[,1])
    ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
    ddm2<- ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
    ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
    temp4<-ddvec/di
    ddAs<-diag(temp4[,1])
    ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
    ddmat<-ddm1-(ddm2+ddm3)
    uad2<-(uad1*dv)+(u_h1*ddv)
    temp4<- (2*ial^3*u_h1)-(2*ial^2*uad1)+(ial*uad2) 
    dia_dd<-diag(temp4[,1])
    Hdd <-rbind(cbind(t(x)%*%ddmat%*%x,t(x)%*%ddmat%*%z),cbind(t(z)%*%ddmat%*%x,t(z)%*%ddmat%*%z+ dia_dd))
    H <- rbind(cbind(t(x)%*%mat%*%x,t(x)%*%mat%*%z),cbind(t(z)%*%mat%*%x,t(z)%*%mat%*%z+U))
    Hinv<-solve(H)
    oq<-matrix(1,q[i],1)
    dp<-digamma(ial)
    ddp<-trigamma(ial)
    k21a<-t(oq)%*%(2*(vv_h1 -uu_h1))+( q[i]*( (-2*log(alpha_h1[i]))+3 -(2*dp)-((1/alpha_h1[i])*ddp))  )
    k21a<-(ial^3)*k21a
    d2halp<--k21a- t(oq)%*%(c_vh*dv)
    adalp<-0.5*sum(diag(-Hinv%*%Hd%*%Hinv%*%Hd+ Hinv%*%Hdd))
    dalp_2<-d2halp+adalp
    se_lam[i]<-sqrt(1/dalp_2)
    }
   }
  if (varfixed==TRUE) {
    for (i in 1:nrand) se_lam[i]<-"NULL"
  }
    u_h1<-exp(v_h1)
    U <- iD%*%diag(u_h1[,1])
    oq<-matrix(1,qcum[nrand+1],1)
    one<-matrix(1,n,1)
    zd<-t(z)%*%del
    eta<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(eta)
    term0<-t(Mi)%*%expeta
    non<-t(one)%*%del
    done<-matrix(1,idx2,1)
    hlike1<-(t(one)%*%(del*eta) )-( t(done)%*%(di*log(term0)) )
    hlike2<-0
    for (i in 1:nrand) {
       oq<-matrix(1,q[i],1)
       index1<-qcum[i]+1
       index2<-qcum[i+1]
       vv_h1<-matrix(0,q[i],1)
       uu_h1<-matrix(0,q[i],1)
       vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
       uu_h1[1:q[i],1]<-u_h1[index1:index2,1]
       i_alp1<-1/alpha_h1[i]
       c_alp1<-0
       if (alpha_h[i]>1e-05) c_alp1<--log(gamma(i_alp1))-(i_alp1*log(alpha_h1[i]))
       if (alpha_h[i]>1e-05) hlike2<- hlike2+t(oq)%*%( (vv_h1-uu_h1)/alpha_h1[i] + c_alp1)
    }
    hlike<-hlike1+hlike2
    pi<-3.14159265359
    H22<-t(z)%*%mat%*%z+U
    zd<-t(z)%*%del
    if (dord==2) {
       temp4<-1/(i_alp1+zd)
       secd<-diag(temp4[,1])
       second<-sum(diag(secd))/12
    } else second<-0
    cc1<-svd(H22/(2*pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i]>100000) cc1$d[i]<-1
    logdet1<-sum(log(abs(cc1$d)))
##    pvhs<-hlike-0.5*log(det(H22/(2*pi)))
    pvhs<-hlike-0.5*logdet1
    svhs<-pvhs+second
    cc1<-svd(Hinv*2*pi)
    for (i in 1:length(cc1$d)) if (cc1$d[i]<0.00001) cc1$d[i]<-1
##    adj1<-( (0.5*(p+qcum[nrand+1]))*log(2*pi)) + (0.5*log(det(Hinv)) )
    logdet1<-sum(log(abs(cc1$d)))
    adj1<-0.5*logdet1
    hpn1<-hlike+ adj1
    hpn2<-pvhs
    hpn3<-hpn1+second
    df1<-sum(diag(Hinv%*%H0))
    res<-list(se_lam,hlike,hpn1,hpn2,hpn3,hlike1,df1)
    return(res)
}

