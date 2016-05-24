PGfrailtyHL <-
function(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,alpha_h0,mord,dord,varfixed=FALSE){
    n<-nrow(x)
    p<-ncol(x)
    nrand <- length(z)
    q <- rep(0, nrand)
    for (i in 1:nrand) q[i] <- dim(z[[i]])[2]
    qcum <- cumsum(c(0, q))
    beta_h<-beta_h0
    v_h<-v_h0
    for (i in 1:nrand) {
       if (alpha_h0[i]<0.000001) alpha_h0[i]<-0.000001
    }
    alpha_h<-alpha_h0
    zz<-z[[1]]
    if (nrand>1) {
        index1<-nrand
        for (i in 2:index1) zz<-cbind(zz, z[[i]])
    }
    z<-matrix(0,n,qcum[nrand+1])
    z[1:n,1:qcum[nrand+1]]<-zz[1:n,1:qcum[nrand+1]]
    muh<-x%*%beta_h0 + z%*%v_h0
    expeta<-exp(muh)
    cla0<-di/(t(Mi)%*%expeta)
    Wi<-diag(expeta[,1])
    Ai<-diag(cla0[,1])
    done<-matrix(1,idx2,1)

########################
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    clam0<-Mi%*%Ai%*%done
    
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    iu_h0<-exp(v_h0) ## gamma frailty
    U<-iD%*%diag(iu_h0[,1]) ## gamma frailty
    dft1<-t(x)%*%(del-Wi%*%clam0)
######################## pv(hp) for beta ########################
    if (mord==1) {
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%As%*%(t(Mi)%*%Wi)
    dinv<-solve(t(z)%*%mat%*%z+U)
    mu0<-exp(x%*%beta_h0 + z%*%v_h0)
    mu<-exp(x%*%beta_h0 + z%*%v_h0)*clam0
    dcla0be<--As%*%(t(Mi)%*%Wi)
    dcla0b<-matrix(0,idx2,p)
    dv_db<-matrix(0,qcum[nrand+1],p)
    xz<-matrix(0,n,p)
    dmu0<-matrix(0,n,p)
    dw_db1<-matrix(0,n,p)
    xk<-matrix(0,n,1)
    ad1<-matrix(0,p,1)
    for (k in 1:p) {
       xk[,1]<-x[,k]
       dv_db[,k] <--dinv%*%(t(z)%*%mat%*%xk)
       xz[,k]<-xk+z%*%(dv_db[,k])
       dcla0b[,k]<-dcla0be%*%(xz[,k])
       dc<-Mi%*%diag(dcla0b[,k])%*%done
       dmu0[,k]<-mu0*(xz[,k])
       dw_db1[,k]<- dc*mu0 + mu*xz[,k]
       temp4<-(2*cla0/di)*(dcla0b[,k])
       dw_db2<-((diag(dmu0[,k]))%*%Mi%*%As%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%diag(temp4[,1])%*%t(Mi)%*%Wi)+(Wi%*%Mi%*%As%*%t(Mi)%*%(diag(dmu0[,k]))) 
       dw_db<-diag(dw_db1[,k])-dw_db2
       ad1[k,1]<-sum(diag(dinv%*%(t(z)%*%dw_db%*%z))) 
    }
    dft1<-dft1-0.5*ad1
   }
#########################################################
    dft2<-t(z)%*%(del-Wi%*%clam0)+(iD%*%oq)-(iD%*%iu_h0) ## gamma frailty
    dft<-rbind(dft1,dft2)
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    Adi<-diag(temp4[,1])
    As<-Adi
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%Adi%*%(t(Mi)%*%Wi)
    U<-iD%*%diag(iu_h0[,1]) ## gamma frailty
    H <- rbind(cbind(t(x)%*%mat%*%x, t(x)%*%mat%*%z), cbind(t(z)%*%mat%*%x, t(z)%*%mat%*%z+U))
    H0 <- rbind(cbind(t(x)%*%mat%*%x, t(x)%*%mat%*%z), cbind(t(z)%*%mat%*%x, t(z)%*%mat%*%z))
    Hinv<-solve(H)
    be_h0<- rbind(beta_h0, v_h0)
    be_h <- be_h0 + (Hinv%*%dft)
    beta_h[1:p,1]<-be_h[1:p,1]
    se_beta_h<-matrix(0,p,1)
    for (i in 1:p) se_beta_h[i,1]<-sqrt(Hinv[i,i])
    index2<-qcum[nrand+1]
    index3<-p+1
    index4<-p+qcum[nrand+1]
    v_h[1:index2,1]<-be_h[index3:index4,1]
################################################
    for (i in 1:nrand) {
        ial<-1/alpha_h[i]
        dp<-digamma(ial)
        ddp<-trigamma(ial)
        oq<-matrix(1,q[i],1)
        one<-matrix(1,n,1)
        u_h<-exp(v_h) 
        eta<-x%*%beta_h + z%*%v_h
        expeta<-exp(eta)
        Wi<-diag(expeta[,1])
        Wei<-(Wi%*%Bi)
        U<-iD%*%diag(u_h[,1])
        term<-(t(z)%*%mat%*%z+U)
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iA<-iD
        iB<-iA%*%C%*%iA
        c_vh<- iB%*%(u_h-1) ## gamma frailty
        invt<-solve(term)
        dv<-invt%*%c_vh
        dexpeta<-expeta*(z%*%dv)
        dcla0<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta)
        temp4<-expeta*(z%*%dv)
        dWi<-diag(temp4[,1])
        dAi<-diag(dcla0[,1])
        temp4<-Mi%*%dAi%*%done
        dBi<-diag(temp4[,1])
        dvec<-2*(cla0*dcla0)
        temp4<-dvec/di
        dAs<-diag(temp4[,1])
        dmat<-(dWi%*%Bi)+(Wi%*%dBi)-(dWi%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi)
        uad1<-u_h*dv
        temp4<--iB%*%u_h+ial*uad1
        dia1<-diag(temp4[,1])
        Hd <- rbind(cbind(t(x)%*%dmat%*%x,t(x)%*%dmat%*%z),cbind(t(z)%*%dmat%*%x,t(z)%*%dmat%*%z+ dia1))
    if (dord==0) {
        temp4<--iD%*%iD%*%u_h
        dia1_k<-diag(temp4[,1])
        zero1<-matrix(0,p,p)
        zero2<-matrix(0,p,qcum[nrand+1])
        Hd_k<-rbind(cbind(zero1,zero2),cbind(t(zero2),dia1_k))
        hinv2<-solve(t(z)%*%mat%*%z+U)
        Hd2<-t(z)%*%dmat%*%z+dia1
        dk2<--0.5*sum(diag(Hinv%*%Hd_k))
        vv_h<-matrix(0,q[i],1)
        index3<-1
        index4<-q[i]
        vv_h[1:q[i],1]<-v_h[index1:index2,1]
        uu_h<-exp(vv_h)
        k2<-(t(oq)%*%(vv_h -uu_h))+( q[i]*(-log(alpha_h[i])+1 - dp) )
        mu<-exp(x%*%beta_h)*clam0
        zd<-t(z)%*%del
        zmu<-t(z)%*%mu
        cor1<-(1+(alpha_h[i]*zd))^(-2)
        dk3<-sum(cor1)/12
        k2<--((alpha_h[i]^-2)*k2)+dk2
        dterm<-t(z)%*%dmat%*%z+dia1
        ddv<-- invt%*%(dterm)%*%invt%*%c_vh+invt%*%(-2*ial^3*(u_h-1)+ial^2*uad1 )
        ddcla0<--dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
        uad2<-(uad1*dv)+(u_h*ddv)
        ddAi<-diag(ddcla0[,1])
        ddclam0<-Mi%*%ddAi%*%done
        ddBi<-diag(ddclam0[,1])
        temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
        ddWi<-diag(temp4[,1])
        ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
        ddm2<-ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
        ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
        temp4<-ddvec/di
        ddAs=diag(temp4[,1])
        ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
        ddmat<-ddm1-(ddm2+ddm3)
        temp4<-(2*ial^3*u_h)-(2*ial^2*uad1)+(ial*uad2)
        dia2<-diag(temp4[,1])
        Hdd2<-t(z)%*%ddmat%*%z+dia2
        k21<-t(oq)%*%(2*(vv_h -uu_h))+ ( q[i]*( (-2*log(alpha_h[i]))+3 - (2*dp)-((1/alpha_h[i])*ddp))  )
        al<-ial^3
        k21<-al*k21
        cor2<-(1+(alpha_h[i]*zd))^(-3)
        cor22<-cor2*zd
        kcor<-sum(cor22)/6
        kcor<-0
        k22<--k21+0.5*(sum(diag(hinv2*Hdd2))-sum(diag(hinv2*Hd2*hinv2*Hd2)))+kcor
        ialp<-1/k22
        if (varfixed==FALSE) alpha_h[i]<-alpha_h[i] + (ialp*k2)
    }
    if (dord==1 | dord==2) {
        hinv2<-solve(t(z)%*%mat%*%z+U)
        Hd2<-t(z)%*%dmat%*%z+dia1
        dk2<--0.5*sum(diag(Hinv%*%Hd))
        vv_h<-matrix(0,q[i],1)
        index3<-1
        index4<-q[i]
        vv_h[1:index4,1]<-v_h[index1:index2,1]
        uu_h<-exp(vv_h)
        k2<-(t(oq)%*%(vv_h -uu_h))+( q[i]*(-log(alpha_h[i])+1 - dp) )
        mu<-exp(x%*%beta_h)*clam0
        zd<-t(z)%*%del
        zmu<-t(z)%*%mu
        cor1<-(1+(alpha_h[i]*zd))^(-2)
        dk3<-sum(cor1)/12
        if(dord==1) k2<--((alpha_h[i]^(-2))*k2)+dk2
        if(dord==2) k2<--((alpha_h[i]^(-2))*k2)+dk2 +dk3
        dterm<-t(z)%*%dmat%*%z+dia1
        ddv<--invt%*%(dterm)%*%invt%*%c_vh+ invt%*%( -2*ial^3*(u_h-1)+ial^2*uad1 )
        ddcla0<--dAs%*%(t(Mi)%*%Wi%*%z)%*%dv-As%*%(t(Mi)%*%dWi%*%z)%*%dv-As%*%(t(Mi)%*%Wi%*%z)%*%ddv 
        uad2<-(uad1*dv)+(u_h*ddv)
        ddAi<-diag(ddcla0[,1])
        ddclam0<-Mi%*%ddAi%*%done
        ddBi<-diag(ddclam0[,1])
        temp4<-expeta*(z%*%dv)*(z%*%dv)+(expeta*(z%*%ddv))
        ddWi<-diag(temp4[,1])
        ddm1<-(ddWi%*%Bi)+(2*dWi%*%dBi) + (Wi%*%ddBi)
        ddm2<-ddWi%*%Mi%*%As%*%t(Mi)%*%Wi+(2*dWi%*%Mi%*%dAs%*%t(Mi)%*%Wi)+(2*dWi%*%Mi%*%As%*%t(Mi)%*%dWi)
        ddvec<-(2*(dcla0^2)) + (2*(cla0*ddcla0))
        temp4<-ddvec/di
        ddAs<-diag(temp4[,1])
        ddm3<-(Wi%*%Mi%*%ddAs%*%t(Mi)%*%Wi)+(2*Wi%*%Mi%*%dAs%*%t(Mi)%*%dWi)+Wi%*%Mi%*%As%*%t(Mi)%*%ddWi
        ddmat<-ddm1-(ddm2+ddm3)
        temp4<-(2*ial^3*u_h)-(2*ial^2*uad1)+(ial*uad2)
        dia2<-diag(temp4[,1])
        Hdd<-rbind(cbind(t(x)%*%ddmat%*%x,t(x)%*%ddmat%*%z),cbind(t(z)%*%ddmat%*%x,t(z)%*%ddmat%*%z+dia2))
        Hdd2<-t(z)%*%ddmat%*%z+dia2
        k21<-t(oq)%*%(2*(vv_h -uu_h))+(q[i]*((-2*log(alpha_h[i]))+3-(2*dp)-((1/alpha_h[i])*ddp)))
        al<-ial^3
        k21<-al*k21
        cor2<-(1+(alpha_h[i]*zd))^(-3)
        cor22<-cor2*zd
        kcor<-sum(cor22)/6
        if(dord==1) kcor<-0
        k22<--k21+0.5*(sum(diag(Hinv%*%Hdd))-sum(diag(Hinv%*%Hd%*%Hinv%*%Hd)))+kcor
        ialp<-1/k22
        if (varfixed==FALSE) alpha_h[i] <- alpha_h[i] + (ialp*k2)
        if (alpha_h[i]<=0.0) alpha_h[i]<-alpha_h0/2
    }
    }
    res<-list(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,beta_h,v_h,alpha_h0,alpha_h,dft,Hinv,clam0,H,mat,se_beta_h,U,H0)
    return(res)
}

