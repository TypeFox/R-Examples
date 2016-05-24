PNfrailtyHL <-
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
    dft1<-t(x)%*%(del-Wi%*%clam0)
######################## pv(hp) for beta ########################
    if (mord==1) {
    U <- iD
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
    dft2<-t(z)%*%(del-Wi%*%clam0)-(iD%*%v_h0)
    dft<-rbind(dft1,dft2)
    Bi<-diag(clam0[,1])
    temp4<-cla0^2/di
    Adi<-diag(temp4[,1])
    As<-Adi
    mat<-(Wi%*%Bi)-(Wi%*%Mi)%*%Adi%*%(t(Mi)%*%Wi)
    U<-iD
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
    if (dord==0) {
        index1<-p+qcum[i]+1
        index2<-p+qcum[i+1]
        gamma<-sum(diag(Hinv[index1:index2,index1:index2]))/alpha_h[i]
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        if (varfixed==FALSE) alpha_h[i]<-sum(v_h[index1:index2,1]^2)/(q[i]-gamma)
    }
    if (dord==1 | dord==2) {
        H22<-solve(t(z)%*%mat%*%z+U)
        ial1<-1/alpha_h[i]
        iA<-iD
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iB1<-iA%*%C%*%iA
        c_vh1<-iB1%*%v_h
        dv1<-H22%*%c_vh1
        dexpeta1<-expeta*(z%*%dv1)
        dcla01<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta1)
        dWi1<-diag(dexpeta1[,1])
        dAi1<-diag(dcla01[,1])
        temp4<-Mi%*%dAi1%*%done
        dBi1<-diag(temp4[,1])
        dvec1<-2*(cla0*dcla01)
        temp4<-dvec1/di
        dAs1<-diag(temp4[,1])
        dmat1<-(dWi1%*%Bi)+(Wi%*%dBi1)-(dWi1%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs1%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi1)
        dia1<--iB1
        Hd1 <- rbind(cbind(t(x)%*%dmat1%*%x,t(x)%*%dmat1%*%z),cbind(t(z)%*%dmat1%*%x,t(z)%*%dmat1%*%z+ dia1))
        gamma1<--alpha_h[i]*sum(diag((Hinv%*%Hd1)))
        if (dord==2) {
        expeta <- exp(x%*%beta_h + z%*%v_h)
        ial1<-1/alpha_h[i]
        muu<-exp(x%*%beta_h)*clam0
        zmuu<-t(z)%*%muu
        u_h<-exp(v_h)
        ude1<-u_h*dv1
        aa1<-(zmuu*u_h)+ial1
        bb1<-(zmuu*ude1)-(ial1^2)
        term11<-((aa1*zmuu*ude1)-(2*zmuu*u_h*bb1))/(aa1^3)
        term21<-((2*aa1*zmuu*zmuu*u_h*ude1)-(3*((zmuu*u_h)^2)*bb1))/(aa1^4)
        term1<-(3*term11)-(5*term21)
        SS1<-diag(term1[,1])
        gamma21<--(alpha_h[i]/12)*sum(diag(SS1))
        }
        if (dord==1) {
            gamma21<-0
        }
        k21<- q[i]- gamma1- gamma21-sum(v_h[index1:index2,1]^2)/(alpha_h[i])
        if (varfixed==FALSE) alpha_h[i]<-sum(v_h[index1:index2,1]^2)/(q[i]-gamma1-gamma21)
    }
    }
    res<-list(x,z,y,del,Mi,idx2,t2,di,beta_h0,v_h0,beta_h,v_h,alpha_h0,alpha_h,dft,Hinv,clam0,H,mat,se_beta_h,U,H0)
    return(res)
}

