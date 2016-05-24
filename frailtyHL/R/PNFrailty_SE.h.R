PNFrailty_SE.h <-
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
U<-res1[21][[1]]
H0<-res1[22][[1]]

################################################
######## SE for frailty parameter ###############
################################################
    n<-nrow(x)
    p<-ncol(x)
    u_h1<-exp(v_h1)
    oq<-matrix(1,qcum[nrand+1],1)
    oq1<-matrix(1,qcum[nrand+1],1)
    for (i in 1:nrand) {
       index1<-qcum[i]+1
       oq1[index1:qcum[i+1]]<-alpha_h1[i]
    } 
    D<-diag(oq1[,1])
    iD<-solve(D)
    iA<-iD
    Bi<-diag(clam0[,1])
    muh<-x%*%beta_h1 + z%*%v_h1
    expeta<-exp(muh)
    cla0<-di/(t(Mi)%*%expeta)
    temp4<-cla0^2/di
    As<-diag(temp4[,1])
    Wi<-diag(expeta[,1])
    done<-matrix(1,idx2,1)
    H22<-solve(t(z)%*%mat%*%z+U)
    Hessian<-matrix(0,nrand,nrand)
  if (varfixed==FALSE) {
    for (i in 1:nrand) {
        C<-matrix(0,qcum[nrand+1],qcum[nrand+1])
        index1<-qcum[i]+1
        index2<-qcum[i+1]
        for (j in index1:index2) C[j,j]<-1
        iB1<-iA%*%C%*%iA
        c_vh1<-iB1%*%v_h1
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
        ddk1<- -0.5*sum(diag(iA%*%C%*%iA%*%C)) +t(v_h1)%*%(iA%*%C%*%iB1)%*%v_h1 -t(dv1)%*%iB1%*%v_h1
        dia11<-(iB1%*%C%*%iA+iA%*%C%*%iB1)
        dv11<--H22%*%((t(z)%*%dmat1%*%z+dia1)%*%dv1-iB1%*%dv1 + dia11%*%v_h1)
        temp4<-(z%*%dv1)*(z%*%dv1)*expeta  +(z%*%dv11)*expeta
        ddW11<-diag(temp4[,1]) 
        ddcla011<- -( dAs1%*%(t(Mi)%*%Wi%*%z)%*%dv1 +As%*%(t(Mi)%*%dWi1%*%z)%*%dv1 +As%*%(t(Mi)%*%Wi%*%z)%*%dv11)  
        temp4<-Mi%*%diag(ddcla011[,1])%*%done
        ddB11<-diag(temp4[,1])
        temp4<-(2*(dcla01^2) + 2*(cla0*ddcla011) ) /di
        ddAs11<-diag(temp4[,1])
        ddm1_11<-(ddW11%*%Bi)+ (2*dWi1%*%dBi1) + (Wi%*%ddB11)
        ddm2_11<- (ddW11%*%Mi%*%As%*%t(Mi)%*%Wi) +(2*dWi1%*%Mi%*%dAs1%*%t(Mi)%*%Wi) +(2*dWi1%*%Mi%*%As%*%t(Mi)%*%dWi1)
        ddm3_11<-(Wi%*%Mi%*%ddAs11%*%t(Mi)%*%Wi) +(2*Wi%*%Mi%*%dAs1%*%t(Mi)%*%dWi1)  +(Wi%*%Mi%*%As%*%t(Mi)%*%ddW11)
        ddmat11<-ddm1_11-(ddm2_11+ddm3_11)
        Hd11 <-rbind(cbind(t(x)%*%ddmat11%*%x,t(x)%*%ddmat11%*%z),cbind(t(z)%*%ddmat11%*%x,t(z)%*%ddmat11%*%z+ dia11 ))
        ddk1<-ddk1 +0.5*sum(diag(-Hinv%*%Hd1%*%Hinv%*%Hd1+ Hinv%*%Hd11))
        Hessian[i,i]<-ddk1
        for (kk in 1:nrand) {
          if (kk>i) {
             D<-matrix(0,qcum[nrand+1],qcum[nrand+1])
             index1<-qcum[kk]+1
             index2<-qcum[kk+1]
             for (j in index1:index2) D[j,j]<-1
             iB2<-iA%*%D%*%iA
             c_vh2<-iB2%*%v_h1
             dv2<-H22%*%c_vh2
             dexpeta2<-expeta*(z%*%dv2)
             dcla02<--(di/((t(Mi)%*%expeta)^2))*(t(Mi)%*%dexpeta2)
             dWi2<-diag(dexpeta2[,1])
             dAi2<-diag(dcla02[,1])
             temp4<-Mi%*%dAi2%*%done
             dBi2<-diag(temp4[,1])
             dvec2<-2*(cla0*dcla02)
             temp4<-dvec2/di
             dAs2<-diag(temp4[,1])
             dd12<--0.5*sum(diag(iA%*%D%*%iA%*%C))+0.5*t(v_h1)%*%(iA%*%D%*%iB1+iA%*%C%*%iB2)%*%v_h1-t(dv1)%*%iB2%*%v_h1
             dia12<-(iB1%*%D%*%iA+iA%*%D%*%iB1)
             dv12<--H22%*%((t(z)%*%dmat1%*%z+dia1)%*%dv2-iB2%*%dv1 + dia12%*%v_h1)
             temp4<-(z%*%dv1)*(z%*%dv2)*expeta + (z%*%dv12)*expeta
             ddW12<-diag(temp4[,1])
             ddcla012<--( dAs2%*%(t(Mi)%*%Wi%*%z)%*%dv1 +As%*%(t(Mi)%*%dWi2%*%z)%*%dv1 +As%*%(t(Mi)%*%Wi%*%z)%*%dv12)  
             temp4<-Mi%*%diag(ddcla012[,1])%*%done
             ddB12<-diag(temp4[,1])
             temp4<-(2*(dcla02*dcla01) + 2*(cla0*ddcla012) )/di
             ddAs12<-diag(temp4[,1])
             ddm1_12<-(ddW12%*%Bi)+ (dWi1%*%dBi2 + dWi2%*%dBi1) + (Wi%*%ddB12)
             ddm2_12<-(ddW12%*%Mi%*%As%*%t(Mi)%*%Wi) +(dWi1%*%Mi%*%dAs2%*%t(Mi)%*%Wi) +(dWi1%*%Mi%*%As%*%t(Mi)%*%dWi2)
             ddm3_12<-(dWi2%*%Mi%*%dAs1%*%t(Mi)%*%Wi) +(Wi%*%Mi%*%ddAs12%*%t(Mi)%*%Wi)  +(Wi%*%Mi%*%dAs1%*%t(Mi)%*%dWi2)
             ddm4_12<-(dWi2%*%Mi%*%As%*%t(Mi)%*%dWi1) +(Wi%*%Mi%*%dAs2%*%t(Mi)%*%dWi1)  +(Wi%*%Mi%*%As%*%t(Mi)%*%ddW12)
             ddmat12<-ddm1_12-(ddm2_12+ddm3_12+ddm4_12)
             dmat2<-(dWi2%*%Bi)+(Wi%*%dBi2)-(dWi2%*%Mi%*%As%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%dAs2%*%t(Mi)%*%Wi)-(Wi%*%Mi%*%As%*%t(Mi)%*%dWi2)
             dia2<--iB2
             Hd2 <- rbind(cbind(t(x)%*%dmat2%*%x,t(x)%*%dmat2%*%z),cbind(t(z)%*%dmat2%*%x,t(z)%*%dmat2%*%z+ dia2))
             Hd12<-rbind(cbind(t(x)%*%ddmat12%*%x,t(x)%*%ddmat12%*%z),cbind(t(z)%*%ddmat12%*%x,t(z)%*%ddmat12%*%z+ dia12))
             dd12<-dd12 +0.5*sum(diag(-Hinv%*%Hd1%*%Hinv%*%Hd2+ Hinv%*%Hd12))
             Hessian[i,kk]<-Hessian[kk,i]<-dd12
         }
       }
     }
     iAp<-solve(Hessian)
     se_lam<-sqrt(diag(iAp))
   }
     if (varfixed==TRUE) {
        se_lam<-rep(0,nrand)
        for (i in 1:nrand) se_lam[i]<-"NULL"
     }
     eta<-x%*%beta_h1 + z%*%v_h1
     expeta<-exp(eta)
     one<-matrix(1,n,1)
     done<-matrix(1,idx2,1)
     oq<-matrix(1,qcum[nrand+1],1)
     pi<-3.14159265359
     term0<-t(Mi)%*%expeta
     hlike1<-(t(one)%*%(del*eta) )-( t(done)%*%(di*log(term0)))
     hlike2<-0
     hlike3<-0
     for (i in 1:nrand) {
         if (alpha_h[i]>1e-05) hlike2<-hlike2-(q[i]/2)*log(2*pi)-( (q[i]/2)*log(alpha_h1[i]))
         index1<-qcum[i]+1
         index2<-qcum[i+1]
         vv_h1<-matrix(0,q[i],1)
         vv_h1[1:q[i],1]<-v_h1[index1:index2,1]
         if (alpha_h[i]>1e-05) hlike3<-hlike3-(t(vv_h1)%*%vv_h1)/(2*alpha_h1[i])
     }
     hliken<-hlike1+hlike2+hlike3
     cc1<-svd(2*pi*Hinv)
     for (i in 1:length(cc1$d)) if (cc1$d[i]<0.00001) cc1$d[i]<-1
     logdet1<-sum(log(abs(cc1$d)))
##     adj1<- 0.5*(p+qcum[nrand+1])*log(2*pi)+0.5*logdet1
     adj1<- 0.5*logdet1
     hpn1<-hliken+ adj1
     muu<-exp(x%*%beta_h1)*clam0
     zmu<-t(z)%*%muu
     u_h1<-exp(v_h1)
     second<-0
     for (i in 1:nrand) {
         ialph1<-1/alpha_h1[i]
         a21<-(zmu*u_h1)+ialph1
         b31<-zmu*u_h1
         S11<-3*(b31/(a21^2))
         S21<-5*((b31^2)/(a21^3))
         temp4<-S11-S21
         S31<-diag(temp4[,1])
         second<-second-sum(diag(S31))/24
    }
    H22<-t(z)%*%mat%*%z+iD
    cc1<-svd(H22/(2*pi))
    for (i in 1:length(cc1$d)) if (cc1$d[i]>100000) cc1$d[i]<-1
    logdet1<-sum(log(abs(cc1$d)))
##    hpn2<-hliken-0.5*log(det(H22/(2*pi)))
    hpn2<-hliken-0.5*logdet1
    hpn3<-hpn1+second
    df1<-sum(diag(Hinv%*%H0))
    res<-list(se_lam,hliken,hpn1,hpn2,hpn3,hlike1,df1)
    return(res)
}

