separ <-
function(siga,vsiga,v,l,fracta,correa){
# separ() - SE's of parameters fract and corre
 
  vfract <- matrix(0, v+1, l, dimnames=dimnames(fracta))
  vcorre<- matrix(0, v+1, l^2, dimnames=dimnames(correa))

# components
  for(iv in 1: v) {   # component iv
    for(il in 1: l) {   # trait il
      # V(fract)
      ib <- (il-1)*l + il    # block no - col of siga for VCii
      vfract[iv,il] <- varz(
                       varlz(vsiga[(ib-1)*v+iv,(ib-1)*v+iv], siga[iv,(il-1)*l+il])
                + varlz(vt(v,l,vsiga,il,il), siga[v+1,(il-1)*l+il]) 
                - 2.0*covlyz(covcit(v,l,iv,vsiga,il,il), siga[iv,(il-1)*l+il],
                                                         siga[v+1,(il-1)*l+il]),
                                                         fracta[iv,il])

      # V(corre)
      covariance <- F  # assume component iv a var unless following tests detect otherwise
      if(rownames(siga)[iv] == "CovG(Ia,Ma)") {
        covariance <- T
        # work out the c1x & c2y nos for Ia and Ma - ie for VarG(Ia) & VarG(Ma)
        c1 <- match("VarG(Ia)", rownames(siga))
        c2 <- match("VarG(Ma)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Ma,Ia)") {
        covariance <- T
        c1 <- match("VarG(Ma)", rownames(siga))
        c2 <- match("VarG(Ia)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Id,Md)") {
        covariance <- T
        c1 <- match("VarG(Id)", rownames(siga))
        c2 <- match("VarG(Md)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Md,Id)") {
        covariance <- T
        c1 <- match("VarG(Md)", rownames(siga))
        c2 <- match("VarG(Id)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Ia:a,Ma:a)") {
        covariance <- T
        c1 <- match("VarG(Ia:a)", rownames(siga))
        c2 <- match("VarG(Ma:a)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Ma:a,Ia:a)") {
        covariance <- T
        c1 <- match("VarG(Ma:a)", rownames(siga))
        c2 <- match("VarG(Ia:a)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Ia:d,Ma:d)") {
        covariance <- T
        c1 <- match("VarG(Ia:d)", rownames(siga))
        c2 <- match("VarG(Ma:d)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Ma:d,Ia:d)") {
        covariance <- T
        c1 <- match("VarG(Ma:d)", rownames(siga))
        c2 <- match("VarG(Ia:d)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Id:d,Md:d)") {
        covariance <- T
        c1 <- match("VarG(Id:d)", rownames(siga))
        c2 <- match("VarG(Md:d)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovG(Md:d,Id:d)") {
        covariance <- T
        c1 <- match("VarG(Md:d)", rownames(siga))
        c2 <- match("VarG(Id:d)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovGs(Ia,Ma)") {
        covariance <- T
        c1 <- match("VarGs(Ia)", rownames(siga))
        c2 <- match("VarGs(Ma)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovGs(Ma,Ia)") {
        covariance <- T
        c1 <- match("VarGs(Ma)", rownames(siga))
        c2 <- match("VarGs(Ia)", rownames(siga))
      }

      if(rownames(siga)[iv] == "CovE(I,M)") {
        covariance <- T
        c1 <- match("VarE(I)", rownames(siga))
        c2 <- match("VarE(M)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovE(M,I)") {
        covariance <- T
        c1 <- match("VarE(M)", rownames(siga))
        c2 <- match("VarE(I)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovE(I,M&!C)") {
        covariance <- T
        c1 <- match("VarE(I)", rownames(siga))
        c2 <- match("VarE(M&!C)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovE(M&!C,I)") {
        covariance <- T
        c1 <- match("VarE(M&!C)", rownames(siga))
        c2 <- match("VarE(I)", rownames(siga))
      } 
      if(rownames(siga)[iv] == "CovE(I,M&C)") {
        covariance <- T
        c1 <- match("VarE(I)", rownames(siga))
        c2 <- match("VarE(M&C)", rownames(siga))
      }
      if(rownames(siga)[iv] == "CovE(M&C,I)") {
        covariance <- T
        c1 <- match("VarE(M&C)", rownames(siga))
        c2 <- match("VarE(I)", rownames(siga))
      } 

      for(jl in 1: l){    # trait jl
        jb <- (jl-1)*l + jl    # block no - col of siga for VCjj
        ijb <- (il-1)*l + jl    # off diag block no - col of siga for COVCij
        if(covariance) {
          if(is.na(c1) || is.na(c2)) {
            vcorre[iv, (il-1)*l+jl] <- NA
            cat("Need Var's as well as Cov to do SE of correlation:\n")
          }
          vcorre[iv, (il-1)*l+jl] <- varz(
             varlz(vsiga[(ijb-1)*v+iv,(ijb-1)*v+iv], siga[iv,(il-1)*l + jl])
            + 0.25*varlz(vsiga[(ib-1)*v+c1,(ib-1)*v+c1], siga[c1,(il-1)*l+il])
            + 0.25*varlz(vsiga[(jb-1)*v+c2,(jb-1)*v+c2], siga[c2,(jl-1)*l+jl])
            - covlyz(vsiga[(ijb-1)*v+iv,(ib-1)*v+c1], siga[iv,(il-1)*l+jl],
                                                    siga[c1,(il-1)*l+il])
            - covlyz(vsiga[(ijb-1)*v+iv,(jb-1)*v+c2], siga[iv,(il-1)*l+jl],
                                                    siga[c2,(jl-1)*l+jl])
            + 0.5*covlyz(vsiga[(ib-1)*v+c1,(jb-1)*v+c2], siga[c1,(il-1)*l+il],
                                                     siga[c2,(jl-1)*l+jl]),
                                                     correa[iv,(il-1)*l+jl])
        }
        else {
          vcorre[iv, (il-1)*l+jl] <- varz(
             varlz(vsiga[(ijb-1)*v+iv,(ijb-1)*v+iv], siga[iv,(il-1)*l + jl])
           + 0.25*varlz(vsiga[(ib-1)*v+iv,(ib-1)*v+iv], siga[iv,(il-1)*l+il])
           + 0.25*varlz(vsiga[(jb-1)*v+iv,(jb-1)*v+iv], siga[iv,(jl-1)*l+jl])
           - covlyz(vsiga[(ijb-1)*v+iv,(ib-1)*v+iv], siga[iv,(il-1)*l+jl],
                                                     siga[iv,(il-1)*l+il])
           - covlyz(vsiga[(ijb-1)*v+iv,(jb-1)*v+iv], siga[iv,(il-1)*l+jl], 
                                                     siga[iv,(jl-1)*l+jl])
           + 0.5*covlyz(vsiga[(ib-1)*v+iv,(jb-1)*v+iv], siga[iv,(il-1)*l+il],
                                                     siga[iv,(jl-1)*l+jl]),
                                                     correa[iv,(il-1)*l+jl])
        }
      }
    }
  }

# phenotypic
  for(il in 1:l) {
    vfract[v+1,il] <- 0.0  # total fract known exactly 1.0
    for(jl in 1:l) {
      vcorre[v+1,(il-1)*l+jl] <- varz(
           varlz(vt(v,l,vsiga,il,jl), siga[v+1,(il-1)*l+jl])
       + 0.25*varlz(vt(v,l,vsiga,il,il), siga[v+1,(il-1)*l+il])
       + 0.25*varlz(vt(v,l,vsiga,jl,jl), siga[v+1,(jl-1)*l+jl])
       - covlyz(covt(v,l,vsiga,il,jl,il), siga[v+1,(il-1)*l+jl], siga[v+1,(il-1)*l+il])
       - covlyz(covt(v,l,vsiga,il,jl,jl), siga[v+1,(il-1)*l+jl], siga[v+1,(jl-1)*l+jl])
       + 0.5*covlyz(covt(v,l,vsiga,il,il,jl), siga[v+1,(il-1)*l+il],
                                              siga[v+1,(jl-1)*l+jl]),
                                              correa[v+1,(il-1)*l+jl])
    }
  }

  outlist <- list(vfract=vfract, vcorre=vcorre)                 
  return(outlist)
}
