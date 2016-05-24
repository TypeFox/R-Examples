comtopar <-
function(v, l, siga,vara, vsiga, sesiga, ctable){
# comtopar() - calls covtopar for each animal var/cov component
#            - to get hsq (fraction) and correlation for component
#	     - if component is a covariance, correlation uses appropriate variance
#              components in denominator - these variance components must be fitted
#            - sigt is total animal var/cov - ie phenotypic var/cov matrix - l x l
#            - siga is all v components var/cov's - v x l^2
#            - also computes Rp - phenotypic correlation matrix
#            - also computes SE's
#            - v = no of components, l = no of traits
#            - in outlist - labelc is vector of component names - v+1
#                         - correa is component correlations - v+1 x l^2
#                         - fracta is component variance fractions - v+1 x l

# adjust vsiga to ensure posdef
# vsiga <- nearPD(vsiga,ensureSymmetry=T)$mat

# calculate sigt - sum of components
  sigt <- matrix(0,l,l, dimnames=dimnames(vara))
  sigt <- matrix(apply(siga,2,sum),l,l,dimnames=dimnames(vara))
# sampling variances and SE's for sigt
  vsigt <- matrix(0,l,l,dimnames=dimnames(vara))
  for (il in 1:l) {
    ib <- (il-1)*l+il
    for(jl in 1:l) {
      jb <- (jl-1)*l+jl
#     vsigt[ib,jb] <- vt(v,l,vsiga,il,jl)
      vsigt[il,jl] <- vt(v,l,vsiga,il,jl)
    }
  }
  sesigt <- sqrt(vsigt)

# components
  correa <- matrix(0,v+1,l^2, dimnames=list(c(rownames(siga),"VarP(I)"), colnames(siga)))
  fracta <- matrix(0,v+1,l,dimnames=list(c(rownames(siga),"VarP(I)"), colnames(sigt)))
  varcomp <- matrix(0,v+1,l^2,dimnames=list(c(rownames(siga),"VarP(I)"), colnames(siga)))
  sevarcomp <- matrix(0,v+1,l^2,dimnames=list(c(rownames(siga),"VarP(I)"), colnames(siga)))
  labelc <- rep(NULL,v+1)
  for(i in 1:v) {
    covi <- matrix(siga[i, ], l, l)
    varcomp[i, ] <- siga[i, ]
    sevarcomp[i, ] <- sesiga[i, ]

    if(any(!is.na(match(rownames(siga)[i],ctable$allvar)))){
      # variance component - use covtopar()
      pari <- covtopar(covi, sigt)
    }
    else {
      # covariance component - use crosseffectcovtopar()
      if(rownames(siga)[i] == "CovG(Ia,Ma)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ia)","VarG(Ma)")
      }
      if(rownames(siga)[i] == "CovG(Ma,Ia)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ma)","VarG(Ia)")
      }
      if(rownames(siga)[i] == "CovG(Id,Md)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Id)","VarG(Md)")
      }
      if(rownames(siga)[i] == "CovG(Md,Id)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Md)","VarG(Id)")
      }
      if(rownames(siga)[i] == "CovG(Ia:a,Ma:a)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ia:a)","VarG(Ma:a)")
      }
      if(rownames(siga)[i] == "CovG(Ma:a,Ia:a)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ma:a)","VarG(Ia:a)")
      }
      if(rownames(siga)[i] == "CovG(Ia:d,Ma:d)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ia:d)","VarG(Ma:d)")
      }
      if(rownames(siga)[i] == "CovG(Ma:d,Ia:d)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Ma:d)","VarG(Ia:d)")
      }
      if(rownames(siga)[i] == "CovG(Id:d,Md:d)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Id:d)","VarG(Md:d)")
      }
      if(rownames(siga)[i] == "CovG(Md:d,Id:d)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarG(Md:d)","VarG(Id:d)")
      }
      if(rownames(siga)[i] == "CovGs(Ia,Ma)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarGs(Ia)","VarGs(Ma)")
      }
      if(rownames(siga)[i] == "CovGs(Ma,Ia)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarGs(Ma)","VarGs(Ia)")
      }

      if(rownames(siga)[i] == "CovE(I,M)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(I)","VarE(M)")
      }
      if(rownames(siga)[i] == "CovE(M,I)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(M)","VarE(I)")
      }
      if(rownames(siga)[i] == "CovE(I,M&!C)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(I)","VarE(M&!C)")
      }
      if(rownames(siga)[i] == "CovE(M&!C,I)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(M&!C)","VarE(I)")
      }
      if(rownames(siga)[i] == "CovE(I,M&C)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(I)","VarE(M&C)")
      }
      if(rownames(siga)[i] == "CovE(M&C,I)") {
        pari <- crosseffectcovtopar(covi,sigt,varcomp,"VarE(M&C)","VarE(I)")
      }
    }
    correa [i, ] <- pari$corre
    fracta[i, ] <- pari$fract
    labelc[i] <- dimnames(siga)[[1]][i]
  }
  
# phenotypic variance component
  pari <- covtopar(sigt, sigt)
  correa[v+1, ] <- pari$corre
  fracta[v+1, ] <- pari$fract
  varcomp[v+1, ] <- matrix(sigt,1,l^2)
  sevarcomp[v+1, ] <- as.vector(sesigt)
  labelc[v+1] <- "VarP(I)"

# SE's
  sep.list <- separ(varcomp, vsiga, v,l, fracta, correa)
# print(sep.list)

  outlist <- list(component=labelc,
          correlation=correa, correlation.variance=sep.list$vcorre,
          correlation.se=sqrt(sep.list$vcorre),
          fraction=fracta, fraction.variance=sep.list$vfract,
          fraction.se=sqrt(sep.list$vfract),
          variance.components=varcomp,variance.components.se=sevarcomp,
          phenotypic.variance=sigt, phenotypic.variance.se=sesigt,
          observed.variance=vara)
  return(outlist)
}
