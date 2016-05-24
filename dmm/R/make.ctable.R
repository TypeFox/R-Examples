make.ctable <-
function(){
#  make.ctable() - setup lists of components
#                - VarGa() and VarGs() added
  cohortvar <- c("VarE(C)","VarE(M&C)","VarE(M&!C)")
  cohortcov <- c("CovE(I,M&!C)","CovE(M&!C,I)")
  cohort <- c(cohortvar,cohortcov)
  evar <- c("VarE(I)","VarE(M)")
  ecov <- c("CovE(I,M)","CovE(M,I)")
  e <- c(evar,ecov)

  sexlinaddgvar <- c("VarGs(Ia)","VarGs(Ma)")
  addgvar <- c("VarG(Ia)","VarG(Ma)",sexlinaddgvar)
  domgvar <- c("VarG(Id)","VarG(Md)")
  epiaddgvar <- c("VarG(Ia:a)","VarG(Ma:a)")
  epidomgvar <- c("VarG(Ia:d)","VarG(Ma:d)","VarG(Id:d)","VarG(Md:d)")

  gvar <- c(addgvar,domgvar,epiaddgvar,epidomgvar)
  allvar <- c(evar,gvar,cohortvar)

  sexlinaddgcov <- c("CovGs(Ia,Ma)","CovGs(Ma,Ia)")
  addgcov <- c("CovG(Ia,Ma)","CovG(Ma,Ia)",sexlinaddgcov)
  domgcov <- c("CovG(Id,Md)","CovG(Md,Id)")
  epiaddgcov <- c("CovG(Ia:a,Ma:a)","CovG(Ma:a,Ia:a)")
  epidomgcov <- c("CovG(Ia:d,Ma:d)","CovG(Ma:d,Ia:d)","CovG(Id:d,Md:d)","CovG(Md:d,Id:d)")

  gcov <- c(addgcov,domgcov,epiaddgcov,epidomgcov)
  g <- c(gvar,gcov)
  sexlinaddg <- c(sexlinaddgvar,sexlinaddgcov)
  allcov <- c(ecov,gcov,cohortcov)
  addg <- c(addgvar,addgcov)
  domg <- c(domgvar,domgcov)
  epiaddg <- c(epiaddgvar,epiaddgcov)
  epidomg <- c(epidomgvar,epidomgcov)

  indvar <- c("VarE(I)","VarG(Ia)","VarG(Id)","VarG(Ia:a)","VarG(Ia:d)","VarG(Id:d)","VarGs(Ia)")
  indcov <- NULL  
  # if indcov is NULL matcov is same as allcov
  ind <- c(indvar,indcov)

  matvar <- c("VarE(M)","VarG(Ma)","VarG(Md)","VarG(Ma:a)","VarG(Ma:d)","VarG(Md:d)","VarE(M&C)","VarE(M&!C)","VarGs(Ma)")
  matcov <- c("CovE(I,M)","CovE(M,I)","CovG(Ia,Ma)","CovG(Ma,Ia)","CovG(Id,Md)","CovG(Md,Id)","CovG(Ia:a,Ma:a)","CovG(Ma:a,Ia:a)","CovG(Ia:d,Ma:d)","CovG(Ma:d,Ia:d)","CovG(Id:d,Md:d)","CovG(Md:d,Id:d)","CovG(Ms,Is)",cohortcov,"CovGs(Ia,Ma)","CovGs(Ma,Ia)")
  mat <- c(matvar,matcov)

  all <- c(allvar,allcov)

  ctable <- list(cohortvar=cohortvar, cohortcov=cohortcov, cohort=cohort, evar = evar, ecov=ecov, e=e, addgvar=addgvar, domgvar=domgvar, epiaddgvar=epiaddgvar, epidomgvar=epidomgvar, sexlinaddgvar=sexlinaddgvar, gvar=gvar, allvar=allvar, addgcov=addgcov, domgcov=domgcov, epiaddgcov=epiaddgcov, epidomgcov=epidomgcov, sexlinaddgcov=sexlinaddgcov, gcov=gcov, g=g,  allcov=allcov, addg=addg, domg=domg, epiaddg=epiaddg, epidomg=epidomg, sexlinaddg=sexlinaddg, indvar=indvar, indcov=indcov, ind=ind,  matvar=matvar, matcov=matcov, mat=mat,  all=all)
  return(ctable)
}
