data(dipper)
# CJS with R
mod=crm(dipper,model="CJS",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)))
mod$results$beta
# CJS with ADMB
prepare_admb()
mod=crm(dipper,model="CJS",use.admb=TRUE,model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)))
mod$results$beta
# CJS with HMM in R
mod=crm(dipper,model="hmmCJS",use.admb=TRUE,model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)))
mod$results$par
# CJS with RMark
summary(mark(dipper,model="CJS",output=FALSE,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time))))$beta

# JS with R
mod=crm(dipper,model="JS",groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)))
mod$results$beta
# JS with RMark
summary(mark(dipper,model="POPAN",output=FALSE,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time))))$beta

#MSCJS with HMM
data(mstrata)
mod=crm(mstrata,model="hmmMSCJS",model.parameters=list(S=list(formula=~stratum),p=list(formula=~stratum),Psi=list(formula=~-1+stratum:tostratum)))
mod$results$beta

# MSCJS with RMark
summary(mark(mstrata,model="Multistrata",output=FALSE,model.parameters=list(S=list(formula=~stratum),p=list(formula=~stratum),Psi=list(formula=~-1+stratum:tostratum))))$beta

