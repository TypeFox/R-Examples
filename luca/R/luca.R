# luca - Likelihood Under Covariate Assumptions
# Copyright (C) 2005 J.Graham,  B.McNeney, Ji-Hyung Shin

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or   
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

########################################################################

"luca" <-
function(pen.model, gLabel, dat, HWP = FALSE, dep.model = NULL) {
require(survival)
require(genetics)

call<-match.call() #so we can report the call to users

n.obs<-nrow(dat)
dLabel<-all.vars(pen.model)[1] #disease status is on LHS of pen. model

if(!is.genotype(dat[,gLabel])) 
  stop(paste(gLabel,"must be a genotype object"))

g.levels<-levels(dat[,gLabel])
## Need g.levels to be a genotype object as well (for 
## operations like tem.g[1:length(tem.g)]<-g.levels[g] later on) and we
## need to keep order of the alleles within genotypes as-is
g.levels<-genotype(g.levels,reorder="no") 

datnames<-dimnames(dat)[[2]]
attrib.ind<-(datnames!=gLabel & datnames!=dLabel)
attribs<-dat[,attrib.ind, drop=FALSE] #Don't coerce to vector if only 1 column

## Set up "covariates", "response" and "stratum" for pseudo-subjects
## TODO: We just rbind pseudo-subjects onto pseudo.covar for 
## now, but Matt says this is inefficient and we should set
## up empty variables that can be filled in.
pseudo.covar<-NULL
pseudo.response<-NULL
pseudo.stratum<-NULL

disease.status<-c(0,1) #Assume 0/1 for now
for(l in 1:2) {
 for(g in 1:length(g.levels)) {
  tem.d<-rep(disease.status[l],n.obs) 
  tem.g<-dat[,gLabel] #Start w/ orig. so tem.g has all the levels/genotypes
  tem.g[1:length(tem.g)]<-g.levels[g] #overwrite w/ level/genotype of interest

  ## Start with disease status variable name hard-coded as d
  ## and genetic factor variable name hard-coded as g
  ## and then use names() to rename columns after tem.dat created
  tem.dat<-data.frame(d=tem.d,g=tem.g,attribs)
  names(tem.dat)[1]<-dLabel
  names(tem.dat)[2]<-gLabel

  pseudo.covar<-rbind( pseudo.covar,
     get.lambda(dat=tem.dat,pen.model=pen.model,gLabel=gLabel,HWP=HWP,
                dep.model=dep.model)  )
  pseudo.response<-c(pseudo.response, 
         as.numeric( (tem.d==dat[,dLabel]) & (tem.g==dat[,gLabel]) ) )
  pseudo.stratum<-c(pseudo.stratum, (1:n.obs) )
 }
}

## Extract the offset variable from the covariates matrix
## -- currently the first column of the matrix

offsets<-pseudo.covar[,1]
pseudo.covar<-pseudo.covar[,-1]

## Now call the cox proportional hazards function. 

## Survival times are all 1, failure/censored is the pseudo.response --
## all "censored" observations end up in the risk set this way
Y<-Surv(rep(1,length(pseudo.response)), pseudo.response) 

## The pseudo.stratum variable stratifies the analysis so that --
## 1. Each match set makes a separate contribution to the likelihood.
## 2. Within each match set the "effected" (real subjects) appears in 
## the numerator and everyone in the match set (all pseudo-subjects) appear
## in the denominator.

## DEBUG
#p.dat<-list(cov=pseudo.covar,resp=pseudo.response,strat=pseudo.stratum)
#return(p.dat)

fit<-coxph.fit(x=pseudo.covar,y=Y,strata=pseudo.stratum, offset=offsets, 
               control=coxph.control(), weights=rep(1,length(Y)),
               method="breslow", rownames=row.names(pseudo.covar)) 
## Notes: method doesn't matter if only one case per matchset (`, "breslow
##        is the default. 

luca.coef<-fit$coefficients[-1] #remove nuisance parameter "delta"
luca.var<-fit$var[-1,-1] #remove elements corresponding to "delta"
out<-list(call=call, coefficients=luca.coef, var=luca.var, iter=fit$iter)
class(out)<-"luca"
return(out)
}

########################################################################

"get.lambda" <-
function(dat,pen.model,gLabel,HWP=FALSE,dep.model=NULL) {
 if( (HWP==TRUE || !is.null(dep.model)) && !is.genotype(dat[,gLabel]) )
     stop(paste(gLabel,"must be a genotype object for models with \n G,A dependence, or when HWP==TRUE"))

 #Construct the "X" part of Lambda corresp. to the penetrance model
 X<-model.matrix(pen.model, dat)
 X<-X[,-1, drop=FALSE] #trim off interecept column -- don't coerce to
                       #vector if only one column left
 dimnames(X)[[2]]<-paste("penmod",dimnames(X)[[2]],sep=".")

 #Construct the "Z" part of Lambda corresp. to the genetic/dependence 
 #model in controls.

 if(!is.null(dep.model)) { #Then dep. btw G and A -- HWP doesn't matter
   Z<-get.Z.dep(dat,gLabel,dep.model)
   offsets<-rep(0,nrow(Z))
 } else if(HWP) { #No dep. model, and geno freqs in controls follow HWP
   Z<-get.Z.HWP(dat,gLabel)
   hvec<-as.numeric(heterozygote(dat[,gLabel])) #1 for heterozygote, 0 o.w.
   offsets<-log(1+hvec) #log 1 if hom, log 2 if het.
 } else { #no dependence, and use a saturated model for P(G|D=0)
   Z<-get.Z.saturated(dat,gLabel)
   offsets<-rep(0,nrow(Z))
 }
 
 #Finally, we need the disease status indicator
 dLabel<-all.vars(pen.model)[1] #disease status is on LHS of pen. model
 d.status<-dat[,dLabel]

 #Now can assemble Lambda
 Lambda<-cbind(offsets,d.status,Z,d.status*X)
 dimnames(Lambda)[[2]][1]<-"offset" #need better name
 dimnames(Lambda)[[2]][2]<-"delta" #need better name
 
 return(Lambda)
}


########################################################################

"get.Z.dep" <-
function(dat,gLabel,dep.model) {
 #attribute part of model
 Z1<-model.matrix(dep.model, dat)
 #genetic variable
 tem.model<-formula( paste("~",gLabel))
 Z2<-model.matrix(tem.model, dat)
 #Now trim off the intercept column
 Z2<-Z2[,-1, drop=FALSE]
 #Need separate intercept and regression parameters for each 
 #non-baseline level of the genetic variable. Make a matrix Z with
 #columns for nu_1,tau_1 , nu_2,tau_2 , etc. 
 Z<-matrix(nrow=nrow(Z2),ncol=(ncol(Z1)*ncol(Z2)))
 Znames<-rep("",(ncol(Z1)*ncol(Z2)))
 for(i in 1:ncol(Z2)) {
  Z[,((i-1)*ncol(Z1)+1):(i*ncol(Z1))]<-Z2[,i]*Z1
  Znames[((i-1)*ncol(Z1)+1):(i*ncol(Z1))]<-
    paste(dimnames(Z2)[[2]][i],dimnames(Z1)[[2]],sep=":")
 }
 dimnames(Z)[[2]]<-paste("covmod",Znames,sep=".")
 return(Z)
}


########################################################################

"get.Z.saturated" <-
function(dat,gLabel) {
 #genetic variable only
 tem.model<-formula( paste("~",gLabel))
 if(!is.factor(dat[,gLabel])) 
   dat[,gLabel]<-factor(dat[,gLabel])
 Z<-model.matrix(tem.model, dat)
 #Now trim off the intercept column
 Z<-Z[,-1, drop=FALSE] #Don't coerce to a vector if only one column left
 Znames<-dimnames(Z)[[2]]
 dimnames(Z)[[2]]<-paste("covmod",Znames,sep=".")
 return(Z)
}


########################################################################

"get.Z.HWP" <-
function(dat,gLabel) {
 if(!is.genotype(dat[,gLabel])) 
   stop(paste(gLabel,"must be a genotype object"))
 Z<-allele.count(dat[,gLabel])
 #Now trim off column of baseline allele counts -- will be most freq allele
 Z<-Z[,-1, drop=FALSE] #Don't coerce to a vector if only one column left
 Znames<-dimnames(Z)[[2]]
 dimnames(Z)[[2]]<-paste("covmod",Znames,sep=".")
 return(Z)
}

########################################################################

"summary.luca" <- 
function(object, ...) {

coefnames<-names(object$coef)
penmodel.ind<-substr(coefnames,start=1,stop=7)=="penmod."
coefnames<-substr(coefnames[penmodel.ind],start=8,stop=1000000)

coef.table<-cbind(object$coef[penmodel.ind],
                  sqrt(diag(object$var[penmodel.ind,penmodel.ind])))
numbeta<-length(object$beta)

#Compute z-scores for regression coefficients
coef.table<-cbind(coef.table,coef.table[,1]/coef.table[,2])

#two-sided p-values
coef.table<-cbind(coef.table,1-pchisq(coef.table[,3]^2,df=1))
dimnames(coef.table)<-list(coefnames,
                           c("Estimate","Std. Error","zscore","Pr(>|z|)"))

sumList <- list(call=object$call, coefficients=coef.table)
class(sumList) <- "summary.luca"
return(sumList)
}

########################################################################


print.summary.luca<-function(x, ...){
  cat("Call:\n")
  print(x[[1]])

  cat("\nCoefficients:\n")
  print(x[[2]])

}
