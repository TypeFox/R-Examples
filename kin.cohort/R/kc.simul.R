kc.simul<-
function(nfam, f, hr, rand = 0, mean.sibs = 2, mean.desc = 1.5, 
         a.age = 8, b.age = 80, a.cancer = 3, b.cancer = 180){

# nfam = number of families
# f = allele frequency
# hr = hazard ratio
# rand = variance of random effect for cancer incidence (ratio of hr)

# ngenes = number of genes (not implemented)
# ndiseases = number of competing diseases (not implemented)

# model parameters
# mean.sibs<- 2       # mean number of sibllings and descendants (Poisson)
# mean.desc<- 1.5
# a.age<- 8           # shape for age
# b.age<- 80          # scale for age
# a.cancer<- 3        # shape for cancer incidence
# b.cancer<- 180      # scale for cancer incidence

p <- f^2 + 2*f*(1-f) # carrier probability

###################################################
# families
# generate a sample of probands and their family characteristics
nsib<-rpois(nfam,mean.sibs)
ndes<-rpois(nfam,mean.desc)
af<-rnorm(nfam, 0, log(hr)*rand) # family specific

###################################################
# Probands

age<-round(rweibull(nfam,a.age,b.age))   # range ~40-100, some asymetry
age[age<20]<-20                          # truncate extremes
age[age>100]<-100

gender<-rbinom(nfam, 1, 0.5)             # not used, but easy

carrier <-rbinom(nfam, 1, p)

agecancer<-round(rweibull(nfam,a.cancer, b.cancer/exp((carrier*log(hr)+af)/a.cancer) ))
cancer <-ifelse(agecancer<=age, 1, 0)    # event if age cancer <= current age
agecancer[cancer==0]<-age[cancer==0]     # censor at current age

###################################################
# Relatives

family<-function(famid,nsib,ndes,af,age,carrier,hr,p ){
# return matrix (relatives x vars) that can be appended
nrel<-2+nsib+ndes
rel<-c(1,1,rep(2,nsib),rep(3,ndes))
ages<-round(c( age+rnorm(2,10,5),
               age+rnorm(nsib, 0, 5),
               age+rnorm(ndes, -20, 5) ))

ages[ages<10]<-10
ages[ages>105]<-105


# depend on carrier, rel, p [cp1exact]
p.carrier<-sapply(rel,FUN=function(r, carrier,p){
   q<-1-p
   ifelse(carrier==0, 			#non carrier
      ifelse(r==2,
         (1/4)*p^2+(1/2)*p*(1+q), 
         p),          			 #carrier
      ifelse(r==2,
         ((1/4)*p^2*(1+p)^2 + p^2*q*(1+p)+p*q*(1+p*q))/(p^2+2*p*q),  
         p*(1+p*q)/(p^2+2*p*q) )
   )}, carrier, p )


rel.carrier<-rbinom(nrel,1,p.carrier)

agecancer<-round(rweibull(nrel,a.cancer,b.cancer/exp((rel.carrier*log(hr)+af)/a.cancer)))
cancer <-ifelse(agecancer<=ages, 1, 0)
agecancer[cancer==0]<-ages[cancer==0]

gender<-c(1,0, rbinom(nrel-2, 1, 0.5) )
rels<-cbind(rep(famid,nrel), rel, ages, gender, agecancer, cancer, rep(carrier,nrel),rel.carrier )
rels
}

###################################################
# families

fam<-lapply(1:nfam, function(id){
  proband<-rbind(c(id, 0, age[id], gender[id], agecancer[id], cancer[id], carrier[id], carrier[id] ))
  relatives<-family(id,nsib[id],ndes[id],af[id],age[id],carrier[id],hr,p )
  t(rbind(proband, relatives)) }
  )

fam<-unlist(fam)
fam<-t(matrix(fam,8,length(fam)/8))
colnames(fam)<- c("famid", "rel", "age", "gender", "agecancer", "cancer", "carrier", "real.carrier")
fam<-as.data.frame(fam)
class(fam)<-c("kin.cohort.sample","data.frame")
fam
}