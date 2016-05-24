rtrios<-function(n,ssize,f0,f1,f2,g=1,tmother,tfather,maf){

#-----------------------------------------
# Simulation of parents genotype
#-----------------------------------------
gen=rbern(n*4,maf)
gf1=gen[1:n]				    #father allele 1
gf2=gen[(n+1):(2*n)]		#father allele 2
gm1=gen[(2*n+1):(3*n)]	#mother allele 1
gm2=gen[(3*n+1):(4*n)]	#mother allele 2
gf=gf1+gf2
gm=gm1+gm2
gp=rbind(gf1,gf2,gf,gm1,gm2,gm)

#--********-----
# clean disk
#--********-----
rm(gen,gf1,gf2,gf,gm1,gm2,gm)

#-----------------------------------------
# Simulation of transmission index
#-----------------------------------------
tf=rbern(n,tfather)	#father
tm=rbern(n,tmother)	#mother

#-----------------------------------------
# Simulation of child genotype received
#-----------------------------------------
gc2=tf
gc2[gp[1,]==gp[2,]]<-gp[1,][gp[1,]==gp[2,]]
gc1=tm
gc1[gp[4,]==gp[5,]]<-gp[4,][gp[4,]==gp[5,]]

gc=gc1+gc2
gt=rbind(gp,gc1,gc2,gc)

# m = mother transmit minor allele for heterozygous child
# m2 = father transmit minor allele for heterozygous child

m=as.numeric((gc1==1)&(gc==1))
m2=as.numeric((gc2==1)&(gc==1))
#--********-----
# clean disk
#--********-----
rm(gc1,gc2,gc,gp,tf,tm)

#-----------------------------------------
# Simulation of child phenotype based on penetrance factors
# of child and maternal factors
#-----------------------------------------
pheno=rep(0,n)

p222=rbern(n,f2)
p212=rbern(n,f2)
p211=rbern(n,f1*g)
p122=rbern(n,f2)
p121=rbern(n,f1)
p201=rbern(n,f1*g)
p021=rbern(n,f1)
p112=rbern(n,f2)
p111m=rbern(n,f1*g)
p111f=rbern(n,f1)
p110=rbern(n,f0)
p101=rbern(n,f1*g)
p100=rbern(n,f0)
p011=rbern(n,f1)
p010=rbern(n,f0)
p000=rbern(n,f0)

pheno[gt[6,]==2&gt[3,]==2&gt[9,]==2&p222==1]=1
pheno[gt[6,]==2&gt[3,]==1&gt[9,]==2&p212==1]=1
pheno[gt[6,]==2&gt[3,]==1&gt[9,]==1&p211==1]=1
pheno[gt[6,]==1&gt[3,]==2&gt[9,]==2&p122==1]=1
pheno[gt[6,]==1&gt[3,]==2&gt[9,]==1&p121==1]=1
pheno[gt[6,]==2&gt[3,]==0&gt[9,]==1&p201==1]=1
pheno[gt[6,]==0&gt[3,]==2&gt[9,]==1&p021==1]=1
pheno[gt[6,]==1&gt[3,]==1&gt[9,]==2&p112==1]=1
pheno[gt[6,]==1&gt[3,]==1&gt[9,]==1&p111m==1&m==1]=1
pheno[gt[6,]==1&gt[3,]==1&gt[9,]==1&p111f==1&m==0]=1
pheno[gt[6,]==1&gt[3,]==1&gt[9,]==0&p110==1]=1
pheno[gt[6,]==1&gt[3,]==0&gt[9,]==1&p101==1]=1
pheno[gt[6,]==1&gt[3,]==0&gt[9,]==0&p100==1]=1
pheno[gt[6,]==0&gt[3,]==1&gt[9,]==1&p011==1]=1
pheno[gt[6,]==0&gt[3,]==1&gt[9,]==0&p010==1]=1
pheno[gt[6,]==0&gt[3,]==0&gt[9,]==0&p000==1]=1

id=c(1:n)
data=t(rbind(gt,pheno))
data=cbind(id,data[,c(10,6,3,9,7,8)])

#--********-----
# clean disk
#--********-----
rm(id,gt,p222,p212,p211,p122,p121,p201,p021)
rm(p112,p111m,p111f,p110,p101,p100,p011,p010,p000,pheno)

#-----------------------------------------
# calculate disease prevalence
#-----------------------------------------
case=data[data[,2]==1,]
ctrl=data[data[,2]==0,]
d=dim(case)[1]/n

#-----------------------------------------
# collect case-trios and control-trios
#-----------------------------------------
case=data[data[,2]==1,][1:ssize,2:5]
ctrl=data[data[,2]==0,][1:ssize,2:5]
idx.case=grep(TRUE, data[,2]==1)[1:ssize]
idx.ctrl=grep(TRUE, data[,2]==0)[1:ssize]

case=data.frame(cbind(case,m[idx.case]))
ctrl=data.frame(cbind(ctrl,m[idx.ctrl]))
names(case)=c('pheno','m','f','c','idx')
names(ctrl)=c('pheno','m','f','c','idx')

out=list(case,ctrl)
names(out)=c('case','ctrl')
out

}
