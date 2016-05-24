`vrais.LDL.add.pere` <-
function(moyenne.pere,alpha.Q,s,CD,perf,PLA,LD.m,LD.chrom1,LD.chrom2,mean.gene)
{
sd=sqrt(s)
stde=sd/CD
dnorm.Q    = dnorm(perf,mean.gene+moyenne.pere+alpha.Q,stde)
dnorm.q    = dnorm(perf,mean.gene+moyenne.pere-alpha.Q,stde)
dnorm.0    = dnorm(perf,mean.gene+moyenne.pere,stde)



PLA.dnormQ = PLA * dnorm.Q
PLA.dnormq = PLA * dnorm.q
PLA.dnorm0 = PLA * dnorm.0

LD.dnormQ  = LD.m * dnorm.Q
LD.dnormq  = LD.m * dnorm.q
LD.dnorm0  = LD.m * dnorm.0

LD.PLA.Q  = PLA * LD.dnormQ
LD.PLA.q  = PLA * LD.dnormq
LD.PLA.0  = PLA * LD.dnorm0

# voir formule simplifiée de la section 1.3.3. 
S1        = LD.dnormQ + dnorm.0 - LD.dnorm0
S2        = LD.dnorm0 + dnorm.q - LD.dnormq
S3        = LD.PLA.Q - 2*LD.PLA.0 + PLA.dnorm0 - PLA.dnormq + LD.PLA.q

#S1[abs(S1)<0.000000001]=0
#S2[abs(S2)<0.000000001]=0
#S3[abs(S3)<0.000000001]=0


total = S1 * LD.chrom2+ S2*(1-LD.chrom2)+S3 *(LD.chrom1-LD.chrom2)

total[total <=0]=-750
total[total>0]=log(total[total>0])

vrais.intra.pere = sum(total)


#write(c(alpha.Q,s,vrais.intra.pere),"res.txt",ncolumns=3,append=T,sep=" ; ")

vrais.intra.pere
}

