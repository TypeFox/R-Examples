## Create data set modeled after the Isavuconazole study
# See Study 0104 and the ITT analysis in: 
# www.fda.gov/downloads/AdvisoryCommittees/
#    CommitteesMeetingMaterials/Drugs/
#    Anti-InfectiveDrugsAdvisoryCommittee/UCM430748.pdf
#
#  Figure 20 shows that data has about 20% failures by day 42 and 30% 
#  failures by day 84, so we use a exponential cure model with 60% cured
#  and exponential mean=60.5932 that matches those two failure rates. 
#  We censor only those that make it until day 100 without an event

R.Version()$version.string
set.seed(1)
n<- 258
y<- rexp(2*n,rate=1/60.5932)
cured<- rbinom(2*n,1,0.60)
y[cured==1]<- 100
status<-rep(1,2*n)
status[cured==1]<-0
z<-rep(c(1,2),each=n)

IsaSimData<-data.frame(y=y,z=z,status=status)

#nicq.test(y,g=nimDiffOR,delta0=0.10,q=0.20, z=z, status=status)


