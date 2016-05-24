# NOTE: This script assumes that the SEER data is in /data/SEER.  Please change
# this line to match your data location if it differs
.seerHome="~/data/SEER" 
rm(list=ls()) # note: dot variables defined above persist through such cleanings
# install.packages("RSQLite")
# install.packages("bbmle")
# install.packages("ggplot2")
# install.packages("plyr")
library(RSQLite)
m=dbDriver("SQLite")
library(plyr)
library(bbmle)
library(ggplot2)
graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=8,height=7)
# windows(height=7,width=8,xpos=-100,ypos=-100) #need to control the device size to
                                    # make things work across computers/screen sizes
# Note: Figure 1 was drawn using Adobe Illustrator so it is not reproduced here.
######################   
# Code 9863 versus 205.1 trends over the decades: Figure 2
con=dbConnect(m,dbname=file.path(.seerHome,"73/all.db"))
# con=dbConnect(m,dbname=file.path(.seerHome,"mrgd/cancDef.db"))
# dbListFields(con,"canc")
#######################################################################
dbListTables(con)
dbListFields(con,"pops")
dbListFields(con,"lymyleuk")
pop=dbGetQuery(con,"SELECT * from pops where popage>5 and popage<19")
head(pop)
(pop<-ddply(pop, .(popage,popsex,popyear), summarise,py=sum(population)))
pop$dec= cut(pop$popyear,breaks=c(1972,1984,1996,2013),labels=c("73to84","85to96","97to11"))
(pop<-ddply(pop, .(popage,popsex,dec), summarise,py=sum(py)))
head(pop,20)

d=dbGetQuery(con, 
#      "SELECT * from lymyleuk where histo2=9863 and seqnum<2 and agerec>5 and agerec<19")
    "SELECT * from lymyleuk where cancer='CML' and seqnum<2")
head(d)
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>5&d$agerec<19,]
d=d[!is.na(d$agerec),]

(d<-ddply(d, .(agerec,sex,yrdx), summarise,cases=length(agerec))) 
d$decade= cut(d$yrdx,breaks=c(1972,1984,1996,2013),labels=c("73to84","85to96","97to11"))
head(d)
(d<-ddply(d, .(agerec,sex,decade), summarise,cases=sum(cases)))
head(cbind(d,pop)) # just to see that they match up
d=cbind(d,py=pop[,"py"]) # only take the non-redundant py column
head(d)

d9=dbGetQuery(con, 
        "SELECT * from lymyleuk where ICD9=2051 and seqnum<2")
d9$agerec=as.numeric(cut(d9$agedx,breaks=c(0,1,seq(5,85,5),130)))
d9=d9[d9$agerec>5&d9$agerec<19,]
d9=d9[!is.na(d9$agerec),]
d9<-ddply(d9, .(sex,agerec,yrdx), summarise,cases=length(agerec))
d9$decade= cut(d9$yrdx,breaks=c(1972,1984,1996,2013),labels=c("73to84","85to96","97to11"))
d9<-ddply(d9, .(agerec,sex,decade), summarise,cases=sum(cases))
d9=cbind(d9,py=pop[,"py"])
head(d9)
(d=cbind(rbind(d,d9),code=gl(2,dim(d9)[1],labels=c("9863","205.1"))) ) 
head(d)
d=transform(d,incid=1e6*cases/py)
d$sex=factor(d$sex,labels=c("Male","Female"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)
names(d)[3]="Decade" # make legend start with capital
(p <- ggplot(d,aes(x=age,y=incid,shape=Decade,col=Decade))+geom_point(size=5)
 + labs(title="CML Incidence",x="Age (years)",
        y=expression(paste("Cases per ",10^6," Person-Years")))    
 + scale_y_log10(limits=c(3,200)) )
(p=p + facet_grid(code ~ sex))
mythem=theme( 
  plot.title = element_text(size = rel(1.5)),
  axis.title = element_text(size = rel(1.4)),
  axis.text = element_text(size = rel(1.2)),
  strip.text = element_text(size = rel(1.2))  )
p=p+mythem
(p=p+theme(legend.position = c(.60, .86),  #           legend.direction = 'vertical',
         legend.title = element_text(size = rel(1.4)) ,
  legend.text = element_text(size = rel(1.4))  ) )  

# ggsave(p,file="/users/radivot/igv/trends.wmf")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/trends.png")

#######################################################################
## CMML Incidence: Figure 3
#######################################################################
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
d=dbGetQuery(con, 
  "SELECT * from lymyleuk where histo3==9945 and seqnum<2")
  # "SELECT * from lymyleuk where histo3==9945 and seqnum<2 and agerec>7 and agerec<19")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>7&d$agerec<19,]
d=d[!is.na(d$agerec),]

(d<-ddply(d, .(sex,agerec), summarise,cases=length(agerec))) 
pops=dbGetQuery(con, "SELECT * from pops where popyear>2000 and popage>7 and popage<19")
(pop<-ddply(pops, .(popage,popsex), summarise,py=sum(population)))
d=cbind(d,py=pop[order(pop$popsex),"py"])
d=transform(d,incid=1e6*cases/py) # make it cases per million person years
d$Sex=factor(d$sex,labels=c("Male","Female"))
d$age=age[d$agerec+1]
head(d,15)

summary(mm1<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                  parameters=list(c~Sex),
                  method="Nelder-Mead",
                  start=list(c=-12,k=0.04),data=d)) 
d$EI=1e6*predict(mm1)/d$py
head(d,15)
(CI2=round(cbind(coef(mm1),confint(mm1)),4)  )
exp(-CI2["c.SexFemale",c(1,3,2)])
mkCIk=function(v) sprintf("k = %4.3f (%4.3f, %4.3f)",v[1],v[2],v[3])
(lab=mkCIk(CI2["k",]))

library(ggplot2)
(p <- ggplot(d,aes(x=age,y=incid,shape=Sex))+geom_point(size=5)
 + labs(title="CMML Incidence: SEER 2000-2011",x="Age (years)",
        y=expression(paste("Cases per ",10^6," Person-Years")))  #  y=expression(frac(Cases,paste(10^6," Person-Years")))) 
 + scale_y_log10(limits=c(0.1,70)) )

(p=p+theme(plot.title = element_text(size = rel(2.4)),
           axis.title = element_text(size = rel(2.4)),
           axis.text = element_text(size = rel(2)))  )
(p=p+theme(legend.position = c(0.8, .2), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(2))  ) )  
(p=p+geom_line(aes(y=EI))) # update just y component of aes

(p=p+annotate("text",x=35,y=70, hjust=0, label = paste("y ~ exp(C+k*age)  C~sex"),size=8 )   )  
(p=p+annotate("text",x=35,y=40,hjust=0,label=lab, size=8))

# ggsave(p,file="/users/Tomas/downloads/CMML.wmf")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/CMML.png")



###########################################################################################
# SEER CML versus Race Figure 4
########################################
pops=dbGetQuery(con, "SELECT * from pops where popyear>2000 and popage>5 and popage<19")
head(pops,20)
summary(as.factor(pops$poprace))
pops$poprace[pops$poprace>2]=3
(pop<-ddply(pops, .(popage,popsex,poprace), summarise,py=sum(population)))
head(pop)
d=dbGetQuery(con, 
# ICDO3 CML codes 9875 = bcr-abl+ and 9876 =bcr-abl neg CML are not in full use yet, i.e.
# many ICD-O3 CML codes are still 9863 carried over from ICD-O2
# "SELECT * from lymyleuk where histo3=9876 and seqnum<2 and race<98 and agerec>5 and agerec<19")
"SELECT * from lymyleuk where cancer='CML' and seqnum<2 and race<98")
# "SELECT * from lymyleuk where histo3=9863 and seqnum<2 and race<98 and agerec>5 and agerec<19")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>5&d$agerec<19,]
d=d[!is.na(d$agerec),]

d$race[d$race>2]=3
head(d)
(ddply(d, .(agerec,sex,race), summarise,cases=length(agerec))) #only 5 of the 16s 
count(d, vars=c("agerec","sex","race"))  # same here, and no options to fix it
(d=ddply(d, .(agerec,sex,race), summarise,cases=length(agerec),.drop=F))# this gets them all
head(cbind(d,pop)) # see alignment
d=cbind(d,py=pop[,"py"])
d=transform(d,incid=1e6*cases/py)
head(d,15)
d$Sex=factor(d$sex,labels=c("Male","Female"))
d$race=factor(d$race,labels=c("White","Black","Asian"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)
# ## playing around with stuff like this first
# (s1=summary(m1<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
#                      parameters=list(c~-1+sex:race,k~-1+race),   
#                      method="L-BFGS-B",
#                      lower=c(rep(-15,6),rep(0.01,3)),
#                      upper=c(rep(-10,6),rep(0.05,3)),
#                      start=list(c=-13,k=0.02),
#                      data=d,
#                      control=list(maxit=500,trace=1,parscale=c(rep(13,6),rep(0.03,3)) )
# ) )  )
# I finally found good initial parameter estimates
(s1=summary(m1<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                     parameters=list(c~-1+Sex:race,k~-1+race),   
                     start=list(c=-13,k=0.02),
                     data=d
) )  )


# This is fine if all we want is k estimates and absolute amplitude params
# But for M/F ratios we need deltaC difference parameters, obtainable as
(sA=summary(mA<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                     parameters=list(c~-1+race+Sex:race,k~-1+race), # = work to figure out  
                     start=list(c=-4,k=0.02),
                     data=d
) )  )

(cA=confint(mA)) #good  
d$EI=1e6*predict(mA)/d$py
head(d,20) #incidence predictions look reasonable
# d$EI=1e6*as.numeric(predict(mA))/d$py 
# here as.numeric clears .Dimnames in the matrix output of predict()  
# If ggplot gives warnings about rownames being the same, this may help


# And for Tf shifts we need 
(sT=summary(mT<-mle2(cases~dpois(lambda=exp(c+k*(age-Tf))*py),
                     parameters=list(c~-1+race,k~-1+race,Tf~-1+Sex:race),   
                    fixed=list("Tf.SexMale:raceWhite"=0,
                               "Tf.SexMale:raceBlack"=0,
                               "Tf.SexMale:raceAsian"=0),
                     start=list(c=-12,k=0.02,Tf=10),
               control=list(maxit=100,trace=0,parscale=c(rep(-12,3),rep(0.03,3),rep(10,3))),
                     data=d
) )  )
(cT=confint(mT))#but unfortunately, we cannot get Tf CI for Asians or Blacks
d$EIT=1e6*as.numeric(predict(mT))/d$py
head(d,20) #predictions of T model are close to the A model 

# # the next block attempts to reparameterize the T model to get Tf CI 
# # This can be collapsed for now, since I'm not sure if anything will come of it
# #######################
# # try moving it all into Tf
# (sT=summary(mT<-mle2(cases~dpois(lambda=exp(k*(age-Tf))*py),
#                      parameters=list(k~-1+race,Tf~-1+race+sex:race),   
#                     start=list(k=0.02,Tf=100),
#              control=list(maxit=100,trace=0,parscale=c(rep(.03,3),rep(300,3),rep(10,3))),
#                      data=d
# ) )  ) # bad fit, want around 470
# 
# # try controlling IC
# (IC=c("k.raceWhite"=0.04,"k.raceBlack"=0.03,"k.raceAsian"=0.02,
# "Tf.raceWhite"=13/0.04,"Tf.raceBlack"=13/0.04,"Tf.raceAsian"=13/0.02,
# "Tf.raceWhite:sexF"=10,"Tf.raceBlack:sexF"=10,"Tf.raceAsian:sexF"=20) )
# 
# (sT=summary(mT<-mle2(cases~dpois(lambda=exp(k*(age-Tf))*py),
#                      parameters=list(k~-1+race,Tf~-1+race+sex:race),   
#                      start=as.list(IC),
# #                      control=list(maxit=100,trace=0,parscale=IC),
#                      data=d
# ) )  )   
# 
# # => An error messsage that I cannot understand. Time is up, I surrender
# ##################

# Back to the basics: reduce the dimensionality by doing the races separately.

dw=subset(d,race=="White")

# first show no significant difference in k values
summary(mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
                       parameters=list(c~Sex,k~sex),
                       start=list(c=-12,k=0.04),data=dw))

(sw2=summary(mw2<-mle2(cases~dpois(lambda=exp(c+k*(age-Tf))*py),
                       method="Nelder-Mead",
                       parameters=list(Tf~Sex),
                       fixed=list("Tf.(Intercept)"=0),
                       start=list(c=-12,k=0.04,Tf=10),data=dw))     )

(cwT=confint(mw2))
cT    # good to see same answer for white Tf CI



db=subset(d,race=="Black")
summary(mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
             parameters=list(c~Sex,k~sex),
             start=list(c=-12,k=0.04),data=db))


(sb2=summary(mb2<-mle2(cases~dpois(lambda=exp(c+k*(age-Tf))*py),
                       method="Nelder-Mead",
                       parameters=list(Tf~Sex),
                       fixed=list("Tf.(Intercept)"=0),
                       start=list(c=-12,k=0.04,Tf=10),data=db))     )
(cbT=confint(mb2))
cT    # adding all those whites didn't tighten black k and c estimates.
# This is because in Poisson regression, the mean = the variance, so no error 
# estimate improvements result from pooling across races or sexes. Joint fits 
# here thus do more harm than good, as they increase the dimensionality of the
# optimization. The only advantage is that it compacts your codes. 
da=subset(d,race=="Asian")
summary(mle2(cases~dpois(lambda=exp(c+k*age)*py),method="Nelder-Mead",
             parameters=list(c~Sex,k~sex),
             start=list(c=-12,k=0.04),data=da))
# => removal of k dependence on sex throughout was OK

(sa2=summary(ma2<-mle2(cases~dpois(lambda=exp(c+k*(age-Tf))*py),
                       method="Nelder-Mead",
                       parameters=list(Tf~Sex),
                       fixed=list("Tf.(Intercept)"=0),
                       start=list(c=-12,k=0.04,Tf=10),data=da))     )
(caT=confint(ma2))
cT # same finding
# Conclusion: try going to lower dimensions when confint() fails on you

mkCIA=function(v) {v=exp(-v); sprintf("M/F = %4.2f (%4.2f, %4.2f)",v[1],v[3],v[2])}
mkCIT=function(v) sprintf("Tf = %3.1f (%3.1f, %3.1f)",v[1],v[2],v[3])
mkCIk=function(v) sprintf("k = %4.3f (%4.3f, %4.3f)",v[1],v[2],v[3])
(CIA=round(cbind(coef(mA),cA),3))
(AS=apply(CIA[4:6,],1,mkCIA))
(kS=apply(CIA[7:9,],1,mkCIk))
(CIT=cbind(c(coef(mw2)[4],coef(mb2)[4],coef(ma2)[4]),rbind(cwT,cbT,caT)[c(3,6,9),]) )
(TS=apply(CIT,1,mkCIT))

(cnts=ddply(d, .(race), summarise,cases=sum(cases)))
mkCases=function(n) sprintf("%d cases",n)
(cS=sapply(cnts[,2],mkCases))
(labS=c(kS,AS,TS,cS))
# note that capital S at the end of variable name => strings in my codes
(tdf=data.frame(race=rep(c("White","Black","Asian"),4),  #panel
           Sex=rep("Male",12),  # dummy, it wants to know sex of the "points"
           age=c(rep(30,9),rep(40,3)),      #x-coordinate of label
           incid=rep(c(2.5,1.67,1.1,90),each=3),     #y-coordinate of label
           Text=labS,stringsAsFactors = F) )           #label text


graphics.off()
# windows(height=6,width=12,xpos=-100,ypos=-100) #need to control the device size to
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=12,height=6)
library(ggplot2)
require(grid) # to avoid "cannot find unit()"
(p <- ggplot(d,aes(x=age,y=incid,col=Sex,shape=Sex)) + geom_point(size=7) 
 + labs(title="CML Incidence by Race, Sex and Age",x="Age",
        y=expression(paste("Cases per ",10^6," Person-Years")))    
 + scale_y_log10(limits=c(1,100)) +xlim(20,90) )
(p=p + facet_grid(. ~ race))
# (p=p+scale_shape_discrete(name = "") )

(p=p+theme(plot.title = element_text(size = rel(2)),
           strip.text = element_text(size = rel(2)),
           axis.title = element_text(size = rel(2)),
           axis.text = element_text(size = rel(2)))  )
(p=p+theme(legend.position = c(0.07, .75),    
            legend.direction = 'vertical',
#            legend.background = element_rect(size=4,color="red"),
#            legend.margin = unit(0, "cm"),
           legend.title = element_text(size = rel(1.5)) ,
           legend.text = element_text(size = rel(1.5))  ) )  
(p=p+geom_line(aes(y=EI))) 
(p=p + geom_text(aes(label=Text), size=6,data=tdf, hjust=0))

# ggsave(p,file="~/race.eps")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/race.png") 
#####################   END Figure 4

######################   
# 9863 k trends sampled every 5 years: Figure 5
con=dbConnect(m,dbname=file.path(.seerHome,"73/all.db"))
#######################################################################
dbListTables(con)
dbListFields(con,"pops")
dbListFields(con,"lymyleuk")
pop=dbGetQuery(con,"SELECT * from pops where popage>5 and popage<19 and poprace=1")
head(pop)
(pop<-ddply(pop, .(popage,popsex,popyear), summarise,py=sum(population)))
yearS=c("73to77","78to82","83to87","88to92","93to99","00to04","05to13")
pop$Years=cut(pop$popyear,breaks=c(1972,1977,1982,1987,1992,1999,2004,2014),labels=yearS)
(pop<-ddply(pop,.(popage,popsex,Years),summarise,py=sum(py)))
head(pop,20)

d=dbGetQuery(con, 
"SELECT * from lymyleuk where race=1 and cancer='CML' and seqnum<2")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>5&d$agerec<19,]
d=d[!is.na(d$agerec),]

# "SELECT * from lymyleuk where race=1 and histo3=9863 and seqnum<2 and agerec>5 and agerec<19")
# "SELECT * from lymyleuk where race=1 and ICD9=2051 and seqnum<2 and agerec>5 and agerec<19")
head(d)
d<-ddply(d, .(agerec,sex,yrdx), summarise,cases=length(agerec))
d$Years=cut(d$yrdx,breaks=c(1972,1977,1982,1987,1992,1999,2004,2014),labels=yearS)
head(d)
d<-ddply(d, .(agerec,sex,Years), summarise,cases=sum(cases))
head(cbind(d,pop)) # just to see that they match up
d=cbind(d,py=pop[,"py"]) # only take the non-redundant py column
head(d)
d$sex=factor(d$sex,labels=c("Male","Female"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)

# summary(mk<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
#                   parameters=list(c~sex,k~-1+Years),
#                   method="Nelder-Mead",
#                   start=list(c=-12,k=rep(0.05,7)),data=d)) 

years=c(1975,1980,1985,1990,1996,2002,2007)
kic=0.05
k=NULL
kL=NULL
kU=NULL

for (Y in yearS) {
  summary(mk<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                   parameters=list(c~sex),
                   method="Nelder-Mead",
                   start=list(c=-12,k=kic),data=subset(d,Years==Y))) 
  k=c(k,coef(mk)[3])
  (CIk=round(confint(mk)[3,],4))
  kL=c(kL,CIk[1])
  kU=c(kU,CIk[2])
}
dfk9=data.frame(years,k,kL,kU,Registries="SEER9")

con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
pop=dbGetQuery(con, 
"SELECT * from pops where poprace=1 and popyear>1999 and popage>5 and popage<19")
pop<-ddply(pop, .(popage,popsex,popyear), summarise,py=sum(population))
head(pop,20)
yearS=c("00to01","02to03","04to05","06to07","08to13")
pop$Years=cut(pop$popyear,breaks=c(1999,2001,2003,2005,2007,2014),labels=yearS)
pop<-ddply(pop,.(popage,popsex,Years),summarise,py=sum(py))
head(pop,20)

d=dbGetQuery(con, 
"SELECT * from lymyleuk where race=1 and cancer='CML' and seqnum<2")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>5&d$agerec<19,]
d=d[!is.na(d$agerec),]

# "SELECT * from lymyleuk where race=1 and histo3=9863 and seqnum<2 and agerec>5 and agerec<19")
# "SELECT * from lymyleuk where race=1 and ICD9=2051 and seqnum<2 and agerec>5 and agerec<19")
d<-ddply(d, .(agerec,sex,yrdx), summarise,cases=length(agerec))
d$Years=cut(d$yrdx,breaks=c(1999,2001,2003,2005,2007,2014),labels=yearS)
d<-ddply(d, .(agerec,sex,Years), summarise,cases=sum(cases))
head(cbind(d,pop)) # just to see that they match up
d=cbind(d,py=pop[,"py"])
head(d,15)
d$Sex=factor(d$sex,labels=c("Male","Female"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)

years=c(2001,2003,2005,2007,2009)
kic=0.05
k=NULL
kL=NULL
kU=NULL

for (Y in yearS) {
  summary(mk<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                   parameters=list(c~sex),
                   method="Nelder-Mead",
                   start=list(c=-12,k=kic),data=subset(d,Years==Y))) 
  k=c(k,coef(mk)[3])
  (CIk=round(confint(mk)[3,],4))
  kL=c(kL,CIk[1])
  kU=c(kU,CIk[2])
}

dfk18=data.frame(years,k,kL,kU,Registries="SEER18")
dfk=rbind(dfk9,dfk18)

if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )
windows(width=9,height=7)
library(ggplot2)
pd <- position_dodge(1) 
(p=ggplot(dfk, aes(x=years, y=k,shape=Registries)) + 
   geom_point(size=6,position=pd))
(p=p+labs(title="SEER9 1973-2011 and SEER18 2000-2013",x="Year", y="k") )     
(p=p+geom_errorbar(aes(ymin=kL, ymax=kU),width=.01,position=pd))
(p=p+theme(plot.title = element_text(size = rel(2)),
           axis.title = element_text(size = rel(2)),
           axis.text = element_text(size = rel(2)))  )
(p=p+theme(legend.position = c(0.2, .2), 
           legend.title = element_text(size = rel(1.5)) ,
           legend.text = element_text(size = rel(1.5))  ) )  
# ggsave(p,file="/users/radivot/igv/ktrends.wmf")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/ktrends.png")

###################### End Fig 5 on k trends



#######################################################################
## AML Incidence: Figure 11
#######################################################################
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
d=dbGetQuery(con, 
       "SELECT * from lymyleuk where icd9==2050 and seqnum<2")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>3&d$agerec<19,]
d=d[!is.na(d$agerec),]

(d<-ddply(d, .(sex,agerec), summarise,cases=length(agerec))) 
pops=dbGetQuery(con, "SELECT * from pops where popyear>2000 and popage>3 and popage<19")
(pop<-ddply(pops, .(popage,popsex), summarise,py=sum(population)))
d=cbind(d,py=pop[order(pop$popsex),"py"])
d=transform(d,incid=1e5*cases/py) # make it cases per million person years
d$Sex=factor(d$sex,labels=c("Male","Female"))
d$age=age[d$agerec+1]
head(d,15)

library(ggplot2)
(p <- ggplot(d,aes(x=age,y=incid,shape=Sex))+geom_point(size=5)
 + labs(title="AML Incidence: SEER 2000-2011",x="Age (years)",
        y=expression(paste("Cases per ",10^5," Person-Years")))  #  y=expression(frac(Cases,paste(10^6," Person-Years")))) 
 + scale_y_log10(limits=c(0.7,30)) )

(p=p+theme(plot.title = element_text(size = rel(2.4)),
           axis.title = element_text(size = rel(2.4)),
           axis.text = element_text(size = rel(2)))  )
(p=p+theme(legend.position = c(0.8, .2), 
           legend.title = element_text(size = rel(2)) ,
           legend.text = element_text(size = rel(2))  ) )  

# ggsave(p,file="/users/radivot/igv/AML.wmf")
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/AML.png")



####################################### End paper
# The following figure was taken out because it was redundant with the race figure
#############################
# CML Incidence in SEER 2000-2009 (parallel lines) 
##############################
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
d=dbGetQuery(con, 
"SELECT * from lymyleuk where cancer='CML' and seqnum<2")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>5&d$agerec<19,]
d=d[!is.na(d$agerec),]

# "SELECT * from lymyleuk where histo3=9863 and seqnum<2 and agerec>5 and agerec<19")
(d<-ddply(d, .(agerec,sex), summarise,cases=length(agerec))) 
pops=dbGetQuery(con, "SELECT * from pops where popyear>2000 and popage>5 and popage<19")
(pop<-ddply(pops, .(popage,popsex), summarise,py=sum(population)))
d=cbind(d,py=pop[,"py"])
d=transform(d,incid=1e6*cases/py)
head(d,15)
d$sex=factor(d$sex,labels=c("M","F"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)

summary(m0<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                 parameters=list(c~sex,k~sex),   # see that slopes do not depend on sex
                 method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=d))

summary(m1<-mle2(cases~dpois(lambda=exp(c+k*age)*py),
                 parameters=list(c~sex),
                 method="Nelder-Mead",
                 start=list(c=-12,k=0.04),data=d)) 
anova(m0,m1)  # => reject null that extra param carries more than noise
d$EI=1e6*predict(m1)/d$py
head(d,15)
(CIA=round(cbind(coef(m1),confint(m1)),4)  )

# reparameterizing the model to get the delay time Tf for females
(s2=summary(m2<-mle2(cases~dpois(lambda=exp(c+k*(age-Tf))*py),
                     method="Nelder-Mead",
                     parameters=list(Tf~sex),
                     fixed=list("Tf.(Intercept)"=0),
                     start=list(c=-12,k=0.04,Tf=10),data=d))     )
coef(s2)
coef(m2)


AICtab(m0,m1,m2)  # should be no difference with m1, it's just a reparameterization
(CIT=round(cbind(coef(s2)[,1],confint(m2)),4)  )

mkCIA=function(v) {v=exp(-v); sprintf("M/F = %4.2f (%4.2f, %4.2f)",v[1],v[3],v[2])}
mkCIT=function(v) sprintf("Tf = %3.1f (%3.1f, %3.1f)",v[1],v[2],v[3])

(p <- ggplot(d,aes(x=age,y=incid,shape=sex))+geom_point(size=5)
 + labs(title="CML Incidence: SEER 2000-2013",x="Age (years)",
        y=expression(frac(Cases,paste(10^6," Person-Years")))) 
 + scale_y_log10(limits=c(3,100)) )
(p=p+theme(title = element_text(size = rel(2)),axis.text = element_text(size = rel(2)))  )
(p=p+geom_line(aes(y=EI))) # update just y component of aes
(p=p+annotate("text",x=30,y=100, hjust=0, label = paste("y ~ exp(c+k*age)  c~sex"),size=7 )   )  
(p=p+annotate("text",x=30,y=80,hjust=0,label=mkCIk(CIA["k",]), size=6))
(p=p+annotate("text",x=30,y=65,hjust=0,label=mkCIA(CIA["c.sexF",]), size=6 )   )  
(p=p+annotate("text",x=30,y=47, hjust=0, label = paste("y ~ exp(c+k*(age-Tf))"),size=7)   )  
(p=p+annotate("text",x=30,y=37,hjust=0, label=mkCIT(CIT["Tf.sexF",]), size=6))
# ggsave(p,file="/users/radivot/case/grants/sachs/HSCepi/figs/CML.wmf")

###########################################################
# This section explores CML ICD-O3 codes in whites as 9863, 9875 (ba+) and 9876 (ba-)
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
pop=dbGetQuery(con, 
     "SELECT * from pops where poprace=1 and popyear>1999 and popage>9 and popage<19")
pop<-ddply(pop, .(popage,popsex,popyear), summarise,py=sum(population),.drop=F)
head(pop,20)
yearS=c("01to03","04to05","06to07","08to13")
pop$Years=cut(pop$popyear,breaks=c(1999,2003,2005,2007,2014),labels=yearS)
pop<-ddply(pop,.(popage,popsex,Years),summarise,py=sum(py),.drop=F)
head(pop,20)
dim(pop)

d=dbGetQuery(con,"SELECT * from lymyleuk where race=1 and ICD9=2051 and seqnum<2")
d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
d=d[d$agerec>9&d$agerec<19,]
d=d[!is.na(d$agerec),]

# "SELECT * from lymyleuk where race=1 and histo2=9863 and
#               seqnum<2 and agerec>9 and agerec<19")
# "SELECT * from lymyleuk where race=1 and 
#              ((histo3=9863) or (histo3=9876) or (histo3=9875)) and
#              seqnum<2 and agerec>9 and agerec<19")
# "SELECT * from lymyleuk where race=1 and 
#              (histo3=9863 or histo3=9875) and
#              seqnum<2 and agerec>9 and agerec<19")
# "SELECT * from lymyleuk where race=1 and histo3=9875 and
# seqnum<2 and agerec>9 and agerec<19")
# "SELECT * from lymyleuk where race=1 and histo3=9945 and
# seqnum<2 and agerec>9 and agerec<19")
head(d)
d<-ddply(d, .(agerec,sex,yrdx), summarise,cases=length(agerec),.drop=F)
sum(d$cases)  #2970 icdo2 9863 cases is same as icdO3 9863+9875+9876
# of the latter, 2425 are 9863, 507 are 9875 and only 38 are 9876  
d$Years=cut(d$yrdx,breaks=c(1999,2003,2005,2007,2014),labels=yearS)
d<-ddply(d, .(agerec,sex,Years), summarise,cases=sum(cases),.drop=F)
head(cbind(d,pop)) # just to see that they match up
d=cbind(d,py=pop[,"py"])
head(d,15)
d$Sex=factor(d$sex,labels=c("Male","Female"))
age=c(0.5,3,seq(7.5,87.5,5))
d$age=age[d$agerec+1]
head(d,15)
years=c(2002,2005,2007,2009)
kic=0.05
c=k=NULL
cL=kL=NULL
cU=kU=NULL
for (Y in yearS) {
  summary(mk<-mle2(cases~dpois(lambda=exp(c+k*(age-65))*py),
                   parameters=list(c~sex),
                   method="Nelder-Mead",
                   start=list(c=-12,k=kic),data=subset(d,Years==Y))) 
  k=c(k,coef(mk)[3])
  c=c(c,coef(mk)[1])
  CI=confint(mk)
  (CIk=round(CI[3,],4))
  (CIc=round(CI[1,],4))
  kL=c(kL,CIk[1])
  kU=c(kU,CIk[2])
  cL=c(cL,CIc[1])
  cU=c(cU,CIc[2])
}
data.frame(years,c,cL,cU)
data.frame(years,k,kL,kU)

# 9876 = ba neg
# > data.frame(years,k,kL,kU)
# years          k      kL     kU
# 1  2002 0.04739624  0.0042 0.0918
# 2  2005 0.11644800  0.0455 0.2098
# 3  2007 0.06525722  0.0131 0.1207
# 4  2009 0.04218057 -0.0077 0.0926

# 9875 = ba + 
# years           k      kL     kU
# 1  2002 0.012260690 -0.0035 0.0276
# 2  2005 0.009064363 -0.0086 0.0234
# 3  2007 0.007569936 -0.0078 0.0245
# 4  2009 0.004961402 -0.0077 0.0165

# so bcr-abl negative CML's do have a larger k/aging rate constant, 
# but since we used icdo2 9863, this cannot explain the ongoing drops in k

