# setup.R
library(SEERaBomb)
# first chunk builds the SEER binaries 
(df=getFields("/Users/radivot/data/SEER")) 
(rdf=pickFields(df)) # this uses default picks of fields
mkSEER(rdf,seerHome="/Users/radivot/data/SEER") 
# the rest fixes tAML and rewrites the fix into cancDef.RData
load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
yearEnd=2013
sc=canc%>%filter(histo3==9920|histo3==9987) #tAML and tMDS  subset of cancs (sc)
table(sc$yrdx,sc$histo3)  # see that they start in 2001
table(sc$ICD9,sc$histo3)  # only 2 9920 are ICD9 9999, so cant reverse it using ICD9
head(sc)
sc%>%group_by(histo3)%>%summarize(age=mean(agedx))# for kicks, see that tMDS cases arrive 5 years later in life
d=sc%>%filter(yrdx>2000)%>%group_by(sex,yrdx,histo3)%>%summarize(cnt=n()) #count rows in sex-year-ICDO3 groups
# Zeros for MDS in 2012 do not show up in d. The next line adds them in just to make the plots look right.
d=rbind(d,data.frame(sex=c("male","female"),yrdx=2012,histo3=9987,cnt=0)) 
library(ggplot2)
theme_update(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             strip.text = element_text(size = rel(1.5)))
qplot(yrdx,cnt,data=d,facets=sex~histo3,ylab="Cases",xlab="Year of Diagnosis")+
  scale_x_continuous(breaks=seq(2002,2010,4))
ggsave("~/Results/amlMDS/tAMLfix.png")  

#fit line to left up to 2009, predict and give right the difference
m=lm(cnt~yrdx,subset(d,histo3==9920&yrdx<2010&sex=="male"))
f=lm(cnt~yrdx,subset(d,histo3==9920&yrdx<2010&sex=="female"))
nd=list(yrdx=2010:yearEnd)
dm=d%>%filter(yrdx>=2010&histo3==9920,sex=="male")
df=d%>%filter(yrdx>=2010&histo3==9920,sex=="female")
(xm=dm$cnt-round(predict(m,nd))) # numbers of male tAMLs that need to become tMDSs in 2010, 2011, and 2012
(xf=df$cnt-round(predict(f,nd))) # same for females
canc$id=1:dim(canc)[1]  # set up an id column to use to pull the random sample
set.seed(8675309)   # Call Jenny to make it reproducible. 
rws=NULL
for (i in 1:3){
  rws=c(rws,sample((canc%>%filter(histo3==9920&sex=="male"&yrdx==(2009+i)))$id,xm[i]))
  rws=c(rws,sample((canc%>%filter(histo3==9920&sex=="female"&yrdx==(2009+i)))$id,xf[i]))
} 
canc[rws,"cancer"]="MDS" # note that we don't need tMDS and tAML as separate cancer types
save(canc,file="~/data/SEER/mrgd/cancDef.RData") # this takes a bit of time because it's writing a big file

