# aveAge.R (plot of average ages in PY and observed cases)
rm(list=ls()) 
graphics.off()
.simpleCap <- function(s) paste(toupper(substring(s, 1, 1)), substring(s, 2),sep = "", collapse = " ")
library(SEERaBomb)
if (0) {  # 1 to run, 0 to comment (this way it can be folded while looking more useful than commented code)
  load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
  load("~/data/SEER/mrgd/popsae.RData") # loads in popsae (extended to ages 85-99)
  canc=canc%>%select(-reg,-COD,-radiatn,-histo3,-ICD9)
  popsa=popsae%>%group_by(db,race,sex,age,year)%>%summarize(py=sum(py)) # sum on regs
  # canc$cancer=as.character(canc$cancer) # need this if merging types into a new name, like lymphoma=NHL+HL
  canc$cancer[canc$cancer=="APL"] ="AML" # overwrite back to AML
  canc$cancer[canc$cancer=="AMLti"] ="AML" # overwrite back to AML
  canc$cancer=factor(canc$cancer) 
  pm=seerSet(canc,popsa,Sex="male",ageStart=0,ageEnd=100) #pooled (races) male seerSet
  head(pm$canc)
  pf=seerSet(canc,popsa,Sex="female",ageStart=0,ageEnd=100) #pooled (races) female seerSet
  head(pf$canc)
  pm=mk2D(pm,secondS=c("AML","MDS"))
  pf=mk2D(pf,secondS=c("AML","MDS"))
  levels(pm$canc$cancer)
  brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12) 
  pm=tsd(pm,brks=brks,trts=c("rad","noRad"))
  pf=tsd(pf,brks=brks,trts=c("rad","noRad"))
  system.time(save(pm,pf,file="~/Results/amlMDS/aveAgeFull.RData")) #~25 seconds 
}

load("~/Results/amlMDS/aveAgeFull.RData")
# load("~/Results/amlMDS/pm.RData")
# load("~/Results/amlMDS/pf.RData")
# debug(mkDF)
dm=mkDF(pm,1)
df=mkDF(pf,1)

names(pf$L[[1]][["rad"]])
head(dm)
library(ggplot2)
library(grid)
library(reshape2)
# glimpse(df)

D=rbind(dm%>%filter(cancer1=="prostate")%>%select(t,Observed=ageO,Expected=age,trt,cancer2)%>%mutate(cancer1="Prostate"),
    df%>%filter(cancer1=="breast")%>%select(t,Observed=ageO,Expected=age,trt,cancer2)%>%mutate(cancer1="Breast"))
D=melt(D,id.vars=c("t","trt","cancer1","cancer2"))
head(D)
D$trt=factor(as.numeric(D$trt),labels=c("Radiation","No Radiation"))


graphics.off()
quartz(width=6,height=6)
theme_set(theme_bw())
# theme_update()
theme_update(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.text=element_text(size=rel(1)),
             strip.text = element_text(size = rel(1.5)))

g=qplot(x=t,y=value,data=D,col=trt,shape=variable,geom=c("line","point"),
        xlab=paste("Years Since First Cancer Dx"),ylab="Age at Second Cancer Dx")
g=g+facet_grid(cancer1~cancer2,scales="free")+geom_abline(intercept=1, slope=0) 
g = g + scale_color_grey(start = 0, end = 0.6)
g+  theme(legend.title = element_blank(),legend.key.height=unit(.45, "cm"),legend.position = c(.22, .5))
ggsave(paste0("~/Results/amlMDS/prosBreastAge.png"))
ggsave(paste0("~/Results/amlMDS/prosBreastAge.eps"))

DL=rbind(dm%>%filter(cancer1=="prostate")%>%select(t,Observed=aoL,Expected=aeL,trt,cancer2)%>%mutate(cancer1="Prostate"),
        df%>%filter(cancer1=="breast")%>%select(t,Observed=aoL,Expected=aeL,trt,cancer2)%>%mutate(cancer1="Breast"))

DU=rbind(dm%>%filter(cancer1=="prostate")%>%select(t,Observed=aoU,Expected=aeU,trt,cancer2)%>%mutate(cancer1="Prostate"),
         df%>%filter(cancer1=="breast")%>%select(t,Observed=aoU,Expected=aeU,trt,cancer2)%>%mutate(cancer1="Breast"))
head(DU)
DU=melt(DU,id.vars=c("t","trt","cancer1","cancer2"))
DL=melt(DL,id.vars=c("t","trt","cancer1","cancer2"))
head(DU,3)
head(DL,3)
D$L=DL$value
D$U=DU$value
head(D,3)
D$t=D$t+0.1*(as.numeric(D$trt)-1)
quartz(width=8,height=6)
g=qplot(x=t,y=value,data=D,col=trt,shape=variable,geom=c("line","point"),
        xlab=paste("Years Since First Cancer Dx"),ylab="Age at Second Cancer Dx")
g=g+facet_grid(cancer1~cancer2,scales="free")+geom_abline(intercept=1, slope=0) 
g = g + scale_color_grey(start = 0, end = 0.6)
g=g+geom_errorbar(aes(ymin=L,ymax=U),width=.15)
g+  theme(legend.title = element_blank(),legend.key.height=unit(.45, "cm"),legend.position = c(.22, .5))
ggsave(paste0("~/Results/amlMDS/prosBreastAgeCI.png"))
ggsave(paste0("~/Results/amlMDS/prosBreastAgeCI.eps"))

summary(pm)
summary(pf)