# PYnMDSvsDBs.R   (Figure 1)
# the following was made earlier using SEERaBomb's mkSEER
load("~/data/SEER/mrgd/popsa.RData")
library(dplyr)
library(ggplot2)
(pys=popsa%>%group_by(db,year)%>%summarize(npy=sum(py)))
(pys2=pys%>%mutate(year=year+0.999)) # to make step function plots hold flat across each year
pys=rbind(pys,pys2)
graphics.off() # start fresh on plots
quartz(width=7,height=5)
# quartz(width=7,height=5,dpi=300)
# quartz(width=3.4,height=2.5,dpi=300)
# theme_set(theme_bw())
theme_update(axis.title = element_text(size = rel(2.3)),
             axis.text = element_text(size = rel(1.5)),
             legend.position = c(.1, .16),
             legend.title=element_text(size=rel(1.2)),
             legend.text=element_text(size=rel(1.2)),
             strip.text = element_text(size = rel(1.5)))
g=ggplot(pys)+geom_ribbon(aes(x=year,ymax=npy/1e6,ymin=0,fill=db))+facet_grid(db~.) #+ theme_bw()
# g = g + scale_fill_grey(start = 0.8, end = 0)
g=g+labs(y=expression(paste(10^6," Person Years")),x="Year") 
g+guides(fill = guide_legend("SEER\nDatabase"))+scale_x_continuous(breaks=c(1973,seq(1980,2010,by=10))) 

ggsave("~/Results/amlMDS/fig1A.png")
ggsave("~/Results/amlMDS/fig1A.eps")

load("~/data/SEER/mrgd/cancDef.RData") #loads in canc
mds=canc%>%filter(cancer=="MDS")
mds=mds%>%filter(yrdx>2000) # filter out 22 scattered cases before 2001 (that may be typos, i.e. bad data)
head(mds)
(D=mds%>%group_by(db,yrdx)%>%summarize(cases=n()))
(D2=D%>%mutate(yrdx=yrdx+0.999))
D=rbind(D,D2)

# quartz(width=7,height=5,dpi=300)
quartz(width=7,height=5)
# quartz(width=3.4,height=2.5)
theme_update(axis.title = element_text(size = rel(2.3)),
             axis.text = element_text(size = rel(1.5)),
             legend.position = "none",
             strip.text = element_text(size = rel(1.5)))
g=ggplot(D)+geom_ribbon(aes(x=yrdx,ymax=cases,ymin=0,fill=db))+facet_grid(db~.)
# g = g + scale_fill_grey(start = 0.8, end = 0)
g=g+labs(y="MDS Cases",x="Year")+scale_x_continuous(breaks=seq(2001,2013,by=2)) 
g+scale_y_continuous(breaks=seq(0,2000,by=1000)) 
ggsave("~/Results/amlMDS/fig1B.png")
ggsave("~/Results/amlMDS/fig1B.eps")



