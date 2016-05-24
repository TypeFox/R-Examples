# specificFirst.R (makes plots in Figure 4)
# first run amlMDScomputing.R, then start here to make time course plots 
rm(list=ls())
.simpleCap <- function(s) paste(toupper(substring(s, 1, 1)), substring(s, 2),sep = "", collapse = " ")
library(dplyr)
library(reshape2)
library(ggplot2)
library(grid)
graphics.off()
quartz(width=7,height=4)
theme_set(theme_bw())
theme_update(legend.position = c(.40, .75))
theme_update(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(0.9)),
             legend.text=element_text(size=rel(0.9)),
             strip.text = element_text(size = rel(1.5)))

system.time(load("~/Results/amlMDS/pm.RData")) 
system.time(load("~/Results/amlMDS/pf.RData")) 
names(pm$L)
hires=TRUE
hires=FALSE
if (hires) {
  FirstS=c("prostate","breast") # first cancers to be plotted 
  brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12) 
} else {
  FirstS=c("MM","HL","NHL","lung") # first cancers to be plotted 
  brks=c(0,1,2,3,6,9,12)  # intermediate level resolution
}

dm=mkDF(pm,brks)
df=mkDF(pf,brks)
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(cancer1%in%FirstS)%>%select(Sex,cancer1,cancer2,trt,RR,rrL,rrU,t,int)
d=d%>%filter(!((cancer1=="breast")&(Sex=="Male"))) # remove male breast cancer cases
d$Radiation=factor(as.numeric(d$trt),labels=c("Yes","No"))
d$trt=NULL
d$int=NULL
head(d)

for (First in FirstS) {
  D=d%>%filter(cancer1==First)
  if (First=="prostate")   theme_update(legend.position = c(.87, .87)) else 
    if (First=="breast")  theme_update(legend.position = c(.87, .87)) else 
      if (First%in%c("lung","MM"))    theme_update(legend.position = c(.58, .87)) else  
        theme_update(legend.position = c(.4, .87))
  g=qplot(x=t,y=RR,data=D,col=Radiation,geom=c("line","point"),#xlim=c(-.1,24),
          xlab=paste("Years Since",.simpleCap(First),"First Cancer"),ylab="Relative Risk")
  if (First=="breast"|First=="prostate") g=g+facet_grid(cancer2~.,scales="free")+geom_abline(intercept=1, slope=0) else
    g=g+facet_grid(Sex~cancer2,scales="free")+geom_abline(intercept=1, slope=0) #g=g+scale_y_continuous(breaks=1:10)
  g = g + scale_color_grey(start = 0, end = 0.6)
    g=g+geom_errorbar(aes(ymin=rrL,ymax=rrU,width=.15))+ theme(legend.key.height=unit(.45, "cm"))
  print(g)
  ggsave(paste0("~/Results/amlMDS/",First,ifelse(hires,"hires",""),".eps"))
  # ggsave(paste0("~/Results/amlMDS/",First,ifelse(hires,"hires",""),".png"))
} # for loop on firsts


