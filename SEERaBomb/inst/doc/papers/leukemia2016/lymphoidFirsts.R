# lymphoidFirsts.R (Figure 4)
# first run amlMDScomputing.R, then start here to make time course plots 
rm(list=ls())
library(dplyr)
library(reshape2)
library(ggplot2)
library(SEERaBomb)
LC=c("CLL","SLL","HCL","NHL","MM","hodgkin") # lymphoid first cancers
system.time(load("~/Results/amlMDS/pm.RData")) # 4 secs to load. 
system.time(load("~/Results/amlMDS/pf.RData")) # 4 secs to load. 
brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12)
df=mkDF(pf,brks)
dm=mkDF(pm,brks)
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(cancer1%in%LC)%>%group_by(cancer2,Sex,int)%>%summarize(O=sum(O),E=sum(E),t=weighted.mean(t,py,na.rm=T))
D=d%>%mutate(RR=O/E, L=qchisq(.025,2*O)/(2*E),U=qchisq(.975,2*O+2)/(2*E))
D[D$cancer2=="MDS","t"]=D[D$cancer2=="MDS","t"]+0.05
graphics.off()
quartz(width=7,height=4)
theme_set(theme_bw())
theme_update(legend.position = c(.92, .825),
             axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             axis.title.x=element_text(size=rel(1.0)),
             legend.title=element_text(size=rel(1)),
             legend.text=element_text(size=rel(1)),
             strip.text = element_text(size = rel(1.5)))
g=qplot(x=t,y=RR,data=D,col=cancer2,geom=c("line","point"),
        xlab="Years Since Dx of Lymphoid First Cancer",ylab="Relative Risk")
g=g+facet_grid(Sex~.,scales="free")+geom_abline(intercept=1, slope=0)
g = g + scale_color_grey(start = 0, end = 0.6)
g1 <- guide_legend("Second\nCancer")
g=g + guides(color=g1) 
g=g+  geom_errorbar(aes(ymin=L,ymax=U,width=.15))+scale_y_continuous(breaks=c(0,5,10,15))
print(g)
ggsave("~/Results/amlMDS/lymphoidFirst.eps")  
ggsave("~/Results/amlMDS/lymphoidFirst.png")  


