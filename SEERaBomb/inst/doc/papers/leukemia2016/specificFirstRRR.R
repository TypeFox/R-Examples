# specificFirstRRR.R  (Radiation Relative Risk variant of specificFirst.R)
rm(list=ls())
library(SEERaBomb)
library(reshape2)
system.time(load("~/Results/amlMDS/pm.RData")) # 4 secs to load. 
system.time(load("~/Results/amlMDS/pf.RData")) # 4 secs to load.
.simpleCap <- function(s) paste(toupper(substring(s, 1, 1)), substring(s, 2),sep = "", collapse = " ")

names(pm$L)
hires=FALSE
hires=TRUE
if (hires) {
  brks=c(0,0.25,0.5,0.75,1,1.5,2,2.5,3,4,5,6,8,10,12)
  FirstS=c("prostate","breast") # first cancers to be plotted 
}   else {  # intermediate level resolution
  FirstS=c("NHL","lung")
  brks=c(0,1,2,3,6,9,12)  
}
dm=mkDF(pm,brks)
df=mkDF(pf,brks)
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(cancer1%in%FirstS)%>%select(Sex,cancer1,cancer2,trt,O,E,t,int)
d=d%>%filter(!((cancer1=="breast")&(Sex=="Male"))) # remove male breast cancer cases
head(d)
tail(d)



mkRRR=function(d,first,Ncol=5000) {
#   first="prostate"
  d0=d%>%filter(cancer1==first,trt=="noRad")
  di=d%>%filter(cancer1==first,trt=="rad")
  O=d0$O
  E=d0$E
  Oi=di$O
  Ei=di$E
  Xi=matrix(0,ncol=Ncol,nrow=length(Oi))
  X=matrix(0,ncol=Ncol,nrow=length(O))
  for (i in 1:Ncol) Xi[,i]=rpois(n=length(Oi),lambda=Oi)
  for (i in 1:Ncol) X[,i]=rpois(n=length(O),lambda=O)
  X=(Xi/X)*(E/Ei)%*%t(rep(1,Ncol))
  T=t(apply(X,1,function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm=TRUE)))
  data.frame(Sex=d0$Sex,cancer1=d0$cancer1,cancer2=d0$cancer2,RR=T[,2],L=T[,1],U=T[,3],t=d0$t)
}

graphics.off()
quartz(width=7,height=4)
theme_set(theme_bw())
theme_update(legend.position = "none")
theme_update(axis.text=element_text(size=rel(1.2)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(0.9)),
             legend.text=element_text(size=rel(0.9)),
             strip.text = element_text(size = rel(1.5)))

for (First in FirstS) {
  (D=mkRRR(d,First))
  g=qplot(x=t,y=RR,data=D,geom=c("line","point"),#xlim=c(-.1,24),
          xlab=paste("Years Since",.simpleCap(First),"First Cancer"),ylab="Radiation Relative Risk Ratio")
  if (First=="breast"|First=="prostate") g=g+facet_grid(cancer2~.,scales="free")+geom_abline(intercept=1, slope=0) else
    g=g+facet_grid(Sex~cancer2,scales="free")+geom_abline(intercept=1, slope=0) #g=g+scale_y_continuous(breaks=1:10)
  g=g+geom_errorbar(aes(ymin=L,ymax=U,width=.15))
  print(g)
  ggsave(paste0("~/Results/amlMDS/",First,ifelse(hires,"hires",""),"RRR.png"))
  ggsave(paste0("~/Results/amlMDS/",First,ifelse(hires,"hires",""),"RRR.eps"))
} # for loop on first types

