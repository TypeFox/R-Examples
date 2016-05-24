# amlAbombAge.R  (Figure 2J)
library(dplyr)
# library(SEERaBomb)
# mkAbomb() #~\data\abomb\lsshempy.csv, lssinc07.csv=> tables heme and solid in abomb.db
db <- src_sqlite("~/data/abomb/abomb.db")
d=collect(tbl(db, sql("SELECT * from heme")))%>%   
  mutate(Dose=cut(D,c(-1,.01,.4,10),labels=c("Low","Medium","High"),include.lowest=TRUE)) %>%
  mutate(agec=cut(age,c(seq(0,80,20),110),labels=seq(10,90,20))) %>%
  group_by(Dose,agec) %>%
  summarise(age=mean(age),py=sum(py)/1e5,AML=sum(AMLtot),
            L=qpois(.025,AML),U=qpois(.975,AML),Incid=AML/py,IL=L/py,IU=U/py)
d
library(ggplot2)
graphics.off()
quartz(width=7,height=5)
theme_set(theme_bw(base_size = 18))
# theme_set(theme_gray(base_size = 18)) 
theme_update(legend.position = c(.7, .2),
             axis.text=element_text(size=rel(1.3)),
             axis.title=element_text(size=rel(1.3)),
             legend.title=element_text(size=rel(1.1)),
             legend.text=element_text(size=rel(1.1)))
g=qplot(x=age+0.5*(as.numeric(Dose)-1),y=Incid,data=d,col=Dose,log="y",ylab="Incidence (Cases/100,000 PY)",
        geom=c("line","point"),xlab="Age")
g = g + scale_color_grey(start = 0.8, end = 0)
g+  geom_errorbar(aes(ymin=IL,ymax=IU,width=.15))
ggsave("~/Results/amlMDS/abombAMLage.eps")  
ggsave("~/Results/amlMDS/abombAMLage.png")  

# Side calculation
(2.5/9)^0.5 # = whole body A-bomb dose equivalent of average radiation therapy AML carcinogenicity


