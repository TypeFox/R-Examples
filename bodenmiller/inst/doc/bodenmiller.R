## ----setup,echo=FALSE,include=FALSE--------------------------------------
library(knitr)
library(bodenmiller)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

knitr::opts_chunk$set(warning=FALSE,
                      fig.keep='high',
                      fig.align='center')

do.fan <- function(x,step=0.01) {
  data.frame(ymin=quantile(x,probs=seq(0,1,step))[-length(seq(0,1,step))],
             ymax=quantile(x,probs=seq(0,1,step))[-1],
             id=seq(1,length(seq(step,1,step))),
             percent=abs(seq(step,1,step)-0.5))
}

scale_fill_fan <- function(...) scale_fill_gradientn(colours=rev(brewer.pal(9,'Oranges')) )

## ----ref_pheno_boxplot---------------------------------------------------
data(refPhenoMat)
refPhenoFrame <- melt(refPhenoMat)
names(refPhenoFrame) <- c('cell_id','channel','value')
ggplot(data=refPhenoFrame,aes(x=channel,y=value))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----annots--------------------------------------------------------------
data('refAnnots')
refPhenoFrame$Cells <- rep(refAnnots$Cells,ncol(refPhenoMat))
cell.colors <- setNames(c('#9CA5D5','#0015C5','#5B6CB4','#BFC5E8','#C79ED0','#850094',
                          '#A567B1','#DBBCE2','#D3C6A1','#5E4500','#BBDEB1','#8A1923',
                          '#B35E62','#CEA191'),
                        c('cd14-hladr-','cd14-hladrhigh','cd14-hladrmid','cd14-surf-',
                          'cd14+hladr-','cd14+hladrhigh','cd14+hladrmid','cd14+surf-',
                          'cd4+','cd8+','dendritic','igm-','igm+','nk'))

## ----ref_pheno_pop_boxplot,fig.width=6,fig.height=4----------------------
cd7.pops <- refPhenoFrame %>% filter(channel=='CD7')
ggplot(data=cd7.pops,
       aes(x=Cells,y=value,fill=Cells))+
  geom_boxplot()+
  scale_fill_manual(values=cell.colors)+
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----ref_pheno_sub_fan---------------------------------------------------
ggplot(refPhenoFrame %>% filter(Cells=='cd4+') %>% group_by(Cells,channel) %>% do(do.fan(.$value)),
                   aes(x=channel,fill=percent,group=id))+
  geom_ribbon(aes(ymin=ymin,ymax=ymax))+
  guides(fill=F)+
  facet_wrap(~Cells)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_fan()

## ----ref_func_sub_fan,fig.width=5,fig.height=3---------------------------
data(refFuncMat)
refFuncFrame <- melt(refFuncMat)
names(refFuncFrame) <- c('cell_id','channel','value')
refFuncFrame$Cells <- rep(refAnnots$Cells,ncol(refFuncMat))
ggplot(refFuncFrame %>% filter(Cells=='cd4+') %>% group_by(Cells,channel) %>% do(do.fan(.$value)),
                   aes(x=channel,fill=percent,group=id))+
  geom_ribbon(aes(ymin=ymin,ymax=ymax))+
  facet_wrap(~Cells)+
  guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_fan()

## ----untreated_func------------------------------------------------------
data('untreatedFuncMat')
data('untreatedAnnots')
untreatedFuncFrame <- melt(untreatedFuncMat)
names(untreatedFuncFrame) <- c('cell_id','channel','value')
untreatedFuncFrame$Cells <- rep(untreatedAnnots$Cells,ncol(untreatedFuncMat))
untreatedFuncFrame$Treatment <- rep(untreatedAnnots$Treatment,ncol(untreatedFuncMat))

## ----un_func_sub_fan,fig.width=6,fig.height=6----------------------------
refFuncLine <- refFuncFrame %>% filter(Cells=='cd4+') %>% group_by(Cells,channel) %>% summarise(value=median(value))
refFuncLine <- do.call(rbind,lapply(seq(1,3),function(x) refFuncLine))
refFuncLine$Treatment <- rep(levels(untreatedFuncFrame$Treatment),each=nlevels(refFuncLine$channel))
refFuncLine$percent <- 0
refFuncLine$id <- 'ref'
ggplot(untreatedFuncFrame %>% filter(Cells=='cd4+') %>% group_by(Treatment,Cells,channel) %>% do(do.fan(.$value)),
                   aes(x=channel,fill=percent,group=id))+
  geom_ribbon(aes(ymin=ymin,ymax=ymax))+
  geom_line(data=refFuncLine,aes(y=value),
            col='black',linetype=4)+
  guides(fill=F)+
  facet_wrap(~Cells*Treatment,ncol=2,scale='free_x')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_fan()

