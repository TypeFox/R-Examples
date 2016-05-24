## ----results="hide",message=FALSE----------------------------------------
require(eyelinker)
require(dplyr)

#Look for file 
fpath <- system.file("extdata/mono500.asc.gz",package="eyelinker")

## ------------------------------------------------------------------------
dat <- read.asc(fpath)

## ------------------------------------------------------------------------
names(dat)

## ------------------------------------------------------------------------
str(dat$info)

## ------------------------------------------------------------------------
raw <- dat$raw
head(raw,3)

## ------------------------------------------------------------------------
dat.bi <- system.file("extdata/bino1000.asc.gz",package="eyelinker") %>% read.asc

head(dat.bi$raw,3)

## ------------------------------------------------------------------------
library(tidyr)

raw.long <- dplyr::select(raw,time,xp,yp,block) %>% gather("coord","pos",xp,yp)
head(raw.long,2)
tail(raw.long,2)

## ----fig.width=5, fig.height=5-------------------------------------------
require(ggplot2)
raw.long <- mutate(raw.long,ts=(time-min(time))/1e3) #let's have time in sec. 
ggplot(raw.long,aes(ts,pos,col=coord))+geom_point()

## ----fig.width=5, fig.height=5-------------------------------------------
ggplot(raw.long,aes(ts,pos,col=coord))+geom_line()+facet_wrap(~ block)

## ------------------------------------------------------------------------
sac <- dat$sac
head(sac,2)

## ------------------------------------------------------------------------
head(dat.bi$sac,3)

## ------------------------------------------------------------------------
Sac <- cbind(sac$stime,sac$etime) #Define a set of intervals with these endpoints
#See also: intervals package
raw <- mutate(raw,saccade=time %In% Sac)
head(raw,3)
mean(raw$saccade)*100 #6% of time samples correspond to saccades

## ----fig.width=5, fig.height=5-------------------------------------------
mutate(raw.long,saccade=time %In% Sac) %>% filter(block==1) %>% ggplot(aes(ts,pos,group=coord,col=saccade))+geom_line()

## ------------------------------------------------------------------------
fix <- dat$fix
head(fix,3)

## ----fig.width=5, fig.height=5-------------------------------------------
Fix <- cbind(fix$stime,fix$etime) #Define a set of intervals 
mutate(raw.long,fixation=time %In% Fix) %>% filter(block==1) %>% ggplot(aes(ts,pos,group=coord,col=fixation))+geom_line()

## ------------------------------------------------------------------------
mutate(raw,fix.index=whichInterval(time,Fix)) %>% head(4)

## ------------------------------------------------------------------------
raw <- mutate(raw,fix.index=whichInterval(time,Fix))
fix.check <- filter(raw,!is.na(fix.index)) %>% group_by(fix.index) %>% summarise(axp=mean(xp),ayp=mean(yp)) %>% ungroup
head(fix.check,3)

## ------------------------------------------------------------------------
all.equal(fix.check$axp,fix$axp)
all.equal(fix.check$ayp,fix$ayp)

## ------------------------------------------------------------------------
fpath <- system.file("extdata/monoRemote500.asc.gz",package="eyelinker")
dat <- read.asc(fpath)
dat$blinks

## ------------------------------------------------------------------------
Blk <- cbind(dat$blinks$stime,dat$blinks$etime) #Define a set of intervals

filter(dat$raw,time %In% Blk) %>% head

## ------------------------------------------------------------------------
require(intervals)
Suspect <- Intervals(Blk) %>% expand(100,"absolute")
Suspect

## ----fig.width=5, fig.height=5-------------------------------------------
raw.long <- dplyr::select(dat$raw,time,xp,yp,block) %>% gather("coord","pos",xp,yp)
raw.long <- mutate(raw.long,ts=(time-min(time))/1e3) #let's have time in sec. 
ex <- mutate(raw.long,suspect=time %In% Suspect) %>% filter(block==2) 
ggplot(ex,aes(ts,pos,group=coord,col=suspect))+geom_line()+coord_cartesian(xlim=c(34,40))+labs(x="time (s)")

## ------------------------------------------------------------------------
head(dat$msg)

## ------------------------------------------------------------------------
library(stringr)
filter(dat$msg,str_detect(text,fixed("blank_screen"))) 

