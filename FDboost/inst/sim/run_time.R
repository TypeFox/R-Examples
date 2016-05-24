###############################################################################
# run the simulations for boosting of functional data
# code based on code by Fabian Scheipl 
# author: Sarah Brockhaus
###############################################################################

##### test gain in time for computation as array model compared to normal computation

print(R.Version()$version.string)

#library(refund)
library(FDboost)
library(splines)
pathResults <- NULL

library(pryr)

source("simfuns.R")

##################################### M=50

set.seed(18102012)

settingsTime <- makeSettings(
  M=c(50, 100, 500),
  ni=c(1),
  Gy=c(30, 100),
  Gx=c(30),
  snrEps=c(2),
  snrE=c(0),
  snrB=c(2),
  scenario=3,
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1:10) # 1:10

length(settingsTime)

usecores <- 10
options(cores=usecores)
res <- try(doSimFDboost(settings=settingsTime, oneRepTime))

save(res, file=paste(pathResults, "resTime.Rdata", sep=""))

res$M <- factor(res$M)
res$Gy <- factor(res$Gy)

res$mod <- factor(res$model-2, levels=c(0,1), labels = c("FLAM", "long"))


#### generate boxplots of errors
library(ggplot2)

pdf("arrayLongTime.pdf", width=9, height=7)
p <- ggplot(res, aes(y=time.elapsed, x=M, fill=mod)) 
p + geom_boxplot(aes(colour = mod), outlier.size=.6) +
  facet_grid( ~ Gy, labeller="label_both") + 
  scale_colour_manual('mod', values = c('FLAM' = 'blue', 'long' = 'black')) +
  scale_fill_manual(values=c("lightblue", "white")) + 
  scale_y_continuous(breaks=c(1, 2, 5, 10, 20, 60, 120, 300, 600, 1200, 2700, 
                              5400, 10800, 21600, 43200), 
                     trans="log10",
                     labels=c("1s", "2s", "5s", "10s", "20s", "1 min", "2 min", "5 min", 
                              "10 min", "20 min", "45 min", "90 min", "3h", "6h", "12h")) +
  labs(title="Computation time", x="N", y="time") +
  theme(text=element_text(size = 20)) #+ theme_bw()
dev.off()


pdf("arrayLongMemory.pdf", width=9, height=7)
p <- ggplot(res, aes(y=memory, x=M, fill=mod)) 
p + geom_boxplot(aes(colour = mod), outlier.size=.6) +
  facet_grid( ~ Gy, labeller="label_both") + 
  scale_colour_manual('mod', values = c('FLAM' = 'blue', 'long' = 'black')) +
  scale_fill_manual(values=c("lightblue", "white")) + 
#   scale_y_continuous(breaks=c(1, 2, 5, 10, 20, 60, 120, 300, 600, 1200, 2700, 
#                               5400, 10800, 21600, 43200), 
#                      trans="log10",
#                      labels=c("1s", "2s", "5s", "10s", "20s", "1 min", "2 min", "5 min", 
#                               "10 min", "20 min", "45 min", "90 min", "3h", "6h", "12h")) +
  labs(title="Memory Consumption", x="N", y="memory") +
  theme(text=element_text(size = 20)) #+ theme_bw()
dev.off()

