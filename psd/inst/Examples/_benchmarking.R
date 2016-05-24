##
##
## test differences in processing time
## between psdcore and spec.mtm
##
## TODO: ensure results are equal, if possible??
##
##
library(ggplot2)
library(plyr)
library(multitaper)
library(rbenchmark)
library(reshape2)
library(psd)
#
reps <- 10
PSD <- psdcore
run.bench <- function(nd, nt=8, nreps=reps){
  set.seed(1234)
  message(c(nd, nt, nreps))
  #
  X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
  #
  psdbench <- benchmark(PSD(X.d, ntaper=nt, plot=FALSE, refresh=TRUE),
                        multitaper::spec.mtm(X.d,k=nt,plot=FALSE),
                        replications=nreps)
  psdbench$num_terms <- nd
  psdbench$num_taps <- nt
  psdbench$test <- c("psd::psdcore","multitaper::spec.mtm")
  return(psdbench)
}
#
nds <- round(10**seq.int(from=1,to=5,by=0.15))
# do the benchmarking
allbench <- lapply(X=nds, FUN=function(x) run.bench(nd=x))
# manipulate into data.frame
allbench.df <- plyr::ldply(allbench)
# save
save(allbench.df, file="current_benchmarking_data.rda")
tests <- unique(allbench.df$test)
# drop some info
allbench.df.drp <- subset(allbench.df, select = c(test, num_terms, user.self, sys.self, elapsed, relative))
# melt
allbench.df.mlt <- melt(allbench.df.drp, id.vars=c("test","num_terms"))
# summarize for each group
head(tmpd <- plyr::ddply(allbench.df.mlt, .(variable,num_terms), summarise, summary="maximums", value=max(value)))
# create data for summary curves
allmeds <- ldply(lapply(X=tests, FUN=function(x,df=tmpd){df$test <- x; return(df)}))
rm(tmpd)

g <- ggplot(data=allbench.df.mlt, 
            aes(x=log10(num_terms), 
                y=log2(value), 
                # 1/sqrt(n) standard errors [assumes N(0,1)]
                ymin=log2(value*(1-1/sqrt(reps))),
                ymax=log2(value*(1+1/sqrt(reps))),
                colour=test, 
                group=test)) + 
  scale_colour_discrete(guide="none") + geom_pointrange() 
g2 <- g + facet_grid(variable~test, scales="free_y")
print(g2 + geom_path(colour="black", data=allmeds, aes(group=test)))

