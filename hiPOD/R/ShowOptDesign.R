ShowOptDesign <-
function(opt.design.results, num.designs=10)
{
is.valid.design<-upper.sample.good<-Xmean.good<-NULL
rm(is.valid.design, upper.sample.good, Xmean.good)	
designs.good <- subset(opt.design.results[[5]], subset=(is.valid.design & upper.sample.good & Xmean.good))
head(designs.good[order(-designs.good$pred.power), ], num.designs)
}

