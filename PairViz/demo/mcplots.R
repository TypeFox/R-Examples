####### Code for constructing figures of Section 3.1 (Hurley and Oldford, 2008) ###################
	 			
library(PairViz)	 			
#############################

data(cancer, package="PairViz")

bx <- with(cancer, split(sqrt(Survival),Organ))
a <-  aov(sqrt(Survival) ~ Organ,data=cancer)

dev.new(height=3, width=3)
par(mar=c(3,4.8,2,1),cex.main=.75,mgp=c(1.5, .5, 0),cex.axis=.6,cex.lab=.75)


plot(TukeyHSD(a,conf.level = 0.95),las=1,tcl = -.3)


abline(v = 0, col = "grey60", lty=3) # should not need
abline(h = 1:10, col = "grey85", lty=1) # should not need


#############################
# library(HH)
# mm <- glht.mmc(a, linfct = mcp(Organ = "Tukey"))
# dev.new(height=4, width=5)




# dev.new(height=6.5, width=9)
# par(mgp=c(2, .5, 0),mar=c(3,4,3,8))

# plot(mm, x.offset=1,main="95% family-wise confidence level",main2.method.phrase="",xlab = "Differences in mean sqrt Survival",col.mca.signif="red")
# par(mgp=c(3, .5, 0))
# title(ylab="Mean sqrt Survival")


#############################


dev.new(height=4.5, width=9.5)
par(cex.axis=.75, cex.main = 1.0, cex.lab=1)
par(mar=c(3,5,3,5))


cols <- c("#FF8080" ,"#8080FF", "#80FFFF" ,"#EEEE77", "#FFDAE2")
mc_plot(bx,a,main="Pairwise comparisons of cancer types", ylab="Sqrt Survival",col=cols)


###########draw in stages ##########


dev.new(height=4.5, width=9.5)
par(cex.axis=.75, cex.main = 1.6, cex.lab=1.4)
par(mar=c(3,5,3,5))

# suppress CIs
mc_plot(bx,a,main="Pairwise comparisons of cancer types", ylab="Sqrt Survival",levels=NULL,col=cols)

# draw CIs

mc_plot(bx,a,main="Pairwise comparisons of cancer types", ylab="Sqrt Survival", sig.col=NULL,col=cols)

# mark significant CIs

mc_plot(bx,a,main="Pairwise comparisons of cancer types", ylab="Sqrt Survival",col=cols)
