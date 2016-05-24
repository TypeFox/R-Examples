##
## demo file for a variety of fan charts based on the th.mcmc object 
##
library("fanplot")

##
##install packages if not already done so (uncomment)
##
# install.packages("colorspace")
# install.packages("RColorBrewer")
# install.packages("zoo")
# install.packages("tsbugs")



# 1. Defaults.
par(mar=c(2,2,2,1.1))
# empty plot
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Defaults")
# add fan
fan(th.mcmc)



# 2. Coarser fan.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Coarser Palette")
fan(th.mcmc, probs = seq(10, 90, 10))



# 3. Fan with no lines and text.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Plain")
fan(th.mcmc, ln = NULL, rlab = NULL)



# 4. Black lines.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Black Contour Lines")
fan(th.mcmc, probs = seq(10, 90, 10), ln.col = "black")



# 5. Transparant fan fill to leave only percentile contor lines.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Transparent Fill, Black Lines")
fan(th.mcmc, alpha=0, ln.col="black")



# 6. Prediction interval.
plot(NULL, xlim = c(1, 985), ylim = range(th.mcmc)*0.85, main="Prediction Interval")
fan(th.mcmc, type = "interval")



# 7. Left labels too.
plot(NULL, xlim = c(-40, 985), ylim = range(th.mcmc)*0.85, main="Include Left Labels")
fan(th.mcmc, type = "interval", llab = TRUE)



# 8. Selected right labels.
plot(NULL, xlim = c(-20, 965), ylim = range(th.mcmc)*0.85, main="User Selected Labels")
fan(th.mcmc, llab=c(0.1,0.5,0.9), rlab=c(0.2,0.5,0.8))



# 9. Change prefixes of labels.
plot(NULL, xlim = c(-50, 995), ylim = range(th.mcmc)*0.85, main="User Named Labels")
fan(th.mcmc, type = "interval", llab = TRUE, rcex = 0.6,
    upplab = "Upp ", lowlab = "Low ", medlab="Med")



# 10. Different colour scheme.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="User Colour Ramp")
fan(th.mcmc, probs = seq(10, 90, 10), 
    fan.col = colorRampPalette(c("royalblue", "grey", "white")))



# 11. colorspace colours.
library("colorspace")
plot(NULL, xlim = c(0, 965), ylim = range(th.mcmc)*0.85, main="colorspace Library: diverge_hcl Colour Ramp")
fan(data = th.mcmc, fan.col = diverge_hcl) 



# 12. colorspace colours 2.
plot(NULL, xlim = c(0, 965), ylim = range(th.mcmc)*0.85, main="colorspace Library: sequential_hcl Colour Ramp")
fan(data = th.mcmc, fan.col = sequential_hcl ) 



# 13. RColorBrewer colours.
library("RColorBrewer")
plot(NULL, xlim = c(0, 965), ylim = range(th.mcmc)*0.85, main="RColorBrewer Library: Accent Colour Ramp")
fan(data = th.mcmc, 
    fan.col = colorRampPalette(colors = brewer.pal(8,"Accent")) ) 



# 14. RColorBrewer colours 2.
plot(NULL, xlim = c(0, 965), ylim = range(th.mcmc)*0.85, main="RColorBrewer Library: Oranges Colour Ramp")
fan(data = th.mcmc, 
    fan.col = colorRampPalette(colors = rev(brewer.pal(9,"Oranges"))) ) 



# 15. RColorBrewer colours 3.
plot(NULL, xlim = c(0, 965), ylim = range(th.mcmc)*0.85, main="RColorBrewer Library: Spectral Colour Ramp")
fan(data = th.mcmc, 
    fan.col = colorRampPalette(colors = brewer.pal(11,"Spectral")) ) 



# 16. Irregular time series.
library("zoo")
library("tsbugs")
th.mcmc2 <- zoo(th.mcmc, order.by=svpdx$date) 
plot(th.mcmc2[,1], type="n", ylim = range(th.mcmc)*0.85, main="X-Axis Based On (Irregular) Dates")
fan(data = th.mcmc2, rcex=0.5)



# 17. Spaghetti plot.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Default Spaghetti Plots")
fan(th.mcmc, style = "spaghetti")



# 18. More transparency in lines.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="Less Transparent Spaghetti Plots")
fan(th.mcmc, style = "spaghetti", alpha = 0.1)



# 19. More lines, red.
plot(NULL, xlim = c(1, 965), ylim = range(th.mcmc)*0.85, main="More Lines, Red")
fan(th.mcmc, style = "spaghetti", ln.col = "red", n.spag = 100, alpha = 0.1)



# 20. Overlay spaghetti on transparent fan.
plot(NULL, xlim = c(-20, 965), ylim = range(th.mcmc)*0.85, main="Spaghetti and Fan")
# transparent fan with visible lines
fan(th.mcmc, ln=c(5, 50, 95), llab=TRUE, alpha=0, ln.col="orange" )
# spaghetti lines
fan(th.mcmc, style="spaghetti")


##
##save as png
##
# dev.copy(png, file = "svplots.png", width=10, height=50, units="in", res=400)
# dev.off()

