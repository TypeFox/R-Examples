
###################################################################### 
## Very simple test Monotonic increasing-decreasing, invest function
######################################################################

# Load data
data(mfidata)
data(ecdata)
ndf <- data_selection(mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="RANTES",], ecdata)
ndf <- ndf[[1]]$standard
ndf$log10_concentration <- log10(ndf$ec)
ndf$log10_mfi <- log10(ndf$mfi)


# Increasing data
model1 <- scluminex("xxx",ndf, NULL, lfct="SSl5", "ignore", fmfi="mfi", fec="ec")
conf1 <- conf_bands(model1, analyte="RANTES", 2.7)
inv1 <- invest(model1,"RANTES", 2.7)           

# Same decreasing data
ndf2 <- ndf
ndf2[,"mfi"] <- 10^(-log10(ndf2[,"mfi"]))
model2 <- scluminex("xxx",ndf2,NULL, lfct="SSl5", "ignore", fmfi="mfi", fec="ec")
conf2 <- conf_bands(model2,"RANTES", -2.7)
inv2 <- invest(model2, "RANTES", -2.7)           

                
# Plots
plot(model1,"sc")
plot(model2,"sc")

loq_cv(model1)
loq_cv(model2)

loq_derivatives(model1)
loq_derivatives(model2)

loq_interval(model1, low.asymp=2)
loq_interval(model2, low.asymp=2)




