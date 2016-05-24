
######################################################
### Test for interval LOQ
### Plots to know if the functions works properly
######################################################

# rm(list=ls(all=TRUE))
# library(devtools)

# dev_mode()
# document()
# load_all()

# Load data
data(ecdata)
data(mfidata)
sdf <- data_selection(mfidata[mfidata$plate=="plate_1" & mfidata$analyte=="FGF",], ecdata)[[1]]

igmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
                      lfct="SSl4", bkg="ignore", fmfi="mfi", verbose=FALSE)

submodels <- scluminex("plate_1",sdf$standard, sdf$background, 
                      lfct="SSl4", bkg="subtract", fmfi="mfi", verbose=FALSE)

incmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
                      lfct="SSl4", bkg="include", fmfi="mfi", verbose=FALSE)

consmodels <- scluminex("plate_1",sdf$standard, sdf$background, 
                      lfct="SSl4", bkg="constrain", fmfi="mfi", verbose=FALSE)


ff <- function(model, title, lp, hp){
  ndf <- model$FGF$data
  
  back <- log10(model$FGF$bkg_mean)
  
  conf <- conf_bands(model,"FGF" ,
                     seq(min(ndf[,"log10_concentration"]),max(ndf[,"log10_concentration"]),0.1),interval="prediction")
  
  with(ndf, plot(log10_concentration, log10_mfi))
  with(conf, lines(xvalue, log10_mfi))
  with(conf, lines(xvalue, log10_mfi.lci,col="grey"))
  with(conf, lines(xvalue, log10_mfi.uci, col="grey"))
  
  fixassymp_loq <- loq_interval(model, low.asymp=lp, high.asymp=hp)
  fixvalue_loq <- loq_interval(model, lowci=2, highci=3)
  
  aux <- fixassymp_loq[[1]]
  abline(v=c(aux$lloq, aux$uloq), col="blue")
  abline(h=c(aux$ly, aux$uy), col="blue", lty=2)
  
  aux <- fixvalue_loq[[1]]
  abline(v=c(aux$lloq, aux$uloq), col="red")
  abline(h=c(aux$ly, aux$uy), col="red",lty=2)
  abline(h=c(2, 3), col="green",lty=2)
  
  abline(h=back)
  title(title)
}


# red: LOQs fix value
# blue: LOQs asymptote
# green: fixed value
# black: background
par(mfrow=c(2,2))
ff(igmodels,"ignore", 2, 3)
ff(submodels,"subtract", 2, 3)
ff(incmodels,"include", 2, 3)
ff(consmodels,"constraint", 2, 2)

####################
##### Another plot
###################
par(mfrow=c(1,1))
ndf <- igmodels$FGF$data
model <- igmodels
back <- log10(model$FGF$bkg_mean)
conf <- conf_bands(model, "FGF",
                   seq(min(ndf[,"log10_concentration"]),max(ndf[,"log10_concentration"]),0.1),interval="prediction")

with(ndf, plot(log10_concentration, log10_mfi, ylim=c(1,4)))
with(conf, lines(xvalue, log10_mfi))
with(conf, lines(xvalue, log10_mfi.lci,col="grey"))
with(conf, lines(xvalue, log10_mfi.uci, col="grey"))

fixassymp_loq <- loq_interval(model, low.asymp="c", high.asymp="d")
fixvalue_loq <- loq_interval(model, lowci=2, highci=3)

confcoef <- confint(model$FGF$model)
abline(v=fixassymp_loq$FGF$lloq)
abline(v=fixassymp_loq$FGF$uloq)

abline(h=fixassymp_loq$FGF$ly, lty=2,col="blue")
abline(h=fixassymp_loq$FGF$uy, lty=2,col="blue")

abline(h=confcoef[2,],lty=2,col="red")
abline(h=confcoef[3,],lty=2,col="red")



