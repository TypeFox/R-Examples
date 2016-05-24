################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
############################## MAIN PROGRAM ####################################
# Session: working with the model
library(ZeBook)
# TODO : change the year of the weather you want (1999 to 2008)
# TODO : change the site in France (1 to 40)
weather=watbal.weather(working.year=2007, working.site=1)
# Define parameter values of the model function
watbal.factors = watbal.define.param()

# input variable describing the soil
# WP : Water content at wilting Point (cm^3.cm^-3)
soil.WP = 0.06
# FC : Water content at field capacity (cm^3.cm^-3)
soil.FC = 0.06+0.13
# WAT0 : Initial Water content (cm^3.cm^-3)
soil.WAT0 = soil.FC*400
pred.WAT = watbal.model(watbal.factors["nominal",],weather, soil.WP,soil.FC, soil.WAT0)

# Write output to a file (to open with notepad or excel)
options(digits=3)
#write.table(format(pred.WAT), file = "pred.WAT.csv", quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)
# Produce graphical output of the state variables
#png("watbal.example.png", width = 8.5, height = 2*8.5,units = "cm", res = 300,pointsize = 12)
par(mfrow=c(3,1), mar=c(4.1,4.1,1.1,0.2))
barplot(pred.WAT$RAIN, xlab = "day", ylab = "Rain (mm)" )
barplot(pred.WAT$ETr, xlab = "day", ylab = "ETr (mm)")
plot(pred.WAT$day,pred.WAT$WATp*100, xlab = "day", ylab = "WATp (%)",type="l",lwd=2,ylim=c(0,soil.FC*110))
abline(h=soil.FC*100, lty=2)
abline(h=soil.WP*100, lty=2)
#dev.off()
# End of file