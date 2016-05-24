################################################################################
# "Working with dynamic models for agriculture"
# R script for pratical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2010-08-09
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
############################## MAIN PROGRAM ####################################
# Session 1: working with the model
library(ZeBook)
# TODO : change the year of the weather you want (1984 to 2011)
# TODO : change the site in France (1 to 40)
weather=maize.weather(working.year=2010, working.site=1)
# Define parameter values of the model function
# Tbase  : the baseline temperature for growth (degC)
Tbase <- 7.0
# RUE : radiation use efficiency (g.MJ-1)
RUEmax <- 1.85
# K : extinction coefficient (-)
K <- 0.7
#alpha : the relative rate of leaf area index increase for small values of leaf area index ((degC.day)-1)
alpha <- 0.00243
#LAImax : maximum leaf area index (m2 leaf/m2 soil)
LAImax <- 7.0
#TTM :  temperature sum for crop maturity (degC.day)
TTM <- 1200
#TTL : temperature sum at the end of leaf area increase (degC.day)
TTL <- 700
# sdate : sowing date
sdate <- 100
# ldate : last date
ldate <- 250
# Running the model
output<-maize.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather,sdate,ldate)

# Write output to a file (to open with notepad or excel)
options(digits=3)
format(output)
#write.table(format(output), file = "output.csv", quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)
# Produce graphical output of the state variables
dev.new()
par(mfcol=c(3,2))
plot(output$day,output$TT, xlab = "day", ylab = "Temperature sum",type="l")
plot(output$day,output$B, xlab = "day", ylab = "Biomass",type="l")
plot(output$day,output$LAI, xlab = "day", ylab = "LAI",type="l")
# Produce graphical output of the input variables
plot(1:365,weather$Tmin, xlab = "day", ylab = "Temperature min")
plot(1:365,weather$Tmax, xlab = "day", ylab = "Temperature max")
plot(1:365,weather$I, xlab = "day", ylab = "Solar Radiation")

################################################################################
## Question 1
res = output[which(output["day"]==115),]
print(res)

res = output[which(output["day"]==250),]
print(res)

#question 2
weather1=maize.weather(working.year=2010, working.site=1)
weather2=maize.weather(working.year=2011, working.site=1)
output1=maize.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather1,sdate,ldate)
output2=maize.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather2,sdate,ldate)

#plotting the Biomass of the 2 years
dev.new()
plot(output1$day,output1$B, xlab = "day", ylab = "Biomass",type="l")
lines(output1$day,output2$B, xlab = "day", ylab = "Biomass",type="l",lty=2)
legend("topleft",c("weather1","weather2"),lty=c(1,2))

 which(output1$TT>TTM)

#computing date of maturity (as simulation day) for the 2 years
print(min(which(output1$TT>TTM)))
print(min(which(output2$TT>TTM)))

# additionnal plot
dev.new()
par(mfrow=c(3,1))
plot(output1$day,output1$TT, xlab = "day", ylab = "TT",type="l")
lines(output2$day,output2$TT, xlab = "TT", ylab = "TT",lty=2)
legend("topleft",c("weather1","weather2"),lty=c(1,2))
plot(output1$TT,output1$LAI, xlab = "TT", ylab = "LAI",type="l")
lines(output2$TT,output2$LAI, xlab = "TT", ylab = "LAI/TT",lty=2)
legend("topleft",c("weather1","weather2"),lty=c(1,2))
plot(output1$LAI,output1$B, xlab = "LAI", ylab = "B",type="l")
lines(output2$LAI,output2$B, xlab = "LAI", ylab = "B",lty=2)
legend("topleft",c("weather1","weather2"),lty=c(1,2))

#Question 4
output=maize_cir.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather1,sdate,ldate)
dev.new()
par(mfrow=c(2,1))
plot(output$day,output$CumInt, xlab = "day", ylab = "CumInt",col="green")
plot(output$CumInt,output$B, xlab = "CumInt", ylab = "B",col="green")
#estimation of RUE
c1 = 1500
c2 = 100 ;
i1 = min(which(output$CumInt>c1))
i2 = min(which(output$CumInt>c2))
RUE_estim = (output$B[i2]-output$B[i1])/(c2 - c1)
print(RUE_estim)

# Question 5 : RUE function

dev.new()
Tpt=seq(0,50,by=0.05)
plot(Tpt,maize.RUEtemp(Tpt,RUEmax,6.2,16.5,33,44),xlab="T",ylab="RUE",typ='l')

#generating results with the 2 models for the 2 years
output1=maize.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather1,sdate,ldate)
output1_rue=maize_cir_rue.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather1,sdate,ldate)
output2=maize.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather2,sdate,ldate)
output2_rue=maize_cir_rue.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather2,sdate,ldate)

#plot comparison
dev.new()
par(mfrow=c(2,2))
plot(output1$day,output1$B, xlab = "day (year 1)", ylab = "B",type="l",col="green");lines(output1_rue$day,output1_rue$B,type="l",col="red")
plot(output2$day,output2$B, xlab = "day (year 2)", ylab = "B",type="l",col="green");lines(output2_rue$day,output2_rue$B, type="l",col="red")
plot(output1$day,output1$LAI, xlab = "day (year 1)", ylab = "LAI",type="l",col="green");lines(output1_rue$day,output1_rue$LAI,type="l",col="red")
plot(output2$day,output2$LAI, xlab = "day (year 2)", ylab = "LAI",type="l",col="green");lines(output2_rue$day,output2_rue$LAI, type="l",col="red")


#Question 6

output_ear1=maize_cir_rue_ear.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather1,sdate,ldate)
output_ear2=maize_cir_rue_ear.model(Tbase,RUEmax,K,alpha,LAImax,TTM,TTL,weather2,sdate,ldate)

# plotting LAI, B and BE for year 1
dev.new()
par(mfrow=c(3,1))
plot(output_ear1$TT,output_ear1$LAI, xlab = "TT", ylab = "MAI",type="l");
plot(output_ear1$TT,output_ear1$BE, xlab = "TT", ylab = "BE",type="l");
plot(output_ear1$TT,output_ear1$B, xlab = "TT", ylab = "B",type="l");

# comparison of some ratios
print(max(output_ear1$BE)/max(output_ear1$B))
print(max(output_ear2$BE)/max(output_ear2$B))

output_ear1750=maize_cir_rue_ear.model(Tbase,RUEmax,K,alpha,LAImax,TTM,750,weather1,sdate,ldate)
print(max(output_ear1750$BE)/max(output_ear1750$B))

output_ear17501300=maize_cir_rue_ear.model(Tbase,RUEmax,K,alpha,LAImax,1300,750,weather1,sdate,ldate)
print(max(output_ear17501300$BE)/max(output_ear17501300$B))
# End