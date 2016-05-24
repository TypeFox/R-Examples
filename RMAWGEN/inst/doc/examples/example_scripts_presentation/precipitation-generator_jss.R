# file precipitation-generator.R
# 
# This file contains a script example with two coupled temperature and precipitation stochastic generations 
#
#
# author: Emanuele Cordano on 12-01-2012
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################




rm(list=ls())
library(RMAWGEN)
set.seed(1222)

data(trentino)


station <- c("T0090","T0083") #,"T0099","T0001")
# Calibration period 
# MONTHLY CLIMATOLOGY

TX_CLIMATE <- NULL #Tx_1961_1990[,station]
TN_CLIMATE <- NULL #Tn_1961_1990[,station]
PREC_CLIMATE <- NULL #prec_1961_1990[,station] # NULL # Adjusts prec_1961_1990 with days!!!! 

# Calibration period 
year_max <- 1990
year_min <- 1961
origin <- "1961-1-1"

# Simulation period (Stochastic Generation)
# MONTHLY CLIMATOLOGY

# specific parameter for model calibration 

n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5

p_test <- 1
p_prec <- 3

exogen <- NULL
exogen_sim <- exogen

generationP03GPCA_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_prec,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
generationP01GPCA_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)


generationP03_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_prec,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
generationP01_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)

# VAR select 

VARselect(generationP03_prec$data_prec,lag.max=20) 
VARselect(generationP03GPCA_prec$var@GPCA_data$final_results,lag.max=20)


normality_test(generationP01_prec$var)
normality_test(generationP03_prec$var)
normality_test(generationP01GPCA_prec$var)
normality_test(generationP03GPCA_prec$var)


serial_test(generationP01_prec$var)
serial_test(generationP03_prec$var)
serial_test(generationP01GPCA_prec$var)
serial_test(generationP03GPCA_prec$var)

# ORGANIZZARE PRECIPITAZIONI COME PER LE TEMPERATURE!!! 


# Collecting the measured and generated time series 

prec_mes <- generationP01_prec$prec_mes


prec_gen <- list(P03GPCA=generationP03GPCA_prec$prec_gen,
		P01GPCA=generationP01GPCA_prec$prec_gen,
		P03=generationP03GPCA_prec$prec_gen,
		P01=generationP01GPCA_prec$prec_gen)
# season 


NDAY <- nrow(prec_mes)	
days <- list()
days$DJF <-  extractmonths(data=1:NDAY,when=c("Dec","Jan","Feb"),origin=origin)
days$MAM <- extractmonths(data=1:NDAY,when=c("Mar","Apr","May"),origin=origin)
days$JJA <- extractmonths(data=1:NDAY,when=c("Jun","Jul","Aug"),origin=origin)
days$SON <- extractmonths(data=1:NDAY,when=c("Sep","Oct","Nov"),origin=origin)	  


# SET THE CORRECT PATH WHERE TO PLOT THE FIGURES 
wpath <- "./"
station00 <- "T0090"
CEX <- 1.4

for (it in names(days)) {
	
	str(it)
	name <- it
	season <- days[[it]]
	lag <- 1
	pdf_prec <- paste(wpath,"/prec_qqplot_",lag,"_",year_min,"_",year_max,"_",it,".pdf",sep="")
	main_prec  <- paste("prec",names(prec_gen),station00,lag,"days",it,sep=" ")
	qqplot_RMAWGEN_prec(prec_mes=prec_mes,prec_gen=prec_gen,main=main_prec,station=station00,when=season,pdf=pdf_prec,lag=lag,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)

	lag <- 2
	pdf_prec <- paste(wpath,"/prec_qqplot_",lag,"_",year_min,"_",year_max,"_",it,".pdf",sep="")
	main_prec  <- paste("prec",names(prec_gen),station00,lag,"days",it,sep=" ")
	qqplot_RMAWGEN_prec(prec_mes=prec_mes,prec_gen=prec_gen,main=main_prec,station=station00,when=season,pdf=pdf_prec,lag=lag,xlim=range(prec_mes)*lag,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	
	
	
}



# ACF Function 
pdf(paste(wpath,"acf_prec_P03GPCA.pdf",sep="/"))
plot(acf(prec_gen$P03GPCA,lag=10),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_prec_mes.pdf",sep="/"))
plot(acf(prec_mes,lag=10))
dev.off()

# GENERATION OF DAILY TEMPERATURE COUPLED WITH PRECIPITATION WITH PREVIOUSLY GENERATED


n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5

p_test <- 1
p_temp <- 10

# Exogenous variables : Daily Precipitation (Observed and P03GPCA)  

exogen <- normalizeGaussian_severalstations(x=prec_mes,data=prec_mes,sample="monthly",origin_x=origin,origin_data=origin,step=0)
exogen_sim <- normalizeGaussian_severalstations(x=prec_gen$P03GPCA,data=prec_gen$P03GPCA,sample="monthly",origin_x=origin,origin_data=origin,step=0)

#exogen <- prec_mes
#exogen_sim <- prec_gen$P03GPCA

generationP10GPCA_temp <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,p=p_temp,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_Tn=TN_CLIMATE,mean_climate_Tx=TX_CLIMATE)
generationP01GPCA_temp <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_Tn=TN_CLIMATE,mean_climate_Tx=TX_CLIMATE)


generationP10_temp <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,p=p_temp,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_Tn=TN_CLIMATE,mean_climate_Tx=TX_CLIMATE)
generationP01_temp <- ComprehensiveTemperatureGenerator(station=station,Tx_all=TEMPERATURE_MAX,Tn_all=TEMPERATURE_MIN,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_Tn=TN_CLIMATE,mean_climate_Tx=TX_CLIMATE)

# VAR select 

VARselect(generationP01_temp$input$data_for_var,lag.max=20) 
VARselect(generationP01GPCA_temp$var@GPCA_data$final_results,lag.max=20)


normality_test(generationP01_temp$var)
normality_test(generationP10_temp$var)
normality_test(generationP01GPCA_temp$var)
normality_test(generationP10GPCA_temp$var)


serial_test(generationP01_temp$var)
serial_test(generationP10_temp$var)
serial_test(generationP01GPCA_temp$var)
serial_test(generationP10GPCA_temp$var)


# Collecting the measured and generated time series 

Tn_mes <- generationP01_temp$input$Tn_mes
Tx_mes <- generationP01_temp$input$Tx_mes
Tx_spline <- generationP01_temp$input$Tx_spline
Tn_spline <- generationP01_temp$input$Tn_spline

Tx_gen <- list(P10GPCA=generationP10GPCA_temp$output$Tx_gen,
		P01GPCA=generationP01GPCA_temp$output$Tx_gen,
		P10=generationP10GPCA_temp$output$Tx_gen,
		P01=generationP01GPCA_temp$output$Tx_gen)



Tn_gen <- list(P10GPCA=generationP10GPCA_temp$output$Tn_gen,
		P01GPCA=generationP01GPCA_temp$output$Tn_gen,
		P10=generationP10GPCA_temp$output$Tn_gen,
		P01=generationP01GPCA_temp$output$Tn_gen)



for (it in names(days)) {
	
	str(it)
	name <- it
	season <- days[[it]]
	pdf_Tx <- paste(wpath,"/tx_qqplot_coupled_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tn <- paste(wpath,"/tn_qqplot_coupled_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_deltaT <- paste(wpath,"/dt_qqplot_coupled_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tx_anom <- paste(wpath,"/tx_anom_qqplot_coupled_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tn_anom <- paste(wpath,"/tn_anom_qqplot_coupled_",year_min,"_",year_max,"_",it,".pdf",sep="")
	main_tx  <- paste("Tx",names(Tx_gen),station00,it,sep=" ")
	main_tn  <- paste("Tn",names(Tn_gen),station00,it,sep=" ")
	main_deltat  <- paste("dT",names(Tx_gen),station00,it,sep=" ")
	main_tx_anom  <- paste("Tx_anom",names(Tx_gen),station00,it,sep=" ")
	main_tn_anom  <- paste("Tn_anom",names(Tn_gen),station00,it,sep=" ")
	
	qqplot_RMAWGEN_Tx(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,main=main_tx,station=station00,when=season,pdf=pdf_Tx,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	qqplot_RMAWGEN_Tn(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,main=main_tn,station=station00,when=season,pdf=pdf_Tn,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	qqplot_RMAWGEN_deltaT(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,main=main_deltat,station=station00,when=season,pdf=pdf_deltaT,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	qqplot_RMAWGEN_Tx(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,Tx_spline=Tx_spline,Tn_spline=Tn_spline,main=main_tx_anom,station=station00,when=season,pdf=pdf_Tx_anom,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	qqplot_RMAWGEN_Tn(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,Tx_spline=Tx_spline,Tn_spline=Tn_spline,main=main_tn_anom,station=station00,when=season,pdf=pdf_Tn_anom,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
	
}

# ACF Function 
pdf(paste(wpath,"acf_coupled_tx_anom_P10GPCA.pdf",sep="/"))

plot(acf(Tx_gen$P10GPCA-Tx_spline,lag=50),xlab="lag [day]",cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
dev.off()
pdf(paste(wpath,"acf_coupled_tx_anom_mes.pdf",sep="/"))
plot(acf(Tx_mes-Tx_spline,lag=50))
dev.off()

pdf(paste(wpath,"acf_coupled_tn_anom_P10GPCA.pdf",sep="/"))
plot(acf(Tn_gen$P10GPCA-Tn_spline,lag=50),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_coupled_tn_anom_mes.pdf",sep="/"))
plot(acf(Tn_mes-Tn_spline,lag=50),cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
dev.off()

pdf(paste(wpath,"acf_coupled_deltat_P10GPCA.pdf",sep="/"))
plot(acf(Tx_gen$P10GPCA-Tn_gen$P10GPCA,lag=50),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_coupled_deltat_mes.pdf",sep="/"))
plot(acf(Tx_mes-Tn_mes,lag=50),cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
dev.off()

print("acf PUT PRECIPITAION!!!!")
dT_mes <- Tx_mes - Tn_mes 
dT_gen <- Tx_gen$P10GPCA-Tn_gen$P10GPCA



names(Tx_gen$P10GPCA) <- paste("Tx",names(Tx_gen$P10GPCA),sep="_")
names(Tx_spline)      <- paste("Tx",names(Tx_spline),sep="_")
names(Tx_mes)         <- paste("Tx",names(Tx_mes),sep="_")


names(Tn_gen$P10GPCA) <- paste("Tn",names(Tn_gen$P10GPCA),sep="_")
names(Tn_spline)      <- paste("Tn",names(Tn_spline),sep="_")
names(Tn_mes)         <- paste("Tn",names(Tn_mes),sep="_")

# plottare auto-covarianze !!!! (VEDI SOPRA) ACF Function 
names(prec_mes)         <- paste("prec",names(prec_mes),sep="_")
names(prec_gen$P03GPCA) <- paste("prec",names(prec_gen$P03GPCA),sep="_")


names(dT_mes)         <- paste("dT",names(dT_mes),sep="_")
names(dT_gen)         <- paste("dT",names(dT_gen),sep="_")

# Auto-Covariance daily thermal range/precipitation 


dT_station <- "dT_T0090"
prec_station <- "prec_T0090"
val_gen <- as.data.frame(cbind(prec_gen$P03GPCA[,prec_station],dT_gen[,dT_station]))		
val_mes <- as.data.frame(cbind(prec_mes[,prec_station],dT_mes[,dT_station]))

names(val_mes) <- c(prec_station,dT_station)
names(val_gen) <- c(prec_station,dT_station)
pdf(paste(wpath,"acf_prec_dt_anom_P10GPCA.pdf",sep="/"))
plot(acf(val_gen,lag=50),xlab="lag [day]",cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
dev.off()

pdf(paste(wpath,"acf_prec_dt_anom_mes.pdf",sep="/"))
plot(acf(val_mes,lag=50),xlab="lag [day]",cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
dev.off()



