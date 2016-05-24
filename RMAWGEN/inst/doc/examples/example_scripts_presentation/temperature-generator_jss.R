# file weather-generator.R
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
p_temp <- 10

exogen <- NULL
exogen_sim <- exogen

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



NDAY <- nrow(Tx_mes)	
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
	pdf_Tx <- paste(wpath,"/tx_qqplot_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tn <- paste(wpath,"/tn_qqplot_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_deltaT <- paste(wpath,"/dt_qqplot_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tx_anom <- paste(wpath,"/tx_anom_qqplot_",year_min,"_",year_max,"_",it,".pdf",sep="")
	pdf_Tn_anom <- paste(wpath,"/tn_anom_qqplot_",year_min,"_",year_max,"_",it,".pdf",sep="")
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
print("acf")
# ACF Function 
pdf(paste(wpath,"acf_tx_anom_P10GPCA.pdf",sep="/"))
plot(acf(Tx_gen$P10GPCA-Tx_spline,lag=50),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_tx_anom_mes.pdf",sep="/"))
plot(acf(Tx_mes-Tx_spline,lag=50))
dev.off()

pdf(paste(wpath,"acf_tn_anom_P10GPCA.pdf",sep="/"))
plot(acf(Tn_gen$P10GPCA-Tn_spline,lag=50),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_tn_anom_mes.pdf",sep="/"))
plot(acf(Tn_mes-Tn_spline,lag=50))
dev.off()

pdf(paste(wpath,"acf_deltat_P10GPCA.pdf",sep="/"))
plot(acf(Tx_gen$P10GPCA-Tn_gen$P10GPCA,lag=50),xlab="lag [day]")
dev.off()
pdf(paste(wpath,"acf_deltat_mes.pdf",sep="/"))
plot(acf(Tx_mes-Tn_mes,lag=50))
dev.off()

# COMPUTING CORRELATION OF GAUSSIANIZED VECTORS .... 

#qqplot(generationP01_temp$input$res_multigen,generationP01_temp$input$data_for_var)

#cor(generationP01_temp$output$res_multigen) #,generationP01_temp$input$data_for_var)
#cor(generationP01_temp$input$data_for_var)

#cor(generationP10_temp$input$data_for_var)
#VARselect(generationP01_temp$input$data_for_var,lag.max=20) 
#VARselect(generationP01GPCA_temp$var@GPCA_data$final_results,lag.max=20)
#.... 

## 
## pdf_Tx <-  '/Users/ecor/Dropbox/IASMA_CRI_DAES/beamerposter/EGU2012_EC/images_article/fig_Tx.pdf'
## pdf_Tn <- '/Users/ecor/Dropbox/IASMA_CRI_DAES/beamerposter/EGU2012_EC/images_article/fig_Tn.pdf'
## pdf_DeltaT <- '/Users/ecor/Dropbox/IASMA_CRI_DAES/beamerposter/EGU2012_EC/images_article/fig_DeltaT.pdf'
## #station0 <- "T0090"
## main <- names(Tx_gen)
## #source('/Users/ecor/Dropbox/iasma/RMAWGENdev/RMAWGEN/inst/doc/private/additional_functions/plot_RMAWGENts.R') 
## #source('/Users/ecor/Dropbox/iasma/RMAWGENdev/RMAWGEN/inst/doc/private/additional_functions/qqplot_RMAWGENts.R') 
## 
## season <- DJF
## 
## qqplot_RMAWGEN_Tn(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,station="T0090",when=season)
## qqplot_RMAWGEN_Tx(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,station="T0090",when=season)
## qqplot_RMAWGEN_DeltaT(Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tx_gen=Tx_gen,Tn_gen=Tn_gen,station="T0090",when=season)
## 
## #qqplot_RMAWGEN_Tx(Tx_mes=Tx_mes,Tx_gen=Tx_gen,Tn_mes=Tn_mes,Tn_gen=Tn_gen,pdf=pdf_Tx,station=station0)	
## #qqplot_RMAWGEN_Tn(Tx_mes=Tx_mes,Tx_gen=Tx_gen,Tn_mes=Tn_mes,Tn_gen=Tn_gen,pdf=pdf_Tn,station=station0)		
## #qqplot_RMAWGEN_deltaT(Tx_mes=Tx_mes,Tx_gen=Tx_gen,Tn_mes=Tn_mes,Tn_gen=Tn_gen,pdf=pdf_DeltaT,station=station0)		
## 