# TODO: Add comment
# 
# Author: ecor
###############################################################################



rm(list=ls())
library(RMAWGEN)
library(RGENERATE)
set.seed(1222)

data(trentino)
## path where to plot the output figures
wpath <-  "/home/ecor/github/RMAWGENCodeCorner/data"
## ADJUST DATASET 


year_min <- 1961
year_max <- 1990

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max
station_prec <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
station_temp <- names(TEMPERATURE_MAX)[!(names(TEMPERATURE_MAX) %in% c("day","month","year"))]
prec_mesx <- PRECIPITATION[period,station_prec]
temp_mesx <- TEMPERATURE_MAX[period,station_temp]

names <- intersect(station_prec,station_temp)

## removing nonworking stations (e.g. time series with NA)
accepted <- array(TRUE,length(names))
names(accepted) <- names
for (it in names) {
	##cond  <- (length(which(!is.na(prec_mesx[,it])))==length(prec_mesx[,it]))
	cond <- TRUE
	accepted[it]  <- (length(which(!is.na(temp_mesx[,it])))==length(temp_mesx[,it])) & cond 
}


station <- names[accepted]
print(station)






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

#p_test <- 1
#p_prec <- 3
#
#exogen <- NULL
#exogen_sim <- exogen
#
#generationP03GPCA_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_prec,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
#generationP01GPCA_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=n_GPCA_iter,n_GPCA_iteration_residuals=n_GPCA_iteration_residuals,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
#
#
#generationP03_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_prec,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
#generationP01_prec <- ComprehensivePrecipitationGenerator(station=station,prec_all=PRECIPITATION,year_min=year_min,year_max=year_max,p=p_test,n_GPCA_iteration=0,n_GPCA_iteration_residuals=0,exogen=exogen,exogen_sim=exogen_sim,sample="monthly",mean_climate_prec=PREC_CLIMATE,no_spline=FALSE)
#
## VAR select 
#
#VARselect(generationP03_prec$data_prec,lag.max=20) 
#VARselect(generationP03GPCA_prec$var@GPCA_data$final_results,lag.max=20)
#
#
#normality_test(generationP01_prec$var)
#normality_test(generationP03_prec$var)
#normality_test(generationP01GPCA_prec$var)
#normality_test(generationP03GPCA_prec$var)
#
#
#serial_test(generationP01_prec$var)
#serial_test(generationP03_prec$var)
#serial_test(generationP01GPCA_prec$var)
#serial_test(generationP03GPCA_prec$var)
#
## ORGANIZZARE PRECIPITAZIONI COME PER LE TEMPERATURE!!! 
#
#
## Collecting the measured and generated time series 
#
#prec_mes <- generationP01_prec$prec_mes
#
#
#prec_gen <- list(P03GPCA=generationP03GPCA_prec$prec_gen,
#		P01GPCA=generationP01GPCA_prec$prec_gen,
#		P03=generationP03GPCA_prec$prec_gen,
#		P01=generationP01GPCA_prec$prec_gen)
## season 
#
#
#NDAY <- nrow(prec_mes)	
#days <- list()
#days$DJF <-  extractmonths(data=1:NDAY,when=c("Dec","Jan","Feb"),origin=origin)
#days$MAM <- extractmonths(data=1:NDAY,when=c("Mar","Apr","May"),origin=origin)
#days$JJA <- extractmonths(data=1:NDAY,when=c("Jun","Jul","Aug"),origin=origin)
#days$SON <- extractmonths(data=1:NDAY,when=c("Sep","Oct","Nov"),origin=origin)	  
#
#
## SET THE CORRECT PATH WHERE TO PLOT THE FIGURES 
#wpath <-  "/Users/ecor/Dropbox/iasma/RMAWGENdev/RMAWGEN/inst/doc/private/elsevier/article/Rscript_2/images_plural_stations" ##  "./"
#station00 <- "T0090"
#CEX <- 1.4
#
#for (it in names(days)) {
#	
#	str(it)
#	name <- it
#	season <- days[[it]]
#	lag <- 1
#	pdf_prec <- paste(wpath,"/prec_qqplot_",lag,"_",year_min,"_",year_max,"_",it,".pdf",sep="")
#	main_prec  <- paste("prec",names(prec_gen),station00,lag,"days",it,sep=" ")
#	qqplot_RMAWGEN_prec(prec_mes=prec_mes,prec_gen=prec_gen,main=main_prec,station=station00,when=season,pdf=pdf_prec,lag=lag,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
#
#	lag <- 2
#	pdf_prec <- paste(wpath,"/prec_qqplot_",lag,"_",year_min,"_",year_max,"_",it,".pdf",sep="")
#	main_prec  <- paste("prec",names(prec_gen),station00,lag,"days",it,sep=" ")
#	qqplot_RMAWGEN_prec(prec_mes=prec_mes,prec_gen=prec_gen,main=main_prec,station=station00,when=season,pdf=pdf_prec,lag=lag,xlim=range(prec_mes)*lag,cex.main=CEX,cex.lab=CEX,cex.axis=CEX)
#	
#	
#	
#}
#
#
#
## ACF Function 
#pdf(paste(wpath,"acf_prec_P03GPCA.pdf",sep="/"))
#plot(acf(prec_gen$P03GPCA,lag=10),xlab="lag [day]")
#dev.off()
#pdf(paste(wpath,"acf_prec_mes.pdf",sep="/"))
#plot(acf(prec_mes,lag=10))
#dev.off()

# GENERATION OF DAILY TEMPERATURE COUPLED WITH PRECIPITATION WITH PREVIOUSLY GENERATED


n_GPCA_iter <- 5
n_GPCA_iteration_residuals <- 5

p_test <- 1
p_temp <- 6

# Exogenous variables : Daily Precipitation (Observed and P03GPCA)  

##exogen <- normalizeGaussian_severalstations(x=prec_mes,data=prec_mes,sample="monthly",origin_x=origin,origin_data=origin,step=0)
##exogen_sim <- normalizeGaussian_severalstations(x=prec_gen$P03GPCA,data=prec_gen$P03GPCA,sample="monthly",origin_x=origin,origin_data=origin,step=0)

exogen <- NULL
exogen_sim <- NULL
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
		P10=generationP10_temp$output$Tx_gen,
		P01=generationP01_temp$output$Tx_gen)



Tn_gen <- list(P10GPCA=generationP10GPCA_temp$output$Tn_gen,
		P01GPCA=generationP01GPCA_temp$output$Tn_gen,
		P10=generationP10_temp$output$Tn_gen,
		P01=generationP01_temp$output$Tn_gen)





## save inpts 
results <- list(Tx_gen=Tx_gen,Tn_gen=Tn_gen,Tx_mes=Tx_mes,Tn_mes=Tn_mes,Tn_spline=Tn_spline,Tx_spline=Tx_spline,prec_mes=prec_mes,prec_gen=prec_gen,year_min=year_min,year_max=year_max,
		data_tempGPCA=generationP10GPCA_temp$var@GPCA_data$final_results,
		data_temp=generationP01_temp$input$data_for_var)
filename <- paste(wpath,"results_uncoupled_temperature_generator_P10.rda",sep="/")
save(results,file=filename)








