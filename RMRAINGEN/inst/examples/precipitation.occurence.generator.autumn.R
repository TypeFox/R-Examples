# TODO: Add comment
# 
# Author: ecor
###############################################################################

rm(list=ls())

set.seed(789)

library(RMRAINGEN)


data(trentino)

####
additional_function_path  <- "/Users/ecor/Dropbox/iasma/RMAWGENdev/RMRAINGEN/inst/doc/examples/additional_functions"
list.files <- list.files(additional_function_path,pattern=".R",full.names=TRUE)
for (it in list.files) source(it)
#####


year_min <- 1961
year_max <- 1990

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max

period <- period & (PRECIPITATION$month %in% c(9,10,11)) ## AUTUMN 

station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
prec_mes <- PRECIPITATION[period,station]

## removing nonworking stations (e.g. time series with NA)
accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
}

prec_mes <- prec_mes[,accepted]
## the dateset is reduced!!!
#prec_mes <- prec_mes[,]

coeff2p <- CoeffYWeq(data=prec_mes,p=2,tolerance=1.e-10)
coeff1p <- CoeffYWeq(data=prec_mes,p=1,tolerance=1.e-10)


generation1p <- generate(coeff1p,n=nrow(prec_mes),names=names(prec_mes))
generation2p <- generate(coeff2p,n=nrow(prec_mes),names=names(prec_mes))

corg <- cor(generation1p)
CCGamma <- CCGamma(prec_mes,lag=0,only.matrix=TRUE,tolerance=1.e-10) ## DISTEMERE CCGAMMA BLOCKMATRIX
CCGammablock  <- CCGammaToBlockmatrix(prec_mes,lag=0,p=2,tolerance=1.e-10)
corg_mes <- CCGamma ##$nooccurence_gcorrelation
## DRAFT COR PLOT 

plot(corg,corg_mes,xlim=c(0.5,1),ylim=c(0.5,1),xlab="generated",ylab="observed")
abline(0,1)
corplot(x=NULL,y=generation1p,corx=corg_mes)


####




