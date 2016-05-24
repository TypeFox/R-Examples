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
origin <- paste(year_min,1,1,sep="-")

period <- PRECIPITATION$year>=year_min & PRECIPITATION$year<=year_max

period <- period 

station <- names(PRECIPITATION)[!(names(PRECIPITATION) %in% c("day","month","year"))]
prec_mes <- PRECIPITATION[period,station]

## removing nonworking stations (e.g. time series with NA)
accepted <- array(TRUE,length(names(prec_mes)))
names(accepted) <- names(prec_mes)
for (it in names(prec_mes)) {
	accepted[it]  <- (length(which(!is.na(prec_mes[,it])))==length(prec_mes[,it]))
}

prec_mes <- prec_mes[,accepted]
station <- names(prec_mes)


## the dateset is reduced!!!
#prec_mes <- prec_mes[,]

coeff2p <- CoeffYWeq(data=prec_mes,p=2,tolerance=1.e-10)
coeff1p <- CoeffYWeq(data=prec_mes,p=1,tolerance=1.e-10)


generation1p <- generate(coeff1p,n=nrow(prec_mes),names=names(prec_mes))
generation2p <- generate(coeff2p,n=nrow(prec_mes),names=names(prec_mes))



##### 


#### DRAFT CORRELATION PLOT OF ALL YEARS 

#corg <- cor(generation1p)
CCGamma <- CCGamma(prec_mes,lag=0,only.matrix=TRUE,tolerance=1.e-10) ## DISTEMERE CCGAMMA BLOCKMATRIX
#CCGammablock  <- CCGammaToBlockmatrix(prec_mes,lag=0,p=2,tolerance=1.e-10)
#corg_mes <- CCGamma ##$nooccurence_gcorrelation

corplot(x=NULL,y=generation1p,corx=corg_mes,origin=origin)

#### 

use <- "pairwise.complete.obs"
method <- "pearson" 


prec_mes_s <- addseason(prec_mes,origin=origin,nodata=TRUE)
prec_mes_sl <- lapply(X=levels(prec_mes_s$season),FUN=function(x,prec_mes_s){ prec_mes_s[prec_mes_s$season==x,names(prec_mes_s)!="season"]},prec_mes_s=prec_mes_s)
names(prec_mes_sl) <- prec_mes_sl
corg_mes_s <- lapply(X=prec_mes_sl,FUN=CCGamma,lag=0,only.matrix=TRUE,tolerance=1.e-10)

###
###
###

corplot(x=NULL,y=generation1p,corx=corg_mes_s,origin=origin,use=use,method=method,season=TRUE)





