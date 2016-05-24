### R code from vignette source 'VIM-Imputation.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: load_VIM (eval = FALSE)
###################################################
## library(VIM)


###################################################
### code chunk number 2: load_VIM_GUI (eval = FALSE)
###################################################
## vmGUImenu()	


###################################################
### code chunk number 3: imp_data_sets (eval = FALSE)
###################################################
## sleep_kNN  <- kNN( sleep, k = 5)
## tao_kNN <- kNN(tao, k = 5)
## chorizon_kNN  <- kNN( chorizonDL ,variable=c("Ca","As","Bi"), dist_var=c("Hg","K","La","Li"), k=5)


###################################################
### code chunk number 4: Aggr_sleep_imp_part (eval = FALSE)
###################################################
## sleep_kNN_part <- kNN(sleep, variable=c("Sleep","Dream","NonD"), dist_var=colnames(sleep), k=5)
## aggr(sleep_kNN_part, delimiter="_imp", numbers=TRUE, prop=c(TRUE,FALSE))


###################################################
### code chunk number 5: Aggr_sleep_custom (eval = FALSE)
###################################################
## aggr(sleep_kNN_part, delimiter="_imp", sortVars=TRUE, numbers=TRUE, prop=c(FALSE, FALSE), only.miss=FALSE)
## aggr(sleep_kNN_part, delimiter="_imp", sortVars=TRUE, numbers=TRUE, prop=c(FALSE, FALSE), only.miss=TRUE)
## aggr(sleep_kNN_part, delimiter="_imp", combined=TRUE, sortVars=FALSE, sortCombs=TRUE, numbers=TRUE, prop=c(FALSE, FALSE), only.miss=TRUE)


###################################################
### code chunk number 6: Hist_sleep_imp (eval = FALSE)
###################################################
## tao_vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## histMiss(tao_kNN[,tao_vars], delimiter="imp", selection="any")
## 
## sleep_vars <- c("Exp","Dream","NonD","Dream_imp","NonD_imp")
## barMiss(sleep_kNN[,sleep_vars], delimiter="_imp", selection="any")


###################################################
### code chunk number 7: Hist_sleep_custom (eval = FALSE)
###################################################
## vars <- c("Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## histMiss(sleep_kNN[,vars], delimiter="imp", selection="any")
## histMiss(sleep_kNN[,vars], delimiter="imp", selection="all")
## histMiss(sleep_kNN[,vars], delimiter="imp", selection="all", only.miss=FALSE)


###################################################
### code chunk number 8: Spine_imp (eval = FALSE)
###################################################
## tao_vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## spineMiss(tao_kNN[,tao_vars], delimiter="imp", selection="any")
## 
## sleep_vars <- c("Exp","Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## spineMiss(sleep_kNN[,sleep_vars], delimiter="_imp")


###################################################
### code chunk number 9: Spine_imp_custom (eval = FALSE)
###################################################
## vars <- c("Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## spineMiss(sleep_kNN[,vars], delimiter="imp", selection="all", only.miss=FALSE)


###################################################
### code chunk number 10: Box_sleep_imp (eval = FALSE)
###################################################
## vars <- c("Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## pbox(sleep_kNN[,vars], delimiter="_imp", selection="any")


###################################################
### code chunk number 11: Box_sleep_custom (eval = FALSE)
###################################################
## vars <- c("Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## pbox(sleep_kNN[,vars], delimiter="_imp", selection="all")


###################################################
### code chunk number 12: Box_sleep_imp (eval = FALSE)
###################################################
## vars <- c("Dream","NonD","Sleep","Dream_imp","NonD_imp","Sleep_imp")
## pbox(sleep_kNN[,vars], delimiter="_imp", selection="none")


###################################################
### code chunk number 13: Marginplot_tao_imp (eval = FALSE)
###################################################
## vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## marginplot(tao_kNN[,vars], delimiter="imp", alpha=0.6)


###################################################
### code chunk number 14: Marginplot_imp_methods (eval = FALSE)
###################################################
## tao_kNN <- kNN(tao, k = 5)
## tao_mean <- as.data.frame(impute(tao, what = "mean"))
## tao_mean <- cbind(tao_mean, tao_kNN[,9:11])
## 
## vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## marginplot(tao_kNN[,vars], delimiter="imp", alpha=0.6, main="kNN")
## marginplot(tao_mean[,vars], delimiter="imp", alpha=0.6, main="mean")


###################################################
### code chunk number 15: Scatterplot_tao_imp (eval = FALSE)
###################################################
## vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## scattMiss(tao_kNN[,vars], delimiter="imp", alpha=0.6)


###################################################
### code chunk number 16: Jitt_tao_imp (eval = FALSE)
###################################################
## vars <- c("Air.Temp","Humidity","Air.Temp_imp","Humidity_imp")
## scattJitt(tao_kNN[,vars], delimiter="imp", alpha=0.6)


###################################################
### code chunk number 17: Marginmatrix_tao_imp (eval = FALSE)
###################################################
## vars <- c("Air.Temp","Humidity","Sea.Surface.Temp","Air.Temp_imp","Humidity_imp","Sea.Surface.Temp_imp")
## marginmatrix(tao_kNN[,vars], delimiter = "_imp", alpha=0.6)


###################################################
### code chunk number 18: Scattmatrix_tao_imp (eval = FALSE)
###################################################
## vars <- c("Air.Temp","Humidity","Sea.Surface.Temp","Air.Temp_imp","Humidity_imp","Sea.Surface.Temp_imp")
## scattmatrixMiss(tao_kNN[,vars], highlight="Air.Temp", delimiter = "_imp", alpha=0.6)


###################################################
### code chunk number 19: Scattmatrix_tao_highlight (eval = FALSE)
###################################################
## vars <- c("Air.Temp","UWind","VWind","Air.Temp_imp")
## scattmatrixMiss(tao_kNN[,vars], highlight="Air.Temp", plotvars=c("UWind","VWind"), delimiter = "_imp", alpha=0.6)


###################################################
### code chunk number 20: Parcoord_chorizon_imp (eval = FALSE)
###################################################
## chorizon_kNN <- kNN(chorizonDL[,c(15,101:110)], k=5)
## parcoordMiss(chorizon_kNN, delimiter = "_imp" , plotvars=2:11)


###################################################
### code chunk number 21: Matrixplot_tao_imp (eval = FALSE)
###################################################
## matrixplot(sleep_kNN, delimiter="_imp", sortby="Span")


###################################################
### code chunk number 22: Mosaic_sleep_imp (eval = FALSE)
###################################################
## mosaicMiss(sleep_kNN, highlight=3, plotvars=8:10, delimiter="_imp", miss.labels=FALSE)


###################################################
### code chunk number 23: Mapmiss_chorizon_imp (eval = FALSE)
###################################################
## vars <- c("As","Bi","Ca","As_imp","Bi_imp")
## coo <- chorizon_kNN[, c("XCOO", "YCOO")]
## mapMiss(chorizon_kNN[,vars], coo, map=kola.background, delimiter="_imp", alpha=0.6)


###################################################
### code chunk number 24: Growdot_chorizon_imp (eval = FALSE)
###################################################
## vars <- c("Ca","As","Bi","As_imp","Bi_imp")
## coo <- chorizon_kNN[, c("XCOO", "YCOO")]
## growdotMiss(chorizon_kNN[,vars], coo, kola.background, delimiter="_imp", alpha=0.6)


