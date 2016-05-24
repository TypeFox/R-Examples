## ----load_gmwm, message=FALSE--------------------------------------------
library("gmwm")

## ----install_imu, message=FALSE------------------------------------------
if(!require("imudata")){
  install_imudata()
  library("imudata")
}

## ----imuData, message=FALSE----------------------------------------------
data(imu6)
head(imu6)

## ----imu_object----------------------------------------------------------
sensor = imu(imu6, 
             gyros = 1:3,
             accels = 4:6, 
             axis = c('X','Y','Z'),
             freq = 100)

## ----imuwvar-------------------------------------------------------------
wv.imu = wvar(sensor)

## ----imuPlot, fig.width=7, fig.height=4.25-------------------------------
plot(wv.imu)

## ----imu_plot_split, fig.width=7, fig.height=4.25------------------------
plot(wv.imu, split = F)

## ----imuModel------------------------------------------------------------
# Define model
TS.mod.imu = 3*AR1()

# Compute GMWM estimator
model.imu = gmwm.imu(TS.mod.imu, data = imu6[,1])

## ----imu_model_summary, fig.width=7, fig.height=4.25---------------------
summary(model.imu)
plot(model.imu)

