calib_swat_ex <-
function(flowgage,rch=3){
print("An example of the commands necessary to calibrate SWAT for a given reach using flowgage measurements.")
print("This could take 30 minutes to several days depending on how fast your computer is.")
change_params=""
rm(change_params)
load(paste(path.package("EcoHydRology"), "data/change_params.rda", sep = "/"))
#data(change_params)
calib_range=c("1999-12-31","2010-12-31")
params_select=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,19,21,23,24,32,33)
calib_params=change_params[params_select,]
setup_swatcal(calib_params)
x=calib_params$current
outDEoptim<-DEoptim(swat_objective_function_rch,calib_params$min,calib_params$max, DEoptim.control(strategy = 6,NP = 16,itermax=200,parallelType = 1, packages = c("SWATmodel")),calib_range,calib_params,flowgage,rch)
x=outDEoptim$optim$bestmem
swat_objective_function_rch(x,calib_range,calib_params,flowgage,rch,save_results=T)
save(outDEoptim,file="outDEoptim")


}
