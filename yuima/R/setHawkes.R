# setHawkes <- function(yuima, counting.var="N", gFun, Kernel,
#                       var.dx, var.dt = "s", lambda.var = "lambda",
#                       lower.var="0", upper.var = "t",
#                       nrow =1 ,ncol=1){
#
#   g <- yuima:::setMaps(func = gFun, yuima = yuima,
#     nrow = nrow, ncol = ncol)
#
#   yuimadum <- yuima
#   yuimadum@time.variable <- var.dt
#
#   HawkesType <- FALSE
#   if(counting.var %in% var.dx){
#     HawkesType <- TRUE
#   }
#   if(!HawkesType){
#   Integral <- yuima:::setIntegral(yuima=yuimadum,
#     integrand = Kernel, var.dx = var.dx,
#     lower.var = lower.var, upper.var = upper.var,
#     out.var = "", nrow = nrow, ncol = ncol)
#   }else{
#     Integral <- yuima:::setIntegral(yuima=yuimadum,
#       integrand = Kernel, var.dx = var.dx,
#       lower.var = lower.var, upper.var = upper.var,
#       out.var = "", nrow = nrow, ncol = ncol, type ="")
#   }
#   if(g@Output@dimension[1]!=Integral$dimIntegrand[1]){
#     yuima.stop("dimension gFun and kernel mismatch")
#   }
#
#
#   allparam <-unique(c(yuima@parameter@all, g@param@allparamMap,
#                Integral$param$IntegrandParam))
#   common <- unique(c(g@param@common, Integral$param$common))
#   paramHawkes <- list(allparam = allparam, common = common,
#                     gFun = g@param@allparamMap,
#                     Kern = Integral$param$IntegrandParam)
#
# #   IntPpr<- yuima:::setIntegral(yuima=yuimadum,
# #     integrand = Kernel, var.dx = "N",
# #     lower.var = lower.var, upper.var = upper.var,
# #     out.var = "", nrow = nrow, ncol = ncol)
#
#   return(list(Count.Proc = counting.var,
#     gFun = list(param=g@param, output=g@Output),
#     Kernel = Integral, paramHawkes = paramHawkes,
#     model = yuima, SelfEx = HawkesType))
#
# }
