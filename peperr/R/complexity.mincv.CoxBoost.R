`complexity.mincv.CoxBoost` <- function(response, x, full.data, ...){
   require(CoxBoost)
   time <- response[,"time"]
   status <- response[,"status"]

   cv.res <- cv.CoxBoost(time=time, status=status, x=x, ...)

   min.cv <- cv.res$optimal.step
   min.cv
}

