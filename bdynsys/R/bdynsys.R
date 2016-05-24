# original code for bdynsys and the functions it calls in Matlab by Shyam Ranganathan 

# before calling bdynsys, load data (loadData.R) and create panel data structure called "datap" (see loadData.R)
# main function to be used in the console for BDYNSYS modeling, calls all other core function
# indnr is the number of indicators to be included in the modeling
# paramnr is the maximum number of polyomial terms to be included in the modeling, recommended is 5-6
# x, y, z, v are the indicators (up to four) that will be included in the modeling from the data 
# define argument rangeX, rangeY, rangeZ, rangeV as following example shows: rangeX = seq(0, 1, by = 0.01)
# visualisation functions are separate from the modeling procedure and can be called on their own

# example of a call of dsm from console: bdynsys(datap, 2, 6, datap$logGDP, datap$EmanzV)

bdynsys <- function(dataset, indnr, paramnr, x, y, z, v)
{
  if (indnr == 2)
  {
    procdata <- preprocess_data(indnr, x, y)
    results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
                      procdata$mx, procdata$my)
    bayesfactor <- bayesfac(indnr, paramnr, results$SEtestx, results$SEtesty, procdata$xs, 
                            procdata$ys, procdata$chXs, procdata$chYs)
   }
  
  if (indnr == 3)
  {
    procdata <- preprocess_data(indnr, x, y, z)
    results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
                       procdata$mx, procdata$my, procdata$zs, procdata$chZs, procdata$mz)
    bayesfactor <- bayesfac(indnr, paramnr, results$SEtestx, results$SEtesty, procdata$xs, 
                            procdata$ys, procdata$chXs, procdata$chYs, procdata$zs, 
                            procdata$chZs, results$SEtestz)
  }
   
  if (indnr == 4)
  {
    procdata <- preprocess_data(indnr, x, y, z, v)
    results <- dysymod(indnr, paramnr, procdata$xs, procdata$ys, procdata$chXs, procdata$chYs, 
                      procdata$mx, procdata$my, procdata$zs, procdata$chZs, procdata$mz,
                      procdata$vs, procdata$chVs, procdata$mv)
    bayesfactor <- bayesfac(indnr, paramnr, results$SEtestx, results$SEtesty, procdata$xs, 
                            procdata$ys, procdata$chXs, procdata$chYs, procdata$zs, procdata$chZs,
                            results$SEtestz, procdata$vs, procdata$chVs, results$SEtestv)
  }
}