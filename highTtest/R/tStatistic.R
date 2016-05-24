#------------------------------------------------------------------------------#
# Two-sample t-statistic                                                       #
#------------------------------------------------------------------------------#
#                                                                              #
#  Inputs :                                                                    #
#                                                                              #
# data1 : matrix of first subset data                                          #
#                                                                              #
# data2 : matrix of second subset data                                         #
#                                                                              #
#  Outputs :                                                                   #
#                                                                              #
#  list object containing t-statistic and p-value for each observation         #
#------------------------------------------------------------------------------#
tStatistic <- function(data1,
                       data2){

  m <- nrow(data1)

  tmp <- sapply(1:m, 
                function(x, data1, data2){
                  tt <- t.test(data2[x,],data1[x,])
                  return(c(tt$statistic,tt$p.value))
                },
                data1 = data1,
                data2 = data2)

  return(list("statistic" = tmp[1,], 
              "pv" = tmp[2,]))
}

