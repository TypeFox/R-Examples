CatDynPred <-
function(x, method)
   {
     fleet.name <- x$Data$Properties$Fleets$Fleet;
     p          <- x$Model[[method]]$Type;
     options(warn=-1)
     if(class(x) != "catdyn")
       {stop("Pass an object of class 'catdyn' to CatDynPred")}
     if(class(method) != "character")
       {stop("method must be a character corresponding to one of the numerical methods passed to the CatDyn function")}  
     if(length(method) != 1)
       {stop("Provide the name of just one the numerical methods used to fit the model")}  
     if(sum(method == names(x$Model)) == 0)
       {stop("The method provided in 'method' was not used to numerically fit the model")}
     if(is.na(x$Model[[method]]$AIC))
       {stop("The selected method failed. Consider trying a different method")}
     if(length(fleet.name) == 1)
       {
        parlist <- list(par=log(as.numeric(x$Model[[method]]$bt.par)), 
                        dates=x$Model[[method]]$Dates, 
                        obscat1=x$Data$Data[[fleet.name]][,5], 
                        obseff1=x$Data$Data[[fleet.name]][,2], 
                        obsmbm1=x$Data$Data[[fleet.name]][,4], 
                        distr=x$Model[[method]]$Distr, 
                        properties=x$Data$Properties,
                        output="predict")
        if(p == 0)
          {
          results <- do.call(.CDMN0P, parlist);
          }
        else if(p == 1)
          {
          results <- do.call(.CDMN1P, parlist);
          }
        else if(p == -1)
          {
          results <- do.call(.CDMNT1P, parlist);
          }
        else if(p == 2)
          {
          results <- do.call(.CDMN2P, parlist);
          }
        else if(p == -2)
          {
          results <- do.call(.CDMNT2P, parlist);
          }
        else if(p == 3)
          {
          results <- do.call(.CDMN3P, parlist);
          }
        else if(p == -3)
          {
          results <- do.call(.CDMNT3P, parlist);
          }
        else if(p == 4)
          {
          results <- do.call(.CDMN4P, parlist);
          }
        else if(p == -4)
          {
          results <- do.call(.CDMNT4P, parlist);
          }
        else if(p == 5)
          {
          results <- do.call(.CDMN5P, parlist);
          }
        else if(p == -5)
          {
          results <- do.call(.CDMNT5P, parlist);
          }
        else if(p == 6)
          {
          results <- do.call(.CDMN6P, parlist);
          }
        else if(p == -6)
          {
          results <- do.call(.CDMNT6P, parlist);
          }
        else if(p == 7)
          {
          results <- do.call(.CDMN7P, parlist);
          }
        else if(p == -7)
          {
          results <- do.call(.CDMNT7P, parlist);
          }
        else if(p == 8)
          {
          results <- do.call(.CDMN8P, parlist);
          }
        else if(p == -8)
          {
          results <- do.call(.CDMNT8P, parlist);
          }
        else if(p == 9)
          {
          results <- do.call(.CDMN9P, parlist);
          }
        else if(p == -9)
          {
          results <- do.call(.CDMNT9P, parlist);
          }
        else if(p == 10)
          {
          results <- do.call(.CDMN10P, parlist);
          }
        else if(p == -10)
          {
          results <- do.call(.CDMNT10P, parlist);
          }
        else if(p == 11)
          {
          results <- do.call(.CDMN11P, parlist);
          }
        else if(p == -11)
          {
          results <- do.call(.CDMNT11P, parlist);
          }
        else if(p == 12)
          {
          results <- do.call(.CDMN12P, parlist);
          }
        else if(p == -12)
          {
          results <- do.call(.CDMNT12P, parlist);
          }
        else if(p == 13)
          {
          results <- do.call(.CDMN13P, parlist);
          }
        else if(p == -13)
          {
          results <- do.call(.CDMNT13P, parlist);
          }
        else if(p == 14)
          {
          results <- do.call(.CDMN14P, parlist);
          }
        else if(p == -14)
          {
          results <- do.call(.CDMNT14P, parlist);
          }
        else if(p == 15)
          {
          results <- do.call(.CDMN15P, parlist);
          }
        else if(p == -15)
          {
          results <- do.call(.CDMNT15P, parlist);
          }
        else if(p == 16)
          {
          results <- do.call(.CDMN16P, parlist);
          }
        else if(p == -16)
          {
          results <- do.call(.CDMNT16P, parlist);
          }
        else if(p == 17)
          {
          results <- do.call(.CDMN17P, parlist);
          }
        else if(p == -17)
          {
          results <- do.call(.CDMNT17P, parlist);
          }
        else if(p == 18)
          {
          results <- do.call(.CDMN18P, parlist);
          }
        else if(p == -18)
          {
          results <- do.call(.CDMNT18P, parlist);
          }
        else if(p == 19)
          {
          results <- do.call(.CDMN19P, parlist);
          }
        else if(p == -19)
          {
          results <- do.call(.CDMNT19P, parlist);
          }
        else if(p == 20)
          {
          results <- do.call(.CDMN20P, parlist);
          }
        else if(p == -20)
          {
          results <- do.call(.CDMNT20P, parlist);
          }
       }
     else if(length(fleet.name) == 2) 
       {
        parlist <- list(par=log(as.numeric(x$Model[[method]]$bt.par)), 
                        dates=x$Model[[method]]$Dates, 
                        obscat1=x$Data$Data[[fleet.name[1]]][,5], 
                        obseff1=x$Data$Data[[fleet.name[1]]][,2], 
                        obsmbm1=x$Data$Data[[fleet.name[1]]][,4],
                        obscat2=x$Data$Data[[fleet.name[2]]][,5], 
                        obseff2=x$Data$Data[[fleet.name[2]]][,2], 
                        obsmbm2=x$Data$Data[[fleet.name[2]]][,4], 
                        distr=x$Model[[method]]$Distr, 
                        properties=x$Data$Properties,
                        output="predict")
        if(sum(p==c(0,0)) == length(p))
          {
          results <- do.call(.CDMN0P0P, parlist);
          }
        else if(sum(p==c(0,1)) == length(p))
          {
          results <- do.call(.CDMN0P1P, parlist);
          }
        else if(sum(p==c(0,2)) == length(p))
          {
          results <- do.call(.CDMN0P2P, parlist);
          }
        else if(sum(p==c(0,3)) == length(p))
          {
          results <- do.call(.CDMN0P3P, parlist);
          }
        else if(sum(p==c(0,4)) == length(p))
          {
          results <- do.call(.CDMN0P4P, parlist);
          }
        else if(sum(p==c(0,5)) == length(p))
          {
          results <- do.call(.CDMN0P5P, parlist);
          }
        else if(sum(p==c(1,1)) == length(p))
          {
          results <- do.call(.CDMN1P1P, parlist);
          }                                                                                                                                            
        else if(sum(p==c(1,2)) == length(p))
          {
          results <- do.call(.CDMN1P2P, parlist);
          }
        else if(sum(p==c(1,3)) == length(p))
          {
          results <- do.call(.CDMN1P3P, parlist);
          }
        else if(sum(p==c(1,4)) == length(p))
          {
          results <- do.call(.CDMN1P4P, parlist);
          }
        else if(sum(p==c(1,5)) == length(p))
          {
          results <- do.call(.CDMN1P5P, parlist);
          }
        else if(sum(p==c(2,2)) == length(p))
          {
          results <- do.call(.CDMN2P2P, parlist);
          }
        else if(sum(p==c(2,3)) == length(p))
          {
          results <- do.call(.CDMN2P3P, parlist);
          }
        else if(sum(p==c(2,4)) == length(p))
          {
          results <- do.call(.CDMN2P4P, parlist);
          }
        else if(sum(p==c(2,5)) == length(p))
          {
          results <- do.call(.CDMN2P5P, parlist);
          }
        else if(sum(p==c(3,3)) == length(p))
          {
          results <- do.call(.CDMN3P3P, parlist);
          }
        else if(sum(p==c(3,4)) == length(p))
          {
          results <- do.call(.CDMN3P4P, parlist);
          }
        else if(sum(p==c(3,5)) == length(p))
          {
          results <- do.call(.CDMN3P5P, parlist);
          }
        else if(sum(p==c(4,4)) == length(p))
          {
          results <- do.call(.CDMN4P4P, parlist);
          }
        else if(sum(p==c(4,5)) == length(p))
          {
          results <- do.call(.CDMN4P5P, parlist);
          }
        else if(sum(p==c(5,5)) == length(p))
          {
          results <- do.call(.CDMN5P5P, parlist);
          }
        else if(sum(p==c(6,6)) == length(p))
          {
          results <- do.call(.CDMN6P6P, parlist);
          }
        else if(sum(p==c(7,7)) == length(p))
          {
          results <- do.call(.CDMN7P7P, parlist);
          }
        else if(sum(p==c(8,8)) == length(p))
          {
          results <- do.call(.CDMN8P8P, parlist);
          }
        else if(sum(p==c(9,9)) == length(p))
          {
          results <- do.call(.CDMN9P9P, parlist);
          }
        else if(sum(p==c(10,10)) == length(p))
          {
          results <- do.call(.CDMN10P10P, parlist);
          }
        else if(sum(p==c(11,11)) == length(p))
          {
          results <- do.call(.CDMN11P11P, parlist);
          }
        else if(sum(p==c(12,12)) == length(p))
          {
          results <- do.call(.CDMN12P12P, parlist);
          }
        else if(sum(p==c(13,13)) == length(p))
          {
          results <- do.call(.CDMN13P13P, parlist);
          }
        else if(sum(p==c(14,14)) == length(p))
          {
          results <- do.call(.CDMN14P14P, parlist);
          }
        else if(sum(p==c(15,15)) == length(p))
          {
          results <- do.call(.CDMN15P15P, parlist);
          }
        else if(sum(p==c(16,16)) == length(p))
          {
          results <- do.call(.CDMN16P16P, parlist);
          }
        else if(sum(p==c(17,17)) == length(p))
          {
          results <- do.call(.CDMN17P17P, parlist);
          }
        else if(sum(p==c(18,18)) == length(p))
          {
          results <- do.call(.CDMN18P18P, parlist);
          }
        else if(sum(p==c(19,19)) == length(p))
          {
          results <- do.call(.CDMN19P19P, parlist);
          }
        else if(sum(p==c(20,20)) == length(p))
          {
          results <- do.call(.CDMN20P20P, parlist);
          }
       }
     results$Model$Method     <- method;
     results$Model$AIC        <- x$Model[[method]]$AIC;
     class(results)           <- "CatDynMod";
     return(results);
    }
