
"lines.frailtyPenal"<-
function (x, type.plot="hazard", conf.bands=TRUE, ...) 
{

    plot.type <- charmatch(type.plot, c("hazard", "survival"), 
        nomatch = 0)
    if (plot.type == 0) {
        stop("estimator must be hazard or survival")
       }

    if(plot.type==1)
       {
         if(conf.bands)
           {
             matlines(x$x1,x$lam,...) 
           }
         else
             lines(x$x1,x$lam[,1],...)   
       }   


    if(plot.type==2)
       {
         if(conf.bands)
           {
             matlines(x$x1,x$surv,...) 
           }
         else
             lines(x$x1,x$surv[,1],...)   
       }
  
  return(invisible())

}



