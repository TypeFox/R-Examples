m.inf.2a <- function(x,
                     which1,
                     which2,
                     dispersion=c('mm', 's', 'se'))
{
  switch(match.arg(dispersion),
         mm = {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which2]],
                                      group=x$model[[which1]]),      
                              function(x) c(mean=mean(x),
                                            min=min(x),
                                            max=max(x)))[,2:3]
         }, s = {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which2]],
                                      group=x$model[[which1]]),      
                              function(x) c(mean=mean(x),
                                            'm - s'=mean(x) - sd(x),
                                            'm + s'=mean(x) + sd(x)))[,2:3]
         }, se= {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which2]],
                                      group=x$model[[which1]]),      
                              function(x) c(mean=mean(x),
                                            se.min=mean(x) - (sd(x) / sqrt(length(x))),
                                            se.max=mean(x) + (sd(x) / sqrt(length(x)))))[,2:3]
         })
}
