m.inf.3a <- function(x,
                     which1,
                     which2,
                     which3,
                     dispersion=c('mm', 's', 'se'))
{
  switch(match.arg(dispersion),
         mm = {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which1]],
                                      group=x$model[[which2]],
                                      group2=x$model[[which3]]),
                              function(x) c(mean=mean(x),
                                            min=min(x),
                                            max=max(x)))[,2:4]
         }, s = {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which1]],
                                      group=x$model[[which2]],
                                      group2=x$model[[which3]]),
                              function(x) c(mean=mean(x),
                                            'm - s'=mean(x) - sd(x),
                                            'm + s'=mean(x) + sd(x)))[,2:4]
         }, se= {
           m.inf <- aggregate(x$model[,1],
                              by=list(x$model[[which1]],
                                      group=x$model[[which2]],
                                      group2=x$model[[which3]]),
                              function(x) c(mean=mean(x),
                                            se.min=mean(x) - (sd(x) / sqrt(length(x))),
                                            se.max=mean(x) + (sd(x) / sqrt(length(x)))))[,2:4]
         })
}
