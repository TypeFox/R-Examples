paired.Scale <-
function(x,y, method="parametric",  conf.level=0.95, alternative="two.sided" )
{

    if (method != "parametric" & method != "nonparametric")
        stop(print("No such method; valid options are 'parametric' and 'nonparametric'"))
 
    var.x <- var(x)
    var.y <- var(y)
    samsize <- length (x)
    xname <- substitute(x)
    yname <- substitute(y)
    if(alternative=="two.sided"){conf.l <- 1-(1-conf.level)/2}
    else if (alternative=="less" | alternative =="greater"){conf.l <- conf.level}
    else {stop (print("No such alternative; valid options are 'two.sided', 'less' and 'greater'"))}

    switch(method, "parametric" =
       {
           correl <-cor(x,y)
           g_term <- 1 + (2*(1-correl^2)*(qt(conf.l,samsize-2))^2)/(samsize-2)
           ratio<- var.x/var.y
           lower_CI <- ratio * (g_term - sqrt(g_term^2-1 ))
           upper_CI <- ratio * (g_term + sqrt(g_term^2-1 ))

           cint <- c(sqrt(lower_CI), sqrt(upper_CI))
           if(alternative=="less"){cint[1] <- 0}
           else if (alternative=="greater"){cint[2] <- Inf}
           attr(cint, "conf.level") <- conf.level
           rval <- list(method="Bonett's parametric confidence interval for the ratio of two paired standard deviations (inversion of Pitman-Morgan test)", data.name=paste(xname, yname, sep=" and "), conf.int=cint, estimate=sqrt(var.x)/sqrt(var.y))
           class(rval) <- "htest"

           return(rval)
       }
           ,


           "nonparametric"=
       {
           mean.x <- mean(x)
           mean.y <- mean(y)
           median.x <- median(x)
           median.y <- median(y)
           madm.x <- sum(abs(x-median.x))/samsize
           madm.y <- sum(abs(y-median.y))/samsize
           ratio <- madm.x/madm.y
           var.log.madm.x <-(var(x)/madm.x^2+((mean.x-median.x)/madm.x)^2-1)/samsize
           var.log.madm.y <-(var(y)/madm.y^2+((mean.y-median.y)/madm.y)^2-1)/samsize
           d.x <- abs(x - median.x)
           d.y <- abs(y - median.y)
           cor.dxdy <- cor(d.x, d.y)
           var.log.madm.xy <- var.log.madm.x+var.log.madm.y-2*cor.dxdy*sqrt(var.log.madm.x*var.log.madm.y)

           lower_CI      <-exp(log(ratio)-qnorm(conf.l)*sqrt(var.log.madm.xy))
           upper_CI      <-exp(log(ratio)+qnorm(conf.l)*sqrt(var.log.madm.xy))

           cint         <- c(lower_CI,upper_CI)
           if(alternative=="less"){cint[1] <- 0}
           else if (alternative=="greater"){cint[2] <- Inf}
           attr(cint, "conf.level") <- conf.level
           rval <- list(method="Asymptotic nonparametric confidence interval for the ratio of two paired mean absolute deviations of the median (Bonett and Seier)", data.name=paste(xname, yname, sep=" and "), conf.int=cint, estimate=ratio)
           class(rval) <- "htest"

           return(rval)

       }
           )
}
