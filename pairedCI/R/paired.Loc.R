paired.Loc <-
function(x,y, method="parametric", exact=FALSE, conf.level=0.95, alternative="two.sided" )
{

    if (method != "parametric" & method != "nonparametric")
        stop(print("No such method; valid options are 'parametric' and 'nonparametric'"))
    


    mean.x <- mean(x)
    mean.y <- mean(y)
    samsize <- length (x)
    xname <- substitute(x)
    yname <- substitute(y)
    if(alternative=="two.sided"){conf.l <- 1-(1-conf.level)/2}
    else if (alternative=="less" | alternative =="greater"){conf.l <- conf.level}
    else {stop (print("No such alternative; valid options are 'two.sided', 'less' and 'greater'"))}

    switch(method, "parametric" =
       {
           test      <-  (samsize*(samsize-1)*mean.y^2)/sum((y-mean.y)^2)
           if (test < (qt(0.95, samsize-1))^2) stop(print("Side condition not fulfilled"))
           else
           {
               ratio     <-  mean.x/mean.y

               A <- samsize*(samsize-1)*(mean.y^2)-(qt(conf.l,samsize-1)^2)*sum((y-mean.y)^2)
               B <- -2*(samsize*(samsize-1)*mean.x*mean.y-(qt(conf.l,samsize-1)^2)*sum((x-mean.x)*(y-mean.y)))
               C <- samsize*(samsize-1)*(mean.x^2)-(qt(conf.l,samsize-1)^2)*sum((x-mean.x)^2)

               lower_CI <- (-B-sqrt(B^2-4*A*C))/(2*A)
               upper_CI <- (-B+sqrt(B^2-4*A*C))/(2*A)
               cint <- c(lower_CI, upper_CI)
               if(alternative=="less"){cint[1] <- 0}
               else if (alternative=="greater"){cint[2] <- Inf}
               attr(cint, "conf.level") <- conf.level
               rval <- list(method="Ogawa's parametric confidence interval for the ratio of two paired means", data.name=paste(xname, yname, sep=" and "), conf.int=cint, estimate=ratio)
               class(rval) <- "htest"

               return(rval)
           }
       },

           "nonparametric"=
       {


           diffs    <- x - y
           L        <- 100000
           rho      <- seq(mean.x/mean.y-20, mean.x/mean.y+20,length = L)
           V        <- rep(NA,L)

           for (i in 1:L)
           {
               z    <- x - rho[i]*y
               absz <- abs(z)
               V[i] <- sum(rank(absz)[z>0])
           }

           mV               <- abs(V - samsize*(samsize+1)/4)
           estimate.nonpara <- median(range(rho[mV== min(mV)]))

           if(exact==TRUE)
           {

               Uq1 <- qsignrank(1-conf.l, samsize, lower.tail = TRUE, log.p = FALSE)
               Uq2 <- qsignrank(conf.l, samsize, lower.tail = TRUE, log.p = FALSE)

               Lalpha <- abs(V - Uq2)
               Ualpha <- abs(V - Uq1)

               lower_CI   <- range(rho[Lalpha== min(Lalpha)])[1]
               upper_CI   <- range(rho[Ualpha== min(Ualpha)])[2]
               cint         <- c(lower_CI,upper_CI)
               if(alternative=="less"){cint[1] <- 0}
               else if (alternative=="greater"){cint[2] <- Inf}
               attr(cint, "conf.level") <- conf.level
               rval <- list(method="Bennett's exact nonparametric confidence interval for the ratio of two paired means", data.name=paste(xname, yname, sep=" and "), conf.int=cint, estimate=estimate.nonpara)
               class(rval) <- "htest"

               return(rval)

           }


           else if (exact==FALSE){
               Wmean   <- 0.25*samsize*(samsize+1)
               Wvar    <- (1/24)*samsize*(samsize+1)*(2*samsize + 1)
               qz      <- qnorm(1-conf.l, lower.tail=FALSE)
               Uaq1    <- Wmean - qz*sqrt(Wvar)
               Uaq2    <- Wmean + qz*sqrt(Wvar)
               Lalpha  <- abs(V - Uaq2)
               Ualpha  <- abs(V - Uaq1)

               lower_CI <- range(rho[Lalpha== min(Lalpha)])[1]
               upper_CI <- range(rho[Ualpha== min(Ualpha)])[2]
               cint         <- c(lower_CI,upper_CI)
               if(alternative=="less"){cint[1] <- 0}
               else if (alternative=="greater"){cint[2] <- Inf}
               attr(cint, "conf.level") <- conf.level
               rval <- list(method="Bennett's asymptotic nonparametric confidence interval for the ratio of two paired means", data.name=paste(xname, yname, sep=" and "), conf.int=cint, estimate=estimate.nonpara)
               class(rval) <- "htest"

               return(rval)


           }
       }
           )

}
