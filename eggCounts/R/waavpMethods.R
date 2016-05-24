######################################
## compute standard FECR
#######################################

fecrtCI <- function(epg1, epg2, paired=FALSE, alpha=0.05, R=1999){
  if (!is.logical(paired))
    stop("The paired argument must be a logical", call.=FALSE)
  if ((R < 1)|(ceiling(R)!=floor(R)))
    stop("'R' should be a positive integer", call.=FALSE)
  if ((alpha<0)||(alpha>1))
    stop("alpha must be between 0 and 1", call.=FALSE)
  if (R > 9999)      cat("NOTE: 'R' seems high\n")
  
   conf.level <- 1-alpha

   if(paired){
      nas <- is.na(epg1) | is.na(epg2)
      epg1 <- epg1[!nas]
      epg2 <- epg2[!nas]
   } else {
      epg1 <- epg1[!is.na(epg1)]
      epg2 <- epg2[!is.na(epg2)]
   }

   ni <- c(length(epg1), length(epg2))
   data <- data.frame(epgs = c(epg1,epg2), grp = as.factor(rep(1:2, ni)))
   meanratio_unpaired <- function(d, i) {
      ind1 <- i[1:ni[1]]
      ind2 <- i[-(1:ni[1])]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
   }
   meanratio_paired <- function(d,i){
      ind1 <- i[1:ni[1]]
      ind2 <- ind1+ni[1]
      x1 <- d[ind1, 1]
      x2 <- d[ind2, 1]
      1-mean(x2)/mean(x1)
   }

   if(paired){
      meanratio <- meanratio_paired
   } else {
      meanratio <- meanratio_unpaired
   }
   boot.out <- boot(data, meanratio,R=R, stype="i", strata=data[,"grp"])
   ## Ratio of means (percentile bootstrap) ##
   if (sum(epg1)!=0){
   if(sd(boot.out$t[!is.na(boot.out$t) & is.finite(boot.out$t)])==0){
      conf.int <- c(NA,NA)
   } else {
      conf.int <- boot.ci(boot.out = boot.out, conf = conf.level, type = c("perc"))$perc[4:5]
   }} else {conf.int <- c(NA,NA)}

   estimate <- 1-mean(epg2)/mean(epg1)

   # approximate CI according to Coles for unpaired situation
   talpha<- qt(1-alpha/2,sum(ni)-2)
   lowerCI <- ifelse(estimate==100, NA, 100*(1-mean(epg2)/mean(epg1)*exp(talpha*sqrt(var(epg1)/ni[1]/mean(epg1)^2+var(epg2)/ni[2]/mean(epg2)^2))) )
   upperCI <- ifelse(estimate==100, NA, 100*(1-mean(epg2)/mean(epg1)*exp(-talpha*sqrt(var(epg1)/ni[1]/mean(epg1)^2+var(epg2)/ni[2]/mean(epg2)^2))) )

   return(list(estimate = estimate*100, bootCI = conf.int*100, approxCI=c(lowerCI,upperCI)))
}
