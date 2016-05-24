
cci.nostrat <- function(data, maxx = NULL, com.est = TRUE, conf.int = FALSE, conf.int.level = NULL, no.iter = NULL, points = NULL, fig = TRUE)
{

# data summary:
summary <- data.frame(cbind(c("The total number of patients","The number of patients who achieved at least first disease remission","The number of patients who died in follow-up","Follow-up (in months)","The maximum number of achieved disease remissions"),c(nrow(data),sum(!is.na(data[,1])),sum(data[,ncol(data)]),floor(max(data[,ncol(data)-1])/(365/12)),ceiling(max(rowSums(!is.na(data[,1:(ncol(data)-2)])))/2) )))
colnames(summary) <- c("","Summary")

# create a matrix and vectors for estimation of the CCI function:
E <- data[,1:(ncol(data)-2)] # a matrix with times to events in days (columns are times from therapy initiation to the achievement of the first disease remission, to the loss of the first remission, to the achievement of the second remission etc.)
LastContact <- data[,ncol(data)-1] # a vector with times from therapy initiation to death or to the last contact with patients
Exitus <- data[,ncol(data)] # a censoring indicator (1..patient died, 0..patient is censored)

# estimate the CCI function:
ccid <- cci.pest(E,LastContact,Exitus,maxx) # the current cumulative incidence function estimate

# compute confidence intervals for the CCI estimates using bootstrap:
if (conf.int) {
  inx=1:(nrow(E)) # a vector of row indexes of the matrix E
  ccib = array(0,c(no.iter,(maxx+1))) # allocation of a matrix in which rows are CCI estimates computed in each iteration
  for (i in 1:no.iter) {
    inxs = sample(inx,replace=TRUE) # take a subsample from the vector of row indexes 
    ccipom <- cci.pest(E[inxs,],LastContact[inxs],Exitus[inxs],maxx) # current cumulative incidence function estimate for the subsample
    ccib[i,]=ccipom$y # save the CCI estimate as the ith row of the matrix ccib
    if ((i%%10)==0) {
      print(paste("Computation of confidence intervals for the CCI - iteration ",i,"/",no.iter,sep="")) # display an information that every 10th iteration has been computed
    }
  }
  CI <- matrix(0,maxx+1,4) # a matrix with CCI estimates and confidence intervals (CI) in each day
  colnames(CI) <- c('day', 'lower CI', 'CCI', 'upper CI') # add column names
  for (i in 1:(maxx+1)) {
    CI[i,]=c(i-1,quantile(ccib[!is.na(ccib[,i]),i],(1-conf.int.level)/2),ccid$y[i],quantile(ccib[!is.na(ccib[,i]),i],(1-conf.int.level)/2+conf.int.level))
  }
  # replace values higher than 1 by 1:
  pci=CI[,2:4] # create matrix pciy with CCI estimates and confidence intervals
  pci[CI[,2:4]>1]=1 # replace values higher than 1 by 1
  CI[,2:4]=pci # replace the original values in CI by modified values from pciy
  # replace values smaller than 0 by 0:
  pci=CI[,2:4] # create matrix pciy with CCI estimates and confidence intervals
  pci[CI[,2:4]<0]=0 # replace values smaller than 0 by 0
  CI[,2:4]=pci # replace the original values in CI by modified values from pciy
}

# estimate the comCI function:
if (com.est) {
  comci <- comci.pest(E[,1],LastContact,Exitus,maxx) # the common cumulative incidence function estimate
  
  # compute confidence intervals for the comCI estimates using bootstrap:
  if (conf.int) {
    inx=1:(nrow(E)) # a vector of row indexes of the matrix E
    ccib = array(0,c(no.iter,(maxx+1))) # allocation of a matrix in which rows are comCI estimates computed in each iteration
    for (i in 1:no.iter) {
      inxs = sample(inx,replace=TRUE) # take a subsample from the vector of row indexes 
      ccipom <- comci.pest(E[inxs,1],LastContact[inxs],Exitus[inxs],maxx) # common cumulative incidence function estimate for the subsample
      ccib[i,]=ccipom$y # save the comCI estimate as the ith row of the matrix ccib
      if ((i%%10)==0) {
        print(paste("Computation of confidence intervals for the comCCI - iteration ",i,"/",no.iter,sep="")) # display an information that every 10th iteration has been computed
      }
    }
    CI2 <- matrix(0,maxx+1,4) # a matrix with comCI estimates and confidence intervals (CI) in each day
    colnames(CI2) <- c('day', 'lower CI', 'comCI', 'upper CI') # add column names 
    for (i in 1:(maxx+1)) {
      CI2[i,]=c(i-1,quantile(ccib[!is.na(ccib[,i]),i],(1-conf.int.level)/2),comci$y[i],quantile(ccib[!is.na(ccib[,i]),i],(1-conf.int.level)/2+conf.int.level))
    }
    # replace values higher than 1 by 1:
    pci=CI2[,2:4] # create matrix pciy with comCI estimates and confidence intervals
    pci[CI2[,2:4]>1]=1 # replace values higher than 1 by 1
    CI2[,2:4]=pci # replace the original values in CI2 by modified values from pciy
    # replace values smaller than 0 by 0:
    pci=CI2[,2:4] # create matrix pciy with comCI estimates and confidence intervals
    pci[CI2[,2:4]<0]=0 # replace values smaller than 0 by 0
    CI2[,2:4]=pci # replace the original values in CI2 by modified values from pciy
  }
}


# create a table with point estimates at the defined time points:
if (com.est) {
  if (conf.int) {
    pest <- cbind(points,CI[1+floor(points*(365/12)),2:4],CI2[1+floor(points*(365/12)),2:4]) # create a matrix with the CCI and comCI estimates accompanied with confidence intervals at the defined time points
    colnames(pest) <- c('Month',paste('CCI_',((1-conf.int.level)/2)*100,'%',sep=""),'CCI',paste('CCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep=""),paste('comCI_',((1-conf.int.level)/2)*100,'%',sep=""),'comCI',paste('comCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest
  } else {
    pest <- cbind(points,ccid$y[1+floor(points*(365/12))],comci$y[1+floor(points*(365/12))]) # create a matrix with the CCI and comCI estimates at the defined time points
    colnames(pest) <- c('Month','CCI','comCI') # change column names of the matrix pest
  }
} else {
  if (conf.int) {
    pest <- cbind(points,CI[1+floor(points*(365/12)),2:4]) # create a matrix with the CCI estimates accompanied with confidence intervals at the defined time points
    colnames(pest) <- c('Month',paste('CCI_',((1-conf.int.level)/2)*100,'%',sep=""),'CCI',paste('CCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest
  } else {
    pest <- cbind(points,ccid$y[1+floor(points*(365/12))]) # create a matrix with the CCI estimates at the defined time points
    colnames(pest) <- c('Month','CCI') # change column names of the matrix pest
  }
}
pest[,2:ncol(pest)] <- round(pest[,2:ncol(pest)],digits=3) # rounds the estimates to 3 decimal places
pest <- rbind(array(0,ncol(pest)),pest) # add estimates at the time point 0


# create a table with point estimates at each day:
if (com.est) {
  if (conf.int) {
    pest.day <- cbind(0:maxx,CI[,2:4],CI2[,2:4]) # create a matrix with the CCI and comCI estimates accompanied with confidence intervals at each day
    colnames(pest.day) <- c('Day',paste('CCI_',((1-conf.int.level)/2)*100,'%',sep=""),'CCI',paste('CCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep=""),paste('comCI_',((1-conf.int.level)/2)*100,'%',sep=""),'comCI',paste('comCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest.day
  } else {
    pest.day <- cbind(0:maxx,ccid$y,comci$y) # create a matrix with the CCI and comCI estimates at each day
    colnames(pest.day) <- c('Day','CCI','comCI') # change column names of the matrix pest.day
  }
} else {
  if (conf.int) {
    pest.day <- cbind(0:maxx,CI[,2:4]) # create a matrix with the CCI estimates accompanied with confidence intervals in each day
    colnames(pest.day) <- c('Day',paste('CCI_',((1-conf.int.level)/2)*100,'%',sep=""),'CCI',paste('CCI_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest.day
  } else {
    pest.day <- cbind(0:maxx,ccid$y) # create a matrix with the CCI estimates at each day
    colnames(pest.day) <- c('Day','CCI') # change column names of the matrix pest.day
  }
}
pest.day[,2:ncol(pest.day)] <- round(pest.day[,2:ncol(pest.day)],digits=3) # rounds the estimates to 3 decimal places


# plot the CCI and comCI estimates and their confidence intervals:
if (fig) {
  x=0:maxx
  yrs <- floor(maxx/365) # a number of years
  plot(0,0,pch='.',cex=0.01,xlab="Years after therapy initiation",ylab="Probability",axes=FALSE,xlim=c(0,maxx),ylim=c(0,1))
  axis(2,at=seq(0,1,0.2)) # set points in which tick-marks are drawn on the y-axis
  axis(1,at=seq(0,((yrs+1)*365),365),labels=seq(0,(yrs+1),1)) # set points at which tick-marks are drawn on the x-axis
  if (com.est) {
    if (conf.int) {
      lines(ccid$x,ccid$y,type="S",lty=1,lwd=2) # plot the CCI estimate
      lines(comci$x,comci$y,type="S",lty=2,lwd=2) # plot the comCI estimate
      lines(x,CI[,2],type="S",lty=1,lwd=1) # plot the lower confidence interval for the CCI estimate
      lines(x,CI[,4],type="S",lty=1,lwd=1) # plot the upper confidence interval for the CCI estimate
      lines(x,CI2[,2],type="S",lty=2,lwd=1) # plot the lower confidence interval for the comCI estimate
      lines(x,CI2[,4],type="S",lty=2,lwd=1) # plot the upper confidence interval for the comCI estimate
      legend("bottomright",legend=c("CCI",paste(conf.int.level*100,"% conf. int.",sep=""),"comCI",paste(conf.int.level*100,"% conf. int.",sep="")),lwd=c(2,1,2,1),lty=c(1,1,2,2),bty="n",cex=0.9)
    } else {
      lines(ccid$x,ccid$y,type="S",lty=1,lwd=1) # plot the CCI estimate
      lines(comci$x,comci$y,type="S",lty=2,lwd=1) # plot the comCI estimate
      legend("bottomright",legend=c("CCI","comCI"),lty=c(1,2),bty="n",cex=0.9)
    }
  } else {
    if (conf.int) {
      lines(ccid$x,ccid$y,type="S",lty=1,lwd=2) # plot the CCI estimate
      lines(x,CI[,2],type="S",lty=1,lwd=1) # plot the lower confidence interval for the CCI estimate
      lines(x,CI[,4],type="S",lty=1,lwd=1) # plot the upper confidence interval for the CCI estimate
      legend("bottomright",legend=c("CCI",paste(conf.int.level*100,"% conf. int.",sep="")),lwd=c(2,1),lty=c(1),bty="n",cex=0.9)
    } else {
      lines(ccid$x,ccid$y,type="S",lty=1,lwd=1) # plot the CCI estimate
      #legend("bottomright",legend=c("CCI"),lwd=c(1),lty=c(1),bty="n",cex=0.9)
    }
  }
}


# numbers at risk:
no.risk <- cbind(c(0,points),array(0,length(points)+1)) # allocation of no.risk
colnames(no.risk) <- c('Month','CCI_Nrisk') # set column names of the matrix no.risk
for (j in 1:(length(points)+1)) {
   no.risk[j,2] <- sum(LastContact>=(1+floor(no.risk[j,1]*(365/12)))) # no.risk for the CCI estimates
}
if (com.est) {
  no.risk <- cbind(no.risk,array(0,length(points)+1)) # add a new column to no.risk
  colnames(no.risk) <- c('Month','CCI_Nrisk','comCI_Nrisk') # change column names of the matrix no.risk
  E1 <- E[,1]
  E1[is.na(E[,1])] <- LastContact[is.na(E[,1])]
  for (j in 1:(length(points)+1)) {
    no.risk[j,3] <- sum(E1>=(1+floor(no.risk[j,1]*(365/12)))) # no.risk for the comCI estimates
  }
}

rownames(pest) <- NULL
rownames(pest.day) <- NULL

cci.nostrat <- list(summary=summary,no.risk=no.risk,pest=pest,pest.day=pest.day)

}
