
clfs.nostrat <- function(data, maxx = NULL, com.est = TRUE, conf.int = FALSE, conf.int.level = NULL, no.iter = NULL, points = NULL, fig = TRUE)
{

# remove patients who did not achieve the first disease remission:
data <- data[!is.na(data[,1]),] 

# create a matrix and vectors for estimation of the CLFS function:
E <- data[,1:(ncol(data)-2)] # a matrix with times to events in days (columns are times from therapy initiation to the achievement of the first disease remission, to the loss of the first remission, to the achievement of the second remission etc.)
LastContact <- data[,ncol(data)-1] # a vector with times from therapy initiation to death or to the last contact with patients
Exitus <- data[,ncol(data)] # a censoring indicator (1..patient died, 0..patient is censored)

# subtract the time to the achievement of the first disease remission from all columns of E and LastContact:
LastContact <- LastContact - E[,1] # subtract the time to the achievement of the first disease remission from LastContact
E <- E - E[,1] # subtract the time to the achievement of the first disease remission from all columns of E
E <- E[,2:ncol(E)] # remove first column from the matrix E

# data summary:
summary <- data.frame(cbind(c("The number of patients who achieved at least first disease remission","The number of patients with loss of the first disease remission","The number of patients who died after achievement of the first disease remission","Follow-up since achievement of the first disease remission (in months)","The maximum number of achieved disease remissions"),c(nrow(E),sum(!is.na(E[,1])),sum(Exitus),floor(max(LastContact)/(365/12)),ceiling(max(rowSums(!is.na(data[,1:(ncol(data)-2)])))/2) )))
colnames(summary) <- c("","Summary")


# estimate the CLFS function:
clfsd <- clfs.pest(E,LastContact,Exitus,maxx) # the current leukaemia-free survival function estimate

# compute confidence intervals for the CLFS estimates using bootstrap:
if (conf.int) {
  inx=1:(nrow(E)) # a vector of row indexes of the matrix E
  clfsb = array(0,c(no.iter,(maxx+1))) # allocation of a matrix in which rows are CLFS estimates computed in each iteration
  for (i in 1:no.iter) {
    inxs = sample(inx,replace=TRUE) # take a subsample from the vector of row indexes 
    clfspom <- clfs.pest(E[inxs,],LastContact[inxs],Exitus[inxs],maxx) # current leukaemia-free survival function estimate for the subsample
    clfsb[i,]=clfspom$y # save the CLFS estimate as the ith row of the matrix clfsb
    if ((i%%10)==0) {
      print(paste("Computation of confidence intervals for the CLFS - iteration ",i,"/",no.iter,sep="")) # display an information that every 10th iteration has been computed
    }
  }
  CI <- matrix(0,maxx+1,4) # a matrix with CLFS estimates and confidence intervals (CI) in each day
  colnames(CI) <- c('day', 'lower CI', 'CLFS', 'upper CI') # add column names  
  for (i in 1:(maxx+1)) {
    CI[i,]=c(i-1,quantile(clfsb[!is.na(clfsb[,i]),i],(1-conf.int.level)/2),clfsd$y[i],quantile(clfsb[!is.na(clfsb[,i]),i],(1-conf.int.level)/2+conf.int.level))
  }
  # replace values higher than 1 by 1:
  pci=CI[,2:4] # create matrix pciy with CLFS estimates and confidence intervals
  pci[CI[,2:4]>1]=1 # replace values higher than 1 by 1
  CI[,2:4]=pci # replace the original values in CI by modified values from pciy
  # replace values smaller than 0 by 0:
  pci=CI[,2:4] # create matrix pciy with CLFS estimates and confidence intervals
  pci[CI[,2:4]<0]=0 # replace values smaller than 0 by 0
  CI[,2:4]=pci # replace the original values in CI by modified values from pciy
}

# estimate the LFS function:
if (com.est) {
  lfs <- lfs.pest(E[,1],LastContact,Exitus,maxx) # the common leukaemia-free survival estimate
  
  # compute confidence intervals for the LFS estimates using bootstrap:
  if (conf.int) {
    inx=1:(nrow(E)) # a vector of row indexes of the matrix E
    clfsb = array(0,c(no.iter,(maxx+1))) # allocation of a matrix in which rows are LFS estimates computed in each iteration
    for (i in 1:no.iter) {
      inxs = sample(inx,replace=TRUE) # take a sample from the vector of row indexes 
      clfspom <- lfs.pest(E[inxs,1],LastContact[inxs],Exitus[inxs],maxx) # common leukaemia-free survival estimate for the subsample
      clfsb[i,]=clfspom$y # save the LFS estimate as the ith row of the matrix clfsb
      if ((i%%10)==0) {
        print(paste("Computation of confidence intervals for the LFS - iteration ",i,"/",no.iter,sep="")) # display an information that every 10th iteration has been computed
      }
    }
    CI2 <- matrix(0,maxx+1,4) # a matrix with LFS estimates and confidence intervals (CI) in each day
    colnames(CI2) <- c('day', 'lower CI', 'LFS', 'upper CI') # add column names
    for (i in 1:(maxx+1)) {
      CI2[i,]=c(i-1,quantile(clfsb[!is.na(clfsb[,i]),i],(1-conf.int.level)/2),lfs$y[i],quantile(clfsb[!is.na(clfsb[,i]),i],(1-conf.int.level)/2+conf.int.level))
    }
    # replace values higher than 1 by 1:
    pci=CI2[,2:4] # create matrix pciy with LFS estimates and confidence intervals
    pci[CI2[,2:4]>1]=1 # replace values higher than 1 by 1
    CI2[,2:4]=pci # replace the original values in CI2 by modified values from pciy
    # replace values smaller than 0 by 0:
    pci=CI2[,2:4] # create matrix pciy with LFS estimates and confidence intervals
    pci[CI2[,2:4]<0]=0 # replace values smaller than 0 by 0
    CI2[,2:4]=pci # replace the original values in CI2 by modified values from pciy
  }
}


# create a table with point estimates at the defined time points:
if (com.est) {
  if (conf.int) {
    pest <- cbind(points,CI[1+floor(points*(365/12)),2:4],CI2[1+floor(points*(365/12)),2:4]) # create a matrix with the CLFS and LFS estimates accompanied with confidence intervals at the defined time points
    colnames(pest) <- c('Month',paste('CLFS_',((1-conf.int.level)/2)*100,'%',sep=""),'CLFS',paste('CLFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep=""),paste('LFS_',((1-conf.int.level)/2)*100,'%',sep=""),'LFS',paste('LFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest
  } else {
    pest <- cbind(points,clfsd$y[1+floor(points*(365/12))],lfs$y[1+floor(points*(365/12))]) # create a matrix with the CLFS and LFS estimates at the defined time points
    colnames(pest) <- c('Month','CLFS','LFS') # change column names of the matrix pest
  }
} else {
  if (conf.int) {
    pest <- cbind(points,CI[1+floor(points*(365/12)),2:4]) # create a matrix with the CLFS estimates accompanied with confidence intervals at the defined time points
    colnames(pest) <- c('Month',paste('CLFS_',((1-conf.int.level)/2)*100,'%',sep=""),'CLFS',paste('CLFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest
  } else {
    pest <- cbind(points,clfsd$y[1+floor(points*(365/12))]) # create a matrix with the CLFS estimates at the defined time points
    colnames(pest) <- c('Month','CLFS') # change column names of the matrix pest
  }
}
pest[,2:ncol(pest)] <- round(pest[,2:ncol(pest)],digits=3) # rounds the estimates to 3 decimal places
pest <- rbind(c(0,array(1,ncol(pest)-1)),pest) # add estimates at the time point 0


# create a table with point estimates at each day:
if (com.est) {
  if (conf.int) {
    pest.day <- cbind(0:maxx,CI[,2:4],CI2[,2:4]) # create a matrix with the CLFS and LFS estimates accompanied with confidence intervals at each day
    colnames(pest.day) <- c('Day',paste('CLFS_',((1-conf.int.level)/2)*100,'%',sep=""),'CLFS',paste('CLFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep=""),paste('LFS_',((1-conf.int.level)/2)*100,'%',sep=""),'LFS',paste('LFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest.day
  } else {
    pest.day <- cbind(0:maxx,clfsd$y,lfs$y) # create a matrix with the CLFS and LFS estimates at each day
    colnames(pest.day) <- c('Day','CLFS','LFS') # change column names of the matrix pest.day
  }
} else {
  if (conf.int) {
    pest.day <- cbind(0:maxx,CI[,2:4]) # create a matrix with the CLFS estimates accompanied with confidence intervals at each day
    colnames(pest.day) <- c('Day',paste('CLFS_',((1-conf.int.level)/2)*100,'%',sep=""),'CLFS',paste('CLFS_',((1-conf.int.level)/2+conf.int.level)*100,'%',sep="")) # change column names of the matrix pest.day
  } else {
    pest.day <- cbind(0:maxx,clfsd$y) # create a matrix with the CLFS estimates at each day
    colnames(pest.day) <- c('Day','CLFS') # change column names of the matrix pest.day
  }
}
pest.day[,2:ncol(pest.day)] <- round(pest.day[,2:ncol(pest.day)],digits=3) # rounds the estimates to 3 decimal places


# plot the CLFS and LFS estimates and their confidence intervals:
if (fig) {
  x=0:maxx
  yrs <- floor(maxx/365) # a number of years
  plot(0,1,pch='.',cex=0.01,xlab="Years after achievement of the first disease remission",ylab="Probability",axes=FALSE,xlim=c(0,maxx),ylim=c(0,1))
  axis(2,at=seq(0,1,0.2)) # set points in which tick-marks are drawn on the y-axis
  axis(1,at=seq(0,((yrs+1)*365),365),labels=seq(0,(yrs+1),1)) # set points in which tick-marks are drawn on the x-axis
  if (com.est) {
    if (conf.int) {
      lines(clfsd$x,clfsd$y,type="S",lty=1,lwd=2) # plot the CLFS estimate
      lines(lfs$x,lfs$y,type="S",lty=2,lwd=2) # plot the LFS estimate
      lines(x,CI[,2],type="S",lty=1,lwd=1) # plot the lower confidence interval for the CLFS estimate
      lines(x,CI[,4],type="S",lty=1,lwd=1) # plot the upper confidence interval for the CLFS estimate
      lines(x,CI2[,2],type="S",lty=2,lwd=1) # plot the lower confidence interval for the LFS estimate
      lines(x,CI2[,4],type="S",lty=2,lwd=1) # plot the upper confidence interval for the LFS estimate
      legend("bottomright",legend=c("CLFS",paste(conf.int.level*100,"% conf. int.",sep=""),"LFS",paste(conf.int.level*100,"% conf. int.",sep="")),lwd=c(2,1,2,1),lty=c(1,1,2,2),bty="n",cex=0.9)
    } else {
      lines(clfsd$x,clfsd$y,type="S",lty=1,lwd=1) # plot the CLFS estimate
      lines(lfs$x,lfs$y,type="S",lty=2,lwd=1) # plot the LFS estimate
      legend("bottomright",legend=c("CLFS","LFS"),lty=c(1,2),bty="n",cex=0.9)
    }
  } else {
    if (conf.int) {
      lines(clfsd$x,clfsd$y,type="S",lty=1,lwd=2) # plot the CLFS estimate
      lines(x,CI[,2],type="S",lty=1,lwd=1) # plot the lower confidence interval for the CLFS estimate
      lines(x,CI[,4],type="S",lty=1,lwd=1) # plot the upper confidence interval for the CLFS estimate
      legend("bottomright",legend=c("CLFS",paste(conf.int.level*100,"% conf. int.",sep="")),lwd=c(2,1),lty=c(1),bty="n",cex=0.9)
    } else {
      lines(clfsd$x,clfsd$y,type="S",lty=1,lwd=1) # plot the CLFS estimate
      #legend("bottomright",legend=c("CLFS"),lwd=c(1),lty=c(1),bty="n",cex=0.9)
    }
  }
}


# numbers at risk:
no.risk <- cbind(c(0,points),array(0,length(points)+1)) # allocation of no.risk
colnames(no.risk) <- c('Month','CLFS_Nrisk') # set column names of the matrix no.risk
for (j in 1:(length(points)+1)) {
   no.risk[j,2] <- sum(LastContact>=(1+floor(no.risk[j,1]*(365/12)))) # no.risk for the CLFS estimates
}
if (com.est) {
  no.risk <- cbind(no.risk,array(0,length(points)+1)) # add a new column to no.risk
  colnames(no.risk) <- c('Month','CLFS_Nrisk','LFS_Nrisk') # change column names of the matrix no.risk
  E1 <- E[,1]
  E1[is.na(E[,1])] <- LastContact[is.na(E[,1])]
  for (j in 1:(length(points)+1)) {
    no.risk[j,3] <- sum(E1>=(1+floor(no.risk[j,1]*(365/12)))) # no.risk for the LFS estimates
  }
}

rownames(pest) <- NULL
rownames(pest.day) <- NULL

clfs.nostrat <- list(summary=summary,no.risk=no.risk,pest=pest,pest.day=pest.day)

}

