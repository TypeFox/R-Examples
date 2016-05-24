as.receiver<-function(x) {
  class(x)<-"receiver"
  return(x)
}
findintersects<-function(x) {
  
  ## Convert compass bearing to standard radian measure (Lenth 1981)
  x$theta = (pi/180*(90-x$Azimuth)) 
  
  ## Calculate additional variables
  x$sin = sin(x$theta)
  x$cos = cos(x$theta)
  x$tan = tan(x$theta)
  
  ## Create data frame to store coordinates of all intersections
  results<-data.frame(GID=NA,X=NA,Y=NA) 
  
  ## Keep track of total number of intersections in the data
  counter=1; 
  for (group in unique(x$GID)) {
    
    ## Create subdata to store data of the single GID
    subdata=data.frame(X=x$Easting[x$GID==group],Y=x$Northing[x$GID==group],TAN=x$tan[x$GID==group])
    
    ## Calculate number of bearings in the single GID
    subdata_length=length(subdata)
    
    ## Cycle through all combinations of bearings in the single GID 
    for (i in 1:subdata_length) {
      for (j in 1:subdata_length) {
        
        ## Do not calcuate intersections of bearings with themselves
        if (i !=j ) {
          
          ## Import GID variable
          results[counter,1]=group
          
          ## Calculate bearing intersect
          results[counter,3]=((subdata[j,'X']-subdata[i,'X'])*subdata[i,'TAN']*subdata[j,'TAN']-subdata[j,'Y']*subdata[i,'TAN']+subdata[i,'Y']*subdata[j,'TAN'])/(subdata[j,'TAN']-subdata[i,'TAN'])
          results[counter,2]=(results[counter,3]-subdata[j,'Y'])/subdata[j,'TAN']+subdata[j,'X']
          
          ## Increment total number of calculated intersects
          counter=counter+1
        }
      }
    }
  }
  
  class(results)<-"intersect"
  return(results)
}
locate<-function(x) {
  
  ## Define function to solve system of MLE equations
  solver<-function(par) {
    
    ## Isolate portion of data corresponding to the unique grouping
    subdata=data.frame(X=x$Easting[x$GID==group],Y=x$Northin[x$GID==group],SIN=x$sin[x$GID==group],COS=x$cos[x$GID==group])
    
    ## Create variables to hold transmitter location estimates
    loc_est = numeric(length(par))
    transmitter_x<-par[1]
    transmitter_y<-par[2]
    
    ## Define system of MLE equations (Lenth 1981)
    loc_est[1] = -sum((transmitter_y-subdata$Y)*(subdata$SIN*(transmitter_x-subdata$X)-subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+(transmitter_y-subdata$Y)^2)^0.5)^3)
    loc_est[2] = sum((transmitter_x-subdata$X)*(subdata$SIN*(transmitter_x-subdata$X)-subdata$COS*(transmitter_y-subdata$Y))/(((transmitter_x-subdata$X)^2+(transmitter_y-subdata$Y)^2)^0.5)^3)
    
    ## Return transmitter location estimate
    loc_est  
  }
  
  ## Convert compass bearing to standard radian measure (Lenth 1981)
  x$theta = (pi/180*(90-x$Azimuth)) 
  
  ## Calculate necessary variables
  x$sin = sin(x$theta)
  x$cos = cos(x$theta)
  x$tan = tan(x$theta)
  
  ## Calculates intersection of all bearings in each grouping
  intersections<-findintersects(x) 
  
  ## Create data frame to store calculated transmitter locations
  transmitter<-data.frame(X=NA,Y=NA,BadPoint=NA,Var_X=NA,Var_Y=NA,Cov_XY=NA,AngleDiff=NA,Date=NA,Time=NA)
  
  ## Declare buffer (in meters) for ensuring location validity  
  buffer<-5
  
  ## Calculate transmitter location for each grouping
  for (group in unique(x$GID)) {
    
    ## Calculate starting point for solver function based on the mean of bearing intersections
    start_point_x<-mean(intersections$X[intersections$GID==group]) 
    start_point_y<-mean(intersections$Y[intersections$GID==group]) 
    
    ## Temporarily store results of solver function
    triangulation_results<-nleqslv(c(start_point_x,start_point_y),solver) 
    
    ## Withdraws transmitter location from list of solver function results
    location<-triangulation_results$x
    
    ## Ensure location validity prior to recording
    valid=TRUE
    if (location[1]<min(intersections$X[intersections$GID==group]-buffer)) valid=FALSE
    if (location[1]>max(intersections$X[intersections$GID==group]+buffer)) valid=FALSE
    if (location[2]<min(intersections$Y[intersections$GID==group]-buffer)) valid=FALSE
    if (location[2]>max(intersections$Y[intersections$GID==group]+buffer)) valid=FALSE
    
    ## Set default color scheme for transmitter ploting
    transmitter[group,3]<-0
    
    ## Transfer calculated locations into results variable
    if (valid) {
      transmitter[group,1]<-location[1]
      transmitter[group,2]<-location[2]
    }
    else {
      warning("Bad point detected in Grouping ",group)
      transmitter[group,1]<-start_point_x
      transmitter[group,2]<-start_point_y
      ## Set color scheme to red for transmitter ploting of bad points
      transmitter[group,3]<-1
    }
    
    ## Create error data frame
    errors=data.frame(X=x$Easting[x$GID==group],Y=x$Northing[x$GID==group],theta=(pi/180*(90-x$Azimuth[x$GID==group])),si=x$sin[x$GID==group],ci=x$cos[x$GID==group])
    errors$X_hat<-transmitter[group,1]
    errors$Y_hat<-transmitter[group,2]
    errors$di<-((errors$X_hat-errors$X)^2+(errors$Y_hat-errors$Y)^2)^(1/2)
    errors$si_star<-(errors$Y_hat-errors$Y)/(errors$di)^3
    errors$ci_star<-(errors$X_hat-errors$X)/(errors$di)^3
    #errors$mew_i<-atan((errors$Y_hat-errors$Y)/(errors$X_hat-errors$X))
    
    ## Calculate azimuth bearings to transmitter location and append to error data frame
    errors$mew_i<-NA
    for (e in 1:nrow(errors)) {
      if ((errors$X[e]==errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {errors$mew_i[e]<-0}
      if ((errors$X[e]==errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {errors$mew_i[e]<-180}
      if ((errors$Y[e]==errors$Y_hat[e])&&(errors$X[e]<errors$X_hat[e])) {errors$mew_i[e]<-90}
      if ((errors$Y[e]==errors$Y_hat[e])&&(errors$X[e]>errors$X_hat[e])) {errors$mew_i[e]<-270}
      if ((errors$X[e]<errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$X_hat[e]-errors$X[e])/abs(errors$Y_hat[e]-errors$Y[e]))
      }
      if ((errors$X[e]<errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$Y_hat[e]-errors$Y[e])/abs(errors$X_hat[e]-errors$X[e]))+90
      }
      if ((errors$X[e]>errors$X_hat[e])&&(errors$Y[e]>errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$X_hat[e]-errors$X[e])/abs(errors$Y_hat[e]-errors$Y[e]))+180
      }
      if ((errors$X[e]>errors$X_hat[e])&&(errors$Y[e]<errors$Y_hat[e])) {
        errors$mew_i[e]<-180/pi*atan(abs(errors$Y_hat[e]-errors$Y[e])/abs(errors$X_hat[e]-errors$X[e]))+270
      }
    }
    
    ## Convert calculated bearings to standard radian measure (Lenth 1981)
    errors$mew_i = (pi/180*(90-errors$mew_i))  
    
    ## Calculate covariance entries as per Lenth (1981)
    C_bar=sum(cos(errors$theta-errors$mew_i)/nrow(errors))
    k_inv=2*(1-C_bar)+(1-C_bar)^2*(.48794-.82905*C_bar-1.3915*C_bar^2)/C_bar
    Q_mat=rbind(c(sum(errors$si_star*errors$si),(-0.5)*sum(errors$si_star*errors$ci+errors$ci_star*errors$si)),c((-0.5)*sum(errors$si_star*errors$ci+errors$ci_star*errors$si),sum(errors$ci_star*errors$ci)))
    Q_hat=k_inv*solve(Q_mat)
    if (transmitter$BadPoint[group]==1) transmitter[group,4]<-0 else transmitter[group,4]<-Q_hat[1,1] #Var_X
    if (transmitter$BadPoint[group]==1) transmitter[group,5]<-0 else transmitter[group,5]<-Q_hat[2,2] #Var_Y
    if (transmitter$BadPoint[group]==1) transmitter[group,6]<-0 else transmitter[group,6]<-Q_hat[1,2] #Cov_XY
    
    ## Calculate average angular difference
    if (transmitter$BadPoint[group]==1) transmitter[group,7]<-0 else transmitter[group,7]=mean(abs(errors$mew_i-errors$theta))*180/pi
    
    ## Attach date and time stamps to results variable
    transmitter[group,8]<-(x$Date[x$GID==group])[1]
    transmitter[group,9]<-round(mean(x$Time[x$GID==group]))
  }
    
  ## Return transmitter locations
  class(transmitter)<-"transmitter"
  return(transmitter)
}

plot.receiver<-function(x,add=FALSE,pch=1,cex=1,col=1,bearings=FALSE,...) {
  if (add) points(x$Easting,x$Northing,pch=pch,cex=cex,col=col)
  if (!add) plot(x$Easting,x$Northing,pch=pch,cex=cex,col=col,...)
  if (bearings) segments(x$Easting,x$Northing,x$Easting+1E10*cos(pi/180*(90-x$Azimuth)),x$Northing+1E10*sin(pi/180*(90-x$Azimuth)),lty=3)
}
plot.intersect<-function(x,add=FALSE,pch="*",cex=1,col=4,...) {
  if (add) points(x$X,x$Y,pch=pch,cex=cex,col=col)
  if (!add) plot(x$X,x$Y,pch=pch,cex=cex,col=col,...)
}
plot.transmitter<-function(x,add=FALSE,errors=TRUE,pch=19,cex=1,col=2,badcolor=FALSE,...) {
  if (add) if (badcolor) points(x$X,x$Y,pch=pch,cex=cex,col=x$BadPoint+1) else points(x$X,x$Y,pch=pch,cex=cex,col=col)
  if (!add) if (badcolor) plot(x$X,x$Y,pch=pch,cex=cex,col=x$BadPoint+1,...) else plot(x$X,x$Y,pch=pch,cex=cex,col=col,...)
  if (errors) {
    for (i in 1:nrow(as.data.frame(x))) {
      if (x$BadPoint[i]==0) lines(ellipse(x$Cov_XY[i]/(sqrt(x$Var_X[i])*sqrt(x$Var_Y[i])),scale=c(sqrt(x$Var_X[i]),sqrt(x$Var_Y[i])),centre=c(x$X[i],x$Y[i])))
    }
  }
}

print.receiver<-function(x,...) {
  class(x)<-"data.frame"
  print(x)
}
print.intersect<-function(x,...) {
  class(x)<-"data.frame"
  print(x)
}
print.transmitter<-function(x,...) {
  class(x)<-"data.frame"
  print(x)
}

as.data.frame.receiver<-function(x,row.names=NULL,optional=FALSE, ...) {  
  class(x)<-"data.frame"
  return(x)
}
as.data.frame.intersect<-function(x,row.names=NULL,optional=FALSE, ...) {
  class(x)<-"data.frame"
  return(x)
}
as.data.frame.transmitter<-function(x,row.names=NULL,optional=FALSE, ...) {
  class(x)<-"data.frame"
  return(x)
}