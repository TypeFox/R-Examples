.intersection_onecolumn <-
function(alpha, SSN1, SSN2, nrow_S1, nrow_S2, Ncol, .quantile_intersection){
  SP_inters <- .quantile_intersection(vec1 = SSN1, vec2 = SSN2, a = alpha/2, b = 1 - alpha/2, n = nrow_S1)
  vol_intersect <- (SP_inters[6])
  vol_B <- (SP_inters[3])
  r <- 1
  if (min(SP_inters[4] == SP_inters[1]) == 0 | min(SP_inters[5] == SP_inters[2]) == 0){
    if (vol_intersect < vol_B) {r <- vol_intersect/vol_B}
    if (vol_B == 0) {r <- 0}
  }
  return(r)
}
.intersection_severalcol <-
function(alpha,SSN1,SSN2,nrow_S1,nrow_S2,Ncol, .quantile_intersection){
  SP_inters<-mapply(.quantile_intersection, vec1 = SSN1, vec2 = SSN2, MoreArgs = list(a=((1-(1-alpha)^(1/ncol(SSN1)))/2),  b=(1-(1-(1-alpha)^(1/ncol(SSN1)))/2),n=nrow_S1))  
  product<-prod(SP_inters[6,])
  vol_intersect<-c(product,mean(SP_inters[6,]),product^{1/length(SP_inters[6,])})  
  product<-prod(SP_inters[3,])
  vol_B<-c(product,mean(SP_inters[3,]),product^{1/length(SP_inters[3,])})  
  r<-c(1,1,1)                                                                                             
  if(min(SP_inters[4,]==SP_inters[1,])==0 | min(SP_inters[5,]==SP_inters[2,])==0){                  
    r[vol_intersect<vol_B]<-vol_intersect[vol_intersect<vol_B]/vol_B[vol_intersect<vol_B]
    r[vol_B==0]<-0                                                                             
  }
  return(r)                                                                                        
}
.portionAinB_coordinates_full <-
function(S1,S2,steps=101 ){ 
  nt<-ncol(S1)
  integral_coord<-rep(0,length=nt)
  portionAinB_function<-function(x1,x2,steps){
    r<-.portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$integral_approx)
  }
  portionAinB_function2<-function(x1,x2,steps){
    r<-.portionAinB_full_onecolumn(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$overlap)
  }
  integral_coord<-mapply(portionAinB_function, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  plot_data_overlap<-mapply(portionAinB_function2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  alpha_grid<-seq(0,1,length=steps)[1:(steps-1)]
  erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_data_overlap=plot_data_overlap)
  return(erg)
}
.portionAinB_full_onecolumn <-
function (S1, S2, steps = 101) {
  S <- rbind(S1, S2)
  Ncol <- ncol(S1)
  alpha_grid <- seq(0, 1, length = steps)
  Spans <- data.frame(trait_nr = 1:Ncol, min = min(S, na.rm = TRUE), max = max(S, na.rm = TRUE))
  ab <- Spans$max - Spans$min
  SSN <- (S - Spans$min)/(ab)
  SSN[, ab == 0] <- 0.5
  nrow_S1 <- nrow(S1)
  nrow_S2 <- nrow(S2)
  z <- mapply(.intersection_onecolumn, alpha = alpha_grid,  MoreArgs = list(SSN1 = SSN[1:nrow_S1, ], SSN2 = SSN[(nrow_S1 + 1):(nrow_S1 + nrow_S2), ], nrow_S1 = nrow_S1, nrow_S2 = nrow_S2, Ncol = Ncol, .quantile_intersection = .quantile_intersection))
  integral_approx <- c(prod = trapz(alpha_grid, z))
  erg <- list(alpha_grid = alpha_grid, overlap = z, integral_approx = integral_approx)
  return(erg)
}
.portionAinB2_full <-
function(S1,S2,steps=101,alpha_grid){
  steps0<-steps
  alpha_grid0<-alpha_grid
  .portionAinB2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha_grid0)
}
.portionAinB2_full_severalcol <-
function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)[1:(steps-1)]){
  S<-rbind(S1,S2)
  alpha_grid<-seq(0,1,length=steps) 
  Ncol<-ncol(S1)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min          
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  nrow_S1<-nrow(S1)
  nrow_S2<-nrow(S2)
  z<-mapply(.intersection_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow_S1,], SSN2=SSN[(nrow_S1+1):(nrow_S1+nrow_S2),],nrow_S1=nrow_S1,nrow_S2=nrow_S2,Ncol=Ncol,.quantile_intersection=.quantile_intersection)) 
  integral_approx<-c(prod=trapz(seq(0,1,length=steps), z[1,]),mean=trapz(seq(0,1,length=steps),z[2,]),gmean=trapz(seq(0,1,length=steps),z[3,])) 
  plot_data_prod<-z[1,]
  erg<-list(alpha_grid=seq(0,1,length=steps)[1:(steps-1)],overlap=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod)                                                                                             
  return(erg)
}
.quantile_intersection <-
function(vec1,vec2,a,b,n){
  x<-c( quantile(vec1, probs = c(a,b), na.rm = TRUE),  quantile(vec2, probs = c(a,b), na.rm = TRUE) )   
  x<-c(x,x[4]-x[3])                                      
  x[5]<-ifelse(x[5]<=0,0,x[5])                            
  r<-c(max(x[c(1,3)]),min(x[c(2,4)]))   
  r<-c(x[3:5],r,r[2]-r[1])   
  r[6]<-ifelse(r[6]<=0,0,r[6])                                        
  return(r)                                                            
}
.volume_onecol <-
function(alpha,SSN1){
  SP1<-quantile(SSN1,probs=c(alpha/2,1-alpha/2),na.rm = TRUE)  
  vol<-(SP1[2]-SP1[1])  
  return(vol)                                                                                      
}
.volume_severalcol <-
function(alpha,SSN1){
  alpha<-(1-(1-alpha)^(1/ncol(SSN1))) 
  a<-alpha/2
  b<-1-alpha/2
  SP1<-mapply(quantile, SSN1, MoreArgs = list(probs=c(a,b),na.rm=TRUE))   
  product<-prod(SP1[2,]-SP1[1,])
  vol<-c(prod=product,mean=mean(SP1[2,]-SP1[1,]),gmean=product^{1/length(SP1[2,]-SP1[1,])})   
  return(vol)                                                                                        
}
.volumeA_coordinates_full <-
function(S1,S2,steps=101 ){ 
  nt<-ncol(S1)
  volumeA_function<-function(x1,x2,steps){
    r<-.volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$integral_approx)
  }
  volumeA_function2<-function(x1,x2,steps){
    r<-.volumeA_full_onecol(data.frame(v1=x1) ,data.frame(v1=x2),steps=steps  )
    return(r$volume)
  }
  integral_coord<-mapply(volumeA_function, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  plot_volume<-mapply(volumeA_function2, x1=S1, x2=S2, MoreArgs =list(steps=steps))
  alpha_grid<-seq(0,1,length=steps)    
  erg<-list(alpha_grid=alpha_grid,integral_coord=integral_coord,plot_volume=plot_volume)
  return(erg)
}
.volumeA_full_onecol <-
function(S1,S2,steps=101){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  alpha_grid<-seq(0,1,length=steps) 
  Spans<-data.frame(trait_nr=1:Ncol,min=min(S,na.rm = TRUE),max=max(S,na.rm = TRUE))      
  ab<-Spans$max-Spans$min          
  SSN<-(S-Spans$min)/(ab)          
  SSN[,ab==0]<-0.5                 
  z<-mapply(.volume_onecol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),])) 
  z<-1/(1-alpha_grid[-length(alpha_grid)])*z[-length(alpha_grid)]
  z <- ifelse(z <= 1, z, 1)
  integral_approx<-c(prod=trapz(alpha_grid[-length(alpha_grid)], z))
  erg<-list(alpha_grid=alpha_grid,volume=z,integral_approx=integral_approx) 
  return(erg)
}
.volumeA2_full <-
function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
  steps0<-steps
  alpha0_grid<-alpha_grid
  .volumeA2_full_severalcol(S1,S2,steps=steps0,alpha_grid=alpha0_grid  )
}
.volumeA2_full_severalcol <-
function(S1,S2,steps=101,alpha_grid=seq(0,1,length=steps)){
  S<-rbind(S1,S2)
  Ncol<-ncol(S1)
  Spans<-data.frame(trait_nr=1:Ncol,min=matrix(mapply(min, S, MoreArgs = list(na.rm = TRUE))),max=matrix(mapply(max, S, MoreArgs = list(na.rm = TRUE))))    
  ab<-Spans$max-Spans$min          
  SSN<-t((t(S)-Spans$min)/(ab))    
  SSN[,ab==0]<-0.5                 
  SSN<-as.data.frame(SSN)          
  z<-mapply(.volume_severalcol, alpha = alpha_grid, MoreArgs = list(SSN1=SSN[1:nrow(S1),]))
  z<-rbind(1/(1-alpha_grid[-length(alpha_grid)])*z[1,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[2,][-length(alpha_grid)],1/(1-alpha_grid[-length(alpha_grid)])*z[3,][-length(alpha_grid)])
  z <- ifelse(z <= 1, z, 1)
  integral_approx<-c(prod=trapz(alpha_grid[-length(alpha_grid)], z[1,]),mean=trapz(alpha_grid[-length(alpha_grid)],z[2,]),gmean=trapz(alpha_grid[-length(alpha_grid)],z[3,]))
  plot_data_prod<-z[1,]
  erg<-list(alpha_grid=seq(0,1,length=steps),volume=z,integral_approx=integral_approx,plot_data_prod=plot_data_prod) 
  return(erg)
}
.trpca <- function(data, va){  # data = dataset as for dynRB, va = how much variance (0-1) should included axes explain
  PCA <- prcomp(data[,-1])
    vars <- apply(PCA$x, 2, var)
    prop <- cumsum(vars / sum(vars))
    k <- which(prop >= va)[1]
    k <- ifelse(k==1, 2, k)
    data1 <- data.frame(data[,1], PCA$x[,1:k])
    colnames(data1)[1] <- "Species"
  return(data1)
}

