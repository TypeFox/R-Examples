ORICC2 <-
function (data,data.col,id.col=NULL,n.rep,n.top=NULL,transform=NULL, name.profile=NULL, 
                   cyclical.profile=NULL, onefile=NULL,plot.format=NULL){

   if (is.null(id.col)){
       id.col <- 1;
       id <- 1:nrow(data);
       data <- cbind(id,data);
       data.col <- data.col+1;
   }
   if (is.null(transform))  {transform <- 0}
   if (is.null(cyclical.profile))  {cyclical.profile <- 0}        
   if (is.null(onefile))  {onefile <- TRUE}
   if (is.null(plot.format)) {plot.format <- "eps"}  
  
   x <- data;
   max.inf <- 1e30;


   if (transform==0)  x[,data.col] <- x[,data.col]
   if (transform==1)  x[,data.col] <- log(x[,data.col])
   if (transform==2)  x[,data.col] <- sqrt(x[,data.col])
   if (transform==3)  x[,data.col] <- x[,data.col]^(1/3)

   qq <- nrow(x)
   n <- sum(n.rep)
   n.time <- length(n.rep)
   uo <- 1
   for(i in 1:n.time){
      uo <- uo+1/i
   }
   vv <- log(n)

   
   n.name.profile <- 0;
   n.profile <- 0;
   if (is.null(name.profile)==FALSE){
       if (name.profile[1]=="all") profile <- 1:(2*n.time-2);
       n.name.profile <- length(name.profile)
       name.pro1 <- paste("up down max at",2:(n.time-1),sep=" ")
       name.pro2 <- paste("down up min at",2:(n.time-1),sep=" ")
       all.pro <- c("decreasing",name.pro1,"increasing",name.pro2)
       if(name.profile[1]!="all"){
            profile <- rep(0,n.name.profile)
            for(i in 1:n.name.profile){
                 for(j in 1:(2*n.time-2)){
                     if (name.profile[i]==all.pro[j]) profile[i] <- j       
                 }
            }
        }
        if (any(profile==0))  stop(" The input of name.profile is wrong!")
         n.profile <- length(profile)
         profile <- sort(profile);
    }

   
  
   
    
    
   cyclical.profile <- as.matrix(cyclical.profile)
   n.cyclical.profile <- 0;
   if(max(cyclical.profile)>0){ 
      if (ncol(cyclical.profile)!=2 )  stop(" There must be two column in cyclical!")
      if (min(cyclical.profile)<2&cyclical.profile[1,1]!=0 )  stop("the minimum time point must be larger than 2 in cyclical !")
      if (any(cyclical.profile[,1]==cyclical.profile[,2]) )  stop("the minimum time point and the maximum time point must not be the same in cyclical!")
      if (max(cyclical.profile)>(n.time-1)) stop("the maximum time point must be smaller than n.time in cyclical !")
      n.cyclical.profile <- nrow(cyclical.profile)
      name.cyclical.profile <- c(0,n.cyclical.profile)
      for (i in 1:n.cyclical.profile){
         name.cyclical.profile[i] <- paste("cyclical min at",cyclical.profile[i,1],"max at",cyclical.profile[i,2] ,sep=" ")
      }

   }



    n.all.profile <- 2
    deta <- matrix(0,n.all.profile,qq)
    bic <- rep(0,n.all.profile)
    uu <- matrix(0,nrow=qq,ncol=n.time)
    for(j in 1:qq){
        bic1 <- rep(0,n.all.profile)
        mu <- matrix(0,n.all.profile,n.time)
        x1 <- x[j,data.col]
        x.mean <- rep(0,n.time)
        len1 <- 1;
        len2 <- 0;
        x1 <- as.numeric(x1)
        for(i in 1:n.time){
            len2 <- len2+n.rep[i]
            x.mean[i] <- mean(x1[len1:len2])
            len1 <- len1+n.rep[i]
        }
        um1 <- complete.profile(x1,x.mean,n.rep)
        mu[1,] <- um1$mu
        h1 <- um1$logelr
        bic1[1] <- -2*(h1)+(1+n.time)*vv;
     
        um1 <- flat.pattern(x1,x.mean,n.rep)
        mu[2,] <- um1$mu
        h1 <- um1$logelr
        bic1[2] <- -2*(h1)+2*vv;
     
        re <- bic1
        mm <- min(re)
        idmin <- which.min(re)
        deta[idmin,j] <- 1;
        uu[j,] <- mu[idmin,]
        rf <- rep(0,n.all.profile)
        rf[re==mm] <- 1
        bic <- bic+rf
       # print(bic)
    }

    xx <- data[deta[n.all.profile,]==0,]
    xx2 <- x[deta[n.all.profile,]==0,]

    qq <- nrow(xx2)
    complete.profile <- 1;
    n.all.profile <- n.profile+n.cyclical.profile+complete.profile;
    deta <- matrix(0,n.all.profile,qq)
    bic <- rep(0,n.all.profile)
    uu <- matrix(0,nrow=qq,ncol=n.time)



    for(j in 1:qq){
        bic1 <- rep(max.inf,n.all.profile)
        mu <- matrix(0,n.all.profile,n.time)

        x1 <- xx2[j,data.col]
        x.mean <- rep(0,n.time)
        len1 <- 1;
        len2 <- 0;
        x1 <- as.numeric(x1)
        for(i in 1:n.time){
            len2 <- len2+n.rep[i]
            x.mean[i] <- mean(x1[len1:len2])
            len1 <- len1+n.rep[i]
        }

        k <- 1;
        if (is.null(name.profile)==FALSE){
             if (any(profile==1)){
                  um1 <- decreasing(x1,x.mean,n.rep)
                  mu[k,] <- um1$mu
                  h1 <- um1$logelr
                  if(mu[k,1]>mu[k,n.time])   bic1[k] <- -2*(h1)+uo*vv;
                  k <- k+1;
              }
              for(i in 2:(n.time-1)){
                  if(any(profile==i )){
                       um1 <- up.down(x1,x.mean,n.rep,i)
                       mu[k,] <- um1$mu
                       h1 <- um1$logelr
                       if((mu[k,1]<mu[k,i])&& (mu[k,i]>mu[k,n.time]) )   bic1[k] <- -2*(h1)+uo*vv;
                       k=k+1;
                  }
              }
              if(any(profile==n.time) ){
                  um1 <- increasing(x1,x.mean,n.rep)
                  mu[k,] <- um1$mu
                  h1 <- um1$logelr
                  if(mu[k,1]<mu[k,n.time])   bic1[k] <- -2*(h1)+uo*vv;
                  k <- k+1;
              }
              for(i in 2:(n.time-1)){
                  if(any(profile==(n.time+i-1)) ){
                       um1 <- down.up(x1,x.mean,n.rep,i)
                       mu[k,] <- um1$mu
                       h1 <- um1$logelr
                       if((mu[k,1]>mu[k,i])&&(mu[k,i]<mu[k,n.time]) )    bic1[k] <- -2*(h1)+uo*vv;
                       k <- k+1
                  }
              }
          }     
     
          if(n.cyclical.profile>0){
               for (i in 1:n.cyclical.profile){
                   if (cyclical.profile[i,1]>cyclical.profile[i,2]){
                        um1 <- cyclical.max.min(x1,x.mean,n.rep,cyclical.profile[i,2],cyclical.profile[i,1])
                        mu[k,] <- um1$mu
                        h1 <- um1$logelr
                        if((mu[k,1]<mu[k,cyclical.profile[i,2]])&&(mu[k,cyclical.profile[i,2]]>mu[k,cyclical.profile[i,1]])&&(mu[k,cyclical.profile[i,1]]<mu[k,n.time])) 
                        bic1[k] <- -2*(h1)+uo*vv;
                        k <- k+1
                   }
                   if (cyclical.profile[i,1]<cyclical.profile[i,2]){
                        um1 <- cyclical.min.max(x1,x.mean,n.rep,cyclical.profile[i,1],cyclical.profile[i,2])
                        mu[k,] <- um1$mu
                        h1 <- um1$logelr
                        if((mu[k,1]>mu[k,cyclical.profile[i,1]])&&(mu[k,cyclical.profile[i,1]]<mu[k,cyclical.profile[i,2]])&&(mu[k,cyclical.profile[i,2]]>mu[k,n.time])) 
                        bic1[k] <- -2*(h1)+uo*vv;
                        k <- k+1
                   }                
              }  
          }     
          if(complete.profile==1){
              um1 <- complete.profile(x1,x.mean,n.rep)
              mu[k,] <- um1$mu
              h1 <- um1$logelr
              bic1[k] <- -2*(h1)+(1+n.time)*vv;
              k <- k+1
          }
          re <- bic1
          mm <- min(re)
          idmin <- which.min(re)
          deta[idmin,j] <- 1;
          uu[j,] <- mu[idmin,]
          rf <- rep(0,n.all.profile)
          rf[re==mm] <- 1
          bic <- bic+rf
        #  print(bic)
    }

   suu <- uu
   inn <- nrow(suu)
    vb <- apply(suu,1,sd)
    vb <- vb^2
    if (is.null(n.top))  {
   n.top <- inn
   ino <- 1:inn

   }    
 if(n.top<inn){

    svb <- rev(sort(vb))
   ino <- rev(order(vb))

}







 
    if(n.top<=inn){
        ino <- ino[1:n.top]; 
        infor1 <- paste("ORICC select", inn, "genes");
        print(infor1)
        infor2 <- paste("retain the top", n.top, "genes");
        print(infor2)
    }
    if(n.top>inn ){
        infor1 <- paste("ORICC select", inn, "genes");
        print(infor1)
        infor2 <- paste("retain the top", inn, "genes");
        print(infor2)
    }

    top.id <- xx[ino,id.col];

    match <- rep(0,length(ino))
    for(j in 1:length(ino)){
        for( i in 1:n.all.profile){
           if(deta[i,ino[j]]==1) match[j] <-  i
        }
    }
 
    cluster <- rep(0,length(ino));
    k.cluster <- 1;
    if (n.profile>0){
       for(j.cluster in 1:n.profile){
         flag <- 0;  
         for(i in 1:length(match)){       
            if(match[i]==j.cluster) {cluster[i] <- k.cluster;flag <- 1;}
         }
         if(flag==1){        
            k.cluster <- k.cluster+1;
                  
         }
      }
   }
   if(n.cyclical.profile>0){
      for(j.cluster in (n.profile+1):(n.profile+n.cyclical.profile)){
         flag <- 0; 
         for(i in 1:length(match)){       
            if(match[i]==j.cluster) {cluster[i] <- k.cluster;flag <- 1;}
         }
          if(flag==1){        
             k.cluster <- k.cluster+1;
                 
         }
      }
   }

   if(complete.profile==1){
      flag <- 0; 
      for(i in 1:length(match)){       
         if(match[i]==n.all.profile) {cluster[i] <- k.cluster;flag <- 1;}
      }
       if(flag==1){             
          k.cluster <- k.cluster+1;       
      }          
   } 

  #### write the data for each cluster into external txt file
        Cluster <- 1;
        name.output.file <- "cluster of raw data.txt"
        this.cluster.mat1 <- c()
        if (n.profile>0){
            for(j.cluster in 1:n.profile){
                this.cluster.mat <- c()  
                for(i in 1:length(match)){       
                    if(match[i]==j.cluster) {this.cluster.mat <- rbind(this.cluster.mat,xx[ino[i],])}
                }
                if(is.null(nrow(this.cluster.mat))==FALSE){ 
                    this.cluster.mat <- cbind(Cluster,this.cluster.mat)    
                    Cluster <- Cluster+1;  
                }
        
                this.cluster.mat1 <- rbind(this.cluster.mat1,this.cluster.mat) 
            }
        }

        if(n.cyclical.profile>0){
             for(j.cluster in (n.profile+1):(n.profile+n.cyclical.profile)){
                 this.cluster.mat <- c()  
                 for(i in 1:length(match)){       
                     if(match[i]==j.cluster) {this.cluster.mat <- rbind(this.cluster.mat,xx[ino[i],])}
                 }
                 if(is.null(nrow(this.cluster.mat))==FALSE){ 
                     this.cluster.mat <- cbind(Cluster,this.cluster.mat)       
                     Cluster <- Cluster+1;
                 }
         
                 this.cluster.mat1 <- rbind(this.cluster.mat1,this.cluster.mat) 
             }
        }

        if(complete.profile==1){
             this.cluster.mat <- c()  
             for(i in 1:length(match)){       
                 if(match[i]==(n.profile+n.cyclical.profile+1)) {this.cluster.mat<- rbind(this.cluster.mat,xx[ino[i],])}
             }
             if(is.null(nrow(this.cluster.mat))==FALSE){ 
                 this.cluster.mat<- cbind(Cluster,this.cluster.mat)           
                 Cluster<- Cluster+1;    
             }
   
             this.cluster.mat1<- rbind(this.cluster.mat1,this.cluster.mat)      
        } 
        write.table(this.cluster.mat1, name.output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  
   


    if (plot.format!="eps" && plot.format!="jpg") stop("the plot.format must be 'eps' or 'jpg' !")  
    #### create the data for making plots
    ### create the  series of each cluster
    name.output.file <- "cluster of fitted mean data.txt"
    k.cluster <- 1;
    data.for.this.cluster1 <- c()
    if (n.profile>0){
        for (j.cluster in 1:n.profile){
             data.for.this.cluster=c()  
             for (i in 1:length(match)){       
                  if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,c(xx[ino[i],id.col],vb[ino[i]],suu[ino[i],]))}
             }
             if(is.null(nrow(data.for.this.cluster))==FALSE){
                  data.for.this.cluster <- cbind(k.cluster,data.for.this.cluster)
                  k.cluster <- k.cluster+1;                 
             }
             data.for.this.cluster1 <- rbind(data.for.this.cluster1,data.for.this.cluster) 

       }
    }

    if(n.cyclical.profile>0){
         for (j.cluster in (n.profile+1):(n.profile+n.cyclical.profile)){
              data.for.this.cluster <- c()  
              for (i in 1:length(match)){       
                   if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,c(xx[ino[i],id.col],vb[ino[i]],suu[ino[i],]))}
              }
              if(is.null(nrow(data.for.this.cluster))==FALSE){
                   data.for.this.cluster <- cbind(k.cluster,data.for.this.cluster)
                   k.cluster <- k.cluster+1; 
              }
              data.for.this.cluster1 <- rbind(data.for.this.cluster1,data.for.this.cluster) 

        }
   }

   if(complete.profile==1){
       data.for.this.cluster=c()  
       for (i in 1:length(match)){       
           if (match[i]==(n.profile+n.cyclical.profile+1)) {data.for.this.cluster <- rbind(data.for.this.cluster,c(xx[ino[i],id.col],vb[ino[i]],suu[ino[i],]))}
       }
       if(is.null(nrow(data.for.this.cluster))==FALSE){
           data.for.this.cluster <- cbind(k.cluster,data.for.this.cluster)
           k.cluster <- k.cluster+1; 
       }
       data.for.this.cluster1 <- rbind(data.for.this.cluster1,data.for.this.cluster) 
   }
   name1 <- paste("Time ",1:n.time,sep="")
   colnames(data.for.this.cluster1) <- c("Cluster","id","Vg" ,name1)
   write.table(data.for.this.cluster1, name.output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)  

      if (onefile==FALSE){  
 
       ### create the  series of each cluster
       k.cluster=1;
       if (n.profile>0){
            for (j.cluster in 1:n.profile){
            #### create the data for making plots
                 data.for.this.cluster <- c()  
                 for (i in 1:length(match)){       
                     if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
                 }
                 ### now produce plots of the series
                 name.title <- paste("cluster ", k.cluster, ": ",all.pro[profile[j.cluster]],sep="")
                 cov.this.cluster <- c(1:n.time)
                 size.this.cluster <- nrow(data.for.this.cluster)
                 if(is.null(size.this.cluster)==FALSE){
                     if(plot.format=="eps"){
                          name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".eps", sep="")      
                          postscript(name.plot.profile)
                     }
                     if(plot.format=="jpg"){
                          name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".jpg", sep="")                         
                          jpeg(name.plot.profile)
                     }
                     k.cluster <- k.cluster+1;
                     plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression",main=name.title )
                     if(size.this.cluster>1) {
                          for(ind.series in 2:size.this.cluster) {
                              points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                          }
                     }             
                     dev.off()
                 }
            }
       }

       if(n.cyclical.profile>0){
            for (j.cluster in (n.profile+1):(n.profile+n.cyclical.profile)){
                 data.for.this.cluster <- c()  
                 for (i in 1:length(match)){       
                      if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
                 }
                 ### now produce plots of the series
                 name.title <- paste("cluster ", k.cluster, ": ",name.cyclical.profile[j.cluster-n.profile],sep="")
                 cov.this.cluster <- c(1:n.time)
                 size.this.cluster <- nrow(data.for.this.cluster)
                 if(is.null(size.this.cluster)==FALSE){
                      if(plot.format=="eps"){
                           name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".eps", sep="")
                           postscript(name.plot.profile)
                      } 
                      if(plot.format=="jpg"){
                           name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".jpg", sep="")
                           jpeg(name.plot.profile)
                     }
                      k.cluster <- k.cluster+1;
                      plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression", main=name.title )
                      if(size.this.cluster>1) {
                           for(ind.series in 2:size.this.cluster) {
                                points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                           }
                      }             
                      dev.off()
                 }
             }
        }

       if(complete.profile==1){
            #### create the data for making plots
            data.for.this.cluster <- c()  
            for (i in 1:length(match)){       
                 if (match[i]==(n.profile+n.cyclical.profile+1)) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
            }
            ### now produce plots of the series
            name.title <- paste("cluster ", k.cluster, ": ", "complete profile", sep="")
            cov.this.cluster <- c(1:n.time)
            size.this.cluster=nrow(data.for.this.cluster)
            if(is.null(size.this.cluster)==FALSE){
                 if(plot.format=="eps"){
                      name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".eps", sep="")
                      postscript(name.plot.profile)
                 }
                if(plot.format=="jpg"){
                      name.plot.profile <- paste("cluster ", k.cluster ," of fitted mean data", ".jpg", sep="")  
                      jpeg(name.plot.profile)
                }
                 k.cluster=k.cluster+1;
                 plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression",main=name.title )
        
                 if(size.this.cluster>1) {
                     for(ind.series in 2:size.this.cluster) {
                          points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                     }
                 }             
                 dev.off()
            }
      }
   }

   if (onefile==TRUE){  
         nn=length(as.numeric(levels(as.factor(match))))
         ### create the  series of each cluster
         if (nn==1)  {a <- 1;b <- 1;}
         if (nn==2)  {a <- 1;b <- 2;}
         if (nn>=3 && nn<=6){  a <- 2; b <- floor(nn/(a+0.1))+1; }
         if (nn>=7 && nn<=12){a <- 3; b <- floor(nn/(a+0.1))+1;}
         if (nn>=13 ){ a <- 4; b <- floor(nn/(a+0.1))+1;}
    
          if(plot.format=="eps"){
              name.plot <- "cluster of fitted mean data.eps";
              postscript(name.plot);
          }
          if(plot.format=="jpg"){
              name.plot <- "cluster of fitted mean data.jpg";
              jpeg(name.plot);
          }

          par(mfrow=c(a,b))

          k.cluster <- 1;
          if (n.profile>0){
               for (j.cluster in 1:n.profile){
                    #### create the data for making plots
                    data.for.this.cluster <- c()  
                    for (i in 1:length(match)){       
                         if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
                    }
                    ### now produce plots of the series   
                    cov.this.cluster <- c(1:n.time)
                    size.this.cluster <- nrow(data.for.this.cluster)
                    if(is.null(size.this.cluster)==FALSE){
                         name.title <- paste("cluster ", k.cluster, ": ", all.pro[profile[j.cluster]], sep="")
                         k.cluster <- k.cluster+1;
                         plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression",main=name.title )
                         if(size.this.cluster>1) {
                              for(ind.series in 2:size.this.cluster) {
                                   points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                              } 
                         }
                    }
                }
           }

           if(n.cyclical.profile>0){
                for (j.cluster in (n.profile+1):(n.profile+n.cyclical.profile)){
                    data.for.this.cluster <- c()  
                    for (i in 1:length(match)){       
                        if (match[i]==j.cluster) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
                    }
                    ### now produce plots of the series       
                     cov.this.cluster <- c(1:n.time)
                     size.this.cluster <- nrow(data.for.this.cluster)
                     if(is.null(size.this.cluster)==FALSE){
                          name.title <- paste("cluster ", k.cluster, ": ", name.cyclical.profile[j.cluster-n.profile], sep="")
                          k.cluster <- k.cluster+1;
                          plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression", main=name.title )
                          if(size.this.cluster>1) {
                               for(ind.series in 2:size.this.cluster) {
                                    points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                               }
                          }    
                     }
                }
           }

          if(complete.profile==1){
               #### create the data for making plots
               data.for.this.cluster <- c()  
               for (i in 1:length(match)){       
                   if (match[i]==(n.profile+n.cyclical.profile+1)) {data.for.this.cluster <- rbind(data.for.this.cluster,suu[ino[i],])}
               }
               ### now produce plots of the series        
               cov.this.cluster <- (1:n.time)
               size.this.cluster <- nrow(data.for.this.cluster)
               if(is.null(size.this.cluster)==FALSE){
                   name.title <- paste("cluster ", k.cluster, ": ", "complete profile", sep="")
                   k.cluster <- k.cluster+1;
                   plot(cov.this.cluster, data.for.this.cluster[1, ], type="l",ylim=range(data.for.this.cluster, na.rm=TRUE), xlab="Time", ylab="Expression",main=name.title )
                   if(size.this.cluster>1) {
                        for(ind.series in 2:size.this.cluster) {
                             points(cov.this.cluster, data.for.this.cluster[ind.series, ], type="l")
                        }
                   }
               }
          }
          dev.off()
      }

       
list( cluster=cluster, top.id=top.id)
}

