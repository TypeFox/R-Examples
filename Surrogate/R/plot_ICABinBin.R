plot.ICA.BinBin <- function(x, R2_H=TRUE, R_H=FALSE, Theta_T=FALSE, Theta_S=FALSE, 
                            Type="Density", Labels=FALSE, Xlab.R2_H, Main.R2_H, Xlab.R_H, Main.R_H,
                            Xlab.Theta_S, Main.Theta_S, Xlab.Theta_T, Main.Theta_T, Cex.Legend=1, 
                            Cex.Position="topright", 
                            col, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), ylim, ...){
  
  if (missing(ylim)==TRUE) {Eigen.Y.Lim = 0}
  if (missing(ylim)==FALSE) {Eigen.Y.Lim = 1}
  Object <- x 
  
#  if (C3==TRUE){
#    dev.new()
#    par=Par 
#    if (missing(Xlab.C3)) {Xlab.C3 <- expression(C[3])}
#    if (missing(col)) {col <- c(8)}
#    if (missing(Main.C3)) {Main.C3=" "}  
    
#    if (x$Monotonicity!="General"){
#      plot(density(x$C3, na.rm = T), xlab=Xlab.C3, ylab="Density", main=Main.C3, lwd=2)
#    }
    
#    if (x$Monotonicity=="General"){
#      resul <- cbind(x$Pi.Vectors, x$C3, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
#      colnames(resul) <-  
#        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
#          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "C3", "R2_H", 
#          "Theta_T", "Theta_S", "H_Delta_T")
#      C3_General <- resul$C3
#      C3_No <- resul$C3[resul$Monotonicity=="No"]
#      C3_Surr <- resul$C3[resul$Monotonicity=="Surr"]
#      C3_True <- resul$C3[resul$Monotonicity=="True"]
#      C3_SurrTrue <- resul$C3[resul$Monotonicity=="SurrTrue"]
#      max_val <- max(max(density(C3_No, na.rm = T)$y),  max(density(C3_Surr, na.rm = T)$y), 
#                     max(density(C3_True, na.rm = T)$y),  max(density(C3_SurrTrue, na.rm = T)$y))   
#      plot(density(x$C3, na.rm = T), xlab=Xlab.C3, ylab="Density", main=Main.C3, lwd=2, ylim = c(0, max_val), col=0)
#      try(lines(density((C3_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
#      try(lines(density((C3_Surr), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
#      try(lines(density((C3_True), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
#      try(lines(density((C3_SurrTrue), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
    
#      legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
#             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
      
#    }
    
    
#    if (Type=="Freq"){
#      h <- hist(Object$C3, ...)
#      h$density <- h$counts/sum(h$counts)
#      cumulMidPoint <- ecdf(x=Object$C3)(h$mids)
#      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
#      if (Labels==FALSE){
#        plot(h,freq=T, xlab=Xlab.C3, ylab="Frequency", col=col, main=Main.C3)
#      }
#      if (Labels==TRUE){
#        plot(h,freq=T, xlab=Xlab.C3, ylab="Frequency", col=col, main=Main.C3, labels=labs)
#      }
#    }
    
#    if (Type=="Percent"){
#      h <- hist(Object$C3, ...)
#      h$density <- h$counts/sum(h$counts)
#      cumulMidPoint <- ecdf(x=Object$C3)(h$mids)
#      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
#      if (Labels==FALSE){
#        plot(h,freq=F, xlab=Xlab.C3, ylab="Percentage", col=col, main=Main.C3)
#      }
#      if (Labels==TRUE){
#        plot(h,freq=F, xlab=Xlab.C3, ylab="Percentage", col=col, main=Main.C3, labels=labs)
#      }
#    }
    
#    if (Type=="CumPerc"){
#      h <- hist(Object$C3, breaks=length(Object$C3), ...)
#      h$density <- h$counts/sum(h$counts)
#      cumulative <- cumsum(h$density)
#      plot(x=h$mids, y=cumulative, xlab=Xlab.C3, ylab="Cumulative percentage", col=0, main=Main.C3)
#      lines(x=h$mids, y=cumulative)
#    }    
#  }
  
  if (R2_H==TRUE){

    dev.new()
    par=Par 
    if (missing(Xlab.R2_H)) {Xlab.R2_H <- expression(R[H]^2)}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.R2_H)) {Main.R2_H <- " "}  
    
    if (Type=="Density"){    

    if (x$Monotonicity!="General"){
      plot(density(x$R2_H, na.rm = T), xlab=Xlab.R2_H, ylab="Density", main=Main.R2_H, lwd=2)
    }
    
    if (x$Monotonicity=="General"){
      resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
      colnames(resul) <-  
        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
          "Theta_T", "Theta_S", "H_Delta_T")
      R2_H_General <- resul$R2_H
      R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
      R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
      R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
      R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
      
      try(max1 <- max(density(R2_H_No, na.rm = T)$y), silent=TRUE)
      if (exists("max1")==FALSE){max1 <- 0}
      try(max2 <- max(density(R2_H_Surr, na.rm = T)$y), silent=TRUE)
      if (exists("max2")==FALSE){max2 <- 0}
      try(max3 <- max(density(R2_H_True, na.rm = T)$y), silent=TRUE)
      if (exists("max3")==FALSE){max3 <- 0}
      try(max4 <- max(density(R2_H_SurrTrue, na.rm = T)$y), silent=TRUE)
      if (exists("max4")==FALSE){max4 <- 0}
      
      try(max_val <- max(max1, max2, max3, max4), silent=TRUE)   
      
      if (exists("max_val")==FALSE){max_val <- max(density(R2_H_General)$y)}
      
      if (Eigen.Y.Lim == 0){
      plot(density(x$R2_H, na.rm = T), xlab=Xlab.R2_H, ylab="Density", main=Main.R2_H, lwd=2, ylim = c(0, max_val), col=0)
      }
      if (Eigen.Y.Lim == 1){
        plot(density(x$R2_H, na.rm = T), xlab=Xlab.R2_H, ylab="Density", main=Main.R2_H, lwd=2, ylim = ylim, col=0)
      }
      try(lines(density((R2_H_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
      try(lines(density((R2_H_Surr), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
      try(lines(density((R2_H_True), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
      try(lines(density((R2_H_SurrTrue), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
      #try(lines(density((R2_H_General), na.rm = T), lty=5, col=6, lwd=2), silent=TRUE)
      
      legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
      }
    }
    
    if (Type=="Freq"){
      
    if (x$Monotonicity!="General"){  
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main=Main.R2_H)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main=Main.R2_H, labels=labs)
      }
    }
    
    if (x$Monotonicity=="General"){
      resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
      colnames(resul) <-  
        c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
          "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
          "Theta_T", "Theta_S", "H_Delta_T")
      R2_H_General <- resul$R2_H
      R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
      R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
      R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
      R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
      
      
      par(mfrow=c(2,2))
      # No
      R_hier <- R2_H_No
      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
      if (exists("h")){
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=R_hier)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="No monotonicity")
        }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
        }
      rm(h)
      }
    
    # Surr
    R_hier <- R2_H_Surr
    try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
    if (exists("h")){
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=R_hier)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S")
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
      }
      rm(h)
    }
        
    # True
    R_hier <- R2_H_True
    try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
    if (exists("h")){
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=R_hier)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity T")
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
      }
      rm(h)
    }
        
    # Surr en True
    R_hier <- R2_H_SurrTrue
    try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
    if (exists("h")){
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=R_hier)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S and T")
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
      }
      rm(h)
    }
    
    par(mfrow=c(1,1))
    
    
    }
    }
    
    if (Type=="Percent"){
      
      if (x$Monotonicity!="General"){ 
      
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h, freq=F, xlab=Xlab.R2_H, ylab="Percentage", col=col, main=Main.R2_H)
      }
      if (Labels==TRUE){
        plot(h, freq=F, xlab=Xlab.R2_H, ylab="Percentage", col=col, main=Main.R2_H, labels=labs)
      }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R2_H_General <- resul$R2_H
        R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
        R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
        R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
        R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(2,2))
        # No
        R_hier <- R2_H_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        R_hier <- R2_H_Surr
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
          }
          rm(h)
        }
        
        # True
        R_hier <- R2_H_True
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity T")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
          }
          rm(h)
        }
        
        # Surr en True
        R_hier <- R2_H_SurrTrue
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S and T")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R2_H, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
          }
          rm(h)
        }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="CumPerc"){
      
      if (x$Monotonicity!="General"){ 
      h <- hist(Object$R2_H, breaks=length(Object$R2_H), plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, main=Main.R2_H)
      lines(x=h$mids, y=cumulative)
      }

      if (x$Monotonicity=="General"){ 
        
        resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R2_H_General <- resul$R2_H
        R2_H_No <- resul$R2_H[resul$Monotonicity=="No"]
        R2_H_Surr <- resul$R2_H[resul$Monotonicity=="Surr"]
        R2_H_True <- resul$R2_H[resul$Monotonicity=="True"]
        R2_H_SurrTrue <- resul$R2_H[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(2,2))
        # No
        R_hier <- R2_H_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, 
               main="No Monotonicity")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # Surr
        R_hier <- R2_H_Surr
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity S")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # True
        R_hier <- R2_H_True
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity T")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        
        # Surr en True
        R_hier <- R2_H_SurrTrue
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R2_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity S and T")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        par(mfrow=c(1,1))
        
        
      }
    
    }    
  }
  
  if (R_H==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.R_H)) {Xlab.R_H <- expression(R[H])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.R_H)) {Main.R_H <- expression(R[H])}  
   
    if (Type=="Density"){    
      
      if (x$Monotonicity!="General"){
        plot(density(sqrt(x$R2_H), na.rm = T), xlab=Xlab.R_H, ylab="Density", main=Main.R_H, lwd=2)
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, sqrt(x$R2_H), x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R_H_General <- resul$R_H
        R_H_No <- resul$R_H[resul$Monotonicity=="No"]
        R_H_Surr <- resul$R_H[resul$Monotonicity=="Surr"]
        R_H_True <- resul$R_H[resul$Monotonicity=="True"]
        R_H_SurrTrue <- resul$R_H[resul$Monotonicity=="SurrTrue"]
        max_val <- max(max(density(R_H_No, na.rm = T)$y),  max(density(R_H_Surr, na.rm = T)$y), 
                       max(density(R_H_True, na.rm = T)$y),  max(density(R_H_SurrTrue, na.rm = T)$y))   
        plot(density(sqrt(x$R2_H), na.rm = T), xlab=Xlab.R_H, ylab="Density", main=Main.R_H, lwd=2, ylim = c(0, max_val), col=0)
        try(lines(density((R_H_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
        try(lines(density((R_H_Surr), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
        try(lines(density((R_H_True), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
        try(lines(density((R_H_SurrTrue), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
        #try(lines(density((R_H_General), na.rm = T), lty=5, col=6, lwd=2), silent=TRUE)
        
        legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
               legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
      }
    }
    
    if (Type=="Freq"){
      
      if (x$Monotonicity!="General"){  
        h <- hist(sqrt(Object$R2_H), plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=sqrt(Object$R2_H))(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main=Main.R_H)
        }
        if (Labels==TRUE){
          plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main=Main.R_H, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, sqrt(x$R2_H), x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R_H_General <- resul$R_H
        R_H_No <- resul$R_H[resul$Monotonicity=="No"]
        R_H_Surr <- resul$R_H[resul$Monotonicity=="Surr"]
        R_H_True <- resul$R_H[resul$Monotonicity=="True"]
        R_H_SurrTrue <- resul$R_H[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(2,2))
        # No
        R_hier <- R_H_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        R_hier <- R_H_Surr
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
          }
          rm(h)
        }
        
        # True
        R_hier <- R_H_True
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity T")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
          }
          rm(h)
        }
        
        # Surr en True
        R_hier <- R_H_SurrTrue
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S and T")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
          }
          rm(h)
        }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="Percent"){
      
      if (x$Monotonicity!="General"){ 
        
        h <- hist(sqrt(Object$R2_H), plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=sqrt(Object$R2_H))(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h, freq=F, xlab=Xlab.R_H, ylab="Percentage", col=col, main=Main.R_H)
        }
        if (Labels==TRUE){
          plot(h, freq=F, xlab=Xlab.R_H, ylab="Percentage", col=col, main=Main.R_H, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, sqrt(x$R2_H), x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R_H_General <- resul$R_H
        R_H_No <- resul$R_H[resul$Monotonicity=="No"]
        R_H_Surr <- resul$R_H[resul$Monotonicity=="Surr"]
        R_H_True <- resul$R_H[resul$Monotonicity=="True"]
        R_H_SurrTrue <- resul$R_H[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(2,2))
        # No
        R_hier <- R_H_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        R_hier <- R_H_Surr
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
          }
          rm(h)
        }
        
        # True
        R_hier <- R_H_True
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity T")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
          }
          rm(h)
        }
        
        # Surr en True
        R_hier <- R_H_SurrTrue
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S and T")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.R_H, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
          }
          rm(h)
        }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="CumPerc"){
      
      if (x$Monotonicity!="General"){ 
        h <- hist(sqrt(Object$R2_H), breaks=length(sqrt(Object$R2_H)), plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulative <- cumsum(h$density)
        plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, main=Main.R_H)
        lines(x=h$mids, y=cumulative)
      }
      
      if (x$Monotonicity=="General"){ 
        
        resul <- cbind(x$Pi.Vectors, sqrt(x$R2_H), x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        R_H_General <- resul$R_H
        R_H_No <- resul$R_H[resul$Monotonicity=="No"]
        R_H_Surr <- resul$R_H[resul$Monotonicity=="Surr"]
        R_H_True <- resul$R_H[resul$Monotonicity=="True"]
        R_H_SurrTrue <- resul$R_H[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(2,2))
        # No
        R_hier <- R_H_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, 
               main="No Monotonicity")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # Surr
        R_hier <- R_H_Surr
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity S")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # True
        R_hier <- R_H_True
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity T")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        
        # Surr en True
        R_hier <- R_H_SurrTrue
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.R_H, ylab="Cumulative percentage", col=0, 
               main="Monotonicity S and T")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        par(mfrow=c(1,1))
        
        
      }
      
    }     
       
  }
  
  if (Theta_S==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.Theta_S)) {Xlab.Theta_S <- expression(theta[S])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.Theta_S)) {Main.Theta_S  <- " "}  
    
    if (Type=="Density"){    
      
      if (x$Monotonicity!="General"){
        plot(density(x$Theta_S, na.rm = T), xlab=Xlab.Theta_S, ylab="Density", main=Main.Theta_S, lwd=2)
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_S_General <- resul$Theta_S
        Theta_S_No <- resul$Theta_S[resul$Monotonicity=="No"]
        #      Theta_S_S <- resul$Theta_S[resul$Monotonicity=="Surr"]
        Theta_S_T <- resul$Theta_S[resul$Monotonicity=="True"]
        #      Theta_S_ST <- resul$Theta_S[resul$Monotonicity=="SurrTrue"]
        max_val <- max(try(max(density(Theta_S_No, na.rm = T)$y), silent=TRUE), 
                       #                     try(max(density(Theta_S_S, na.rm = T)$y), silent=TRUE),
                       try(max(density(Theta_S_T, na.rm = T)$y), silent=TRUE) 
                       #                     try(max(density(Theta_S_ST, na.rm = T)$y), silent=TRUE)
        )   
        plot(density(x$Theta_S, na.rm = T), xlab=Xlab.Theta_S, ylab="Density", main=Main.Theta_S, lwd=2, ylim = c(0, max_val), col=0)
        try(lines(density((Theta_S_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
#        try(lines(density((Theta_S_S), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
        try(lines(density((Theta_S_T), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
#        try(lines(density((Theta_S_ST), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
        #try(lines(density((Theta_S_General), na.rm = T), lty=5, col=6, lwd=2), silent=TRUE)
        
        legend(Cex.Position, lwd=c(3, 3), col=c(1, 3), lty=c(3, 3), cex = Cex.Legend,
               legend=c("No monotonicity", "Monotonicity T"))
        
        
        #      legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
        #             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
      }
    }
    
    if (Type=="Freq"){
      
      if (x$Monotonicity!="General"){  
        h <- hist(Object$Theta_S, plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$Theta_S)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main=Main.Theta_S)
        }
        if (Labels==TRUE){
          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main=Main.Theta_S, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_S_General <- resul$Theta_S
        Theta_S_No <- resul$Theta_S[resul$Monotonicity=="No"]
        #      Theta_S_S <- resul$Theta_S[resul$Monotonicity=="Surr"]
        Theta_S_T <- resul$Theta_S[resul$Monotonicity=="True"]
        #      Theta_S_ST <- resul$Theta_S[resul$Monotonicity=="SurrTrue"]
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_S_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        #      R_hier <- Theta_S_S
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        # True
        R_hier <- Theta_S_T
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity T")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
          }
          rm(h)
        }
        
        # Surr en True
        #      R_hier <- Theta_S_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S and T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=T, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="Percent"){
      
      if (x$Monotonicity!="General"){ 
        
        h <- hist(Object$Theta_S, plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$Theta_S)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h, freq=F, xlab=Xlab.Theta_S, ylab="Percentage", col=col, main=Main.Theta_S)
        }
        if (Labels==TRUE){
          plot(h, freq=F, xlab=Xlab.Theta_S, ylab="Percentage", col=col, main=Main.Theta_S, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_S_General <- resul$Theta_S
        Theta_S_No <- resul$Theta_S[resul$Monotonicity=="No"]
        #      Theta_S_S <- resul$Theta_S[resul$Monotonicity=="Surr"]
        Theta_S_T <- resul$Theta_S[resul$Monotonicity=="True"]
        #      Theta_S_ST <- resul$Theta_S[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_S_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        #      R_hier <- Theta_S_S
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        # True
        R_hier <- Theta_S_T
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity T")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
          }
          rm(h)
        }
        
        # Surr en True
        #      R_hier <- Theta_S_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S and T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=F, xlab=Xlab.Theta_S, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="CumPerc"){
      
      if (x$Monotonicity!="General"){ 
        h <- hist(Object$Theta_S, breaks=length(Object$Theta_S), plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulative <- cumsum(h$density)
        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, main=Main.Theta_S)
        lines(x=h$mids, y=cumulative)
      }
      
      if (x$Monotonicity=="General"){ 
        
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_S_General <- resul$Theta_S
        Theta_S_No <- resul$Theta_S[resul$Monotonicity=="No"]
        #      Theta_S_S <- resul$Theta_S[resul$Monotonicity=="Surr"]
        Theta_S_T <- resul$Theta_S[resul$Monotonicity=="True"]
        #      Theta_S_ST <- resul$Theta_S[resul$Monotonicity=="SurrTrue"]
        
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_S_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, 
               main="No Monotonicity")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # Surr
        #      R_hier <- Theta_S_S
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulative <- cumsum(h$density)
        #        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, 
        #             main="Monotonicity S")
        #        lines(x=h$mids, y=cumulative)
        #        rm(h)
        #      }
        
        # True
        R_hier <- Theta_S_T
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, 
               main="Monotonicity T")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        
        # Surr en True
        #      R_hier <- Theta_S_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulative <- cumsum(h$density)
        #        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_S, ylab="Cumulative percentage", col=0, 
        #             main="Monotonicity S and T")
        #        lines(x=h$mids, y=cumulative)
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
        
      }
      
    }
    
  }
  
  if (Theta_T==TRUE){
    dev.new()
    par=Par 
    if (missing(Xlab.Theta_T)) {Xlab.Theta_T <- expression(theta[T])}
    if (missing(col)) {col <- c(8)}
    if (missing(Main.Theta_T)) {Main.Theta_T  <- " "}  
    
    if (Type=="Density"){    
      
      if (x$Monotonicity!="General"){
        plot(density(x$Theta_T, na.rm = T), xlab=Xlab.Theta_T, ylab="Density", main=Main.Theta_T, lwd=2)
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$R2_H, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_T_General <- resul$Theta_T
        Theta_T_No <- resul$Theta_T[resul$Monotonicity=="No"]
        Theta_T_S <- resul$Theta_T[resul$Monotonicity=="Surr"]
        #      Theta_T_T <- resul$Theta_T[resul$Monotonicity=="True"]
        #      Theta_T_ST <- resul$Theta_T[resul$Monotonicity=="SurrTrue"]
        max_val <- max(try(max(density(Theta_T_No, na.rm = T)$y), silent=TRUE), 
                       try(max(density(Theta_T_S, na.rm = T)$y), silent=TRUE)
                       #                     try(max(density(Theta_T_T, na.rm = T)$y), silent=TRUE), 
                       #                     try(max(density(Theta_T_ST, na.rm = T)$y), silent=TRUE)
        )   
        plot(density(x$Theta_T, na.rm = T), xlab=Xlab.Theta_T, ylab="Density", main=Main.Theta_T, lwd=2, ylim = c(0, max_val), col=0)
        try(lines(density((Theta_T_No), na.rm = T), lty=1, col=1, lwd=3), silent=TRUE) 
        try(lines(density((Theta_T_S), na.rm = T), lty=2, col=2, lwd=3), silent=TRUE)
        #      try(lines(density((Theta_T_T), na.rm = T), lty=3, col=3, lwd=3), silent=TRUE)
        #      try(lines(density((Theta_T_ST), na.rm = T), lty=4, col=4, lwd=2), silent=TRUE)
        #try(lines(density((Theta_T_General), na.rm = T), lty=5, col=6, lwd=2), silent=TRUE)
        
        legend(Cex.Position, lwd=c(3, 3), col=c(1, 2), lty=c(1, 2), cex = Cex.Legend,
               legend=c("No monotonicity", "Monotonicity S"))
        
        
        #      legend(Cex.Position, lwd=c(3, 3, 3, 3), col=c(1, 2, 3, 4), lty=c(1, 2, 3, 4), cex = Cex.Legend,
        #             legend=c("No monotonicity", "Monotonicity S", "Monotonicity T", "Monotonicity S and T"))
      }
    }
    
    if (Type=="Freq"){
      
      if (x$Monotonicity!="General"){  
        h <- hist(Object$Theta_T, plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$Theta_T)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main=Main.Theta_T)
        }
        if (Labels==TRUE){
          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main=Main.Theta_T, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_T_General <- resul$Theta_T
        Theta_T_No <- resul$Theta_T[resul$Monotonicity=="No"]
        Theta_T_S <- resul$Theta_T[resul$Monotonicity=="Surr"]
        #      Theta_T_T <- resul$Theta_T[resul$Monotonicity=="True"]
        #      Theta_T_ST <- resul$Theta_T[resul$Monotonicity=="SurrTrue"]
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_T_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        R_hier <- Theta_T_S
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S")
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
          }
          rm(h)
        }
        
        # True
        #      R_hier <- Theta_T_T
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        # Surr en True
        #      R_hier <- Theta_T_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S and T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=T, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="Percent"){
      
      if (x$Monotonicity!="General"){ 
        
        h <- hist(Object$Theta_T, plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulMidPoint <- ecdf(x=Object$Theta_T)(h$mids)
        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        
        if (Labels==FALSE){
          plot(h, freq=F, xlab=Xlab.Theta_T, ylab="Percentage", col=col, main=Main.Theta_T)
        }
        if (Labels==TRUE){
          plot(h, freq=F, xlab=Xlab.Theta_T, ylab="Percentage", col=col, main=Main.Theta_T, labels=labs)
        }
      }
      
      if (x$Monotonicity=="General"){
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_T_General <- resul$Theta_T
        Theta_T_No <- resul$Theta_T[resul$Monotonicity=="No"]
        Theta_T_S <- resul$Theta_T[resul$Monotonicity=="Surr"]
        #      Theta_T_T <- resul$Theta_T[resul$Monotonicity=="True"]
        #      Theta_T_ST <- resul$Theta_T[resul$Monotonicity=="SurrTrue"]
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_T_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="No monotonicity")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="No monotonicity", labels=labs)
          }
          rm(h)
        }
        
        # Surr
        R_hier <- Theta_T_S
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=R_hier)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S")
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S", labels=labs)
          }
          rm(h)
        }
        
        # True
        #      R_hier <- Theta_T_T
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        # Surr en True
        #      R_hier <- Theta_T_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulMidPoint <- ecdf(x=R_hier)(h$mids)
        #        labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
        #        
        #        if (Labels==FALSE){
        #          plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S and T")
        #        }
        #        if (Labels==TRUE){
        #          plot(h,freq=F, xlab=Xlab.Theta_T, ylab="Frequency", col=col, main="Monotonicity S and T", labels=labs)
        #        }
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
        
      }
    }
    
    if (Type=="CumPerc"){
      
      if (x$Monotonicity!="General"){ 
        h <- hist(Object$Theta_T, breaks=length(Object$Theta_T), plot = FALSE, ...)
        h$density <- h$counts/sum(h$counts)
        cumulative <- cumsum(h$density)
        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, main=Main.Theta_T)
        lines(x=h$mids, y=cumulative)
      }
      
      if (x$Monotonicity=="General"){ 
        
        resul <- cbind(x$Pi.Vectors, x$Theta_S, x$Theta_T, x$Theta_S, x$H_Delta_T)
        colnames(resul) <-  
          c("Pi_0000", "Pi_0100", "Pi_0010", "Pi_0001", "Pi_0101", "Pi_1000", "Pi_1010", "Pi_1001", "Pi_1110", "Pi_1101", "Pi_1011", 
            "Pi_1111", "Pi_0110", "Pi_0011", "Pi_0111", "Pi_1100", "Sum.Pi.f", "Monotonicity", "R2_H", 
            "Theta_T", "Theta_S", "H_Delta_T")
        Theta_T_General <- resul$Theta_T
        Theta_T_No <- resul$Theta_T[resul$Monotonicity=="No"]
        Theta_T_S <- resul$Theta_T[resul$Monotonicity=="Surr"]
        #      Theta_T_T <- resul$Theta_T[resul$Monotonicity=="True"]
        #      Theta_T_ST <- resul$Theta_T[resul$Monotonicity=="SurrTrue"]
        
        par(mfrow=c(1,2))
        # No
        R_hier <- Theta_T_No
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, 
               main="No Monotonicity")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # Surr
        R_hier <- Theta_T_S
        try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        if (exists("h")){
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, 
               main="Monotonicity S")
          lines(x=h$mids, y=cumulative)
          rm(h)
        }
        
        # True
        #      R_hier <- Theta_T_T
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulative <- cumsum(h$density)
        #        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, 
        #             main="Monotonicity T")
        #        lines(x=h$mids, y=cumulative)
        #        rm(h)
        #      }
        
        
        # Surr en True
        #      R_hier <- Theta_T_ST
        #      try(h <- hist(R_hier, plot = FALSE, ...), silent=TRUE)
        #      if (exists("h")){
        #        h$density <- h$counts/sum(h$counts)
        #        cumulative <- cumsum(h$density)
        #        plot(x=h$mids, y=cumulative, xlab=Xlab.Theta_T, ylab="Cumulative percentage", col=0, 
        #             main="Monotonicity S and T")
        #        lines(x=h$mids, y=cumulative)
        #        rm(h)
        #      }
        
        par(mfrow=c(1,1))
        
      }
      
    }
  
  }
}
