score.GGUM <-
  function(itemPars,                      
           respData,
           numCats,
           recode = FALSE,
           scoreMissing = TRUE,
           removePersons = NULL,
           respCutoff = NULL,
           numQuads = 30,
           plotTheta = TRUE,
           outFile = NULL) {
    
    # Store item parameter estimates in "pars" and response data in "data"
    
    pars <- as.matrix(itemPars)
    data <- as.matrix(respData)
    
    # Check length of numCats and respCutoff
    
    if (length(numCats)==1) {numCats <- rep(numCats,nrow(pars))} 
    if (length(respCutoff)==1) {respCutoff <- rep(respCutoff,nrow(pars))}
    
    # Recode response values from 1:C to 0:(C-1) for all items
    
    if (recode==TRUE) { data[,2:ncol(data)] <- data[,2:ncol(data)] - 1 }
    
    # Add column names to response file
    
    data.names <- character(ncol(data))   
    data.names[1] <- "SUBJECT"    
    data.names[2:ncol(data)] <- paste(c("ITEM"), 1:(ncol(data)-1), sep="")     
    colnames(data) <- data.names
    
    # Remove persons from response matrix that match removePersons vector
    
    if (!is.null(removePersons)) {
      remove <- numeric(0)
      for (j in 1:nrow(data)) {
        if (data[j,1] %in% removePersons) {remove <- c(remove,j)}
      }
      remove <- as.vector(remove)
      removePersons <- remove*-1
      if (length(removePersons)>0) { data <- data[removePersons,] }
    }
    
    # Remove persons that do not meet minimum agreement response cutoff
    
    if (!is.null(respCutoff)) {  
      meets <- rep(FALSE,(ncol(data)-1))
      remove <- numeric(0)
      for (j in 1:nrow(data)) {
        for (i in 2:ncol(data)) {
          if (!is.na(data[j,i]) & (data[j,i]>=respCutoff[i-1])) {meets[i-1] <- TRUE}
        }
      if (!any(meets)) {remove <- c(remove,j)}
      }
      removePersons <- remove*-1
      if (length(removePersons)>0) { data <- data[removePersons,] }
    }
    
    # Stop program and print error if input values are not correct dimensions
    
    if (nrow(pars)!=(ncol(data)-1)){
      stop("ERROR: number of rows in itemPars must be equal to number of columns minus one in respData")
    }
    
    if ((length(numCats)!=1) & (length(numCats)!=nrow(pars))) {
      print(length(numCats))
      stop("ERROR: length of numCats not equal to 1 or to the number of items")
    }
    
    if (!is.null(respCutoff)) {
      if ((length(respCutoff)!=1) & (length(respCutoff)!=nrow(pars))) {
        stop("ERROR: length of respCutoff not equal to 1 or to the number of items")
      }
    }
    
    # Stop program and print error if input values are out of range
    
    if ((min(numCats)<2) | (max(numCats)>(ncol(pars)-3))) {
      stop("ERROR: number of response categories is out of range
           \n:number of response categories must be between 2 and C for each item")
    }
    
    if (!is.null(respCutoff)) {
      if ((min(respCutoff)<0) | (max(respCutoff)>(ncol(pars)-4))) {
        stop("ERROR: response cutoffs are out of range
           \n:response cutoffs must be between 0 and C-1")
      }
    }
    
    # Calculate quadrature points and normal densities based on user-specified number of quadrature points
    
    point <- seq(-4,4,8/(numQuads-1))
    density <- dnorm(point)
    density.r <- density/sum(density)
    quads <- matrix(c(point,density),ncol=2)
    colnames(quads)[1:2] <- c("point","density")
    
    # Create empty columns for storing response probabilities and estimates
    
    M = numCats*2 - 1
    calc <- matrix(NA,nrow=nrow(pars),ncol=3)
    colnames(calc) <- c("num","den","prob")
    pars <- cbind(pars,M,calc)

    EAP <- matrix(NA,nrow=numQuads,ncol=5)
    colnames(EAP) <- c("EAP_num","EAP_den","like","std_num","std_den")
    estimates <- matrix(NA,nrow=nrow(data),ncol=3)
    colnames(estimates) <- c("Subject","Theta_j","Std_j")
    
    # Create labels for TAU parameters
    
    thresh.labels <- character(10)
    for (k in 0:9) {
      thresh.labels[k+1] <- paste("TAU",k,sep="")
    }
    
    for (k in 0:(max(numCats)-1)) {
      x <- pars[,paste("TAU",k,sep="")]
      pars <- cbind(pars,x)
      colnames(pars)[ncol(pars)] <- paste("TAUd",k,sep="")
    }
      
    # Function for calculating numerator of GGUM response probability function 
      
    num <- function(a,numCats,b,scoreMissing=TRUE) { 
        x <- (exp(a["ALPHA"]*(a["z"]*(b["point"]-a["DELTA"]) - sum(a[4:(4+numCats-1)],na.rm=scoreMissing))) +      
                exp(a["ALPHA"]*((a["M"]-a["z"])*(b["point"]-a["DELTA"]) - sum(a[4:(4+numCats-1)],na.rm=scoreMissing))))   
        return(x)
    }
      
    # Function for calculating denominator of GGUM response probability function
      
    den <- function(a,numCats,b,scoreMissing=TRUE) {
        denom.sum <- 0
        for (w in 0:((a["M"]-1)/2)) {
          x <- (exp(a["ALPHA"]*(w*(b["point"]-a["DELTA"]) - sum(a[(8+numCats):(8+numCats+w)],
                                                                     na.rm=scoreMissing))) +      
                  exp(a["ALPHA"]*((a["M"]-w)*(b["point"]-a["DELTA"]) - sum(a[(8+numCats):(8+numCats+w)],
                                                                                na.rm=scoreMissing))))
          denom.sum <- denom.sum + x
        }
        return(denom.sum)
    }
      
    # Function for calculating GGUM response probabilities at each quadrature point
    
    quad.fun <- function(quads,a,numCats,scoreMissing=TRUE) { 
        b <- quads
        a[,"num"] <- apply(a,1,num,numCats=max(numCats),b,scoreMissing) 
        a[,"den"] <- apply(a,1,den,numCats=max(numCats),b,scoreMissing)     
        a[,"prob"] <- a[,"num"]/a[,"den"]            
        for (na.replace in 1:nrow(a)) {
          if (is.na(a[na.replace,"z"])) { a[na.replace,"den"] <- NA } 
        }       
        like <- prod(a[,"prob"],na.rm=scoreMissing)
        EAP_num <- quads["point"]*like*quads["density"]
        EAP_den <- like*quads["density"]        
        EAP.row <- numeric(3)
        EAP.row <- cbind(EAP_num,EAP_den,like)       
        return(EAP.row) 
    }
    
    # Function for calculating EAP estimate and associated standard deviation for each person
    
    person.score <- function(data,pars,numCats,EAP,numQuads=30,scoreMissing=TRUE) {
      z <- matrix(data[2:length(data)],ncol=1)
      pars <- cbind(pars,z)
      colnames(pars)[ncol(pars)] <- "z"
      a <- pars
      for (k in 0:(max(numCats)-2)) {
        pars[pars[,"z"]==k,thresh.labels[(k+2):(max(numCats))]] <- 0 
      }
      EAP.temp <- numeric(3*numQuads)
      EAP.temp <- apply(quads,1,quad.fun,a,numCats,scoreMissing) 
      EAP[,1:3] <- t(EAP.temp)      
      theta_j <- sum(EAP[,"EAP_num"])/sum(EAP[,"EAP_den"])      
      EAP <- cbind(EAP,quads)
      colnames(EAP)[(ncol(EAP)-1):ncol(EAP)] <- c("point","density")        
      EAP[,"std_num"] <- ((EAP[,"point"] - theta_j)^2)*EAP[,"like"]*EAP[,"density"]      
      EAP[,"std_den"] <- EAP[,"like"]*EAP[,"density"]     
      std_j <- sqrt(sum(EAP[,"std_num"])/sum(EAP[,"std_den"]))     
      estimates_j <- cbind(data["SUBJECT"],theta_j,std_j)      
      return(estimates_j)
    }
    
    # Calculate estimates and store in J x 3 matrix
    
    estimates.temp <- numeric(3*nrow(data))
    estimates.temp <- apply(data,1,person.score,pars,numCats,EAP,numQuads,scoreMissing) 
    estimates[,1:3] <- t(estimates.temp)
    
    # Plot distribution of EAP estimates
    
    if (plotTheta) {
      hist(estimates[,"Theta_j"], main = paste("Distribution of EAP Estimates"), 
           xlab = "Theta", ylab = "Frequency")
    }
    
    # Save file of subject numbers, EAP estimates, and associated posterior standard deviations in working directory
    
    if (!is.null(outFile)) {
      write.table(estimates,outFile,sep="\t",row.names=FALSE,quote=FALSE)
    }
    
    # Return a J x 3 matrix of subject numbers, EAP estimates and associated posterior standard deviations
    
    estimates <- round(estimates,8)
    return(estimates)
  }
