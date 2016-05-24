mpci<- function(index = c("shah", "taam", "pan", "wangw", "wang", "xeke", "wangw"),x ,
                LSL, USL, Target, npc, alpha, Method, perc, graphic, xlab = "Var 1", ylab = "Var 2", ...){   

     index <- match.arg(index)

     ###Variables
     p <- ncol(x) # number of quality characteristics
     m <- nrow(x) # sample number
     Xmv <- colMeans(x) # mean vector
     S <- cov(x) # covariance matrix
	 rr <- cor(x) # correlation matrix
     Sinv <- solve(S)
     s3<-matrix(0,p-1,p-1)
	 ###
     if (missing(index)){ #used as default index
         index<-"shah"
	 }
	 
     if(missing(Target)) { # Estimating Target 
         Target <- LSL + (USL - LSL) / 2
     }

     ###Plotting 
     if(missing(alpha)) { # Setting alpha 
         alpha <- 0.0027
     } 
     if(missing(graphic)) { 
         graphic <- FALSE
     }      
	 
	 
	 if(p == 2 && graphic == TRUE ){
         LPL <- UPL <- matrix(0,nrow = p,ncol = 1)

         for(i in 1:p){
             s3 <- matrix(Sinv[-i,-i])
             LPL[i,1] <- Xmv[i] - sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
             UPL[i,1] <- Xmv[i] + sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
             }
         
		 plot(x,xlim=c(min(LSL[1],LPL[1]),max(USL[1],UPL[1]) ), ylim=c(min(LSL[2],LPL[2]), max(USL[2],UPL[2])),xlab = xlab, ylab = ylab)
         rect( xleft <- LSL[1], xright <- LSL[2], yleft <- USL[1], yright <- USL[2], border = 2)
         points(Target[1],Target[2], pch=3, col=2, cex=0.7)
		 points(Xmv[1],Xmv[2], pch=16, col=1, cex=0.7)
		 
		 ### Making the confidence ellipsoid
	     Ue<-eigen(S, symmetric=TRUE)$vectors    # eigenvectors
         DDe<-eigen(S, symmetric=TRUE)$values

         angle <- seq(0, 2 * pi, length.out = 200)
         ch <- cbind(sqrt(qchisq(1 - alpha,2)) * cos(angle),sqrt(qchisq(1 - alpha,2)) * sin(angle))
         lines(t(Xmv - ((Ue %*% diag(sqrt(DDe))) %*% t(ch))),type = "l")
		 
		 if(index=="shah"){
			rect( xleft <- LPL[1], xright <- LPL[2], yleft <- UPL[1], yright <- UPL[2],lty = 2, border = 4) 
			cur <- "Modified Process Region"
		 }
		 
		 #### Making the largest ellipsoid within the tolerance region centered at Target
         if(index=="taam"){
		 hi <- (USL[1] - LSL[1]) / 2
         lo <- (USL[2] - LSL[2]) / 2
         Xm <- colMeans(rbind(LSL, USL))
         nn <- 201

         d2 <- (hi-lo) * (hi+lo)                 
         ang <- 2 * pi * seq(0,1, len = 201)
         r <- lo * hi / sqrt(lo^2 + d2 * sin(ang)^2)
         pro <- r * cbind(cos(ang), sin(ang))
         al <- alpha * pi/180
         xx <- pro %*% rbind(c(cos(al), sin(al)), c(-sin(al), cos(al))) + cbind(rep(Xm[1], nn),rep(Xm[2], nn))
	     points(xx, type="l", lty = 2, col=4)
		 cur <- "Modified Tolerance Region"
		}
         
	 
		 ### Making ellipsoid from Pan index	
	     if(index=="pan"){
			A <- matrix(0,p,p)
		
			 for (i in 1 :nrow(rr)){
				 for (j in 1 :ncol(rr)){
					 A[i,j] <- rr[i,j] * ((USL[i] - LSL[i])/(2 * sqrt(qchisq((1 - alpha),p)))) * ((USL[j] - LSL[j])/(2 * sqrt(qchisq((1 - alpha), p))))
		   		  }         
			  }
		 	 
		 Ue<-eigen(A, symmetric=TRUE)$vectors    # eigenvectors
         DDe<-eigen(A, symmetric=TRUE)$values

         angle <- seq(0, 2 * pi, length.out = 200)
         ch <- cbind(sqrt(qchisq(1 - alpha,2)) * cos(angle),sqrt(qchisq(1 - alpha,2)) * sin(angle))
         lines(t(Target - ((Ue %*% diag(sqrt(DDe))) %*% t(ch))),type = "l",lty = 2, col=4)
		 cur <- "Modified Tolerance Region"
		 }
		 
		
	     ###
	     
         legend(locator(1),cex=0.85,c("Process Region","Tolerance Region",paste(cur),"Target","Process Mean"),lty=c(1,1,2,NA,NA),pch=c(NA,NA,NA,3,16),col=c(1,2,4,2,1))

         #print(list("CpM index of Shahriari et al. (1995) is the ratio of the Tolerance Region and the Modified Process Region"))
         #print(list("MCpm index of Taam et al. (1993) is the ratio of the ellipsoids: Modified Tolerance Region and the Process Region "))
 
     }	 
	 
     if (index=="shah") {
	   
         if(missing(alpha)) { # Setting alpha 
             alpha <- 0.0027
         } 

         ###First Component
         LPL <- UPL <- matrix(0,nrow = p,ncol = 1)

         for(i in 1:p){
             if (p == 2) {s3 <- matrix(Sinv[-i,-i])}
			 else {s3 <- Sinv[-i,-i]}
             LPL[i,1] <- Xmv[i] - sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
             UPL[i,1] <- Xmv[i] + sqrt((det(s3) * (qchisq((1 - alpha),df = p))) / det(Sinv))
         }

         spec <- cbind(LSL,USL)
	     proc <- cbind(LPL,UPL)
	  
	     CpM <- (prod(USL - LSL) / prod(UPL - LPL)) ^ (1 / p)  
	 
         ###Second Component
         PV <- 1 - pf((t(Target - Xmv) %*% solve(S) %*% (Target - Xmv)) * ((m * (m - p))/(p * (m - 1))), p, m - p)

         ###Third Component
	 
	     if( any(spec[,1] > proc[,1]) || any(spec[,2] < proc[,2])) LI <- 0 else LI <- 1

         ###Vector
         return(list ("Shahriari et al. (1995) Multivariate Capability Vector","CpM" = CpM,"PV" = PV,"LI" = LI))
     } 
 
 
     if (index=="taam") {	 
	 
         if(missing(alpha)) { # Setting alpha 
             alpha <- 0.0027
         }

         LPL <- UPL <- matrix(0,nrow = p,ncol = 1)

         VTR <- 2 * (prod((USL - LSL) / 2)) * pi ^ (p / 2)/(p * gamma(p / 2)) # Vol. Tolerance Region   
         VE <- det(S) ^ 0.5 * ((pi * qchisq(1 - alpha, p)) ^ (p / 2))/gamma(p / 2 + 1) # Vol. Estimated 99.73% Process Region
         Cp <- VTR / VE
         D <- (1 + m / (m - 1) * (t(Target - Xmv) %*% solve(S) %*% (Target - Xmv))) ^ 0.5
         MCpm <- Cp / D

         return(list ("Taam et al. (1993) Multivariate Capability Index (MCpm)","MCpm"=MCpm))
     }
	 
     if (index=="pan") {
	     A <- matrix(0,p,p)
	 	 for (i in 1 :nrow(rr)){
			for (j in 1 :ncol(rr)){
             A[i,j] <- rr[i,j]*((USL[i] - LSL[i]) / (2 * sqrt(qchisq((1 - alpha), p)))) * ((USL[j]-LSL[j])/(2 * sqrt(qchisq((1 - alpha),p))))
			}         
		}
	 	 
	     VTR <- det(A) ^ 0.5 * ((pi * qchisq(1 - alpha, p)) ^ (p / 2))/gamma(p / 2 + 1)# Vol. Tolerance Region   
         VE <- det(S) ^ 0.5 * ((pi * qchisq(1 - alpha, p)) ^ (p / 2))/gamma(p / 2 + 1) # Vol. Estimated 99.73% Process Region
         NMCp <- VTR / VE
         D <- (1 + m / (m - 1) * (t(Target - Xmv) %*% solve(S) %*% (Target - Xmv))) ^ 0.5
         NMCpm <- NMCp / D
		 
		 return(list ("Pan and Lee (2010) Multivariate Capability Index (NMCpm)","NMCpm"=NMCpm))
	 
	 }
     if (index == "wang" || index == "xeke"|| index == "wangw") {	 
	 
         spec <- cbind(LSL,USL) # matrix of specifications
     
         if(missing(perc)){
	         perc<-0.8
	     }
 
         Ue<-eigen(S, symmetric=TRUE)$vectors    # eigenvectors
         DDe<-eigen(S, symmetric=TRUE)$values     # eigenvalues

         if(!missing(npc)) { 
             if(npc<=0 || npc>p || !is.numeric(npc) || npc != as.integer(npc) || length(npc) > 1){
                 stop("Attention: the number of principal components (npc) must be a integer between 1 
                 and the number of quality characteristics")
             }
         }

         if(missing(npc)) {    #number of principal components
             #Modified Algorithm of Rencher,A.C.(2002) Methods of Multivariate Analysis. John Wiley and Sons.
             #12.6 DECIDING HOW MANY COMPONENTS TO RETAIN
  
             if(missing(alpha)) { # Setting alpha 
                 alpha<-0.05
             }
 
             if(!missing(alpha)) {  
                 if(alpha<0 || alpha>1){
                     stop("Attention: the significance level (alpha) must be between 0 and 1")
                 }
             } 

             if(missing(Method)) { # Method to select the number of principal components
                 Method<-"Percentage"   
             }
     
             if(!missing(Method)) { 
                 if((Method<=0 || Method>5) && (Method != "Percentage" && Method != "Average" && Method != "Scree" && Method != "Bartlett.test" && Method!="Anderson.test")){
                     stop("Attention: the Method must be a integer between 1 and 5 or one of the followings:
                     Percentage, Average, Scree, Bartlett.test or Anderson.test")
                 }
             }
     
             ###Methods for npc
             if (Method == "Percentage" || Method == 1) {
                 npc <- which(cumsum(DDe) / sum(DDe) > perc)[1]
             } 

             ###
             if (Method == "Average" || Method == 2){
		         npc <- length(which(DDe > mean(DDe)))
		     }
 
             ###
             if (Method == "Scree" || Method == 3){ 
                 plot(DDe,main = "Scree graph for eigenvalue",xlab="Eigenvalue number",ylab="Eigenvalue size",type="o",pch=1,lty=1)
                 cat("\n","Enter the number of principal component(npc) according to the scree graph:: ","\n")
                 npc <- scan(n = 1) 
         
                 if(npc <= 0 || npc > p || !is.numeric(npc) || npc != as.integer(npc) || length(npc) > 1){
                     stop("Attention: the number of principal components (npc) must be a integer between 1 
                     and the number of quality characteristics")
                 }
             }

             ###
             if (Method == "Bartlett.test" || Method==4){
                chi.t <- chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value
         
		         for (i in seq_along(DDe)){
                     DD<-DDe[i:p]
                     chi.p[i]<-(m-(2*p+11)/6)*((p-i+1)*log(mean(DD))-sum(log(DD)))
			         chi.t[i]<- qchisq(1-alpha,((length(DD)-1)*(length(DD)+2)/2)) 
                 }

                 npc <- length(which(chi.p>chi.t))
		 
                 if (npc == 0){
                     stop("There are no difference between principal components according to Bartlett's Test.
                     Please use another method (1,2,3 or 5)")}
       	         }

             ####
             if (Method == "Anderson.test" || Method==5){
     
                 chi.t <- chi.p<-matrix(0,1,(length(DDe))) # (practical) chi-squared value

                 for (i in seq_along(DDe)){
                     DD <- DDe[i:p]
                     chi.p[i] <- (m-1) * length(DD) * log(sum(DD) / length(DD)) - (m-1) * sum(log(DD))
                     chi.t[i] <- qchisq(1 - alpha, ((length(DD) - 1) * (length(DD) + 2) / 2)) 
                 }

                 npc <- length(which(chi.p > chi.t))
                 if (npc == 0){
                     stop("There are no difference between principal components according to Anderson's Test.
                     Please use another method (1,2,3 or 4)")
				 }
             }

         }
     
	 
         SLpce <- t(Ue) %*% spec  #matrix of spec. limits for each PC (rows:LSLPC, USLPC)
         Xmvpce <- t(Ue) %*% Xmv ##vector of sample means for each PC
         Targetpce <- t(Ue) %*% Target ##vector of target values for each PC
         numCpe <- abs(SLpce[,2] - SLpce[,1])  ##COMPUTING MATRIX OF Cp 
         MatCpe <- numCpe / (6 * sqrt(DDe)) ##this is matrix of the Cppc
         MatCpme <- numCpe / (6 * sqrt(DDe + (Xmvpce - Targetpce) ^ 2)) ##COMPUTING MATRIX OF Cpm 

         ###COMPUTING MATRIX OF Cpk

         ##this is matrix of the UPPER
         Cppcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe)) 
         ##this is matrix of the LOWER
         Cppcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe)) 
         appo1e<-cbind(Cppcle,Cppcue)
         MatCpke<-matrix(NA, nrow=p, ncol=1)

         ###COMPUTING MATRIX OF Cpmk 

         ##this is matrix of the UPPER
         Cpmpcue<-abs((SLpce[,2]-Xmvpce))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
         #this is matrix of the LOWER
         Cpmpcle<-abs((Xmvpce-SLpce[,1]))/(3*sqrt(DDe+(Xmvpce-Targetpce)^2)) 
         appo2e<-cbind(Cpmpcle,Cpmpcue)
         MatCpmke<-matrix(NA, nrow=p, ncol=1)

         ##  MatCpke MatCpmke
         MatCpke<-apply(appo1e,1,min)
	     MatCpmke<-apply(appo2e,1,min)
	 }
	 
	 
	 ###WANG-CHEN
     if (index == "wang") {
	     MCpe <- (prod(MatCpe[1 : npc])) ^ (1 / npc)
         MCpke <- (prod(MatCpke[1 : npc])) ^ (1 / npc)
         MCpme <- (prod(MatCpme[1 : npc])) ^ (1 / npc)
         MCpmke <- (prod(MatCpmke[1 : npc])) ^ (1 / npc)

         return(list ("Wang and Chen (1998) Multivariate Process Capability Indices(PCI) based on PCA",
         "number of principal components"=npc,"MCp"=MCpe,"MCpk"=MCpke,"MCpm"=MCpme,"MCpmk"=MCpmke))
     } 

     ###Xekalaki-Perakis
     if (index == "xeke") {
		
         su_we <- sum(DDe[1 : npc])

         vMXCpe <- MatCpe * DDe
         MXCpe <- (sum(vMXCpe[1 : npc])) / su_we

         vMXCpke <- MatCpke * DDe  
         MXCpke <- (sum(vMXCpke[1 : npc])) / su_we

         vMXCpme <- MatCpme * DDe
         MXCpme <- (sum(vMXCpme[1 : npc])) / su_we

         vMXCpmke <- MatCpmke * DDe
         MXCpmke <- (sum(vMXCpmke[1 : npc])) / su_we

         return(list ("Xekalaki and Perakis (2002) Multivariate Process Capability Indices(PCI) based on PCA",
         "number of principal components"=npc,"MCp"=MXCpe,"MCpk"=MXCpke,"MCpm"=MXCpme,"MCpmk"=MXCpmke))

     }

     ###CH Wang
     if (index == "wangw") {
         ##SUM OF WEIGHTS
         su_we <- sum(DDe[1 : npc])

    
         vMWCpe <- (MatCpe) ^ (DDe)
         MWCpe <- (prod(vMWCpe[1 : npc])) ^ (1 / su_we)

         vMWCpke<-(MatCpke) ^ (DDe)
         MWCpke <- (prod(vMWCpke[1 : npc])) ^ (1 / su_we)

         vMWCpme <- (MatCpme) ^ (DDe)
         MWCpme <- (prod(vMWCpme[1 : npc])) ^ (1 / su_we)

         vMWCpmke <- (MatCpmke) ^ (DDe)
         MWCpmke <- (prod(vMWCpmke[1 : npc])) ^ (1 / su_we)

         return(list ("Wang(2005) Multivariate Process Capability Indices(PCI) based on PCA","number of principal components" = npc,
         "MCp" = MWCpe,"MCpk" = MWCpke,"MCpm" = MWCpme,"MCpmk" = MWCpmke))
     }  
	 
}
