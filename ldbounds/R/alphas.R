"alphas" <-
function(iuse,asf,alpha,phi,side,t){
   tol <- 10^(-13)
   if (iuse==1){
      pe <- 2*(1-pnorm(qnorm(1-(alpha/side)/2)/sqrt(t)))
      spend <- "O'Brien-Fleming"
    }
   else if (iuse==2){
      pe <- (alpha/side)*log(1+(exp(1)-1)*t)
      spend <- "Pocock"
    }
   else if (iuse==3){
      pe <- (alpha/side)*t^phi
      spend <- "Power Family: alpha * t^phi"
    }
   else if (iuse==4){
      pe <- (alpha/side)*(1-exp(phi*t))/(1-exp(-phi))
      spend <- "Hwang-Shih-DeCani Family"
    }
   else if (iuse==5){
      checkalpha <- asf((c(1:51)-1)/50)
      if ({checkalpha[1]!=0}|{checkalpha[51]!=1}|{sum(checkalpha[-1]-checkalpha[-51]<=0)>=1}){
         stop("Alpha spending function error")
       }
      spend <- "User-specified spending function"
      pe <- (alpha/side)*asf(t)
    }
   else stop("Must choose 1, 2, 3, 4, or 5 as spending function.")
   pe <- side*pe
   pd <- pe-c(0,pe[-length(pe)])
   if (sum(as.integer({pd<0}|{pd>1})) >= 1){
      warning("Spending function error")
      pd <- min(1,pd)
      pd <- max(0,pd)
    }
   for (j in 1:length(pd)){
      if (pd[j] < tol){
         warning("Type I error spent too small for analysis #",j,"\n",
                 "Zero used as approximation for ",pd[j])
         pd[j] <- 0
       }
    }
   ans <- list(pe=pe,pd=pd,spend=spend)
   return(ans)
 }

