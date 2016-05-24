
"summary.SimpleTable" <- function(object, estimand=c("ATE", "ATT", "ATC",
                                       "RR", "RRT", "RRC",
                                       "logRR", "logRRT", "logRRC"),
                                  percent=95, ...){
  S <- object
  
  estimand <- match.arg(estimand)

  S.estimand.pf <- paste("S$", estimand, ".pf", sep="")
  S.estimand.sens <- paste("S$", estimand, ".sens", sep="")
  S.estimand.prior <- paste("S$", estimand, ".prior", sep="")

  stuff.pf <- eval(parse(text=S.estimand.pf))
  stuff.sens <- eval(parse(text=S.estimand.sens))
  stuff.prior <- eval(parse(text=S.estimand.prior))
  
  S.estimand.min <- paste("S$", estimand, ".min", sep="")
  S.estimand.max <- paste("S$", estimand, ".max", sep="")

  stuff.min <- eval(parse(text=S.estimand.min))
  stuff.max <- eval(parse(text=S.estimand.max))
  stuff.bounds <- paste("[", round(stuff.min, 5), ", ",
                        round(stuff.max, 5), "]", sep="")
  
  if (estimand %in% c("ATE", "ATT", "ATC")){
    fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                     ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
    fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                       renorm=TRUE,
                       ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
    fit.prior <- locfit(~lp(stuff.prior, nn=.3, deg=3), renorm=TRUE,
                        ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
    
    x <- seq(from=-1+1e-6, to=1-1e-6, by=.0005)
    
    
    y.pf <- predict(fit.pf, newdata=x)
    y.sens <- predict(fit.sens, newdata=x)
    y.prior <- predict(fit.prior, newdata=x)

    hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)
    hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    hdr.prior  <- hdr(den=list(x=x, y=y.prior), all.modes=FALSE, prob=percent)

    
    mode.pf <- hdr.pf$mode
    mode.sens <- hdr.sens$mode
    mode.prior <- hdr.prior$mode
    
    hdr.pf <- hdrStringList(hdr.pf, x, y.pf)
    hdr.sens <- hdrStringList(hdr.sens, x, y.sens)
    hdr.prior <- hdrStringList(hdr.prior, x, y.prior)


    printSummary(mean.pf=mean(stuff.pf), mode.pf=mode.pf,
                 sd.pf=sd(stuff.pf), hpd.pf=hdr.pf,
                 mean.sens=mean(stuff.sens), mode.sens=mode.sens,
                 sd.sens=sd(stuff.sens), hpd.sens=hdr.sens,
                 mean.prior=mean(stuff.prior), mode.prior=mode.prior,
                 sd.prior=sd(stuff.prior), hpd.prior=hdr.prior,
                 bounds=stuff.bounds, estimand=estimand, hpd.percent=percent)
  }
  if (estimand %in% c("RR", "RRT", "RRC")){
    fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                     ev=rbox(cut=.3), xlim=c(0, 1e10), maxk=500)
    fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                       renorm=TRUE,
                       ev=rbox(cut=.3), xlim=c(0, 1e10), maxk=500)
   
    x <- seq(from=min(c(fit.pf$box, fit.sens$box)),
             to=max(c(quantile(stuff.pf, .999),
               quantile(stuff.sens, .999))),
             length.out=4001)
    
    y.pf <- predict(fit.pf, newdata=x)
    y.sens <- predict(fit.sens, newdata=x)

    hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)
    hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    
    mode.pf <- hdr.pf$mode
    mode.sens <- hdr.sens$mode
    
    hdr.pf <- hdrStringList(hdr.pf, x, y.pf)
    hdr.sens <- hdrStringList(hdr.sens, x, y.sens)


    printSummary(mean.pf=mean(stuff.pf), mode.pf=mode.pf,
                 sd.pf=sd(stuff.pf), hpd.pf=hdr.pf,
                 mean.sens=mean(stuff.sens), mode.sens=mode.sens,
                 sd.sens=sd(stuff.sens), hpd.sens=hdr.sens,
                 mean.prior=NULL, mode.prior=NULL,
                 sd.prior=NULL, hpd.prior=NULL,
                 bounds=stuff.bounds, estimand=estimand, hpd.percent=percent)  
  }
  if (estimand %in% c("logRR", "logRRT", "logRRC")){
    fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                     ev=rbox(cut=.3), maxk=500)
    fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                       renorm=TRUE,
                       ev=rbox(cut=.3), maxk=500)

    x <- seq(from=min(c(stuff.pf, stuff.sens)),
             to=max(c(quantile(stuff.pf), .999),
               quantile(stuff.sens, .999)),
             length.out=4001)
    
    y.pf <- predict(fit.pf, newdata=x)
    y.sens <- predict(fit.sens, newdata=x)

    hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)
    hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    
    mode.pf <- hdr.pf$mode
    mode.sens <- hdr.sens$mode
    
    hdr.pf <- hdrStringList(hdr.pf, x, y.pf)
    hdr.sens <- hdrStringList(hdr.sens, x, y.sens)


    printSummary(mean.pf=mean(stuff.pf), mode.pf=mode.pf,
                 sd.pf=sd(stuff.pf), hpd.pf=hdr.pf,
                 mean.sens=mean(stuff.sens), mode.sens=mode.sens,
                 sd.sens=sd(stuff.sens), hpd.sens=hdr.sens,
                 mean.prior=NULL, mode.prior=NULL,
                 sd.prior=NULL, hpd.prior=NULL,
                 bounds=stuff.bounds, estimand=estimand, hpd.percent=percent)  
        
  }
  
  


} ## end summary.SimpleTable





  
## function to take hdr output and return a nicely formated list of strings
## with the hdr regions
"hdrStringList" <- function(H,x,y){

  H <- hdrCheckInterval(H, x, y)
    
  n.hdr <- ncol(H$hdr)/2

  hdrlist <- vector(mode="list", length=max(n.hdr, 1))
  
  if (n.hdr >= 1){
    hdrlist[[1]] <- paste("[", round(H$hdr[,1], 3), ", ",
                        round(H$hdr[,2], 3), "]", sep="")
  }
  if (n.hdr > 1){
    for (count in 2:n.hdr){
      hdrlist[[count]] <- paste("U [", round(H$hdr[,count*2-1], 3), ", ",
                                round(H$hdr[,count*2], 3), "]", sep="")
    }
  }

  return(hdrlist)
  
}


## cleans up H$hdr so that endpoints are dealt with properly
"hdrCheckInterval" <- function(H, x, y){

  if (y[which.min(x)] > H$falpha){
    H$hdr <- matrix(c(min(x), H$hdr[1,]), 1, ncol(H$hdr)+1) 
  }
  if (y[which.max(x)] > H$falpha){
    H$hdr <- matrix(c(H$hdr[1,], max(x)), 1, ncol(H$hdr)+1)
  }

  return(H)

}




## function to print the summary quantities
"printSummary" <- function(mean.pf, mode.pf, sd.pf, hpd.pf,
                           mean.sens, mode.sens, sd.sens, hpd.sens,
                           mean.prior, mode.prior, sd.prior, hpd.prior,
                           bounds, estimand, hpd.percent){

  n.hpd.pf <- length(hpd.pf)
  n.hpd.sens <- length(hpd.sens)
  n.hpd.prior <- length(hpd.prior)
  n.hpd.max <- max(n.hpd.pf, n.hpd.sens, n.hpd.prior)

  if (estimand == "ATE"){
    supplement = " (Average Treatment Effect)"
  }
  if (estimand == "ATT"){
    supplement = " (Average Treatment Effect Within the Treated)"
  }
  if (estimand == "ATC"){
    supplement = " (Average Treatment Effect Within Controls)"
  }
  if (estimand == "RR"){
    supplement = " (Relative Risk)"
  }
  if (estimand == "RRT"){
    supplement = " (Relative Risk Within the Treated)"
  }
  if (estimand == "RRC"){
    supplement = " (Relative Risk Within Controls)"
  }
  if (estimand == "logRR"){
    supplement = " (Log Relative Risk)"
  }
  if (estimand == "logRRT"){
    supplement = " (Log Relative Risk Within the Treated)"
  }
  if (estimand == "logRRC"){
    supplement = " (Log Relative Risk Within Controls)"
  }
  


  
  if (estimand %in% c("ATE", "ATT", "ATC")){
    cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    cat("Estimand:", estimand, supplement, "\n")
    cat("-------------------------------------------------------------------------------\n")
    cat("                  Prima Facie          Sensitivity Analysis          Prior\n")  
    cat("-------------------------------------------------------------------------------\n")
    nchar.pf <- nchar(as.character(round(mode.pf, 5)))
    nchar.sens <- nchar(as.character(round(mode.sens, 5)))
    nchar.prior <- nchar(as.character(round(mode.prior, 5)))
    mode.string <- paste("Mode:", spaces(24-5-nchar.pf), round(mode.pf, 5),
                         spaces(49-24-nchar.sens), round(mode.sens, 5),
                         spaces(65-45-nchar.prior), round(mode.prior, 5),
                         "\n")
    cat(mode.string)
    
    
    nchar.pf <- nchar(as.character(round(mean.pf, 5)))
    nchar.sens <- nchar(as.character(round(mean.sens, 5)))
    nchar.prior <- nchar(as.character(round(mean.prior, 5)))
    mean.string <- paste("Mean:", spaces(24-5-nchar.pf), round(mean.pf, 5),
                         spaces(49-24-nchar.sens), round(mean.sens, 5),
                         spaces(65-45-nchar.prior), round(mean.prior, 5),
                         "\n")
    cat(mean.string)
    
    nchar.pf <- nchar(as.character(round(sd.pf, 5)))
    nchar.sens <- nchar(as.character(round(sd.sens, 5)))
    nchar.prior <- nchar(as.character(round(sd.prior, 5)))
    sd.string <- paste("Std. Dev.:", spaces(24-10-nchar.pf), round(sd.pf, 5),
                       spaces(49-24-nchar.sens), round(sd.sens, 5),
                       spaces(65-45-nchar.prior), round(sd.prior, 5),
                       "\n")
    cat(sd.string)
    
    
    nchar.pf <- nchar(hpd.pf[[1]])
    nchar.sens <- nchar(hpd.sens[[1]])
    nchar.prior <- nchar(hpd.prior[[1]])
    hpd.string <- paste(round(hpd.percent,2), "% HPD Region:",
                        spaces(32-15-nchar.pf), hpd.pf[[1]],
                        spaces(57-32-nchar.sens),
                        hpd.sens[[1]],
                        spaces(79-57-nchar.prior),
                        hpd.prior[[1]], "\n", sep="")
    cat(hpd.string)
    if (n.hpd.max > 1){
      for (i in 2:n.hpd.max){
        
        if (length(hpd.pf) < i){
          hpd.pf[[i]] <- "  "
        }
        if (length(hpd.sens) < i){
          hpd.sens[[i]] <- "  "
        }
        if (length(hpd.prior) < i){
          hpd.prior[[i]] <- "  "
        }
        
        nchar.pf <- nchar(hpd.pf[[i]])
        nchar.sens <- nchar(hpd.sens[[i]])
        nchar.prior <- nchar(hpd.prior[[i]])
        hpd.string <- paste(spaces(32-nchar.pf), hpd.pf[[i]],
                            spaces(57-32-nchar.sens),
                            hpd.sens[[i]],
                            spaces(79-57-nchar.prior),
                            hpd.prior[[i]], "\n", sep="")
        cat(hpd.string)    
      }
    }
    
    cat("-------------------------------------------------------------------------------\n")
    
    cat("Large Sample Nonparametric Bounds:", bounds)
    cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")
  }
  else{ ## a relative risk like quantity is the estimand
    cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n")
    cat("Estimand:", estimand, supplement, "\n")
    cat("-------------------------------------------------------------------------------\n")
    cat("                  Prima Facie          Sensitivity Analysis         \n")  
    cat("-------------------------------------------------------------------------------\n")
    nchar.pf <- nchar(as.character(round(mode.pf, 5)))
    nchar.sens <- nchar(as.character(round(mode.sens, 5)))
    mode.string <- paste("Mode:", spaces(24-5-nchar.pf), round(mode.pf, 5),
                         spaces(49-24-nchar.sens), round(mode.sens, 5),
                         "\n")
    cat(mode.string)
    
    
    nchar.pf <- nchar(as.character(round(mean.pf, 5)))
    nchar.sens <- nchar(as.character(round(mean.sens, 5)))
    mean.string <- paste("Mean:", spaces(24-5-nchar.pf), round(mean.pf, 5),
                         spaces(49-24-nchar.sens), round(mean.sens, 5),
                         "\n")
    cat(mean.string)
    
    nchar.pf <- nchar(as.character(round(sd.pf, 5)))
    nchar.sens <- nchar(as.character(round(sd.sens, 5)))
    sd.string <- paste("Std. Dev.:", spaces(24-10-nchar.pf), round(sd.pf, 5),
                       spaces(49-24-nchar.sens), round(sd.sens, 5),
                       "\n")
    cat(sd.string)
    
    
    nchar.pf <- nchar(hpd.pf[[1]])
    nchar.sens <- nchar(hpd.sens[[1]])
    hpd.string <- paste(round(hpd.percent,2), "% HPD Region:",
                        spaces(32-15-nchar.pf), hpd.pf[[1]],
                        spaces(57-32-nchar.sens),
                        hpd.sens[[1]],
                        "\n", sep="")
    cat(hpd.string)
    if (n.hpd.max > 1){
      for (i in 2:n.hpd.max){
        
        if (length(hpd.pf) < i){
          hpd.pf[[i]] <- "  "
        }
        if (length(hpd.sens) < i){
          hpd.sens[[i]] <- "  "
        }
        
        nchar.pf <- nchar(hpd.pf[[i]])
        nchar.sens <- nchar(hpd.sens[[i]])
        hpd.string <- paste(spaces(32-nchar.pf), hpd.pf[[i]],
                            spaces(57-32-nchar.sens),
                            hpd.sens[[i]], "\n", sep="")
        cat(hpd.string)    
      }
    }
    
    cat("-------------------------------------------------------------------------------\n")
    
    cat("Large Sample Nonparametric Bounds:", bounds)
    cat("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n")
  }

  
}

      
      
## creates a string with n blank spaces
"spaces" <- function(n){
  output <- NULL
  for (i in 1:n){
    output <- paste(output, " ", sep="")
  }
  return(output)
}



