### redist.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Nov 28 2015 (10:30) 
## Version: 
## last-updated: Nov 28 2015 (10:35) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' Calculation of Efron's re-distribution to the right algorithm to obtain the
##' Kaplan-Meier estimate.
##'
#' @param time A numeric vector of event times.
#' @param status The event status vector takes the value \code{1} for observed events and 
#' the value \code{0} for right censored times.  
##' @return Calculations needed to 
##' @seealso prodlim 
##' @examples
##' redist(time=c(.35,0.4,.51,.51,.7,.73),status=c(0,1,1,0,0,1))
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
redist <- function(time,status){
    library(prodlim)
    cat("\nKaplan-Meier estimate via re-distribution to the right algorithm:\n")
    order <- order(time,-status)
    time <- time[order]
    status <- status[order]
    N <- length(time)
    mass <- as.list(rep(1/N,N))
    fractions <- as.list(rep(paste("1/",N,sep=""),N))
    names(mass) <- paste("subject",1:N)
    for (i in 1:N) names(mass[[i]]) <- "own"
    for (i in 1:N) names(fractions[[i]]) <- "own contribution"
    surv <- 1
    for (i in 1:N) {
        cat("\nSubject ",i,":\n---------------------------\nSurvival before = ",round(surv*100,2),"%\n",sep="")
        if (status[i]==0){
            if (i==N){
                cat("Last subject lost to follow-up event free at time = ",time[i],"\n",sep="")
            }
            else{
                cat("No event until time = ",time[i],"\nRe-distribute mass ",signif(sum(mass[[i]]),2)," to remaining ",N-i,ifelse(N-i==1," subject"," subjects"),"\n",sep="")
                for (j in ((i+1):N)){
                    mass[[j]] <- c(mass[[j]],mass[[i]]/(N-i))
                    fractions[[j]] <- c(fractions[[j]],paste(fractions[[i]],"*1/",(N-i),sep=""))
                    names(fractions[[j]])[length(fractions[[j]])-length(mass[[i]])+1] <- paste("from subject ",i,sep="")
                    names(mass[[j]])[length(mass[[j]])] <- paste("from subject ",i,sep="")
                }
            }
            cat("Survival after = ",round(surv*100,2),"%\n",sep="")
        } else{
            cat("Event at time = ",time[i],"\nContribution to Kaplan-Meier estimate:\n\n",sep="")
            contr <- rbind(fractions[[i]],format(mass[[i]],digits=4,nsmall=4))
            rownames(contr) <- c("fractions","decimal")
            contr <- rbind(t(contr),c("sum",format(sum(mass[[i]]),digits=4,nsmall=4)))
            print(contr,quote=FALSE)
            surv.before <- surv
            surv <- surv-sum(mass[[i]])
            cat("\nSurvival after = ",round(100*surv.before,2),"% - (",paste(fractions[[i]],collapse=" + ") ,")",
                "\n               = ",round(100*surv.before,2),"% - ",round(100*sum(mass[[i]]),2) ,"% = ",round(surv*100,2),"%\n",sep="")
        }
    }
    table <- summary(f <- prodlim(Hist(time,status)~1,data=data.frame(time,status)),times=c(0,time),percent=TRUE)
    cat("\nSummary table:\n\n")
    tab <- table$table[,c("time","n.risk","n.event","n.lost","surv")]
    print(tab)
    out <- list(fit=f,table=tab)
    invisible(out)
}



#----------------------------------------------------------------------
### redist.R ends here
