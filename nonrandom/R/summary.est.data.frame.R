summary.est.data.frame <- function(object,
                                   ...)
{
  if( object$family == "binomial" ){
    meas <- "odds ratios"
    col.n <- c("or","SE[log[or]]","[95%-CI[or]]")
  }else{
    meas <- "effects"
    col.n <- c("effect","SE[effect]","[95%-CI[effect]]")
  }
  
  if ( !is.list(object$lr.estimation) ){
    lr.eff.c <- lr.se.c <- 
      lr.eff.m <- lr.se.m <- 
        lr.ci.c <- lr.ci.m <- NULL
    
  }else{
    if( object$family == "binomial" ){
      lr.eff.c <- round(object$lr.estimation$effect,3)
      lr.se.c  <- round(object$lr.estimation$se,4)
      lr.eff.m <- round(object$lr.estimation$effect.marg,3)
      lr.se.m  <- round(object$lr.estimation$se.marg,4)

      lr.ci.c <- round(c(exp(log(lr.eff.c) - qnorm(0.975)*lr.se.c),
                         exp(log(lr.eff.c) + qnorm(0.975)*lr.se.c)),3)

      lr.ci.m <- round(c(exp(log(lr.eff.m) - qnorm(0.975)*lr.se.m),
                         exp(log(lr.eff.m) + qnorm(0.975)*lr.se.m)),3)
    }else{
      lr.eff.c <- round(object$lr.estimation$effect,3)
      lr.se.c  <- round(object$lr.estimation$se,4)
      
      lr.ci.c <- round(c(lr.eff.c - qnorm(0.975)*lr.se.c,
                         lr.eff.c + qnorm(0.975)*lr.se.c),3)
    }
  }
    
  if ( !is.list(object$ps.estimation$adj) ){    
    ps.adj.eff <-
      ps.adj.se <-
        ps.adj.ci <- NULL
    
    if (!is.null(object$stratum.index)){ ps.adj.str <- NULL }
    
  }else{
    ps.adj.eff <- round(object$ps.estimation$adj$effect,3)
    ps.adj.se  <- round(object$ps.estimation$adj$se,4)

    if (!is.null(object$stratum.index)){
      ps.adj.str <- round(object$ps.estimation$adj$effect.str,3)
    }else{
      ps.adj.str <- NULL      
    }
    
    if( object$family == "binomial" ){
      ps.adj.ci <- round(c(exp(log(ps.adj.eff) - qnorm(0.975)*ps.adj.se),
                           exp(log(ps.adj.eff) + qnorm(0.975)*ps.adj.se)),3)     
    }else{
      ps.adj.ci <- round(c(ps.adj.eff - qnorm(0.975)*ps.adj.se,
                           ps.adj.eff + qnorm(0.975)*ps.adj.se),3)     
    }
  }

  
  if ( object$family=="binomial" ){

    crude <- round(object$ps.estimation$crude$effect,3)
    crude.se <- round(object$ps.estimation$crude$se,4)
    crude.ci <- round(c(exp(log(crude)-qnorm(0.975)*crude.se),
                        exp(log(crude)+qnorm(0.975)*crude.se)),3)

    if (!is.null(object$stratum.index)){ ## binomial, stratified
    
      rr <- round(object$ps.estimation$unadj$effect,3)
      rr.se <- round(object$ps.estimation$unadj$se,4)  
      mh <- round(object$ps.estimation$unadj$effect.mh,3)
      mh.se <- round(object$ps.estimation$unadj$se.mh,4)
      
      rr.ci <- round(c(exp(log(rr)-qnorm(0.975)*rr.se),
                       exp(log(rr)+qnorm(0.975)*rr.se)),3)
      mh.ci <- round(c(exp(log(mh)-qnorm(0.975)*mh.se),
                       exp(log(mh)+qnorm(0.975)*mh.se)),3)
      
      str <-
        data.frame(rbind(rep("-----",times=length(levels(object$stratum.index))),
                         round(object$ps.estimation$unadj$p0.str,2),
                         round(object$ps.estimation$unadj$p1.str,2),
                         round(object$ps.estimation$unadj$odds.str,2)),
                   row.names=c("",
                     "  outcome rates 'p0'  ",
                     "  outcome rates 'p1'  ",
                     "  odds ratio"))      
      colnames(str) <-
        c(paste("S", seq(1:length(levels(object$stratum.index))),sep=""))

    }else{ ## binomial, matched

      crude <- round(object$ps.estimation$crude$effect,3)
      crude.se <- round(object$ps.estimation$crude$se,4)
      crude.ci <- round(c(crude-qnorm(0.975)*crude.se,
                          crude+qnorm(0.975)*crude.se),3)
      
      diff <- round(object$ps.estimation$unadj$effect,3)
      diff.se <- round(object$ps.estimation$unadj$se,4)  
      
      diff.ci <- round(c(exp(log(diff)-qnorm(0.975)*diff.se),
                         exp(log(diff)+qnorm(0.975)*diff.se)),3)  
      
    }
    
  }else{ ## gaussian
    
    crude <- round(object$ps.estimation$crude$effect,3)
    crude.se <- round(object$ps.estimation$crude$se,4)
    crude.ci <- round(c(crude-qnorm(0.975)*crude.se,
                        crude+qnorm(0.975)*crude.se),3)
    
    diff <- round(object$ps.estimation$unadj$effect,3)
    diff.se <- round(object$ps.estimation$unadj$se,4)  
    
    diff.ci <- round(c(diff-qnorm(0.975)*diff.se,
                       diff+qnorm(0.975)*diff.se),3)      

    if (!is.null(object$stratum.index)){ ## gaussian, stratified
      str <-
        data.frame(rbind(rep("-----",times=length(levels(object$stratum.index))),
                         round(object$ps.estimation$unadj$effect.str,3)),
                   row.names=c("","Unadjusted effect"))
      
      colnames(str) <-
        c(paste("S", seq(1:length(levels(object$stratum.index))),sep=""))
    }
  }
  
  
  if ( object$family=="binomial" ){ ## binomial
    
    if (!is.null(object$stratum.index)){ ## binomial, stratified
    
      eff <- c(" -----",
               paste("",crude,sep=""),
               "",
               paste("",rr,sep=""),
               paste("",mh,sep=""),
               paste("",ps.adj.eff,sep=""),
               "",
               paste("",lr.eff.c,sep=""),
               paste("",lr.eff.m,sep=""),
               "")
      
      se <- c(" -----------",
              paste("",crude.se,"", sep=""),
              "",
              paste("",rr.se,"", sep=""),
              paste("",mh.se,"", sep=""),
              paste("",ps.adj.se,"", sep=""),
              "",
              paste("",lr.se.c,"", sep=""),
              paste("",lr.se.m,"", sep=""),
              "")
      
      ci <- c(" ------------",
              paste("[",crude.ci[1],",",crude.ci[2], "]",sep=""),
              "",
              paste("[",rr.ci[1],",",rr.ci[2],"]",sep=""),
              paste("[",mh.ci[1],",",mh.ci[2],"]",sep=""),
              paste("[",ps.adj.ci[1],",",ps.adj.ci[2],"]",sep=""),
              "",
              paste("[",lr.ci.c[1],",",lr.ci.c[2],"]",sep=""),
              paste("[",lr.ci.m[1],",",lr.ci.m[2],"]",sep=""),
              "")
      
    df <- data.frame(cbind(eff, se, ci),
                     row.names=c("",
                       "Crude",
                       "Stratification", "  Outcome rates", "  MH", "  Adjusted",
                       "Regression", "  Conditional", "  Marginal", " "))

    }else{ ## binomial, matched

      eff <- c(" ------",
               paste("",crude,sep=""),
               "",
               paste("",diff,sep=""),
               paste("",ps.adj.eff,sep=""),
               "",
               paste("",lr.eff.c,sep=""),
               paste("",lr.eff.m,sep=""),
               "")
      
      se <- c(" ----------",
              paste("",crude.se,"", sep=""),
              "",
              paste("",diff.se,"", sep=""),
              paste("",ps.adj.se,"", sep=""),
              "",
              paste("",lr.se.c,"", sep=""),
              paste("",lr.se.m,"", sep=""),
              "")
      
      ci <- c(" ----------------",
              paste("[",crude.ci[1],",",crude.ci[2], "]",sep=""),
              "",
              paste("[",diff.ci[1],",",diff.ci[2],"]",sep=""),
              paste("[",ps.adj.ci[1],",",ps.adj.ci[2],"]",sep=""),
              "",
              paste("[",lr.ci.c[1],",",lr.ci.c[2],"]",sep=""),
              paste("[",lr.ci.m[1],",",lr.ci.m[2],"]",sep=""),
              "")
      
      df <- data.frame(cbind(eff, se, ci),
                       row.names=c("",
                         "Crude",
                         "Matching", " Unadjusted", " Adjusted",
                         "Regression", "  Conditional", "  Marginal", " "))
      
    }
    colnames(df) <- col.n
    
  }else{ ## gaussian
    
    eff <- c(" ------",
             paste("",crude,sep=""),
             "",
             paste("",diff,sep=""),
             paste("",ps.adj.eff,sep=""),
             paste("",lr.eff.c,sep=""),
             "")
    
    se <- c(" ----------",
            paste("",crude.se,"", sep=""),
            "",
            paste("",diff.se,"", sep=""),
            paste("",ps.adj.se,"", sep=""),
            paste("",lr.se.c,"", sep=""),
            "")
    
    ci <- c(" ----------------",
            paste("[",crude.ci[1],",",crude.ci[2], "]",sep=""),
            "",
            paste("[",diff.ci[1],",",diff.ci[2],"]",sep=""),
            paste("[",ps.adj.ci[1],",",ps.adj.ci[2],"]",sep=""),
            paste("[",lr.ci.c[1],",",lr.ci.c[2],"]",sep=""),
            "")

    if (!is.null(object$stratum.index)){
    
      df <- data.frame(cbind(eff, se, ci),
                       row.names=c("",
                         "Crude",
                         "Stratification", " Unadjusted", " Adjusted",
                         "Regression", " "))
    }else{
      df <- data.frame(cbind(eff, se, ci),
                       row.names=c("",
                         "Crude",
                         "Matching", " Unadjusted", " Adjusted",
                         "Regression", " "))
    }
    colnames(df) <- col.n

  }

  if (!is.null(object$stratum.index)){ ## stratified
    
    sum.eff <- list(treat   = object$name.treat,
                    outcome = object$name.resp,
                    measure = meas,
                    eff.tab = df,
                    weights = round(object$ps.estimation$weights.str,2),                 
                    stratum = str,
                    adjust  = ps.adj.str,
                    type    = "str")

  }else{ ## matched
    sum.eff <- list(treat   = object$name.treat,
                    outcome = object$name.resp,
                    measure = meas,
                    eff.tab = df,
                    type    = "mat")
  }

  class(sum.eff) <- "summary.est.data.frame"
  
  sum.eff
 
}
