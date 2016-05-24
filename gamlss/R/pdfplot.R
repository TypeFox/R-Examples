#------------------------------------------------------------------------------------------
#                                   pdf.plot
#                  Kalliope Akantziliotou, Mikis Stasinopoulos
#                 last       Thursday, February 6, 2003 at 12:04
#------------------------------------------------------------------------------------------
pdf.plot<- function (obj = NULL, 
                     obs = c(1), 
                  family = NO(),  
                      mu = NULL, 
                   sigma = NULL,  
                      nu = NULL,  
                     tau = NULL, 
                     min = NULL, 
                     max = NULL,
                    step = NULL, 
                allinone = FALSE, 
                no.title = FALSE, 
                     ...) 
{
     # gamlss.bi.list<-c("BI", "Binomial", "BB", "Beta Binomial") # binomial denominators
      if (is.null(min))   stop(paste("The min has not been specified"))
      if (is.null(max))   stop(paste("The max has not been specified"))
      if (is.null(step))  stop(paste("The step has not been specified"))
     interval <- seq(from=min, to=max, by=step)     
    # if the a gamlss object exist  
    if(!is.null(obj))
    {      
       if (!is.gamlss(obj))  stop(paste("This is not an gamlss object", "\n", ""))
          lobs <- length(obs) 
      if (lobs > 9)   stop(paste("Use up to eight observations for ploting"))
         plots <- list(c(1,1),c(2,1),c(3,1),c(2,2),c(3,2),c(3,2),c(4,2),c(4,2))
           vvv <- unlist(plots[lobs])     
      if(obj$type=="Discrete")   {typelh <- "h"; step <- 1 } else   {typelh <- "l" }          
        linter <- length(interval)
        family <- as.gamlss.family(obj$family[1])
       if(any(family%in%.gamlss.bi.list))  
       {
       stop(paste("This function is not working for BI and BB distibutions", "\n", ""))  
       }
         y.var <- seq(from=min, to=max, by=step)
    pdfunction <- rep(0, length(interval))   
    G.dev.expr <- body(family$G.dev.inc)  
        # nopar <- length(obj$parameters)
        # param <- obj$parameters
      if ("mu"%in%names(family$parameters))
             { mu.var <- fitted(obj)[obs]
               if (!family$mu.valid(mu.var))  stop( "`mu' parameter out of range")
             }
      if ("sigma"%in%names(family$parameters))
             {  sigma.var <- fitted(obj,"sigma")[obs]
               if (!family$sigma.valid(sigma.var))  stop( "`sigma' parameter out of range")
             }
      if ("nu"%in%names(family$parameters)) 
             {  nu.var <- fitted(obj,"nu")[obs]
               if (!family$nu.valid(nu.var))  stop( "`nu' parameter out of range")
             }
      if ("tau"%in%names(family$parameters)) 
             { tau.var <- fitted(obj,"tau")[obs]
               if (!family$tau.valid(tau.var))  stop( "`tau' parameter out of range") 
             }   
             
           pdfArr <- array(0:0, c(linter,lobs))
      title.label <- array(0:0,c(1,lobs))
    t.title.label <- array(0:0,c(1,lobs))
     op <- par(mfrow=vvv, mar=par("mar")+c(0,1,0,0), 
                    col.axis="blue4",col.main="blue4",col.lab="blue4", cex=0.6, font=2) 
     for (j in 1:lobs)
     {
        for (i in 1:linter)
        {
          if ("mu"%in%names(family$parameters))     assign("mu", mu.var[j])
          if ("sigma"%in%names(family$parameters))  assign("sigma", sigma.var[j])
          if ("nu"%in%names(family$parameters))     assign("nu", nu.var[j])
          if ("tau"%in%names(family$parameters))    assign("tau", tau.var[j])  
            assign("y", y.var[i])
            pdfunction[i] <- exp(-eval(G.dev.expr)/2)  
        }
        pdfArr[,j] <- pdfunction 
        if (!no.title)
        {
        mu.title <- if (!is.null(family$parameters$mu)) paste(" mu = ", format(signif(mu.var[j],4)),
                        ifelse(!is.null(family$parameters$sigma),",", " "), sep="")
        sigma.title <- if (!is.null(family$parameters$sigma)) paste(" sigma = ", 
                           format(signif(sigma.var[j],4)),
                        ifelse(!is.null(family$parameters$nu),",", " "), sep="")
        nu.title <- if (!is.null(family$parameters$nu)) paste(" nu = ", format(signif(nu.var[j],4)),
                        ifelse(!is.null(family$parameters$tau),",", " "), sep="")
        tau.title <- if (!is.null(family$parameters$tau)) paste(" tau = ", format(signif(tau.var[j],4)),
                         "", sep="")
        if (!allinone ) long.dist.name <- paste(family$family[2], ", ", 
                            family$family[1], "\n\n", sep="") else long.dist.name <- ""
        title.label[,j] <- paste(long.dist.name,family$family[1],"(",mu.title, 
                                sigma.title, nu.title, tau.title, ")", sep="")
        }
        else
            title.label[,j] <- ""
       if(obj$type=="Discrete")  y.title <- "pf, f(y)" else y.title <- "pdf, f(y)"   
       plot(y.var , pdfArr[,j], xlab="y",  
         ylab=y.title, 
         main=title.label[,j],
         col="darkgreen", 
         frame.plot = TRUE,
         type=typelh, lty=1, lwd=1)   
     }
 
    par(op)         
    }
    else # if only family
    {    
          lobs <- max(c(length(mu),length(sigma),length(nu),length(tau)))
          if (lobs > 9)   stop(paste("Use up to eight observations for ploting"))
         plots <- list(c(1,1),c(2,1),c(3,1),c(2,2),c(3,2),c(3,2),c(4,2),c(4,2))
           vvv <- unlist(plots[lobs])  
         y.var <- seq(from=min, to=max, by=step)   
        family <- as.gamlss.family(family)
         if(any(family$family[1]%in%.gamlss.bi.list)) bd <- max     
    pdfunction <- rep(0, length(interval))   
    G.dev.expr <- body(family$G.dev.inc) 
         if(family$type=="Discrete")   {typelh <- "h" } else   {typelh <- "l" }    
         if ("mu"%in%names(family$parameters))
             {  mu.var <- rep(mu,lobs)[1:lobs]
               if (!family$mu.valid(mu.var))  stop( "`mu' parameter out of range")
             }
         if ("sigma"%in%names(family$parameters))
             {  sigma.var <- rep(sigma,lobs)[1:lobs]
               if (!family$sigma.valid(sigma.var))  stop( "`sigma' parameter out of range")
             }
         if ("nu"%in%names(family$parameters)) 
             {  nu.var <- rep(nu,lobs)[1:lobs]
               if (!family$nu.valid(nu.var))  stop( "`nu' parameter out of range")
             }
         if ("tau"%in%names(family$parameters)) 
             {  tau.var <- rep(tau,lobs)[1:lobs]
               if (!family$tau.valid(tau.var))  stop( "`tau' parameter out of range") 
             }    
        if (!family$y.valid(y.var))  stop( "response variable out of range")
  
        pdfArr <- array(0:0, c(length(interval),lobs))
   title.label <- array(0:0,c(1,lobs))
 t.title.label <- array(0:0,c(1,lobs))
            op <- par(mfrow=vvv, mar=par("mar")+c(0,1,0,0), 
                    col.axis="blue4",col.main="blue4",col.lab="blue4", cex=0.6, font=2) 
    for (j in 1:lobs)
    {
        for (i in 1:length(interval))
        {
          if ("mu"%in%names(family$parameters))     assign("mu", mu.var[j])
          if ("sigma"%in%names(family$parameters))  assign("sigma", sigma.var[j])
          if ("nu"%in%names(family$parameters))     assign("nu", nu.var[j])
          if ("tau"%in%names(family$parameters))    assign("tau", tau.var[j])
          if(any(family%in%.gamlss.bi.list))       assign("bd",bd)   
          assign("y", y.var[i])
     pdfunction[i] <- exp(-eval(G.dev.expr)/2)  
        }
        pdfArr[,j] <- pdfunction 
         if (!no.title)
        {
        mu.title <- if (!is.null(family$parameters$mu)) paste(" mu = ", format(signif(mu.var[j],4)),
                        ifelse(!is.null(family$parameters$sigma),",", " "), sep="")
        sigma.title <- if (!is.null(family$parameters$sigma)) paste(" sigma = ", 
                           format(signif(sigma.var[j],4)),
                        ifelse(!is.null(family$parameters$nu),",", " "), sep="")
        nu.title <- if (!is.null(family$parameters$nu)) paste(" nu = ", format(signif(nu.var[j],4)),
                        ifelse(!is.null(family$parameters$tau),",", " "), sep="")
        tau.title <- if (!is.null(family$parameters$tau)) paste(" tau = ", format(signif(tau.var[j],4)),
                         "", sep="")
        if (!allinone) long.dist.name <- paste(family$family[2], ", ", 
                            family$family[1], "\n\n", sep="") else long.dist.name <- ""
        title.label[,j] <- paste(long.dist.name,family$family[1],"(",mu.title, 
                                sigma.title, nu.title, tau.title, ")", sep="")
        }
        else
            title.label[,j] <- ""     
       if(family$type=="Discrete")  y.title <- "pf, f(y)" else y.title <- "pdf, f(y)"   
       plot(y.var , pdfArr[,j], xlab="y",  
         ylab=y.title, 
         main=title.label[,j],
         col="darkgreen", 
         frame.plot = TRUE,
         type=typelh, lty=1, lwd=1)   
     }
   par(op)  
   }  
  par(op)  
}
