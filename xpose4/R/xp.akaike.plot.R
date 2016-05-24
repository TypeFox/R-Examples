# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

"xp.akaike.plot"<-
  function(title = NULL,
           xlb = "Akaike value",
           ylb="Models",
           gamobj=NULL,
           ...) {

    if(is.null(gamobj)){
      gamobj <- check.gamobj()
      if(is.null(gamobj)){
        return()
      } else {
      }
    } else {
      c1 <- call("assign",pos=1, "current.gam", gamobj,immediate=T)
      eval(c1)
    }
    ##eval(parse(text=paste("current.gam","$steppit",sep="")))
    ##if(is.null(current.gam$steppit)) {
    if(is.null(eval(parse(text=paste("current.gam","$steppit",sep=""))))) {
      cat("This plot is not applicable without stepwise covariate selection.\n")
      return()
    }
    
    
    keep <- eval(parse(text=paste("current.gam","$keep",sep=""))) #current.gam$keep
    aic <- apply(keep, 2, function(x)
                 return(x$AIC))
    df.resid <- apply(keep, 2, function(x)
                      return(x$df.resid))
    term <- apply(keep, 2, function(x)
                  return(x$term))
    pdata <- data.frame(aic, df.resid, term)
    aic.ord <- order(pdata$aic)
    pdata <- pdata[aic.ord,  ]

    ##
    ## Select the 30 models with lowest AIC
    ##
    if(dim(pdata)[1] > 30){
      pdata1 <- pdata[1:30,  ]
      pdata2 <- pdata[1:30,  ]
    } else {
      pdata1 <- pdata
      pdata2 <- pdata
    }
    pdata1$term <- unclass(pdata1$term)
    pdata1$term <- reorder(as.factor(pdata1$term), pdata1$aic)
    names(pdata1$term) <- pdata2$term
    
    if(is.null(title)) {
      title <- paste("AIC values from stepwise GAM search on ",
                     eval(parse(text=paste("current.gam","$pars",sep=""))),
                     #current.gam$pars,
                     " (Run ",
                     eval(parse(text=paste("current.gam","$runno",sep=""))),
                     #current.gam$runno,
                     ")",sep="")
    }
        
    xplot <- dotplot(term~aic,
                     pdata1,
                     main=title,
                     xlab=xlb,
                     ylab=ylb,
                     scales=list(cex=0.7,
                       tck=-0.01,
                       y=list(labels=pdata2$term,cex=0.6 )
                       ),
                     ...
                     )

    #print(xplot)
    return(xplot)
    #invisible()
    
    
  }
