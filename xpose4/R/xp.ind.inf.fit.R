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

"xp.ind.inf.fit" <-
  function(plot.ids=TRUE,
           idscex=0.7,
           ptscex=0.7,
           title = NULL,
           recur = FALSE,
           xlb = NULL,
           ylb = NULL,
           gamobj=NULL,
           ...){

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
      
      
      sd <- sqrt(eval(parse(text=paste("current.gam","$deviance",sep="")))/eval(parse(text=paste("current.gam","$df.residual",sep=""))))

      #sd   <- sqrt(current.gam$deviance/current.gam$df.residual)
      

      #pear <- residuals(current.gam,type="pearson")/sd
      #h    <- lm.influence(current.gam)$hat
      #p    <- current.gam$rank

      pear <- residuals(eval(parse(text="current.gam")),type="pearson")/sd
      h    <- lm.influence(eval(parse(text="current.gam")))$hat
      p    <- eval(parse(text=paste("current.gam","$rank",sep="")))

      rp   <- pear/sqrt(1-h)
                                        #for (i in 1:length(rp)){
    #  if(is.na(rp[i])){
    #    rp[i] <- pear[i]
    #  }
    #}
    cook <- (h*rp^2)/((1-h)*p)

    #for (i in 1:length(cook)){
    #  if(is.na(cook[i])){
    #    cook[i] <- h[i]*rp[i]^2
    #  }
    #}

      
      #n    <- p + current.gam$df.residual
      n    <- p + eval(parse(text=paste("current.gam","$df.residual",sep="")))
    cooky <- 8/(n-2*p)
    ##    hy   <- (2*p)/(n-2*p)
    hy   <- (2*p)/n
    leve <- h/(1-h)
    #for (i in 1:length(leve)){
    #  if(!is.finite(leve[i])){
    #    leve[i] <- h[i]
    #  }
    #}


    if(is.null(xlb))
      xlb <- "Leverage (h/(1-h))"
    if(is.null(ylb))
      ylb <- "Cooks distance"
      
      


      if(is.null(title)) {
      title <- paste("Individual influence on the GAM fit for ",
                     eval(parse(text=paste("current.gam","$pars",sep=""))),
                     " (run ",
                     #current.gam$runno,
                     eval(parse(text="current.gam$runno")),
                     ")",sep="")
    }

    ## Get the idlabs
    if(any(is.null(eval(parse(text="current.gam$data$ID"))))){
      ids <- "n"
    } else {
      ids <- eval(parse(text="current.gam$data$ID"))
    }

    ## inform about NaN and INF values
    for (i in 1:length(cook)){
      if(is.na(cook[i])||!is.finite(cook[i])||
         is.na(leve[i])||!is.finite(leve[i])){
        cat("\nFor ID ",ids[i], ":\n", sep="")
        cat("   Cook distance is ", cook[i],"\n",sep="")
        cat("   Leverage is ", leve[i],"\n",sep="")
        cat("   => the point is not included in the plot\n")
      }
    }

    xplot <- xyplot(cook~leve,
                    ylab=ylb,
                    xlab=xlb,
                    main=title,
                    aspect=1,
                    cooky=cooky,
                    hy=hy,
                    scales = list(cex=0.7,tck=-0.01),
                    ids = eval(parse(text="current.gam$data[,1]")),
                    panel=
                    function(x,y,ids,...) {
                      if(!any(ids == "n")&& plot.ids==TRUE) {
                        addid(x,y,ids=ids,
                              idsmode=TRUE,
                              idsext =0.05,
                              idscex = idscex,
                              idsdir = "both")
                      } else {
                        panel.xyplot(x,y,cex=ptscex,col="black")
                      }
                    }
                    )

    return(xplot)

}

