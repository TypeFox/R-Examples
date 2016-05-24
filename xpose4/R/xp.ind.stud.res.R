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

"xp.ind.stud.res" <-
  function(title = NULL,
           recur = FALSE,
           xlb = NULL,
           ylb = NULL,
           gamobj=NULL){

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

    if(eval(parse(text="current.gam$family$family")) == "gaussian"){
      sd <- sqrt(eval(parse(text="current.gam$deviance"))/eval(parse(text="current.gam$df.residual")))
    } else {
      sd <- 1
    }

    dev  <- residuals(eval(parse(text="current.gam")), type = "deviance")/sd
    pear <- residuals(eval(parse(text="current.gam")), type = "pearson")/sd
    h    <- lm.influence(eval(parse(text="current.gam")))$hat
    rp   <- pear/sqrt(1 - h)
    #for (i in 1:length(rp)){
    #  if(is.na(rp[i])){
    #    rp[i] <- pear[i]
    #  }
    #}

    sgn <- function(x)
      ifelse(x > 0, 1, ifelse(x < 0, -1, 0))

    res  <- sgn(dev) * sqrt(dev^2 + h * rp^2)

    pdata <- data.frame(cbind(eval(parse(text="current.gam$data[,1]")),res))
    names(pdata) <- c("ID","studres")
    studres.ord <- order(pdata$studres)
    pdata <- pdata[studres.ord,  ]
    pdata$ID <- reorder(as.factor(pdata$ID),pdata$studres)

    if(is.null(xlb))
      xlb <- "Studentized residual"
    if(is.null(ylb))
      ylb <- "ID"

    if(is.null(title)) {
      title <- paste("Studentized residual of the GAM fit for ", eval(parse(text="current.gam$pars"))," (Run ",
                     eval(parse(text="current.gam$runno")), ")",sep="")
    }

#    xplot <- dotchart(pdata$studres,labels=as.character(pdata$ID),
#                      main = title, xlab = xlb,ylab=ylb,cex=0.6)

    xplot <- dotplot(ID~studres,
                     pdata,
                     main=title,
                     xlab=xlb,
                     ylab=ylb,
                     scales=list(
                       tck=-0.01,
                       y=list(cex=0.6 )
                       )
                     )
    #print(xplot)
    return(xplot)
    #invisible()

  }
