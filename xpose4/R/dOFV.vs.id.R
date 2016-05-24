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

dOFV.vs.id <-
  function(xpdb1,
           xpdb2,
           sig.drop=-3.84,
           decrease.label.number=3,
           increase.label.number=3,
           id.lab.cex=0.6,
           id.lab.pos=2,
           type="o",
           xlb="Number of subjects removed",
           ylb=expression(paste(Delta,"OFV")),
           main="Default",
           sig.line.col = "red",
           sig.line.lty = "dotted",
           tot.line.col = "grey",
           tot.line.lty = "dashed",
           key=list(#title = expression(paste("Individual influence on ",Delta,"OFV")),
             columns = 1,
             lines = list(pch = c(super.sym$pch[1:2],NA,NA),
                          type=list("o","o","l","l"),
                          col = c(super.sym$col[1:2],sig.line.col,tot.line.col),
                          lty = c(super.sym$lty[1:2],sig.line.lty,tot.line.lty)
             ),
             ## points = list(pch = c(super.sym$pch[1:2],NA),
             ## col = c(super.sym$col[1:2],"red")),
             text = list(c(
               expression(paste(Delta, OFV[i] < 0)),
               expression(paste(Delta, OFV[i] > 0)),
               expression(paste("Significant  ",Delta, OFV)),
               expression(paste("Total  ",Delta, OFV))
               ##"Individuals with OFV drop",
               ##"Individuals with OFV increase",
               ##"Significant drop in ")
             )),
             ##space="right",
             corner=c(0.95,0.5),border=T
           ),
           ...)
{
    
    
    if(is.null(xpdb2)){
      cat("Comparison database needed for this plot!")
      return(NULL)
    }
    
    #require(graphics)
    
    iv1 <- xpdb1@Data.firstonly
    ##str(iv2)
    
    iv2 <- xpdb2@Data.firstonly
    ##str(iv2)
    
    if(!all(iv1$ID == iv2$ID)){
      cat("All ID labels for both databases must match\n")
      return(NULL)
    }
    
    comp.frame <- data.frame(id=iv1$ID,obj1=iv1$OBJ,obj2=iv2$OBJ)
    
    # fix for R version 3.1.0
    if(is.factor(comp.frame$obj1)) comp.frame$obj1 <- as.numeric(levels(comp.frame$obj1))[comp.frame$obj1] 
    if(is.factor(comp.frame$obj2)) comp.frame$obj2 <- as.numeric(levels(comp.frame$obj2))[comp.frame$obj2] 
    comp.frame$d.obj <- comp.frame$obj2 - comp.frame$obj1
    ##str(comp.frame)
    
    ##hist(comp.frame$d.obj)
    
    ofv1 <- sum(comp.frame$obj1)
    ofv2 <- sum(comp.frame$obj2)
    d.ofv <- sum(comp.frame$d.obj)
    ##ofv1
    ##ofv2
    ##d.ofv
    ##ofv2-ofv1
    
    
    comp.frame$rnk <- rank(comp.frame$d.obj,ties.method="random")
    ##comp.frame[comp.frame$rnk==74,]
    ##comp.frame[comp.frame$rnk==1,]
    
    max.rank <- max(comp.frame$rnk)
    id.increase <- c(NA)
    id.decrease<- c(NA)
    increase.sum <- c(d.ofv)
    decrease.sum <- c(d.ofv)
    
    for(i in 1:max.rank){
      ##i=1
      
      ## compute largest decrease
      tmp <- comp.frame[comp.frame$rnk==i,]
      if(tmp$d.obj<=0){
        id.decrease <- c(id.decrease,tmp$id)
        decrease.sum <- c(decrease.sum,sum(comp.frame[comp.frame$rnk > i,"d.obj"]))
      }
      
      ## compute largest increase
      tmp <- comp.frame[comp.frame$rnk==(max.rank-i+1),]
      if(tmp$d.obj>0){
        id.increase <- c(id.increase,tmp$id)
        increase.sum <- c(increase.sum,
                          sum(comp.frame[comp.frame$rnk < (max.rank - i +1),"d.obj"]))
      }
    }
    
    ##length(id.increase)
    ##length(id.decrease)
    
    n.removed <- c(1:max(length(id.increase),length(id.decrease)))-1
    ##n.removed
    
    if(length(id.increase) < length(n.removed)){
      n.na <- length(n.removed) - length(id.increase)
      id.increase <- c(id.increase,rep(NA,times=n.na))
      increase.sum <- c(increase.sum,rep(NA,times=n.na))
    }
    if(length(id.decrease) < length(n.removed)){
      n.na <- length(n.removed) - length(id.decrease)
      id.decrease <- c(id.decrease,rep(NA,times=n.na))
      decrease.sum <- c(decrease.sum,rep(NA,times=n.na))
    }
    
    plot.frame <- data.frame(n.removed,
                             id.increase,
                             increase.sum,
                             id.decrease,
                             decrease.sum)
    
    
    
    ##plot.frame
    
    default.plot.title <- paste("Individual influence on change in OFV\n",
                                "(Run",xpdb2@Runno," - Run",xpdb1@Runno,")",sep="" )
    plotTitle <- xpose.multiple.plot.title(object=xpdb1,
                                           plot.text = default.plot.title,
                                           main=main,
                                           no.runno=T,
                                           ...)
    
    ##id label
    idlab.decrease <- subset(plot.frame,n.removed > 0 & n.removed <= decrease.label.number)
    idlab.increase <- subset(plot.frame,n.removed > 0 & n.removed <= increase.label.number)
    
    super.sym <- trellis.par.get("superpose.symbol")
    xplot <- xyplot(decrease.sum+increase.sum~n.removed,data=plot.frame,
                    type=type,
                    xlab=xlb,
                    ylab=ylb,
                    sig.line.col = sig.line.col,
                    sig.line.lty = sig.line.lty,
                    tot.line.col = tot.line.col,
                    tot.line.lty = tot.line.lty,
                    id.lab.cex = id.lab.cex,
                    id.lab.pos = id.lab.pos,
                    panel=function(x,y,
                                   sig.line.col,sig.line.lty,
                                   tot.line.col,tot.line.lty,
                                   id.lab.cex,id.lab.pos,
                                   ...){
                      panel.xyplot(x,y,...)
                      panel.abline(h=sig.drop,col=sig.line.col,lty=sig.line.lty,...)
                      panel.abline(h=y[1],col=tot.line.col,lty=tot.line.lty,...)
                      ltext(x=idlab.decrease$n.removed,y=idlab.decrease$decrease.sum,
                            label=idlab.decrease$id.decrease,
                            cex=id.lab.cex,
                            pos=id.lab.pos,
                            #adj=c(2,1),
                            ...)
                      ltext(x=idlab.increase$n.removed,y=idlab.increase$increase.sum,
                            label=idlab.increase$id.increase,
                            cex=id.lab.cex,
                            pos=id.lab.pos,
                            #adj=c(2,1),
                            ...)
                    },
                    ##auto.key =T
                    main= plotTitle,
                    key = key,
                    ...
    )
    return(xplot)
  }
