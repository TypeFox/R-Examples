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

"parm.summary" <- function(object,
                          onlyfirst=TRUE,
                          subset=xsubset(object),
                          inclZeroWRES=FALSE,
                          out.file=".screen", # can be ".ask" ".graph" or a file name,
                          #out.file.sep=",",
                          main="Default",
                          fill = "gray",
                          values.to.use=xvardef("parms",object),
                          value.name="Parameter",
                           max.plots.per.page=1,
                          ...){

  ##  if(is.null(object@Prefs@Xvardef$parms)) 
  if(is.null(values.to.use)){
    cat("The current database has no", value.name,"defined!\n")
    return()
  }


  data <- Data(object,onlyfirst=onlyfirst,subset=subset,inclZeroWRES=inclZeroWRES)

  #if(any(is.null(data))) {
  #  return("The subset expression is invalid.")
  #}
  
  #data <- object@Data

  ##parnams <- object@Prefs@Xvardef$parms
  parnams <- values.to.use

  cats  <- NULL
  conts <- NULL
  
  for(parm in parnams) {
    if(is.factor(data[[parm]])) {
      cats <- c(cats,parm)
    } else {
      conts <- c(conts,parm)
    }
  } 
  
  if(!is.null(cats)) {
    cat.mat <- categorical.table(object, cats, onlyfirst=onlyfirst, subset=subset,
       inclZeroWRES=inclZeroWRES)
  }
  if(!is.null(conts)) {
    con.mat <- continuous.table(object, conts, onlyfirst=onlyfirst, subset=subset,
       inclZeroWRES=inclZeroWRES)
  }

  if (out.file==".ask"){
    cat("Would you like to export the table(s) as a text file? n(y)\n")
    out.to.text <- readline()
    if(out.to.text == "y") {
      cat("Please type a filename (excluding the .csv extension):\n")
      out.file <- readline()
    } else {
      cat("Would you like the table to be output as a graph? n(y)\n")
      out.to.text <- readline()
      if(out.to.text == "y") {
        out.file <- ".graph"
      } else {
        out.file <- ".screen"
      }
    }
  }
  
  if (out.file==".screen" | out.file==".graph"){
    
    if (out.file==".screen"){
      if(!is.null(cats))
        print.char.matrix(cat.mat)
      cat("\n")
      if(!is.null(conts))
        print.char.matrix(con.mat)
    }
    if (out.file==".graph"){
      table.list=list()
      iii <- 1
      if(!is.null(conts)) table.list[[iii]] <- con.mat ; iii <- iii+1
      if(!is.null(cats)) {
        ## find max height for row
        num.rows <- dim(cat.mat)[2]
        num.cols <- dim(cat.mat)[1]
        for(jjj in 1:num.cols){
          max.lines <- 1
          num.lines <- rep(1,num.rows)
          cell.ht <- gregexpr("\n",cat.mat[jjj,])
          for(k in 1:num.rows){
            if(all(cell.ht[[k]]==-1)) {
              num.lines[k]=1
            } else {
              num.lines[k]=length(cell.ht[[k]])+1
            }
          }
          
          max.lines <- max(num.lines)
          line.diff <- max.lines-num.lines
          
          for(kk in 1:num.rows){
            tmp <- paste(rep("\n",line.diff[kk]),sep="",collapse="")
            cat.mat[jjj,kk] <- paste(cat.mat[jjj,kk],tmp,sep="")
          }
        }
        table.list[[iii]] <- cat.mat 
        iii <- iii+1
      }

      #table.list <- list(con.mat,cat.mat)
      num.tables <- length(table.list)
      plotList <- vector("list",num.tables)
      vp1 <- viewport(x=0, y=1, just=c("left","top"),
                      width=1, height=0.9,
                      gp=gpar(#lineheight=1.0,
                        cex=0.9#txt.cex,
                        ##  font=0.01#txt.font
                        ),
                      name="vp1")
      for(jj in 1:num.tables){
        psobj <- table.list[[jj]]
        ##iter  <- 7 * length(xvardef("parms", object))
        ##iter  <- 7 * length(conts)
        if(is.null(psobj)) break
        
        cols <- psobj[1,]
        
        nr <- dim(psobj)[1]
        nc <- dim(psobj)[2]
        
                                        #grid.newpage()
                                        #      xpose.multiple.plot.default(list(1),plotTitle=plotTitle,...)
        
        textColumnList <- vector("list",nc) 
        for(ii in 1:nc){
          textColumnList[[ii]] <- psobj[-1,ii]
        }

        xpose.table <- add.grid.table(textColumnList,col.nams=cols,ystart=unit(1,"npc"),
                                      vp=list(vp1),cell.padding=1,center.table=TRUE,
                                      fill.type="both",
                                      v.space.before=0.25,
                                      v.space.after=0.5,
                                      draw=FALSE,
                                      use.rect=TRUE,...)
        
        plotList[[jj]] <- xpose.table$xpose.table

        
      }
      default.plot.title <- paste(value.name,"Summary",sep=" ")
      plotTitle <- xpose.multiple.plot.title(object=object,
                                             plot.text = default.plot.title,
                                             main=main,
                                             subset=subset,
                                             ...)
      obj <- xpose.multiple.plot(plotList,plotTitle,...)
      return(obj)

    }
  } else {
    if(!is.null(cats)) {
      print(cat.mat, file = paste(out.file, ".csv", sep = ""),
            hsep=",",vsep="",csep="",top.border=FALSE,left.border=FALSE)
      #write.table(cat.mat, file = paste(out.file, ".csv", sep = ""),
      #            append = FALSE, quote = FALSE, sep = ",",
      #            row.names = FALSE,
      #            col.names = FALSE)
    }
    if(!is.null(conts)){
      write.table(con.mat, file = paste(out.file, ".csv", sep = ""),
                  append = TRUE, quote = FALSE, sep = ",",
                  row.names = FALSE,
                  col.names = FALSE)
      # print(con.mat, file = paste(out.file, ".csv", sep = ""),
      #       hsep=",",vsep=NULL,csep=NULL)
    }
    invisible()
  }

  #invisible()
  #return(cat(""))
  #return()

}

