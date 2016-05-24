#############################################################################################
## Project: PKgraph
## File: ggobi.R
## Author: Xiaoyong Sun
## Date: 11/17/2009
## Goal: PKgraph
##        - handle ggobi function
## Notes:
#############################################################################################

ggobiMap <- function()
{
    ## check data exist
    if (!checkDataExist())
    {
        ErrorMessage("No data exist!")
        return(invisible(NULL))
    }



}

ggobi.data <- function()
{
    labelMessage1 <- "Please move data for diagnostics from left TABLE to right TABLE."
    labelMessage2 <- "After choosing data, click to interactive diagnostics."
    winTitle <- "Configure datasets"
    statusMessage <- "Data is ready for interactive diagnostics."

    selectDataDialog(winTitle, labelMessage1, labelMessage2, statusMessage, menuOption=2)
}

ggobi.map <- function(Message)
{
   # check data exist
    choose.dataname <- getItDataName()
    if (is.null(choose.dataname))
    {
        ErrorMessage("Please choose specific data first using submenu above.")
        return(invisible(NULL))
    }
    
    ggobi.map = gwindow("Map data", horizontal=FALSE)

    gtgroup1 = ggroup(cont=ggobi.map, horizontal=FALSE)

    gf1 <- gframe(text = "", markup = FALSE, pos = 0, horizontal=TRUE, container = gtgroup1)
    tbl <- glayout(cont=gf1)
    
      cline <- 0
      tbl.list <- list()
      for (i in 1: length(choose.dataname))
      {
          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "Data Name"
          tbl[cline, 2, anchor = c(-1,-1)] = choose.dataname[i]

          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "Mapping key"
          tbl.list[[i]] = gdroplist(items=colnames(getCurrentData(choose.dataname[i])))
          tbl[cline, 2, anchor = c(-1,-1)] = tbl.list[[i]]

      }

    #tbl2 <- glayout(cont=gtgroup2)

    gb1 = gbutton(text="Map data", horizontal=FALSE )
    addhandlerclicked(gb1, function(h,...)
                    {
                        key <- rep("", length(choose.dataname))
                        sapply(1:length(tbl.list), function(i) key[i] <<- svalue(tbl.list[[i]]))
                        setItMap(key)
                        dispose(ggobi.map)
                        svalue(pmg.statusBar) <- "Data is mapped for interactive graphics."

                    })
    cline <- cline + 1
    tbl[cline, 2, anchor = c(-1,-1)] = gb1
}

ggobi.diagose <- function()
{
    #check data exist
    key <- getItMap()
    if (is.null(key))
    {
        ErrorMessage("Please choose data and map data first!")
        return(invisible(NULL))
    }
    
    gwin <- gwindow(parent=NULL, height=getSubHeight(), width=getSubWidth())

    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, expand=FALSE) # expand

    #size(g1) =   c(rightWidth*0.3,mainHeight*0.3)
    glabel("Step 1: Map data", cont=g1)
    gf = gframe(horizontal=FALSE, expand=FALSE, cont=g1) # expand
          size(gf) <- c(getSubWidth()*0.5, getSubHeight()*0.8)
      tbl <- glayout(cont=gf)

    data.name <- getItDataName()
    plot.list <- list()

      cline <- 0
      for (i in 1: length(data.name))
      {
          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "Data Name:"
          plot.list[[data.name[i]]] <- glabel(text=data.name[i])
          tbl[cline, 2, anchor = c(-1,-1)] = plot.list[[data.name[i]]]

          cline <- cline + 1
          tbl[cline, 1, anchor = c(0,0)] = "x"
          xname <- paste(data.name[i], "_x", sep="")
          plot.list[[xname]] <- gdroplist(items=c("", colnames(getCurrentData(data.name[i]))))
          tbl[cline, 2, anchor = c(0,0)] = plot.list[[xname]]

          cline <- cline + 1
          tbl[cline, 1, anchor = c(0,0)] = "y"
          yname <- paste(data.name[i], "_y", sep="")
          plot.list[[yname]] <- gdroplist(items=c("", colnames(getCurrentData(data.name[i]))))
          tbl[cline, 2, anchor = c(0,0)] = plot.list[[yname]]
      }

    bg1 = ggroup(cont=g1)

     addSpace(bg1, getSubWidth()/10, horizontal=TRUE)
     mapButton = gbutton("Map Data", cont=bg1,
        handler = NULL)

     addSpace(bg1, getSubWidth()/10, horizontal=TRUE)
     removeButton = gbutton("Remove selection", cont=bg1,
        handler = NULL)
        
    
    g2 = ggroup(horizontal=FALSE, expand=FALSE) # expand
    #size(g2) =   c(rightWidth*0.3,mainHeight*0.3)
    glabel("Step 2: Configure plots", cont=g2)
    g22 = ggroup(horizontal=FALSE, expand=FALSE, cont=g2) # expand
    gt <- gtable(items=c(letters[1:24]), cont=g22, expand=TRUE)
    size(gt) <- c(getSubWidth()*0.5, getSubHeight()*0.8)
    gt[] <- c()

    bg2 = ggroup(cont=g2)
        
     addSpace(bg2, getSubWidth()/15, horizontal=TRUE)
     plotButton = gbutton("Plot selection", cont=bg2,
        handler = NULL)

     addSpace(bg2, getSubWidth()/15, horizontal=TRUE)
     allButton = gbutton("Plot all", cont=bg2,
        handler = NULL)

     addSpace(bg2, getSubWidth()/15, horizontal=TRUE)
     closeButton = gbutton("Close plots", cont=bg2,
        handler = NULL)
        
    leftpane =gpanedgroup(g1)
    rightpane =gpanedgroup(g2)

    part.pg = gpanedgroup(leftpane, rightpane)


    add(gwin, part.pg)
    
    addhandlerclicked(mapButton, handler=function(h,...)
                      {
                                key <- rep("", length(data.name))
                                mycheck <- sapply(1:length(data.name), function(i)
                                          {
                                              if (svalue(plot.list[[i*3-1]])=="")
                                              {
                                                  ErrorMessage(paste(data.name[i], ": No x value is selected!", sep=" "))
                                                  return(TRUE)
                                              }
                                              else
                                              {
                                                  if (svalue(plot.list[[i*3]])=="")
                                                  {
                                                      key[i] <<- paste(svalue(plot.list[[i*3-2]]), svalue(plot.list[[i*3-1]]), sep=" ")
                                                  }
                                                  else
                                                  {
                                                     key[i] <<- paste(svalue(plot.list[[i*3-2]]), svalue(plot.list[[i*3]]), "vs", svalue(plot.list[[i*3-1]]), sep=" ")
                                                  }
                                              }
                                              return(FALSE)
                                          })
                                if (any(mycheck))
                                {
                                    return(invisible(NULL))
                                }

                            old <- gt[]
                  
                            if (length(old)==1 && is.na(old))
                            {
                                gt[] <- key
                            }
                            else
                            {
                                #gt[] <- c(old, key)
                                
                                if (any(key %in% old))
                                {
                                    delete.key <- which(key %in% old)
                                    gt[] <- c(old, key[-delete.key])
                                }
                                else
                                {
                                    gt[] <- c(old,key)
                                }
                            }
                                          




                      })


    addhandlerclicked(removeButton, handler=function(h,...)
                      {
                          select.value <- svalue(gt)
                          all.value <- gt[]
                          
                          if (length(select.value)==0)
                          {
                              ErrorMessage("Please select value first in the right table!")
                              return(invisible(NULL))

                          }
                          
                          select.ind <- which(gt[]==select.value)
                          gt[] <<- all.value[-select.ind[1]]
                          
                      })    

    addhandlerclicked(plotButton, handler=function(h,...)
                      {
                          select.value <- svalue(gt)
                          #all.value <- gt[]

                          if (length(select.value)==0)
                          {
                              ErrorMessage("Please select value first in the right table!")
                              return(invisible(NULL))

                          }

                          ## data manage
                          select.data <- data.manage()
                          select.names <- names(select.data)
                          ggobi.object <- ggobi(select.data[1])
                          for (i in 2: length(select.data))
                          {
                              ggobi.object[i] <- select.data[i]
                          }

                          split.value <- unlist(strsplit(select.value, " "))

                          ## split =2, hist
                          ## split =4, scatterplot
                          ## split.value: dataname, y, vs, x
                          if (length(split.value)==2)
                          {
                              data.index <- which(split.value[1]==select.names)
                              myx <- paste("X", split.value[1], ".",split.value[2], sep="")
			if (Sys.info()[["sysname"]] != "Windows")
			{
				myx <- gsub("/", ".", myx)
			}
                              display(ggobi.object[data.index], pmode="Barchart", vars=list(X= myx))
                          }
                          else
                          {

                              data.index <- which(split.value[1]==select.names)
                              myx <- paste("X", split.value[1], ".",split.value[4], sep="")
			if (Sys.info()[["sysname"]] != "Windows")
			{
				myx <- gsub("/", ".", myx)
			}
                              myy <- paste("X", split.value[1], ".",split.value[2], sep="")
			if (Sys.info()[["sysname"]] != "Windows")
			{
				myy <- gsub("/", ".", myy)
			}
                              display(ggobi.object[data.index], pmode="Scatterplot Display", vars=list(X= myx, Y= myy))
                          }

                      })
                      
    addhandlerclicked(allButton, handler=function(h,...)
                      {
                          select.value <- gt[]
                          #all.value <- gt[]

                          if (length(select.value)==0)
                          {
                              ErrorMessage("Please select value first in the right table!")
                              return(invisible(NULL))

                          }

                          ## data manage
                          select.data <- data.manage()
                          select.names <- names(select.data)
                          ggobi.object <- ggobi(select.data[[1]])
                          for (i in 2: length(select.data))
                          {
                              ggobi.object[i] <- select.data[[i]]
                          }

                          for ( i in 1: length(select.value))
                          {
                              split.value <- unlist(strsplit(select.value[i], " "))

                              ## split =2, hist
                              ## split =4, scatterplot
                              ## split.value: dataname, y, vs, x
                              if (length(split.value)==2)
                              {
                                  data.index <- which(split.value[1]==select.names)
                                  myx <-  split.value[4]
                        			    if (Sys.info()[["sysname"]] != "Windows")
                        			    {
                        				myx <- gsub("/", ".", myx)
                        			    }
                                  display(ggobi.object[data.index], pmode="Barchart", vars=list(X= myx))
                              }
                              else
                              {

                                  data.index <- which(split.value[1]==select.names)
                                  myx <-  split.value[4]
                        			    if (Sys.info()[["sysname"]] != "Windows")
                        			    {
                        				      myx <- gsub("/", ".", myx)
                        			    }
                                  myy <- split.value[2]
                        			    if (Sys.info()[["sysname"]] != "Windows")
                        			    {
                              				myy <- gsub("/", ".", myy)
                        			    }
                                  display(ggobi.object[data.index], pmode="Scatterplot Display", vars=list(X= myx, Y= myy))
                              }

                           }
                      })
                      
    addhandlerclicked(closeButton, handler=function(h,...)
                      {
                          close(ggobi_get())
                      })

}

## map data with key id
data.manage <- function()
{
    map.dataset.names <- getItDataName()
    map.ids <- getItMap()
    map.data <- list()
    
    sapply(1:length(map.dataset.names), function(i)
           {
              tmp.data <- getCurrentData(map.dataset.names[i])
              map.col <- colnames(tmp.data)
              map.col <- gsub(map.ids[i], map.ids[1], map.col)
              colnames(tmp.data) <- map.col
              map.data[[map.dataset.names[i]]] <<- tmp.data

           })
              
    return(map.data)
    
}
