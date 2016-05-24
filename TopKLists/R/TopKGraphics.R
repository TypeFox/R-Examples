TopKListsGUI <- function(lists, autorange.delta = FALSE, override.errors = TRUE,  venndiag.pdf.size = c(7, 7), venndiag.size = c(380, 420), gui.size = c(900, 810), directory = NULL, venndiag.res = 70, aggmap.res = 100) {
  requireNamespace("RGtk2")
  requireNamespace("gWidgets2")
  requireNamespace("gWidgetsRGtk2")
  options("guiToolkit"="RGtk2")

  ##setting up the directory
  if(is.null(directory)) {
    directory <- paste(tempdir(), "/TopKLists-temp",sep="")
    if(!file.exists(directory)) dir.create(directory)
  }
  message(paste("Writing files to", directory))

  
                              #delta_symbol = substitute(delta)
                              #nu_symbol = expression(nu)

                              #check if given lists are of type 'data.frame' and contains data
  if (!is.data.frame(lists)) {
    gWidgets2::gmessage("ERROR:\nGiven list-object is not of type data.frame!\nProgram will exit.", title = "Incorrect input", icon = "error")
    return("TopKLists GUI Error: Given list-object is not of type data.frame!")
  } else if ((dim(lists)[1] < 2) | (dim(lists)[2] < 2)) {
    gWidgets2::gmessage("ERROR:\nGiven list-object has eighter only one feature or only one list!\nProgram will exit.", title = "Incorrect input", icon = "error")
    return("TopKLists GUI Error: Given list-object has either only one feature or only one list!")
  }

  win <- gWidgets2::gwindow("TopKLists GUI", visible = FALSE, width = gui.size[1], height = gui.size[2]) #main window
  maingroup <- gWidgets2::ggroup(horizontal = FALSE, container = win, expand = TRUE) #main container
  wid <- environment() #contains all the widgets
  current.arguments <- environment() #contains all the arguments (N,L,v,....) of the current calculation

                                        #create input area for arguments
  calcframe <- gWidgets2::gframe("ARGUMENTS", pos = 0.5, container = maingroup)
  calcgroup <- gWidgets2::ggroup(horizontal = FALSE, container = calcframe, expand = TRUE)
  agroup <- gWidgets2::ggroup(horizontal = TRUE, container = calcgroup, expand = TRUE)
  a1frame <- gWidgets2::gframe("Multiple ranked lists input", container = agroup, expand = TRUE)
  a2frame <- gWidgets2::gframe("Select nu and choose range of delta", container = agroup, expand = TRUE)
  a3frame <- gWidgets2::gframe("Threshold", container = agroup, expand = TRUE)

  wid$calculate <- gWidgets2::gbutton("Calculate", container = calcgroup, handler = function(h,...) { calculate.all.truncationlists() })	
  gWidgets2::svalue(calcgroup) <- 10

                                        #layout of the first argument input area
  tbl1 <- gWidgets2::glayout(spacing = 10,container = a1frame)
  tbl1[2,2, anchor = c(1,0)] <- gWidgets2::glabel("N", container = tbl1)
  tbl1[2,3] <- wid$N <- gWidgets2::gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = dim(lists)[1], container = tbl1)
  tbl1[2,4, anchor = c(-1,0)] <- gWidgets2::glabel("number of objects" ,container = tbl1)
  tbl1[3,2, anchor = c(1,0)] <- gWidgets2::glabel("L ", container = tbl1)
  tbl1[3,3, anchor = c(1,0)] <- gWidgets2::glabel(dim(lists)[2], container = tbl1)
  tbl1[3,4, anchor = c(-1,0)] <- gWidgets2::glabel("number of lists", container = tbl1)

                                        #layout of the second argument input area
  tbl2 <- gWidgets2::glayout(spacing = 10, container = a2frame)
  tbl2[2,2, anchor = c(1,0)] <- gWidgets2::glabel(expression(nu), container = tbl2)
  tbl2[2,3] <- wid$v <- gWidgets2::gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 10, container = tbl2, handler = function(h,...) { update.deltarange() })
  tbl2[2,4:5, anchor = c(1,0)] <- gWidgets2::glabel("for list truncation", container = tbl2)
  tbl2[3,2, anchor = c(1,0)] = gWidgets2::glabel(expression(delta), container = tbl2)
  tbl2[3,3] <- wid$delta.start <- gWidgets2::gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 0, container = tbl2)
  if (autorange.delta) {
    gWidgets2::enabled(wid$delta.start) <- FALSE
  }
  tbl2[3,4, anchor = c(0,0)] = gWidgets2::glabel("<- between ->", container = tbl2)
  tbl2[3,5] <- wid$delta.stop <- gWidgets2::gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 10, container = tbl2)
  tbl2[3,6, anchor = c(0,0)] = gWidgets2::glabel("by", container = tbl2)
  tbl2[3,7] <- wid$delta.by <- gWidgets2::gspinbutton(from = 0, to = dim(lists)[1], by = 1, value = 1, container = tbl2)############
  if (autorange.delta) {
    gWidgets2::enabled(wid$delta.stop) <- FALSE
  }

                                        #layout of the third argument input area
  tbl3 <- gWidgets2::glayout(spacing = 10, container = a3frame)
  tbl3[2,2, anchor = c(1,0)] <- gWidgets2::glabel("min", container = tbl3)
  tbl3[2,3] <- wid$threshold <- gWidgets2::gedit(text = "50", width = 3, container = tbl3)
  tbl3[2,4, anchor = c(-1,0)] <- gWidgets2::glabel("%" ,container = tbl3)
  tbl3[3,2:4, anchor = c(1,0)] <- gWidgets2::glabel("for gray-shading of objects", container = tbl3)

                                        #create the delta slider for choosing the desired delta value (and to show the calculated data for the selected delta value)
  sldrframe <- gWidgets2::gframe("DELTA-SLIDER", pos = 0.5, container = maingroup)
                                        #the parameters 'from' and 'to' are dummy values because the slider will be updated after the calculation
  wid$delta.slider <- gWidgets2::gslider(from = 1, to = 2, by = 1, container = sldrframe, expand = TRUE, handler = function(h,...) { load.data() })
  gWidgets2::enabled(wid$delta.slider) <- FALSE

                                        #create four tabs for presenting the results of the calculation
  nb <- gWidgets2::gnotebook(container = maingroup, expand = TRUE)
  aggmapg <- gWidgets2::ggroup(horizontal = FALSE, container = nb, label = "Aggregation map", use.scrollwindow = TRUE)
  summtblg <- gWidgets2::ggroup(horizontal = FALSE, container = nb, label = "Summary table")
  venng <- gWidgets2::ggroup(horizontal = TRUE, container = nb, label = "Venn-diagram & Venn-table")
  gWidgets2::svalue(nb) <- 1 #set the first tab as the selected tab

                                        #create a status bar to inform the user of the current status
  wid$status <- RGtk2::gtkStatusbar()
  gWidgets2::add(maingroup, wid$status)
  info <- wid$status$getContextId("info")
  wid$status$push(info, "GUI started")

                                        #create a progress bar to show the progress of the calculation
  wid$progress <- RGtk2::gtkProgressBar(show = TRUE)
  gWidgets2::add(maingroup, wid$progress)	

                                        #update the allowed range for delta (based on the nu-value)
  update.deltarange <- function() {
    requireNamespace("RGtk2")
    requireNamespace("gWidgets2")
    requireNamespace("gWidgetsRGtk2")
    
                                        #check if allowed range for delta should be updated
    if (autorange.delta) {
      deltarange <- .determineDeltaRange(lists, as.numeric(gWidgets2::svalue(wid$N)), dim(lists)[2], as.numeric(gWidgets2::svalue(wid$v)))
      gWidgets2::svalue(wid$delta.start) <- deltarange[1]
      gWidgets2::svalue(wid$delta.stop) <- deltarange[2]
      gWidgets2::enabled(wid$delta.start) <- TRUE
      gWidgets2::enabled(wid$delta.stop) <- TRUE
      wid$status$push(info, expression(paste("Recalculated allowed range for delta with ",nu," =", ss), list(ss=gWidgets2::svalue(wid$v))))
    }
  }

                                        #update the calculated data set (aggmap, summary-table and Venn) in the three tabs
  load.data <- function(truncated.lists) {
    requireNamespace("RGtk2")
    requireNamespace("gWidgets2")
    requireNamespace("gWidgetsRGtk2")
    
                                        #update is only performed when delta slider is enabled
    if (gWidgets2::enabled(wid$delta.slider)) {
                                        #check if calculated data set can be shown in the GUI (no file = error in the calculation for this delta-value)
      if (file.exists(paste(directory, "/N" , current.arguments$val.N, "_L", current.arguments$val.L, "_delta", gWidgets2::svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))) {
                                        #load the calculated data set
        load(paste(directory, "/N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", gWidgets2::svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".Rdata", sep = ""))

                                        #first tab: show aggmap
        if (!is.null(wid$aggmap)) {
					#delete current aggmap image
          gWidgets2::delete(aggmapg, wid$aggmap)
        }
        wid$aggmap <- gWidgets2::gimage(paste(directory, "/aggmap_N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", gWidgets2::svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = aggmapg)


                                       #second tab: show summary-table
        if (!is.null(wid$sumtable)) {
					#delete current summary-table
          gWidgets2::delete(summtblg, wid$sumtable)
        }
        wid$sumtable <- gWidgets2::gtable(truncated.lists$summarytable, container = summtblg, expand = TRUE)#########################

                                        #third tab: show Venn-diagram and Venn-table (only if L <= 4)
        if (current.arguments$val.L <= 4) {			
          if (!is.null(wid$venndiag) & !is.null(wid$venntbl)) {
                                        #delete current Venn-diagram and Venn-table
            gWidgets2::delete(venng, wid$venndiag)
            gWidgets2::delete(venng, wid$venntbl)
          }
          if (!is.null(wid$nodraw) & !is.null(wid$nodrawimage)) {
                                        #delete the icon and message
            gWidgets2::delete(venng, wid$nodrawimage)
            gWidgets2::delete(venng, wid$nodraw)
          }
          wid$venndiag <- gWidgets2::gimage(paste(directory, "/venn_N", current.arguments$val.N, "_L", current.arguments$val.L, "_delta", gWidgets2::svalue(wid$delta.slider), "_v", current.arguments$val.v, "_thrshld", current.arguments$val.threshold, ".png", sep = ""), container = venng)
          mydata <- as.data.frame(truncated.lists$venntable)
          mydata <- matrix(unlist(mydata), byrow=FALSE,ncol=2)
          wid$venntbl <- gWidgets2::gtable(mydata, multiple=TRUE, container = venng, expand = TRUE, drop=FALSE, index=TRUE)

        } else {
					#output message that Venn-diagram and Venn-table are not shown for L > 4
          if (is.null(wid$nodrawimage) & is.null(wid$nodraw)) {
            wid$nodrawimage <- gWidgets2::gimage("info", dirname = "stock", container = venng)
            wid$nodraw <- gWidgets2::glabel("There are more than four input lists.\nVenn-diagram and Venn-table are only drawn to a maximum of four input lists.", container = venng)
          }
        }# end for third tab if
  
      } else {
        gWidgets2::gmessage("There is no data to show in the GUI!\nThis is caused by an error in the calculation.\nYou may have to adjust the arguments!", title = "Loading data failed", icon = "error")
      }
    }
  }
  
                                        #calculates aggmap, summary-table and Venn for each delta in the given range for delta
  calculate.all.truncationlists <- function() {
    requireNamespace("RGtk2")
    requireNamespace("gWidgets2")
    requireNamespace("gWidgetsRGtk2")
    
                                        #check if entered values are valid
    if ((as.numeric(gWidgets2::svalue(wid$delta.start)) < 0) || (as.numeric(gWidgets2::svalue(wid$delta.stop)) < 0) || (gWidgets2::svalue(wid$delta.start) > gWidgets2::svalue(wid$delta.stop))) {
      gWidgets2::gmessage("Invalid range for delta entered.\nPlease check your input.", title = "Invalid range", icon = "error")
      return()
    }
    temp.threshold <- as.numeric(gWidgets2::svalue(wid$threshold))
    if (is.na(temp.threshold) || (temp.threshold < 0) || (temp.threshold > 100)) {
      gWidgets2::gmessage("The entered threshold is not a valid number!\nPlease enter a threshold between 0 and 100.", title = "Threshold invalid", icon = "error")
      return()
    }

    wid$status$push(info, "Begin calculation for the selected range of delta")
    RGtk2::gtkProgressBarSetFraction(wid$progress, 0) #set progress-bar to start-value
    RGtk2::gtkMainIterationDo(FALSE)

                                        #disable delta-slider to prevent sliding while calculating
    gWidgets2::enabled(wid$delta.slider) <- FALSE

                                        #calculate data set for each delta in the given range for delta
    temp.tocalc <- length(seq(as.numeric(gWidgets2::svalue(wid$delta.start)),as.numeric(gWidgets2::svalue(wid$delta.stop)),by=as.numeric(gWidgets2::svalue(wid$delta.by))))
    temp.loopcounter <- 0
    wid$error <- FALSE

	
       					 #draw the Delta-plot and save it as pdf-image in the destination directory
        Mdelta = suppressWarnings(deltaplot(lists, directory=directory))

        save(Mdelta, file=paste(directory, "/Mdelta.rdata", sep=""))
        for (aname in names(Mdelta))
        {
	write.table(Mdelta[[aname]], file=paste(directory, "/Mdelta",aname,".txt", sep=""), sep="\t", row.names=FALSE)
        }

        rm(Mdelta)

    for (i in seq(as.numeric(gWidgets2::svalue(wid$delta.start)),as.numeric(gWidgets2::svalue(wid$delta.stop)),by=as.numeric(gWidgets2::svalue(wid$delta.by)))) {
      RGtk2::gtkMainIterationDo(FALSE)
                                        #try to calculate data set for current delta and selected N
      tryCatch({
        truncated.lists <- calculate.maxK(lists[1:as.numeric(gWidgets2::svalue(wid$N)),], dim(lists)[2], i, as.numeric(gWidgets2::svalue(wid$v)), as.numeric(gWidgets2::svalue(wid$threshold)))
                                        #save calculated truncated lists for current delta in the destination directory
        save(truncated.lists, file = paste(directory, "/N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".Rdata", sep = ""))

        #draw the aggmap and save it as pdf-image in the destination directory
	      if(is.null(truncated.lists)){max.length=10}else{max.length <- min(length(truncated.lists$comparedLists[[1]]), as.numeric(gWidgets2::svalue(wid$N)))}

      	png(filename = paste(directory, "/aggmap_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".png", sep = ""), width = 870, height = 60+28*max.length, res=aggmap.res, type="cairo-png")
        aggmap(truncated.lists)
        dev.off()
        
        pdf(file = paste(directory, "/aggmap_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".pdf", sep = ""), width = 9, height = 0.8+0.29*max.length)
        aggmap(truncated.lists)
        dev.off()
        

                                        #draw the Venn-diagram and save it as pdf-image in the destination directory (only if L is between 2 and 4 & if there is a j0 estimated)
        if((dim(lists)[2] >= 2) & (dim(lists)[2] <= 4) & !is.null(truncated.lists)) {
          pdf(file = paste(directory, "/venn_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".pdf", sep = ""), width = venndiag.pdf.size[1], height = venndiag.pdf.size[2], bg = "transparent")
          venn(truncated.lists$vennlists)
          dev.off()


          png(filename = paste(directory, "/venn_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res=venndiag.res)
          venn(truncated.lists$vennlists)
          dev.off()
        }else{
          pdf(file = paste(directory, "/venn_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".pdf", sep = ""), width = venndiag.pdf.size[1], height = venndiag.pdf.size[2], bg = "transparent")
          plot(1,1, type="n", axes=FALSE)
          text(1,1,"No overlap on selected parameters")
          dev.off()

 	      png(filename = paste(directory, "/venn_N", as.numeric(gWidgets2::svalue(wid$N)), "_L", dim(lists)[2], "_delta", i, "_v", as.numeric(gWidgets2::svalue(wid$v)), "_thrshld", as.numeric(gWidgets2::svalue(wid$threshold)), ".png", sep = ""), width = venndiag.size[1], height = venndiag.size[2], bg = "transparent", res=venndiag.res)
          plot(1,1, type="n", axes=FALSE)
          text(1,1,"No overlap on selected parameters")
          dev.off()
        }
      }, error = function(e) {
          print(e)
          traceback()
                                        #catch errors while calculating the data set
        wid$status$push(info, paste("Error in the current calculation for delta =", i))
        wid$error <- TRUE
        if(!override.errors) {
          return()
        }
      })# end tryCatch
      if (wid$error & !override.errors) {
                                        #exit the loop
        break()
      }

                                        #update progress bar
      temp.loopcounter <- temp.loopcounter + 1
      temp.fraction <- (1 / temp.tocalc) * temp.loopcounter
      RGtk2::gtkProgressBarSetFraction(wid$progress, temp.fraction)
      RGtk2::gtkMainIterationDo(FALSE)
    }

                                        #if errors are overridden, the calculation continues (but there is no data to view in the GUI for this delta value)
    if (wid$error) {
      if(override.errors) {
        wid$status$push(info, "Finished calculation of data sets for range of delta (with errors occurred)")
      } else {
        gWidgets2::enabled(wid$delta.slider) <- FALSE
        return()
      }
      wid$error <- FALSE
    } else {		
      wid$status$push(info, "Finished calculation of data sets for range of delta")
    }

                                        #update the range of the delta-slider
    wid$delta.slider[] <- seq(as.numeric(gWidgets2::svalue(wid$delta.start)), as.numeric(gWidgets2::svalue(wid$delta.stop)), by = as.numeric(gWidgets2::svalue(wid$delta.by)))

                                        #set the start position of the slider and enable slider
    if (temp.tocalc > 2) {
      gWidgets2::svalue(wid$delta.slider) <- as.numeric(gWidgets2::svalue(wid$delta.start)) + round((temp.tocalc - 1) / 2)
    } else {
      gWidgets2::svalue(wid$delta.slider) <- as.numeric(gWidgets2::svalue(wid$delta.start))
    }
    gWidgets2::enabled(wid$delta.slider) <- TRUE

                                        #save the current arguments
                                        #(this prevents errors when sliding through the delta-values - otherwise changing the arguments in the GUI could cause errors when showing the calculated data)
    current.arguments$val.N <- as.numeric(gWidgets2::svalue(wid$N))
    current.arguments$val.L <- dim(lists)[2]
    current.arguments$val.v <- as.numeric(gWidgets2::svalue(wid$v))
    current.arguments$val.threshold <- as.numeric(gWidgets2::svalue(wid$threshold))
    
                                        #show the data in the tabs
    load.data()
  }#end calculate.all.truncationlists
  
                                        #show the GUI
  gWidgets2::visible(win) <- TRUE
}


.calculateVennValues <- function(genetable, L) {
                                        #check if L is between 2 and 4 (only in this range a Venn-table will be calculated)
  if ((L >= 2) & (L <= 4)) {
                                        #create new lists with the object-symbols or names as list-entries (currently only the distances are in the columns because the summary-table is used)
    newlists <- list()
    for (x in 2:(L + 1)) {
      currentlist <- c()
      for (y in 1:length(genetable[,1])) {
                                        #check if object-symbol is in the current list (if the distance is NA, the object is not in the current list)
        if (!is.na(genetable[,x][y])) {
          currentlist <- append(currentlist, genetable[,1][y], after = y)
        }
      }
      newlists[[colnames(genetable)[x]]] <- currentlist
    }
    
    if(L == 2) {
                                        #create Venn-table for two input lists
      venntable <- matrix(ncol = 2, nrow = 1)
      colnames(venntable) <- c("intersection", "objects")
      venntable[1,1] <- paste(names(newlists)[1], names(newlists)[2], sep = "_")
      venntable[1,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
      return(list("vennlists" = newlists, "venntable" = venntable))
    }else if(L == 3) {
                                        #create Venn-table for three input lists
      venntable <- matrix(ncol = 2, nrow = 4)
      temp.L1_L2 <- intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]]) #to calculate L1_L2_L3
      temp.rownames <- c(paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[2], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[3], sep = "_"),
                         paste(names(newlists)[2], names(newlists)[3], sep = "_"))
      venntable[,1] <- temp.rownames
      colnames(venntable) <- c("intersection", "objects")
      venntable[1,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[3]]])) #L1_L2_L3
      venntable[2,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
      venntable[3,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[3]]])) #L1_L3
      venntable[4,2] <- toString(intersect(newlists[[names(newlists)[2]]], newlists[[names(newlists)[3]]])) #L2_L3
      return(list("vennlists" = newlists, "venntable" = venntable))
    }else if(L == 4) {
                                        #create Venn-table for four input lists
      venntable <- matrix(ncol = 2, nrow = 9)
      temp.L1_L2 <- intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])
      temp.L3_L4 <- intersect(newlists[[names(newlists)[3]]], newlists[[names(newlists)[4]]])
      temp.rownames <- c(paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], names(newlists)[4], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[2], names(newlists)[3], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[2], names(newlists)[4], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[3], names(newlists)[4], sep = "_"),
                         paste(names(newlists)[2], names(newlists)[3], names(newlists)[4], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[2], sep = "_"),
                         paste(names(newlists)[2], names(newlists)[3], sep = "_"),
                         paste(names(newlists)[3], names(newlists)[4], sep = "_"),
                         paste(names(newlists)[1], names(newlists)[4], sep = "_"))
      venntable[,1] <- temp.rownames
      colnames(venntable) <- c("intersection", "objects")
      venntable[1,2] <- toString(intersect(temp.L1_L2, temp.L3_L4)) #L1_L2_L3_L4
      venntable[2,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[3]]])) #L1_L2_L3
      venntable[3,2] <- toString(intersect(temp.L1_L2, newlists[[names(newlists)[4]]])) #L1_L2_L4
      venntable[4,2] <- toString(intersect(newlists[[names(newlists)[1]]], temp.L3_L4)) #L1_L3_L4
      venntable[5,2] <- toString(intersect(newlists[[names(newlists)[2]]], temp.L3_L4)) #L2_L3_L4
      venntable[6,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[2]]])) #L1_L2
      venntable[7,2] <- toString(intersect(newlists[[names(newlists)[2]]], newlists[[names(newlists)[3]]])) #L2_L3
      venntable[8,2] <- toString(intersect(newlists[[names(newlists)[3]]], newlists[[names(newlists)[4]]])) #L3_L4
      venntable[9,2] <- toString(intersect(newlists[[names(newlists)[1]]], newlists[[names(newlists)[4]]])) #L1_L4
      return(list("vennlists" = newlists, "venntable" = venntable))
    }			
  } else {
    return(list("vennlists" = NA, "venntable" = NA))
  }
}





.determineDeltaRange <- function(lists, N, L, v) {
  start.delta <- 0
  stop.delta <- 0
  start.delta.found <- FALSE
  
                                        #iterate through all possible delta-values
  for (d in 1:N) {
                                        #check if truncated lists are calculatable with current delta
                                        #if not calculatable, then the current delta value cannot be applied on the lists (aggmap cannot be drawn)
    if (!.isCalculatable(lists, d, L, v)) {
      if (start.delta.found) {
        return(c(start.delta, d-1))
      }
    } else {
      if (!start.delta.found) {
        start.delta <- d
        start.delta.found <- TRUE
      } else {
        stop.delta <- d
      }
    }
  }
  
  return(c(start.delta, stop.delta))
}

.prepareDeltaplot <- function(lists, i, j, deltas, Mdelta,xxs)
{		
	
	## zero count calculation and deltaplot for LiLj (in one window)
    Mdelta.temp = data.frame(Object=c(as.character(lists[,i]), "#zeros"), L1=c(c(1:nrow(lists)), NA), L2 = c(match(lists[,i],lists[,j]), NA))
    names(Mdelta.temp)[2:3] = names(lists)[c(i,j)]
    xx = c()
    for (d in deltas)
    {	
        a = prepare.idata(lists[,c(i,j)],d=d)
        x = table(as.numeric(a$Idata))['0']
        xx = c(xx,x)
        Mdelta.temp[,paste("delta_",d)] = c(a$Idata, x)
    }# end for d
	Mdelta[[paste(names(lists)[i],"_",names(lists)[j], sep="")]] = Mdelta.temp
    xxs[[paste(names(lists)[i],"_",names(lists)[j], sep="")]] = xx  ##saving xx for plotting single deltaplot with subplot in the corner
	par(mar=c(5,5,1,1))
    plot(deltas,xx, xlab=expression(delta), ylab="# of 0's", las=1,cex.axis=0.7, main=paste(names(lists)[i]," vs ",names(lists)[j], sep=""))
	
	MdeltaXxs = list(Mdelta=Mdelta, xxs=xxs)
	return(MdeltaXxs)
}

###function that generates Delta-plot and Delta-matrix
##if subset.plotted is NA no subplots are created
deltaplot<-function(lists, deltas=NULL, subset.lists=NULL, subplot = FALSE, perc.subplot=50, directory = NULL)
{

  if (is.null(deltas)) {
    warning(paste("The vector of delta values was not specified, using range from 0 to ", floor(nrow(lists)*0.25), "by 1"), "\n")
    deltas=c(0:(floor(nrow(lists)*0.25)))
  }
  if (max(deltas)>nrow(lists)) {
    warning(paste("The maximum for delta you specified is larger than the number of objects in your lists. The vector of deltas was truncated,","\n","so that max(delta)<",nrow(lists)*0.25), "\n")
    deltas=deltas[which(deltas<nrow(lists)*0.25)]
  }
  if(is.null(names(lists)) | any(names(lists)=="")){
    names(lists) <- paste("L",1:ncol(lists),sep="")
    warning(paste("List names not given or incorrect. Replaced by 'L1', 'L2',... L",ncol(lists),'\n'))
  }

  if(!is.null(subset.lists)){
    lists = lists[1:subset.lists,] ## takes only specified subset of input list  
    if (max(deltas)>nrow(lists)) {
    	cat(paste("The maximum for delta you specified is larger than the number of objects to be analyzed in your lists.","\n","The vector of deltas was truncated so that max(delta)<",nrow(lists)*0.25), "\n")
    	deltas=deltas[which(deltas<nrow(lists)*0.25)]
     }

  }
  if(!is.null(directory) && (!file.exists(directory) || !file.info(directory)$isdir)){
    stop("The given directory does not exist")
  }
  
  ## calculate zero count for each pair of lists  
  MdeltaXxs=list(Mdelta=list(), xxs = list())
  
  for (i in 1:(ncol(lists)-1))
    {
      for (j in (i+1):ncol(lists))
        {
            ##create a pdf
			if (!is.null(directory)) {
				pdf(file=paste(directory,'/deltaplotL',i,"-" ,j,'.pdf',sep=''), width=16)
				par(mfrow=c(1,2))
				MdeltaXxs=.prepareDeltaplot(lists,i,j,deltas,Mdelta=MdeltaXxs$Mdelta, xxs=MdeltaXxs$xxs)
				MdeltaXxs=.prepareDeltaplot(lists,j,i,deltas,Mdelta=MdeltaXxs$Mdelta, xxs=MdeltaXxs$xxs)
				dev.off()
                              }else {
                                    ##the next line need to be fixed in order to be more generic (dev.new did not work so far)
##                                if( .Platform$OS.type =="unix" ) { X11() } else { windows() }
                                dev.new()
                                par(mfrow=c(1,2))
                                MdeltaXxs=.prepareDeltaplot(lists,i,j,deltas,Mdelta=MdeltaXxs$Mdelta, xxs=MdeltaXxs$xxs)
                                MdeltaXxs=.prepareDeltaplot(lists,j,i,deltas,Mdelta=MdeltaXxs$Mdelta, xxs=MdeltaXxs$xxs)
				}#end plotpdf
		}# end for j
    }# end for i

  if(subplot){
    ## deltaplot with subplot in the top right corner:
	xxs = MdeltaXxs$xxs
    for (i in 1:ncol(lists))
      {
        for (j in 1:ncol(lists))
          {
            if (i!=j){
              ##the next line need to be fixed in order to be more generic (dev.new did not work so far)
                ##if( .Platform$OS.type =="unix" ) { X11() } else { windows() }
                ##print pdf
                if(!is.null(directory)){
                pdf(file=paste(directory,'/subplotL',i,"-" ,j,'.pdf',sep=''), width=16)
              } else {
                dev.new()
              }
              par(mar=c(5,5,1,1))
              plot(deltas,as.vector(xxs[[paste(names(lists)[i],"_",names(lists)[j], sep="")]]), xlab=expression(delta), ylab="# of 0's", las=1,cex.axis=0.7, main=paste(names(lists)[i]," vs ",names(lists)[j], sep=""))			
              extremes = par("usr")
              dimen = par("pin")					   
              subplot(plot(deltas[1:((perc.subplot/100)*length(deltas))],as.vector(xxs[[paste(names(lists)[i],"_",names(lists)[j], sep="")]])[1:((perc.subplot/100)*length(deltas))], xlab="", ylab="", las=1, cex.axis=0.7) , extremes[2], extremes[4], size = c(dimen[1]*0.5, dimen[2]*0.4),hadj=1, vadj=1, pars=list(col="black", mar=c(5,5,1,1)))
              if(!is.null(directory))
                dev.off()
            }
          }
      }
  }
  return(MdeltaXxs$Mdelta) # returns vectors of zeros and ones

}#end deltaplot

## Aggregation map (aggmap) for the graphical representation of more than two top-k lists obtained from the method by Hall and Schimek (2012) - it applies Paul Murrell's grid package
## Input variables 
# L - number of lists
# N - maximal number of objects
# v - value of tuning parameter nu used for j0 estimation
# res.temp - data frame, result from the Hall and Schimek (2012) algorithm for all pairings of lists, 4 columns, the first and second column contain the list names that were compared, the third column contains the selected nu, and the fourth column the estimated j0
# The lists data frame contains the complete set of L ordered input lists - column names are obligatory

aggmap <- function(truncated.lists) {
N=truncated.lists$N
v=truncated.lists$v
d=truncated.lists$d
lists=truncated.lists$lists
threshold=truncated.lists$threshold
L=truncated.lists$L
                                       #delta_symbol = substitute(delta)
                                        #nu_symbol = expression(nu)
      wid <- 1
      hei <- 1
                                        #create list of background colours to shade the distance-polygons  
      red.orange.white = c( "#BC0000", "#FF3700", "#FF6E00", "#FFA500", "#FFC300", "#FFE100", "#FFFF00", "#FFFF54", "#FFFFA9", "#FFFFFF")
 
   if(is.null(truncated.lists)){max.length=10}else{max.length <- min(length(truncated.lists$comparedLists[[1]]), N)}


      grid.newpage()

      #pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = c(0.1, 0.9)), name="basic"))# creating the basic viewport with two rows and one column

      pushViewport(viewport(layout = grid.layout(2, 1, widths = c(1, 1), heights = unit(c(40,max.length*20), "points")), name="basic"))# creating the basic viewport with two rows and one column


                         #output the header (the parameters and the colour gradient)

      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, name="header")) #writing in first row of the parent region basic

 	     pushViewport(viewport(layout = grid.layout(2,2, widths = c(1, 1,1,1), heights = c(1, 1,1,1)), name="HeaderSplit")) # in this first row and first column of the parent region basic we are defining another split with two rows and two columns
		      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1, name="headerLeft1")) #writing in first column and first row of the header region 
		    	  grid.text(substitute(paste(hat(k)[max],"= ", kmax,", N = ", N, ", L = ", L, sep = ""), list(N=N, L=L, kmax=max.length)), x = 0.5, y = 0.4)
		      upViewport() # coming back to the HeaderSplit

		      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1, name="headerLeft2")) #writing in first column and second row of the header region 
 		      	grid.text(substitute(paste(delta, "= ", a, ", ",nu," = ", b, ", threshold = ", f, "%", sep = ""), list(a=d, b=v, f=threshold)), x = 0.5, y = 0.6)
    		      upViewport() # coming back to the HeaderSplit


     		      pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2, name="headerRight1")) #writing in the second column and first row of the header region 
      			pushViewport(plotViewport(margins = c(1, 1, 1, 1))) #plotting in the headerRight
				grid.points(seq(0.2, 0.8, length =  10), rep(1, 10), pch = 15, gp = gpar(col = red.orange.white))
      			upViewport() # coming back to headerRight1 level
		      upViewport() # coming back to the HeaderSplit
      
     		      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2, name="headerRight2")) #writing in the second column and second row of the header region 
      			grid.text("close", x = 0.15, y = 0.8, gp = gpar(fontsize = 8, col = "black"))
     			grid.text("RANK POSITION", x = 0.5, y = 0.5, gp = gpar(fontsize = 7, col = "black"))
    			grid.text("distant", x = 0.85, y = 0.8, gp = gpar(fontsize = 8, col = "black"))
      		      upViewport() #coming back to the HeaderSplit 
     	     upViewport() #coming back to the header
      upViewport() #coming back to the basic

  if (!is.null(truncated.lists)) 
    {
                                         #get the necessary data
      compared.lists <- truncated.lists$comparedLists
      info <- truncated.lists$info
      grayshaded.lists <- truncated.lists$grayshadedLists

                                        #set the layout of the heatmap part
      pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))	 #writing in the second row of the parent region basic
 
      pushViewport(viewport(layout = grid.layout(max.length + 1, 1+truncated.lists$Ntoplot, widths = rep(1, 1+truncated.lists$Ntoplot) * (max.length), heights = rep(1, 1 + truncated.lists$Ntoplot) * (max.length)))) #dividing the second row into rows and columns
      
                                        #output row-numbers
      for (g in c(1:max.length)) {
        pushViewport(viewport(width = wid, height = hei, layout.pos.col = 1, layout.pos.row = g + 1))
        grid.text(g, gp = gpar(fontsize = 8, col = "black"))
        popViewport()
      }
      
      h.sizes <- c(L:2)
      
                                        #output the lists
      for (i in 1:length(compared.lists)) {
                                        #get the name of the list and output the list name
        listname <- info[2, i]
        pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = 1))
        grid.text(listname, gp = gpar(fontsize = 8, col = "black"))
        popViewport()
        
                                        #output the list (check output for reference or truncated list)
        if (info[3,i] == "R") {
                                        #output reference list
          for (y in 1:min(length(compared.lists[[i]]), max.length)) {
            pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = y + 1))				
                                          #check if the object name has to be gray-shaded
            bg <- "white"
            if (!is.na(grayshaded.lists[[i]][y])) { if(grayshaded.lists[[i]][y]==1) {
              bg <- "grey"
            }}# end for y
            grid.rect(gp = gpar(col = "black", fill = bg))
            grid.text(compared.lists[[i]][y], gp = gpar(fontsize = 8, col = "black", fontface = "bold"))
            popViewport()
          }# end if info
        } else if (info[3,i] == "T") {
            #output truncated list                            
          for (y in 1:min(length(compared.lists[[i]]), max.length)) {
            distance <- compared.lists[[i]][y]
            pushViewport(viewport(width = wid, height = hei, layout.pos.col = i + 1, layout.pos.row = y + 1))
            if (!is.na(distance)) {
              grid.polygon(x = c(1, 1, 0), y = c(0, 1, 0), gp = gpar(col = "black", fill = red.orange.white[abs(distance)/N*10 + 1]))
              grid.text(-distance, x = 0.9, y = 0.4, gp = gpar(fontsize = 6)) #  cex = 0.7
            } else {
              grid.polygon(x = c(1, 1, 0), y = c(0, 1, 0), gp = gpar(col = "black", fill = "white"))
              grid.text("NA", x = 0.9, y = 0.4, gp = gpar(fontsize = 6)) #cex = 0.7
            }
            bg <- "white"
            if (!is.na(grayshaded.lists[[i]][y])) { if(grayshaded.lists[[i]][y]==1) {
              bg <- "grey"
            }}#end for else if
            grid.polygon(x = c(0, 1,0), y = c(1, 1, 0), gp = gpar(col = "black", fill = bg))
            popViewport()				
          }
        }
      }#end for i
      upViewport()# returning to the heatmap part
      upViewport()# returning to the basic viewport
    } else {                                        #set the layout of the heatmap
            pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))	
            pushViewport(viewport(layout = grid.layout(1, 2, widths = c(1), heights = c(1,1))))
            grid.text(substitute(paste("On selected ", delta ,", no overlap found", sep="")))
            upViewport()
            upViewport()
            upViewport()  } #end for is.null 
}#end for aggmap


.isCalculatable <- function(lists, d, L, v) {
  res.j0.temp <- j0.multi(lists, d, v)
  res.temp <- as.matrix(res.j0.temp$L)
  
  temp2 = c()
                                        #calculate the reference lists
  for (i in c(1:L)) {
    if (i > 1) {
      res.temp.temp <- res.temp[-c(which(!is.na(match(res.temp[,1], names(temp2)))), which(!is.na(match(res.temp[,2], names(temp2))))),]
    } else {
      res.temp.temp <- res.temp
    }
    temp <- tapply(as.numeric(res.temp.temp[,4]), res.temp.temp[,1], function(x) max(x, na.rm = TRUE))
    which.max.temp <- which.max(temp)
    temp2 <- c(temp2, temp[which.max.temp])
  }
  
  lists.to.remove <- ""
  truncation.points <- c()
  
                                        #calculate the truncated lists for each block (for each reference list)
  for (first.list.name in names(temp2)) {
    res.temp.selection <- res.temp[which(res.temp[,1] == first.list.name), ]
    res.temp.selection <- res.temp.selection[which(is.na(match(res.temp.selection[,2], lists.to.remove))),]
    if (is.matrix(res.temp.selection)) {
      ranked.list.names <- c(first.list.name, res.temp.selection[order(res.temp.selection[,4], decreasing = TRUE),2])
      lists.rank <- rank(ranked.list.names)
    } else {
      ranked.list.names <- c(first.list.name, res.temp.selection[2])
      lists.rank <- rank(ranked.list.names)
    }
    
                                        #get the truncation points for the truncated lists in the current block
    for(currennt.truncated.list in ranked.list.names[2:length(ranked.list.names)]) {
      truncation.points <- c(truncation.points, as.numeric(res.temp[res.temp[,1] == first.list.name & res.temp[,2] == currennt.truncated.list,4]))
    }
    
    lists.to.remove <- c(lists.to.remove, first.list.name)
  }
  
                                        #check if a truncation point is missing - if true, then the delta-value cannot be applied on the lists (aggmap is not drawn)
  if(NA %in% truncation.points) {
    return(FALSE)
  }else {
    return(TRUE)
  }
}

