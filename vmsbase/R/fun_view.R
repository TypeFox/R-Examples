

bathyview <- function(bathy, xrange, yrange)
{
  isob <- c(-20, -50, -100, -200, -500, -1000, -1500)
  lon <- unique(as.numeric(rownames(bathy$data)))
  lat <- unique(as.numeric(colnames(bathy$data)))
  icol <- rgb(0,0, seq(250,20, len = length(isob)), maxColorValue = 255)
  contour(lon, lat, bathy$data, add=TRUE,
          #lwd=0.5,
          col=icol,
          ylim = yrange, xlim = xrange,
          levels = c(-20, -50, -100, -200, -500, -1000, -1500),
          drawlabel=TRUE) 
}

plot_comp <- function(sel_opt, obj_lst){
  if(length(obj_lst) == 0){
    plot.new()
  }else{
    switch(sel_opt ,
           Absolute = plot_coun(obj_lst),
           Relative = plot_dens(obj_lst),
           IDs = plot_IDs(obj_lst),
           Stats = plot_statot(obj_lst),
           Box = plot_box(obj_lst))
  }
}


plot_coun <- function(obj_list){
  op <- par(no.readonly = TRUE)
  plot.new()
  in_dat <- data.frame(date = character(), sorg = character(), stringsAsFactors = FALSE)
  for(j in 1:length(obj_list)){
    tmp_dat <- data.frame(date = as.numeric(obj_list[[j]]$vessDate), sorg = obj_list[[j]]$sourceName, stringsAsFactors = FALSE)
    if(nrow(tmp_dat) > 0){
      in_dat <- rbind(in_dat, tmp_dat)
    }
  }
  the_plot <- ggplot(data = as.data.frame(in_dat), aes_string(x = "date", fill = "sorg")) + xlab("Date") + ylab("Number of Pings") + scale_x_chron(format = "%m - %Y") + stat_bin(binwidth = 1) + facet_grid(sorg ~ .) + ggtitle("Absolute Frequency of Signal") 
  print(the_plot)
  par(op)
}


plot_statot <- function(obj_list){
  op <- par(no.readonly = TRUE)
  plot.new()
  in_dat <- data.frame(typ = character(), val = character(), sorg = character(), stringsAsFactors = FALSE)
  for(j in 1:length(obj_list)){
    tmp_dat <- data.frame(typ = c("Pings","Days","Vessels"), val = c(as.numeric(obj_list[[j]]$n_ping),  as.numeric(obj_list[[j]]$n_day), as.numeric(obj_list[[j]]$n_vess)), sorg = obj_list[[j]]$sourceName, stringsAsFactors = FALSE)
    if(nrow(tmp_dat) > 0){
      in_dat <- rbind(in_dat, tmp_dat)
    }
  }
  the_plot <- ggplot(data = as.data.frame(in_dat), aes_string(x = as.factor("typ"), y = "val", colour = "sorg")) + xlab("") + ylab("")
  the_plot <- the_plot + geom_text(aes_string(label = "val", ymax = "val"), size = 4, position = position_dodge(width=0.9), vjust=0)
  the_plot <- the_plot + theme_minimal()
  the_plot <- the_plot + facet_wrap( ~ typ, scales = "free")
  the_plot <- the_plot + geom_bar(aes_string(fill = "sorg"), position = "dodge", stat = "identity")
  print(the_plot)
  par(op)
}


plot_dens <- function(obj_list){
  op <- par(no.readonly = TRUE)
  plot.new()
  in_dat <- data.frame(date = character(), sorg = character(), stringsAsFactors = FALSE)
  for(j in 1:length(obj_list)){
    tmp_dat <- data.frame(date = as.numeric(obj_list[[j]]$vessDate), sorg = obj_list[[j]]$sourceName, stringsAsFactors = FALSE)
    if(nrow(tmp_dat) > 0){
      in_dat <- rbind(in_dat, tmp_dat)
    }
  }
  print(ggplot(data = as.data.frame(in_dat), aes_string(x = "date", colour = "sorg")) + geom_histogram(aes_string(y = "..density..", fill = "sorg"), alpha=.3, position="dodge", binwidth = 7) + xlab("Date") + ylab("Relative Frequency") + scale_x_chron(format = "%m - %Y") + ggtitle("Relative Frequency of Signal"))
  #   the_plot <- ggplot(data = as.data.frame(in_dat), aes(x = date, colour = sorg)) + geom_histogram(aes(y = ..density.., fill=sorg), alpha=.3, position="dodge", binwidth = 7) + xlab("Date") + ylab("Relative Frequency") + scale_x_chron(format = "%m - %Y") + ggtitle("Relative Frequency of Signal")
  #   print(the_plot)
  par(op)
}


plot_IDs <- function(obj_list){
  grid.newpage()
  in_dat <- list()
  if(length(obj_list) < 6){
    for(j in 1:length(obj_list)){
      in_dat <- c(in_dat, list(obj_list[[j]]$vessIds))
      names(in_dat)[j] <- obj_list[[j]]$sourceName}
    #   sto_zz <- venn.diagram(in_dat, filename = NULL, col = rainbow(length(obj_list)), alpha = 0.5, sub.cex = 0.7, main.cex = 0.8, cat.pos = 0)
    sto_zz <- VennDiagram::venn.diagram(in_dat, filename = NULL, alpha = 0.5, sub.cex = 0.7, main.cex = 0.8, cat.pos = 0)
    grid.draw(sto_zz)
  }else{
    plot.new()
    title("ID Diagram not available with more than 5 sources")
  }
}


plot_box <- function(obj_list){
  op <- par(no.readonly = TRUE)
  plot.new()
  tmp_dat <- list()
  for(j in 1:length(obj_list)){
    tmp_lst <- list(obj_list[[j]]$n_pida)
    names(tmp_lst) <- obj_list[[j]]$sourceName
    tmp_dat <- c(tmp_dat, tmp_lst)
  }
  boxplot(tmp_dat, outline=F, ylab = "Pings per day", main = "Distributions for each source")
  par(op)
}
########################################
########################################


########################################
#            Ping Viewer
########################################

pingview <- function(vessel, bathy, xrange, yrange)
{
  pal_col <- add.alpha(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"), 0.5)
  leg_col <- add.alpha(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF"), 0.9)
  span <- 0.25
  if(diff(xrange) == 0){
    xrange <- c(vessel[1,"LON"]-0.01, vessel[1,"LON"]+0.01)
  }
  if(diff(yrange) == 0){
    yrange <- c(vessel[1,"LAT"]-0.01, vessel[1,"LAT"]+0.01)
  }
  map2plot <- try(maps::map("worldHires",ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
  if(class(map2plot) != "try-error")
  {
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    if("SOU" %in% colnames(vessel)){
      points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = pal_col[as.numeric(as.factor(vessel[,"SOU"]))], pch = 20, cex = 0.7)
      legend("bottomright", sort(unique(vessel[,"SOU"])), col = leg_col[1:length(unique(vessel[,"SOU"]))], pch = 20, cex = 0.9, title = "Source", bg = "white")
    }else{
      points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = pal_col[1], pch = 20, cex = 0.7)
    }
    maps::map.axes()
    title(main = paste("Vessel ", vessel[1,1]," with ", nrow(vessel)," pings" , sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
  }else{
    while(class(map2plot) == "try-error")
    {
      span <- span + 0.25
      xrange <- extendrange(x = vessel["LON"], f = span)
      yrange <- extendrange(x = vessel["LAT"], f = span)
      map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
      maps::map.axes()
    }
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    if("SOU" %in% colnames(vessel)){
      points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = pal_col[as.numeric(as.factor(vessel[,"SOU"]))], pch = 20, cex = 0.7)
      legend("bottomright", sort(unique(vessel[,"SOU"])), col = leg_col[1:length(unique(vessel[,"SOU"]))], pch = 20, cex = 0.9, title = "Source", bg = "white")
    }else{
      points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = pal_col[1], pch = 20, cex = 0.7)
    }
    maps::map.axes()
    title(main = paste("Vessel ", vessel[1,1]," with ", nrow(vessel)," pings" , sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    span <- 0.25
  }
}
########################################
########################################



########################################
#            Track Viewer
########################################
trackview <- function(vessel, bathy = "", xrange, yrange)
{
  span <- 0.25
  if(diff(xrange) == 0){
    xrange <- c(vessel[1,"LON"]-0.01, vessel[1,"LON"]+0.01)
  }
  if(diff(yrange) == 0){
    yrange <- c(vessel[1,"LAT"]-0.01, vessel[1,"LAT"]+0.01)
  }
  map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
  if(class(map2plot) != "try-error")
  {
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel ", vessel[1,1]," with ", length(unique(vessel[,"T_NUM"]))," tracks" , sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    colsca <- rainbow(max(vessel[, "T_NUM"]), alpha = 0.5)
    textn <- 0
    for (i in 1:(nrow(vessel)-1))
    {
      trackn <- vessel[i,"T_NUM"]
      if(trackn != 0)
      {
        if(trackn == vessel[i+1,"T_NUM"])
        {
          lines(cbind(c(vessel[i,"LON"],vessel[i+1,"LON"]), c(vessel[i,"LAT"], vessel[i+1,"LAT"])), col = colsca[vessel[i,"T_NUM"]], lty = "dashed", cex = 1.5)
          if((vessel[i,"W_HARB"] != 1) & (trackn != textn))
          {
            text(vessel[i,"LON"], vessel[i,"LAT"], col = colsca[vessel[i,"T_NUM"]], labels = vessel[i,"T_NUM"])
            textn <- trackn
          }
        }
      }
    }
    points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = colsca[vessel[,"T_NUM"]], pch = 20, cex = 0.6)
    points(cbind(vessel[c(1, nrow(vessel)),"LON"], vessel[c(1, nrow(vessel)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
  }else{
    map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
    while(class(map2plot) == "try-error")
    {
      span <- span + 0.25
      xrange <- extendrange(x = vessel["LON"], f = span)
      yrange <- extendrange(x = vessel["LAT"], f = span)
      map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(0, 0, 0, 0), plot = F), silent = T)
    }
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel ", vessel[1,1]," with ", length(unique(vessel[,"T_NUM"]))," tracks" , sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    colsca <- rainbow(length(unique(vessel[, "T_NUM"])), alpha = 0.5)
    textn <- 0
    for (i in 1:(nrow(vessel)-1))
    {
      trackn <- vessel[i,"T_NUM"]
      if(trackn != 0)
      {
        if(trackn == vessel[i+1,"T_NUM"])
        {
          lines(cbind(c(vessel[i,"LON"],vessel[i+1,"LON"]), c(vessel[i,"LAT"], vessel[i+1,"LAT"])), col = colsca[vessel[i,"T_NUM"]], lty = "dashed", cex = 1.5)
          if((vessel[i,"W_HARB"] != 1) & (trackn != textn))
          {
            text(vessel[i,"LON"], vessel[i,"LAT"], col = colsca[vessel[i,"T_NUM"]], labels = vessel[i,"T_NUM"])
            textn <- trackn
          }
        }
      }
    }            
    points(cbind(vessel[,"LON"], vessel[,"LAT"]), col = colsca[vessel[,"T_NUM"]], pch = 20, cex = 0.6)
    points(cbind(vessel[c(1, nrow(vessel)),"LON"], vessel[c(1, nrow(vessel)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
    span <- 0.25
  }
}
########################################
########################################



########################################
#           Track Viewer 2
########################################
trackview2 <- function(track, trackn, bathy = "", xrange, yrange)
{
  span <- 0.25
  if(diff(xrange) == 0){
    xrange <- c(track[1,"LON"]-0.01, track[1,"LON"]+0.01)
  }
  if(diff(yrange) == 0){
    yrange <- c(track[1,"LAT"]-0.01, track[1,"LAT"]+0.01)
  }
  map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
  if(class(map2plot) != "try-error")
  {
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel: ", track[1,1],"  -  Track: ",trackn, sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    textn <- 0
    for (i in 1:(nrow(track)-1))
    {
      lines(cbind(c(track[i,"LON"], track[i+1,"LON"]), c(track[i,"LAT"], track[i+1,"LAT"])), lty = "dashed", cex = 1.5)
      if((track[i,"W_HARB"] != 1) & (trackn != textn))
      {
        text(track[i,"LON"], track[i,"LAT"], labels = track[i,"T_NUM"])
        textn <- trackn
      }
    }
    points(cbind(track[,"LON"], track[,"LAT"]), pch = 20, cex = 0.6)
    points(cbind(track[c(1, nrow(track)),"LON"], track[c(1, nrow(track)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
  }else{
    while(class(map2plot) == "try-error")
    {
      span <- span + 0.25
      xrange <- extendrange(x = track["LON"], f = span)
      yrange <- extendrange(x = track["LAT"], f = span)
      map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
    }
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel: ", track[1,1],"  -  Track: ",trackn, sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    textn <- 0
    for (i in 1:(nrow(track)-1))
    {
      lines(cbind(c(track[i,"LON"], track[i+1,"LON"]), c(track[i,"LAT"], track[i+1,"LAT"])), lty = "dashed", cex = 1.5)
      if((track[i,"W_HARB"] != 1) & (trackn != textn))
      {
        text(track[i,"LON"], track[i,"LAT"], labels = track[i,"T_NUM"])
        textn <- trackn
      }
    }            
    points(cbind(track[,"LON"], track[,"LAT"]), pch = 20, cex = 0.6)
    points(cbind(track[c(1, nrow(track)),"LON"], track[c(1, nrow(track)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
    span <- 0.25
  }
}
########################################
########################################



########################################
#       Interpolation Viewer 2
########################################
intrpview2 <- function(track, trackn, bathy = "", fishi = FALSE, xrange, yrange)
{
  span <- 0.25
  if(diff(xrange) == 0){
    xrange <- c(track[1,"LON"]-0.01, track[1,"LON"]+0.01)
  }
  if(diff(yrange) == 0){
    yrange <- c(track[1,"LAT"]-0.01, track[1,"LAT"]+0.01)
  }
  map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
  if(class(map2plot) != "try-error")
  {
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel: ", track[1,1],"  -  Track: ",trackn, sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    textn <- 0
    for (i in 1:(nrow(track)-1))
    {
      lines(cbind(c(track[i,"LON"], track[i+1,"LON"]), c(track[i,"LAT"], track[i+1,"LAT"])), lty = "dashed", cex = 1.5)
      if((track[i,"W_HARB"] != 1) & (trackn != textn))
      {
        text(track[i,"LON"], track[i,"LAT"], labels = track[i,"T_NUM"])
        textn <- trackn
      }
      if(track[i,"P_INT"] == 0)
      {
        text(track[i,"LON"], track[i,"LAT"], labels = i, pos = 2, cex = 0.8)
      }
    }
    points(cbind(track[,"LON"], track[,"LAT"]), pch = 20, cex = 0.6)
    if(fishi == TRUE) {points(cbind(track[which(track["FISH"] == 1),"LON"], track[which(track["FISH"] == 1),"LAT"]), pch = 19, cex = 0.5, col = "goldenrod1")}
    points(cbind(track[c(1, nrow(track)),"LON"], track[c(1, nrow(track)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
  }else{
    while(class(map2plot) == "try-error")
    {
      span <- span + 0.25
      xrange <- extendrange(x = track["LON"], f = span)
      yrange <- extendrange(x = track["LAT"], f = span)
      map2plot <- try(maps::map("worldHires", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), plot = F), silent = T)
    }
    maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0))
    if(bathy$path != ""){
      bathyview(bathy, xrange, yrange)
      maps::map("worldHires", fill=T,col="honeydew3",bg="lightsteelblue1", ylim = yrange, xlim = xrange, mar = c(6, 6, 0, 0), add = TRUE)
    }
    maps::map.axes()
    title(main = paste("Vessel: ", track[1,1],"  -  Track: ",trackn, sep = ""), line = 0.3)
    title(xlab = "Lon", ylab = "Lat", line = 2)
    textn <- 0
    for (i in 1:(nrow(track)-1))
    {
      lines(cbind(c(track[i,"LON"], track[i+1,"LON"]), c(track[i,"LAT"], track[i+1,"LAT"])), lty = "dashed", cex = 2)
      if((track[i,"W_HARB"] != 1) & (trackn != textn))
      {
        text(track[i,"LON"], track[i,"LAT"], labels = track[i,"T_NUM"])
        textn <- trackn
      }
      if(track[i,"P_INT"] == 0)
      {
        text(track[i,"LON"], track[i,"LAT"], labels = i, pos = 2, cex = 0.8)
      }
    }
  }
  points(cbind(track[,"LON"], track[,"LAT"]), pch = 20, cex = 0.6)
  if(fishi == TRUE) {points(cbind(track[which(track["FISH"] == 1),"LON"], track[which(track["FISH"] == 1),"LAT"]), pch = 19, cex = 0.5, col = "goldenrod1")}
  points(cbind(track[c(1, nrow(track)),"LON"], track[c(1, nrow(track)),"LAT"]), col = "yellow", bg = "blue", pch = c(24, 25), cex = 0.8)
  span <- 0.25
}
########################################
########################################





add.alpha <- function(col, alpha=1){
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}