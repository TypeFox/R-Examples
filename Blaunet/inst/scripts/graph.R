showgraph <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    n <- network(adj)
    assign("g1",0,envir=.GlobalEnv)
    assign("g2",0,envir=.GlobalEnv)
    assign("g3",0,envir=.GlobalEnv)
    assign("g4",0,envir=.GlobalEnv)
    assign("g5","fruchtermanreingold",envir=.GlobalEnv)
    toplevel <- gwindow("Plot Graph", width = 400, height = 220, parent = window)
    cg <- ggroup(cont = toplevel, use.scrollwindow=T, horizontal = FALSE)
    tbl <- glayout(cont = cg)
    i <- 1
    tbl[i,1] <- gcheckboxgroup("Show vertex name",handler = function(h,...){
       assign("g1",svalue(h$obj),envir=.GlobalEnv)
    })
    i <- i + 1
    tbl[i,1] <- "Vertex color by"
    tbl[i,2] <- gcombobox(c("None",names(cov)), selected = 1, cont = tbl, 
      handler = function(h,...){
        if (svalue(h$obj) %in% names(cov)) {
          assign("g2",which(names(cov) %in% svalue(h$obj)),envir=.GlobalEnv)
        } else {assign("g2",0,envir=.GlobalEnv)}
    })
    i <- i + 1
    tbl[i,1] <- "Vertex side by"
    tbl[i,2] <- gcombobox(c("None",names(cov)), selected = 1, cont = tbl, 
      handler = function(h,...){
        if (svalue(h$obj) %in% names(cov)) {
          assign("g3",which(names(cov) %in% svalue(h$obj)),envir=.GlobalEnv)
        } else {assign("g3",0,envir=.GlobalEnv)}
    })
    i <- i + 1
    tbl[i,1] <- "Vertex size by"
    tbl[i,2] <- gcombobox(c("None",names(cov)), selected = 1, cont = tbl, 
      handler = function(h,...){
        if (svalue(h$obj) %in% names(cov)) {
          assign("g4",which(names(cov) %in% svalue(h$obj)),envir=.GlobalEnv)
        } else {assign("g4",0,envir=.GlobalEnv)}
    })
    i <- i + 1
    tbl[i,1] <- "Layout"
    tbl[i,2] <- gcombobox(c("fruchtermanreingold","kamadakawai","spring","circle","eigen","hall","mds","princoord","target","random"), cont = tbl, 
       handler = function(h,...)assign("g5",svalue(h$obj),envir=.GlobalEnv)
    )
    button <- gbutton("Plot Graph", cont = cg, handler = function(h, ...) {
      add_legend <- function(...) {
        opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
        mar=c(0, 0, 0, 0), new=TRUE)
        on.exit(par(opar))
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
        legend(...)
      }
      net_plot <- function() {
        par(mar = c(5, 4, 1.4, 0.2))
        if (g1==0 & g2==0 & g3==0 & g4==0) gplot(n,mode=g5)
        if (g1=="Show vertex name" & g2==0 & g3==0 & g4==0) gplot(n, label=network.vertex.names(n),mode=g5)
        if (g1=="Show vertex name" & g2==0 & g3>0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.sides=temp6,vertex.cex=1.5,mode=g5)
        if (g1!="Show vertex name" & g2>0 & g3==0 & g4==0) gplot(n, vertex.col=temp3,mode=g5)
        if (g1!="Show vertex name" & g2>0 & g3>0 & g4==0) gplot(n, vertex.col=temp3,vertex.sides=temp6,vertex.cex=1.5,mode=g5)
        if (g1!="Show vertex name" & g2==0 & g3==0 & g4>0) gplot(n, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
        if (g1!="Show vertex name" & g2==0 & g3>0 & g4>0) gplot(n, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
        if (g1!="Show vertex name" & g2>0 & g3==0 & g4>0) gplot(n, vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
        if (g1!="Show vertex name" & g2>0 & g3>0 & g4>0) gplot(n, vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
        if (g1=="Show vertex name" & g2>0 & g3==0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.col=temp3,mode=g5)
        if (g1=="Show vertex name" & g2>0 & g3>0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.col=temp3,vertex.sides=temp6,vertex.cex=1.5,mode=g5)
        if (g1=="Show vertex name" & g2==0 & g3==0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
        if (g1=="Show vertex name" & g2==0 & g3>0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
        if (g1=="Show vertex name" & g2>0 & g3==0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
        if (g1=="Show vertex name" & g2>0 & g3>0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
        if (g2>0) {
          add_legend("topleft",temp2,pch=16,col=rainbow(length(temp1)),bty='n',title = names(cov)[g2])
        }
        if (g3>0) {
          add_legend("bottomleft",temp8,bty='n',title = names(cov)[g3])
        }
        par(mar=c(5, 4, 4, 2) + 0.1)
      }
      if (g2>0) {
        temp1 <- sort(unique(cov[,g2]))
        temp2 <- as.character(temp1)
        temp3 <- rep("",nrow(cov))
        for (i in 1:nrow(cov)) temp3[i] <- rainbow(length(temp1))[which(temp1 %in% cov[i,g2])]
      }
      if (g3>0) {
        temp4 <- sort(unique(cov[,g3]))
        temp5 <- as.character(temp4)
        temp6 <- rep(0,nrow(cov))
        for (i in 1:nrow(cov)) temp6[i] <- which(temp4 %in% cov[i,g3])+2
        temp7 <- c("triangle","diamond","pentagon","hexagon","heptagon","octagon","nonagon","decagon","hendecagon","dodecagon")[1:length(temp5)]
        temp8 <- rep("",length(temp5))
        for (i in 1:length(temp5)) temp8[i] <- paste(temp7[i],":",temp5[i],sep="")
      }
      sndlevel <- gwindow("Network Graph", width = 800, height = 800)
      sg <- ggroup(cont = sndlevel, horizontal = F, expand=T)
      assign("o1",600,envir=.GlobalEnv)
      assign("o2",600,envir=.GlobalEnv)
      assign("o3",100,envir=.GlobalEnv)
      assign("o4","png",envir=.GlobalEnv)
      tbl <- glayout(cont = sg)
      i <- 1
      tbl[i,i] <- "width"
      tbl[i,2] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
        handler = function(h,...){assign("o1",svalue(h$obj),envir=.GlobalEnv)})
      tbl[i,3] <- "height"
      tbl[i,4] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
        handler = function(h,...){assign("o2",svalue(h$obj),envir=.GlobalEnv)})
      tbl[i,5] <- "DPI"
      tbl[i,6] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
        handler = function(h,...){assign("o3",svalue(h$obj),envir=.GlobalEnv)})
      tbl[i,7] <- "format"
      tbl[i,8] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
        handler = function(h,...){assign("o4",svalue(h$obj),envir=.GlobalEnv)})
      button <- gbutton("Save Plot", cont = sg, handler = function(h, ...) {
        if (o4=="png") {
          png("network-graph.png", width = o1, height = o2, res = o3)
          net_plot()
          dev.off()}
        if (o4=="jpg") {
          jpeg("network-graph.jpg", width = o1, height = o2, res = o3)
          net_plot()
          dev.off()}
        if (o4=="bmp") {
          bmp("network-graph.bmp", width = o1, height = o2, res = o3)
          net_plot()
          dev.off()}
        if (o4=="tiff") {
          tiff("network-graph.tiff", width = o1, height = o2, res = o3)
          net_plot()
          dev.off()}
        if (o4=="pdf") {
          pdf("network-graph.pdf", width = o1/100, height = o2/100)
          net_plot()
          dev.off()}
      })
      plot1 <- ggraphics(container=sg)
      dev1 <- dev.cur();
      Sys.sleep(0.5)
      par(mar = c(5, 4, 1.4, 0.2))
      if (g1==0 & g2==0 & g3==0 & g4==0) gplot(n,mode=g5)
      if (g1=="Show vertex name" & g2==0 & g3==0 & g4==0) gplot(n, label=network.vertex.names(n),mode=g5)
      if (g1=="Show vertex name" & g2==0 & g3>0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.sides=temp6,vertex.cex=1.5,mode=g5)
      if (g1!="Show vertex name" & g2>0 & g3==0 & g4==0) gplot(n, vertex.col=temp3,mode=g5)
      if (g1!="Show vertex name" & g2>0 & g3>0 & g4==0) gplot(n, vertex.col=temp3,vertex.sides=temp6,vertex.cex=1.5,mode=g5)
      if (g1!="Show vertex name" & g2==0 & g3==0 & g4>0) gplot(n, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
      if (g1!="Show vertex name" & g2==0 & g3>0 & g4>0) gplot(n, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
      if (g1!="Show vertex name" & g2>0 & g3==0 & g4>0) gplot(n, vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
      if (g1!="Show vertex name" & g2>0 & g3>0 & g4>0) gplot(n, vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
      if (g1=="Show vertex name" & g2>0 & g3==0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.col=temp3,mode=g5)
      if (g1=="Show vertex name" & g2>0 & g3>0 & g4==0) gplot(n, label=network.vertex.names(n), vertex.col=temp3,vertex.sides=temp6,vertex.cex=1.5,mode=g5)
      if (g1=="Show vertex name" & g2==0 & g3==0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
      if (g1=="Show vertex name" & g2==0 & g3>0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
      if (g1=="Show vertex name" & g2>0 & g3==0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),mode=g5)
      if (g1=="Show vertex name" & g2>0 & g3>0 & g4>0) gplot(n, label=network.vertex.names(n), vertex.col=temp3, vertex.cex=cov[,g4]/mean(cov[,g4]),vertex.sides=temp6,mode=g5)
      if (g2>0) {
        add_legend("topleft",temp2,text.col=rainbow(length(temp1)),bty='n',title = names(cov)[g2])
      }
      if (g3>0) {
        add_legend("bottomleft",temp8,bty='n',title = names(cov)[g3])
      }
      par(mar=c(5, 4, 4, 2) + 0.1)
      dispose(toplevel)
    })
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}

showhoutdegree <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    outdegree <- degree(network(adj),cmode="outdegree")
    sndlevel <- gwindow("Out-degree Distribution", width = 600, height = 600)
    sg <- ggroup(cont = sndlevel, horizontal = F, expand=T)
    assign("o1",600,envir=.GlobalEnv)
    assign("o2",600,envir=.GlobalEnv)
    assign("o3",100,envir=.GlobalEnv)
    assign("o4","png",envir=.GlobalEnv)
    tbl <- glayout(cont = sg)
    i <- 1
    tbl[i,i] <- "width"
    tbl[i,2] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o1",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,3] <- "height"
    tbl[i,4] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o2",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,5] <- "DPI"
    tbl[i,6] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o3",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,7] <- "format"
    tbl[i,8] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o4",svalue(h$obj),envir=.GlobalEnv)})
    button <- gbutton("Save Plot", cont = sg, handler = function(h, ...) {
      if (o4=="png") {
        png("out-degree.png", width = o1, height = o2, res = o3)
        hist(outdegree)
        dev.off()}
      if (o4=="jpg") {
        jpeg("out-degree.jpg", width = o1, height = o2, res = o3)
        hist(outdegree)
        dev.off()}
      if (o4=="bmp") {
        bmp("out-degree.bmp", width = o1, height = o2, res = o3)
        hist(outdegree)
        dev.off()}
      if (o4=="tiff") {
        tiff("out-degree.tiff", width = o1, height = o2, res = o3)
        hist(outdegree)
        dev.off()}
      if (o4=="pdf") {
        pdf("out-degree.pdf", width = o1/100, height = o2/100)
        hist(outdegree)
        dev.off()}
    })
    plot1 <- ggraphics(container=sg)
    dev1 <- dev.cur();
    Sys.sleep(0.5)
    hist(outdegree)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showhindegree <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    indegree <- degree(network(adj),cmode="indegree")
    sndlevel <- gwindow("In-degree Distribution", width = 600, height = 600)
    sg <- ggroup(cont = sndlevel, horizontal = F, expand=T)
    assign("o1",600,envir=.GlobalEnv)
    assign("o2",600,envir=.GlobalEnv)
    assign("o3",100,envir=.GlobalEnv)
    assign("o4","png",envir=.GlobalEnv)
    tbl <- glayout(cont = sg)
    i <- 1
    tbl[i,i] <- "width"
    tbl[i,2] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o1",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,3] <- "height"
    tbl[i,4] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o2",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,5] <- "DPI"
    tbl[i,6] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o3",svalue(h$obj),envir=.GlobalEnv)})
    tbl[i,7] <- "format"
    tbl[i,8] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
      handler = function(h,...){assign("o4",svalue(h$obj),envir=.GlobalEnv)})
    button <- gbutton("Save Plot", cont = sg, handler = function(h, ...) {
      if (o4=="png") {
        png("in-degree.png", width = o1, height = o2, res = o3)
        hist(indegree)
        dev.off()}
      if (o4=="jpg") {
        jpeg("in-degree.jpg", width = o1, height = o2, res = o3)
        hist(indegree)
        dev.off()}
      if (o4=="bmp") {
        bmp("in-degree.bmp", width = o1, height = o2, res = o3)
        hist(indegree)
        dev.off()}
      if (o4=="tiff") {
        tiff("in-degree.tiff", width = o1, height = o2, res = o3)
        hist(indegree)
        dev.off()}
      if (o4=="pdf") {
        pdf("in-degree.pdf", width = o1/100, height = o2/100)
        hist(indegree)
        dev.off()}
    })
    plot1 <- ggraphics(container=sg)
    dev1 <- dev.cur();
    Sys.sleep(0.5)
    hist(indegree)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}

