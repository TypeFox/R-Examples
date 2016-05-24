nicheplot <- function(h,...) {
  if ("cov" %in% ls(envir=.GlobalEnv)==FALSE) {gmessage("Sorry! Attribute file is not loaded.", parent = window)} else {
  assign("d1",character(0),envir=.GlobalEnv)
  assign("d2",character(0),envir=.GlobalEnv)
  assign("d3","",envir=.GlobalEnv)
  assign("d4",character(0),envir=.GlobalEnv)
  assign("d5",character(0),envir=.GlobalEnv)
  assign("d6",1.5,envir=.GlobalEnv)
  assign("d7",character(0),envir=.GlobalEnv)
  assign("d8","FALSE",envir=.GlobalEnv)
  assign("d9",character(0),envir=.GlobalEnv)
  assign("m1",names(cov),envir=.GlobalEnv)
  assign("m3",names(cov),envir=.GlobalEnv)
  toplevel <- gwindow("Niche Plot", width=800, height=800, parent = window)
  tbl0 <- gtable(m1,expand=TRUE,multiple=TRUE,cont=toplevel)
  cg <- gpanedgroup(width=100, cont = toplevel)
  cg1 <- ggroup(horizontal = FALSE, cont = cg)
  addSpace(cg1,20,horizontal = FALSE)
  gbplus1 <- gbutton("+", cont=cg1)
  gbminus1 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus2 <- gbutton("+", cont=cg1)
  gbminus2 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus4 <- gbutton("+", cont=cg1)
  gbminus4 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus5 <- gbutton("+", cont=cg1)
  gbminus5 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus7 <- gbutton("+", cont=cg1)
  gbminus7 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  cg2 <- gframe("Options", horizontal=FALSE, cont=cg)
  d1temp <- data.frame(Node.ids="",stringsAsFactors=FALSE)
  d2temp <- data.frame(Ecology.ids="",stringsAsFactors=FALSE)
  d4temp <- data.frame(Dimensions=rep("",length(m3)),stringsAsFactors=FALSE)
  d5temp <- data.frame(Groups=rep("",length(m3)),stringsAsFactors=FALSE)
  d7temp <- data.frame(Weights="",stringsAsFactors=FALSE)
  tbl1 <- gtable(d1temp,expand=TRUE,multiple=FALSE,cont=cg2)
  tbl2 <- gtable(d2temp,expand=TRUE,multiple=FALSE,cont=cg2)
  if (('el' %in% ls(envir=.GlobalEnv))==FALSE) {
    gcheckboxgroup("Show nodes",cont=cg2,handler = function(h,...) assign("d3",svalue(h$obj),envir=.GlobalEnv))
  }
  if ('el' %in% ls(envir=.GlobalEnv)) {
    gradio(c("Do not show nodes and network ties","Show nodes","Show nodes and network ties"),selected=1,cont=cg2,handler = function(h,...) assign("d3",svalue(h$obj),envir=.GlobalEnv))
  }
  tbl4 <- gtable(d4temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl5 <- gtable(d5temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl7 <- gtable(d7temp,expand=TRUE,multiple=FALSE,cont=cg2)
  glabel("Dev.range",cont=cg2)
  gslider(from = 0, to = 5, by = .05, value = 1.5, cont=cg2, handler = function(h,...) assign("d6",svalue(h$obj),envir=.GlobalEnv))
  glabel("Complete.cases",cont=cg2)
  gradio(c("TRUE","FALSE"), selected = 2, cont = cg2,  handler = function(h,...) assign("d8",svalue(h$obj),envir=.GlobalEnv))
  addHandlerClicked(gbplus1, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl1[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d1",m3[which(m3 %in% union(d1,temp))],envir=.GlobalEnv) 
      tbl1[1] <- d1
    }
  })
  addHandlerClicked(gbminus1, handler = function(h,...) {
    temp <- svalue(tbl1)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d1",character(0),envir=.GlobalEnv) 
      tbl1[1] <- ""
    }
  })
  addHandlerClicked(gbplus2, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl2[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d2",m3[which(m3 %in% union(d2,temp))],envir=.GlobalEnv) 
      tbl2[1] <- d2
    }
  })
  addHandlerClicked(gbminus2, handler = function(h,...) {
    temp <- svalue(tbl2)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d2",character(0),envir=.GlobalEnv) 
      tbl2[1] <- ""
    }
  })
  addHandlerClicked(gbplus4, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d4",m3[which(m3 %in% union(d4,temp))],envir=.GlobalEnv) 
      kd4 <- c(d4,rep("",length(m3)-length(d4)))
      for (j in 1:length(m3)) tbl4[j] <- kd4[j]
    }
  })
  addHandlerClicked(gbminus4, handler = function(h,...) {
    temp <- svalue(tbl4)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d4",m3[which(m3 %in% setdiff(d4,temp))],envir=.GlobalEnv) 
      kd4 <- c(d4,rep("",length(m3)-length(d4)+1))
      for (j in 1:length(m3)) tbl4[j] <- kd4[j]
    }
  })
  addHandlerClicked(gbplus5, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d5",m3[which(m3 %in% union(d5,temp))],envir=.GlobalEnv) 
      kd5 <- c(d5,rep("",length(m3)-length(d5)))
      for (j in 1:length(m3)) tbl5[j] <- kd5[j]
    }
  })
  addHandlerClicked(gbminus5, handler = function(h,...) {
    temp <- svalue(tbl5)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d5",m3[which(m3 %in% setdiff(d5,temp))],envir=.GlobalEnv) 
      kd5 <- c(d5,rep("",length(m3)-length(d5)+1))
      for (j in 1:length(m3)) tbl5[j] <- kd5[j]
    }
  })
  addHandlerClicked(gbplus7, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl7[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d7",m3[which(m3 %in% union(d7,temp))],envir=.GlobalEnv) 
      tbl7[1] <- d7
    }
  })
  addHandlerClicked(gbminus7, handler = function(h,...) {
    temp <- svalue(tbl7)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("d7",character(0),envir=.GlobalEnv) 
      tbl7[1] <- ""
    }
  })
  gbutton("Continue", expand = FALSE, cont = cg2, handler = function(h, ...) {
    if (length(d1)==0 | length(d4)==0 | length(d5)==0) {gmessage("Missing required information.", parent = toplevel)} else {
    dispose(toplevel)
    if (length(d7)>0) tmpweight <- d7 else tmpweight <- NULL  
    bcolors <- rep(c(26,116,142,47,8,12,31,32,33,41,51,53,62,139,151,175,153,85,450,477),5)
    if (length(d2)>0) {
      extralevel <- gwindow("Niche Plot Options", width=600, height=600, parent = window)
      extrag <- ggroup(cont = extralevel, use.scrollwindow=T, horizontal = FALSE)
      tbl <- glayout(cont = extrag)
      j <- 1
      tbl[j,1] <- "Please select which ecology you want to make niche plot:"
      tbl[j,2] <- gcheckboxgroup(c("all",as.matrix(unique(cov[which(colnames(cov)==d2)]))), handler = function(h,...){
        assign("d9",svalue(h$obj),envir=.GlobalEnv)
      })
      button <- gbutton("Continue", expand = FALSE, cont = extrag, handler = function(h, ...) {
        dispose(extralevel)
        if ("all" %in% d9) {
          cov1 <- cov
          if ('el' %in% ls(envir=.GlobalEnv)) el1 <- el
        } else {
          cov1 <- cov[which(cov[,which(colnames(cov)==d2)] %in% d9),] 
          if ('el' %in% ls(envir=.GlobalEnv)) {
            el1 <- el[which(el[,1] %in% unique(cov1[,which(colnames(cov1)==d1)]) & el[,2] %in% unique(cov1[,which(colnames(cov1)==d1)])),]
          }
        }
        if (('el' %in% ls(envir=.GlobalEnv))==FALSE) {
          b <- blau(cov1, node.ids=d1, ecology.ids=d2, dimension=d4, memberships=d5,weights=tmpweight,complete.cases=d8)
          b <- niches(b, dev.range = d6)
        }
        if ('el' %in% ls(envir=.GlobalEnv)) {
          b <- blau(cov1, node.ids=d1, ecology.ids=d2, graph = el1, dimension=d4, memberships=d5,weights=tmpweight,complete.cases=d8 )
          b <- niches(b, dev.range = d6)
        }
        assign("bobj",b,envir=.GlobalEnv)
        if (length(d4)==2) {
          lowb <- b$lowbounds
          topb <- b$topbounds
          k1 <- unique(b$ids$ecologyId)
          k2 <- c(which(names(cov1)==d4[1]),which(names(cov1)==d4[2]))
          niche2d1 <- function() {
            add_legend <- function(...) {
              opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
              on.exit(par(opar))
              plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
              legend(...)
            }
            plot(-10000,-10000,xlim=c(min(lowb[,1],cov1[,k2[1]],na.rm = TRUE),max(topb[,1],cov1[,k2[1]],na.rm = TRUE)),ylim=c(min(lowb[,2],cov1[,k2[2]],na.rm = TRUE),max(topb[,2],cov1[,k2[2]],na.rm = TRUE)), xlab=d4[1], ylab=d4[2])
            if (length(grep('Show nodes',d3))==1) {
              plot(cov1[,k2[1]],cov1[,k2[2]],xlim=c(min(lowb[,1],cov1[,k2[1]],na.rm = TRUE),max(topb[,1],cov1[,k2[1]],na.rm = TRUE)),ylim=c(min(lowb[,2],cov1[,k2[2]],na.rm = TRUE),max(topb[,2],cov1[,k2[2]],na.rm = TRUE)), xlab=d4[1], ylab=d4[2], col="red")
            }
            bcol <- 1
            for (i in 1:length(k1)) {
              for (j in 1:length(d5)) {
                k <- (i-1)*length(d5)+j
                rect(lowb[k,1], lowb[k,2], topb[k,1], topb[k,2], border=colors()[bcolors[bcol]])
                bcol <- bcol+1
              }
            }        
            if (d3=="Show nodes and network ties") {
              x1 <- y1 <- x2 <- y2 <- rep(0,nrow(el))
              for (i in 1:nrow(el1)) {
                s <- which(as.character(el1[i,1])==as.character(cov[,d1]))
                x1[i] <- cov[s,k2[1]]
                y1[i] <- cov[s,k2[2]]
                t <- which(as.character(el1[i,2])==as.character(cov[,d1]))
                x2[i] <- cov[t,k2[1]]
                y2[i] <- cov[t,k2[2]]
              }
              segments(x1,y1,x2,y2,col="red")
            }
            add_legend("topleft",paste(rep(k1,each=length(d5)),"_",rownames(lowb),sep=""),text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
          }
          sndlevel <- gwindow("2D Niche Plot", width = 600, height = 600)
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
              png("niche2d.png", width = o1, height = o2, res = o3)
              niche2d1()
              dev.off()}
            if (o4=="jpg") {
              jpeg("niche2d.jpg", width = o1, height = o2, res = o3)
              niche2d1()
              dev.off()}
            if (o4=="bmp") {
              bmp("niche2d.bmp", width = o1, height = o2, res = o3)
              niche2d1()
              dev.off()}
            if (o4=="tiff") {
              tiff("niche2d.tiff", width = o1, height = o2, res = o3)
              niche2d1()
              dev.off()}
            if (o4=="pdf") {
              pdf("niche2d.pdf", width = o1/100, height = o2/100)
              niche2d1()
              dev.off()}
          })
          plot1 <- ggraphics(container=sg)
          dev1 <- dev.cur();
          Sys.sleep(0.5)
          niche2d1()
        } else if (length(d4) == 3) {
          lowb <- b$lowbounds
          topb <- b$topbounds
          k1 <- unique(b$ids$ecologyId)
          k2 <- c(which(names(cov)==d4[1]),which(names(cov)==d4[2]),which(names(cov)==d4[3]))
          niche3d <- function(phi1,theta1) {
            add_legend <- function(...) {
              opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
              on.exit(par(opar))
              plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
              legend(...)
            }
            box3D(lowb[,1], lowb[,2], lowb[,3],
                  topb[,1], topb[,2], topb[,3], phi=phi1, theta = theta1,
                  xlim=c(min(lowb[,1],cov1[,k2[1]],na.rm = TRUE),max(topb[,1],cov1[,k2[1]],na.rm = TRUE)),
                  ylim=c(min(lowb[,2],cov1[,k2[2]],na.rm = TRUE),max(topb[,2],cov1[,k2[2]],na.rm = TRUE)),
                  zlim=c(min(lowb[,3],cov1[,k2[3]],na.rm = TRUE),max(topb[,3],cov1[,k2[3]],na.rm = TRUE)),
                  xlab=d4[1],ylab=d4[2],zlab=d4[3],bty = "f", ticktype = "detailed",
                  col = colors()[bcolors[1:nrow(lowb)]], alpha = 0.3,
                  border = "black", lwd = 2)
            if (length(grep('Show nodes',d3))==1) {
              scatter3D(cov1[,k2[1]],cov1[,k2[2]],cov1[,k2[3]], pch = 16, 
                    add = TRUE,
                    colkey = FALSE, col="red")
            }
            if (d3=="Show nodes and network ties") {
              x1 <- y1 <- z1 <- x2 <- y2 <- z2 <- rep(0,nrow(el1))
              for (i in 1:nrow(el1)) {
                s <- which(as.character(el1[i,1])==as.character(cov1[,d1]))
                x1[i] <- cov1[s,k2[1]]
                y1[i] <- cov1[s,k2[2]]
                z1[i] <- cov1[s,k2[3]]
                t <- which(as.character(el1[i,2])==as.character(cov1[,d1]))
                x2[i] <- cov1[t,k2[1]]
                y2[i] <- cov1[t,k2[2]]
                z2[i] <- cov1[t,k2[3]]
              }
              segments3D(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
            }
            add_legend("topleft",paste(rep(k1,each=length(d5)),"_",rownames(lowb),sep=""),text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
          } 
          niche3drgl1 <- function() {
            scatter3Drgl(-10000,-10000,-10000, 
                   xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),
                   ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)),
                   zlim=c(min(lowb[,3],cov[,k2[3]],na.rm = TRUE),max(topb[,3],cov[,k2[3]],na.rm = TRUE)),
                   xlab=d4[1],ylab=d4[2],zlab=d4[3],
                   colkey = FALSE, col="red")
            par3d(windowRect = c(900, 50, 1500, 650))
            if (length(grep('Show nodes',d3))==1) {
              rgl.close()
              scatter3Drgl(cov[,k2[1]],cov[,k2[2]],cov[,k2[3]], 
                    xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),
                    ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)),
                    zlim=c(min(lowb[,3],cov[,k2[3]],na.rm = TRUE),max(topb[,3],cov[,k2[3]],na.rm = TRUE)),
                    xlab=d4[1],ylab=d4[2],zlab=d4[3],
                    colkey = FALSE, col="red")
              par3d(windowRect = c(900, 50, 1500, 650))
            }
            axes3d()
            box3Drgl(lowb[,1], lowb[,2], lowb[,3],
                  topb[,1], topb[,2], topb[,3], 
                  add = TRUE,
                  col = colors()[bcolors[1:nrow(lowb)]], alpha = 0.5,
                  border = "black", lwd = 2)
            if (d3=="Show nodes and network ties") {
              x1 <- y1 <- z1 <- x2 <- y2 <- z2 <- rep(0,nrow(el1))
              for (i in 1:nrow(el1)) {
                s <- which(as.character(el1[i,1])==as.character(cov1[,d1]))
                x1[i] <- cov1[s,k2[1]]
                y1[i] <- cov1[s,k2[2]]
                z1[i] <- cov1[s,k2[3]]
                t <- which(as.character(el1[i,2])==as.character(cov1[,d1]))
                x2[i] <- cov1[t,k2[1]]
                y2[i] <- cov1[t,k2[2]]
                z2[i] <- cov1[t,k2[3]]
              }
              segments3Drgl(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
            }
            legend3d("topleft",paste(rep(k1,each=length(d5)),"_",rownames(lowb),sep=""),text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
          } 
          niche3drgl1()
          fthlevel <- gwindow("3D Niche Plot", width = 800, height = 800)
          fg <- ggroup(cont = fthlevel, horizontal = F, expand=T)
          ang1 <- ang2 <- "0"
          assign("n3d1",600,envir=.GlobalEnv)
          assign("n3d2",600,envir=.GlobalEnv)
          assign("n3d3",100,envir=.GlobalEnv)
          assign("n3d4","png",envir=.GlobalEnv)
          tbl <- glayout(cont = fg)
          i <- 1
          tbl[i,i] <- "horizonal angle:"
          tbl[i,2] <- gedit("0") 
          tbl[i,3] <- "vertical angle"
          tbl[i,4] <- gedit("0") 
          tbl[i,5] <- "width"
          tbl[i,6] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
            handler = function(h,...){assign("n3d1",svalue(h$obj),envir=.GlobalEnv)})
          tbl[i,7] <- "height"
          tbl[i,8] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
            handler = function(h,...){assign("n3d2",svalue(h$obj),envir=.GlobalEnv)})
          tbl[i,9] <- "DPI"
          tbl[i,10] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
            handler = function(h,...){assign("n3d3",svalue(h$obj),envir=.GlobalEnv)})
          tbl[i,11] <- "format"
          tbl[i,12] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
            handler = function(h,...){assign("n3d4",svalue(h$obj),envir=.GlobalEnv)})
          button <- gbutton("Refresh Plot", cont = fg, handler = function(h, ...) {
             niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
          })
          button <- gbutton("Save Plot", cont = fg, handler = function(h, ...) {
            if (n3d4=="png") {
              png("niche3d.png", width = n3d1, height = n3d2, res = n3d3)
              niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
              dev.off()}
            if (n3d4=="jpg") {
              jpeg("niche3d.jpg", width = n3d1, height = n3d2, res = n3d3)
              niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
              dev.off()}
            if (n3d4=="bmp") {
              bmp("niche3d.bmp", width = n3d1, height = n3d2, res = n3d3)
              niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
              dev.off()}
            if (n3d4=="tiff") {
              tiff("niche3d.tiff", width = n3d1, height = n3d2, res = n3d3)
              niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
              dev.off()}
            if (n3d4=="pdf") {
              pdf("niche3d.pdf", width = n3d1/100, height = n3d2/100)
              niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
              dev.off()}
          })
          plot1 <- ggraphics(container=fg)
          dev1 <- dev.cur();
          Sys.sleep(0.5)
          niche3d(as.numeric(ang1),as.numeric(ang2))
        } else {
          dispose(toplevel)
          gmessage("You must select 2  or 3 dimensions.")
        }
      })  
    } else {
      if (('el' %in% ls(envir=.GlobalEnv))==FALSE) {
        b <- blau(cov, node.ids=d1, dimension=d4, memberships=d5,weights=tmpweight,complete.cases=d8)
        b <- niches(b, dev.range = d6)
      }
      if ('el' %in% ls(envir=.GlobalEnv)) {
        b <- blau(cov, node.ids=d1, graph = el, dimension=d4, memberships=d5,weights=tmpweight,complete.cases=d8)
        b <- niches(b, dev.range = d6)
      }
      assign("bobj",b,envir=.GlobalEnv)
      if (length(d4)==2) {
        lowb <- b$lowbounds
        topb <- b$topbounds
        k2 <- c(which(names(cov)==d4[1]),which(names(cov)==d4[2]))
        niche2d2 <- function() {
          add_legend <- function(...) {
            opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
            mar=c(0, 0, 0, 0), new=TRUE)
            on.exit(par(opar))
            plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
            legend(...)
          }
          plot(-10000,-10000,xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)), xlab=d4[1], ylab=d4[2])
          if (length(grep('Show nodes',d3))==1) {
            plot(cov[,k2[1]],cov[,k2[2]],xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)), xlab=d4[1], ylab=d4[2],col="red")
          }
          for (i in 1:length(d5)) {
            rect(lowb[i,1], lowb[i,2], topb[i,1], topb[i,2], border=colors()[bcolors[i]])
          }
          if (d3=="Show nodes and network ties") {
            r <- which(d1==names(cov))
            x1 <- y1 <- x2 <- y2 <- rep(0,nrow(el))
            for (i in 1:nrow(el)) {
              s <- which(as.character(el[i,1])==as.character(cov[,r]))
              x1[i] <- cov[s,k2[1]]
              y1[i] <- cov[s,k2[2]]
              t <- which(as.character(el[i,2])==as.character(cov[,r]))
              x2[i] <- cov[t,k2[1]]
              y2[i] <- cov[t,k2[2]]
            }
            segments(x1,y1,x2,y2,col="red")
          }  
          add_legend("topleft",d5,text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
        }
        sndlevel <- gwindow("2D Niche Plot", width = 600, height = 600)
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
            png("niche2d.png", width = o1, height = o2, res = o3)
            niche2d2()
            dev.off()}
          if (o4=="jpg") {
            jpeg("niche2d.jpg", width = o1, height = o2, res = o3)
            niche2d2() 
            dev.off()}
          if (o4=="bmp") {
            bmp("niche2d.bmp", width = o1, height = o2, res = o3)
            niche2d2()
            dev.off()}
          if (o4=="tiff") {
            tiff("niche2d.tiff", width = o1, height = o2, res = o3)
            niche2d2() 
            dev.off()}
          if (o4=="pdf") {
            pdf("niche2d.pdf", width = o1/100, height = o2/100)
            niche2d2() 
            dev.off()}
        })
        plot1 <- ggraphics(container=sg)
        dev1 <- dev.cur();
        Sys.sleep(0.5)
        niche2d2()       
      } else if (length(d4) == 3) {
        lowb <- b$lowbounds
        topb <- b$topbounds
        k2 <- c(which(names(cov)==d4[1]),which(names(cov)==d4[2]),which(names(cov)==d4[3]))
        niche3d <- function(phi1,theta1) {
          add_legend <- function(...) {
            opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
            mar=c(0, 0, 0, 0), new=TRUE)
            on.exit(par(opar))
            plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
            legend(...)
          }
          box3D(lowb[,1], lowb[,2], lowb[,3],
                topb[,1], topb[,2], topb[,3], phi=phi1, theta = theta1,
                xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),
                ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)),
                zlim=c(min(lowb[,3],cov[,k2[3]],na.rm = TRUE),max(topb[,3],cov[,k2[3]],na.rm = TRUE)),
                xlab=d4[1],ylab=d4[2],zlab=d4[3],bty = "f", ticktype = "detailed",
                col = colors()[bcolors[1:nrow(lowb)]], alpha = 0.3,
                border = "black", lwd = 2)
          if (length(grep('Show nodes',d3))==1) {
            scatter3D(cov[,k2[1]],cov[,k2[2]],cov[,k2[3]], pch = 16, 
                  add = TRUE,
                  colkey = FALSE, col="red")
          }
          if (d3=="Show nodes and network ties") {
            r <- which(d1==names(cov))
            x1 <- y1 <- z1 <- x2 <- y2 <- z2 <- rep(0,nrow(el))
            for (i in 1:nrow(el)) {
              s <- which(as.character(el[i,1])==as.character(cov[,r]))
              x1[i] <- cov[s,k2[1]]
              y1[i] <- cov[s,k2[2]]
              z1[i] <- cov[s,k2[3]]
              t <- which(as.character(el[i,2])==as.character(cov[,r]))
              x2[i] <- cov[t,k2[1]]
              y2[i] <- cov[t,k2[2]]
              z2[i] <- cov[t,k2[3]]
            }
            segments3D(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
          }
          add_legend("topleft",d5,text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
        }
        niche3drgl2 <- function() {
          scatter3Drgl(-10000,-10000,-10000, 
                 xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),
                 ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)),
                 zlim=c(min(lowb[,3],cov[,k2[3]],na.rm = TRUE),max(topb[,3],cov[,k2[3]],na.rm = TRUE)),
                 xlab=d4[1],ylab=d4[2],zlab=d4[3],
                 colkey = FALSE, col="red")
          par3d(windowRect = c(900, 50, 1500, 650))
          if (length(grep('Show nodes',d3))==1) {
            rgl.close()
            scatter3Drgl(cov[,k2[1]],cov[,k2[2]],cov[,k2[3]], 
                  xlim=c(min(lowb[,1],cov[,k2[1]],na.rm = TRUE),max(topb[,1],cov[,k2[1]],na.rm = TRUE)),
                  ylim=c(min(lowb[,2],cov[,k2[2]],na.rm = TRUE),max(topb[,2],cov[,k2[2]],na.rm = TRUE)),
                  zlim=c(min(lowb[,3],cov[,k2[3]],na.rm = TRUE),max(topb[,3],cov[,k2[3]],na.rm = TRUE)),
                  xlab=d4[1],ylab=d4[2],zlab=d4[3],
                  colkey = FALSE, col="red")
            par3d(windowRect = c(900, 50, 1500, 650))
          }
          axes3d()
          box3Drgl(lowb[,1], lowb[,2], lowb[,3],
                topb[,1], topb[,2], topb[,3], 
                add = TRUE,
                col = colors()[bcolors[1:nrow(lowb)]], alpha = 0.5,
                border = "black", lwd = 2)
          if (d3=="Show nodes and network ties") {
            r <- which(d1==names(cov))
            x1 <- y1 <- z1 <- x2 <- y2 <- z2 <- rep(0,nrow(el))
            for (i in 1:nrow(el)) {
              s <- which(as.character(el[i,1])==as.character(cov[,r]))
              x1[i] <- cov[s,k2[1]]
              y1[i] <- cov[s,k2[2]]
              z1[i] <- cov[s,k2[3]]
              t <- which(as.character(el[i,2])==as.character(cov[,r]))
              x2[i] <- cov[t,k2[1]]
              y2[i] <- cov[t,k2[2]]
              z2[i] <- cov[t,k2[3]]
            }
            segments3Drgl(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
          }
          legend3d("topleft",d5,text.col=colors()[bcolors[1:nrow(lowb)]],bty='n')
        }
        niche3drgl2()
        fthlevel <- gwindow("3D Niche Plot", width = 800, height = 800)
        fg <- ggroup(cont = fthlevel, horizontal = F, expand=T)
        ang1 <- ang2 <- "0"
        assign("n3d1",600,envir=.GlobalEnv)
        assign("n3d2",600,envir=.GlobalEnv)
        assign("n3d3",100,envir=.GlobalEnv)
        assign("n3d4","png",envir=.GlobalEnv)
        tbl <- glayout(cont = fg)
        i <- 1
        tbl[i,i] <- "horizonal angle:"
        tbl[i,2] <- gedit("0") 
        tbl[i,3] <- "vertical angle"
        tbl[i,4] <- gedit("0") 
        tbl[i,5] <- "width"
        tbl[i,6] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
          handler = function(h,...){assign("n3d1",svalue(h$obj),envir=.GlobalEnv)})
        tbl[i,7] <- "height"
        tbl[i,8] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
          handler = function(h,...){assign("n3d2",svalue(h$obj),envir=.GlobalEnv)})
        tbl[i,9] <- "DPI"
        tbl[i,10] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
          handler = function(h,...){assign("n3d3",svalue(h$obj),envir=.GlobalEnv)})
        tbl[i,11] <- "format"
        tbl[i,12] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
          handler = function(h,...){assign("n3d4",svalue(h$obj),envir=.GlobalEnv)})
        button <- gbutton("Refresh Plot", cont = fg, handler = function(h, ...) {
           niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
        })
        button <- gbutton("Save Plot", cont = fg, handler = function(h, ...) {
          if (n3d4=="png") {
            png("niche3d.png", width = n3d1, height = n3d2, res = n3d3)
            niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            dev.off()}
          if (n3d4=="jpg") {
            jpeg("niche3d.jpg", width = n3d1, height = n3d2, res = n3d3)
            niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            dev.off()}
          if (n3d4=="bmp") {
            bmp("niche3d.bmp", width = n3d1, height = n3d2, res = n3d3)
            niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            dev.off()}
          if (n3d4=="tiff") {
            tiff("niche3d.tiff", width = n3d1, height = n3d2, res = n3d3)
            niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            dev.off()}
          if (n3d4=="pdf") {
            pdf("niche3d.pdf", width = n3d1/100, height = n3d2/100)
            niche3d(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            dev.off()}
        })
        plot1 <- ggraphics(container=fg)
        dev1 <- dev.cur();
        Sys.sleep(0.5)
        niche3d(as.numeric(ang1),as.numeric(ang2))
      } else {
        gmessage("You must select 2  or 3 dimensions.")
      }
    }
  }})}
}
