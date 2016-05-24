showdynamics <- function(h,...) {
  if ("cov" %in% ls(envir=.GlobalEnv)==FALSE) {gmessage("Sorry! Attribute file is not loaded.", parent = window)} else {
  assign("dy1",character(0),envir=.GlobalEnv)
  assign("dy2",character(0),envir=.GlobalEnv)
  assign("dy3","",envir=.GlobalEnv)
  assign("dy4",character(0),envir=.GlobalEnv)
  assign("dy5",character(0),envir=.GlobalEnv)
  assign("dy6",1.5,envir=.GlobalEnv)
  assign("dy7",character(0),envir=.GlobalEnv)
  assign("dy8","FALSE",envir=.GlobalEnv)
  assign("dy9","all",envir=.GlobalEnv)
  assign("m1",names(cov),envir=.GlobalEnv)
  assign("m3",names(cov),envir=.GlobalEnv)
  toplevel <- gwindow("Niche Dynamics", width=800, height=800, parent = window)
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
  dy1temp <- data.frame(Node.ids="",stringsAsFactors=FALSE)
  dy2temp <- data.frame(Ecology.ids="",stringsAsFactors=FALSE)
  dy4temp <- data.frame(Dimensions=rep("",length(m3)),stringsAsFactors=FALSE)
  dy5temp <- data.frame(Groups=rep("",length(m3)),stringsAsFactors=FALSE)
  dy7temp <- data.frame(Weights="",stringsAsFactors=FALSE)
  tbl1 <- gtable(dy1temp,expand=TRUE,multiple=FALSE,cont=cg2)
  tbl2 <- gtable(dy2temp,expand=TRUE,multiple=FALSE,cont=cg2)
  if ('el' %in% ls(envir=.GlobalEnv)) {
    gcheckboxgroup("Network included",cont=cg2,handler = function(h,...) assign("dy3",svalue(h$obj),envir=.GlobalEnv))
  }
  tbl4 <- gtable(dy4temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl5 <- gtable(dy5temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl7 <- gtable(dy7temp,expand=TRUE,multiple=FALSE,cont=cg2)
  glabel("Dev.range",cont=cg2)
  gslider(from = 0, to = 5, by = .05, value = 1.5, cont=cg2, handler = function(h,...) assign("dy6",svalue(h$obj),envir=.GlobalEnv))
  glabel("Complete.cases",cont=cg2)
  gradio(c("TRUE","FALSE"), selected = 2, cont = cg2,  handler = function(h,...) assign("dy8",svalue(h$obj),envir=.GlobalEnv))
  addHandlerClicked(gbplus1, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl1[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy1",m3[which(m3 %in% union(dy1,temp))],envir=.GlobalEnv) 
      tbl1[1] <- dy1
    }
  })
  addHandlerClicked(gbminus1, handler = function(h,...) {
    temp <- svalue(tbl1)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy1",character(0),envir=.GlobalEnv) 
      tbl1[1] <- ""
    }
  })
  addHandlerClicked(gbplus2, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl2[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy2",m3[which(m3 %in% union(dy2,temp))],envir=.GlobalEnv) 
      tbl2[1] <- dy2
    }
  })
  addHandlerClicked(gbminus2, handler = function(h,...) {
    temp <- svalue(tbl2)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy2",character(0),envir=.GlobalEnv) 
      tbl2[1] <- ""
    }
  })
  addHandlerClicked(gbplus4, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy4",m3[which(m3 %in% union(dy4,temp))],envir=.GlobalEnv) 
      kd4 <- c(dy4,rep("",length(m3)-length(dy4)))
      for (j in 1:length(m3)) tbl4[j] <- kd4[j]
    }
  })
  addHandlerClicked(gbminus4, handler = function(h,...) {
    temp <- svalue(tbl4)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy4",m3[which(m3 %in% setdiff(dy4,temp))],envir=.GlobalEnv) 
      kd4 <- c(dy4,rep("",length(m3)-length(dy4)+1))
      for (j in 1:length(m3)) tbl4[j] <- kd4[j]
    }
  })
  addHandlerClicked(gbplus5, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy5",m3[which(m3 %in% union(dy5,temp))],envir=.GlobalEnv) 
      kd5 <- c(dy5,rep("",length(m3)-length(dy5)))
      for (j in 1:length(m3)) tbl5[j] <- kd5[j]
    }
  })
  addHandlerClicked(gbminus5, handler = function(h,...) {
    temp <- svalue(tbl5)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy5",m3[which(m3 %in% setdiff(dy5,temp))],envir=.GlobalEnv) 
      kd5 <- c(dy5,rep("",length(m3)-length(dy5)+1))
      for (j in 1:length(m3)) tbl5[j] <- kd5[j]
    }
  })
  addHandlerClicked(gbplus7, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl7[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy7",m3[which(m3 %in% union(dy7,temp))],envir=.GlobalEnv) 
      tbl7[1] <- dy7
    }
  })
  addHandlerClicked(gbminus7, handler = function(h,...) {
    temp <- svalue(tbl7)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("dy7",character(0),envir=.GlobalEnv) 
      tbl7[1] <- ""
    }
  })
  gbutton("Continue", cont = cg2, width=20, handler = function(h, ...) {
    if (length(dy1)==0 | length(dy4)==0 | length(dy5)==0) {gmessage("Missing required information.", parent = toplevel)} else {
    dispose(toplevel)
    if (length(dy7)>0) tmpweight <- dy7 else tmpweight <- NULL   
    cc_calc <- function(x,y) {
      cc <- data.frame(matrix(rep(0,m1*n1*(length(y)+3)),nrow=m1*n1))
      names(cc)[1:2] <- dy4
      names(cc)[3:(length(y)+2)] <- paste("cc",y,"_niche",sep="") 
      names(cc)[length(y)+3] <- "meancc"
      cc[,1] <- rep(1:m1,each=n1)
      cc[,2] <- rep(1:n1,m1)
      for (i in 1:length(y)) {
        k1 <- x[which(x[,i+3]==1),c(1:3,i+3)]
        k2 <- table(k1[,2],k1[,3])
        for (f in 1:(nrow(k2))) {
          for (g in 1:(ncol(k2))) {
            o <- which(cc[,1]==as.numeric(rownames(k2)[f]) & cc[,2]==as.numeric(colnames(k2)[g]))
            cc[o,(i+2)] <- k2[f,g]/sum(k2)
          }
        }
      }
      cc[,(length(y)+3)] <- rowMeans(cc[,3:(length(y)+2)])
      return(cc)
    }
    mr_calc <- function(x,y) {
      mr <- data.frame(matrix(rep(0,m1*n1*(length(y)+3)),nrow=m1*n1))
      names(mr)[1:2] <- dy4
      names(mr)[3:(length(y)+2)] <- paste("mr",y,sep="") 
      names(mr)[length(y)+3] <- "meanmr"
      mr[,1] <- rep(1:m1,each=n1)
      mr[,2] <- rep(1:n1,m1)
      for (i in 1:length(y)) {
        k1 <- x[which(x[,i+length(y)+3]==1),c(1:3,i+length(y)+3)]
        k2 <- table(k1[,2],k1[,3])
        for (f in 1:(nrow(k2))) {
          for (g in 1:(ncol(k2))) {
            o <- which(mr[,1]==as.numeric(rownames(k2)[f]) & mr[,2]==as.numeric(colnames(k2)[g]))
            mr[o,(i+2)] <- k2[f,g]/sum(k2)
          }
        }
      }
      mr[,(length(y)+3)] <- rowMeans(mr[,3:(length(y)+2)])
      return(mr)
    }
    if (length(dy4)==2) {
      single <- function(attr) {
        b <- blau(attr, node.ids=dy1, dimension=dy4, memberships=dy5,weights=tmpweight,complete.cases=dy8)
        b <- niches(b, dev.range = dy6)
        assign("bobj",b,envir=.GlobalEnv)
        k <- data.frame(cbind(b$ids$nodeId,b$dimensions,b$isInNiche,b$memberships))
        names(k)[1] <- "nodeId"
        ccw <- gwindow("Dimension Category Selection", width=600, height=600, parent = window)
        cg1 <- ggroup(cont = ccw, use.scrollwindow=T, horizontal = FALSE)
        dd1 <- data.frame(t(table(k[,2])))[,2:3]
        colnames(dd1)[1] <- dy4[1]
        dd2 <- data.frame(t(table(k[,3])))[,2:3]
        colnames(dd2)[1] <- dy4[2]
        glabel(paste("This is the frequency table for ",dy4[1],".",sep=""),cont=cg1)
        dim1 <- gtable(dd1, expand = TRUE, cont = cg1)
        glabel(paste("Please slide the bar below to set the categories for dimension ",dy4[1],"?",sep=""),cont=cg1)
        assign("m1",1,envir=.GlobalEnv)
        assign("n1",1,envir=.GlobalEnv)
        gslider(from = 1, to = max(k[,2],na.rm = TRUE)-min(k[,2],na.rm = TRUE)+1, by = 1, value = 1, cont=cg1, handler = function(h,...){
          assign("m1",svalue(h$obj),envir=.GlobalEnv)
        })
        gseparator(cont = cg1)
        gdim2 <- glabel(paste("This is the frequency table for ",dy4[2],".",sep=""),cont=cg1)
        dim2 <- gtable(dd2, expand = TRUE, cont = cg1)
        glabel(paste("Please slide the bar below to set the categories for dimension ",dy4[2],"?",sep=""),cont=cg1)
        gslider(from = 1, to = max(k[,3],na.rm = TRUE)-min(k[,3],na.rm = TRUE)+1, by = 1, value = 1, cont=cg1, handler = function(h,...){
          assign("n1",svalue(h$obj),envir=.GlobalEnv)
        })
        gseparator(cont = cg1)
        gbutton("Continue", cont = cg1, width=20, handler = function(h, ...){
          dispose(ccw)  
          x <- seq(min(k[,2],na.rm = TRUE),max(k[,2],na.rm = TRUE),by = (max(k[,2],na.rm = TRUE)-min(k[,2],na.rm = TRUE))/(m1-1))
          y <- seq(min(k[,3],na.rm = TRUE),max(k[,3],na.rm = TRUE),by = (max(k[,3],na.rm = TRUE)-min(k[,3],na.rm = TRUE))/(n1-1)) 
          k[,2] <- cut(k[,2],b=m1,label=c(1:m1)) 
          k[,3] <- cut(k[,3],b=n1,label=c(1:n1)) 
          cc <- cc_calc(k,dy5)
          mr <- mr_calc(k,dy5)
          ie <- cbind(cc[,1],cc[,2],mr[,(length(dy5)+3)]-cc[,(length(dy5)+3)])
          ie_poly <- lm(ie[,3]~ie[,1]+I(ie[,1]^2)+I(ie[,1]^3)+ie[,2]+I(ie[,2]^2)+I(ie[,2]^3)+ie[,1]*ie[,2])
          ie <- cbind(ie,predict(ie_poly))
          final <- data.frame(cbind(cc,mr[,3:ncol(mr)],ie[,3:4]))
          colnames(final)[ncol(final)-1] <- "ie"
          colnames(final)[ncol(final)] <- "ie_poly"
          trdlevel <- gwindow("Plot...", width=20)
          g2 <- gpanedgroup(cont = trdlevel, horizontal = TRUE)
          cg2 <- ggroup(cont = g2, horizontal = TRUE)
          button1 <- gbutton("Plot Carrying Capacity", cont = cg2, width=20, handler = function(h, ...) {
            x1 <- x2 <- y1 <- y2 <- z1 <- z2 <- rep(0,length(dy5))
            for (i in 1:length(dy5)) {
              x1[i] <- x2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[1]],na.rm = TRUE)
              y1[i] <- y2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[2]],na.rm = TRUE)
              z1[i] <- min(cc[,(3+length(dy5))])
              z2[i] <- max(cc[,(3+length(dy5))])
            }
            z <- matrix(cc[,(3+length(dy5))],nrow=length(x))
            cc3D <- function(phi1,theta1) {
              persp3D(x,y,z, facets=F, phi=phi1, theta=theta1,
                             xlab = dy4[1], ylab = dy4[2], zlab = "Carrying capacity",
                             col = ramp.col(c("blue", "red")), colkey = FALSE,
                             bty = "f", ticktype = "detailed")
              scatter3D(x2,y2,z2,add = TRUE, col = "red")
              scatter3D(x1,y1,z1,add = TRUE, col = "red")
              segments3D(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3D(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            cc3Drgl <- function() {
              nbcol = m1*n1
              color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
              zcol  = cut(cc[,(3+length(dy5))], nbcol)
              persp3d(x,y,z, 
                         xlab = dy4[1], ylab = dy4[2], zlab = "Carrying capacity",
                         col=color[zcol])
              par3d(windowRect = c(900, 50, 1500, 650))
              scatter3Drgl(x2,y2,z2,dev = rgl.cur(), col = "red")
              rgl.close()
              scatter3Drgl(x2,y2,z2,add = TRUE, col = "red")
              scatter3Drgl(x1,y1,z1,add = TRUE, col = "red")
              segments3Drgl(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3Drgl(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            cc3Drgl()
            gthlevel <- gwindow("Carrying Capacity Plot", width = 800, height = 800)
            gg <- ggroup(cont = gthlevel, horizontal = F, expand=T)
            ang1 <- ang2 <- "0"
            assign("cc3d1",600,envir=.GlobalEnv)
            assign("cc3d2",600,envir=.GlobalEnv)
            assign("cc3d3",100,envir=.GlobalEnv)
            assign("cc3d4","png",envir=.GlobalEnv)
            tbl <- glayout(cont = gg)
            i <- 1
            tbl[i,i] <- "horizonal angle:"
            tbl[i,2] <- gedit("0") 
            tbl[i,3] <- "vertical angle"
            tbl[i,4] <- gedit("0") 
            tbl[i,5] <- "width"
            tbl[i,6] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("cc3d1",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,7] <- "height"
            tbl[i,8] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("cc3d2",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,9] <- "DPI"
            tbl[i,10] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
              handler = function(h,...){assign("cc3d3",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,11] <- "format"
            tbl[i,12] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
              handler = function(h,...){assign("cc3d4",svalue(h$obj),envir=.GlobalEnv)})
            button <- gbutton("Refresh Plot", cont = gg, handler = function(h, ...) {
              Sys.sleep(0.5)
              cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            })
            button <- gbutton("Save Plot", cont = gg, handler = function(h, ...) {
              if (cc3d4=="png") {
                png("cc3d.png", width = cc3d1, height = cc3d2, res = cc3d3)
                cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (cc3d4=="jpg") {
                jpeg("cc3d.jpg", width = cc3d1, height = cc3d2, res = cc3d3)
                cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (cc3d4=="bmp") {
                bmp("cc3d.bmp", width = cc3d1, height = cc3d2, res = cc3d3)
                cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (cc3d4=="tiff") {
                tiff("cc3d.tiff", width = cc3d1, height = cc3d2, res = cc3d3)
                cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (cc3d4=="pdf") {
                pdf("cc3d.pdf", width = cc3d1/100, height = cc3d2/100)
                cc3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
            })
            plot1 <- ggraphics(container=gg)
            dev1 <- dev.cur();
            Sys.sleep(0.5)
            cc3D(as.numeric(ang1),as.numeric(ang2))
          })
          addSpace(cg2, 60) 
          button2 <- gbutton("Plot Membership Rate", cont = cg2, width=20, handler = function(h, ...) {
            x1 <- x2 <- y1 <- y2 <- z1 <- z2 <- rep(0,length(dy5))
            for (i in 1:length(dy5)) {
              x1[i] <- x2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[1]],na.rm = TRUE)
              y1[i] <- y2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[2]],na.rm = TRUE)
              z1[i] <- min(mr[,(3+length(dy5))])
              z2[i] <- max(mr[,(3+length(dy5))])
            }
            z <- matrix(mr[,(3+length(dy5))],nrow=length(x))
            mr3D <- function(phi1,theta1) {
              persp3D(x,y,z, facets=F, phi=phi1, theta=theta1,
                             xlab = dy4[1], ylab = dy4[2], zlab = "Membership rate",
                           col = ramp.col(c("blue", "red")), colkey = FALSE,
                           bty = "f", ticktype = "detailed")
              scatter3D(x2,y2,z2,add = TRUE, col = "red")
              scatter3D(x1,y1,z1,add = TRUE, col = "red")
              segments3D(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3D(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            mr3Drgl <- function() {
              nbcol = m1*n1
              color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
              zcol  = cut(mr[,(3+length(dy5))], nbcol)
              persp3d(x,y,z, 
                      col=color[zcol], 
                      xlab = dy4[1], ylab = dy4[2], zlab = "Membership rate")
              par3d(windowRect = c(900, 50, 1500, 650))
              scatter3Drgl(x2,y2,z2,dev = rgl.cur(), col = "red")
              rgl.close()
              scatter3Drgl(x2,y2,z2,add = TRUE, col = "red")
              scatter3Drgl(x1,y1,z1,add = TRUE, col = "red")
              segments3Drgl(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3Drgl(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            mr3Drgl()
            hthlevel <- gwindow("Membership Rate Plot", width = 800, height = 800)
            hg <- ggroup(cont = hthlevel, horizontal = F, expand=T)
            ang1 <- ang2 <- "0"
            assign("mr3d1",600,envir=.GlobalEnv)
            assign("mr3d2",600,envir=.GlobalEnv)
            assign("mr3d3",100,envir=.GlobalEnv)
            assign("mr3d4","png",envir=.GlobalEnv)
            tbl <- glayout(cont = hg)
            i <- 1
            tbl[i,i] <- "horizonal angle:"
            tbl[i,2] <- gedit("0") 
            tbl[i,3] <- "vertical angle"
            tbl[i,4] <- gedit("0") 
            tbl[i,5] <- "width"
            tbl[i,6] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("mr3d1",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,7] <- "height"
            tbl[i,8] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("mr3d2",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,9] <- "DPI"
            tbl[i,10] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
              handler = function(h,...){assign("mr3d3",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,11] <- "format"
            tbl[i,12] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
              handler = function(h,...){assign("mr3d4",svalue(h$obj),envir=.GlobalEnv)})
            button <- gbutton("Refresh Plot", cont = hg, handler = function(h, ...) {
              Sys.sleep(0.5)
              mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            })
            button <- gbutton("Save Plot", cont = hg, handler = function(h, ...) {
              if (mr3d4=="png") {
                png("mr3d.png", width = mr3d1, height = mr3d2, res = mr3d3)
                mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (mr3d4=="jpg") {
                jpeg("mr3d.jpg", width = mr3d1, height = mr3d2, res = mr3d3)
                mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (mr3d4=="bmp") {
                bmp("mr3d.bmp", width = mr3d1, height = mr3d2, res = mr3d3)
                mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (mr3d4=="tiff") {
                tiff("mr3d.tiff", width = mr3d1, height = mr3d2, res = mr3d3)
                mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (mr3d4=="pdf") {
                pdf("mr3d.pdf", width = mr3d1/100, height = mr3d2/100)
                mr3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
            })
            plot1 <- ggraphics(container=hg)
            dev1 <- dev.cur();
            Sys.sleep(0.5)
            mr3D(as.numeric(ang1),as.numeric(ang2))
          })
          addSpace(cg2, 60) 
          button3 <- gbutton("Plot Intensity of Exploitation", cont = cg2, width=20, handler = function(h, ...) {
            x1 <- x2 <- y1 <- y2 <- z1 <- z2 <- rep(0,length(dy5))
            for (i in 1:length(dy5)) {
              x1[i] <- x2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[1]],na.rm = TRUE)
              y1[i] <- y2[i] <- mean(attr[which(attr[,dy5[i]]==1),dy4[2]],na.rm = TRUE)
              z1[i] <- min(ie[,4])
              z2[i] <- max(ie[,4])
            }
            z <- matrix(ie[,4],nrow=length(x))
            ie3D <- function(phi1,theta1) {
              persp3D(x,y,z, facets=F, phi=phi1, theta=theta1,
                             xlab = dy4[1], ylab = dy4[2], zlab = "Intensity of exploitation",
                             col = ramp.col(c("blue", "red")), colkey = FALSE,
                             bty = "f", ticktype = "detailed")
              scatter3D(x2,y2,z2,add = TRUE, col = "red")
              scatter3D(x1,y1,z1,add = TRUE, col = "red")
              segments3D(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3D(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            ie3Drgl <- function() {
              nbcol = m1*n1
              color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
              zcol  = cut(ie[,3], nbcol)
              persp3d(x,y,z,
                      col=color[zcol], 
                      xlab = dy4[1], ylab = dy4[2], zlab = "Intensity of exploitation")
              par3d(windowRect = c(900, 50, 1500, 650))
              scatter3Drgl(x2,y2,z2,dev = rgl.cur(), col = "red")
              rgl.close()
              scatter3Drgl(x2,y2,z2,add = TRUE, col = "red")
              scatter3Drgl(x1,y1,z1,add = TRUE, col = "red")
              segments3Drgl(x1,y1,z1,x2,y2,z2,add = TRUE,col="red")
              text3Drgl(x2, y2, z2, labels=dy5,add = TRUE, col = "red", cex=0.9)
            }
            ie3Drgl()
            ithlevel <- gwindow("Intensity of Exploitation Plot", width = 800, height = 800)
            ig <- ggroup(cont = ithlevel, horizontal = F, expand=T)
            ang1 <- ang2 <- "0"
            assign("ie3d1",600,envir=.GlobalEnv)
            assign("ie3d2",600,envir=.GlobalEnv)
            assign("ie3d3",100,envir=.GlobalEnv)
            assign("ie3d4","png",envir=.GlobalEnv)
            tbl <- glayout(cont = ig)
            i <- 1
            tbl[i,i] <- "horizonal angle:"
            tbl[i,2] <- gedit("0") 
            tbl[i,3] <- "vertical angle"
            tbl[i,4] <- gedit("0") 
            tbl[i,5] <- "width"
            tbl[i,6] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("ie3d1",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,7] <- "height"
            tbl[i,8] <- gcombobox(c(600,800,1000,1200,1400,1600,1800,2000), selected = 1, cont = tbl, 
              handler = function(h,...){assign("ie3d2",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,9] <- "DPI"
            tbl[i,10] <- gcombobox(c(100,200,300,400,500,600), selected = 1, cont = tbl, 
              handler = function(h,...){assign("ie3d3",svalue(h$obj),envir=.GlobalEnv)})
            tbl[i,11] <- "format"
            tbl[i,12] <- gcombobox(c("png","jpg","bmp","tiff","pdf"), selected = 1, cont = tbl, 
              handler = function(h,...){assign("ie3d4",svalue(h$obj),envir=.GlobalEnv)})
            button <- gbutton("Refresh Plot", cont = ig, handler = function(h, ...) {
              Sys.sleep(0.5)
              ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
            })
            button <- gbutton("Save Plot", cont = ig, handler = function(h, ...) {
              if (ie3d4=="png") {
                png("ie3d.png", width = ie3d1, height = ie3d2, res = ie3d3)
                ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (ie3d4=="jpg") {
                jpeg("ie3d.jpg", width = ie3d1, height = ie3d2, res = ie3d3)
                ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (ie3d4=="bmp") {
                bmp("ie3d.bmp", width = ie3d1, height = ie3d2, res = ie3d3)
                ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (ie3d4=="tiff") {
                tiff("ie3d.tiff", width = ie3d1, height = ie3d2, res = ie3d3)
                ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
              if (ie3d4=="pdf") {
                pdf("ie3d.pdf", width = ie3d1/100, height = ie3d2/100)
                ie3D(as.numeric(svalue(tbl[i,2])),as.numeric(svalue(tbl[i,4])))
                dev.off()}
            })
            plot1 <- ggraphics(container=ig)
            dev1 <- dev.cur();
            Sys.sleep(0.5)
            ie3D(as.numeric(ang1),as.numeric(ang2))
            gseparator(cont = ig)
            polyf <- data.frame(matrix(ie_poly$coefficients,nrow=1))
            names(polyf) <- c("Intercept",dy4[1],paste(dy4[1],"^2",sep=""),paste(dy4[1],"^3",sep=""),
              dy4[2],paste(dy4[2],"^2",sep=""),paste(dy4[2],"^3",sep=""),paste(dy4[1],"*",dy4[2],sep=""))
            gtable(polyf,cont=ig)
          })
          addSpace(cg2, 60) 
          button4 <- gbutton("Show Table", cont = cg2, width=20, handler = function(h, ...) {
            fourthlevel <- gwindow("Table for Carrying Capacity, Membership Rate, & Intensity of Exploitation",width = 800, height = 600)
            ng4 <- ggroup(horizontal = FALSE, cont = fourthlevel)
            button1 <- gbutton("Save as csv file: cc&mr&ie.csv", expand = FALSE, cont = ng4, handler = function(h, ...) {
              write.table(final, "cc&mr&ie.csv", row.names=F, col.names=T, sep=",")
            })
            button2 <- gbutton("Save as R file: cc&mr&ie.Rdata", expand = FALSE, cont = ng4, handler = function(h, ...) {
              save(final, file="cc&mr&ie.Rdata")
            })
            button3 <- gbutton("Save as SAS file: cc&mr&ie.txt & cc&mr&ie.sas", expand = FALSE, cont = ng4, handler = function(h, ...) {
              write.foreign(final, "cc&mr&ie.txt", "cc&mr&ie.sas",   package="SAS")
            })
            button4 <- gbutton("Save as Stata file: cc&mr&ie.dta", expand = FALSE, cont = ng4, handler = function(h, ...) {
              write.dta(final, ("cc&mr&ie.dta"))
            })
            button5 <- gbutton("Save as SPSS file: cc&mr&ie.txt & cc&mr&ie.sps", expand = FALSE, cont = ng4, handler = function(h, ...) {
              write.foreign(final, "cc&mr&ie.txt", "cc&mr&ie.sps",   package="SPSS")
            })
            gseparator(cont = ng4)
            vars <- gtable(final, expand = TRUE, cont = ng4)
          })
        })
      }
      if (length(dy2)==0) {
        single(cov)
      }
      if (length(dy2)>0) {
        if (nrow(unique(cov[which(colnames(cov)==dy2)]))==1) {
          single(cov)
        } else if (nrow(unique(cov[which(colnames(cov)==dy2)]))>1) {
          dylevel <- gwindow("Niche Dynamics Options", width=600, height=600, parent = window)
          dyg <- ggroup(cont = dylevel, use.scrollwindow=T, horizontal = FALSE)
          tbl <- glayout(cont = dyg)
          j <- 1
          tbl[j,1] <- "Please select which ecology you want to test niche dynamics:"
          tbl[j,2] <- gradio(c("all",as.matrix(unique(cov[which(colnames(cov)==dy2)]))), check=1, handler = function(h,...){
            assign("dy9",svalue(h$obj),envir=.GlobalEnv)
          })
          buttondy <- gbutton("Continue", expand = FALSE, cont = dyg, handler = function(h, ...) {
            dispose(dylevel)
            if (dy9!="all") {
              cov1 <- cov[which(cov[,which(names(cov)==dy2)]==dy9),]
              single(cov1)
            } else {
              z1 <- which(names(cov)==dy1)
              z2 <- which(names(cov)==dy2)
              z3 <- which(names(cov)==dy4[1])
              z4 <- which(names(cov)==dy4[2])
              t <- table(cov[z2])
              k <- data.frame(cbind(cov[z1],cov[z3],cov[z4]))
              ccw <- gwindow("Dimension Category Selection", width=600, height=600, parent = window)
              cg1 <- ggroup(cont = ccw, use.scrollwindow=T, horizontal = FALSE)
              dd1 <- data.frame(t(table(k[,2])))[,2:3]
              colnames(dd1)[1] <- dy4[1]
              dd2 <- data.frame(t(table(k[,3])))[,2:3]
              colnames(dd2)[1] <- dy4[2]
              glabel(paste("This is the frequency table for ",dy4[1],".",sep=""),cont=cg1)
              dim1 <- gtable(dd1, expand = TRUE, cont = cg1)
              glabel(paste("Please slide the bar below to set the categories for dimension ",dy4[1],"?",sep=""),cont=cg1)
              assign("m1",1,envir=.GlobalEnv)
              assign("n1",1,envir=.GlobalEnv)
              gslider(from = 1, to = max(k[,2],na.rm = TRUE)-min(k[,2],na.rm = TRUE)+1, by = 1, value = 1, cont=cg1, handler = function(h,...){
                assign("m1",svalue(h$obj),envir=.GlobalEnv)
              })
              gseparator(cont = cg1)
              gdim2 <- glabel(paste("This is the frequency table for ",dy4[2],".",sep=""),cont=cg1)
              dim2 <- gtable(dd2, expand = TRUE, cont = cg1)
              glabel(paste("Please slide the bar below to set the categories for dimension ",dy4[2],"?",sep=""),cont=cg1)
              gslider(from = 1, to = max(k[,3],na.rm = TRUE)-min(k[,3],na.rm = TRUE)+1, by = 1, value = 1, cont=cg1, handler = function(h,...){
                assign("n1",svalue(h$obj),envir=.GlobalEnv)
              })
              gseparator(cont = cg1)
              gbutton("Continue", cont = cg1, width=20, handler = function(h, ...){
                dispose(ccw) 
                nichem1 <- nichem2 <- matrix(rep(0,(length(t)-1)*length(dy5)),nrow=length(t)-1)
                netop1 <- netop2 <- rightc1 <- rightc2 <- leftc1 <- leftc2 <- matrix(rep(0,length(t)*length(dy5)),nrow=length(t))
                for (o in 1:(length(t)-1)) {
                  a1 <- cov[which(cov[,which(names(cov)==dy2)]==names(t)[o]),]
                  a1[,z3] <- as.numeric(cut(a1[,z3],b=m1,label=c(1:m1))) 
                  a1[,z4] <- as.numeric(cut(a1[,z4],b=n1,label=c(1:n1))) 
                  a2 <- cov[which(cov[,which(names(cov)==dy2)]==names(t)[o+1]),]
                  a2[,z3] <- as.numeric(cut(a2[,z3],b=m1,label=c(1:m1))) 
                  a2[,z4] <- as.numeric(cut(a2[,z4],b=n1,label=c(1:n1))) 
                  for(p in 1:length(dy5)) {
                    nichem1[o,p] <- mean(a2[which(a2[,p+3]==1),z3],na.rm=T)-mean(a1[which(a1[,p+3]==1),z3],na.rm=T)
                    nichem2[o,p] <- mean(a2[which(a2[,p+3]==1),z4],na.rm=T)-mean(a1[which(a1[,p+3]==1),z4],na.rm=T) 
                  }
                }
                for (o in 1:length(t)) {
                   a <- cov[which(cov[,which(names(cov)==dy2)]==names(t)[o]),]
                   bo <- blau(a, node.ids=dy1, dimension=dy4, memberships=dy5,weights=tmpweight,complete.cases=dy8)
                   bo <- niches(bo, dev.range = dy6)
                   ko <- data.frame(cbind(bo$ids$nodeId,bo$dimensions,bo$isInNiche,bo$memberships))
                   cco <- cc_calc(ko,dy5)
                   mro <- mr_calc(ko,dy5)
                   ieo <- data.frame(cbind(cco[,1],cco[,2],mro[,(length(dy5)+3)]-cco[,(length(dy5)+3)]))
                   names(ieo) <- c("x1","x2","y")
                   ie_polyo <- lm(y~x1+I(x1^2)+I(x1^3)+x2+I(x2^2)+I(x2^3)+x1*x2,data=ieo)
                   cuto <- seq(from = dy6/10, to = dy6, by = dy6/10)
                   for(p in 1:length(dy5)) {
                     mean1 <- mean(a[which(a[,which(names(a)==dy5[p])]==1),z3],na.rm=T)
                     sd1 <- sd(a[which(a[,which(names(a)==dy5[p])]==1),z3],na.rm=T)
                     mean2 <- mean(a[which(a[,which(names(a)==dy5[p])]==1),z4],na.rm=T)
                     sd2 <- sd(a[which(a[,which(names(a)==dy5[p])]==1),z4],na.rm=T)
                     r1 <- mean1+cuto*sd1
                     iep <- data.frame(r1,rep(mean2,length(r1)))
                     names(iep) <- c("x1","x2")
                     rightc1[o,p] <- sum(predict(ie_polyo,iep))
                     l1 <- mean1-cuto*sd1
                     iep <- data.frame(l1,rep(mean2,length(l1)))
                     names(iep) <- c("x1","x2")
                     leftc1[o,p] <- sum(predict(ie_polyo,iep))     
                     netop1[o,p] <- leftc1[o,p]-rightc1[o,p]
                     r2 <- mean2+cuto*sd2
                     iep <- data.frame(rep(mean1,length(r2)),r2)
                     names(iep) <- c("x1","x2")
                     rightc2[o,p] <- sum(predict(ie_polyo,iep))
                     l2 <- mean2-cuto*sd2
                     iep <- data.frame(rep(mean1,length(l2)),l2)
                     names(iep) <- c("x1","x2")
                     leftc2[o,p] <- sum(predict(ie_polyo,iep))     
                     netop2[o,p] <- leftc2[o,p]-rightc2[o,p]
                   }
                 }
                 niche_movement <- data.frame(cbind(c(nichem1),c(netop1[1:(length(t)-1),]),c(nichem2),c(netop2[1:(length(t)-1),])))
                 names(niche_movement)[1] <- "Niche_movement_dim1"
                 names(niche_movement)[2] <- "Niche_opportunity_dim1"
                 names(niche_movement)[3] <- "Niche_movement_dim2"
                 names(niche_movement)[4] <- "Niche_opportunity_dim2"
                 assign("niche_movement",niche_movement,envir=.GlobalEnv)
                 eq3 <- lm(Niche_movement_dim1~Niche_opportunity_dim1,data=niche_movement)
                 eq4 <- lm(Niche_movement_dim2~Niche_opportunity_dim2,data=niche_movement)
                 outeq3 <- paste(capture.output(summary(eq3)), collapse="\n")
                 outeq4 <- paste(capture.output(summary(eq4)), collapse="\n")
                 eqw <- gwindow("Predicted Niche Movement", width=800, height=750, parent = window)
                 eqg1 <- ggroup(cont = eqw, use.scrollwindow=T, horizontal = FALSE)
                 button1 <- gbutton("Save as csv file: niche_movement.csv", expand = FALSE, cont = eqg1, handler = function(h, ...) {
                   write.table(niche_movement, "niche_movement.csv", row.names=F, col.names=T, sep=",")
                 })
                 button2 <- gbutton("Save as R file: niche_movement.Rdata", expand = FALSE, cont = eqg1, handler = function(h, ...) {
                   save(niche_movement, file="niche_movements.Rdata")
                 })
                 button3 <- gbutton("Save as SAS file: niche_movement.txt & niche_movement.sas", expand = FALSE, cont = eqg1, handler = function(h, ...) {
                   write.foreign(niche_movement, "niche_movements.txt", "niche_movements.sas",   package="SAS")
                 })
                 button4 <- gbutton("Save as Stata file: niche_movement.dta", expand = FALSE, cont = eqg1, handler = function(h, ...) {
                   write.dta(niche_movement, ("niche_movements.dta"))
                 })
                 button5 <- gbutton("Save as SPSS file: niche_movement.txt & niche_movement.sps", expand = FALSE, cont = eqg1, handler = function(h, ...) {
                   write.foreign(niche_movement, "niche_movement.txt", "niche_movements.sps",   package="SPSS")
                 })
                 gseparator(cont = eqg1)
                 glabel(paste("dim1 = ",dy4[1],sep=""),cont=eqg1)
                 gtext(outeq3, cont=eqg1, expand=TRUE, font.attr=c(family="monospace"))
                 glabel(paste("dim2 = ",dy4[2],sep=""),cont=eqg1)
                 gtext(outeq4, cont=eqg1, expand=TRUE, font.attr=c(family="monospace"))
              })
            }
          })
        }
      }
    } else {
      gmessage("You must select 2 dimensions.")
    }
  }})}
}

