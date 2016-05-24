showanalysis <- function(h,...) {
  if ("cov" %in% ls(envir=.GlobalEnv)==FALSE) {gmessage("Sorry! Attribute file is not loaded.", parent = window)} else {
  assign("b1",character(0),envir=.GlobalEnv)
  assign("b2",character(0),envir=.GlobalEnv)
  assign("b3","",envir=.GlobalEnv)
  assign("b4",character(0),envir=.GlobalEnv)
  assign("b5",character(0),envir=.GlobalEnv)
  assign("b6",1.5,envir=.GlobalEnv)
  assign("b7",character(0),envir=.GlobalEnv)
  assign("b8","FALSE",envir=.GlobalEnv)
  assign("m1",names(cov),envir=.GlobalEnv)
  assign("m3",names(cov),envir=.GlobalEnv)
  toplevel <- gwindow("Blaunet Analysis", width=800, height=800, parent = window)
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
  b1temp <- data.frame(Node.ids="",stringsAsFactors=FALSE)
  b2temp <- data.frame(Ecology.ids="",stringsAsFactors=FALSE)
  b4temp <- data.frame(Dimensions=rep("",length(m3)),stringsAsFactors=FALSE)
  b5temp <- data.frame(Groups=rep("",length(m3)),stringsAsFactors=FALSE)
  b7temp <- data.frame(Weights="",stringsAsFactors=FALSE)
  tbl1 <- gtable(b1temp,expand=TRUE,multiple=FALSE,cont=cg2)
  tbl2 <- gtable(b2temp,expand=TRUE,multiple=FALSE,cont=cg2)
  if ('el' %in% ls(envir=.GlobalEnv)) {
    gcheckboxgroup("Network included",cont=cg2,handler = function(h,...) assign("b3",svalue(h$obj),envir=.GlobalEnv))
  }
  tbl4 <- gtable(b4temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl5 <- gtable(b5temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl7 <- gtable(b7temp,expand=TRUE,multiple=FALSE,cont=cg2)
  glabel("Dev.range",cont=cg2)
  gslider(from = 0, to = 5, by = .05, value = 1.5, cont=cg2, handler = function(h,...) assign("b6",svalue(h$obj),envir=.GlobalEnv))
  glabel("Complete.cases",cont=cg2)
  gradio(c("TRUE","FALSE"), selected = 2, cont = cg2,  handler = function(h,...) assign("b8",svalue(h$obj),envir=.GlobalEnv))
  addHandlerClicked(gbplus1, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl1[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b1",m3[which(m3 %in% union(b1,temp))],envir=.GlobalEnv) 
      tbl1[1] <- b1
    }
  })
  addHandlerClicked(gbminus1, handler = function(h,...) {
    temp <- svalue(tbl1)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b1",character(0),envir=.GlobalEnv) 
      tbl1[1] <- ""
    }
  })
  addHandlerClicked(gbplus2, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl2[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b2",m3[which(m3 %in% union(b2,temp))],envir=.GlobalEnv) 
      tbl2[1] <- b2
    }
  })
  addHandlerClicked(gbminus2, handler = function(h,...) {
    temp <- svalue(tbl2)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b2",character(0),envir=.GlobalEnv) 
      tbl2[1] <- ""
    }
  })
  addHandlerClicked(gbplus4, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b4",m3[which(m3 %in% union(b4,temp))],envir=.GlobalEnv) 
      kb4 <- c(b4,rep("",length(m3)-length(b4)))
      for (j in 1:length(m3)) tbl4[j] <- kb4[j]
    }
  })
  addHandlerClicked(gbminus4, handler = function(h,...) {
    temp <- svalue(tbl4)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b4",m3[which(m3 %in% setdiff(b4,temp))],envir=.GlobalEnv) 
      kb4 <- c(b4,rep("",length(m3)-length(b4)+1))
      for (j in 1:length(m3)) tbl4[j] <- kb4[j]
    }
  })
  addHandlerClicked(gbplus5, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b5",m3[which(m3 %in% union(b5,temp))],envir=.GlobalEnv) 
      kb5 <- c(b5,rep("",length(m3)-length(b5)))
      for (j in 1:length(m3)) tbl5[j] <- kb5[j]
    }
  })
  addHandlerClicked(gbminus5, handler = function(h,...) {
    temp <- svalue(tbl5)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b5",m3[which(m3 %in% setdiff(b5,temp))],envir=.GlobalEnv) 
      kb5 <- c(b5,rep("",length(m3)-length(b5)+1))
      for (j in 1:length(m3)) tbl5[j] <- kb5[j]
    }
  })
  addHandlerClicked(gbplus7, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl7[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b7",m3[which(m3 %in% union(b7,temp))],envir=.GlobalEnv) 
      tbl7[1] <- b7
    }
  })
  addHandlerClicked(gbminus7, handler = function(h,...) {
    temp <- svalue(tbl7)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("b7",character(0),envir=.GlobalEnv) 
      tbl7[1] <- ""
    }
  })
  button <- gbutton("Continue", cont = cg2, handler = function(h, ...) {
    if (length(b1)==0 | length(b4)==0 | length(b5)==0) {gmessage("Missing required information.", parent = toplevel)} else {
    if (length(b7)>0) tmpweight <- b7 else tmpweight <- NULL
    if (length(b2)==0 & b3!="Network included") {
      b <- blau(cov, node.ids=b1, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8)
      b <- niches(b, dev.range = b6)
    }
    if (length(b2)==0 & b3=="Network included") {
      b <- blau(cov, node.ids=b1, graph = el, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8)
      b <- niches(b, dev.range = b6)
    }
    if (length(b2)!=0 & b3!="Network included") {
      b <- blau(cov, node.ids=b1, ecology.ids=b2, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8 )
      b <- niches(b, dev.range = b6)
    }
    if (length(b2)!=0 & b3=="Network included") {
      b <- blau(cov, node.ids=b1, ecology.ids=b2, graph = el, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8 )
      b <- niches(b, dev.range = b6)
    }
    if (b3=="Network included") {
      b <- nodal.network(b)
    }
    assign("bobj",b,envir=.GlobalEnv)
    dispose(toplevel)
    sndlevel <- gwindow("Niche Analysis Results", width=20)
    g <- gpanedgroup(cont = sndlevel, horizontal = TRUE)
    cg <- ggroup(cont = g, horizontal = TRUE)
    button <- gbutton("Show Object", cont = cg, handler = function(h, ...) {
      blauobj <- data.frame(capture.output(b))
      nw1 <- gwindow("Object Information", width = 800, height = 600, scrollwindow=TRUE)
      ng1 <- ggroup(horizontal = FALSE, cont = nw1 )
      button <- gbutton("Save Blau object: blauobj.Rdata", expand = FALSE, cont = ng1, handler = function(h, ...) {
        save(blauobj, file="blauobj.Rdata")
      })
      gseparator(cont = ng1)
      gdf(blauobj, expand = TRUE, cont = ng1)
    })
    addSpace(cg, 60) 
    button <- gbutton("Nodal Result", cont = cg, handler = function(h, ...) {
      b <- nodal.global(b)
      Nodalstatus <- b$nodalGlobal[,2]
      nstemp1 <- nstemp2 <- nstemp3 <- rep(0,nrow(cov))
      nstemp1[as.numeric(b$nodalGlobal[,2])==0] <- 1
      nstemp2[as.numeric(b$nodalGlobal[,2])==1] <- 1
      nstemp3[as.numeric(b$nodalGlobal[,2])>1] <- 1
      b0 <- data.frame(matrix(rep(NA,length(b$memberships)),nrow=nrow(b$memberships)))
      for (i in 1:nrow(b$memberships)) {
        for (j in 1:length(b5)) {
          if (is.na(b$memberships[i,j])==FALSE & is.na(b$isInNiche[i,j])==FALSE) {
            if (b$memberships[i,j]==1 & b$isInNiche[i,j]==1) b0[i,j] <- "Member & nicher" 
            if (b$memberships[i,j]==1 & b$isInNiche[i,j]==0) b0[i,j] <- "Member not in niche" 
            if (b$memberships[i,j]==0 & b$isInNiche[i,j]==1) b0[i,j] <- "Non-member but in niche"
            if (b$memberships[i,j]==0 & b$isInNiche[i,j]==0) b0[i,j] <- "Neither member nor nicher"
          }     
        }
      }
      colnames(b0) <- paste(b5,"_mem&niche",sep="")
      if (b3!="Network included") {
        nodaloutput <- data.frame(cbind(b$ids,b$nodalGlobal,nstemp1,nstemp2,nstemp3,b$memberships,b$isInNiche[,1:(ncol(b$isInNiche)-1)],b0))
      } else {
        nodaloutput <- data.frame(cbind(b$ids,b$nodalGlobal,nstemp1,nstemp2,nstemp3,b$memberships,b$isInNiche[,1:(ncol(b$isInNiche)-1)],b0,b$nodalNetwork))
      }
      names(nodaloutput)[4] <- "Niches"
      names(nodaloutput)[6] <- "Outsider"
      names(nodaloutput)[7] <- "Insider_Exclusive"
      names(nodaloutput)[8] <- "Insider_Manifolder"
      nw2 <- gwindow("Nodal Result",width = 800, height = 600)
      ng2 <- ggroup(horizontal = FALSE, cont = nw2)
      button1 <- gbutton("Save as csv file: nodaloutput.csv", expand = FALSE, cont = ng2, handler = function(h, ...) {
        write.table(nodaloutput, "nodaloutput.csv", row.names=F, col.names=T, sep=",")
      })
      button2 <- gbutton("Save as R file: nodaloutput.Rdata", expand = FALSE, cont = ng2, handler = function(h, ...) {
        save(nodaloutput, file="nodaloutput.Rdata")
      })
      button3 <- gbutton("Save as SAS file: nodaloutput.txt & nodaloutput.sas", expand = FALSE, cont = ng2, handler = function(h, ...) {
        write.foreign(nodaloutput, "nodaloutput.txt", "nodaloutput.sas",   package="SAS")
      })
      button4 <- gbutton("Save as Stata file: nodaloutput.dta", expand = FALSE, cont = ng2, handler = function(h, ...) {
        write.dta(nodaloutput, ("nodaloutput.dta"))
      })
      button5 <- gbutton("Save as SPSS file: nodaloutput.txt & nodaloutput.sps", expand = FALSE, cont = ng2, handler = function(h, ...) {
        write.foreign(nodaloutput, "nodaloutput.txt", "nodaloutput.sps",   package="SPSS")
      })
      gseparator(cont = ng2)
      vars <- gtable(nodaloutput, expand = TRUE, cont = ng2)
    })
    addSpace(cg, 60) 
    button <- gbutton("Niche Breadth Summary", cont = cg, handler = function(h, ...) {
      nbtemp <- data.frame(bobj$lowbounds,stringsAsFactors = FALSE)
      for (o in 1:nrow(nbtemp)) {
        for (p in 1:length(b4)) {
           nbtemp[o,p] <- paste(format(round(bobj$lowbounds[o,p],2),nsmall=2),"-",format(round(bobj$topbounds[o,p],2),nsmall=2),sep="")
        }
      }
      nichebreadth <- cbind(rownames(nbtemp),nbtemp)
      names(nichebreadth)[1] <- "GROUPS"
      nw13 <- gwindow("Niche Breadth Summary",width = 800, height = 600)
      ng13 <- ggroup(horizontal = FALSE, cont = nw13)
      button1 <- gbutton("Save as csv file: nichebreadth.csv", expand = FALSE, cont = ng13, handler = function(h, ...) {
        write.table(nichebreadth, "nichebreadth.csv", row.names=F, col.names=T, sep=",")
      })
      button2 <- gbutton("Save as R file: nichebreadth.Rdata", expand = FALSE, cont = ng13, handler = function(h, ...) {
        save(nichebreadth, file="nichebreadth.Rdata")
      })
      button3 <- gbutton("Save as SAS file: nichebreadth.txt & nichebreadth.sas", expand = FALSE, cont = ng13, handler = function(h, ...) {
        write.foreign(nichebreadth, "nichebreadth.txt", "nichebreadth.sas",   package="SAS")
      })
      button4 <- gbutton("Save as Stata file: nichebreadth.dta", expand = FALSE, cont = ng13, handler = function(h, ...) {
        write.dta(nichebreadth, ("nichebreadth.dta"))
      })
      button5 <- gbutton("Save as SPSS file: nichebreadth.txt & nichebreadth.sps", expand = FALSE, cont = ng13, handler = function(h, ...) {
        write.foreign(nichebreadth, "nichebreadth.txt", "nichebreadth.sps",   package="SPSS")
      })
      gseparator(cont = ng13)
      vars <- gtable(nichebreadth, expand = TRUE, cont = ng13)
    })
    addSpace(cg, 60)
    button <- gbutton("Focal Niche Summary", cont = cg, handler = function(h, ...) {
      fns <- data.frame(niche.summary(b),stringsAsFactors = FALSE)
      es2 <- data.frame(ecology.summary(b,percent=TRUE),stringsAsFactors = FALSE)
      if (length(b2)==0) {
        pm <- c()
        fnstemp <- matrix(as.numeric(fns[,3]),ncol=1)
        es2temp <- es2[,3:(2+length(b5))]
        for (p in 1:length(b5)) {
          pm <- c(pm,round(matrix(as.numeric(es2temp[p,]),nrow=1) %*% fnstemp))
        }
      }
      if (length(b2)==1) {
        pm <- c()
        for (o in unique(es2$Ecology)) {
          fnstemp <- matrix(as.numeric(fns[which(fns$Ecology==o),3]),ncol=1)
          es2temp <- es2[which(es2$Ecology==o),3:(2+length(b5))]
          for (p in 1:length(b5)) {
            pm <- c(pm,round(matrix(as.numeric(es2temp[p,]),nrow=1) %*% fnstemp))
          }
        }
      }
      nichesummary <- cbind(fns[,1:4],pm,fns[,5:7],round(as.numeric(fns[,5])/as.numeric(fns[,4]),4)*100)
      names(nichesummary)[5] <- "PredictedNicheMem"
      names(nichesummary)[9] <- "ExclusivePercent"
      nw5 <- gwindow("Focal Niche Summary",width = 800, height = 600)
      ng5 <- ggroup(horizontal = FALSE, cont = nw5)
      button1 <- gbutton("Save as csv file: nichesummary.csv", expand = FALSE, cont = ng5, handler = function(h, ...) {
        write.table(nichesummary, "nichesummary.csv", row.names=F, col.names=T, sep=",")
      })
      button2 <- gbutton("Save as R file: nichesummary.Rdata", expand = FALSE, cont = ng5, handler = function(h, ...) {
        save(nichesummary, file="nichesummary.Rdata")
      })
      button3 <- gbutton("Save as SAS file: nichesummary.txt & nichesummary.sas", expand = FALSE, cont = ng5, handler = function(h, ...) {
        write.foreign(nichesummary, "nichesummary.txt", "nichesummary.sas",   package="SAS")
      })
      button4 <- gbutton("Save as Stata file: nichesummary.dta", expand = FALSE, cont = ng5, handler = function(h, ...) {
        write.dta(nichesummary, ("nichesummary.dta"))
      })
      button5 <- gbutton("Save as SPSS file: nichesummary.txt & nichesummary.sps", expand = FALSE, cont = ng5, handler = function(h, ...) {
        write.foreign(nichesummary, "nichesummary.txt", "nichesummary.sps",   package="SPSS")
      })
      gseparator(cont = ng5)
      vars <- gtable(nichesummary, expand = TRUE, cont = ng5)
    })
    addSpace(cg, 60) 
    button <- gbutton("Niche by Niche Summary", cont = cg, handler = function(h, ...) {
      es1 <- data.frame(ecology.summary(b),stringsAsFactors = FALSE)
      es2 <- data.frame(ecology.summary(b,percent=TRUE),stringsAsFactors = FALSE)
      names(es2)[3:(2+length(b5))] <- paste(names(es2)[3:(2+length(b5))],"_CC",sep="")
      if (length(b2)==0) {
        AVG_CC <- mean(as.numeric(as.matrix(es2[,3:(2+length(b5))])),na.rm=TRUE)
        STD_CC <- sd(as.numeric(as.matrix(es2[,3:(2+length(b5))])),na.rm=TRUE)
        AVG_CC <- rep(AVG_CC,length(b5))
        STD_CC <- rep(STD_CC,length(b5)) 
        ecologysummary <- cbind(es1,es2[,3:(2+length(b5))],AVG_CC,STD_CC)
      }
      if (length(b2)==1) {
        AVG_CC <- STD_CC <- c()
        for (o in unique(es2$Ecology)) {
          es2temp <- es2[which(es2$Ecology==o),]
          AVG_CC <- c(AVG_CC,mean(as.numeric(as.matrix(es2temp[,3:(2+length(b5))])),na.rm=TRUE))
          STD_CC <- c(STD_CC,sd(as.numeric(as.matrix(es2temp[,3:(2+length(b5))])),na.rm=TRUE))
        }
        AVG_CC <- rep(AVG_CC,each=length(b5))
        STD_CC <- rep(STD_CC,each=length(b5))
        ecologysummary <- cbind(es1,es2[,3:(2+length(b5))],AVG_CC,STD_CC)
      }
      for (o in 1:nrow(ecologysummary)) {
        for (p in (3+length(b5)):ncol(ecologysummary)) {
           ecologysummary[o,p] <- format(round(as.numeric(ecologysummary[o,p]),2),nsmall=2)
        }
      }
      nw4 <- gwindow("Niche by Niche Summary",width = 800, height = 600)
      ng4 <- ggroup(horizontal = FALSE, cont = nw4)
      button1 <- gbutton("Save as csv file: ecologysummary.csv", expand = FALSE, cont = ng4, handler = function(h, ...) {
        write.table(ecologysummary, "ecologysummary.csv", row.names=F, col.names=T, sep=",")
      })
      button2 <- gbutton("Save as R file: ecologysummary.Rdata", expand = FALSE, cont = ng4, handler = function(h, ...) {
        save(ecologysummary, file="ecologysummary.Rdata")
      })
      button3 <- gbutton("Save as SAS file: ecologysummary.txt & ecologysummary.sas", expand = FALSE, cont = ng4, handler = function(h, ...) {
        write.foreign(ecologysummary, "ecologysummary.txt", "ecologysummary.sas",   package="SAS")
      })
      button4 <- gbutton("Save as Stata file: ecologysummary.dta", expand = FALSE, cont = ng4, handler = function(h, ...) {
        write.dta(ecologysummary, ("ecologysummary.dta"))
      })
      button5 <- gbutton("Save as SPSS file: ecologysummary.txt & ecologysummary.sps", expand = FALSE, cont = ng4, handler = function(h, ...) {
        write.foreign(ecologysummary, "ecologysummary.txt", "ecologysummary.sps",   package="SPSS")
      })
      gseparator(cont = ng4)
      vars <- gtable(ecologysummary, expand = TRUE, cont = ng4)
    })
    addSpace(cg, 60)
    if (b3=="Network included") {
      button <- gbutton("Dyadic Result", cont = cg, handler = function(h, ...) {
        b <- dyadic(b)
        dyadicoutput <- data.frame(b$dyadic)
        nw6 <- gwindow("Dyadic Result",width = 800, height = 600)
        ng6 <- ggroup(horizontal = FALSE, cont = nw6)
        button1 <- gbutton("Save as csv file: dyadicoutput.csv", expand = FALSE, cont = ng6, handler = function(h, ...) {
          write.table(dyadicoutput, "dyadicoutput.csv", row.names=F, col.names=T, sep=",")
        })
        button2 <- gbutton("Save as R file: dyadicoutput.Rdata", expand = FALSE, cont = ng6, handler = function(h, ...) {
          save(dyadicoutput, file="dyadicoutput.Rdata")
        })
        button3 <- gbutton("Save as SAS file: dyadicoutput.txt & dyadicoutput.sas", expand = FALSE, cont = ng6, handler = function(h, ...) {
          write.foreign(dyadicoutput, "dyadicoutput.txt", "dyadicoutput.sas",   package="SAS")
        })
        button4 <- gbutton("Save as Stata file: dyadicoutput.dta", expand = FALSE, cont = ng6, handler = function(h, ...) {
          write.dta(dyadicoutput, "dyadicoutput.dta")
        })
        button5 <- gbutton("Save as SPSS file: dyadicoutput.txt & dyadicoutput.sps", expand = FALSE, cont = ng6, handler = function(h, ...) {
          write.foreign(dyadicoutput, "dyadicoutput.txt", "dyadicoutput.sps",   package="SPSS")
        })
        gseparator(cont = ng6)
        vars <- gtable(dyadicoutput, expand = TRUE, cont = ng6)
      })
      addSpace(cg, 60) 
    }
    if (b3!="Network included") {
      button <- gbutton("Correlation Matrix", cont = cg, handler = function(h, ...) {
        if (length(b2)==0) {
          b <- blau(cov, node.ids=b1, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8)
          b <- niches(b, dev.range = b6)
        }
        if (length(b2)!=0) {
          b <- blau(cov, node.ids=b1, ecology.ids=b2, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8 )
          b <- niches(b, dev.range = b6)
        }
        b <- nodal.global(b)
        y1 <- cbind(data.frame(rownames(b$nodalGlobal)),data.frame(matrix(as.numeric(b$nodalGlobal[,1:2]),ncol=2)))
        names(y1)[1] <- "nodeId"
        names(y1)[2] <- "TotalOrgs"
        names(y1)[3] <- "Niches"
        y2 <- cbind(cov[,which(b1==names(cov))],cov[,b4[1:length(b4)]])
        names(y2)[1] <- "nodeId"
        x <- merge(y1,y2,by="nodeId")
        x <- x[,-1]
        x <- cor(x)
        coroutput <- data.frame(cbind(rownames(x),x))
        names(coroutput)[1] <- "Names"
        nw7 <- gwindow("Correlation Matrix",width = 800, height = 600)
        ng7 <- ggroup(horizontal = FALSE, cont = nw7)
        button1 <- gbutton("Save as csv file: correlation.csv", expand = FALSE, cont = ng7, handler = function(h, ...) {
          write.table(coroutput, "correlation.csv", row.names=F, col.names=T, sep=",")
        })
        button2 <- gbutton("Save as R file: correlation.Rdata", expand = FALSE, cont = ng7, handler = function(h, ...) {
          save(coroutput, file="correlation.Rdata")
        })
        button3 <- gbutton("Save as SAS file: correlation.txt & correlation.sas", expand = FALSE, cont = ng7, handler = function(h, ...) {
          write.foreign(coroutput, "correlation.txt", "correlation.sas",   package="SAS")
        })
        button4 <- gbutton("Save as Stata file: correlation.dta", expand = FALSE, cont = ng7, handler = function(h, ...) {
          write.dta(coroutput, ("correlation.dta"))
        })
        button5 <- gbutton("Save as SPSS file: correlation.txt & correlation.sps", expand = FALSE, cont = ng7, handler = function(h, ...) {
          write.foreign(coroutput, "correlation.txt", "correlation.sps",   package="SPSS")
        })
        gseparator(cont = ng7)
        vars <- gtable(coroutput, expand = TRUE, cont = ng7)
      })
    }
    if (b3=="Network included") {
      button <- gbutton("Correlation Matrix", cont = cg, handler = function(h, ...) {
        if (length(b2)==0) {
          b <- blau(cov, node.ids=b1, graph = el, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8)
          b <- niches(b, dev.range = b6)
        }
        if (length(b2)!=0) {
          b <- blau(cov, node.ids=b1, ecology.ids=b2, graph = el, dimension=b4, memberships=b5,weights=tmpweight,complete.cases = b8 )
          b <- niches(b, dev.range = b6)
        }
        b <- nodal.global(b)
        y1 <- cbind(data.frame(rownames(b$nodalGlobal)),data.frame(matrix(as.numeric(b$nodalGlobal[,1:2]),ncol=2)))
        names(y1)[1] <- "nodeId"
        names(y1)[2] <- "TotalOrgs"
        names(y1)[3] <- "Niches"
        y2 <- cbind(cov[,which(b1==names(cov))],cov[,b4[1:length(b4)]])
        names(y2)[1] <- "nodeId"
        x <- merge(y1,y2,by="nodeId")
        n <- network(adj)
        outdegree <- degree(n,cmode="outdegree")
        indegree <- degree(n,cmode="indegree")
        betweenness <- betweenness(n)
        closeness <- closeness(n)
        eigenvector <- round(evcent(n),4)
        y3 <- cbind(data.frame(rownames(b$nodalGlobal)),data.frame(cbind(outdegree,indegree,betweenness,closeness,eigenvector)))
        names(y3)[1] <- "nodeId"
        x <- merge(x,y3,by="nodeId")
        x <- x[,-1]
        x <- cor(x)
        coroutput <- data.frame(cbind(rownames(x),x))
        names(coroutput)[1] <- "Names"
        nw8 <- gwindow("Correlation Matrix",width = 800, height = 600)
        ng8 <- ggroup(horizontal = FALSE, cont = nw8)
        button1 <- gbutton("Save as csv file: correlation.csv", expand = FALSE, cont = ng8, handler = function(h, ...) {
          write.table(coroutput, "correlation.csv", row.names=F, col.names=T, sep=",")
        })
        button2 <- gbutton("Save as R file: correlation.Rdata", expand = FALSE, cont = ng8, handler = function(h, ...) {
          save(coroutput, file="correlation.Rdata")
        })
        button3 <- gbutton("Save as SAS file: correlation.txt & correlation.sas", expand = FALSE, cont = ng8, handler = function(h, ...) {
          write.foreign(coroutput, "correlation.txt", "correlation.sas",   package="SAS")
        })
        button4 <- gbutton("Save as Stata file: correlation.dta", expand = FALSE, cont = ng8, handler = function(h, ...) {
          write.dta(coroutput, ("correlation.dta"))
        })
        button5 <- gbutton("Save as SPSS file: correlation.txt & correlation.sps", expand = FALSE, cont = ng8, handler = function(h, ...) {
          write.foreign(coroutput, "correlation.txt", "correlation.sps",   package="SPSS")
        })
        gseparator(cont = ng8)
        vars <- gtable(coroutput, expand = TRUE, cont = ng8)
      })
    }
  }})}
}

