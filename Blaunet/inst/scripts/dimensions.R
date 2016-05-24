showdimensions <- function(h,...) {
  if ("cov" %in% ls(envir=.GlobalEnv)==FALSE) {gmessage("Sorry! Attribute file is not loaded.", parent = window)} else {
  assign("c1",character(0),envir=.GlobalEnv)
  assign("c2",character(0),envir=.GlobalEnv)
  assign("c3","all",envir=.GlobalEnv)
  assign("c4",character(0),envir=.GlobalEnv)
  assign("c5",character(0),envir=.GlobalEnv)
  assign("c6",0.05,envir=.GlobalEnv)
  assign("c7",character(0),envir=.GlobalEnv)
  assign("c8",character(0),envir=.GlobalEnv)
  assign("c9",character(0),envir=.GlobalEnv)
  assign("m1",names(cov),envir=.GlobalEnv)
  assign("m3",names(cov),envir=.GlobalEnv)
  toplevel <- gwindow("Salient Dimensions", width=600, height=600, parent = window)
  tbl0 <- gtable(m1,expand=TRUE,multiple=TRUE,cont=toplevel)
  cg <- gpanedgroup(width=100, cont = toplevel)
  cg1 <- ggroup(horizontal = FALSE, cont = cg)
  addSpace(cg1,20,horizontal = FALSE)
  gbplus1 <- gbutton("+", cont=cg1)
  gbminus1 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus2 <- gbutton("+", cont=cg1)
  gbminus2 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus4 <- gbutton("+", cont=cg1)
  gbminus4 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  gbplus5 <- gbutton("+", cont=cg1)
  gbminus5 <- gbutton("-", cont=cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  addSpring(cg1)
  cg2 <- gframe("Options", horizontal=FALSE, cont=cg)
  c1temp <- data.frame(Node.ids="",stringsAsFactors=FALSE)
  c2temp <- data.frame(Ecology.ids="",stringsAsFactors=FALSE)
  c4temp <- data.frame(Dimensions=rep("",length(m3)),stringsAsFactors=FALSE)
  c5temp <- data.frame(Groups=rep("",length(m3)),stringsAsFactors=FALSE)
  tbl1 <- gtable(c1temp,expand=TRUE,multiple=FALSE,cont=cg2)
  tbl2 <- gtable(c2temp,expand=TRUE,multiple=FALSE,cont=cg2)
  tbl4 <- gtable(c4temp,expand=TRUE,multiple=TRUE,cont=cg2)
  tbl5 <- gtable(c5temp,expand=TRUE,multiple=TRUE,cont=cg2)
  glabel("Alpha (0.05 by default)",cont=cg2)
  gcombobox(c(0.05,0.001,0.01,0.1), cont=cg2, handler = function(h,...) assign("c6",svalue(h$obj),envir=.GlobalEnv))
  addHandlerClicked(gbplus1, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl1[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c1",m3[which(m3 %in% union(c1,temp))],envir=.GlobalEnv) 
      tbl1[1] <- c1
    }
  })
  addHandlerClicked(gbminus1, handler = function(h,...) {
    temp <- svalue(tbl1)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c1",character(0),envir=.GlobalEnv) 
      tbl1[1] <- ""
    }
  })
  addHandlerClicked(gbplus2, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & tbl2[1]=="" & length(temp)==1) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c2",m3[which(m3 %in% union(c2,temp))],envir=.GlobalEnv) 
      tbl2[1] <- c2
    }
  })
  addHandlerClicked(gbminus2, handler = function(h,...) {
    temp <- svalue(tbl2)
    if ("" %in% temp==FALSE & length(temp)==1) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c2",character(0),envir=.GlobalEnv) 
      tbl2[1] <- ""
    }
  })
  addHandlerClicked(gbplus4, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c4",m3[which(m3 %in% union(c4,temp))],envir=.GlobalEnv) 
      kc4 <- c(c4,rep("",length(m3)-length(c4)))
      for (j in 1:length(m3)) tbl4[j] <- kc4[j]
    }
  })
  addHandlerClicked(gbminus4, handler = function(h,...) {
    temp <- svalue(tbl4)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c4",m3[which(m3 %in% setdiff(c4,temp))],envir=.GlobalEnv) 
      kc4 <- c(c4,rep("",length(m3)-length(c4)+1))
      for (j in 1:length(m3)) tbl4[j] <- kc4[j]
    }
  })
  addHandlerClicked(gbplus5, handler = function(h,...) {
    temp <- svalue(tbl0)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% setdiff(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c5",m3[which(m3 %in% union(c5,temp))],envir=.GlobalEnv) 
      kc5 <- c(c5,rep("",length(m3)-length(c5)))
      for (j in 1:length(m3)) tbl5[j] <- kc5[j]
    }
  })
  addHandlerClicked(gbminus5, handler = function(h,...) {
    temp <- svalue(tbl5)
    if ("" %in% temp==FALSE & length(temp)>0) {
      assign("m1",m3[which(m3 %in% union(m1,temp))],envir=.GlobalEnv) 
      km1 <- c(m1,rep("",length(m3)-length(m1)+1))
      for (i in 1:length(m3)) tbl0[i] <- km1[i]
      assign("c5",m3[which(m3 %in% setdiff(c5,temp))],envir=.GlobalEnv) 
      kc5 <- c(c5,rep("",length(m3)-length(c5)+1))
      for (j in 1:length(m3)) tbl5[j] <- kc5[j]
    }
  })
  button <- gbutton("Continue", expand = FALSE, cont = cg2, handler = function(h, ...) {
    if (length(c1)==0 | length(c4)==0) {gmessage("Missing required information.", parent = toplevel)} else {
    dispose(toplevel)
    if ('el' %in% ls(envir=.GlobalEnv)==FALSE & length(c5)==0) {
      gmessage("You must select group options.")
    } else {
      if (length(c5)>0) {
        if (length(c2)!=0) {
          seclevel <- gwindow("Salient Dimension Options", width=600, height=600, parent = window)
          secg <- ggroup(cont = seclevel, use.scrollwindow=T, horizontal = FALSE)
          tbl <- glayout(cont = secg)
          j <- 1
          tbl[j,1] <- "Please select which ecology you want to identify salient dimensions:"
          tbl[j,2] <- gcheckboxgroup(c("all",as.matrix(unique(cov[which(colnames(cov)==c2)]))), handler = function(h,...){
            assign("c3",svalue(h$obj),envir=.GlobalEnv)
          })
          j <- j + 1
          tbl[j,1] <- "Please identify categorical variables:"
          tbl[j,2] <- gcheckboxgroup(c4, handler = function(h,...){
            assign("c7",svalue(h$obj),envir=.GlobalEnv)
          })
          j <- j + 1        
          rm(j)
          button <- gbutton("Continue", expand = FALSE, cont = secg, handler = function(h, ...) {
            dispose(seclevel)
            if ("all" %in% c3) {
              cov1 <- cov
              if ('el' %in% ls(envir=.GlobalEnv)) adj1 <- adj
            } else {
              cov1 <- cov[which(cov[,which(colnames(cov)==c2)] %in% c3),]
              if ('el' %in% ls(envir=.GlobalEnv)) adj1 <- adj[which(cov[,which(colnames(cov)==c2)] %in% c3),which(cov[,which(colnames(cov)==c2)] %in% c3)]
            }
            if ('el' %in% ls(envir=.GlobalEnv)) {
              tempadj <- symmetrize(adj1,rule="weak")
              deg <- degree(network(tempadj), cmode="outdegree")
              bet <- betweenness(network(tempadj))
              eig <- evcent(network(tempadj))
              tempel <- as.matrix(network(tempadj),matrix.type="edgelist")
              lclu <- rep(0,nrow(cov1))
              for (i in 1:nrow(cov1)) {
                k <- tempel[which(tempel[,1]==i),2]
                if (length(k)>=2) {
                  m <- rbind(tempel, t(combn(k,2)))
                  lclu[i] <- nrow(m[duplicated(m), , drop = FALSE])/nrow(t(combn(k,2)))
                } else {lclu[i] <- NA}
              }
              rm(i,k,m)
              k <- c()
              for (i in 1:length(c5)) {
                alteringroup <- rep(0,nrow(cov1))
                for (j in 1:nrow(cov1)) {
                  if (length(tempel[which(tempel[,1]==j),2])>=1) alteringroup[j] <- sum(cov1[tempel[which(tempel[,1]==j),2],which(colnames(cov1)==c5[i])])
                  else alteringroup[j] <- 0
                }
                temp <- cbind(cov1[which(colnames(cov1)==c5[i])],cov1[which(colnames(cov1) %in% c5[-i])],cov1[which(colnames(cov1) %in% c4)],deg,bet,eig,lclu,alteringroup)
                fit1 <- glm(temp, family=binomial)
                k <- unique(c(k,names(which(summary(fit1)$coefficients[,4]<c6))))
              }
              assign("c8",c4[which(c4 %in% k)],envir=.GlobalEnv) 
              if (length(c8)>0) {
                g <- c()
                gm <- list()
                for (i in 1:length(c5)) {
                  g <- c(g,which(names(cov)==c5[i]))
                  gm[[i]] <- matrix(rep(0,nrow(cov)*nrow(cov)),nrow=nrow(cov))
                }
                for (l in 1:length(c5)) {
                  tempgm <- cbind(c(1:nrow(cov)),cov1[,g[l]])
                  tempgm <- tempgm[which(tempgm[,2]==1),]
                  tempgm <- t(combn(tempgm[,1],2))
                  tempgm <- network(tempgm,directed=FALSE)
                  tempgm1 <- nrow(cov)-network.size(tempgm)
                  if (tempgm1>0) add.vertices(tempgm,tempgm1)
                  gm[[l]] <- as.matrix(tempgm)
                }
                rm(i,l)
                for (l in 1:length(c5)) {
                  n <- network(gm[[l]])
                  emf <- "n~"
                  for (m in 1:length(c8)) {
                    k1 <- cov1[,which(names(cov1)==c8[m])]
                    k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                    n %v% c8[m] <- k1
                    if (c8[m] %in% c7) {
                      emf <- paste(emf,"+nodematch(","\'",c8[m],"\'",")",sep="")                    
                    } else {
                      emf <- paste(emf,"+absdiff(","\'",c8[m],"\'",")",sep="")                    
                    }
                  }
                  o <- setdiff(union(l,1:length(c5)),intersect(l,1:length(c5)))
                  if (nrow(cov1)<=500) {
                    for (m in 1:length(o)) emf <- paste(emf,"+edgecov(gm[[",o[m],"]])",sep="")
                    emf <- paste(emf,"+edgecov(tempadj)",sep="")
                  }
                  emf <- paste(substr(emf, 1, 2), substr(emf, 4, nchar(emf)), sep='')
                  em <- ergm(as.formula(emf))
                  c9 <- union(c9,c8[which(summary(em)$coefs[1:length(c8),4]<c6)])
                }
                assign("c9",c9,envir=.GlobalEnv) 
              } else {assign("c9","",envir=.GlobalEnv)}
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"  
              thdlevel <- gwindow("Salient Dimensions",width = 200, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg)  
            } else {
              k <- c()
              for (j in 1:length(c5)) {
                temp <- cbind(cov1[which(colnames(cov1)==c5[j])],cov1[which(colnames(cov1) %in% c5[-j])],cov1[which(colnames(cov1) %in% c4)])
                fit1 <- glm(temp, family=binomial)
                k <- unique(c(k,names(which(summary(fit1)$coefficients[,4]<c6))))
              }
              assign("c8",c4[which(c4 %in% k)],envir=.GlobalEnv) 
              if (length(c8)>0) {
                g <- c()
                gm <- list()
                for (i in 1:length(c5)) {
                  g <- c(g,which(names(cov1)==c5[i]))
                  gm[[i]] <- matrix(rep(0,nrow(cov1)*nrow(cov1)),nrow=nrow(cov1))
                }
                for (l in 1:length(c5)) {
                  tempgm <- cbind(c(1:nrow(cov)),cov1[,g[l]])
                  tempgm <- tempgm[which(tempgm[,2]==1),]
                  tempgm <- t(combn(tempgm[,1],2))
                  tempgm <- network(tempgm,directed=FALSE)
                  tempgm1 <- nrow(cov)-network.size(tempgm)
                  if (tempgm1>0) add.vertices(tempgm,tempgm1)
                  gm[[l]] <- as.matrix(tempgm)
                }
                rm(i,l)
                for (l in 1:length(c5)) {
                  n <- network(gm[[l]])
                  emf <- "n~"
                  for (m in 1:length(c8)) {
                    k1 <- cov1[,which(names(cov1)==c8[m])]
                    k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                    n %v% c8[m] <- k1
                    if (c8[m] %in% c7) {
                      emf <- paste(emf,"+nodematch(","\'",c8[m],"\'",")",sep="")                    
                    } else {
                      emf <- paste(emf,"+absdiff(","\'",c8[m],"\'",")",sep="")                    
                    }
                  }
                  o <- setdiff(union(l,1:length(c5)),intersect(l,1:length(c5)))
                  if (nrow(cov1)<=500) {
                    for (m in 1:length(o)) emf <- paste(emf,"+edgecov(gm[[",o[m],"]])",sep="")
                  }
                  emf <- paste(substr(emf, 1, 2), substr(emf, 4, nchar(emf)), sep='')
                  em <- ergm(as.formula(emf))
                  c9 <- union(c9,c8[which(summary(em)$coefs[1:length(c8),4]<c6)])
                }
                assign("c9",c9,envir=.GlobalEnv) 
              } else {assign("c9","",envir=.GlobalEnv)}
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"  
              thdlevel <- gwindow("Salient Dimensions",width = 200, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg)   
            }
          })
        } else {
          seclevel <- gwindow("Salient Dimension Options", width=600, height=400, parent = window)
          secg <- ggroup(cont = seclevel, use.scrollwindow=T, horizontal = FALSE)
          tbl <- glayout(cont = secg)
          j <- 1
          tbl[j,1] <- "Please identify categorical variables:"
          tbl[j,2] <- gcheckboxgroup(c4, handler = function(h,...){
            assign("c7",svalue(h$obj),envir=.GlobalEnv)
          })
          rm(j) 
          button <- gbutton("Continue", expand = FALSE, cont = secg, handler = function(h, ...) {
            dispose(seclevel)
            cov1 <- cov
            if ('el' %in% ls(envir=.GlobalEnv)) {
              adj1 <- adj
              tempadj <- symmetrize(adj1,rule="weak")
              deg <- degree(network(tempadj), cmode="outdegree")
              bet <- betweenness(network(tempadj))
              eig <- evcent(network(tempadj))
              tempel <- as.matrix(network(tempadj),matrix.type="edgelist")
              lclu <- rep(0,nrow(cov1))
              for (i in 1:nrow(cov1)) {
                k <- tempel[which(tempel[,1]==i),2]
                if (length(k)>=2) {
                  m <- rbind(tempel, t(combn(k,2)))
                  lclu[i] <- nrow(m[duplicated(m), , drop = FALSE])/nrow(t(combn(k,2)))
                } else {lclu[i] <- NA}
              }
              rm(i,k,m)
              k <- c()
              for (i in 1:length(c5)) {
                alteringroup <- rep(0,nrow(cov1))
                for (j in 1:nrow(cov1)) {
                  if (length(tempel[which(tempel[,1]==j),2])>=1) alteringroup[j] <- sum(cov1[tempel[which(tempel[,1]==j),2],which(colnames(cov1)==c5[i])])
                  else alteringroup[j] <- 0
                }
                temp <- cbind(cov1[which(colnames(cov1)==c5[i])],cov1[which(colnames(cov1) %in% c5[-i])],cov1[which(colnames(cov1) %in% c4)],deg,bet,eig,lclu,alteringroup)
                fit1 <- glm(temp, family=binomial)
                k <- unique(c(k,names(which(summary(fit1)$coefficients[,4]<c6))))
              }
              assign("c8",c4[which(c4 %in% k)],envir=.GlobalEnv)
              if (length(c8)>0) { 
                g <- c()
                gm <- list()
                for (i in 1:length(c5)) {
                  g <- c(g,which(names(cov)==c5[i]))
                  gm[[i]] <- matrix(rep(0,nrow(cov)*nrow(cov)),nrow=nrow(cov))
                }
                for (l in 1:length(c5)) {
                  tempgm <- cbind(c(1:nrow(cov)),cov1[,g[l]])
                  tempgm <- tempgm[which(tempgm[,2]==1),]
                  tempgm <- t(combn(tempgm[,1],2))
                  tempgm <- network(tempgm,directed=FALSE)
                  tempgm1 <- nrow(cov)-network.size(tempgm)
                  if (tempgm1>0) add.vertices(tempgm,tempgm1)
                  gm[[l]] <- as.matrix(tempgm)
                }
                rm(i,l)
                for (l in 1:length(c5)) {
                  n <- network(gm[[l]])
                  emf <- "n~"
                  for (m in 1:length(c8)) {
                    k1 <- cov1[,which(names(cov1)==c8[m])]
                    k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                    n %v% c8[m] <- k1
                    if (c8[m] %in% c7) {
                      emf <- paste(emf,"+nodematch(","\'",c8[m],"\'",")",sep="")                    
                    } else {
                      emf <- paste(emf,"+absdiff(","\'",c8[m],"\'",")",sep="")                    
                    }
                  }
                  o <- setdiff(union(l,1:length(c5)),intersect(l,1:length(c5)))
                  if (nrow(cov1)<=500) {
                    for (m in 1:length(o)) emf <- paste(emf,"+edgecov(gm[[",o[m],"]])",sep="")
                    emf <- paste(emf,"+edgecov(as.matrix(network(tempel)))",sep="")
                  }
                  emf <- paste(substr(emf, 1, 2), substr(emf, 4, nchar(emf)), sep='')
                  em <- ergm(as.formula(emf))
                  c9 <- union(c9,c8[which(summary(em)$coefs[1:length(c8),4]<c6)])
                }
                assign("c9",c9,envir=.GlobalEnv) 
              } else {assign("c9",character(0),envir=.GlobalEnv)}
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"          
              thdlevel <- gwindow("Salient Dimensions",width = 300, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg)    
            } else {
              k <- c()
              for (j in 1:length(c5)) {
                temp <- cbind(cov1[which(colnames(cov1)==c5[j])],cov1[which(colnames(cov1) %in% c5[-j])],cov1[which(colnames(cov1) %in% c4)])
                fit1 <- glm(temp, family=binomial)
                k <- unique(c(k,names(which(summary(fit1)$coefficients[,4]<c6))))
              }
              assign("c8",c4[which(c4 %in% k)],envir=.GlobalEnv) 
              if (length(c8)>0) {
                g <- c()
                gm <- list()
                for (i in 1:length(c5)) {
                  g <- c(g,which(names(cov)==c5[i]))
                  gm[[i]] <- matrix(rep(0,nrow(cov)*nrow(cov)),nrow=nrow(cov))
                }
                for (l in 1:length(c5)) {
                  tempgm <- cbind(c(1:nrow(cov)),cov1[,g[l]])
                  tempgm <- tempgm[which(tempgm[,2]==1),]
                  tempgm <- t(combn(tempgm[,1],2))
                  tempgm <- network(tempgm,directed=FALSE)
                  tempgm1 <- nrow(cov)-network.size(tempgm)
                  if (tempgm1>0) add.vertices(tempgm,tempgm1)
                  gm[[l]] <- as.matrix(tempgm)
                }
                rm(i,l)
                for (l in 1:length(c5)) {
                  n <- network(gm[[l]])
                  emf <- "n~"
                  for (m in 1:length(c8)) {
                    k1 <- cov1[,which(names(cov1)==c8[m])]
                    k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                    n %v% c8[m] <- k1
                    if (c8[m] %in% c7) {
                      emf <- paste(emf,"+nodematch(","\'",c8[m],"\'",")",sep="")                    
                    } else {
                      emf <- paste(emf,"+absdiff(","\'",c8[m],"\'",")",sep="")                    
                    }
                  }
                  o <- setdiff(union(l,1:length(c5)),intersect(l,1:length(c5)))
                  if (nrow(cov1)<=500) {
                    for (m in 1:length(o)) emf <- paste(emf,"+edgecov(gm[[",o[m],"]])",sep="")
                  }
                  emf <- paste(substr(emf, 1, 2), substr(emf, 4, nchar(emf)), sep='')
                  em <- ergm(as.formula(emf))
                  c9 <- union(c9,c8[which(summary(em)$coefs[1:length(c8),4]<c6)])
                }
                assign("c9",c9,envir=.GlobalEnv) 
              } else {assign("c9",character(0),envir=.GlobalEnv)}
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"  
              thdlevel <- gwindow("Salient Dimensions",width = 300, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg)       
            }
          })
        } 
	} else {
        if (length(c2)!=0) {
          seclevel <- gwindow("Salient Dimension Options", width=600, height=600, parent = window)
          secg <- ggroup(cont = seclevel, use.scrollwindow=T, horizontal = FALSE)
          tbl <- glayout(cont = secg)
          j <- 1
          tbl[j,1] <- "Please select which ecology you want to identify salient dimensions:"
          tbl[j,2] <- gradio(c("all",as.matrix(unique(cov[which(colnames(cov)==c2)]))), selected = 1, handler = function(h,...){
            assign("c3",svalue(h$obj),envir=.GlobalEnv)
          })
          j <- j + 1
          tbl[j,1] <- "Please identify categorical variables:"
          tbl[j,2] <- gcheckboxgroup(c4, handler = function(h,...){
            assign("c7",svalue(h$obj),envir=.GlobalEnv)
          })
          j <- j + 1        
          rm(j)
          button <- gbutton("Continue", expand = FALSE, cont = secg, handler = function(h, ...) {
            dispose(seclevel)
            if (c3=="all") {
              cov1 <- cov
              if ('adj' %in% ls(envir=.GlobalEnv)) adj1 <- adj
            } else {
              cov1 <- cov[which(cov[which(colnames(cov)==c2)]==c3),]
              if ('adj' %in% ls(envir=.GlobalEnv)) adj1 <- adj[which(cov[which(colnames(cov)==c2)]==c3),which(cov[which(colnames(cov)==c2)]==c3)]
            }
            if ('el' %in% ls(envir=.GlobalEnv)) {
              n <- network(adj1)
              emf <- "n~edges+mutual"
              for (m in 1:length(c4)) {
                k1 <- cov1[,which(names(cov1)==c4[m])]
                k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                n %v% c4[m] <- k1
                if (c4[m] %in% c7) {
                  emf <- paste(emf,"+nodematch(","\'",c4[m],"\'",")",sep="")                    
                } else {
                  emf <- paste(emf,"+absdiff(","\'",c4[m],"\'",")",sep="")                    
                }	
              }				
              em <- ergm(as.formula(emf))
              c9 <- c4[which(summary(em)$coefs[3:(length(c4)+2),4]<c6)]
              assign("c9",c9,envir=.GlobalEnv) 
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"  
              thdlevel <- gwindow("Salient Dimensions",width = 300, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg) 
            }
          })
        } else {
          seclevel <- gwindow("Salient Dimension Options", width=600, height=600, parent = window)
          secg <- ggroup(cont = seclevel, use.scrollwindow=T, horizontal = FALSE)
          tbl <- glayout(cont = secg)
          j <- 1
          tbl[j,1] <- "Please identify categorical variables:"
          tbl[j,2] <- gcheckboxgroup(c4, handler = function(h,...){
            assign("c7",svalue(h$obj),envir=.GlobalEnv)
          })
          j <- j + 1        
          rm(j)
          button <- gbutton("Continue", expand = FALSE, cont = secg, handler = function(h, ...) {
            dispose(seclevel)
            if ('el' %in% ls(envir=.GlobalEnv)) {
              adj1 <- adj
              n <- network(adj1)
              emf <- "n~edges+mutual"
              for (m in 1:length(c4)) {
                k1 <- cov1[,which(names(cov1)==c4[m])]
                k1[which(is.na(k1)==TRUE)] <- round(mean(k1,na.rm=TRUE))
                n %v% c4[m] <- k1
                if (c4[m] %in% c7) {
                  emf <- paste(emf,"+nodematch(","\'",c4[m],"\'",")",sep="")                    
                } else {
                  emf <- paste(emf,"+absdiff(","\'",c4[m],"\'",")",sep="")                    
                }	
              }				
              em <- ergm(as.formula(emf))
              c9 <- c4[which(summary(em)$coefs[3:(length(c4)+2),4]<c6)]
              assign("c9",c9,envir=.GlobalEnv) 
              k <- data.frame(c9)
              colnames(k) <- "Dimensions"  
              thdlevel <- gwindow("Salient Dimensions",width = 300, height = 400)
              tg <- ggroup(horizontal = FALSE, cont = thdlevel)
              dims <- gtable(k, expand = TRUE, cont = tg) 
            }
          })
	  }    
      }  
    }
  }})}
}

