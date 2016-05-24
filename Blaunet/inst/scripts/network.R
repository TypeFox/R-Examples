
showinfo <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    out <- capture.output(summary(network(adj)))[1:19]
    window <- gwindow("Network Information", visible = FALSE)
    exp_group <- gexpandgroup("Summary", cont = window)
    label <- glabel(out, cont = exp_group)
    visible(exp_group) <- TRUE
    visible(window) <- TRUE} else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showdensity <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    gmessage(paste("The network density is ",network.density(network(adj)),".",sep=""), parent = window)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showcentrality <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    n <- network(adj)
    outd <- degree(n,cmode="outdegree")
    ind <- degree(n,cmode="indegree")
    bet <- betweenness(n)
    clo <- closeness(n)
    eig <- round(evcent(n),4)
    output <- data.frame(cbind(network.vertex.names(n),outd,ind,bet,clo,eig))
    names(output) <- c("id","outdegree","indegree","betweenness","closeness","eigenvector")
    nw <- gwindow("Centrality", width = 600, height = 400)
    group <- ggroup(horizontal = FALSE, cont = nw)
    button1 <- gbutton("Save as csv file: centrality.csv", expand = FALSE, cont = group, handler = function(h, ...) {
      write.table(output, "centrality.csv", row.names=F, col.names=T, sep=",")
    })
    button2 <- gbutton("Save as R file: centrality.Rdata", expand = FALSE, cont = group, handler = function(h, ...) {
      save(output, file="centrality.Rdata")
    })
    button3 <- gbutton("Save as SAS file: centrality.txt & centrality.sas", expand = FALSE, cont = group, handler = function(h, ...) {
      write.foreign(output, "centrality", "centrality.sas",   package="SAS")
    })
    button4 <- gbutton("Save as Stata file: centrality.dta", expand = FALSE, cont = group, handler = function(h, ...) {
      write.dta(output, ("centrality.dta"))
    })
    button5 <- gbutton("Save as SPSS file: centrality.txt & centrality.sps", expand = FALSE, cont = group, handler = function(h, ...) {
      write.foreign(output, "centrality.txt", "centrality.sps",   package="SPSS")
    })
    gseparator(cont = group)
    vars <- gtable(output, expand = TRUE, cont = group)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showdcensus <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    out <- data.frame(dyad.census(network(adj)))
    names(out) <- c("Mutual","Asymmetric","Null")
    window <- gwindow("Dyad Census", width = 300, height = 100)
    vars <- gtable(out, cont = window)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showreciprocityindex <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    n <- (network(adj))
    ri <- (dyad.census(n)[1]*2)/(dyad.census(n)[1]*2+dyad.census(n)[2]) 
    gmessage(paste("The reciprocity index is ",ri,".",sep=""), parent = window)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showtcensus <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    out <- data.frame(triad.census(network(adj)))
    names(out) <- gsub("X","",names(out))
    window <- gwindow("Triad Census", width = 800, height = 100)
    vars <- gtable(out, cont = window)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showglobalcustering <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    t<-triad.census(network(adj))
    num <- 3*(t[9]+t[10]+t[12]+t[13]+t[14]+t[15]+t[16])
    dem <- num+t[4]+t[5]+t[6]+t[7]+t[8]+t[11]
    gmessage(paste("The global clustering coefficient is ",num/dem,".",sep=""), parent = window)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}
showlocalcustering <- function(h,...) {
  if ("adj" %in% ls(envir=.GlobalEnv)) {
    m1 <- as.matrix(network(adj),matrix.type='edgelist')
    m2 <- rbind(m1,cbind(m1[,2],m1[,1]))
    m2 <- m2[!duplicated(m2), ]
    m9 <- matrix(rep(0,2*attr(m1,"n")),ncol=2)
    for (i in 1:attr(m1,"n")) {
      m9[i,1] <- attr(m1,"vnames")[i]
      k <- m2[which(m2[,1]==i),2]
      if (length(k)>=2) {
        m3 <- t(combn(k,2))
        m4 <- rbind(m2, m3)
        m9[i,2] <- nrow(m4[duplicated(m4), , drop = FALSE])/nrow(m3)
      } else {
        m9[i,2] <- NA
      }
    }
    m9 <- data.frame(m9)
    names(m9) <- c("nodeId","Local clustering coefficient")
    nw <- gwindow("Local Clustering Coefficient",width = 800, height = 600)
    group <- ggroup(horizontal = FALSE, cont = nw)
    button1 <- gbutton("Save as csv file: localclustering.csv", expand = FALSE, cont = group, handler = function(h, ...) {
      write.table(m9, "localclustering.csv", row.names=F, col.names=T, sep=",")
    })
    button2 <- gbutton("Save as R file: localclustering.Rdata", expand = FALSE, cont = group, handler = function(h, ...) {
      save(m9, file="localclustering.Rdata")
    })
    button3 <- gbutton("Save as SAS file: localclustering.txt & localclustering.sas", expand = FALSE, cont = group, handler = function(h, ...) {
      write.foreign(m9, "localclustering.txt", "localclustering.sas",   package="SAS")
    })
    button4 <- gbutton("Save as Stata file: localclustering.dta", expand = FALSE, cont = group, handler = function(h, ...) {
      write.dta(m9, ("localclustering.dta"))
    })
    button5 <- gbutton("Save as SPSS file: localclustering.txt & localclustering.sps", expand = FALSE, cont = group, handler = function(h, ...) {
      write.foreign(m9, "localclustering.txt", "localclustering.sps",   package="SPSS")
    })
    gseparator(cont = group)
    vars <- gtable(m9, expand = TRUE, cont = group)
  } else gmessage("Sorry! Network file is not loaded.", parent = window)
}

