#'Plotting functions used in the manuscript
#'
#'@rdname plot_kernscr
#'@keywords internal
#'@importFrom graphics layout par axis plot mtext text points abline legend
#'@export
plot_kernscr_methodsplit <- function(raw_melted, adj_melted, kernel, method,
                                     pathway_names=TRUE, title=NULL, raw_lower_threshold=round(log10(1/10000), 1), adj_lower_threshold=round(log10(1/10000*70), 1)){
  range <- seq(70,1)

  if(pathway_names){
    nplot <- length(method) + 1
  }else{
    nplot <- length(method)
  }

  tmpy <- log10(raw_melted[raw_melted$kern==kernel & raw_melted$variable==method[1], "value"])
  raw_tmpy <- pmax(tmpy, raw_lower_threshold)


  graphics::layout(matrix(1:nplot, nrow=1), c(2.5, rep(2, nplot-1)), rep(1,nplot))

  graphics::par(mar=c(5,0,1,0))
  if(pathway_names){
    graphics::plot(raw_tmpy,range,type="n",pch=15,ylab="",xlab="",
                   yaxt="n",xaxt="n",xlim=c(-3,0),bty="n",cex=0.7)
    txt.pathway <- adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[1], "pathway_name"]
    for(jj in range){
      graphics::text(0,jj, txt.pathway[70 - jj + 1], srt = 0, xpd = TRUE, pos = 2, cex=0.8)
    }
    method_legend <- method
    if(sum(method_legend=="McDm")>0){
      method_legend[method_legend=="McDm"] <- "SCR"
    }
    graphics::legend("bottom", pch=c(4,15), rev(method_legend))
  }

  graphics::par(mar=c(5,0,1,4))
  graphics::plot(raw_tmpy, range, pch=15,ylab="", xlab="log10 raw p-value",
                 yaxt="n", xaxt="n", xlim=c(-4,0), bty="n", cex=0.7, cex.lab=1.1)

  breaks_raw <- floor(raw_lower_threshold):0
  breaks_raw[1] <- raw_lower_threshold
  graphics::axis(side=1, at = breaks_raw, labels = c("", "-3", "-2", "-1", "0"))
  graphics::axis(side=1, at = breaks_raw[1], labels=c(paste("<", breaks_raw[1])), hadj=1)

  raw_tmpy <- pmax(log10(raw_melted[raw_melted$kern==kernel & raw_melted$variable==method[2], "value"]), -4)
  graphics::points(raw_tmpy, range, pch = 4, cex = 0.8)
  graphics::abline(v = log10(0.05), col = "gray");
  #graphics::mtext(side=3, text="Raw p-values", cex=0.95,line=0.7, adj=0.6)

  tmpy <- log10(adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[1], "value"])
  adj_tmpy = pmax(tmpy, adj_lower_threshold)
  graphics::plot(adj_tmpy, range, pch = 15, ylab = "", xlab = "log10 adjusted p-value",
                 yaxt = "n", xaxt = "n", xlim = c(adj_lower_threshold, 0), bty = "n", cex = 0.7, cex.lab = 1.1)
  breaks_adj <- floor(adj_lower_threshold):0
  breaks_adj[1] <- adj_lower_threshold
  graphics::axis(side = 1, at = breaks_adj, labels = c("", as.character(breaks_adj[-1])), tick=TRUE)
  graphics::axis(side = 1, at = breaks_adj[1],labels = c(paste("<", breaks_adj[1])), hadj = 1)

  adj_tmpy <- pmax(log10(adj_melted[adj_melted$kern==kernel & adj_melted$variable==method[2], "value"]), adj_lower_threshold)
  graphics::points(adj_tmpy, range, pch = 4, cex = 0.8)
  graphics::abline(v = log10(0.05), col = "gray")
  #graphics::mtext(side=3, text="Ajusted p-values", cex=0.95,line=0.7, adj=0.6)
}

#'@rdname plot_kernscr
#'@keywords internal
#'@export
plot_kernscr_kernelsplit <- function(raw_melted, adj_melted, kernel, method,
                                     pathway_names = TRUE, title = NULL, raw_lower_threshold = round(log10(1/10000), 1),
                                     adj_lower_threshold = NULL){

  if(!is.null(adj_lower_threshold)){
    warning("'adj_lower_threshold' argument is ignored")
  }
  range <- seq(70,1)

  if(pathway_names){
    nplot <- length(kernel)+1
  }else{
    nplot <- length(kernel)
  }

  if(is.null(title)){
    title <- method
  }

  graphics::layout(matrix(1:nplot, nrow=1), c(3, rep(2, nplot-1)), rep(1, nplot))
  if(title==""){
    graphics::par(mar=c(5,0,1,0))
  }else{
    graphics::par(mar=c(5,0,5,0))
  }

  dolegend <- TRUE

  for(k in kernel){

    tmpy <- log10(adj_melted[adj_melted$kern==k & adj_melted$variable==method, "value"])
    tmpy <- pmax(tmpy, adj_lower_threshold)

    if(pathway_names){
      graphics::plot(tmpy,range,type="n",pch=15,ylab="",xlab="",
                     yaxt="n",xaxt="n",xlim=c(-3,0),bty="n",cex=0.7)
      txt.pathway <- as.character(adj_melted[adj_melted$kern==k & adj_melted$variable==method, "pathway_name"])
      #Cai 2011 paper order:
      #ordralpha <- order(txt.pathway)
      #txt.pathway <- txt.pathway[ordralpha]
      #ordr <- rev(order(sapply(txt.pathway, nchar)))
      #txt.pathway <- txt.pathway[ordr]
      #txt.pathway <- txt.pathway[c(1:13, 15:16, 14, 17:19, 21:22, 20, 23:24, 29:31, 25:28, 34:35, 32:33, 36:70)]
      for(jj in range){
        graphics::text(0,jj, txt.pathway[70 - jj + 1], srt = 0, xpd = TRUE, pos = 2, cex=0.57)
      }
      pathway_names=FALSE
    }

    if(dolegend){
      graphics::legend("bottom", pch=c(4,15), c("raw p-val.", "adj. p-val."), cex=0.85)
      dolegend <- FALSE
    }

    if(title==""){
      par(mar=c(5,0,1,4))
    }else{
      par(mar=c(5,0,5,4))
    }
    graphics::plot(tmpy,range,pch=15,ylab="",xlab="log10 Pvalue",
                   yaxt="n",xaxt="n",xlim=c(adj_lower_threshold,0),bty="n",cex=0.)
    graphics::axis(side=1, at = c(adj_lower_threshold, -3, -2, -1, 0),
                   labels=c(paste0("< ", adj_lower_threshold),"-3", "-2","-1","0"))
    tmpy <- pmax(log10(raw_melted[raw_melted$kern==k & raw_melted$variable==method, "value"]),
                 adj_lower_threshold)
    graphics::points(tmpy,range,pch=4,cex=0.7)
    tmpy <- pmax(log10(adj_melted[adj_melted$kern==k & raw_melted$variable==method, "value"]),
                 adj_lower_threshold)
    graphics::points(tmpy,range,pch=15,cex=0.62)
    graphics::abline(v = log10(0.05), col="gray")
    if(title != ""){
      graphics::mtext(side=3, text=paste(paste0(toupper(substring(k,1,1)), substring(k,2)), "kernel"),
            cex=0.95, line=0.7, adj=0.6)
    }
  }
  if(title!=""){
    graphics::mtext(side=3, text=title, cex=1, line=3, at=-5)
  }
  #graphics::mtext(side=1, text="log10 Pvalue", cex=0.9, line=3.5, adj=0.5, at=-5)
}