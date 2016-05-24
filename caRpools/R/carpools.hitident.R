carpools.hitident = function(data, type="deseq2", title="DESeq2 plot",print.names=FALSE, cutoff=c(0,0,0,0), inches=0.1, offsetplot=1.2, plot.p=0.01, sgRNA.top=1, separate=FALSE)
{
  
  # Prepare dataset depending on analysis method
  # wilcox
  if(type=="wilcox")
   { # data from wilcox or riger
    res.short=as.data.frame(
      t(
        apply(
          as.data.frame(data),
          1,
          function(x) c(
            as.numeric(x["p.value"]),
            abs(log2(as.numeric(x["foldchange"]))),
            sign(log2(as.numeric(x["foldchange"])))
          )
        )
      )
    )}
   #  DESeqV2
    else if (type=="deseq2")
    {
      res.short=as.data.frame(
        t(
          apply(
            as.data.frame(data$genes),
            1,
            function(x) c(
              if(is.finite(as.numeric(x["padj"]))) {as.numeric(x["padj"])},
              abs(as.numeric(x["log2FoldChange"])),
              sign(as.numeric(x["log2FoldChange"])),
              as.numeric(x["sgRNA"])
            )
          )
        )
      )
    
    }
  
  
  # plotting only for wilcox, DESEQ
  if(type=="wilcox" || type == "deseq2")
  {
    if(type=="deseq2")
    {
      names(res.short)=c("pval","log2fc","sign","sgRNAs")
    }
    else
    {
      names(res.short)=c("pval","log2fc","sign")
    }
    
    #colorplot = c(rgb(1, 0, 0, alpha = 0.65),rgb(0, 0, 0, alpha = 0.65),rgb(0, 0, 1, alpha = 0.65))
    colors = c(rgb(217,35,35, 255, maxColorValue=255),rgb(46,98,166, 255, maxColorValue=255),"black")
    status = c("enriched","depleted",paste("pval threshold:", plot.p, sep="\n"))
    
    # color direction of genes
    # add colors to dataframe
    res.short$color = apply(res.short,1, function(x)
      if(is.na(x["sign"]))
      {
        #return white as color
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
      else if(as.numeric(x["sign"])<0 && as.numeric(x["pval"]) <= plot.p)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(46,98,166, 255, maxColorValue=255))
      }
      else if(as.numeric(x["sign"])>0 && as.numeric(x["pval"]) <= plot.p)
      {
        # majority of sgRNAs is enriched -> color in red and use 
        return(rgb(217,35,35, 255, maxColorValue=255))
      }
      else if(as.numeric(x["sign"])==0)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
      else
      {
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
    )
    
    
    # remove 0 in pval to circumvent problems with log
    res.short$pval = res.short$pval+0.0000001
    
    # create radius for circle plotting according to the p-value
    radius=(sqrt(res.short$log2fc)/2/pi)
    res.short$radius = radius
    
    #if(identical(separate, TRUE)) {
      old.par=par
      par(mfrow=c(2,1)) 
      par(mar=c(4,4,2,0.5))
      
      
      # split datasets into positive and negative ones
      data.split = split(res.short, res.short$sign >= 0)
      data.enriched = data.split[[2]]
      data.depleted = data.split[[1]]
      
      # for plot normalization, copy min and max of all datasets to each and make it invisible    
        ylim = if(-log10(plot.p) > max(-log10(data.enriched$pval),na.rm=TRUE)) 
        { c(0.1,-log10(plot.p)*offsetplot)}
        else
        {c(0.1,max(-log10(data.enriched$pval),na.rm=TRUE)*offsetplot)}
    
      #c(0,max(data.enriched$AUC, na.rm=TRUE)*1.1)
        xlim = c(1,nrow(data.enriched))
  
      
      data.max = as.vector(res.short[res.short$pval == max(res.short$pval, na.rm=TRUE),])
      
      data.max$color = rgb(1, 1, 1, alpha = 0)
  
      data.min = res.short[res.short$pval == min(res.short$pval, na.rm=TRUE),]
      data.min$color = rgb(1, 1, 1, alpha = 0)
      rownames(data.max)= NULL
      rownames(data.min)= NULL
  
      data.enriched = rbind(data.enriched, data.max)
      data.enriched = rbind(data.enriched, data.min)
  
      data.depleted = rbind(data.depleted, data.max)
      data.depleted = rbind(data.depleted, data.min)
      
      
      # plot upper plot
      symbols(seq(1,nrow(data.enriched)),(-log10(data.enriched$pval))*data.enriched$sign,circles=data.enriched$radius,inches=inches,
              fg=data.enriched$color,bg=data.enriched$color,xaxt="n",xlab="Genes",ylab="-log10 p-value", main=title, xlim=xlim, ylim=ylim)
      
      
      lines(seq(1,nrow(data.enriched)), rep(-log10(plot.p), nrow(data.enriched)), col="black", lty=2)
      
      legend("topright",status,,cex=0.8, bty="n", text.col=colors)
      
      if(print.names && cutoff[1] > 0){
        tophits=order(data.enriched[data.enriched$color !="#FFFFFF00","pval"],decreasing = FALSE)[1:cutoff[1]]
        text(seq(1,nrow(data.enriched))[tophits],
             (data.enriched$sign*(-log10(data.enriched$pval)))[tophits],
             row.names(data.enriched)[tophits]
        )
      }
      
      # plot negative plot
      ylim = 
        if(log10(plot.p) < min(log10(data.depleted$pval), na.rm=TRUE))
        {
          c(log10(plot.p)*offsetplot, -0.1)
        }
        else
        {c(min(log10(data.depleted[data.depleted$color!="#FFFFFF00","pval"]), na.rm=TRUE)*offsetplot, -0.1)}
      xlim = c(1,nrow(data.depleted))
      symbols(
        seq(1,nrow(data.depleted)),
        data.depleted$sign*(-log10(data.depleted$pval)),
        circles=data.depleted$radius,
        inches=inches,
        fg=data.depleted$color,
        bg=data.depleted$color,
        xaxt="n",
        xlab="Genes",
        ylab="-log10 p-value",
        main=title, 
        xlim=xlim, 
        ylim=ylim
      )
      # print p-value
      lines(seq(1,nrow(data.depleted)), rep(log10(plot.p), nrow(data.depleted)), col="black", lty=2)
      
      if(print.names && cutoff[3] > 0){
        tophits=order(data.depleted[data.depleted$color !="#FFFFFF00", "pval"],decreasing = FALSE)[1:cutoff[3]]
        text(seq(1,nrow(data.depleted))[tophits],
             (data.depleted$sign*(-log10(data.depleted$pval)))[tophits],
             row.names(data.depleted)[tophits]
        )
      }
      # set plotting par to default
      par=old.par
      
#     }
#     else
#     {
#       #if(is.null(ylim))
#       #{
#         c(-max(-log10(res.short[res.short[,"sign"]==-1,1]))*offsetplot, max(-log10(res.short$pval))*offsetplot)
#       #}
#       symbols(
#         seq(1,nrow(res.short)),
#         res.short$sign*(-log10((res.short$pval))),
#         circles=c(radius),
#         inches=inches,
#         fg=res.short$color,
#         bg=res.short$color,
#         xaxt="n",
#         xlab="Genes",
#         ylab="-log10 p-value",
#         main=title, 
#         xlim=xlim,
#         ylim=ylim
#       )
#       
#       legend("topright",status,,cex=0.8, bty="n", text.col=colors)
#       
#       # plot p-value estimation
#       lines(seq(1,nrow(res.short)), rep(-log10(plot.p), nrow(res.short)), col="black", lty=2)
#       lines(seq(1,nrow(res.short)), rep(log10(plot.p), nrow(res.short)), col="black", lty=2)
#       
#       
#       if(print.names){
#         # top
#         tophitstop=order(res.short[res.short$sign==1,"pval"],decreasing = FALSE)[1:cutoff[1]]
#         text(seq(1,nrow(res.short))[tophitstop],
#              (res.short$sign*(-log(res.short$pval)))[tophitstop],
#              row.names(res.short)[tophitstop]
#         )
#         
#         # bottom
#         tophits=order(res.short[res.short$sign==-1,"pval"],decreasing = FALSE)[1:cutoff[3]]
#         text(seq(1,nrow(res.short))[tophits],
#              (res.short$sign*(-log(res.short$pval)))[tophits],
#              row.names(res.short)[tophits]
#         )
#       }
#       
#     }
    
#     # FOR DESEq2 go for #sgRNA based plotting
#     if(type=="deseq2")
#     {
#       
#       # ########## Now go for sgRNA based plotting
#       
#       
#       # data from stat.deseq2
#       # 8 columns: baseMean, log2FoldCHange, lfcSE, stat, pvalue, padj, genes, sgRNAs
#       colors = c(rgb(217,35,35, 255, maxColorValue=255),rgb(46,98,166, 255, maxColorValue=255))
#       status = c("enriched","depleted")
#       
#       old.par=par
#       par(mfrow=c(2,1)) 
#       #par(mar=c(3,4,2,0.5))
#       par(mar=c(4,4,2,0.5))
#       
#       data.enriched = res.short[res.short$sign > 0,]
#       data.enriched = data.enriched[order(data.enriched$pval, decreasing = FALSE, na.last = TRUE),]
#       
#       #data.enriched = data.enriched[1:sgRNA.top,]
#       #order by name
#       #data.enriched = data.enriched[order(data.enriched$genes, decreasing = FALSE, na.last = TRUE),]
#       data.depleted = res.short[res.short$sign < 0,]
#       data.depleted = data.depleted[order(data.depleted$pval, decreasing = FALSE, na.last = TRUE),]
#       
#       #data.depleted = data.depleted[1:sgRNA.top,]
#       #order by name
#       #data.depleted= data.depleted[order(data.depleted$genes, decreasing = FALSE, na.last = TRUE),]
#       
#       
#       # create radius for circle plotting according to the p-value
#       data.enriched$radius=asinh(abs(data.enriched$log2fc))
#       data.depleted$radius=asinh(abs(data.depleted$log2fc))
#       
#       data.enriched$color = apply(data.enriched,1, function(x)
#         if(as.numeric(x["pval"]) <= plot.p)
#         {
#           # majority of sgRNAs is depleted -> color in blue and use 
#           return(rgb(217,35,35, 255, maxColorValue=255))
#         }
#         else
#         {
#           return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
#         }
#       )
#       
#       data.depleted$color = apply(data.depleted,1, function(x)
#         if(as.numeric(x["pval"]) <= plot.p)
#         {
#           # majority of sgRNAs is depleted -> color in blue and use 
#           return(rgb(46,98,166, 255, maxColorValue=255))
#         }
#         else
#         {
#           return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
#         }
#       )
#       
#       xlim = c(0,-log10(min(data.enriched$pval))*offsetplot)
#       
#       ylim = c(0, max(as.numeric(data.enriched$sgRNAs), na.rm=TRUE)*offsetplot)
#       
#       
#       # new plot 
#       # Plot sgRNA on Y-axis and pval (adjusted) on x-axis
#       # Ball size = change
#       
#       symbols(-log10(data.enriched$pval),data.enriched$sgRNAs,circles=data.enriched$radius,inches=inches,fg=data.enriched$color,bg=data.enriched$color, xlab="Genes", ylab="# significant sgRNAs", main=title, xlim=xlim, ylim=ylim, yaxt="n")
#       axis(2, at = seq(1, max(data.enriched$sgRNAs), by = 1), las=2)
#       
#       
#       legend("topright",status,,cex=0.8, bty="n", text.col=colors)
#       #legend("topleft",c("1","2","3","4"),, pch=16, cex=c(1,2,3,4), bty="n", text.col=colors)
#       
#       # create labels
#       if(print.names){
#         # pval
#         #tophitstop.pval=data.enriched[order(data.enriched[,"pval"],decreasing = FALSE, na.last=TRUE),]
#         #tophitstop.pval = tophitstop.pval[1:cutoff[1],]
#         
#         #text(-log10(tophitstop.pval$pval),
#              
#         #     tophitstop.pval$sgRNAs,
#        #      row.names(tophitstop.pval))
#         
#         #sgRNA
#         tophitstop.sgrna=data.enriched[order(data.enriched[,"sgRNAs"],decreasing = TRUE, na.last=TRUE),]
#         tophitstop.sgrna = tophitstop.sgrna[1:sgRNA.top,]
#         
#         text(-log10(tophitstop.sgrna$pval),
#              
#              tophitstop.sgrna$sgRNAs,
#              row.names(tophitstop.sgrna))
#       } 
#       
#       
#       
#       
#       # plot negative plot
#       
#       ylim = c(max(data.depleted$sgRNAs, na.rm=TRUE)*-offsetplot,0)
#       #xlim = c(1,nrow(data.depleted))
#       xlim = c(0,-log10(min(data.depleted$pval)))
#       
#       # plot normalization, see enriched
#       
#       symbols(-log10(data.depleted$pval),data.depleted$sgRNAs*(-1),circles=data.depleted$radius,inches=inches,fg=data.depleted$color,bg=data.depleted$color, xlab="Genes", ylab="# significant sgRNAs", main=title, xlim=xlim, ylim=ylim, yaxt="n")
#       axis(2, at = seq(-1, -max(data.depleted$sgRNAs), by = -1), las=2)
#       
#       # create labels
#       if(print.names){
#         # pval
#         tophitstop.pval=data.depleted[order(data.depleted[,"pval"],decreasing = FALSE, na.last=TRUE),]
#         #tophitstop.pval = tophitstop.pval[1:cutoff[1],]
#         
#         #text(-log10(tophitstop.pval$pval),
#              
#         #     tophitstop.pval$sgRNAs*(-1),
#         #     row.names(tophitstop.pval))
#         
#         #sgRNA
#         tophitstop.sgrna=data.depleted[order(data.depleted[,"sgRNAs"],decreasing = TRUE, na.last=TRUE),]
#         tophitstop.sgrna = tophitstop.sgrna[1:sgRNA.top,]
#         
#         text(-log10(tophitstop.sgrna$pval),
#              
#              tophitstop.sgrna$sgRNAs*(-1),
#              row.names(tophitstop.sgrna))
#         
#       }
#       
#       
#       
#       # set par back to normal
#       
#       par(mfrow=c(1,1))
#       par=old.par
#       
#       
#     }
#     
    
    return(res.short)
    
  } # end of plotting for p-value based output
  
  
  
  # plotting only for MAGECK pvalue based
  if(type=="mageck")
  {
    # input
    data = as.data.frame(data$genes)
    
    #genes      neg           rank.neg      pos           rank.pos sgrna.neg sgrna.pos
    #CASP8    0.000707        1             0.999891      258        18        82
    
    #colorplot = c(rgb(1, 0, 0, alpha = 0.65),rgb(0, 0, 0, alpha = 0.65),rgb(0, 0, 1, alpha = 0.65))
    colors = c(rgb(217,35,35, 255, maxColorValue=255),rgb(46,98,166, 255, maxColorValue=255),"black")
    status = c("enriched","depleted",paste("threshold:", plot.p, sep="\n"))
    
    #Direction of genes is not known, can only be retrieved by the list
    

    
    # remove 0 in pval to circumvent problems with log
    data$neg = data$neg+0.0000001
    data$pos = data$pos+0.0000001

      old.par=par
      par(mfrow=c(2,1)) 
      par(mar=c(4,4,2,0.5))
      
    #print(data[1:10,])  
    
      # split datasets into positive and negative ones
      data.enriched = data[order(data$rank.pos, decreasing = FALSE, na.last = TRUE),]
      #data.enriched = data.enriched[1:sgRNA.top,]
      #order by name
      data.enriched = data.enriched[order(data.enriched$genes, decreasing = FALSE, na.last = TRUE),]
   
      data.depleted = data[order(data$rank.neg, decreasing = FALSE, na.last = TRUE),]
      #data.depleted = data.depleted[1:sgRNA.top,]
      #order by name
      data.depleted= data.depleted[order(data.depleted$genes, decreasing = FALSE, na.last = TRUE),]
    
      # set radius according to # of sgRNAs
      data.enriched$radius = sqrt(data.enriched$sgrna.pos.good)/2/pi
      data.depleted$radius = sqrt(data.depleted$sgrna.neg.good)/2/pi
    
    
      data.enriched$color = apply(data.enriched,1, function(x)
      if(as.numeric(x["pos"]) <= plot.p)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(217,35,35, 255, maxColorValue=255))
      }
      else
      {
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
      )
    
    data.depleted$color = apply(data.depleted,1, function(x)
      if(as.numeric(x["neg"]) <= plot.p)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(46,98,166, 255, maxColorValue=255))
      }
      else
      {
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
      )
      
      # for plot normalization, copy min and max of all datasets to each and make it invisible    
      ylim = if(-log10(plot.p) > max(-log10(data.enriched$pos),na.rm=TRUE)) 
      { c(0.1,-log10(plot.p)*offsetplot)}
      else
      {c(0.1,max(-log10(data.enriched$pos),na.rm=TRUE)*offsetplot)}
      
      #c(0,max(data.enriched$AUC, na.rm=TRUE)*1.1)
      xlim = c(1,nrow(data.enriched))

      #print(data.enriched[1:10,])
    
      # plot upper plot
      symbols(seq(1,nrow(data.enriched)),(-log10(data.enriched$pos)),inches=inches, circles= data.enriched$radius,
              fg=data.enriched$color,bg=data.enriched$color,xaxt="n",xlab="Genes",ylab="-log10 p-value" ,main=title, xlim=xlim, ylim=ylim)
      
      
      lines(seq(1,nrow(data.enriched)), rep(-log10(plot.p), nrow(data.enriched)), col="black", lty=2)
      
      legend("topright",status,,cex=0.8, bty="n", text.col=colors)
      

      if(print.names){
        tophits=order(data.enriched[,"rank.pos"],decreasing = FALSE)[1:cutoff[1]]
        text(seq(1,nrow(data.enriched))[tophits],
             (-log10(data.enriched$pos))[tophits],
             row.names(data.enriched)[tophits]
        )
      }  
      # plot negative plot
      ylim = 
        if(log10(plot.p) < min(log10(data.depleted$neg), na.rm=TRUE))
        {
          c(log10(plot.p)*offsetplot, -0.1)
        }
      else
      {c(min(log10(data.depleted[,"neg"]), na.rm=TRUE)*offsetplot, -0.1)}
      xlim = c(1,nrow(data.depleted))

      symbols(
        seq(1,nrow(data.depleted)),
        (log10(data.depleted$neg)),
        inches=inches,
        circles = data.depleted$radius,
        fg=data.depleted$color,
        bg=data.depleted$color,
        xaxt="n",
        xlab="Genes",
        ylab="-log10 p-value",
        main=title, 
        xlim=xlim, 
        ylim=ylim
      )
      # print p-value
      lines(seq(1,nrow(data.depleted)), rep(log10(plot.p), nrow(data.depleted)), col="black", lty=2)
      
      if(print.names){
        tophits=order(data.depleted[,"rank.neg"],decreasing = FALSE)[1:cutoff[3]]
        text(seq(1,nrow(data.depleted))[tophits],
             (log10(data.depleted$neg))[tophits],
             row.names(data.depleted)[tophits]
        )
      } 
      # set plotting par to default
      par=old.par
      
    
    
    # ########## Now go for sgRNA based plotting
    # get cutoff
    #cutoffenriched.sgRNA = cutoff[2]
    #cutoffdepleted.sgRNA = cutoff[4]
    
    
    colors = c(rgb(217,35,35, 255, maxColorValue=255),rgb(46,98,166, 255, maxColorValue=255))
    status = c("enriched","depleted")
    
      old.par=par
      par(mfrow=c(2,1)) 
      par(mar=c(4,4,2,0.5))
      
 
    data.enriched = data[order(data$pos, decreasing = FALSE, na.last = TRUE),]
    
    #data.enriched = data.enriched[1:sgRNA.top,]
    #order by name
    data.enriched = data.enriched[order(data.enriched$genes, decreasing = FALSE, na.last = TRUE),]
    
    data.depleted = data[order(data$neg, decreasing = FALSE, na.last = TRUE),]
    #data.depleted = data.depleted[1:sgRNA.top,]
    #order by name
    data.depleted= data.depleted[order(data.depleted$genes, decreasing = FALSE, na.last = TRUE),]
      
    # create radius for circle plotting according to the p-value
    data.enriched$radius=asinh(data.enriched$sgrna.pos.good)
    data.depleted$radius=asinh(data.depleted$sgrna.neg.good)
    
    data.enriched$color = apply(data.enriched,1, function(x)
      if(as.numeric(x["pos"]) <= plot.p)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(217,35,35, 255, maxColorValue=255))
      }
      else
      {
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
    )
    
    data.depleted$color = apply(data.depleted,1, function(x)
      if(as.numeric(x["neg"]) <= plot.p)
      {
        # majority of sgRNAs is depleted -> color in blue and use 
        return(rgb(46,98,166, 255, maxColorValue=255))
      }
      else
      {
        return(rgb(211, 211, 211, alpha = 255, maxColorValue=255))
      }
    )
    
        xlim = c(0,-log10(min(data.enriched$pos))*offsetplot)
      
        ylim = c(0, max(as.numeric(data.enriched$sgrna.pos.good), na.rm=TRUE)*offsetplot)


        # new plot 
        # Plot sgRNA on Y-axis and asinh(AUC) on x-axis
        # Ball size = change
        
        symbols(-log10(data.enriched$pos),data.enriched$sgrna.pos.good,circles=data.enriched$radius,inches=inches,fg=data.enriched$color,bg=data.enriched$color, xlab="Genes", ylab="# significant sgRNAs", main=title, xlim=xlim, ylim=ylim)
        # print p-value
        abline(v=abs(log10(plot.p)), col="black", lty=2)
          
        legend("topright",status,,cex=0.8, bty="n", text.col=colors)
        #legend("topleft",c("1","2","3","4"),, pch=16, cex=c(1,2,3,4), bty="n", text.col=colors)
        
        # create labels
        if(print.names){
          # pval
          #tophitstop.pval=data.enriched[order(data.enriched[,"rank.pos"],decreasing = FALSE, na.last=TRUE),]
          #tophitstop.pval = tophitstop.pval[1:sgRNA.top,]
          
         # text(-log10(tophitstop.pval$pos),
          
          #     tophitstop.pval$sgrna.pos,
          #     row.names(tophitstop.pval))
          
          #sgRNA
          tophitstop.sgrna=data.enriched[order(data.enriched[,"sgrna.neg.good"],decreasing = TRUE, na.last=TRUE),]
          tophitstop.sgrna = tophitstop.sgrna[1:sgRNA.top,]
          
          text(-log10(tophitstop.sgrna$pos),
               
               tophitstop.sgrna$sgrna.pos.good,
               row.names(tophitstop.sgrna))
        } 
        


        
        # plot negative plot

        ylim = c(max(data.depleted$sgrna.neg.good, na.rm=TRUE)*-offsetplot,0)
        #xlim = c(1,nrow(data.depleted))
        xlim = c(0,-log10(min(data.depleted$neg)))
        
        # plot normalization, see enriched

        symbols(-log10(data.depleted$neg),data.depleted$sgrna.neg.good*(-1),circles=data.depleted$radius,inches=inches,fg=data.depleted$color,bg=data.depleted$color, xlab="Genes", ylab="# significant sgRNAs", main=title, xlim=xlim, ylim=ylim)
        abline(v=abs(log10(plot.p)), col="black", lty=2)
        
        # create labels
        if(print.names){
          # pval
          #tophitstop.pval=data.depleted[order(data.depleted[,"rank.neg"],decreasing = FALSE, na.last=TRUE),]
          #tophitstop.pval = tophitstop.pval[1:cutoff[1],]
          
          #text(-log10(tophitstop.pval$neg),
               
          #     tophitstop.pval$sgrna.neg.good*(-1),
          #     row.names(tophitstop.pval))
          
          #sgRNA
          tophitstop.sgrna=data.depleted[order(data.depleted[,"sgrna.neg.good"],decreasing = TRUE, na.last=TRUE),]
          tophitstop.sgrna = tophitstop.sgrna[1:sgRNA.top,]
          
          text(-log10(tophitstop.sgrna$neg),
               
               tophitstop.sgrna$sgrna.neg.good*(-1),
               row.names(tophitstop.sgrna))
          
        }
        

      
      # set par back to normal
      
      par(mfrow=c(1,1))
      par=old.par
      
    
    return(data)
    

   
    
  } # end of plotting for  MAGECK
  
  

  
  
}