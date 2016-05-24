#' @export
"phxplot" <-
  function(phxd=NULL,plot.type,cat.cov,cont.cov,
           forest.ci=c(0.025,0.5,0.975),multip=TRUE,outpdf=TRUE,
           scale=NULL,sel.ID,sparname){
    
    #########################################################################
    # Read in results
    #########################################################################
    
    phxdflag = 1
    
    if(is.null(phxd)){
      phxd = phxdata()
      phxdflag = NULL
    }
    
    #########################################################################
    # Function for multiplot
    #########################################################################
    # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
    # - cols:   Number of columns in layout
    # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
    #
    # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
    # then plot 1 will go in the upper left, 2 will go in the upper right, and
    # 3 will go all the way across the bottom.
    #
    multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
        
       # If layout is NULL, then use 'cols' to determine layout
       if (is.null(layout)) {
         # Make the panel
         # ncol: Number of columns of plots
         # nrow: Number of rows needed, calculated from # of cols
         layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                          ncol = cols, nrow = ceiling(numPlots/cols))
       }
        
       if (numPlots==1) {
         print(plots[[1]])
         
       } else {
         # Set up the page
         grid.newpage()
         pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
          
         # Make each plot, in the correct location
         for (i in 1:numPlots) {
           # Get the i,j matrix positions of the regions that contain this subplot
           matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
           print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
         }
       }
     }
    
    ###############################################################################
    
    ## Plotting starts
    
    ## Parameters correlation scatterplot matrix 
    if(plot.type=="correlation"){
      
      panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
      {
        usr <- par("usr"); on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- abs(cor(x, y))
        txt <- format(c(r, 0.123456789), digits=digits)[1]
        txt <- paste(prefix, txt, sep="")
        if(missing(cex.cor)) cex.cor <- 0.9/strwidth(txt)
        text(0.5, 0.5, txt, cex = cex.cor * r)
      }
      
      if(outpdf){
        pdf("correlation.plot.pdf")
        pairs(phxd$etatable[1:ncol(phxd$etatable)],
              lower.panel=panel.smooth, upper.panel=panel.cor, 
              pch=20)
        dev.off()
      }else{
        pairs(phxd$etatable[1:ncol(phxd$etatable)],
              lower.panel=panel.smooth, upper.panel=panel.cor, 
              pch=20)
      }
    }
        
    ## Obs vs pred plots [which(out.txt$PRED>0),]
    if(plot.type=="obs.pred"){
      #IPRED vs OBS   
      p1 = ggplot(data=phxd$out.txt, aes_string(x="IPRED",y="DV")) + 
        geom_point(shape=1,col="blue") +
        geom_abline(intercept=0,slope=1,col="red") 
      
      # PRED vs OBS 
      p2 = ggplot(data=phxd$out.txt, aes_string(x="PRED",y="DV")) + 
        geom_point(shape=1,col="blue") +
        geom_abline(intercept=0,slope=1,col="red") 
      
      if(outpdf){
        pdf("obs.pred.plot.pdf")
        
        if(is.null(scale)){
          print(p1 + geom_smooth(method=loess))
        }  else{
          print(p1 + scale_y_log10() + scale_x_log10())
        }
        
        if(is.null(scale)){
          print(p2 + geom_smooth(method=loess))
        }  else{
          print(p2 + scale_y_log10() + scale_x_log10())
        }        
        dev.off()    
      }else{
        if(is.null(scale)){
          print(p1 + geom_smooth(method=loess))
        }  else{
          print(p1 + scale_y_log10() + scale_x_log10())
        }
        
        if(is.null(scale)){
          print(p2 + geom_smooth(method=loess))
        }  else{
          print(p2 + scale_y_log10() + scale_x_log10())
        }        
      }
    }   
    
    ## Residual plots
    if(plot.type=="residual.scatter"){
      # WRES vs PRED 
      p3 = ggplot(data=phxd$out.txt, aes_string(x="PRED",y="WRES")) +
        geom_point(shape=1,col="blue") +
        geom_hline(yintercept=0,color="red")+
        geom_hline(yintercept=c(-2,2),linetype=2)
      
      # CWRES vs PRED 
      p4 = ggplot(data=phxd$out.txt, aes_string(x="PRED",y="CWRES")) +
        geom_point(shape=1,col="blue") +
        geom_hline(yintercept=0,color="red")+
        geom_hline(yintercept=c(-2,2),linetype=2)
      
      # CWRES vs TAD 
      p5 = ggplot(data=phxd$out.txt, aes_string(x="TAD",y="CWRES")) +
        geom_point(shape=1,col="blue") +
        geom_hline(yintercept=0,color="red")+
        geom_hline(yintercept=c(-2,2),linetype=2)
      
      # WRES vs TAD
      p6 = ggplot(data=phxd$out.txt, aes_string(x="TAD",y="WRES")) +
        geom_point(shape=1,col="blue") +
        geom_hline(yintercept=0,color="red")+
        geom_hline(yintercept=c(-2,2),linetype=2)
      
      if(outpdf){
        pdf("residual.plot.pdf")
        print (p3)
        print (p4)
        print (p5)
        print (p6)
        dev.off()    
      }else{
        print (p3)
        print (p4)
        print (p5)
        print (p6)
      }  
    }
    
    ## Parameters vs covariate (cat)
    if(plot.type=="param.catcov"){
      
      # parameters list
      parname = names(phxd$param[,-(1:3)]) #discard id, time and X..repl columns
      
      # categorical
      if(outpdf){
        pdf("param.catcov.plot.pdf")
        for (ii in 1:length(parname)){
          for(jj in 1:length(cat.cov)){
            phxd$phxnlme.data.first[,cat.cov[jj]] =  as.factor(phxd$phxnlme.data.first[,cat.cov[jj]])
            p7 = ggplot(data=phxd$phxnlme.data.first, aes_string(x=cat.cov[jj], y=parname[ii])) +
              geom_boxplot() +
              labs(x = cat.cov[jj], y = parname[ii]) 
            
            print (p7)
          }
        }
        dev.off()    
      }else{
        for (ii in 1:length(parname)) {
          for(jj in 1:length(cat.cov)){
            phxd$phxnlme.data.first[,cat.cov[jj]] =  as.factor(phxd$phxnlme.data.first[,cat.cov[jj]])
            p7 = ggplot(data=phxd$phxnlme.data.first, aes_string(x=cat.cov[jj], y=parname[ii])) +
              geom_boxplot() +
              labs(x = cat.cov[jj], y = parname[ii]) 
            
            print (p7)                
          }
        }
      }  
    }
    
    ## Parameters vs covariate (cont)
    if(plot.type=="param.contcov"){
      
      # parameters list
      parname = names(phxd$param[,-(1:3)])
      
      ## Parameters vs continous covariate
      if(outpdf){
        pdf("param.contcov.plot.pdf")
        for (ii in 1:length(parname)){
          for(jj in 1:length(cont.cov)){
            p8 = ggplot(data=phxd$phxnlme.data.first, aes_string(x=cont.cov[jj],y=parname[ii])) +
              geom_point(shape=1,col="blue") +
              geom_smooth(method="loess") +
              labs(x = cont.cov[jj], y = parname[ii]) 
            
            print(p8)
          }
        }
        dev.off()
      }else{
        for (ii in 1:length(parname)) {
          for(jj in 1:length(cont.cov)){
            p8 = ggplot(data=phxd$phxnlme.data.first, aes_string(x=cont.cov[jj],y=parname[ii])) +
              geom_point(shape=1,col="blue") +
              geom_smooth(method="loess") +
              labs(x = cont.cov[jj], y = parname[ii]) 
            
            print(p8)
          }        
        }
      }
    }
    
    ## Parameters histograms
    if(plot.type=="param"){
      
      # parameters list
      parname = names(phxd$param[,-(1:2)])
      
      p9 = list()
      
      for (ii in seq_along(parname)) {
        p9[[ii]]  = ggplot(data=phxd$phxnlme.data, aes_string(x=parname[ii])) +
          geom_histogram() +
          labs(x = parname[ii]) 
        
        # Add shrinkage #note that eta has to be named n followed by name of parameter
        r1 = suppressMessages(print(p9[[ii]]))
        
        for(jj in seq_along(phxd$dmp.txt$eta[,"shrinkage"])){
          rename=substring(names(phxd$dmp.txt$eta[,"shrinkage"])[jj],2)
          if(parname[ii]==rename){
            p9[[ii]] = p9[[ii]] + geom_text(y=max(r1$panel$ranges[[1]]$y.range)*0.8,
                                            x=max(phxd$phxnlme.data[,parname[ii]])*0.8,label=paste("Shrinkage=",phxd$dmp.txt$eta[,"shrinkage"][jj],sep=""))
          }
        }  
      }
      
      if(outpdf){
        pdf("param.histograms.plot.pdf")
        if(multip){
          suppressMessages(multiplot(plotlist=p9,cols=2))
          #suppressMessages(do.call(grid.arrange, p9))    
        }else{
          suppressMessages(print(p9))      
        }
        dev.off()    
      } else{
        if(multip){
          suppressMessages(multiplot(plotlist=p9,cols=2))
          #suppressMessages(do.call(grid.arrange, p9))    
        }else{
          suppressMessages(print(p9))      
        }
      }
      
    }  
    
    ## Forest plots
    if(plot.type=="forest"){
      
      # parameters list
      parname = sparname 
      
      # categorical
      cat.cov=cat.cov
      Group.1=NULL
      outsum=data.frame(Group.1=NULL,x=NULL,Par=NULL)
      for(ii in seq_along(parname)) {
        for(jj in seq_along(cat.cov)){
          tsum = aggregate(phxd$phxnlme.data[,parname[ii]],list(phxd$phxnlme.data[,cat.cov[jj]]),FUN = function(x) { c(q = quantile(x,probs=forest.ci)) })
          tsum[[1]] = paste(cat.cov[jj],tsum[[1]],sep="=")
          tsum$Par=parname[ii]
          outsum = rbind(outsum,tsum)
        }
      }
      ylo = NULL
      y = NULL
      yhi = NULL
      outsum$ylo = outsum[2][[1]][,1]
      outsum$y = outsum[2][[1]][,2]
      outsum$yhi = outsum[2][[1]][,3]
      outsum = merge(outsum,t(as.matrix(phxd$dmp.txt$coefficients$fixed)))  
      
      for(ii in unique(outsum$Par)){
        outsum[outsum$Par==ii,c("ylo","y","yhi")] = outsum[outsum$Par==ii,c("ylo","y","yhi")]/ unique((outsum[,paste("tv",ii,sep="")]))
      }
      
      p10 = ggplot(data=outsum, aes(x=Group.1,y=y)) +
        geom_point() +
        geom_errorbar(aes(ymin=ylo,ymax=yhi),width=0.2) +
        geom_hline(yintercept=1,col="blue") +
        geom_hline(yintercept=c(0.8,1.25),linetype=2) + 
        facet_wrap(~Par) #, scales = "free") 
      
      if(outpdf){
        pdf("forest.plot.pdf")
        print (p10 + coord_flip() + labs(title = "", x = "", y = "Relative to typical value"))         
        dev.off()    
      }else{
        print (p10 + coord_flip() + labs(title = "", x = "", y = "Relative to typical value"))         
      }        
    }
    
    ## QQ plots of etas
    if(plot.type=="qq"){
      
      if(outpdf){pdf("qqnorms.plot.pdf")}
      p11 = par(mfrow=c(1,1))
      for(ii in unique(names(phxd$etatable))){
        qqnorm(phxd$etatable[,ii],main=paste(ii, "Normal Q-Q Plot",sep=" "))
        qqline(phxd$etatable[,ii])   
      }
      if(outpdf){dev.off()}
    }
    
    ## Individual DV/IPRED/PRED vs IDV
    if(plot.type=="ind"){
      
      if(outpdf){pdf("ind.plot.pdf")}
      par(mfrow=c(3,3))
      rout.txt = phxd$out.txt[order(phxd$out.txt$ID5,phxd$out.txt$TAD),]
      
      ulog="y"
      if(is.null(scale)){
        ulog=""
      }
      
      for(ii in unique(rout.txt$ID5)){
        temp <- rout.txt[rout.txt$ID5==ii,]
        x.max=max(temp$TAD)*1.5
        plot(temp$TAD,temp$DV, ylim=c(min(temp$DV)*0.5,max(temp$DV)*1.5),xlim=c(0,x.max),
             xlab="TAD",ylab="DV",type="p",log=ulog) 
        lines(x=temp$TAD,y=temp$IPRED,  col='red',type="l") 
        lines(x=temp$TAD,y=temp$PRED, col='blue',type="l",lty=2) 
        title(paste("ID=",ii))
      }
      
      if(outpdf){dev.off()}  
      
    }
    
    ## Individual DV/IPRED/PRED vs IDV
    if(plot.type=="ind.dynamic"){
      
      rout.txt = phxd$out.txt[order(phxd$out.txt$ID5,phxd$out.txt$TAD),]
      temp <- rout.txt[rout.txt$ID5==sel.ID,]
      x.max=max(temp$TAD)*1.5
      manipulate(  
        plot(temp$TAD,temp$DV, ylim=c(min(temp$DV)*0.5,max(temp$DV)*1.5),xlim=c(0,x.max),
             xlab="TAD",ylab="DV",type="p") +
          lines(x=temp$TAD,y=temp$IPRED,  col='red',type="l") +
          lines(x=temp$TAD,y=temp$PRED, col='blue',type="l",lty=2),x.max=slider(0,x.max) )
      title(paste("ID=",sel.ID))
    }
    
    if(is.null(phxdflag)){
      ### Create results folder and copy files in
      path = getwd()
      if(file.exists("Results")){
        pplots = Sys.glob("*plot.pdf")
        res = paste(path,pplots,sep="/")
        file.copy(res,paste(path,"Results",sep="/"),overwrite=TRUE)
      }else{
        dir.create("Results")
        pplots = Sys.glob("*plot.pdf")
        res = paste(path,pplots,sep="/")
        file.copy(res,paste(path,"Results",sep="/"),overwrite=TRUE)
      }
    }
  }


