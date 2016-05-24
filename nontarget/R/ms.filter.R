ms.filter <-
function(component,x="mz",y="dm",xlim=FALSE,ylim=FALSE,
                    rm.comp=FALSE,plot.comp=TRUE,rm.noncomp=TRUE,
                    select.polygon="inside",res=100,filter.for="raw"){

    ############################################################################
    # check inputs #############################################################
    if(x==y){stop("x should differ from y")}
    if(any(x==c("int","mz","rt","dm"))==FALSE){stop("x must be int,mz,rt or dm")}
    if(any(y==c("int","mz","rt","dm"))==FALSE){stop("x must be int,mz,rt or dm")}
    if(any(select.polygon!=c("inside","outside"))==FALSE){stop("select.polygon must be either inside or outside")}
    if(xlim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
    if(ylim[1]!=FALSE){if(length(xlim)>2){stop("xlim not correct!")}}
    if(rm.comp!=TRUE & rm.comp!=FALSE){stop("rm.components: TRUE or FALSE")}
    if(all(names(component)==c("Components","pattern peak list","adduct peak list","homologue list","Peaks in components","Summary","Parameters"))==FALSE)
    {stop("invalid component argument")}
    if(is.numeric(res)==FALSE){stop("res must be numeric")}
    if(res<10){stop("res must be >10")}
    if(plot.comp!=FALSE & plot.comp!=TRUE){stop("plot.comp must be either TRUE or FALSE")}
    if(filter.for!="raw" & filter.for!="pattern" & filter.for!="adduct"){stop("filter.for must be raw, pattern or adduct")}
    if(rm.noncomp==FALSE & rm.comp==FALSE){cat("both rm.comp and rm.noncomp = FALSE -> nothing will be selected!\n")}
    if(filter.for!="pattern" & length(component[[2]])<1){stop("Cannot filter.for pattern, since no pattern groups provided!")}
    if(filter.for!="adduct" & length(component[[2]])<1){stop("Cannot filter.for adduct, since no adduct groups provided!")}
    ############################################################################
    ############################################################################
    # generate inputs of x and y & plot ########################################
    # x ########################################################################
    if(length(component[[2]])>1){
      if(x=="int"){xplot<-component[[2]][,2]};
      if(x=="mz"){xplot<-component[[2]][,1]};
      if(x=="rt"){xplot<-component[[2]][,3]};
      if(x=="dm"){
        xplot<-c()
        for(i in 1:length(component[[2]][,1])){
            n<-as.numeric(paste("0.",strsplit(as.character((component[[2]][i,1])),".",fixed=TRUE)[[1]][2],sep=""))
            if(n>0.5){n<-c((1-n)*-1);}
            xplot<-c(xplot,n);
        };
      };
    }else{
    if(x=="int"){xplot<-component[[3]][,2]};
      if(x=="mz"){xplot<-component[[3]][,1]};
      if(x=="rt"){xplot<-component[[3]][,3]};
      if(x=="dm"){
        xplot<-c()
        for(i in 1:length(component[[3]][,1])){
            n<-as.numeric(paste("0.",strsplit(as.character((component[[3]][i,1])),".",fixed=TRUE)[[1]][2],sep=""))
            if(n>0.5){n<-c((1-n)*-1);}
            xplot<-c(xplot,n);
        };
      };    
    }
    # y ########################################################################
    if(length(component[[2]])>1){    
        if(y=="int"){yplot<-component[[2]][,2]};
        if(y=="mz"){yplot<-component[[2]][,1]};
        if(y=="rt"){yplot<-component[[2]][,3]};
        if(y=="dm"){
          yplot<-c()
          for(i in 1:length(component[[2]][,1])){
              n<-as.numeric(paste("0.",strsplit(as.character((component[[2]][i,1])),".",fixed=TRUE)[[1]][2],sep=""))
              if(n>0.5){n<-c((1-n)*-1);}
              yplot<-c(yplot,n);
          };
        };
    }else{
        if(y=="int"){yplot<-component[[3]][,2]};
        if(y=="mz"){yplot<-component[[3]][,1]};
        if(y=="rt"){yplot<-component[[3]][,3]};
        if(y=="dm"){
          yplot<-c()
          for(i in 1:length(component[[3]][,1])){
              n<-as.numeric(paste("0.",strsplit(as.character((component[[3]][i,1])),".",fixed=TRUE)[[1]][2],sep=""))
              if(n>0.5){n<-c((1-n)*-1);}
              yplot<-c(yplot,n);
          };
        };
    }
    # plot all data ############################################################
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
    plot.new()
    if(is.numeric(ylim[1])==FALSE){ylim=c(min(yplot),max(yplot))};
    if(is.numeric(xlim[1])==FALSE){xlim=c(min(xplot),max(xplot))};
    plot.window(xlim=xlim,ylim=ylim);
    box();axis(1);axis(2);
    if(y=="int"){title(ylab="Intensity")}
    if(y=="mz"){title(ylab="m/z")}
    if(y=="rt"){title(ylab="Retention time")}
    if(y=="dm"){title(ylab="Mass defect")}
    if(x=="int"){title(xlab="Intensity")}
    if(x=="mz"){title(xlab="m/z")}
    if(x=="rt"){title(xlab="Retention time")}
    if(x=="dm"){title(xlab="Mass defect")}
    points(xplot,yplot,col="darkgrey",bg="darkgrey",pch=22,cex=0.2)
    # add components ###########################################################
    if(length(component[[2]])>1){
      getit<-c(seq(1,length(component[[2]][,1]),1)[component[[5]]]);
    }else{
      getit<-c(seq(1,length(component[[3]][,1]),1)[component[[5]]]);
    }
    if(plot.comp==TRUE){
      points(xplot[component[[5]]],yplot[component[[5]]],col="red",bg="red",pch=22,cex=0.3)
    };
    ############################################################################
    ############################################################################
    title(main="Use left mouse key to draw polygon \n Use right mouse key to finish (stopp)")
    options(locatorBell=FALSE);
    cat("\nNow, select points in plot for polygon vertices. Select STOPP (R Graphic device toolbar) when finished.\n\n");
    these<-locator(type="o");
    lines(c(these$x[1],these$x[length(these$x)]),c(these$y[1],these$y[length(these$y)]));
    cat("Please wait...\n")
    x1<-seq(min(these$x),max(these$x),(max(these$x)-min(these$x))/res);
    dx1<-abs(x1[2]-x1[1]);
    y1<-seq(max(these$y),min(these$y),-(max(these$y)-min(these$y))/res);
    dy1<-abs(y1[2]-y1[1]);
    x3<-c(these$x[1]);
    y3<-c(these$y[1]);
    count<-c(1);
    for(i in 1:(length(these$x))){
      if(i!=length(these$x)){
           #points(these$x[i],these$y[i],col="red",pch=19,cex=1);
           #points(these$x[i+1],these$y[i+1],col="red",pch=19,cex=1);
           dx<-c(these$x[i+1]-these$x[i]);
           dy<-c(these$y[i+1]-these$y[i]);
           if(dx!=0 || dy!=0){
               itx<-abs(dx/(dx1/2));
               ity<-abs(dy/(dy1/2));
               it<-ceiling(max(itx,ity));
               dx<-c(dx/it);
               dy<-c(dy/it);
               for(j in 1:it){
                  x3<-c(x3,x3[count]+dx);
                  y3<-c(y3,y3[count]+dy);
                  #points(x3[count],y3[count],col="red",pch=19,cex=0.3);
                  count<-c(count+1);
               }
           }
           x3<-c(x3,these$x[i+1]);
           y3<-c(y3,these$y[i+1]);
           #points(x3[count],y3[count],col="red",pch=19,cex=0.3);
           count<-c(count+1);                      
      }else{
           #points(these$x[i],these$y[i],col="red",pch=19,cex=1);
           #points(these$x[1],these$y[1],col="red",pch=19,cex=1);
           dx<-c(these$x[1]-these$x[i]);
           dy<-c(these$y[1]-these$y[i]);
           if(dx!=0 || dy!=0){
                 itx<-c(abs(dx/(dx1/2)));
                 ity<-c(abs(dy/(dy1/2)));
                 it<-ceiling(max(itx,ity));
                 dx<-c(dx/it);
                 dy<-c(dy/it);
                 for(j in 1:it){
                    x3<-c(x3,x3[count]+dx);
                    y3<-c(y3,y3[count]+dy);
                    #points(x3[count],y3[count],col="red",pch=19,cex=0.3);
                    count<-c(count+1);
                 }
           }
           x3<-c(x3,these$x[1]);
           y3<-c(y3,these$y[1]);
           #points(x3[count],y3[count],col="red",pch=19,cex=0.3);
           count<-c(count+1);           
      }
    }
    # generate assignment matrix ###############################################
    mat<-matrix(nrow=length(y1),ncol=length(x1),0);
    dx2<-dx1/1.8;
    dy2<-abs(dy1/1.8);
    for(i in 1:length(x1)){
      for(j in 1:length(y1)){
       if(any( 
          x3>=(x1[i]-dx2) & 
          x3<=(x1[i]+dx2) &
          y3>=(y1[j]-dy2) & 
          y3<=(y1[j]+dy2)      
          )){mat[j,i]<-1;}
      }
    }
    # check if matrix is correct for given resolution: matrix tight? ###########
    # inside ...
    for(i in 2:(length(mat[,1])-1)){
      for(j in 2:(length(mat[1,])-1)){
        if(mat[i,j]==1){
            that<-c(  mat[i+1,j]+mat[i-1,j]+mat[i,j+1]+mat[i,j-1]+
                      mat[i-1,j+1]+mat[i-1,j-1]+mat[i+1,j+1]+mat[i+1,j-1])
            if(that<1){stop("increase res!")
            }
        }
      }
    }
    # ... and on the margins ... skipped
    # fill with 2 ##############################################################
    # rigth to left
    for(i in 1:length(mat[,1])){j<-c(1);while(mat[i,j]!=1 & j<=length(mat[,1])){mat[i,j]<-2;j<-c(j+1);}}
    # left to right
    for(i in 1:length(mat[,1])){j<-length(mat[,1]);while(mat[i,j]!=1 & j>0){mat[i,j]<-2;j<-c(j-1);}}
    # from top to bottom
    for(i in 1:length(mat[1,])){j<-c(1);while(mat[j,i]!=1 & j<=length(mat[1,])){mat[j,i]<-2;j<-c(j+1);}}
    # from bottom to top
    for(i in 1:length(mat[1,])){j<-length(mat[1,]);while(mat[j,i]!=1 & j>0){mat[j,i]<-2;j<-c(j-1);}}
    # fill shadowed regions
    itis<-TRUE;
    while(itis==TRUE){
      itis<-FALSE
      for(i in 2:(length(mat[,1])-1)){
        for(j in 2:(length(mat[1,])-1)){
          if(mat[i,j]==2){
            if(mat[(i-1),(j)]==0){mat[(i-1),(j)]<-2;itis<-TRUE;}
            if(mat[(i+1),(j)]==0){mat[(i+1),(j)]<-2;itis<-TRUE;}
            if(mat[(i),(j-1)]==0){mat[(i),(j-1)]<-2;itis<-TRUE;}
            if(mat[(i),(j+1)]==0){mat[(i),(j+1)]<-2;itis<-TRUE;}
    }}}}
    mat[mat==0]<-1;
    # assign values inside / outside that polygon ##############################
    if(length(component[[2]])>1){
      not<-rep(FALSE,length(component[[2]][,1]));
    }else{
      not<-rep(FALSE,length(component[[3]][,1]));
    }
    if(rm.noncomp==TRUE || rm.comp==TRUE){  # remove anthing at all ? ##########
    # intside? #################################################################
    if(select.polygon=="inside"){
    for(i in 1:length(mat[,1])){            # per matrix row
      for(j in 1:length(mat[1,])){          # per matrix column
      if(mat[i,j]==1){
        # remove components AND non-components #################################
        if(rm.noncomp==TRUE & rm.comp==TRUE){
              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) ] <- TRUE
        }
        # remove only non-components ###########################################
        if(rm.noncomp==TRUE & rm.comp==FALSE){
              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) ] <- TRUE
              not[getit]<-FALSE
        }
        # remove only components ###############################################
        if(rm.noncomp==FALSE & rm.comp==TRUE){

              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) &
                    component[[5]]==1                ] <-TRUE
        }
        ########################################################################      
    }}}
    } # inside?
    # outside? #################################################################
    if(select.polygon=="outside"){
    # mark anything outside that matrix area! ##################################
    # remove components AND non-components #####################################
    if(rm.noncomp==TRUE & rm.comp==TRUE){
        not[xplot<min(x3)]<-TRUE
        not[xplot>max(x3)]<-TRUE
        not[yplot<min(y3)]<-TRUE
        not[yplot>max(y3)]<-TRUE
    }
    # remove only non-components ###############################################
    if(rm.noncomp==TRUE & rm.comp==FALSE){
        not[xplot<min(x3)]<-TRUE
        not[xplot>max(x3)]<-TRUE
        not[yplot<min(y3)]<-TRUE
        not[yplot>max(y3)]<-TRUE
        not[getit]<-FALSE
    }
    # remove only components ###################################################
    if(rm.noncomp==FALSE & rm.comp==TRUE){
        not[xplot<min(x3) & component[[5]]==1]<-TRUE
        not[xplot>max(x3) & component[[5]]==1]<-TRUE
        not[yplot<min(y3) & component[[5]]==1]<-TRUE
        not[yplot>max(y3) & component[[5]]==1]<-TRUE
    }
    # then, mark anything intside that matrix area! ############################
    for(i in 1:length(mat[,1])){            # per matrix row
      for(j in 1:length(mat[1,])){          # per matrix column
      if(mat[i,j]==2){
        # remove components AND non-components #################################
        if(rm.noncomp==TRUE & rm.comp==TRUE){
              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) ] <- TRUE
        }
        # remove only non-components ###########################################
        if(rm.noncomp==TRUE & rm.comp==FALSE){
              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) ] <- TRUE
              not[getit]<-FALSE
        }
        # remove only components ###############################################
        if(rm.noncomp==FALSE & rm.comp==TRUE){

              not[  (xplot<=(x1[j]+(dx1/1.8))) &  
                    (xplot>=(x1[j]-(dx1/1.8))) &
                    (yplot<=(y1[i]+(dy1/1.8))) &
                    (yplot>=(y1[i]-(dy1/1.8))) &
                    component[[5]]==1                ] <-TRUE
        }
        ########################################################################      
    }}}
    } # outside?
    } # if: do that at all?
    ############################################################################
    # plot results! ############################################################
    plot.new()
    plot.window(xlim=xlim,ylim=ylim);    box();axis(1);axis(2);
    if(y=="int"){title(ylab="Intensity")}
    if(y=="mz"){title(ylab="m/z")}
    if(y=="rt"){title(ylab="Retention time")}
    if(y=="dm"){title(ylab="Mass defect")}
    if(x=="int"){title(xlab="Intensity")}
    if(x=="mz"){title(xlab="m/z")}
    if(x=="rt"){title(xlab="Retention time")}
    if(x=="dm"){title(xlab="Mass defect")}
    points(xplot[not==FALSE],yplot[not==FALSE],col="darkgrey",bg="darkgrey",pch=22,cex=0.2)
    if(plot.comp==TRUE){
      if(rm.comp==FALSE){
          points(xplot[getit],yplot[getit],col="red",bg="red",pch=22,cex=0.3);
      }else{
          getit<-getit[not[getit]==FALSE];
          points(xplot[getit],yplot[getit],col="red",bg="red",pch=22,cex=0.3);
      }
    }
    if(length(component[[2]])>1){
      cat(paste("Number of data points removed: ",length(not[not==TRUE])," of ",length(component[[2]][,1]),"\n",sep=""))
    }else{
      cat(paste("Number of data points removed: ",length(not[not==TRUE])," of ",length(component[[3]][,1]),"\n",sep=""))    
    }
    par<-def.par #- reset to default
    ############################################################################
    ############################################################################
    # select data to return ####################################################
    if(filter.for=="raw"){
      if(length(component[[2]])>1){
        peak.sel<-data.frame(component[[2]][not==FALSE,c(1,2,3)]);
      }else{
        peak.sel<-data.frame(component[[3]][not==FALSE,c(1,2,3)]);
      }
    }
    if(filter.for=="pattern"){peak.sel<-data.frame(component[[2]][not==FALSE,])};
    if(filter.for=="adduct"){peak.sel<-data.frame(component[[3]][not==FALSE,])};
    return(peak.sel);
    ############################################################################
}
