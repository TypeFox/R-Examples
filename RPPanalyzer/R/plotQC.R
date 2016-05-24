`plotQC` <-
function(x,file="target_vs_blank.pdf",arrays2rm=c("protein")){

         x <- remove.arrays(x,param="target",arrays2rm=arrays2rm)
         
         #inc.run <- levels(as.factor(x[[3]]["incubation_run",]))
         inc.run <- unique(x[[3]]["incubation_run",])

         ctrl.lines <- which ( x[[4]][,"sample_type"]=="control")
         
         controls <- unique(x[[4]][ctrl.lines,"sample"])

          pdf(file=file)
          
    for (j in seq(along=inc.run)){
    
         ## define columns of a single incubation run
          run.cols <- which((x[[3]]["incubation_run",])==inc.run[j])
          
          ## select data of a single incubation run
           temp.dat <- x[[1]][,run.cols]
            temp.arrays <- x[[3]] [,run.cols]
            blk.col <- grep("lank$",temp.arrays["target",])
            
     for (i in 1:ncol(temp.dat)){
     
          #ctrl.lines <- which ( x[[4]][,"sample_type"]=="control")
          temp.dat <- x[[1]][ctrl.lines,run.cols]
          temp.anno <- x[[4]][ctrl.lines,]
          
          par(lwd=2)
          plot(x=NULL,y=NULL
               , main=c(paste("target: ",temp.arrays["target",i],"  Ab: ",temp.arrays["AB_ID",i])
                        ,paste("incubation run: ",inc.run[j]))
               ,xlab="concentration",ylab="signal intensity [a.u.]"
               ,ylim=c(0,max(c(temp.dat[,i],temp.dat[,blk.col])))
               ,xlim=c(0,max(x[[4]][ctrl.lines,"concentration"])))

          for (k in seq(along=controls)){
          
             ctrl.lines.ii <- which( temp.anno[,"sample"]==controls[k])
             
             points(x=temp.anno[ctrl.lines.ii,"concentration"],y=temp.dat[ctrl.lines.ii,i]
                     ,pch=k+1,col=k+1)
          

          }
          
          points ( x[[4]][ctrl.lines,"concentration"],temp.dat[,blk.col],
                  col=1,pch=1)
          legend("topleft",col=c(1,(seq(along=controls)+1)),pch=c(1,(seq(along=controls)+1))
                  ,legend=c("blank",as.character(controls)))

  }
  }
  dev.off()
}

