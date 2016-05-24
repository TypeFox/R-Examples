slider.stem.leaf<-function(x,main=main){
  main<-if(missing(main)) deparse(substitute(x)) else main
  slider(obj.name="main",
         obj.value=paste("stem and leaf of",main))
  slider(obj.name="rule.line",obj.value="Dixon")
  repl<-function(...){
    rule.line<-slider(obj.name="rule.line")
    main<-c(slider(obj.name="main"), paste("rule:", rule.line))
    cex<-slider(no=1); m<-slider(no=2)

    res<-if(m>0) stem.leaf(x,m=m,rule.line=rule.line) else 
                 stem.leaf(x,rule.line=rule.line)
    stem<-rev(res$stem); stem <- paste(rev(res$depth), stem, rev(res$leaves))

    # x.max<-max(nchar(res$stem));  y.max<-length(res$stem)
    x.max<-max(nchar(stem));  y.max<-length(stem) #140922

    opar<-par(mai=c(0.5,0.2,0.5,0),oma=c(0,0,0,0), omi=c(0,0,0,0))
    plot(0,type="n",xlim=c(1,x.max),
         ylim=c(-1.5,y.max+2),ann=FALSE,axes=FALSE)

    text(1,-2,paste(res$info,collapse=","),cex=cex,adj=0,xpd=NA)
    if(!is.null(res$upper)){
      text(1,     -0.5,paste(res$upp,collapse=","),cex=cex,adj=0,xpd=NA)
    }
    if(!is.null(res$lower)) 
      text(1,y.max+1.5,paste(res$low,collapse=","),cex=cex,adj=0,xpd=NA)
    title(main,cex.main=cex)
    for(i in 1:y.max){
      n<-nchar(h<-stem[i]); h<-substring(h,1:n,1:n)
      text(1:n,i,h,cex=cex)
    }
    par(opar)
  } # end of repl
  b1<-function(...){
    slider(obj.name="rule.line",obj.value="Dixon")
    repl(); print(slider(obj.name="rule.line"))
  }
  b2<-function(...){
    slider(obj.name="rule.line",obj.value="Velleman")
    repl(); print(slider(obj.name="rule.line"))
  }
  b3<-function(...){
    slider(obj.name="rule.line",obj.value="Sturges")
    repl(); print(slider(obj.name="rule.line"))
  }
  sp<-function(...){
    din<-par()$din
    filename<-tkgetSaveFile(filetypes="{{Paper Files} {.ps|.jpg}}",
                          title="name of ps-/jpg-file?", initialdir=getwd(),
                          defaultextension=".ps", initialfile="stem-and-leaf")
    filename<-sub("\\.ps$","",sub("\\.jpg$","",filename))
    filename<-paste(filename,".ps",sep="")  
    dev.copy(postscript,filename,width=din[1],height=din[2],
             horizontal=FALSE); dev.off()
    filename<-sub(".ps",".jpg",filename) 
    dev.copy(jpeg,filename,units="in",res=95,width=din[1],height=din[2])
    dev.off()
  }
  nt <- slider(repl, # 140925
         c("character size","rows per unit, 0=by rule"),
         c(0.5,0),c(5,10),c(.1,1),c(1,0),
         c(b1,b2,b3,sp),c("Dixon rule","Velleman rule","Sturges rule",
          "save plot")     
  )
  tkwm.minsize(nt, "450", "110") # set width, height to prevent to small sizes
  repl()
  "choose layout by control widget"
} # slider.stem.leaf(c(co2,100,1000))

