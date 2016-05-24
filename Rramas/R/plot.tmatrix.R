plot.tmatrix <-
function(x, ...){
   #require(diagram)
   main <- deparse(substitute(x))
   di <- dim(x)[1]
   cosa <- summary(x)
   m.names<- cosa$m.names
   col<- rainbow(di)
   
   # set graphical parameters for the diagram
   self.arrpos <- NULL
   self.cex<-0.5
   
   if(di>6){ pos <-NULL 
      
      #self.arrpos <- NULL
      self.shiftx <- NULL
      self.shifty <- NULL
    }
   
    if (di==3){
       pos<-matrix(c( 0.25, 
               0.75,  0.50,
               0.25,
               0.25,  0.85),3,2)
       self.arrpos=c(pi,3*pi,pi)/2
       self.shiftx <- c(-0.08,.1,.13)
       self.shifty <- c(.1,-0.12,0.05)
    }
    if (di==4){ 
      pos <- matrix(c(0.15, 0.50, 0.85, 0.50,
                      0.5, 0.15, 0.5, 0.85)
                    ,4,2)
       self.arrpos=c(pi,3*pi,3*pi,pi,pi)/2
       self.shiftx <- c(0,   .13,  .075, .13)
       self.shifty <- c(0.13,-0.05,.12,  0.05)
    }
    if(di==5){
      pos<-matrix(c(0.15, 0.25, 
                    0.75, 0.85, 0.50,
                    0.60, 0.25,
                    0.25, 0.60, 0.85),5,2)
      self.arrpos=c(pi,3*pi,3*pi,pi,pi)/2
      self.shiftx <- c(0,0,.1,.075,.13)
      self.shifty <- c(0.13,-0.15,-0.12,.12,0.05)
    }
    if(di==6){ 
    pos<-matrix(c(0.15, 0.15, 0.50,
                    0.85, 0.85, 0.50,
                    0.65, 0.35, 0.15,
                    0.35, 0.65, 0.85),6,2)
    self.arrpos=c(pi,3*pi,3*pi,3*pi,pi,pi)/2
    self.shiftx <- c(0,0,-.13,.1,.075,.13)
    self.shifty <- c(0.13,-0.13,-0.05,-0.12,.12,0.05)
    }
    # end graphical parameters for the diagram
   
   old.par <- par(no.readonly = TRUE)
   on.exit(par(old.par)) 
   dev.new()
   plotmat(x, pos =pos, self.cex=self.cex,
            self.shiftx=self.shiftx,self.shifty=self.shifty,
             box.col=col, main=main, name=m.names)
   
   dev.new(width=7,height=3.5);par(mfrow=c(1,2))
   barplot(cosa$stable.stage.distribution, names.arg=m.names, col=col,
           main=paste(main, "\nstable age d."),las=2)
   barplot(cosa$reproductive.value, names.arg=m.names, col=col,
            main=paste(main, "\nreproductive value"),las=2)
   }

