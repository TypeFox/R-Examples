spin3R <- function(x, alpha=1, delay=.015, na.rm=FALSE){
 #################################################################
 # spin3R: simple spin function to rotate a 3-dim cloud of points#
 # pwolf 070831                                              #
 #                                                               #
 # arguments:                                                    #
 #                                                               #
 #  x             (nx3)-matrix of points                         #
 #  alpha         arc of rotation                                #
 #  delay         sleeping time between rotations                #
 #                                                               #
 #################################################################
if(ncol(x)!=3) { print("Error: data matrix must have 3 columns"); return() }
 # require(tcltk) # 131104
 Rot <-tclVar("relax");bw <- 4
 topl<-tktoplevel();   tkwm.geometry(topl,"+100+100")
 f1 <- tkframe(topl);f2 <- tkframe(topl);f3 <- tkframe(topl)
 f4 <- tkframe(topl);f5 <- tkframe(topl);tkpack(f1,f2,f3,f4,f5)

 b12 <- tkbutton(f1, relief="ridge",  width=bw, text="up")
 b21 <- tkbutton(f2, relief="ridge",  width=bw, text="left")
 b22 <- tklabel(f2,  relief="flat",   width=bw)
 b23 <- tkbutton(f2, relief="ridge",  width=bw, text="right")
 b32 <- tkbutton(f3, relief="ridge",  width=bw, text="down")
 b41 <- tkbutton(f4, relief="ridge",  width=bw, text="clock")
 b42 <- tklabel(f4,  relief="flat",   width=bw)
 b43 <- tkbutton(f4, relief="ridge",  width=bw, text="cclock")
 b51 <- tkbutton(f5, relief="raised", width=bw, text="reset")
 b52 <- tklabel(f5,  relief="flat",   width=bw)
 b53 <- tkbutton(f5, relief="raised", width=bw, text="exit")
 tkpack(b12,b32)
 tkpack(b21,b22,b41,b42,b51,b52,side="left")
 tkpack(b23,b43,b53,side="right")

 alpha<-alpha/360*2*pi; ca<-cos(alpha); sa<-sin(alpha)
 rot<-matrix(c(ca,-sa,sa,ca),2,2)

 for(type in c("12","21","23","32","41","43")){
  b<-eval(parse(text=paste("b",type,sep="")))
  tkbind(b, "<Enter>",
      eval(parse(text=paste("function()tclvalue(Rot)<-\"",type,"\"",sep=""))))
  tkbind(b, "<Leave>",function() tclvalue(Rot) <- "relax")
 }
 tkconfigure(b51,command=function() tclvalue(Rot) <- "reset" )
 tkconfigure(b53,command=function() tclvalue(Rot) <- "exit" )

 n <- nrow(x)
 if(any(is.na(x))){
   if(na.rm){ x<-x[!apply(is.na(x),1,any),,drop=FALSE]
     print("Warning: NA elements have been removed!!")
   }else{
     xy.means<-colMeans(x,na.rm=TRUE)
     for(j in 1:ncol(x)) x[is.na(x[,j]),j]<-xy.means[j]
     print("Warning: NA elements have been exchanged by mean values!!")
   }  
 }
 x <- x - matrix(apply(x,2,min),n,3,TRUE)
 x.o<-x<-x / matrix(apply(x,2,max),n,3,TRUE) - 0.5;                xa <- x[,2:3]
 A.o<-A<-0.5*matrix(c(1,0,0, 0,0,0, 0,1,0, 0,0,0, 0,0,1),5,3,TRUE);Aa <- A[,2:3]
 plot(xa, xlim=.7*c(-1,1), ylim=.7*c(-1,1),
           pch=20, xlab="",ylab="",xaxt="n",yaxt="n")
 lines(Aa)


 i <- 0             #  ; i.max<-100
 cat("exit by button Exit\n")
 if(delay < 0.015) delay <- 0.015
 repeat{
   Sys.sleep(delay)
   choice <- tclvalue(Rot)
   if(choice=="exit"
                    #  || ((i<-i+1)>i.max)
                     ){ break }
   if(choice=="relax") next
   if(choice=="reset") {
     points(xa, pch=20, col="white"); lines(Aa, col="white")
     x <- x.o; A <- A.o; xa<-x[,2:3]; Aa<-A[,2:3]
     points(xa, pch=20, col="black"); lines(Aa, col="black")
     tclvalue(Rot)<-"relax"; next
   }
   switch(choice,
    "12" = ind<-c(1,3), "21" = ind<-c(2,1), "23" = ind<-c(1,2),
    "32" = ind<-c(3,1), "41" = ind<-c(3,2), "43" = ind<-c(2,3)
   )
   x[,ind] <- x[,ind]%*%rot; A[,ind] <- A[,ind]%*%rot
   points(xa, pch=20, col="white"); lines(Aa, col="white")
   xa<-x[,2:3]; Aa<-A[,2:3]
   points(xa, pch=20, col="black"); lines(Aa, col="black")
 }

 tkdestroy(topl)
 "control widget closed"

}

