
## function for eliciting priors for psi
ElicitPsi <- function(C00, C01, C10, C11, maxvalue=100,
                      a00=.25, a01=.25, a10=.25, a11=.25, nsamp=50000,
                      output.object="output.SimpleTable"){

  
  b00tcl  <- tclVar(1)
  c00tcl  <- tclVar(1)

  b01tcl  <- tclVar(1)
  c01tcl  <- tclVar(1)

  b10tcl  <- tclVar(1)
  c10tcl  <- tclVar(1)

  b11tcl  <- tclVar(1)
  c11tcl  <- tclVar(1)
  
  
  bw.sav <- 1 # in case replot.maybe is called too early
  
  
  
  replot.maybe <- function(...){
    b00 <<- as.numeric(tclObj(b00tcl))
    c00 <<- as.numeric(tclObj(c00tcl))

    b01 <<- as.numeric(tclObj(b01tcl))
    c01 <<- as.numeric(tclObj(c01tcl))

    b10 <<- as.numeric(tclObj(b10tcl))
    c10 <<- as.numeric(tclObj(c10tcl))
    
    b11 <<- as.numeric(tclObj(b11tcl))
    c11 <<- as.numeric(tclObj(c11tcl))

    gridval <- seq(from=.001, to=.999, by=.001)

    y00 <- dbeta(gridval, b00, c00)
    z00 <- dbeta(gridval, c00, b00)

    y01 <- dbeta(gridval, b01, c01)
    z01 <- dbeta(gridval, c01, b01)

    y10 <- dbeta(gridval, b10, c10)
    z10 <- dbeta(gridval, c10, b10)

    y11 <- dbeta(gridval, b11, c11)
    z11 <- dbeta(gridval, c11, b11)

    par(mfrow=c(2,4))

    
    if (b00+c00-2 <= C00){
      mycol00="black"
    }
    else{
      mycol00="red"
    }

    if (b01+c01-2 <= C01){
      mycol01="black"
    }
    else{
      mycol01="red"
    }

    if (b10+c10-2 <= C10){
      mycol10="black"
    }
    else{
      mycol10="red"
    }

        if (b11+c11-2 <= C11){
      mycol11="black"
    }
    else{
      mycol11="red"
    }

    
    plot(gridval,
         y00,
         type="l", ylab="Density", mex=1.5,
         xlab=expression(paste(psi["00"],
             ": Fraction of Helped in (X=0, Y=0)", sep="")),
         col.lab=mycol00, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,y00,0), col="deepskyblue4")

    plot(gridval,
         y01,
         type="l", ylab="Density",
         xlab=expression(paste(psi["01"],
             ": Fraction of Always Succeed in (X=0, Y=1)", sep="")),
         col.lab=mycol01, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,y01,0), col="orangered2")
        
    plot(gridval,
         y10,
         type="l", ylab="Density",
         xlab=expression(paste(psi["10"],
             ": Fraction of Hurt in (X=1, Y=0)", sep="")),
         col.lab=mycol10, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,y10,0), col="forestgreen")

    plot(gridval,
         y11,
         type="l", ylab="Density",
         xlab=expression(paste(psi["11"],
             ": Fraction of Always Succeed in (X=1, Y=1)", sep="")),
         col.lab=mycol11, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,y11,0), col="purple4")


    
    plot(gridval,
         z00,
         type="l", ylab="Density",
         xlab=expression(paste(1-psi["00"],
             ": Fraction of Never Succeed in (X=0, Y=0)", sep="")),
         col.lab=mycol00, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,z00,0), col="deepskyblue4")

    plot(gridval,
         z01,
         type="l", ylab="Density",
         xlab=expression(paste(1-psi["01"],
             ": Fraction of Hurt in (X=0, Y=1)", sep="")),
         col.lab=mycol01, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,z01,0), col="orangered2")


    plot(gridval,
         z10,
         type="l", ylab="Density",
         xlab=expression(paste(1-psi["10"],
             ": Fraction of Never Succeed in (X=1, Y=0)", sep="")),
         col.lab=mycol10, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,z10,0), col="forestgreen")

    plot(gridval,
         z11,
         type="l", ylab="Density",
         xlab=expression(paste(1-psi["11"],
             ": Fraction of Helped in (X=1, Y=1)", sep="")),
         col.lab=mycol11, cex.lab=1.5)
    polygon(c(0,gridval,1), c(0,z11,0), col="purple4")

    
  }









  replot.effect <- function(...){
    b00 <<- as.numeric(tclObj(b00tcl))
    c00 <<- as.numeric(tclObj(c00tcl))

    b01 <<- as.numeric(tclObj(b01tcl))
    c01 <<- as.numeric(tclObj(c01tcl))

    b10 <<- as.numeric(tclObj(b10tcl))
    c10 <<- as.numeric(tclObj(c10tcl))
    
    b11 <<- as.numeric(tclObj(b11tcl))
    c11 <<- as.numeric(tclObj(c11tcl))

    gridval <- seq(from=.001, to=.999, by=.001)

    y00 <- dbeta(gridval, b00, c00)
    z00 <- dbeta(gridval, c00, b00)

    y01 <- dbeta(gridval, b01, c01)
    z01 <- dbeta(gridval, c01, b01)

    y10 <- dbeta(gridval, b10, c10)
    z10 <- dbeta(gridval, c10, b10)

    y11 <- dbeta(gridval, b11, c11)
    z11 <- dbeta(gridval, c11, b11)

    par(mfrow=c(2,4))

    plot(gridval,
         y00,
         type="l", ylab="Density", mex=1.5,
         xlab="Psi00: Fraction of Helped in (X=0, Y=0)")
    polygon(c(0,gridval,1), c(0,y00,0), col="deepskyblue4")

    plot(gridval,
         y01,
         type="l", ylab="Density",
         xlab="Psi01: Fraction of Always Succeed in (X=0, Y=1)")
    polygon(c(0,gridval,1), c(0,y01,0), col="orangered2")
        
    plot(gridval,
         y10,
         type="l", ylab="Density",
         xlab="Psi10: Fraction of Hurt in (X=1, Y=0)")
    polygon(c(0,gridval,1), c(0,y10,0), col="forestgreen")

    plot(gridval,
         y11,
         type="l", ylab="Density",
         xlab="Psi11: Fraction of Always Succeed in (X=1, Y=1)")
    polygon(c(0,gridval,1), c(0,y11,0), col="purple4")


    
    plot(gridval,
         z00,
         type="l", ylab="Density",
         xlab="1-Psi00: Fraction of Never Succeed in (X=0, Y=0)")
    polygon(c(0,gridval,1), c(0,z00,0), col="deepskyblue4")

    plot(gridval,
         z01,
         type="l", ylab="Density",
         xlab="1-Psi01: Fraction of Hurt in (X=0, Y=1)")
    polygon(c(0,gridval,1), c(0,z01,0), col="orangered2")


    plot(gridval,
         z10,
         type="l", ylab="Density",
         xlab="1-Psi10: Fraction of Never Succeed in (X=1, Y=0)")
    polygon(c(0,gridval,1), c(0,z10,0), col="forestgreen")

    plot(gridval,
         z11,
         type="l", ylab="Density",
         xlab="1-Psi11: Fraction of Helped in (X=1, Y=1)")
    polygon(c(0,gridval,1), c(0,z11,0), col="purple4")

    output <- analyze2x2(C00=C00, C01=C01, C10=C10, C11=C11,
                         a00=a00, a01=a01, a10=a10, a11=a11,
                         b00=b00, b01=b01, b10=b10, b11=b11,
                         c00=c00, c01=c01, c10=c10, c11=c11,
                         nsamp=nsamp)
    

    outstring <- paste(output.object, "<<- output")
    eval(parse(text=outstring))
    
  }






  
    
  base <- tktoplevel()
  tkwm.title(base, "Prior Parameters")
  
  spec.frm <- tkframe(base,borderwidth=2)
  farL.frm <- tkframe(spec.frm, borderwidth=6, relief="groove")
  nearL.frm <- tkframe(spec.frm, borderwidth=6, relief="groove")
  farR.frm <- tkframe(spec.frm, borderwidth=6, relief="groove")
  nearR.frm <- tkframe(spec.frm, borderwidth=6, relief="groove")
  

  frame1 <-tkframe(farL.frm, relief="groove", borderwidth=2)
  tkpack(tklabel(frame1,
                  text="  b00: Number of Pseudo Helped + 1 \n in (X=0, Y=0)  "))
  tkpack(tkscale(frame1, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=b00tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))

  frame2 <-tkframe(farL.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame2,
                  text="  c00: Number of Pseudo Never Succeed + 1 \n in (X=0, Y=0)  "))
  tkpack(tkscale(frame2, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=c00tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  
  
  
  frame3 <-tkframe(nearL.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame3,
                  text="  b01: Number of Pseudo Always Succeed + 1 \n in (X=0, Y=1)  "))
  tkpack(tkscale(frame3, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=b01tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  
  frame4 <-tkframe(nearL.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame4,
                  text="  c01: Number of Pseudo Hurt + 1 \n in (X=0, Y=1)  "))
  tkpack(tkscale(frame4, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=c01tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  
  
  frame5 <-tkframe(nearR.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame5,
                  text="  b10: Number of Pseudo Hurt + 1 \n in (X=1, Y=0)  "))
  tkpack(tkscale(frame5, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=b10tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))

  frame6 <-tkframe(nearR.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame6,
                  text="  c10: Number of Pseudo Never Succeed + 1 \n in (X=1, Y=0)  "))
  tkpack(tkscale(frame6, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=c10tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  
  
  
  frame7 <-tkframe(farR.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame7,
                  text="  b11: Number of Pseudo Always Succeed + 1 \n in (X=1, Y=1)  "))
  tkpack(tkscale(frame7, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=b11tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  
  frame8 <-tkframe(farR.frm, relief="groove", borderwidth=2)
  tkpack(tklabel (frame8,
                  text="  c11: Number of Pseudo Helped + 1 \n in (X=1, Y=1)  "))
  tkpack(tkscale(frame8, command=replot.maybe, from=(maxvalue/5000), to=maxvalue,
                 showvalue=T, variable=c11tcl,
                 resolution=(maxvalue/5000), orient="horiz", length=290))
  

  
  tkpack(frame1, frame2, side="top", fill="x")
  tkpack(frame3, frame4, side="top", fill="x")
  tkpack(frame5, frame6, side="top", fill="x")
  tkpack(frame7, frame8, side="top", fill="x")
  
  tkpack(farL.frm, nearL.frm, nearR.frm, farR.frm, side="left", anchor="n")
  

  but.frm <- tkframe(base, borderwidth=2)
  
  q.but <- tkbutton(but.frm, text="Quit",
                    command=function()tkdestroy(base))

  effect.but <- tkbutton(but.frm, text="Calculate Effects",
                         command=replot.effect)

  
  tkpack(effect.but, q.but, side="left", anchor="n")
  tkpack(spec.frm, but.frm)
    
} ## end elicitPsi
