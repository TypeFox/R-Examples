fagan.plot<-function(probs.pre.test, LR, test.result="+") {
    opar <- par(no.readonly = T)
    on.exit(par(opar))
    par(mar = c(1.5, 6, 2, 6))
    stato <- ifelse(test.result == "+", "disease", "no disease")
    if (probs.pre.test > 1 | probs.pre.test < 0 | LR < 0 | is.infinite(LR) |
        is.nan(LR) | test.result %in% c("+", "-") == F) {
      stop("wrong values !!")
    } else {
      logits <- function(p) log(p/(1 - p))
    }
    logits.pre <- logits(probs.pre.test)
    logits.post <- log(LR) + logits.pre
    probs.post.test <- exp(logits.post)/(1 + exp(logits.post))
    compl.logit.pre <- logits(1 - probs.pre.test)
    LR.vec <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2,
        0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    prob.vec <- c(0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
        0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
        0.8, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.993, 0.995,
        0.997, 0.998, 0.999)
    plot(0, 0, type = "n", ylim = range(logits(prob.vec)), axes = F,
        xlab = "", ylab = "")
    axis(2, rev(logits(prob.vec)), 100 * prob.vec, pos = -1,
        las = 1, cex.axis = 0.7)
    axis(2, rev(logits(prob.vec)), 100 * prob.vec, pos = -1,
        tck = 0.03, labels = F)
    axis(4, logits(prob.vec), 100 * prob.vec, pos = 1, las = 1,
        cex.axis = 0.7)
    axis(4, logits(prob.vec), 100 * prob.vec, pos = 1, tck = 0.03,
        labels = F)
    axis(2, log(LR.vec[1:10])/2, LR.vec[1:10], pos = 0, las = 1,
        cex.axis = 0.7)
    axis(2, log(LR.vec[1:10])/2, LR.vec[1:10], pos = 0, tck = 0.03,
        labels = F)
    axis(4, log(LR.vec[10:19])/2, LR.vec[10:19], pos = 0, las = 1,
        cex.axis = 0.7)
    axis(4, log(LR.vec[10:19])/2, LR.vec[10:19], pos = 0, tck = 0.03,
        labels = F)
    text(0, 4.5, "Likelihood ratio", cex = 1.2)
    segments(-1, compl.logit.pre, 1, logits.post, lwd = 1.5,
        col = 2)
    mtext(side = 2, text = "Pre test probability(%)", line = 2,
        cex = 1.2)
    mtext(side = 4, text = "Post test probability(%)", line = 2,
        cex = 1.2, las = 3)
    title(main = "Fagan's nomogram")
    text(0, -6.3, paste("Pre test prob. of disease =", round(100 *
        probs.pre.test, 2), "% \n", "Likelihood ratio ", "=",
        round(LR, 2), "\n", "Post test prob. of", stato, "=",
        ifelse(test.result == "+", round(100 * probs.post.test,
            2), round(100 * (1 - probs.post.test), 2)), "%"),
        cex = 0.7)
}



plotFagan.old<-function(){
  refresh.code <- function(...) {
    probs.pre.test <- slider(no = 1)
    LR <- slider(no = 2)
    test.result <- slider(obj.name = "test.result")
    fagan.plot(probs.pre.test, LR, test.result)
  }
  slider(refresh.code, sl.names =
         c("pre test probability", "Likelihood Ratio"),
         sl.mins = c(0, 0.01), sl.maxs = c(1, 100),
         title = "Bayes nomogram",
         sl.defaults = c(0.5, 1),
         sl.deltas = c(0.01, 0.01), but.functions = list(function(...) {
           slider(obj.name = "test.result", obj.value = "+")
           refresh.code()
         }, function(...) {
           slider(obj.name = "test.result", obj.value = "-")
           refresh.code()
         }), but.names = c("positive result", "negative result"))
  slider(obj.name = "test.result", obj.value = "+")
  invisible(NULL)
}


plotFagan2.old<-function(){
  refresh.code <- function(...) {
    probs.pre.test <- slider(no = 1)
    LR <- slider(no=2)/(1-slider(no=3))
    test.result <- slider(obj.name = "test.result")
    fagan.plot(probs.pre.test, LR, test.result)
  }
  slider(refresh.code, sl.names = c("pre test probability",
                         "Sensitivity","Specificity"), sl.mins = c(0, 0, 0),
         sl.maxs = c(1, 1, 1), title = "Bayes nomogram",
         sl.defaults = c(0.5, 0.95, 0.95),
         sl.deltas = c(0.01, 0.001, 0.001),
         but.functions = list(function(...) {
           slider(obj.name = "test.result", obj.value = "+")
           refresh.code()
         }, function(...) {
           slider(obj.name = "test.result", obj.value = "-")
           refresh.code()
         }), but.names = c("positive result", "negative result"))
  slider(obj.name = "test.result", obj.value = "+")
  invisible(NULL)
}


plotFagan <- function(hscale=1.5, vscale=1.5, wait=FALSE) {

  if( !requireNamespace('tkrplot', quietly=TRUE) ) stop('This function depends on the tkrplot package being available')

  ppt <- tcltk::tclVar()
  tcltk::tclvalue(ppt) <- 0.5
  lr <- tcltk::tclVar()
  tcltk::tclvalue(lr) <- 1
  tr <- tcltk::tclVar()
  tcltk::tclvalue(tr) <- '+'

  hsc <- tcltk::tclVar()
  tcltk::tclvalue(hsc) <- hscale
  vsc <- tcltk::tclVar()
  tcltk::tclvalue(vsc) <- hscale

  replot <- function(...) {
    probs.pre.test <- as.numeric(tcltk::tclvalue(ppt))
    LR <- as.numeric(tcltk::tclvalue(lr))
    test.result <- tcltk::tclvalue(tr)
    fagan.plot(probs.pre.test, LR, test.result)
  }

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, "Fagan Plot Demo")

  img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
  tcltk::tkpack(img, side='top')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Pre Test Probability: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=ppt, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=0, to=1, resolution=.01), side='right')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Likelihood Ratio: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=lr, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=0.01, to=100, resolution=.01), side='right')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tkcheckbutton(fr, text='Positive Test Result', variable=tr,
                       onvalue='+', offvalue='-',
                       command=function(...) tkrplot::tkrreplot(img,
                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                         vscale=as.numeric(tcltk::tclvalue(vsc)) )
                       ),
         side='left')

  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tkbutton(tfr, text="Refresh", command=function() tkrplot::tkrreplot(img,
                                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                                         vscale=as.numeric(tcltk::tclvalue(vsc)) ) ),
                  side='left',anchor='s')

  tcltk::tkpack(tcltk::tkbutton(tfr, text="Exit", command=function()tcltk::tkdestroy(tt)),
             side='right',anchor='s')


  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tklabel(tfr,text="Hscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=hsc,width=6), side='left')
  tcltk::tkpack(tcltk::tklabel(tfr,text="      Vscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=vsc,width=6), side='left')

  if(wait) {
    tcltk::tkwait.window(tt)
    return( list(ppt = as.numeric(tcltk::tclvalue(ppt)),
                 lr = as.numeric(tcltk::tclvalue(lr)),
                 tr = tcltk::tclvalue(tr)
                 ))
  } else {
    return(invisible(NULL))
  }
}





plotFagan2 <- function(hscale=1.5, vscale=1.5, wait=FALSE) {

  if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

  ppt <- tcltk::tclVar()
  tcltk::tclvalue(ppt) <- 0.5
  sens <- tcltk::tclVar()
  tcltk::tclvalue(sens) <- 0.5
  spec <- tcltk::tclVar()
  tcltk::tclvalue(spec) <- 0.5
  tr <- tcltk::tclVar()
  tcltk::tclvalue(tr) <- '+'

  hsc <- tcltk::tclVar()
  tcltk::tclvalue(hsc) <- hscale
  vsc <- tcltk::tclVar()
  tcltk::tclvalue(vsc) <- hscale

  replot <- function(...) {
    probs.pre.test <- as.numeric(tcltk::tclvalue(ppt))
    sns <- as.numeric(tcltk::tclvalue(sens))
    spc <- as.numeric(tcltk::tclvalue(spec))
    test.result <- tcltk::tclvalue(tr)
    fagan.plot(probs.pre.test, sns/(1-spc), test.result)
  }

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, "Fagan Plot Demo")

  img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
  tcltk::tkpack(img, side='top')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Pre Test Probability: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=ppt, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=0, to=1, resolution=.01), side='right')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Sensitivity: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=sens, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=0, to=1, resolution=.01), side='right')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Specificity: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=spec, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=0, to=1, resolution=.01), side='right')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tkcheckbutton(fr, text='Positive Test Result', variable=tr,
                       onvalue='+', offvalue='-',
                       command=function(...) tkrplot::tkrreplot(img,
                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                         vscale=as.numeric(tcltk::tclvalue(vsc)) )
                       ),
         side='left')

  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tkbutton(tfr, text="Refresh", command=function() tkrplot::tkrreplot(img,
                                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                                         vscale=as.numeric(tcltk::tclvalue(vsc)) ) ),
                  side='left',anchor='s')

  tcltk::tkpack(tcltk::tkbutton(tfr, text="Exit", command=function()tcltk::tkdestroy(tt)),
             side='right',anchor='s')


  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tklabel(tfr,text="Hscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=hsc,width=6), side='left')
  tcltk::tkpack(tcltk::tklabel(tfr,text="      Vscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=vsc,width=6), side='left')

  if(wait) {
    tcltk::tkwait.window(tt)
    return( list(ppt = as.numeric(tcltk::tclvalue(ppt)),
                 sens = as.numeric(tcltk::tclvalue(sens)),
                 spec = as.numeric(tcltk::tclvalue(spec)),
                 tr = tcltk::tclvalue(tr)
                 ))
  } else {
    return(invisible(NULL))
  }
}




