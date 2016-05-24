samp.dist.snap<-function(parent = NULL, parent2 = NULL, biv.parent = NULL, stat = mean,stat2 = NULL, stat3 = NULL, stat4 = NULL, s.size=c(1,3,6,10,20,50), s.size2 = NULL, R=1000, 
func = NULL, xlab = expression(bar(x)),show.SE = TRUE, fits = NULL, show.fits = TRUE, xlim = NULL, ylim = NULL,...){
if(!is.null(s.size2)&(length(s.size)!=length(s.size2))) stop("length of s.size must equal length of size2")
L <- length(s.size)
if(L>12) {stop("L must be <= 12")} else
if(L==1) {par(mfrow=c(1,1),mar=c(5,4,1,1.5))} else
if(L==2) {par(mfrow=c(2,1),mar=c(5,4,1,1.5))} else
if(L==3) {par(mfrow=c(1,3),mar=c(5,4,1.5))} else
if(L==4) {par(mfrow=c(2,2),mar=c(5,4,1,1.5))} else
if(L==5|L==6) {par(mfrow=c(3,2),mar=c(5,4,2,1.5))} else
if(L==7|L==8|L==9) {par(mfrow=c(3,3),mar=c(5,4,2,1.0))}else
if(L==10|L==11|L==12) {par(mfrow=c(4,3),mar=c(5,4,2,1.0))}
if(L>12)stop("s.size vectors must have length <= 12", call. = FALSE) 
for(i in 1:L){
if(is.null(xlim)&is.null(ylim)){                     
samp.dist(parent = parent, parent2 = parent2, biv.parent = biv.parent, stat = stat, stat2 = stat2, stat3 = stat3, stat4 = stat4, s.size = s.size[i], 
s.size2 = s.size2[i], func = func, R = R, xlab = xlab, show.SE = show.SE, anim = FALSE, ...)}
if(!is.null(xlim)&is.null(ylim)){                     
samp.dist(parent = parent, parent2 = parent2, biv.parent = biv.parent, stat = stat, stat2 = stat2, stat3 = stat3, stat4 = stat4, s.size = s.size[i], 
s.size2 = s.size2[i], func = func, R = R, xlab = xlab, show.SE = show.SE, anim = FALSE, xlim = xlim, ...)}
if(!is.null(xlim)&!is.null(ylim)){                     
samp.dist(parent = parent, parent2 = parent2, biv.parent = biv.parent, stat = stat, stat2 = stat2, stat3 = stat3, stat4 = stat4, s.size = s.size[i], 
s.size2 = s.size2[i], func = func, R = R, xlab = xlab, show.SE = show.SE, anim = FALSE, xlim = xlim, ylim = ylim,...)}
if(show.fits == TRUE){
if(!is.null(fits))fits(s.size[i], s.size2[i])}
}
}

#--------------------------------- Simples stats -------------------------------------#


samp.dist.snap.tck1<-function(statc = "mean"){

local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkcheckbutton <- ttkcheckbutton
        tkentry <- ttkentry
        tkframe <- ttkframe
        tklabel <- ttklabel
        tkradiobutton <- ttkradiobutton
    }

tclServiceMode(FALSE)              
        dialog.sd <- function(){
        tt <- tktoplevel()
        tkwm.title(tt,"Sampling distributions")
        biv.parent.entry <- tkentry(tt, textvariable=Biv.parent, width = 16)
        parent.entry1 <- tkentry(tt, textvariable=Parent1, width = 16)
        parent.entry2 <- tkentry(tt, textvariable=Parent2, width = 16)
        s.size.entry<-tkentry(tt, textvariable=SS, width = 16)
        s.size.entry2<-tkentry(tt, textvariable=SS2, width = 16)
        stat.entry1<-tkentry(tt, textvariable=Stat1, width = 16)
        stat.entry2<-tkentry(tt, textvariable=Stat2, width = 16)
        stat.entry3<-tkentry(tt, textvariable=Stat3, width = 16)
        stat.entry4<-tkentry(tt, textvariable=Stat4, width = 16)
        func.entry<-tkentry(tt, textvariable=Func, width = 16)
        R.entry<-tkentry(tt, textvariable=Rep, width = 16)
        x.entry<-tkentry(tt, textvariable=Xlab, width = 16)
        fits.entry <- tkentry(tt, textvariable=fits, width = 16)
        x.lim.entry<-tkentry(tt, textvariable=xlim, width = 16)
        y.lim.entry<-tkentry(tt, textvariable=ylim, width = 16)
        
 
done <- tclVar(0)
show.SE<-tclVar(1)
show.fits <- tclVar(1)
        
        reset <- function()
        {
            tclvalue(Biv.parent)<-"NULL"
            tclvalue(Parent1)<-"NULL"
            tclvalue(Parent2)<-"NULL"
            tclvalue(SS)<-"0"
            tclvalue(SS2)<-"0"
            tclvalue(Stat1)<-"NULL"
            tclvalue(Stat2)<-"NULL"
            tclvalue(Stat3)<-"NULL"
            tclvalue(Stat4)<-"NULL"
            tclvalue(Func)<-"NULL"
            tclvalue(Rep)<-"NULL"
            tclvalue(Xlab)<-"NULL"
            tclvalue(xlim)<-"NULL"
            tclvalue(ylim)<-"NULL"
            tclvalue(fits)<-"NULL"
            }

        reset.but <- tkbutton(tt, text="Reset", command=reset)
        submit.but <- tkbutton(tt, text="Submit",command=function()tclvalue(done)<-1)
       
 
tw <- function(){
        tkdestroy(tt)
        samp.dist.snap.tck2(statc = statc)
        }       
               
        add.button <- tkbutton(tt, text="Add 2nd sample",command=tw)
        
        build <- function()
        {
            parent <-parse(text=tclvalue(Parent1))[[1]]
            s.size <-parse(text=tclvalue(SS))[[1]]
            stat <-parse(text=tclvalue(Stat1))[[1]]
            R <-tclvalue(Rep)
            x<-parse(text=tclvalue(Xlab))[[1]]
            R <-tclvalue(Rep)
            se <- as.logical(tclObj(show.SE))
            xlim<-parse(text=tclvalue(xlim))[[1]]
            ylim<-parse(text=tclvalue(ylim))[[1]]
            fits <- parse(text=tclvalue(fits))[[1]]
            show.fits <- as.logical(tclObj(show.fits))
            substitute(samp.dist.snap(parent = parent, s.size = s.size, stat = stat, R = as.numeric(R), xlab = x,show.SE = se, xlim = xlim, ylim = ylim, fits = fits, show.fits = show.fits))
           
             }
        


        se.cbut <- tkcheckbutton(tt, text="Show SE", variable=show.SE)
        fits.cbut <- tkcheckbutton(tt, text="Show fits", variable=show.fits)
        tkgrid(tklabel(tt,text="Sampling distribution snapshots"), columnspan = 2)
        tkgrid(tklabel(tt,text=""))
        tkgrid(tklabel(tt,text="Parent",font=c("Helvetica","9","bold"), width = 12),parent.entry1) 
        tkgrid(tklabel(tt,text="Sample sizes ",font=c("Helvetica","9","bold"), width = 12),s.size.entry)
        tkgrid(tklabel(tt,text="Stat",font=c("Helvetica","9","bold"), width = 12),stat.entry1)
        tkgrid(tklabel(tt,text="Iterations", width = 12),  R.entry)
        tkgrid(tklabel(tt,text="X-axis label", width = 12), x.entry)
        tkgrid(tklabel(tt,text="Xlim", width = 12), x.lim.entry)
        tkgrid(tklabel(tt,text="Ylim", width = 12), y.lim.entry)
        tkgrid(tklabel(tt,text=""), columnspan = 2)
        tkgrid(se.cbut, sticky = "w")
        tkgrid(fits.cbut,sticky="w")
        tkgrid(tklabel(tt,text="Fit(s)"), fits.entry, sticky = "w")
        tkgrid(tklabel(tt,text=""), columnspan = 2)
        tkgrid(add.button, columnspan = 2)
        tkgrid(submit.but, reset.but)
        
        tkbind(tt, "<Destroy>", function()tclvalue(done)<-2)
        tkwait.variable(done)

        if(tclvalue(done)=="2") stop("aborted")
        tkdestroy(tt)
        cmd <- build()
        eval.parent(cmd)
        tclServiceMode(TRUE)        
}
     if(statc == "mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(1,3,7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("10000")
      Xlab<- tclVar("expression(bar(x))")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL") 
      }
      
      if(statc == "median"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(1,3,7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("median")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(Median)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL") 
      }

if(statc == "trimmed mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      tr.mean <- function(x)mean(x,trim = 0.2)
      Stat1<-tclVar("tr.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(Trimmed.mean)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")}
          
if(statc == "Winsorized mean"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      win.mean <- function(x)mean(win(x))
      Stat1<-tclVar("win.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(Winsorized.mean)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")}
          
if(statc == "Huber estimator"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("huber.mu")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(Huber.estimator)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")}
           
if(statc == "H-L estimator"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rexp(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("HL.mean")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(H-L.estimator)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")}
            
      if(statc == "sd"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(3,7,10,20)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("sd")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(S)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL") 
      }
      
      if(statc == "var"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(3,7,10,20)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("var")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(S^2)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")  
      }
      
if(statc == "MAD"){

      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("mad")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(MAD)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")  
      }
          
if(statc == "IQR"){
      Biv.parent<-tclVar("NULL")
      Parent1<-tclVar("expression(rnorm(s.size))")
      Parent2<-tclVar("NULL")
      SS<-tclVar("c(7,10,20,50)")
      SS2<-tclVar("NULL")
      Stat1<-tclVar("IQR")
      Stat2<-tclVar("NULL")
      Stat3<-tclVar("NULL")
      Stat4<-tclVar("NULL")
      Func<-tclVar("NULL")
      Rep<-tclVar("1000")
      Xlab<- tclVar("expression(IQR)")
      ylim <- tclVar("NULL")
      xlim <- tclVar("NULL")
      fits <- tclVar("NULL")  
      }        
             
dialog.sd()
})
}



#-------------------------------- Complex stats -------------------------------------------#



samp.dist.snap.tck2<-function (statc = "mean") 
{
    local({
        have_ttk <- as.character(tcl("info", "tclversion")) >= 
            "8.5"
        if (have_ttk) {
            tkbutton <- ttkbutton
            tkcheckbutton <- ttkcheckbutton
            tkentry <- ttkentry
            tkframe <- ttkframe
            tklabel <- ttklabel
            tkradiobutton <- ttkradiobutton
        }
        tclServiceMode(FALSE)
        dialog.sd <- function() {
            tt <- tktoplevel()
            tkwm.title(tt, "Sampling distributions")
            biv.parent.entry <- tkentry(tt, textvariable = Biv.parent, 
                width = 16)
            parent.entry1 <- tkentry(tt, textvariable = Parent1, 
                width = 16)
            parent.entry2 <- tkentry(tt, textvariable = Parent2, 
                width = 16)
            s.size.entry <- tkentry(tt, textvariable = SS, width = 16)
            s.size.entry2 <- tkentry(tt, textvariable = SS2, 
                width = 16)
            stat.entry1 <- tkentry(tt, textvariable = Stat1, 
                width = 16)
            stat.entry2 <- tkentry(tt, textvariable = Stat2, 
                width = 16)
            stat.entry3 <- tkentry(tt, textvariable = Stat3, 
                width = 16)
            stat.entry4 <- tkentry(tt, textvariable = Stat4, 
                width = 16)
            func.entry <- tkentry(tt, textvariable = Func, width = 16)
            R.entry <- tkentry(tt, textvariable = Rep, width = 16)
            x.entry <- tkentry(tt, textvariable = Xlab, width = 16)
            fits.entry <- tkentry(tt, textvariable = fits, width = 16)
            x.lim.entry <- tkentry(tt, textvariable = xlim, width = 16)
            y.lim.entry <- tkentry(tt, textvariable = ylim, width = 16)
            done <- tclVar(0)
            show.SE <- tclVar(1)
            show.fits <- tclVar(1)
            reset <- function() {
                tclvalue(Biv.parent) <- "NULL"
                tclvalue(Parent1) <- "NULL"
                tclvalue(Parent2) <- "NULL"
                tclvalue(SS) <- "0"
                tclvalue(SS2) <- "0"
                tclvalue(Stat1) <- "NULL"
                tclvalue(Stat2) <- "NULL"
                tclvalue(Stat3) <- "NULL"
                tclvalue(Stat4) <- "NULL"
                tclvalue(Func) <- "NULL"
                tclvalue(Rep) <- "NULL"
                tclvalue(Xlab) <- "NULL"
                tclvalue(xlim) <- "NULL"
                tclvalue(ylim) <- "NULL"
                tclvalue(fits) <- "NULL"
            }
            reset.but <- tkbutton(tt, text = "Reset", command = reset)
            submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)
            build <- function() {
                biv.parent <- parse(text = tclvalue(Biv.parent))[[1]]
                parent <- parse(text = tclvalue(Parent1))[[1]]
                parent2 <- parse(text = tclvalue(Parent2))[[1]]
                s.size <- parse(text = tclvalue(SS))[[1]]
                stat <- parse(text = tclvalue(Stat1))[[1]]
                s.size2 <- parse(text = tclvalue(SS2))[[1]]
                stat2 <- parse(text = tclvalue(Stat2))[[1]]
                stat3 <- parse(text = tclvalue(Stat3))[[1]]
                stat4 <- parse(text = tclvalue(Stat4))[[1]]
                R <- tclvalue(Rep)
                x <- parse(text = tclvalue(Xlab))[[1]]
                func <- parse(text = tclvalue(Func))[[1]]
                se <- as.logical(tclObj(show.SE))
                xlim <- parse(text = tclvalue(xlim))[[1]]
                ylim <- parse(text = tclvalue(ylim))[[1]]
                fits <- parse(text = tclvalue(fits))[[1]]
                show.fits <- as.logical(tclObj(show.fits))
                substitute(samp.dist.snap(parent = parent, parent2 = parent2, 
                  biv.parent = biv.parent, s.size = s.size, s.size2 = s.size2, 
                  stat = stat, func = func, stat2 = stat2, stat3 = stat3, 
                  stat4 = stat4, R = as.numeric(R), xlab = x, 
                  show.SE = se, xlim = xlim, ylim = ylim, show.fits = show.fits, 
                  fits = fits))
            }
            se.cbut <- tkcheckbutton(tt, text = "Show SE", variable = show.SE)
            fits.cbut <- tkcheckbutton(tt, text = "Show fits", 
                variable = show.fits)
            tkgrid(tklabel(tt, text = "Sampling distribution snapshots"), 
                columnspan = 4)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "Bivariate parent"), 
                biv.parent.entry)
            tkgrid(tklabel(tt, text = "Parent 1", font = c("Helvetica", 
                "9", "bold")), parent.entry1, tklabel(tt, text = "Parent 2"), 
                parent.entry2)
            tkgrid(tklabel(tt, text = "Sample sizes ", font = c("Helvetica", 
                "9", "bold")), s.size.entry, tklabel(tt, text = "Sample sizes 2"), 
                s.size.entry2)
            tkgrid(tklabel(tt, text = "Stat", font = c("Helvetica", 
                "9", "bold")), stat.entry1, tklabel(tt, text = "Stat 2"), 
                stat.entry2)
            tkgrid(tklabel(tt, text = "Stat 3"), stat.entry3, 
                tklabel(tt, text = "Stat 4"), stat.entry4)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "Iterations"), 
                R.entry)
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "Function"), 
                func.entry)
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "X-axis label"), 
                x.entry)
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "Xlim"), 
                x.lim.entry)
            tkgrid(tklabel(tt, text = ""), tklabel(tt, text = "Ylim"), 
                y.lim.entry)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(se.cbut)
            tkgrid(fits.cbut, tklabel(tt, text = "Fit(s)"), fits.entry)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = ""), submit.but, reset.but)
            tkgrid(tklabel(tt, text = ""))
            tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
            tkwait.variable(done)
            if (tclvalue(done) == "2") 
                stop("aborted")
            tkdestroy(tt)
            cmd <- build()
            eval.parent(cmd)
            tclServiceMode(TRUE)
        }
        if (statc == "custom") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("NULL")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("NULL")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("NULL")
            Xlab <- tclVar("NULL")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "mean") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(1,3,7,10,20,50)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("10000")
            Xlab <- tclVar("expression(bar(x))")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "median") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(1,3,7,10,20,50)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("median")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("10000")
            Xlab <- tclVar("expression(Median)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "trimmed mean") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(7,10,20,50)")
            SS2 <- tclVar("NULL")
            tr.mean <- function(x) mean(x, trim = 0.2)
            Stat1 <- tclVar("tr.mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(Trimmed.mean)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "Winsorized mean") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar(("c(7,10,20,50)"))
            SS2 <- tclVar("NULL")
            win.mean <- function(x) mean(win(x))
            Stat1 <- tclVar("win.mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(Winsorized.mean)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "Huberestimator") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar(("c(7,10,20,50)"))
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("huber.mu")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(Huber.estimator)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "H-L estimator") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rexp(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar(("c(7,10,20,50)"))
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("HL.mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(H-L.estimator)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "sd") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,7,10,20)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("sd")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(S)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "var") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,7,10,20)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("var")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(S^2)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "MAD") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar(("c(7,10,20,50)"))
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("mad")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(MAD)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "IQR") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(7,10,20,50)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("IQR")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("NULL")
            Rep <- tclVar("1000")
            Xlab <- tclVar("expression(IQR)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "(n-1)S^2/sigma^2") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,7,10,20)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("var")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("function(s.dist, s.dist2, s.size, s.size2, sigma.sq = 1)((s.size - 1) * s.dist)/sigma.sq")
            Rep <- tclVar("10000")
            Xlab <- tclVar("expression((n - 1)*S^2/sigma^2)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("c(0,40)")
            x <- NULL
            suppressWarnings(rm(x))
            fits <- tclVar("function(s.size, s.size2)curve(dchisq(x, s.size - 1),from = 0, to = 40, add = TRUE, lwd = 2, col = gray(.3))")
        }
        if (statc == "F*") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("expression(rnorm(s.size2))")
            SS <- tclVar("c(5,7,12,20)")
            SS2 <- tclVar("c(7,5,10,15)")
            Stat1 <- tclVar("var")
            Stat2 <- tclVar("var")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("function(s.dist, s.dist2, s.size, s.size2) s.dist/s.dist2")
            Rep <- tclVar("2000")
            Xlab <- tclVar("expression(F.star)")
            ylim <- tclVar("c(0, 0.8)")
            xlim <- tclVar("c(0,20)")
            x <- NULL
            suppressWarnings(rm(x))
            fits <- tclVar("function(s.size, s.size2)curve(df(x, s.size - 1, s.size2 - 1),from = 0, to = 20, add = TRUE, col = gray(.3), lwd = 2)")
        }
        if (statc == "t* (1 sample)") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,7,10,20)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("mean")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("var")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("function(s.dist, s.dist3, s.size, s.size2)s.dist/sqrt(s.dist3/s.size)")
            Rep <- tclVar("10000")
            Xlab <- tclVar("expression(t.star)")
            ylim <- tclVar("c(0, 0.5)")
            xlim <- tclVar("c(-5, 5)")
            x <- NULL
            suppressWarnings(rm(x))
            fits <- tclVar("function(s.size, s.size2)({curve(dnorm(x),from = -10, to = 10, add = TRUE, lty = 2, lwd = 2, col = gray(.3)); curve(dt(x, s.size-1),from = -10,to = 10, add = TRUE, lty = 1, lwd = 2, col = gray(.6))})")
        }
        if (statc == "t* (2 sample)") {
            Biv.parent <- tclVar("NULL")
            Parent1 <- tclVar("expression(rnorm(s.size))")
            Parent2 <- tclVar("expression(rnorm(s.size))")
            SS <- tclVar("c(3,7,10,20)")
            SS2 <- tclVar("c(5,8,10,30)")
            Stat1 <- tclVar("mean")
            Stat2 <- tclVar("mean")
            Stat3 <- tclVar("var")
            Stat4 <- tclVar("var")
            Func <- tclVar("function(s.dist1, s.dist2, s.dist3, s.dist4, s.size = 6, s.size2 = s.size2)({MSE<-(((s.size - 1) * s.dist3) + ((s.size2 - 1) * s.dist4))/((s.size + s.size2) - 2);(s.dist1 - s.dist2)/(sqrt(MSE) * sqrt((1/s.size) + (1/s.size2)))})")
            Rep <- tclVar("10000")
            Xlab <- tclVar("expression(t.star)")
            ylim <- tclVar("c(0, 0.5)")
            xlim <- tclVar("c(-5, 5)")
            x <- NULL
            suppressWarnings(rm(x))
            fits <- tclVar("function(s.size, s.size2)({curve(dnorm(x),from = -10, to = 10, add = TRUE, lty = 2, lwd = 2, col = gray(.3)); curve(dt(x, (s.size + s.size2) - 2),from = -10,to = 10, add = TRUE, lty = 1, lwd = 2, col = gray(.6))})")
        }
        if (statc == "Pearson correlation") {

            Biv.parent <- tclVar("expression(rmvnorm(s.size, c(0,0), sigma = matrix(nrow=2,ncol=2, data =c(1,0,0,1))))")
            Parent1 <- tclVar("NULL")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,5,7,10,20,50)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("NULL")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("cor")
            Rep <- tclVar("2000")
            Xlab <- tclVar("expression(r)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        if (statc == "covariance") {

            Biv.parent <- tclVar("expression(rmvnorm(s.size, c(0,0), sigma = matrix(nrow=2,ncol=2, data =c(1,0,0,1))))")
            Parent1 <- tclVar("NULL")
            Parent2 <- tclVar("NULL")
            SS <- tclVar("c(3,5,7,10,20,50)")
            SS2 <- tclVar("NULL")
            Stat1 <- tclVar("NULL")
            Stat2 <- tclVar("NULL")
            Stat3 <- tclVar("NULL")
            Stat4 <- tclVar("NULL")
            Func <- tclVar("cov")
            Rep <- tclVar("2000")
            Xlab <- tclVar("expression(Covariance)")
            ylim <- tclVar("NULL")
            xlim <- tclVar("NULL")
            fits <- tclVar("NULL")
        }
        dialog.sd()
    })
}


