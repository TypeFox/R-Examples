anm.die.tck<-function () 
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
            tkwm.title(tt, "Frequentist probability and die throws")
            p.entry <- tkentry(tt, textvariable = P, width = 10)
            int.entry <- tkentry(tt, textvariable = Int, width = 10)
            throw.entry <- tkentry(tt, textvariable = Throws, width = 10)
            done <- tclVar(0)
            show.die<-tclVar(1)
            color<-tclVar(1)
            reset <- function() {
                tclvalue(P) <- "c(1/6,1/6,1/6,1/6,1/6,1/6)"
                tclvalue(Int) <- "0.01"
                tclvalue(Throw) <- "1000"
                tclvalue(show.die)<-"0"
                tclvalue(color)<-"0"
            }
            reset.but <- tkbutton(tt, text = "Reset", command = reset)
            submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)
            build <- function() {
                p <-parse(text=tclvalue(P))[[1]]
                interval <- tclvalue(Int)
                reps <- tclvalue(Throws)
                show.die <- as.logical(tclObj(show.die))
                color <- as.logical(tclObj(color))
                substitute(anm.die(rep=as.numeric(reps),p=p, 
                  interval = as.numeric(interval), show.die=show.die, cl = color))
            }
            nc.cbut <- tkcheckbutton(tt, text="Show die", variable=show.die)
            col.cbut <- tkcheckbutton(tt, text="Color", variable=color)
            tkgrid(tklabel(tt, text = "Frequentist probability\nand die throws", justify = "center"), 
                columnspan = 2)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = "P(Throw)"), p.entry)
            tkgrid(tklabel(tt, text = "Throws"), throw.entry)
            tkgrid(tklabel(tt, text = "Anim. int."), 
                int.entry)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(nc.cbut, col.cbut)
            tkgrid(submit.but, reset.but,sticky="w")
            tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
            tkwait.variable(done)
            if (tclvalue(done) == "2") 
                stop("aborted")
            tkdestroy(tt)
            cmd <- build()
            eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
        }
        P <- tclVar("rep(1/6,6)")
        Throws <- tclVar("1000")
        Int <- tclVar("0.01")
        dialog.sd()
    })
    }
