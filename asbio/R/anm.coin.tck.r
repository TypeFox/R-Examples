anm.coin.tck<-function () 
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
            tkwm.title(tt, "Frequentist probability and coin flips")
            p.entry <- tkentry(tt, textvariable = P, width = 10)
            int.entry <- tkentry(tt, textvariable = Int, width = 10)
            flip.entry <- tkentry(tt, textvariable = Flip, width = 10)
            done <- tclVar(0)
            show.coin<-tclVar(1)
            reset <- function() {
                tclvalue(P) <- "0.5"
                tclvalue(Int) <- "0.01"
                tclvalue(Flip) <- "1000"
                tclvalue(show.coin)<-"0"
            }
            tclServiceMode(FALSE)
            reset.but <- tkbutton(tt, text = "Reset", command = reset)
            submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)
            build <- function() {
                p.head <- tclvalue(P)
                interval <- tclvalue(Int)
                flips <- tclvalue(Flip)
                show.coin <- as.logical(tclObj(show.coin))
                substitute(anm.coin(p.head = as.numeric(p.head), 
                  interval = as.numeric(interval), flips = as.numeric(flips),show.coin=show.coin))
            }
            nc.cbut <- tkcheckbutton(tt, text="Show coin", variable=show.coin)
            tkgrid(tklabel(tt, text = "Frequentist probability\nand coin flips", justify = "center"), columnspan = 2)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = "P(Head)"), p.entry)
            tkgrid(tklabel(tt, text = "Flips"), flip.entry)
            tkgrid(tklabel(tt, text = "Anim. int."), int.entry)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(nc.cbut)
            tkgrid(submit.but, reset.but, sticky ="w")
            tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
            tkwait.variable(done)
            if (tclvalue(done) == "2") 
                stop("aborted")
            tkdestroy(tt)
            cmd <- build()
            eval.parent(cmd)
        invisible(tclServiceMode(TRUE))
        }
        P <- tclVar("0.5")
        Flip <- tclVar("1000")
        Int <- tclVar("0.01")
        dialog.sd()
    })
    }
