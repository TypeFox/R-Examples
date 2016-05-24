anm.mc.bvn.tck <- function(){

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
        dialog.ci <- function() {
            tt <- tktoplevel()
            tkwm.title(tt, "MCMC random walk")
            mu.entry <- tkentry(tt, textvariable = Mu, 
                width = 10)
            sigma.entry <- tkentry(tt, textvariable = Sigma, 
                width = 10)
            start.entry <- tkentry(tt, textvariable = Start, 
                width = 10)
            jump.entry <- tkentry(tt, textvariable = Jump, width = 10)
            xlim.entry <- tkentry(tt, textvariable = Xlim, 
                width = 10)
            ylim.entry <- tkentry(tt, textvariable = Ylim, 
                width = 10)
            length.entry <- tkentry(tt, textvariable = Length, 
                width = 10)
            int.entry <- tkentry(tt, textvariable = Int, width = 10)
            done <- tclVar(0)
            Par <- tclVar("Gibbs")
            reset <- function() {
                tclvalue(Mu) <- "c(0, 0)"
                tclvalue(Sigma) <- "matrix(2, 2, data = c(1, 0, 0, 1))"
                tclvalue(Start) <- "c(-4, 4)"
                tclvalue(Jump) <- "0.2"
                tclvalue(Xlim) <- "c(-4, 4)"
                tclvalue(Ylim) <- "c(-4, 4)"
                tclvalue(Length) <- "1000"
                tclvalue(Int) <- "0.01"
                tclvalue(Par) <- "Gibbs"
            }
            reset.but <- tkbutton(tt, text = "Reset", command = reset)
            submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 0)
            build <- function() {
                mu <- parse(text = tclvalue(Mu))[[1]]
                sigma <- parse(text = tclvalue(Sigma))[[1]]
                start <- parse(text = tclvalue(Start))[[1]] 
                jump <- tclvalue(Jump)
                part <- tclvalue(Par)
                if(part == "Gibbs") par.type <- "G"
                if(part == "Metropolis") par.type <- "M"
                if(part == "Metropolis-Hastings") par.type <- "MH"
                xlim <- parse(text = tclvalue(Xlim))[[1]]
                ylim <- parse(text = tclvalue(Ylim))[[1]]
                interval <- tclvalue(Int)
                length <- tclvalue(Length)
                substitute(anm.mc.bvn(mu = mu, sigma = sigma, jump.kernel = as.numeric(jump), 
                  sim = par.type, interval = as.numeric(interval), length = as.numeric(length), xlim = xlim, ylim = ylim))
            }
            alt.rbuts <- tkframe(tt)
            tkpack(tklabel(alt.rbuts, text = "Simulation"))
            for (i in c("Gibbs", "Metropolis", "Metropolis-Hastings")) {
                tmp <- tkradiobutton(alt.rbuts, text = i, variable = Par, 
                  value = i)
                tkpack(tmp, anchor = "w")
            }
            tkgrid(tklabel(tt, text = "MCMC random walk"), 
                columnspan = 2)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(tklabel(tt, text = "\u03BC", font=c("Helvetica","9","bold")), mu.entry)
            tkgrid(tklabel(tt, text = "\u03C3\u00b2", font=c("Helvetica","9","bold")), sigma.entry)
            tkgrid(tklabel(tt, text = "Start"), start.entry)
            tkgrid(tklabel(tt, text = "Jump kernel"), jump.entry)
            tkgrid(tklabel(tt, text = "Length"), length.entry)
            tkgrid(tklabel(tt, text = "Xlim"), xlim.entry)
            tkgrid(tklabel(tt, text = "Ylim"), ylim.entry)
            tkgrid(tklabel(tt, text = "Anim. interval"), 
                int.entry)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(alt.rbuts, columnspan = 2)
            tkgrid(tklabel(tt, text = ""))
            tkgrid(submit.but, reset.but, sticky = "w")
            tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
            tkwait.variable(done)
            if (tclvalue(done) == "2") 
                stop("aborted")
            tkdestroy(tt)
            cmd <- build()
            eval.parent(cmd)
            invisible(tclServiceMode(TRUE))
        }
        Mu <- tclVar("c(0, 0)")
        Sigma <- tclVar("matrix(2, 2, data = c(1, 0, 0, 1))")
        Start <- tclVar("c(-4, 4)")
        Jump <- tclVar("0.2")
        Length <- tclVar("1000")
        Xlim <- tclVar("c(-4, 4)")
        Ylim <- tclVar("c(-4, 4)")
        Int <- tclVar("0.01")
        dialog.ci()
    })
}