.eyefit.prevR = function (vario, silent = FALSE) 
{
  ###############################################################################################
  # Cette fonction est extraite du package Geor    
  # Elle permet de realiser un ajustement du semi variograme a vue de nez
  # Elle a ete quelque peut modifiee. Dans cette version seuls les modeles "Exp", "Sph", "Gau", "Mat" sont autorises 
  # Cette fonction renvoie les parametres des semi variograms. Ces parametres sont utilises
  #     comme valeurs initiales du programme d'ajustement fit.variogram appele par la fonction krige quand on est en mode manual
  # Elle renvoie un objt contenant le modele et les parametres d'ajustement 
  ###############################################################################################
  requireNamespace("tcltk", quietly = TRUE) || stop("The package tcltk is required to use manual fit. Please install it.", domain="R-prevR")
    geterrmessage()
    done <- tcltk::tclVar(0)
    eyefit.env <- new.env()
    assign("eyefit.tmp", list(), envir = eyefit.env)
    dmax <- max(vario$u)
    kappa1 <- tcltk::tclVar("0.5")
    kappa2 <- tcltk::tclVar("1.0")
    kernel <- tcltk::tclVar("exponential")
    mdist <- tcltk::tclVar(max(vario$u))
    nugget <- tcltk::tclVar(0.1 * (max(vario$v)))
    sill <- tcltk::tclVar(0.8 * (max(vario$v)))
    range <- tcltk::tclVar(dmax/3)
    replot <- function(...) {
        k <- as.character(tcltk::tclObj(kernel))
        kp1 <- as.numeric(tcltk::tclObj(kappa1))
        kp2 <- as.numeric(tcltk::tclObj(kappa2))
        maxdist <- as.numeric(tcltk::tclObj(mdist))
        sigmasq <- as.numeric(tcltk::tclObj(sill))
        phi <- as.numeric(tcltk::tclObj(range))
        tausq <- as.numeric(tcltk::tclObj(nugget))
        eval(substitute(plot(vario)))
        fit <- get("eyefit.tmp", envir = eyefit.env)
        lapply(fit, function(x) geoR::lines.variomodel(seq(0, maxdist, 
            l = 100), cov.model = x$cov.model, kappa = x$kappa, 
            cov.pars = x$cov.pars, nug = x$nug, max.dist = x$max.dist))
        if (k == "gneiting.matern" | k == "gencauchy") {
            geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k, 
                kappa = c(kp1, kp2), cov.pars = c(sigmasq, phi), 
                nug = tausq, max.dist = maxdist)
        }
        else if (k == "powered.exponential" || k == "cauchy" || 
            k == "matern") {
            geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k, 
                kappa = kp1, cov.pars = c(sigmasq, phi), nug = tausq, 
                max.dist = maxdist)
        }
        else geoR::lines.variomodel(x = seq(0, maxdist, l = 100), cov.model = k, 
            cov.pars = c(sigmasq, phi), nug = tausq, max.dist = maxdist)
    }
    redraw <- function(...) {
        var <- as.character(tcltk::tclObj(kernel))
        if (var == "gneiting.matern" | var == "gencauchy") {
            tcltk::tkconfigure(entry.kappa1, state = "normal")
            tcltk::tkconfigure(ts5, state = "normal")
            tcltk::tkfocus(entry.kappa1)
            tcltk::tkconfigure(entry.kappa2, state = "normal")
            #tcltk::tkconfigure(ts6, state = "normal")
        }
        else if (var == "powered.exponential" || var == "cauchy" || 
            var == "matern") {
            tcltk::tkconfigure(entry.kappa1, state = "normal")
            tcltk::tkconfigure(ts5, state = "normal")
            tcltk::tkfocus(entry.kappa1)
            tcltk::tkconfigure(entry.kappa2, state = "disabled")
            #tcltk::tkconfigure(ts6, state = "disabled")
        }
        else {
            tcltk::tkconfigure(ts5, state = "disabled")
            #tcltk::tkconfigure(ts6, state = "disabled")
            tcltk::tkconfigure(entry.kappa1, state = "disabled")
            tcltk::tkconfigure(entry.kappa2, state = "disabled")
        }
        replot()
    }
    base <- tcltk::tktoplevel()
    tcltk::tkwm.title(base, "Eyefit 1.0")
    spec.frm <- tcltk::tkframe(base, borderwidth = 2)
    left.frm <- tcltk::tkframe(spec.frm)
    right.frm <- tcltk::tkframe(spec.frm)
    frame1 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame1, text = gettext("Parameters",domain="R-prevR")), fill = "both", 
        side = "top")
    entry.mdist <- tcltk::tkentry(frame1, width = "4", textvariable = mdist)
    tcltk::tkpack(ts1 <- tcltk::tkscale(frame1, label = gettext("Max. Distance",domain="R-prevR"), command = replot, 
        from = 0, to = dmax, showvalue = 0, variable = mdist, 
        resolution = 0.01, orient = "horiz", relief = "groove"), 
        fill = "both", expand = 1, padx = 3, pady = 2, ipadx = 3, 
        ipady = 2, side = "left")
    tcltk::tkpack(entry.mdist, fill = "none", side = "right")
    frame3 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame3, text = gettext("Cov. Parameters",domain="R-prevR")), fill = "both", 
        side = "top")
    entry.sill <- tcltk::tkentry(frame3, width = "4", textvariable = sill)
    tcltk::tkpack(ts2 <- tcltk::tkscale(frame3, label = "Sill (sigmasq):", 
        command = replot, from = 0, to = 2 * max(vario$v), showvalue = 0, 
        variable = sill, resolution = 0.01, orient = "horiz", 
        relief = "groove"), fill = "none", expand = 1, padx = 3, 
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tcltk::tkpack(entry.sill, side = "right")
    frame4 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    entry.range <- tcltk::tkentry(frame4, width = "4", textvariable = range)
    tcltk::tkpack(ts3 <- tcltk::tkscale(frame4, label = "Range (phi):", command = replot, 
        from = 0, to = 2 * dmax, showvalue = 1, variable = range, 
        resolution = 0.01, orient = "horiz", relief = "groove"), 
        fill = "x", expand = 1, padx = 3, pady = 2, ipadx = 3, 
        ipady = 2, side = "left")
    tcltk::tkpack(entry.range, side = "right")
    frame5 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame5, text = "Nugget"), fill = "both", side = "top")
    entry.nugget <- tcltk::tkentry(frame5, width = "4", textvariable = nugget)
    tcltk::tkpack(ts4 <- tcltk::tkscale(frame5, label = "Nugget (tausq):", 
        command = replot, from = 0, to = 2 * max(vario$v), showvalue = 1, 
        variable = nugget, resolution = 0.01, orient = "horiz", 
        relief = "groove"), fill = "x", expand = 1, padx = 3, 
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tcltk::tkpack(entry.nugget, side = "right")
    frame6 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame6, text = "Kappa"), fill = "both", side = "top")
    entry.kappa1 <- tcltk::tkentry(frame6, width = "4", textvariable = kappa1, 
        state = "disabled")
    tcltk::tkpack(ts5 <- tcltk::tkscale(frame6, label = "Kappa 1:", command = replot, 
        from = 0, to = 10, showvalue = 1, variable = kappa1, 
        state = "disabled", resolution = 0.01, orient = "horiz", 
        relief = "groove"), fill = "x", expand = 1, padx = 3, 
        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tcltk::tkpack(entry.kappa1, side = "right", fill = "none")
    frame7 <- tcltk::tkframe(left.frm, relief = "groove", borderwidth = 2)
    entry.kappa2 <- tcltk::tkentry(frame7, width = "4", textvariable = kappa2, 
        state = "disabled")
#    tcltk::tkpack(ts6 <- tcltk::tkscale(frame7, label = "Kappa 2:", command = replot, 
#        from = 0, to = 10, showvalue = 1, variable = kappa2, 
#        state = "disabled", resolution = 0.01, orient = "horiz", 
#        relief = "groove"), fill = "x", expand = 1, padx = 3, 
#        pady = 2, ipadx = 3, ipady = 2, side = "left")
    tcltk::tkpack(entry.kappa2, side = "right", fill = "none")
    frame2 <- tcltk::tkframe(right.frm, relief = "groove", borderwidth = 2)
    tcltk::tkpack(tcltk::tklabel(frame2, text = gettext("Function",domain="R-prevR")))
    for (i in c("exponential", "gaussian",  "matern", "power", "spherical")) {
        tmp <- tcltk::tkradiobutton(frame2, command = redraw, text = i, 
            value = i, variable = kernel)
        tcltk::tkpack(tmp, anchor = "w")
    }
    OnOK <- function() {
        replot()
    }
    OnQuit <- function() {
        k <- as.character(tcltk::tclObj(kernel))
        kp1 <- as.numeric(tcltk::tclObj(kappa1))
        if (k == "gneiting.matern") 
            kp2 <- as.numeric(tcltk::tclObj(kappa2))
        else kp2 <- NULL
        maxdist <- as.numeric(tcltk::tclObj(mdist))
        sigmasq <- as.numeric(tcltk::tclObj(sill))
        phi <- as.numeric(tcltk::tclObj(range))
        tausq <- as.numeric(tcltk::tclObj(nugget))
        aux <- list(cov.model = k, cov.pars = c(sigmasq, phi), 
            nugget = tausq, kappa = c(kp1, kp2), lambda = vario$lambda, 
            trend = vario$trend, practicalRange = geoR::practicalRange(cov.model = k, 
                phi = phi, kappa = kp1), max.dist = maxdist)
        oldClass(aux) <- "variomodel"
        assign("eyefit.tmp", c(get("eyefit.tmp", envir = eyefit.env), 
            list(aux)), envir = eyefit.env)
        tcltk::tclvalue(done) <- 2
    }
    OnClear <- function(aux = vario) {
        assign("eyefit.tmp", list(), envir = eyefit.env)
        plot(aux)
    }
    OnSave <- function() {
        k <- as.character(tcltk::tclObj(kernel))
        kp1 <- as.numeric(tcltk::tclObj(kappa1))
        if (k == "gneiting.matern") 
            kp2 <- as.numeric(tcltk::tclObj(kappa2))
        else kp2 <- NULL
        maxdist <- as.numeric(tcltk::tclObj(mdist))
        sigmasq <- as.numeric(tcltk::tclObj(sill))
        phi <- as.numeric(tcltk::tclObj(range))
        tausq <- as.numeric(tcltk::tclObj(nugget))
        aux <- list(cov.model = k, cov.pars = c(sigmasq, phi), 
            nugget = tausq, kappa = c(kp1, kp2), lambda = vario$lambda, 
            trend = vario$trend, practicalRange = geoR::practicalRange(cov.model = k, 
                phi = phi, kappa = kp1), max.dist = maxdist)
        oldClass(aux) <- "variomodel"
        assign("eyefit.tmp", c(get("eyefit.tmp", envir = eyefit.env), 
            list(aux)), envir = eyefit.env)
        replot()
    }
    tcltk::tkpack(frame1, frame3, frame4, frame5, frame6, frame7, fill = "x")
    tcltk::tkpack(frame2, fill = "x")
    tcltk::tkpack(left.frm, right.frm, side = "left", anchor = "n")
    c.but <- tcltk::tkbutton(base, text = "Clear", command = function() {
        OnClear(vario)
    })
    q.but <- tcltk::tkbutton(base, text = gettext("Choose this variogram",domain="R-prevR"), command = OnQuit)
    save.but <- tcltk::tkbutton(base, text = "Save", command = OnSave)
    tcltk::tkpack(spec.frm)
    tcltk::tkpack(q.but, side = "right")
    #tcltk::tkpack(c.but, side = "left")
    #tcltk::tkpack(save.but, side = "right")
    replot()
    tcltk::tkbind(entry.kappa1, "<Return>", function() {
        replot()
    })
    tcltk::tkbind(entry.kappa2, "<Return>", function() {
        replot()
    })
    tcltk::tkbind(base, "<Destroy>", function() tcltk::tclvalue(done) <- 2)
    tcltk::tkwait.variable(done)
    tcltk::tkdestroy(base)
    if (!silent) {
        fit <- get("eyefit.tmp", envir = eyefit.env)
        oldClass(fit) <- "eyefit"
        return(fit)
    }
    else return(invisible())
}
#out        = eyefit.prevR(varioGeoR)
#print(as.vgm.variomodel(out[[length(out)]]))