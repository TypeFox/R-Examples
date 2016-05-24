## illustrateCLT_tcl.R produces plots for the subseqent laws of T_n
## 
### 

require(tcltk) || stop("tcltk support is absent")
require(graphics); require(stats); require(distrTeach)
options("newDevice"=TRUE)

local({

    k    <- tclVar(1)
    X  <- tclVar("Norm()")


    k.sav <- kv <- 1
    X.sav  <- Xv <- "Norm()"
    renewK <- TRUE
    renewX <- TRUE

    replot <- function(...) {
        if(is.null(k)) kv <- 1
        if(renewK)
           {k.sav <<- kv <<- as.numeric(tclObj(k)); renewk <<- FALSE}
        if(renewX)
           {X.sav <<- Xv <<- tclvalue(X); 
            Distr0 <<- eval(parse(text=Xv)); renewX <<- FALSE}
        illustrateCLT.tcl(Distr = Distr0, k = kv, Distrname = Xv)
    }

    replot.maybe <- function(...)
    {
        if (as.numeric(tclObj(k)) != k.sav)
            renewK <<- TRUE
        if (tclvalue(X) != X.sav)
            renewX <<- TRUE 
            replot()
    }




    base <- tktoplevel()
    tkwm.title(base, gettext("Demo for CLT"))

    spec.frm <- tkframe(base,borderwidth=2)
    top.frm <- tkframe(spec.frm)
    bottom.frm <- tkframe(spec.frm)

    ## 1 top frame:
    frame1 <-tkframe(top.frm, relief="groove", borderwidth=2)
    
    tkpack(tklabel (frame1, text=gettext("# of summands")))
    tkpack(tkscale(frame1, command=replot.maybe, from=1, to=100,
                   showvalue=TRUE, variable=k,
                   resolution=1, orient="horiz"))

    ## 1 bottom frame:
    frame2 <-tkframe(bottom.frm, relief="groove", borderwidth=2)
  
    tkpack(tklabel(frame2, text=gettext("Summation variable")))
    X.entry <- tkentry(frame2, textvariable=X, width="20")
    X.entry.Scroll <- tkscrollbar(frame2, orient="horizontal",
                repeatinterval=5, command = 
                function(...) tkxview(X.entry, ...))    
    tkconfigure(X.entry, xscrollcommand = 
                function(...) tkset(X.entry.Scroll, ...))

    buttonsFrame <- tkframe(bottom.frm, relief="groove", borderwidth=2)
    onOK <- function() {Xv <- tclvalue(X); replot.maybe()}
    OKbutton <- tkbutton(buttonsFrame, text=gettext("OK"), 
                        fg="darkgreen", command=onOK, default="active")
    
    tkpack(X.entry,X.entry.Scroll,OKbutton, fill="x")

    tkpack(frame1, fill="x")
    tkpack(frame2, buttonsFrame, fill="x")
    tkpack(top.frm,bottom.frm,side="top", anchor="n")

    ## `Bottom frame' (on base):
    q.but <- tkbutton(base,text=gettext("Quit"),
                      command=function()tkdestroy(base))

    tkpack(spec.frm, q.but)
})


