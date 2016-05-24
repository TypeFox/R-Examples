#require(detrendeR)
#detrender()

y = function(x, first.year, last.year){
if (last.year<first.year) { temp<-first.year; first.year<-last.year ; last.year<-temp}
subset(x, as.integer(rownames(x))>=first.year & as.integer(rownames(x))<=last.year)
}

s = function (x){
tk.select.list = function (choices, preselect = NULL, multiple = FALSE, title = NULL) 
{

fontFixedWidth <- tkfont.create(family="courier",size=9)

    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if (!have_ttk) 
        ttkbutton <- tkbutton
    lvar <- tclVar()
    tclObj(lvar) <- choices
    oldmode <- tclServiceMode(FALSE)
    dlg <- tktoplevel()
	tkwm.resizable(dlg, 0, 0)
    tkwm.title(dlg, title)
    tkwm.deiconify(dlg)
    tkgrab.set(dlg)
    tkfocus(dlg)
   # if (!is.null(title) && nzchar(title)) {
   #     lab <- if (have_ttk) 
   #         ttklabel(dlg, text = title, foreground = "blue")
   #     else tklabel(dlg, text = title, fg = "blue")
   #     tkpack(lab, side = "top")
   # }
	 lab1<-ttklabel(dlg, text = "Series       First   Last   Span   ",font= fontFixedWidth)
	 tkpack(lab1, side = "top")
	 
    onOK <- function() {
        res <- 1L + as.integer(tkcurselection(box))
        ans.select_list <<- choices[res]
        tkgrab.release(dlg)
        tkdestroy(dlg)
    }
    onCancel <- function() {
        tkgrab.release(dlg)
        tkdestroy(dlg)
    }
    buttons <- tkframe(dlg)
    tkpack(buttons, side = "bottom")
    OK <- ttkbutton(buttons, text = gettext("OK"), width = 6, 
        command = onOK)
    Cancel <- ttkbutton(buttons, text = gettext("Cancel"), command = onCancel)
    tkpack(OK, Cancel, side = "left", fill = "x", padx = "2m")
    scht <- as.numeric(tclvalue(tkwinfo("screenheight", dlg))) - 
       200L
    ht <- min(length(choices), scht%/%20)
    box <- tklistbox(dlg, height = ht, listvariable = lvar, bg = "white", 
        setgrid = 1, selectmode = ifelse(multiple, "multiple", 
            "single"))
    tmp <- tcl("font", "metrics", tkcget(box, font = NULL))
    tmp <- as.numeric(sub(".*linespace ([0-9]+) .*", "\\1", tclvalue(tmp))) + 
        3
    ht <- min(length(choices), scht%/%tmp)
    tkdestroy(box)
    if (ht < length(choices)) {
        scr <- if (have_ttk) 
            ttkscrollbar(dlg, command = function(...) tkyview(box, 
                ...))
        else tkscrollbar(dlg, repeatinterval = 5, command = function(...) tkyview(box, 
            ...))
        box <- tklistbox(dlg, height = ht, width = 0, listvariable = lvar, 
            bg = "white", setgrid = 1, selectmode = ifelse(multiple, 
                "multiple", "single"),font=fontFixedWidth, yscrollcommand = function(...) tkset(scr, 
                ...))
        tkpack(box, side = "left", fill = "both", expand = TRUE)
        tkpack(scr, side = "right", fill = "y")
    }
    else {
        box <- tklistbox(dlg, height = ht, width = 0, listvariable = lvar, 
            bg = "white",font=fontFixedWidth, selectmode = ifelse(multiple, "multiple", 
                "single"))
        tkpack(box, side = "left", fill = "both")
    }
    preselect <- match(preselect, choices)
    preselect <- preselect[preselect > 0L] - 1L
    if (length(preselect)) {
        for (i in preselect) tkselection.set(box, i)
        tkyview(box, preselect[1L])
    }
    ans.select_list <- character()
	
	all = function () for (i in 1:length(choices)) tkselection.set(box, i-1)
	none = function () for (i in 1:length(choices)) tkselection.clear(box, i-1)
	
    tkbind(dlg, "<Destroy>", onCancel)
    tkbind(dlg, "<Double-ButtonPress-1>", onOK)
	tkbind(box, "<Control-a>", all)
    tkbind(box, "<Control-x>", none)
    tkfocus(box)
    tclServiceMode(oldmode)
    tkwait.window(dlg)
    Sys.sleep(0.1)
    if (!multiple && !length(ans.select_list)) 
        ans.select_list <- ""
    ans.select_list
}




yr.range = function(x) {
        yr.vec = as.numeric(names(x))
        mask = !is.na(x)
        range(yr.vec[mask])
    }
info.fun = function(x) {
		first<-yr.range(x)[1]
		last<-yr.range(x)[2]
        paste(format(first,width=6, justify="right"), format(last,width=6), format(last-first+1,width=6), " ")
    }

series<-paste( "",format(colnames(x), width=10),apply(x,2,info.fun), sep=" ")
selected.series<-tk.select.list(series, multiple=TRUE,preselect=series,title="Select the series to keep")
x[,series%in%selected.series]->temp

if(sum(series%in%selected.series)==0) return (invisible())
if(sum(series%in%selected.series)==1) {data.frame(x[,series%in%selected.series])->temp
colnames(temp)<-colnames(x)[series%in%selected.series]
rownames(temp)<-rownames(x)
apply(temp,1,sum, na.rm=T)->years
data.frame(temp[years>0,])->TEMP
colnames(TEMP)<-colnames(temp)
rownames(TEMP)<-rownames(temp)[years>0]
return(TEMP)
}

apply(temp,1,sum, na.rm=T)->years
temp[years>0,]
}
#data(co021)
#trimCol(co021)
