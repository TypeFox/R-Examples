about <- function () {

ab.tt <- tktoplevel()
tktitle(ab.tt) <- c("About")
tkwm.resizable(ab.tt, 0,0)
ico <- tk2ico.load(res = "information")
tk2ico.set(ab.tt, ico)

invisibleframe <- tkwidget(ab.tt, "labelframe", text= "", padx=10, pady=10)
tkconfigure(invisibleframe, borderwidth=0, background="white")
tkgrid(invisibleframe, column=0, row=0)

tk2ico.set(ab.tt, Envir$ico)

#Affichage du logo FBdataM
TTAimg <- tclVar()
tcl("image", "create", "photo", TTAimg, file=file.path(path.package("TTAinterfaceTrendAnalysis"),"aide","TTAinterface.gif",fsep=.Platform$file.sep))
imgAsLabel11 <- tklabel(invisibleframe, image=TTAimg)
tkconfigure(imgAsLabel11, background="white")
tkgrid(imgAsLabel11, column=0, row=1, sticky="w", rowspan=7)

ab.text <- tklabel(invisibleframe, text= paste0("TTAinterface - GUI for R Temporal Trend Analysis"), justify="left")
tkgrid(ab.text, column=1, row=1, sticky="w")
tkconfigure(ab.text, font=tkfont.create(size= 10, weight="bold"), background="white" )

##A usage personnel (doit etre connecte)
#Envir$cran.version <- paste0("", available.packages()["TTAinterfaceTrendAnalysis", ][2])
#if ( Envir$pversion >= Envir$cran.version ) { v.text <- tklabel(invisibleframe, text=paste0("\n\n", "Congratulation your version is up to date!"))
#tkconfigure(v.text, background="white")
#tkgrid(v.text, column=1, row=2, sticky="wn") } else {
#v.text <- tklabel(invisibleframe, text=paste0("\n\n", "You should update the package..."))
#tkconfigure(v.text, background="white")
#tkgrid(v.text, column=1, row=2, sticky="wn") }
#ab.text6 <- tklabel(invisibleframe, text= paste0("Version installed: ", Envir$pversion, " build 2107", "\n", "Version CRAN:     ", Envir$cran.version), justify="left")

##A usage CRAN (diffusion)
ab.text6 <- tklabel(invisibleframe, text= paste0("Version ", Envir$pversion, " (", utils::packageDescription("TTAinterfaceTrendAnalysis")[5], ")"), justify="left")

tkconfigure(ab.text6, background="white")
tkgrid(ab.text6, column=1, row=2, sticky="wn")

emptyframe <- tkwidget(invisibleframe, "labelframe", text= "", padx=10, pady=10, height=15)
tkconfigure(emptyframe, borderwidth=0, background="white")
tkgrid(emptyframe, column=1, row=4)


ab.text3 <- tklabel(invisibleframe, text= paste0("\n", "David Devreker & Alain Lefebvre"), justify="left")
tkconfigure(ab.text3, background="white" )
tkgrid(ab.text3, column=1, row=5, sticky="wn")

ab.text2 <- tklabel(invisibleframe, text= paste0("CREATED BY :"), justify="left")
tkgrid(ab.text2, column=1, row=5, sticky="wn")
tkconfigure(ab.text2, font=tkfont.create(size= 9, weight="bold"), background="white" )


ab.text5 <- tklabel(invisibleframe, text= paste0("\n", "<david.devreker@ifremer.fr> or <alain.lefebvre@ifremer.fr>"), justify="left")
tkconfigure(ab.text5, background="white" )
tkgrid(ab.text5, column=1, row=6, sticky="wn")

ab.text4 <- tklabel(invisibleframe, text= paste0("CONTACTS :"), justify="left")
tkgrid(ab.text4, column=1, row=6, sticky="wn")
tkconfigure(ab.text4, font=tkfont.create(size= 9, weight="bold"), background="white" )

tkgrid(tk2label(ab.tt, text=""), column=0, row=1)
bouton.ok <- tk2button(ab.tt, text="OK", command=function() tkdestroy(ab.tt) )
tkgrid(bouton.ok, column=0, row=2)
tkgrid(tk2label(ab.tt, text=""), column=0, row=3)

ab.text7 <- tklabel(ab.tt, text= paste0("Created with R 3.2  "), justify="right")
tkconfigure(ab.text7, font=tkfont.create(size=7), foreground="grey40")
tkgrid(ab.text7, column=0, row=3, sticky="e")

}
#__________________________________________________________________________________________________________________________
#Fin du code
