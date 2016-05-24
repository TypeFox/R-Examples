Sel_procedure <-
function(){
tt <- tktoplevel() # on cr\'ee une fenetre
tktitle(tt) <- "Number of raters"

textEntryVarTcl <- tclVar("2")
textEntryWidget <- tkentry(tt, width = 1, textvariable = textEntryVarTcl)
tkgrid(tklabel(tt, text = "Enter the number of raters: ", font="bold"), textEntryWidget)
tkgrid(tklabel(tt,text="    ")) # ligne vide pour aerer
OK.but <- tkbutton(tt,text="OK",command=function() tkdestroy(tt), relief="groove",borderwidth=3,width=7,bg="green4",fg="white")
tkgrid(OK.but)
tkfocus(tt)
tkwait.window(OK.but) # TRES IMPORTANT : l'ex\'ecution s'arr\^ete ici tant que l'utilisateur n'a pas cliqu\'e sur OK 

if (as.numeric(tclvalue(textEntryVarTcl))==2) {
 tt <- tktoplevel() # on cr\'ee une fenetre
 tktitle(tt) <- "Weighting"

 rb1 <- tkradiobutton(tt) # premier bouton radio
 rb2 <- tkradiobutton(tt) # second bouton radio
 rb3 <- tkradiobutton(tt) # troisieme bouton radio
 rbValue <- tclVar("unweighted") # cette variable contiendra au final le choix retenu par l'utilisateur via son bouton
 tkconfigure(rb1,variable=rbValue,value="unweighted")
 tkconfigure(rb2,variable=rbValue,value="equal")
 tkconfigure(rb3,variable=rbValue,value="squared")
 tkgrid(tklabel(tt,text="Weighting scheme:", font="bold"))
 tkgrid(tklabel(tt,text="Unweighted Kappa"),rb1)
 tkgrid(tklabel(tt,text="Linear weighting"),rb2)
 tkgrid(tklabel(tt,text="Quadratic Weighting"),rb3)

 tkgrid(tklabel(tt,text="    ")) # ligne vide pour aerer
 OK.but <- tkbutton(tt,text="OK",command=function() tkdestroy(tt), relief="groove",borderwidth=3,width=7,bg="green4",fg="white")
 tkgrid(OK.but)
 tkfocus(tt)
 tkwait.window(OK.but) # TRES IMPORTANT : l'ex\'ecution s'arr\^ete ici tant que l'utilisateur n'a pas cliqu\'e sur OK 

 res = list(as.numeric(tclvalue(textEntryVarTcl)), as.character(tclvalue(rbValue)))
 names(res) = c("Nb_raters", "Weighting")
 return(res)
} else {
 res = list(as.numeric(tclvalue(textEntryVarTcl)), "no")
 names(res) = c("Nb_raters", "Weighting")
 return(res)
}
}
