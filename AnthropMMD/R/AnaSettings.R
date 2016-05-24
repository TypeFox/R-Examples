AnaSettings <-
function(tab, type.data="raw_data") {
# tab = le tableau de donn\'ees
# type_data = raw_data ou summarized_data

.provenv = new.env()

tt <- tktoplevel(width=600, height=420) # on cr\'ee une fenetre de taille fix\'ee 
 tkpack.propagate(tt, FALSE) 
 tktitle(tt) <- "Analysis settings" # Le titre qui s'affichera dans la barre de fenetre 
 fontHeading <- tkfont.create(family="arial",size=12,weight="bold") # un style de fonte 

 frameGroupes <- tkframe(tt, relief="groove", borderwidth=3, width=550, height=130) # declaration d'une frame 
 frameProtocole <- tkframe(tt, relief="groove", borderwidth=3, width=550, height=150) # declaration d'une frame 
 frameOutputs <- tkframe(tt, relief="groove", borderwidth=3, width=550, height=90) # declaration d'une frame 

 tkpack(frameGroupes, frameProtocole, frameOutputs) 
 tkpack.propagate(frameGroupes, "0") 
 tkpack.propagate(frameProtocole, "0")
 tkpack.propagate(frameOutputs, "0")
 tkpack(tklabel(frameGroupes,text="Select the active groups (by default, all groups are active)", font=fontHeading))
 tkpack(tklabel(frameProtocole,text="Variable selection and other parameters", font=fontHeading))
 tkpack(tklabel(frameOutputs,text="Outputs", font=fontHeading))

########## FRAME GROUPE
frameListe <- tkframe(frameGroupes, relief="flat", borderwidth=1, width=300, height=50) # declaration d'une frame 
scr <- tkscrollbar(frameListe, repeatinterval=5, command=function(...)tkyview(maliste,...))
maliste<-tklistbox(frameListe, height=5, selectmode="multiple", yscrollcommand=function(...)tkset(scr,...),background="white")
tkpack(frameListe)
tkpack(maliste, scr, side="left")

# on remplit la liste
if (type.data=="raw_data") {
 groupes <- levels(factor(tab[,1]))
 for (i in (1:length(groupes))) { tkinsert(maliste, "end", groupes[i]) } 
} else if (type.data=="summarized_data") {
 groupes <- rownames(tab)[1:(nrow(tab)/2)]
 for (i in (1:length(groupes))) { tkinsert(maliste, "end", groupes[i]) } 
}
tkselection.set(maliste, 0, nlevels(factor(tab[,1]))-1)

######## FRAME PROTOCOLE
 frameAllTraits <- tkframe(frameProtocole, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 cb1 <- tkcheckbutton(frameAllTraits) # premi\`ere tickbox
 var_all <- tclVar("0")
 tkconfigure(cb1,variable=var_all)

 frameAllBOTraits <- tkframe(frameProtocole, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 cb2 <- tkcheckbutton(frameAllBOTraits) # seconde tickbox
 var_abo <- tclVar("0")
 tkconfigure(cb2,variable=var_abo)

 frameMinIndiv <- tkframe(frameProtocole, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 textEntryVarTcl <- tclVar("10")
 textEntryWidget <- tkentry(frameMinIndiv, width = 2, bg="white", textvariable = textEntryVarTcl)
 
 frameFTCor <- tkframe(frameProtocole, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 cb3 <- tkcheckbutton(frameFTCor) # troisi\`eme tickbox
 var_ftcor <- tclVar("1")
 tkconfigure(cb3,variable=var_ftcor)
 
 frameFormula <- tkframe(frameProtocole, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 rb1 <- tkradiobutton(frameFormula)
 rb2 <- tkradiobutton(frameFormula)
 var_formula <- tclVar("Anscombe") # le choix par d\'efaut
 tkconfigure(rb1,variable=var_formula,value="Anscombe", text="Anscombe formula ; or ")
 tkconfigure(rb2,variable=var_formula,value="FreemanTukey", text="Freeman and Tukey formula")

 tkpack(frameAllTraits)
 tkpack(frameAllBOTraits)
 tkpack(tklabel(frameAllTraits,text="Keep the traits with the same value for all the individuals of the dataset"), cb1, side="left")
 tkpack(tklabel(frameAllBOTraits,text="Keep the traits with the same value for all the individuals, except one    "), cb2, side="left")
 tkpack(frameMinIndiv)
 tkpack(tklabel(frameMinIndiv, text = "Keep only the traits with this minimal number of individuals per group : "), textEntryWidget, side="left")
 tkpack(frameFTCor)
 tkpack(tklabel(frameFTCor,text="Apply Freeman-Tukey correction for small samples                                 "), cb3, side="left")
 tkpack(frameFormula)
 tkpack(tklabel(frameFormula, text="Use :"), rb1, rb2, side="left")
 
 
######### FRAME OUTPUTS
 frameMDSGraph <- tkframe(frameOutputs, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 cb4 <- tkcheckbutton(frameMDSGraph) # premi\`ere tickbox
 var_mds <- tclVar("1")
 tkconfigure(cb4,variable=var_mds)

 frameDiscriPow <- tkframe(frameOutputs, relief="flat", borderwidth=1, width=540, height=50) # declaration d'une frame 
 cb5 <- tkcheckbutton(frameDiscriPow) # seconde tickbox
 var_dp <- tclVar("0")
 tkconfigure(cb5,variable=var_dp)

 tkpack(frameMDSGraph)
 tkpack(frameDiscriPow)
 tkpack(tklabel(frameMDSGraph,text="Classical Multidimensional scaling plot (MDS)                 "), cb4, side="left")
 tkpack(tklabel(frameDiscriPow,text="Save individual measure of divergence for each variable"), cb5, side="left")
 
 
##### OK BUTTON :
 OnOK <- function() {
   #groupes_retenus<<-as.numeric(tkcurselection(maliste))+1
   assign("groupes_retenus", as.numeric(tkcurselection(maliste))+1, envir=.provenv)
  tkdestroy(tt)
 }
 OK.but <- tkbutton(tt,text="OK",command=OnOK, relief="groove",borderwidth=3,width=7,bg="green4",fg="white")
 tkpack(OK.but)
 tkfocus(tt)
 tkwait.window(OK.but) # TRES IMPORTANT : l'ex\'ecution s'arrete ici tant que l'utilisateur n'a pas cliqu\'e sur OK 

 groupes_retenus = get("groupes_retenus", envir=.provenv)
 user_choice = list(groupes_retenus, as.logical(as.numeric(tclvalue(var_all))), as.logical(as.numeric(tclvalue(var_abo))), as.logical(as.numeric(tclvalue(var_ftcor))), as.numeric(tclvalue(textEntryVarTcl)), tclvalue(var_formula), as.logical(as.numeric(tclvalue(var_mds))), as.logical(as.numeric(tclvalue(var_dp))))
 names(user_choice) = c("Groupes_retenus", "All_vars", "All_BO_vars", "FTcor", "Nb_min_indiv", "Formula", "MDS", "IDP")
 return(user_choice)
}

