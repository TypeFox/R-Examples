`selectgraph` <-
function (listnomvar, listgraph) 
{
####################################################
# création de la première boite de dialogue
####################################################

  choixgraphfunc <- function(title, question, liste) 
    {
      tt1 <- tktoplevel()
      tkfocus(tt1)
      tkwm.title(tt1, title)
      graphChoice <- ""
      lstbox <- tklistbox(tt1, height = 3, selectmode = "single",background = "white")
      tkgrid(tklabel(tt1, text = question))
      tkgrid(lstbox)
      graph <- liste

      for (i in (1:3)) 
       {
        tkinsert(lstbox, "end", graph[i])
       }

      OnOK <- function() 
       {
        graphChoice <<- graph[as.numeric(tkcurselection(lstbox)) +1]
        tkdestroy(tt1)
       }

      OK.but <- tkbutton(tt1, text = "   OK   ", command = OnOK)
      tkgrid(OK.but)
   
      tkbind(tt1, "<Destroy>", function() {tkdestroy(tt1)})
      tkfocus(tt1)
      tkwait.window(tt1)
      return(graphChoice)
   }


####################################################
# création de la seconde boite de dialogue
####################################################

choixvarfunc <- function(title, question, liste) 
   {
     tt2 <- tktoplevel()
     scr <- tkscrollbar(tt2, repeatinterval = 5, command = function(...) tkyview(lstbox,...))
     lstbox <- tklistbox(tt2, height = 4, selectmode = "single",yscrollcommand = function(...) tkset(scr, ...), background = "white")
     tkgrid(tklabel(tt2, text = question))
     tkfocus(tt2)
     tkwm.title(tt2, title)
     varChoice <- ""
     tkgrid(lstbox, scr)
     tkgrid.configure(scr, rowspan = 4, sticky = "nsw")
     var <- liste
       
      for (i in 1:length(var)) 
       {
        tkinsert(lstbox, "end", var[i])
       }

     OnOK <- function() 
      {
       varChoice <<- var[as.numeric(tkcurselection(lstbox)) + 1]
       tkdestroy(tt2)
      }

    OK.but <- tkbutton(tt2, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
       
    tkbind(tt2, "<Destroy>", function() {tkdestroy(tt2)})
    tkfocus(tt2)
    tkwait.window(tt2)
    return(varChoice)
   }



####################################################
# programme principal
####################################################


 graphChoice <- choixgraphfunc("Choice of graphic", "Kind of graphic desired",listgraph)

  if (length(graphChoice) == 0) 
   {
    varChoice1 <- ""
    varChoice2 <- ""
    }
   else 
    {
      if (graphChoice == "Histogram") 
       {
        varChoice1 <- choixvarfunc("Choice of variable","Which numeric variable do you want to choose ?",listnomvar)
         if (length(varChoice1) == 0)  
          {
           graphChoice <- ""
          }
         varChoice2 <- ""
        }
       else if (graphChoice == "Barplot") 
           {
            varChoice1 <- choixvarfunc("Choice of variable","Which factor do you want to choose ?",listnomvar)
            if (length(varChoice1) == 0)  
             {
              graphChoice <- ""
             }
            varChoice2 <- ""
           }

          else
          {
           varChoice1 <- choixvarfunc("Choice of 1st variable","Which numeric variable do you want in x-absciss ?",listnomvar)
            if (length(varChoice1) == 0)  
             {
              graphChoice <- ""
              varChoice2 <- ""
             }
            else 
             {
              varChoice2 <- choixvarfunc("Choice of 2nd variable","Which numeric variable do you want in y-abcsis ?",listnomvar)
                if (length(varChoice2) == 0)    
                {
                  graphChoice <- ""
                  varChoice1 <- ""
                }
             } 
         }
    }
    return(list(varChoice1=varChoice1,varChoice2=varChoice2,graphChoice=graphChoice))
}

