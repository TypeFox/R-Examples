see.logic<-function(){

local({
    have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
    if(have_ttk) {
        tkbutton <- ttkbutton
        tkframe <- ttkframe
        tklabel <- ttklabel
    }
   tclServiceMode(TRUE)
   top <- tktoplevel() 
   tkwm.geometry(top, "+500+4")
   tktitle(top) <- "Logical versus fallacious arguments"
    canvas <- tkcanvas(top, relief="raised", width=800, height=600)
    tkpack(canvas, side="top", fill="x")

    plotFont <- "Helvetica 9"
    plotFont2 <- "Calibri 11"
    
#------------------------ Message boxes ----------------------#

I.T<-function(){
 tkmessageBox(message="Correct.  This is a logical form of argument called \U0022modus ponens\U0022 \nor \U0022 affirming the antecedent\U0022.")}
    
I.F<-function(){
 tkmessageBox(message="Incorrect.  This is a logical form of argument called \U0022modus ponens\U0022 \nor \U0022 affirming the antecedent\U0022.", icon="error")}

II.T<-function(){
 tkmessageBox(message="Correct.  This is a logical form of argument called \U0022modus tollens\U0022 or\n\U0022 denying the consequent\U0022.")}
 
II.F<-function(){
 tkmessageBox(message="Incorrect.  This is a logical form of argument called \U0022modus tollens\U0022 or\n\U0022 denying the consequent.\U0022", icon="error")}

III.T<-function(){
 tkmessageBox(message="Correct.  This is another logical example of \U0022modus tollens\U0022 or\n\U0022 denying the consequent\U0022.")}
 
III.F<-function(){
 tkmessageBox(message="Incorrect.  This is another logical example of \U0022modus tollens\U0022 or\n\U0022 denying the consequent\U0022.", icon="error")}

IV.T<-function(){
 tkmessageBox(message="Incorrect.  This is an example of a fallacious type of argument called \n\U0022 affirming the consequent\U0022.", icon="error")}
 
IV.F<-function(){
 tkmessageBox(message="Correct.  This is an example of a fallacious type of argument called \n\U0022 affirming the consequent\U0022.")}
 
V.T<-function(){
 tkmessageBox(message="Incorrect.  This is another example of a fallacious type of argument \ncalled \U0022 affirming the consequent\U0022.", icon="error")}
 
V.F<-function(){
 tkmessageBox(message="Correct.  This is another example of a fallacious type of argument called \n\U0022 affirming the consequent\U0022.")} 
 
VI.T<-function(){
 tkmessageBox(message="Incorrect.  No causality has been established. This is an example of a fallacious type of argument called \U0022post hoc ergo propter hoc\U0022.  Translated this means: \U0022 after this; therefore, because of this\U0022.", icon="error")}
 
VI.F<-function(){
 tkmessageBox(message="Correct.  No causality has been established.  This is an example of a fallacious type of argument called \U0022post hoc ergo propter hoc\U0022.   Translated this means: \U0022 after this; therefore, because of this\U0022.")} 
 
VII.T<-function(){
 tkmessageBox(message="Incorrect.  This is a confusion of cause and effect.", icon="error")}
 
VII.F<-function(){
 tkmessageBox(message="Correct.  This is a confusion of cause and effect.")} 
 
VIII.T<-function(){
 tkmessageBox(message="Correct.  This is a logical form of an argument called \U0022reductio ad absurdum\U0022.  This is closely related to denying the consequent arguments.")}
 
VIII.F<-function(){
 tkmessageBox(message="Incorrect.  This is a logical form of an argument called \U0022reductio ad absurdum\U0022.  This is closely related to denying the consequent arguments.",icon="error")} 

#------------------------------------------------------------------------------#

    tkcreate(canvas, "polygon", 10, 10, 10, 590, 790, 590, 790, 10, width=0, fill = "white")
    tkcreate(canvas, "line", 10, 10, 790, 10, width=2)
    tkcreate(canvas, "line", 10, 590, 790, 590, width=2)
    tkcreate(canvas, "line", 10, 10, 10, 590, width=2)
    tkcreate(canvas, "line", 10, 300, 790, 300, width=2)
    tkcreate(canvas, "line", 10, 155, 790, 155, width=2)
    tkcreate(canvas, "line", 10, 445, 790, 445, width=2)
    tkcreate(canvas, "line", 400, 10, 400, 590, width=2)
    tkcreate(canvas, "line", 790, 10, 790, 590, width=2)
    
 #----------------------------- Argument 1 -------------------------#
    tkcreate(canvas, "text", 205, 20, text="Argument 1",
             font=plotFont2)
    tkcreate(canvas, "text", 20, 50, text="If Joe fails his ecology final he fails the class.   }Premise 1",
             font=plotFont,anchor="w")
    tkcreate(canvas, "text", 20, 75, text="Joe fails his ecology final.                                      }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 20, 85, 350, 85, width=1)
    tkcreate(canvas, "text", 20, 100, text="Joe fails the class.                                                   }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 125, 140, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 255, 140, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemI.T <- tkcreate(canvas, "oval", 110, 135, 120, 145,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointI.T", "withtag", itemI.T)
     itemI.F <- tkcreate(canvas, "oval", 240, 135, 250, 145,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointI.F", "withtag", itemI.F)       

#----------------------------- Argument 2 -------------------------#    

tkcreate(canvas, "text", 205, 165, text="Argument 2",
             font=plotFont2)
    tkcreate(canvas, "text", 20, 190, text="Only if Joe fails his ecology final does he fail the class.  }Premise 1",
             font=plotFont,anchor="w")
    tkcreate(canvas, "text", 20, 215, text="Joe does not fail his ecology final.                                        }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 20, 225, 350, 225, width=1)
    tkcreate(canvas, "text", 20, 240, text="Joe passes the class.                                                             }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 125, 280, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 255, 280, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemII.T <- tkcreate(canvas, "oval", 110, 275, 120, 285,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointII.T", "withtag", itemII.T)
     itemII.F <- tkcreate(canvas, "oval", 240, 275, 250, 285,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointII.F", "withtag", itemII.F)       
    
#----------------------------- Argument 3 -------------------------#    

 tkcreate(canvas, "text", 205, 310, text="Argument 3",
             font=plotFont2)
    tkcreate(canvas, "text", 20, 340, text="If A then not B.   }Premise 1",
             font=plotFont,anchor="w")
    tkcreate(canvas, "text", 20, 365, text="B.                         }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 20, 375, 350, 375, width=1)
    tkcreate(canvas, "text", 20, 390, text="not A.                   }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 125, 430, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 255, 430, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemIII.T <- tkcreate(canvas, "oval", 110, 425, 120, 435,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointIII.T", "withtag", itemIII.T)
     itemIII.F <- tkcreate(canvas, "oval", 240, 425, 250, 435,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointIII.F", "withtag", itemIII.F)
              
#----------------------------- Argument 4 -------------------------#    

tkcreate(canvas, "text", 205, 455, text="Argument 4",
             font=plotFont2)
    tkcreate(canvas, "text", 20, 485, text="If A then B.         }Premise 1",
             font=plotFont,anchor="w")                     
    tkcreate(canvas, "text", 20, 510, text="Not A.                  }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 20, 520, 350, 520, width=1)
    tkcreate(canvas, "text", 20, 535, text="Not B.                  }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 125, 575, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 255, 575, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemIV.T <- tkcreate(canvas, "oval", 110, 570, 120, 580,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointIV.T", "withtag", itemIV.T)
     itemIV.F <- tkcreate(canvas, "oval", 240, 570, 250, 580,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointIV.F", "withtag", itemIV.F)  


#----------------------------- Argument 5 -------------------------#
    tkcreate(canvas, "text", 595, 20, text="Argument 5",
             font=plotFont2)
    tkcreate(canvas, "text", 410, 50, text="If Pocatello is in ID, then it is in the Pacific NW.   }Premise 1",
             font=plotFont,anchor="w")
    tkcreate(canvas, "text", 410, 75, text="Pocatello is in the Pacific NW.                                 }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 410, 85, 740, 85, width=1)
    tkcreate(canvas, "text", 410, 100, text="Pocatello is in ID.                                                       }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 510, 140, text="Logically corect",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 640, 140, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemV.T <- tkcreate(canvas, "oval", 495, 135, 505, 145,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointV.T", "withtag", itemV.T)
     itemV.F <- tkcreate(canvas, "oval", 625, 135, 635, 145,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointV.F", "withtag", itemV.F)     

#----------------------------- Argument 6 -------------------------#    

tkcreate(canvas, "text", 595, 165, text="Argument 6",
             font=plotFont2)
    tkcreate(canvas, "text", 410, 190, text="Sarah ate vitamins during her pregnancy.   }Premise",
             font=plotFont,anchor="w")
    tkcreate(canvas, "line", 410, 200, 740, 200, width=1)
    tkcreate(canvas, "text", 410, 215, text="Sarah's baby will be a genius.                       }Conclusion",
             font=plotFont, anchor="w")
       
    tkcreate(canvas, "text", 510, 280, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 640, 280, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemVI.T <- tkcreate(canvas, "oval", 495, 275, 505, 285,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVI.T", "withtag", itemVI.T)
     itemVI.F <- tkcreate(canvas, "oval", 625, 275, 635, 285,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVI.F", "withtag", itemVI.F) 

#----------------------------- Argument 7 -------------------------#    

 tkcreate(canvas, "text", 595, 310, text="Argument 7",
             font=plotFont2)
    tkcreate(canvas, "text", 410, 340, text="Sober, industrious students own at least one cow.       }Premise",
             font=plotFont,anchor="w")
    tkcreate(canvas, "line", 410, 350, 740, 350, width=1)
    tkcreate(canvas, "text", 410, 365, text="Drunk, lazy students should be given a cow \nto make them sober and industrious.",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 705, 370, text="}Conclusion",
             font=plotFont, anchor="w")   
    tkcreate(canvas, "text", 510, 430, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 640, 430, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemVII.T <- tkcreate(canvas, "oval", 495, 425, 505, 435,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVII.T", "withtag", itemVII.T)
     itemVII.F <- tkcreate(canvas, "oval", 625, 425, 635, 435,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVII.F", "withtag", itemVII.F)

#----------------------------- Argument 8 -------------------------#    

tkcreate(canvas, "text", 595, 455, text="Argument 8",
             font=plotFont2)
    tkcreate(canvas, "text", 410, 485, text="To prove A, we assume not A.         }Premise 1",
             font=plotFont,anchor="w")                     
    tkcreate(canvas, "text", 410, 510, text="We disprove not A                              }Premise 2",
             font=plotFont, anchor="w")
    tkcreate(canvas, "line", 410, 520, 740, 520, width=1)
    tkcreate(canvas, "text", 410, 535, text="conclude A.                                          }Conclusion",
             font=plotFont, anchor="w")
    tkcreate(canvas, "text", 510, 575, text="Logically correct",
             font=plotFont2, anchor="w")
    tkcreate(canvas, "text", 640, 575, text="Fallacious",
             font=plotFont2, anchor="w")         
     itemVIII.T <- tkcreate(canvas, "oval", 495, 570, 505, 580,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVIII.T", "withtag", itemVIII.T)
     itemVIII.F <- tkcreate(canvas, "oval", 625, 570, 635, 580,
                                 width=1, outline="black",
                                 fill="SkyBlue2")        
            tkaddtag(canvas, "pointVIII.F", "withtag", itemVIII.F)  


#------------------------ Drive widgets ----------------------#             

        tkitembind(canvas, "pointI.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointI.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointI.T", "<1>", substitute(I.T()))


        tkitembind(canvas, "pointI.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointI.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointI.F", "<1>", substitute(I.F()))
        
#---------------------------------------------------------------#        
        
          tkitembind(canvas, "pointII.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointII.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointII.T", "<1>", substitute(II.T()))


        tkitembind(canvas, "pointII.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointII.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointII.F", "<1>", substitute(II.F()))
        
#---------------------------------------------------------------#         

        
          tkitembind(canvas, "pointIII.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointIII.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointIII.T", "<1>", substitute(III.T()))


        tkitembind(canvas, "pointIII.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointIII.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointIII.F", "<1>", substitute(III.F()))
        
#---------------------------------------------------------------#  
        
          tkitembind(canvas, "pointIV.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointIV.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointIV.T", "<1>", substitute(IV.T()))


        tkitembind(canvas, "pointIV.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointIV.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointIV.F", "<1>", substitute(IV.F()))
        
#---------------------------------------------------------------#  

   tkitembind(canvas, "pointV.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointV.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointV.T", "<1>", substitute(V.T()))


        tkitembind(canvas, "pointV.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointV.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointV.F", "<1>", substitute(V.F()))
        
#---------------------------------------------------------------#  


   tkitembind(canvas, "pointVI.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVI.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVI.T", "<1>", substitute(VI.T()))


        tkitembind(canvas, "pointVI.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVI.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVI.F", "<1>", substitute(VI.F()))
        
#---------------------------------------------------------------#  

 
   tkitembind(canvas, "pointVII.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVII.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVII.T", "<1>", substitute(VII.T()))


        tkitembind(canvas, "pointVII.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVII.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVII.F", "<1>", substitute(VII.F()))
        
#---------------------------------------------------------------#  

 tkitembind(canvas, "pointVIII.T", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVIII.T", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVIII.T", "<1>", substitute(VIII.T()))


        tkitembind(canvas, "pointVIII.F", "<Any-Enter>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="red"))
        tkitembind(canvas, "pointVIII.F", "<Any-Leave>",
                   function() tkitemconfigure(canvas, "current",
                                              fill="SkyBlue2"))
        tkitembind(canvas, "pointVIII.F", "<1>", substitute(VIII.F()))
        
#---------------------------------------------------------------#  
tkpack(tkbutton(top, text = "Exit", command = function()tkdestroy(top)),side="right")
})
}


    