"guiOutputTask4" <-
function(K,confidenceLevel,secondTimeScaleIsUsed,t,t2,t2max,lowerBounds,upperBounds,BoundsSymmetry,
         enterBoundsManually,alpha,phi,confidenceIntervall,spendingFunctionUsed,Zvalue, taskWindow)
{
  # Initializing
  FunctionNames=NULL;
  FunctionNamesUpper=NULL;
  FunctionNamesLower=NULL;

  #Set Toplevel
  outTask4Toplevel <- tktoplevel(taskWindow)
  tkwm.title(outTask4Toplevel,paste("-4-    K =",K,", Z-value =",Zvalue))

  #Define main Frame
  OutputTask4 <- tkframe(outTask4Toplevel)

  #Define subframes
  parametersFrame <- tkframe(OutputTask4,relief="groove",borderwidth=0)

  #create labels with parameter values:
  tkgrid( tklabel(parametersFrame, text=paste("K=",K)),sticky="w")

  if(!enterBoundsManually)
  {
    if(!BoundsSymmetry==3)
    {
      tkgrid( tklabel(parametersFrame, text=paste("alpha =",alpha[1])),sticky="w")
    }
    else
    {
      tkgrid( tklabel(parametersFrame, text=paste("alpha - Upper Bounds =",alpha[1])),sticky="w")
      tkgrid( tklabel(parametersFrame, text=paste("alpha - Lower Bounds =",alpha[2])),sticky="w")

    }
  }

  if(!BoundsSymmetry==3)
  {
    if(enterBoundsManually)
    {
      tkgrid( tklabel(parametersFrame, text="manually entered Bounds"),sticky="w")
    }
    else
    {
      ##names of spending functions that could have been used
      FunctionNames <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[1],sep=""),
                       paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")

      # substitute according funtion in output
      tkgrid( tklabel(parametersFrame, text=paste("Function: ",FunctionNames[[spendingFunctionUsed[1]]])),sticky="w")
    }
  }
  else
  {
    if(enterBoundsManually)
    {
      tkgrid( tklabel(parametersFrame, text="Upper Bounds: manually entered "),sticky="w")
      tkgrid( tklabel(parametersFrame, text="Lower Bounds: manually entered "),sticky="w")
    }
    else
    {
      ##names of spending functions that could have been used
      FunctionNamesUpper <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[1],sep=""),
                            paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")
      FunctionNamesLower <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[2],sep=""),
                            paste("Hwang-Shih-DeCani Family ( phi =",phi[2],")"),"Exact Pocock Bounds")

      # substitute according funtion in output
      tkgrid( tklabel(parametersFrame, text=paste("Function - Upper Bounds: ",FunctionNamesUpper[spendingFunctionUsed[1]])),sticky="w")
      tkgrid( tklabel(parametersFrame, text=paste("Function - Lower Bounds: ",FunctionNamesLower[spendingFunctionUsed[2]])),sticky="w")
    }
  }

  #confidence- level and intervall
  tkgrid( tklabel(parametersFrame, text=""),sticky="w") #blank line
  tkgrid( tklabel(parametersFrame, text=paste("Z-value at last analysis= ",Zvalue)),sticky="w")
  tkgrid( tklabel(parametersFrame, text=paste("Confidence Level= ",confidenceLevel,"%")),sticky="w")
  tkgrid( tklabel(parametersFrame, text=""),sticky="w") #blank line
  tkgrid( tklabel(parametersFrame, text=paste("Confidence Intervall= <",round(confidenceIntervall[1],digits=5)," , ",
                                                                      round(confidenceIntervall[2],digits=5),">")),sticky="w")
  tkgrid( tklabel(parametersFrame, text="Drift is equal to the expectation of the Z statistic when time=1."),sticky="w")
  if(t2max!=0)
  {
    tkgrid( tklabel(parametersFrame, text=paste("Maximum Information=",t2max)),sticky="w")
  }

  tkgrid( tklabel(parametersFrame, text=""),sticky="w") #blank line
  tkgrid(parametersFrame,sticky="w")


  ###########################################################################
  # function handles click onto button to show results of bounds in a graph #
  ###########################################################################
  onShowGraph <- function()
  {
    if(enterBoundsManually)
    {
      ## if one-Sided-Test we won't see negative Z-Values
      if(BoundsSymmetry==1)
      {
        xCoordinate<-t
        yCoordinate<-upperBounds

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main=paste("-3-  K=",K,", Z=",Zvalue,
             ", confidence level = ",round(confidenceLevel,digits=5), ", confidence intervall= < ",
             round(confidenceIntervall[1],digits=5),",",round(confidenceIntervall[2],digits=5)," >",
             "\n Bounds manually entered", sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(0,4))

        ##...then add lines between them
        lines(t,upperBounds,col="blue")
      }

      else
      {
        xCoordinate<-c(t,t)
        yCoordinate<-c(lowerBounds,upperBounds)

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main=paste("-3-  K=",K,", Z=",Zvalue,
             ", confidence level = ",round(confidenceLevel,digits=5), ", confidence intervall= < ",
             round(confidenceIntervall[1],digits=5),",",round(confidenceIntervall[2],digits=5)," >",
             "\n Bounds manually entered", sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(-4,4))

        ##...then add lines between them
        lines(t,lowerBounds,col="blue")
        lines(t,upperBounds,col="blue")
      }
    }# endif 'if(enterBoundsManually)'

    else # spending function was used
    {
    ## if one-Sided-Test we won't see negative Z-Values
      if(BoundsSymmetry==1)
      {
        xCoordinate<-t
        yCoordinate<-upperBounds

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main=paste("-3-  K=",K,", Z=",Zvalue,
             ", confidence level = ",round(confidenceLevel,digits=5), ", confidence intervall= < ",
             round(confidenceIntervall[1],digits=5),",",round(confidenceIntervall[2],digits=5)," >",
             "\n Function:", FunctionNames[spendingFunctionUsed[1]],", alpha=",round(alpha[1],digits=5), sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(0,4))

        ##...then add lines between them
        lines(t,upperBounds,col="blue")
      }

      else
      {
        xCoordinate<-c(t,t)
        yCoordinate<-c(lowerBounds,upperBounds)

        if(BoundsSymmetry==2)
        {
          ## first plotting bounds as points...
         plot(xCoordinate,yCoordinate,main=paste("-3-  K=",K,", Z=",Zvalue,
             ", confidence level = ",round(confidenceLevel,digits=5), ", confidence intervall= < ",
             round(confidenceIntervall[1],digits=5),",",round(confidenceIntervall[2],digits=5)," >",
             "\n Function:", FunctionNames[spendingFunctionUsed[1]],", alpha=",round(alpha[1],digits=5), sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(-4,4))
        }
        else
        {
         ## first plotting bounds as points...
         plot(xCoordinate,yCoordinate,main=paste("-3-  K=",K,", Z=",Zvalue,
             ", confidence level = ",round(confidenceLevel,digits=5), ", confidence intervall= < ",
             round(confidenceIntervall[1],digits=5),",",round(confidenceIntervall[2],digits=5)," >",
             "\n upper Function:", FunctionNamesUpper[spendingFunctionUsed[1]],", alpha=",round(alpha[1],digits=5),
             "\n lower Function:", FunctionNamesLower[spendingFunctionUsed[2]],", alpha=",round(alpha[2],digits=5),sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(-4,4))
        }

        ##...then add lines between them
        lines(t,lowerBounds,col="blue")
        lines(t,upperBounds,col="blue")
      }
    }
  }

  ##################################################################
  ## function handles click onto button to save results in a file ##
  ##################################################################
  onSave <- function()
  {
     #create file variable
     fileName <- tclvalue(tkgetSaveFile(initialfile=".html",filetypes="{{html Files} {.html}} {{All files} *}"))
     if (fileName=="") return;

     #open file
     zz <- file(fileName,"w")

     #output will be writed in HTML
     cat("<html> <body> \n",file = zz)

     #output K
     cat("K=",K,"<br> \n",file = zz)

     ##ouput alpha
     if(!enterBoundsManually)
     {
       if(!BoundsSymmetry==3)
       {
         cat("&alpha; =",alpha[1],"<br>\n",file = zz)
       }
       else
       {
         cat("Upper &alpha; = ",alpha[1],"<br>\n",file = zz)
         cat("Lower &alpha; = ",alpha[2],"<br>\n",file = zz)
       }
     }
     if(!BoundsSymmetry==3)
     {
       if(enterBoundsManually)
       {
         cat("<b>manually entered</b> Bounds <br>\n",file = zz)
       }
       else
       {
         ##output names of spending functions that were used
         FunctionNames <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: &alpha;&sdot;t<sup>",phi[1],"</sup>"),
                          paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")

         # substitute according funtion in output
         cat("<b>",FunctionNames[spendingFunctionUsed[1]],"</b>"," was used as spending Function.","<br>\n",file = zz)
       }
     }
     else
     {
       if(enterBoundsManually)
       {
         cat("UPPER Bounds: <b>manually entered</b> <br>\n",file = zz)
         cat("LOWER Bounds: <b>manually entered</b> <br>\n",file = zz)
       }
       else
       {
         ##names of spending functions that could have been used
         FunctionNamesUpper <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: &alpha;&sdot;t<sup>",phi[1],"</sup>"),
                          paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")
         FunctionNamesLower <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: &alpha;&sdot;t<sup>",phi[2],"</sup>"),
                          paste("Hwang-Shih-DeCani Family ( phi =",phi[2],")"),"Exact Pocock Bounds")

         # substitute according funtion in output
         cat("Spending Function for UPPER Bound:","<b>",FunctionNamesUpper[spendingFunctionUsed[1]],"</b>","<br>\n",file = zz)
         cat("Spending Function for LOWER Bound:","<b>",FunctionNamesLower[spendingFunctionUsed[2]],"</b>","<br>\n",file = zz)
       }
     }

     #confidence- level and intervall
     cat("<br>Confidence Level= ",confidenceLevel,"% <br>\n",file = zz)
     cat("Confidence Intervall= <",round(confidenceIntervall[1],digits=5)," , ",
                                   round(confidenceIntervall[2],digits=5),"> <br>\n",file = zz)
     cat("Drift is equal to the expectation of the Z statistic when time=1.<br> \n",file=zz)
     if(t2max!=0)
     {
       cat("Maximum Information=",t2max,"<br> \n",file=zz)
     }

     cat("</body> </html> \n",file = zz)
     close(zz)
  }


  ###########################################################################
  ##function handles click onto button to Cancel i.e. close current window ##
  ###########################################################################
  onCancel <- function()
  {
   tkdestroy(outTask4Toplevel)
  }

  #frame for the buttons
  buttonFrame<-tkframe(OutputTask4,relief="groove",borderwidth=0)

  #button to show graphic
  showGraph.button <-tkbutton(buttonFrame,text="  Show Graph  ",command=onShowGraph)

  #button to save in file
  save.button <-tkbutton(buttonFrame,text="  Save to File  ",command=onSave)

  #button to cancel i.e. close current window
  cancel.button <-tkbutton(buttonFrame,text="  Cancel   ",command=onCancel)

  #grid buttons
  tkgrid( tklabel(buttonFrame, text=""))   #blank line
  tkgrid(showGraph.button,tklabel(buttonFrame, text="      "),
         save.button, tklabel(buttonFrame,     text="      "),
         cancel.button, sticky="we")
  tkgrid(buttonFrame)
  tkgrid( tklabel(buttonFrame, text=""))   #blank line

  #grid allover frame and focus
  tkgrid(OutputTask4,sticky="w")
  tkfocus(outTask4Toplevel)



}#end <--*function(...)*
