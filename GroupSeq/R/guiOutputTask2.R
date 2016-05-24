"guiOutputTask2" <-
function(K,probTotal,drift,expectedStoppingTime,secondTimeScaleIsUsed,t,t2,t2max,
         lowerBounds,upperBounds,probStopping,probExceedingUpper,probExceedingLower, confidenceLevel,
         BoundsSymmetry,enterBoundsManually,alpha,phi,spendingFunctionUsed, taskWindow)

{
  ###INITIALiZE VARIABLES###
  resultExitProb<-0 #exit probability
  cumulativeExitProb<-0 #cumulative exit probability
  FunctionNames=NULL;
  FunctionNamesUpper=NULL;
  FunctionNamesLower=NULL;

  ##compute exit probability and cumulative exit probability
  for(i in 1:K)
  {
    resultExitProb[i] <- probExceedingUpper[i]+probExceedingLower[i]

    if(i==1)
    {
      cumulativeExitProb[i] <- resultExitProb[i]
    }
    else
    {
      cumulativeExitProb[i] <- cumulativeExitProb[i-1] + resultExitProb[i]
    }
  }


  #Set Toplevel
  outTask2Toplevel <- tktoplevel(taskWindow)
  tkwm.title(outTask2Toplevel,paste("-2-   K =",K,", drift =",round(drift,digits=5),", Power=",confidenceLevel))

  #Define main Frame
  OutputTask2 <- tkframe(outTask2Toplevel)

  #Define subframes
  staticFrame <- tkframe(OutputTask2,relief="groove",borderwidth=0)
  dynamicFrame <- tkframe(OutputTask2,relief="groove",borderwidth=0)
  parametersFrame <- tkframe(staticFrame,relief="groove",borderwidth=0)
  numbersFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  timesFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  secondTimesFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  lowerBoundsFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  upperBoundsFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  exitProbabilityFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  cumulativeExitProbFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)


  #create label with parameter values:
  tkgrid( tklabel(parametersFrame, text=paste("K =",K)),sticky="w")
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
      tkgrid( tklabel(parametersFrame, text=paste("Function: ",FunctionNames[spendingFunctionUsed[1]])),sticky="w")
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
      ## names of spending functions that could have been used
      FunctionNamesUpper <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[1],sep=""),
                              paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")
      FunctionNamesLower <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[2],sep=""),
                              paste("Hwang-Shih-DeCani Family ( phi =",phi[2],")"),"Exact Pocock Bounds")

      # substitute according funtion in output
      tkgrid( tklabel(parametersFrame, text=paste("Function - Upper Bounds: ",FunctionNamesUpper[spendingFunctionUsed[1]])),sticky="w")
      tkgrid( tklabel(parametersFrame, text=paste("Function - Lower Bounds: ",FunctionNamesLower[spendingFunctionUsed[2]])),sticky="w")

    }
  }

  #drift parameters
  tkgrid( tklabel(parametersFrame, text=""),sticky="w") #blank line
  tkgrid( tklabel(parametersFrame, text=paste("drift=",round(drift,digits=5))),sticky="w")
  if(t2max!=0)
  {
    tkgrid( tklabel(parametersFrame, text=paste("Maximum Information =",t2max)),sticky="w")
  }
  tkgrid( tklabel(parametersFrame, text=paste("Power =",confidenceLevel)),sticky="w")
  tkgrid( tklabel(parametersFrame, text=""),sticky="w") #blank line


  #create head labels
  tkgrid( tklabel(numbersFrame, text="k    "),sticky="w")
  tkgrid( tklabel(timesFrame, text="Times   "),sticky="w")
  tkgrid( tklabel(secondTimesFrame, text="2nd Time Scale   "),sticky="w")
  tkgrid( tklabel(lowerBoundsFrame, text="Lower Bounds  "),sticky="w")
  tkgrid( tklabel(upperBoundsFrame, text="Upper Bounds  "),sticky="w")
  tkgrid( tklabel(exitProbabilityFrame, text="Exit Probability  "),sticky="w")
  tkgrid( tklabel(cumulativeExitProbFrame, text="Cumulative Exit Prob.  "),sticky="w")

  #create labels with results
  for(i in 1:K)
  {
    tkgrid( tklabel(numbersFrame, text=as.character(i)),sticky="w")
    tkgrid( tklabel(timesFrame, text=as.character(round(t[i],digits=3))),sticky="w")
    tkgrid( tklabel(secondTimesFrame, text=as.character(round(t2[i],digits=3))),sticky="w")
    tkgrid( tklabel(lowerBoundsFrame, text=as.character(round(lowerBounds[i],digits=4))),sticky="w")
    tkgrid( tklabel(upperBoundsFrame, text=as.character(round(upperBounds[i],digits=4))),sticky="w")
    tkgrid( tklabel(exitProbabilityFrame, text=as.character(round(resultExitProb[i],digits=10))),sticky="w")
    tkgrid( tklabel(cumulativeExitProbFrame, text=as.character(round(cumulativeExitProb[i],digits=10))),sticky="w")
  }

  #put frames together
  tkgrid(parametersFrame,sticky="w")
  ##if second information time is used output with information time
  if(secondTimeScaleIsUsed)
  {
    tkgrid(numbersFrame, timesFrame,secondTimesFrame, lowerBoundsFrame, upperBoundsFrame, exitProbabilityFrame, cumulativeExitProbFrame, sticky="w")
  }
  else
  {
    tkgrid(numbersFrame, timesFrame,                  lowerBoundsFrame, upperBoundsFrame, exitProbabilityFrame, cumulativeExitProbFrame, sticky="w")
  }
  tkgrid(staticFrame,sticky="w")
  tkgrid(dynamicFrame,sticky="w")


  ###########################################################################
 ## function handles click onto button to show results of bounds in a graph ##
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
        plot(xCoordinate,yCoordinate,main=paste("-2-  K=",K,", drift=",round(drift,digits=5),
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
        plot(xCoordinate,yCoordinate,main=paste("-2-  K=",K,", drift=",round(drift,digits=5),
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
        plot(xCoordinate,yCoordinate,main=paste("-2-  K=",K,", drift=",round(drift,digits=5),
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
           plot(xCoordinate,yCoordinate,main=paste("-2-  K=",K,", drift=",round(drift,digits=5),
           "\n Function:", FunctionNames[spendingFunctionUsed[1]],", alpha=",round(alpha[1],digits=5), sep=""),
           pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
           xlab="Times",ylab="Standardized Z-Value",ylim=c(-4,4))
         }
         else
         {
           ## first plotting bounds as points...
           plot(xCoordinate,yCoordinate,main=paste("-2-  K=",K,", drift=",round(drift,digits=5),
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

     #drift parameters
     cat("drift=",round(drift,digits=5),"<br> \n",file=zz)
     if(t2max!=0)
     {
       cat("Maximum Information=",t2max,"<br> \n",file=zz)
     }
     cat("Power=",confidenceLevel,"<br><br> \n",file=zz)


     ##output the bounds
     #labels
     cat("<br>\n",file = zz)
     cat("<table border=\"3\"> \n",file = zz)
     cat("<tr> \n",file = zz)
     cat("<td>Times &#160</td> \n" ,file=zz)
     if(secondTimeScaleIsUsed) { cat("<td>2nd Time Scale &#160</td> \n" ,file=zz) }
     cat("<td>Lower Bounds &#160</td>  <td>Upper Bounds &#160</td> \n",file = zz)
     cat("<td>Exit Probability &#160</td>  <td>Cumulative Exit Prob. &#160</td> \n",file = zz)
     cat("</tr> \n",file = zz)

     #values
     for(i in 1:K)
     {
       cat("<tr> \n",file = zz)
       cat("<td>",round(t[i],digits=3),"</td> \n",file = zz )
       if(secondTimeScaleIsUsed) { cat("<td>",round(t2[i],digits=3),"</td> \n",file = zz ) }
       cat("<td>",round(lowerBounds[i],digits=4),"</td>",   "<td>",round(upperBounds[i],digits=4),"</td>",
           "<td>",round(resultExitProb[i],digits=10),"</td>",   "<td>",round(cumulativeExitProb[i],digits=10),"</td> \n",file = zz )
       cat("</tr> \n",file = zz)
     }

     cat("</table> \n",file = zz)
     cat("</body> </html> \n",file = zz)
     close(zz)
  }

  ###########################################################################
  ##function handles click onto button to Cancel i.e. close current window ##
  ###########################################################################
  onCancel <- function()
  {
   tkdestroy(outTask2Toplevel)
  }

  #frame for the buttons
  buttonFrame<-tkframe(OutputTask2,relief="groove",borderwidth=0)

  #button to show graphic
  showGraph.button <-tkbutton(buttonFrame,text="  Show Graph  ",command=onShowGraph)

  #button to save in file
  save.button <-tkbutton(buttonFrame,text="  Save to File  ",command=onSave)

  #button to cancel i.e. close current window
  cancel.button <-tkbutton(buttonFrame,text="  Cancel   ",command=onCancel)

  #grid buttons
  tkgrid( tklabel(buttonFrame, text=""))   #blank line
  tkgrid(showGraph.button,tklabel(buttonFrame, text="            "),
         save.button, tklabel(buttonFrame, text="            "),
         cancel.button, sticky="we")
  tkgrid(buttonFrame)
  tkgrid( tklabel(buttonFrame, text=""))   #blank line

  #grid allover frame and focus
  tkgrid(OutputTask2,sticky="w")
  tkfocus(outTask2Toplevel)


}#end <--*function(...)*
