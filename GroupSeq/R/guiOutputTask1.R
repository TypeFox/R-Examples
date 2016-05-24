"guiOutputTask1" <-
function(K,alpha,phi,t,lowerBounds,upperBounds,probDifference,probExit,
         BoundsSymmetry,spendingFunctionUsed, taskWindow)
{
  # Initializing
  FunctionNames=NULL;
  FunctionNamesUpper=NULL;
  FunctionNamesLower=NULL;

  #Set Toplevel
  outTask1Toplevel <- tktoplevel(taskWindow)

  if(!BoundsSymmetry==3)
  {
    tkwm.title(outTask1Toplevel,paste("-1-    K=",K,  ", alpha=",alpha[1]))
  }
  else
  {
    tkwm.title(outTask1Toplevel,paste("-1-    ,K=",K,  ", alphaUPPER=",alpha[1],", alphaLOWER=",alpha[2]))
  }

  #Define main Frame
  OutputTask1 <- tkframe(outTask1Toplevel)

  #Define subframes
  staticFrame <- tkframe(OutputTask1,relief="groove",borderwidth=0)
  dynamicFrame <- tkframe(OutputTask1,relief="groove",borderwidth=0)
  parametersFrame <- tkframe(staticFrame,relief="groove",borderwidth=0)
  numbersFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  timesFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  lowerBoundsFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  upperBoundsFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  probDiffFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)
  probExitFrame <- tkframe(dynamicFrame,relief="groove",borderwidth=2)

  #create label with parameter values:
  tkgrid( tklabel(parametersFrame, text=paste("K =",K)),sticky="w")
  if(!BoundsSymmetry==3)
  {
    tkgrid( tklabel(parametersFrame, text=paste("alpha =",alpha[1])),sticky="w")
  }
  else
  {
    tkgrid( tklabel(parametersFrame, text=paste("alpha - Upper Bounds =",alpha[1])),sticky="w")
    tkgrid( tklabel(parametersFrame, text=paste("alpha - Lower Bounds =",alpha[2])),sticky="w")

  }

  if(!BoundsSymmetry==3)
  {
    ##names of spending functions that could have been used
    FunctionNames <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: alpha*t^",phi[1],sep=""),
                       paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")

    # substitute according funtion in output
    tkgrid( tklabel(parametersFrame, text=paste("Function: ",FunctionNames[[spendingFunctionUsed[1]]])),sticky="w")
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

  #create head labels
  tkgrid( tklabel(numbersFrame, text="k    "),sticky="w")
  tkgrid( tklabel(timesFrame, text="Times   "),sticky="w")
  tkgrid( tklabel(lowerBoundsFrame, text="Lower Bounds  "),sticky="w")
  tkgrid( tklabel(upperBoundsFrame, text="Upper Bounds  "),sticky="w")
  tkgrid( tklabel(probDiffFrame, text="alpha[i]-alpha[i-1]  "),sticky="w")
  tkgrid( tklabel(probExitFrame, text="cumulative alpha  "),sticky="w")

  #create labels with results
  for(i in 1:K)
  {
    tkgrid( tklabel(numbersFrame, text=as.character(i)),sticky="w")
    tkgrid( tklabel(timesFrame, text=as.character(round(t[i],digits=3))),sticky="w")
    tkgrid( tklabel(lowerBoundsFrame, text=as.character(round(lowerBounds[i],digits=4))),sticky="w")
    tkgrid( tklabel(upperBoundsFrame, text=as.character(round(upperBounds[i],digits=4))),sticky="w")
    tkgrid( tklabel(probDiffFrame, text=as.character(round(probDifference[i],digits=10))),sticky="w")
    tkgrid( tklabel(probExitFrame, text=as.character(round(probExit[i],digits=10))),sticky="w")
  }
  tkgrid( tklabel(dynamicFrame, text=""))   #blank line

  #put frames together
  tkgrid(parametersFrame,sticky="w")
  tkgrid(numbersFrame, timesFrame, lowerBoundsFrame, upperBoundsFrame, probDiffFrame, probExitFrame, sticky="w")
  tkgrid(staticFrame,sticky="w")
  tkgrid(dynamicFrame,sticky="w")


  ###########################################################################
  ##function handles click onto button to show results of bounds in a graph##
  ###########################################################################
  onShowGraph <- function()
  {
    ## if one-Sided-Test we won't see negative Z-Values
      if(BoundsSymmetry==1)
      {
        xCoordinate<-t
        yCoordinate<-upperBounds

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main=paste("-1-  K=",K,
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
          plot(xCoordinate,yCoordinate,main=paste("-1-  K=",K,
             "\n Function:", FunctionNames[spendingFunctionUsed[1]],", alpha=",round(alpha[1],digits=5), sep=""),
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4, cex.main=0.9,
             xlab="Times",ylab="Standardized Z-Value",ylim=c(-4,4))
        }
        else
        {
        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main=paste("-1-  K=",K,
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

     #output will be written in HTML
     cat("<html> <body> \n",file = zz)

     #output K
     cat("K=",K,"<br> \n",file = zz)

     ##ouput alpha
     if(!BoundsSymmetry==3)
     {
       cat("&alpha; =",alpha[1],"<br>\n",file = zz)
     }
     else
     {
       cat("Upper &alpha; = ",alpha[1],"<br>\n",file = zz)
       cat("Lower &alpha; = ",alpha[2],"<br>\n",file = zz)
     }

     if(BoundsSymmetry!=3)
     {
       ##output names of spending functions that were used
       FunctionNames <- c("O'Brien-Fleming Type","Pocock Type",paste("Power Family: &alpha;&sdot;t<sup>",phi[1],"</sup>"),
                        paste("Hwang-Shih-DeCani Family ( phi =",phi[1],")"),"Exact Pocock Bounds")

       # substitute according funtion in output
       cat("<b>",FunctionNames[spendingFunctionUsed[1]],"</b>"," was used as spending Function.","<br>\n",file = zz)
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

     ##output the bounds
     cat("<br>\n",file = zz)
     cat("<table border=\"3\"> \n",file = zz)
     cat("<tr> \n",file = zz)
     cat("<td>Times &#160</td>  <td>Lower Bounds &#160</td>  <td>Upper Bounds &#160</td> \n",file = zz)
     cat("<td>&alpha;[i]-&alpha;[i-1] &#160</td>  <td>cumulative &alpha; &#160</td> \n",file = zz)
     cat("</tr> \n",file = zz)

     for(i in 1:K)
     {
       cat("<tr> \n",file = zz)
       cat("<td>",round(t[i],digits=3),"</td>",   "<td>",round(lowerBounds[i],digits=4),"</td>",   "<td>",round(upperBounds[i],digits=4),"</td>",
           "<td>",round(probDifference[i],digits=10),"</td>",   "<td>",round(probExit[i],digits=10),"</td> \n",file = zz )
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
   tkdestroy(outTask1Toplevel)
  }

  #frame for the buttons
  buttonFrame<-tkframe(OutputTask1,relief="groove",borderwidth=0)

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
  tkgrid(OutputTask1,sticky="w")
  tkfocus(outTask1Toplevel)

}#end <--*function(...)*
