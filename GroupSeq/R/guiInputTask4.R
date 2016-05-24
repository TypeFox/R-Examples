"guiInputTask4" <-
function(taskWindow)
{

  ### Initialize Variables ###
  #default inputs
  nMax<-25 # Number of interim analyses is limited to 25
  n<-1 # number of interim analyses
  confidenceLevel<-0.95 # desired power respective confidence level
  alphaInput<-0.05 # 'alphaInput' is input of the desired overall size.
  t<-1  # vector containing interim analyses
  t2<-t # vector containing second time scale - by default t2==t
  t3<-t2 # t3 is t2 divided by the maximum information, if needed
  t2max<-0 # maximum value of t2
  upperBounds<-0 # vector containing upper bounds
  lowerBounds<-0 # vector containing lower bounds
  BoundsSymmetry<-1 # BoundsSymmetry==1 means one-sided bounds, BoundsSymmetry==2 means two-sided symmetric and BoundsSymmetry==3 means asymmetric bounds.
  functionInput<-1 # indicates type I error spending rate function e.g. the function(s) the user choosed
  Zvalue<-0 # standardized statistic (Z value) at the last analysis
  TruncateBoundsInput<-8 ## here 8 equals infinity

  #some status variables (names are self-explanatory)
  nInputEdited<-FALSE
  equallySpacedTimesCheckBoxClicked<-FALSE
  secondTimeScaleCheckBoxClicked<-FALSE
  equallySpacedTimesInput<-TRUE
  secondTimeScaleIsUsedInput<-FALSE
  enterBoundsManually<-FALSE
  truncateBoundsYesNo<-FALSE # default is no truncating of bounds

  #some lists (names are self-explanatory)
  listOfTimePointLabel.unequalTimes<-list() #
  listOfInputFields.unequalTimes<-list()
  listOfEntries.unequalTimes<-list()

  listOfTimePointLabel.secondTimes<-list()
  listOfInputFields.secondTimes<-list()
  listOfEntries.secondTimes<-list()

  listOfTimePointLabel.boundsUPPER<-list() #
  listOfInputFields.boundsUPPER<-list()
  listOfEntries.boundsUPPER<-list()

  listOfTimePointLabel.boundsLOWER<-list() #
  listOfInputFields.boundsLOWER<-list()
  listOfEntries.boundsLOWER<-list()

  #operating variables
  nBackup<-n # backup n
  boundBackup<-BoundsSymmetry # backup BoundsSymmetry
  alpha1<- alphaInput # set alpha
  alpha2<- 0 # alpha2 is need in case of asymmetric bounds
  function1<-functionInput # set function1
  function2<-functionInput # function2 is need in case of asymmetric bounds
  phi1<-1 # optional Parameter referring to Power Family
  phi2<-1 # phi2 is need in case of asymmetric bounds both using Power family

  #Define some Fonts
  fontItalic <- tkfont.create(family="times",size=10,slant="italic")
  fontBig <- tkfont.create(family="times",size=10,weight="bold")


####################################################################################################
#################------------------------------------------------------------#######################
##################  FUNCTIONS THAT HANDLE EVENTS TAKING PLACE IN THE WINDOW ########################
#################------------------------------------------------------------#######################
####################################################################################################

  #########################################################
  # function handles change on number of interim analyses #
  #########################################################
  onChangeInterimAnalyses<-function(nValue)
  {
    #First check whether n has changed
    if(n==nValue)
    {
      #nothing changed - so do nothing
    }
    else #interim times changed
    {
      #set new n
      n<<-nValue

      #we'll have to recompute the lists t,t2 new e.g
      #they got default values equally spaced and also set the bounds to some default values
      t<<-1
      for(i in 1:n)
      {
        t[i]<<-i/n
        upperBounds[i]<<-i
      }
      t2<<-t
      lowerBounds<<- -upperBounds


      #update n in menu bar
      tkdelete(topMenu,0,1)
      tkadd(topMenu,"cascade",label=paste("#Interim Times: K=",as.character(n)),menu=nMenu)

      ### equally or unequally spaced times? get it from the checkbox ###
      equallySpacedTimesInput <- as.logical(as.numeric(tclvalue(equallySpacedTimesCheckBoxValue)))

      # check case unequally checkboxes - grid input fields with number of interim analyses into the frames
      if(!equallySpacedTimesInput)
      {
        #first remove "old" labels and input fields - old n is stored in nBackup
        for(i in 1:nBackup)
        {
  #remove labels and input fields
          tkgrid.remove(listOfTimePointLabel.unequalTimes[[i]],listOfEntries.unequalTimes[[i]])
        }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfInputFields.unequalTimes<<-list()
        listOfTimePointLabel.unequalTimes<<-list()
        listOfEntries.unequalTimes<<-list()

        #create new labels by calling function 'onCheckBoxEquallySpacedTimes()'
        onCheckBoxEquallySpacedTimes()

      }#end <--*if(!equallySpacedTimesInput)*


      ### Second Time scale - will it be used? ###
      secondTimeScaleIsUsedInput <- as.logical(as.numeric(tclvalue(secondTimeScaleIsUsedCheckBoxValue)))

      #check case second times scale checkbox is activated
      if(secondTimeScaleIsUsedInput)
      {
        #first remove "old" labels and input fields - old n is stored in nBackup
        for(i in 1:nBackup)
        {
  #remove labels and input fields
          tkgrid.remove(listOfTimePointLabel.secondTimes[[i]],listOfEntries.secondTimes[[i]])
        }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfTimePointLabel.secondTimes<<-list() #
        listOfInputFields.secondTimes<<-list()
        listOfEntries.secondTimes<<-list()

        #create new labels and input fields by calling function 'onCheckBoxSecondTimeScale()'
        onCheckBoxSecondTimeScale()

      }#end <--*if(secondTimeScaleIsUsedInput)*

      ### if user enters bounds manually we have to adapt appropriate input fields
      if(enterBoundsManually)
      {

        #check on symmetric or asymmetric bounds used
        if(!BoundsSymmetry==3)
        {
          #symmetric
          #first remove "old" labels and input fields - old n is stored in nBackup
          for(i in 1:nBackup)
          {
    #remove labels and input fields
            tkgrid.remove(listOfTimePointLabel.boundsUPPER[[i]],listOfEntries.boundsUPPER[[i]])
          }
          #set the lists to NULL otherwise we would duplicate entries in a next loop
          listOfInputFields.boundsUPPER<<-list()
          listOfTimePointLabel.boundsUPPER<<-list()
          listOfEntries.boundsUPPER<<-list()

          #create new labels and input fields by calling function 'onManualBounds()'
          onManualBounds()

        }#end <--*if(!BoundsSymmetry==3)*

        else #asymmetric bounds - do the same as with symmetric for both frames
        {
          for(i in 1:nBackup)
          {
    #remove labels and input fields
            tkgrid.remove(listOfTimePointLabel.boundsUPPER[[i]],listOfEntries.boundsUPPER[[i]])
            tkgrid.remove(listOfTimePointLabel.boundsLOWER[[i]],listOfEntries.boundsLOWER[[i]])
          }
          #set the lists to NULL otherwise we would duplicate entries in a next loop
          listOfInputFields.boundsUPPER<<-list()
          listOfTimePointLabel.boundsUPPER<<-list()
          listOfEntries.boundsUPPER<<-list()
          listOfInputFields.boundsLOWER<<-list()
          listOfTimePointLabel.boundsLOWER<<-list()
          listOfEntries.boundsLOWER<<-list()

          #create new labels and input fields by calling function 'onManualBounds()'
          onManualBounds()
        }#end <--*else #asymmetric bounds - do the same as with symmetric for both frames*

      }#end <--*if(enterBoundsManually)*


    }#end <--*else #interim times changed*

    #update nBackup
    nBackup<<-n
  }#end <--*onChangeInterimAnalyses<-function(nValue)*


  ###################################################################################
  # function handles a click on checkbox for equally/unequally spaced interim times #
  ###################################################################################
  onCheckBoxEquallySpacedTimes <- function()
  {
    #equally or unequally spaced times? get it from the checkbox
    equallySpacedTimesInput <- as.logical(as.numeric(tclvalue(equallySpacedTimesCheckBoxValue)))

    # case unequally checkboxes - grid input fields with number of interim analyses into the frames
    if(!equallySpacedTimesInput)
    {
      for(i in 1:n)
      {
        #create label in a list thus we can dynamically change number of input fields
        listOfTimePointLabel.unequalTimes<<-c(listOfTimePointLabel.unequalTimes,list(tklabel(unEquallyDynamicFrame, text=paste("time",as.character(i)))))

        #We need a list of Input Fields to be able to save the dynamic created tclVar's
        listOfInputFields.unequalTimes<<-c(listOfInputFields.unequalTimes,list(tclVar(as.character(t[i]))))
        listOfEntries.unequalTimes<<-c(listOfEntries.unequalTimes, list(tkentry(unEquallyDynamicFrame,width="11",textvariable=as.character(listOfInputFields.unequalTimes[[i]]))))

        #put label with Input field
        tkgrid(listOfTimePointLabel.unequalTimes[[i]],listOfEntries.unequalTimes[[i]])
        tkgrid.configure(listOfTimePointLabel.unequalTimes[[i]],sticky="nw")
        tkgrid.configure(listOfEntries.unequalTimes[[i]],sticky="nw")
      }#end <--*for*
    #put frame
    tkgrid(unEquallyDynamicFrame)
    }#end <--*if*

    else #equally spaced - remove all input fields cause they should disappear in the window
    {
      #fade out frame
      tkgrid.forget(unEquallyDynamicFrame)

      for(i in 1:n)
      {
        #remove labels and input fields
tkgrid.remove(listOfTimePointLabel.unequalTimes[[i]],listOfEntries.unequalTimes[[i]])
      }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfInputFields.unequalTimes<<-list()
        listOfTimePointLabel.unequalTimes<<-list()
        listOfEntries.unequalTimes<<-list()
    }
  }#end <--*function()*

  ##############################################################
  # function handles a click on checkbox for second time scale
  #
  # ATTENTION: this feature of a second time scale is currently NOT used!
  #
  ##############################################################
  onCheckBoxSecondTimeScale <- function()
  {
    #second time scale used?
    secondTimeScaleIsUsedInput <- as.logical(as.numeric(tclvalue(secondTimeScaleIsUsedCheckBoxValue)))

    # case unequally checkboxes - grid input fields with number of interim analyses into the frames
    if(secondTimeScaleIsUsedInput)
    {
      for(i in 1:n)
      {
        #create label in a list thus we can dynamically change number of input fields
        listOfTimePointLabel.secondTimes<<-c(listOfTimePointLabel.secondTimes,list(tklabel(secondTimesDynamicFrame, text=paste("time",as.character(i)))))

        #We need a list of Input Fields to be able to save the dynamic created tclVar's
        listOfInputFields.secondTimes<<-c(listOfInputFields.secondTimes,list(tclVar(as.character(t2[i]))))
        listOfEntries.secondTimes<<-c(listOfEntries.secondTimes, list(tkentry(secondTimesDynamicFrame,width="11",textvariable=as.character(listOfInputFields.secondTimes[[i]]))))

        #put label with Input field
        tkgrid(listOfTimePointLabel.secondTimes[[i]],listOfEntries.secondTimes[[i]])
        tkgrid.configure(listOfTimePointLabel.secondTimes[[i]],sticky="nw")
        tkgrid.configure(listOfEntries.secondTimes[[i]],sticky="nw")
      }#end <--*for*
    #put frame
    tkgrid(secondTimesDynamicFrame)
    }#end <--*if*

    else #equally spaced - remove all input fields cause they should disappear in the window
    {
      #fade out frame
      tkgrid.forget(secondTimesDynamicFrame)

      for(i in 1:n)
      {
#remove labels and input fields
tkgrid.remove(listOfTimePointLabel.secondTimes[[i]],listOfEntries.secondTimes[[i]])
      }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfInputFields.secondTimes<<-list()
        listOfTimePointLabel.secondTimes<<-list()
        listOfEntries.secondTimes<<-list()
    }

  }#end <--*onCheckBoxSecondTimeScale <- function()*




  #################################################################
  # function handles a click on checkbox to enter bounds manually #
  #################################################################
  onManualBounds <-function()
  {
    #get checkbox value
    enterBoundsManually<<-as.numeric(tclvalue(manualBoundsCheckBoxValue))


    ### user wants to enter bounds manually ###
    if(enterBoundsManually)
    {
      #fade out frame with choice of functions and truncating bounds checkbox
      tkgrid.forget(computedBoundsFrame)
      tkgrid.forget(TruncateBoundsFrame)

      #symmetric or asymmetric bounds? get checkbox value and check it out
      BoundsSymmetry <<- as.numeric(tclvalue(SymmetryValue))


      #at least we need one input field
      for(i in 1:n)
      {
        #create label in a list thus we can dynamically change number of input fields
        listOfTimePointLabel.boundsUPPER<<-c(listOfTimePointLabel.boundsUPPER,list(tklabel(manualBoundsUPPERframe.InputFields, text=paste("time",as.character(i)))))

        #We need a list of Input Fields to be able to save the dynamic created tclVar's
        listOfInputFields.boundsUPPER<<-c(listOfInputFields.boundsUPPER,list(tclVar(as.character(upperBounds[i]))))
        listOfEntries.boundsUPPER<<-c(listOfEntries.boundsUPPER, list(tkentry(manualBoundsUPPERframe.InputFields,width="11",textvariable=as.character(listOfInputFields.boundsUPPER[[i]]))))

        #put label with Input field
        tkgrid(listOfTimePointLabel.boundsUPPER[[i]],listOfEntries.boundsUPPER[[i]])
        tkgrid.configure(listOfTimePointLabel.boundsUPPER[[i]],sticky="nw")
        tkgrid.configure(listOfEntries.boundsUPPER[[i]],sticky="nw")
      }#end <--*for*

      #if asymmetric bounds we need an additional second input field
      if(BoundsSymmetry==3)
      {
        for(i in 1:n)
        {
          #create label in a list thus we can dynamically change number of input fields
          listOfTimePointLabel.boundsLOWER<<-c(listOfTimePointLabel.boundsLOWER,list(tklabel(manualBoundsLOWERframe.InputFields, text=paste("time",as.character(i)))))

          #We need a list of Input Fields to be able to save the dynamic created tclVar's
          listOfInputFields.boundsLOWER<<-c(listOfInputFields.boundsLOWER,list(tclVar(as.character(lowerBounds[i]))))
          listOfEntries.boundsLOWER<<-c(listOfEntries.boundsLOWER, list(tkentry(manualBoundsLOWERframe.InputFields,width="11",textvariable=as.character(listOfInputFields.boundsLOWER[[i]]))))

          #put label with Input field
          tkgrid(listOfTimePointLabel.boundsLOWER[[i]],listOfEntries.boundsLOWER[[i]])
          tkgrid.configure(listOfTimePointLabel.boundsLOWER[[i]],sticky="nw")
          tkgrid.configure(listOfEntries.boundsLOWER[[i]],sticky="nw")
        }#end <--*for*
        tkgrid(manualBoundsUPPERframe,manualBoundsLOWERframe,sticky="nw")
      }
      else
      {
        tkgrid(manualBoundsUPPERframe,sticky="nw")
      }

      #put whole frame
      tkgrid(manualBoundsFrame,sticky="w")

    }#end <--*if(enterBoundsManually)*

    else #user deactivated checkbox to enter bounds manual
    {
      #fade out frame containing input fields for manual bounds and fade in truncating bounds checkbox
      tkgrid.forget(manualBoundsUPPERframe)
      tkgrid.forget(manualBoundsLOWERframe)
      tkgrid.forget(manualBoundsFrame)
      tkgrid(TruncateBoundsFrame)

      #check on symmetric or asymmetric bounds used
      if(!BoundsSymmetry==3)
      {
        #symmetric
        #remove labels and input fields and clear lists
        for(i in 1:n)
        {
          #remove labels and input fields
          tkgrid.remove(listOfTimePointLabel.boundsUPPER[[i]],listOfEntries.boundsUPPER[[i]])
        }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfInputFields.boundsUPPER<<-list()
        listOfTimePointLabel.boundsUPPER<<-list()
        listOfEntries.boundsUPPER<<-list()

      }#end <--*if(!BoundsSymmetry==3)*

      else #asymmetric bounds - do the same as with symmetric for both frames
      {
        for(i in 1:n)
        {
          #remove labels and input fields and clear lists
          tkgrid.remove(listOfTimePointLabel.boundsUPPER[[i]],listOfEntries.boundsUPPER[[i]])
          tkgrid.remove(listOfTimePointLabel.boundsLOWER[[i]],listOfEntries.boundsLOWER[[i]])
        }
        #set the lists to NULL otherwise we would duplicate entries in a next loop
        listOfInputFields.boundsUPPER<<-list()
        listOfTimePointLabel.boundsUPPER<<-list()
        listOfEntries.boundsUPPER<<-list()
        listOfInputFields.boundsLOWER<<-list()
        listOfTimePointLabel.boundsLOWER<<-list()
        listOfEntries.boundsLOWER<<-list()

      }#end <--*else #asymmetric bounds - do the same as with symmetric for both frames*

      #call onBoundsChosen which will fade in the needed frames
      onBoundsChosen()
      tkgrid(computedBoundsFrame)

    }#end <--*else #user deactivated checkbox to enter bounds manual*

  }#end <--*onManualBounds <-function()*


  ####################################################################
  # functions handle a click on CONFIRM-Buttons to choose a function #
  ####################################################################
  onConfirmFunction1 <-function()
  {
    #check whether we got symmetric or asymmetric bounds
    BoundsSymmetry <<- as.numeric(tclvalue(SymmetryValue))

    #get value from listbox and ask for additional parameters if necessary
    if( BoundsSymmetry==1 || BoundsSymmetry==2)
    {
      function1<<-as.numeric(tkcurselection(listBoxFunction1of1))+1
    }
    else
    {
      function1<<-as.numeric(tkcurselection(listBoxFunction1of2))+1
    }

    #check whether user selected a function
    if( length(function1)==0 )
    {
      tkmessageBox(message="You must select a function!",icon="error",type="ok")
    }

    else # handle select
    {
      ### ASYMMETRIC ###
      if(BoundsSymmetry==3)
      {
        #first of all remove earlier input field
        tkgrid.remove(phiLabel1of2,entry.functionParameter1of2)

        #case Power Family
        if (function1==3)
        {
          #set new label and input field
          phiLabel1of2<<-tklabel(additionalParametersFrame1of2,text="Enter Paramter phi>0:")
          entry.functionParameter1of2 <<-tkentry(additionalParametersFrame1of2,width="6",textvariable=phi1of2)
          tkgrid(phiLabel1of2,sticky="w")
          tkgrid(entry.functionParameter1of2,sticky="w")
        }

        #case Hwang-Shih-DeCani family
        else if (function1==4)
             {
               #set new label and input field
               phiLabel1of2<<-tklabel(additionalParametersFrame1of2,text="Enter Parameter phi=/=0:")
               entry.functionParameter1of2 <<-tkentry(additionalParametersFrame1of2,width="6",textvariable=phi1of2)
               tkgrid(phiLabel1of2,sticky="w")
               tkgrid(entry.functionParameter1of2,sticky="w")
             }
             #case no additional parameters needed
             else
             {
               #do nothing else
             }
      }#end <--*if(BoundsSymmetry==3)*

      else ### SYMMETRIC - do the same in other frame###
      {
        #first of all remove earlier input field
        tkgrid.remove(phiLabel1of1,entry.functionParameter1of1)

        #get value from listbox and ask for additional parameters if necessary
        functionInput <<- as.numeric(tkcurselection(listBoxFunction1of1))+1

        #case Power Family
        if (function1==3)
        {
          #set new label and input field
          phiLabel1of1<<-tklabel(additionalParametersFrame1of1,text="Enter Paramter phi>0:")
          entry.functionParameter1of1 <<-tkentry(additionalParametersFrame1of1,width="6",textvariable=phi1of1)
          tkgrid(phiLabel1of1,sticky="w")
          tkgrid(entry.functionParameter1of1,sticky="w")
        }

        #case Hwang-Shih-DeCani family
        else if (function1==4)
             {
               #set new label and input field
               phiLabel1of1<<-tklabel(additionalParametersFrame1of1,text="Enter Parameter phi=/=0:")
               entry.functionParameter1of1 <<-tkentry(additionalParametersFrame1of1,width="6",textvariable=phi1of1)
               tkgrid(phiLabel1of1,sticky="w")
               tkgrid(entry.functionParameter1of1,sticky="w")
             }
             #case no additional parameters needed
             else
             {
               #do nothing else
             }
      }#end <--*else ### SYMMETRIC - do the same in other frame###     *
    }#end <--*else # handle select*
  }#end <--*onConfirmFunction1 <-function()*




  onConfirmFunction2 <-function()
  {
    #get value from listbox and ask for additional parameters if necessary
    function2<<-as.numeric(tkcurselection(listBoxFunction2of2))+1

    #check whether user selected a function
    if( length(function2)==0 )
    {
      tkmessageBox(message="You must have select a function!",icon="error",type="ok")
    }

    else # handle select
    {
      #first of all remove earlier input field
      tkgrid.remove(phiLabel2of2,entry.functionParameter2of2)

      #case Power Family
      if (function2==3)
      {
        #set new label and input field
        phiLabel2of2<<-tklabel(additionalParametersFrame2of2,text="Enter Paramter phi>0:")
        entry.functionParameter2of2 <<-tkentry(additionalParametersFrame2of2,width="6",textvariable=phi2of2)
        tkgrid(phiLabel2of2,sticky="w")
        tkgrid(entry.functionParameter2of2,sticky="w")
      }

      #case Hwang-Shih-DeCani family
      else if (function2==4)
           {
             #set new label and input field
             phiLabel2of2<<-tklabel(additionalParametersFrame2of2,text="Enter Parameter phi=/=0:")
             entry.functionParameter2of2 <<-tkentry(additionalParametersFrame2of2,width="6",textvariable=phi2of2)
             tkgrid(phiLabel2of2,sticky="w")
             tkgrid(entry.functionParameter2of2,sticky="w")
           }
           #case no additional parameters needed
           else
           {
             #do nothing else
           }
    }#end <--*else # handle select*
  }#end <--*onConfirmFunction2 <-function()*



  ##################################################################
  # function handles a click on Radio Button for ASYMMETRIC BOUNDS #
  ##################################################################
  onBoundsChosen <- function()
  {
    #check whether we got asymmetric bounds
    BoundsSymmetry <<- as.numeric(tclvalue(SymmetryValue))

    #check whether bounds are computed or entered manually by user
    #get checkbox value
    enterBoundsManually<-as.numeric(tclvalue(manualBoundsCheckBoxValue))



    #Asymmetric Bounds!
    if( (BoundsSymmetry==1 || BoundsSymmetry==2) & boundBackup!=3)
    {
      #nothing to change
    }
    else if( (BoundsSymmetry==1 || BoundsSymmetry==2) & boundBackup==3)
         {
           #if users enters bounds manually update frames, if necessary
           if(enterBoundsManually)
           {
             tkgrid.forget(manualBoundsUPPERframe)
             tkgrid.forget(manualBoundsLOWERframe)
             tkgrid.forget(manualBoundsFrame)
             onManualBounds()
           }
           #exchange frames
           tkgrid.remove(nonSymmetricBoundsFrame)
           tkgrid(symmetricBoundsFrame,sticky="nw")
         }
         else if(BoundsSymmetry==3 & boundBackup!=3)
              {
                #if users enters bounds manually update frames, if necessary
                if(enterBoundsManually)
                {
                  tkgrid.forget(manualBoundsUPPERframe)
                  tkgrid.forget(manualBoundsLOWERframe)
                  tkgrid.forget(manualBoundsFrame)
                  onManualBounds()
                }
                #exchange frames
                tkgrid.remove(symmetricBoundsFrame)
                tkgrid(nonSymmetricBoundsFrame,sticky="nw")

              }
  #update boundBackup
  boundBackup<<-BoundsSymmetry
  }


  #######################################################
  # function handles a click on Trunate Bounds Checkbox #
  #######################################################
  onTruncateCheckbox <- function()
  {
    #checkbox activated?
    truncateBoundsYesNo <<- as.logical(as.numeric(tclvalue(TruncateBoundsCheckBoxValue)))

    if(truncateBoundsYesNo) #activated
    {
      ##grid edit box
      tkgrid(TruncateDynamicFrame)
    }
    else #deactivated
    {
      #ungrid edit box
      tkgrid.forget(TruncateDynamicFrame)
    }
  }


  ##################################################
  # function handles a click on 'CALCULATE'-Button #
  ##################################################
  OnCalculateInputTask1 <- function()
  {
    readyForCalculate <- TRUE

    #get values from checkboxes, listboxes and radiobuttons
    equallySpacedTimesInput <<- as.logical(as.numeric(tclvalue(equallySpacedTimesCheckBoxValue)))
    secondTimeScaleIsUsedInput <<- as.logical(as.numeric(tclvalue(secondTimeScaleIsUsedCheckBoxValue)))
    BoundsSymmetry <<- as.numeric(tclvalue(SymmetryValue))
    truncateBoundsYesNo <<- as.logical(as.numeric(tclvalue(TruncateBoundsCheckBoxValue)))
    Zvalue<-as.numeric(tclvalue(ZvalueTclVar))

    #get and check power input to be in (0,1]
    confidenceLevel<<-as.numeric(tclvalue(powerTclVar))
    if( !( confidenceLevel>0 & confidenceLevel<1) )
    {
      readyForCalculate<-FALSE
      tkmessageBox(message="Incorrect Power value entered - please correct!",icon="error",type="ok")
    }


    #truncation point set?
    if(truncateBoundsYesNo)
    {
      TruncateBoundsInput <<- abs(as.numeric(tclvalue(boundsTruncation)))
    }
    else
    {
      TruncateBoundsInput<-8
    }

    #evaluate whether function is used to compute bounds or user entered them manually
    if(enterBoundsManually)
    {
      #manually entered bounds
      for(i in 1:n)
      {
        upperBounds[i]<<- as.numeric(tclvalue(listOfInputFields.boundsUPPER[[i]]))
      }

      #elicit lower bounds
      if(BoundsSymmetry==1)
      {
        #one-sided => lower Bounds == -8 (that is -infinity)
        lowerBounds <<- seq(-8,-8,length=n)
      }
      else if(BoundsSymmetry==2)
           {
             #two-sided symmetric
             lowerBounds <<- -upperBounds
           }
           else
           {
             #asymmetric
             for(i in 1:n)
             {
               lowerBounds[i]<<- as.numeric(tclvalue(listOfInputFields.boundsLOWER[[i]]))
             }
           }
    }#end <--*if(enterBoundsManually)*

    else #bounds are computed
    {
      #get chosen function(s) and check alpha
      ###################
      # case asymmetric #
      ###################
      if(BoundsSymmetry==3)
      {
        alpha1 <<- as.numeric(tclvalue(alpha1of2))
        alpha2 <<- as.numeric(tclvalue(alpha2of2))

        #check alpha
        alphaAll<- alpha1 + alpha2
        if( !(alpha1>=0 & alpha1<=1 & alpha2>=0 & alpha2<=1 & alphaAll<=1) )
        {
          readyForCalculate<-FALSE
          tkmessageBox(message="Alpha out of range! Correct it and try again.",icon="error",type="ok")
        }

        #check phi if entered as parameter
        phi1 <<- as.numeric(tclvalue(phi1of2))
        phi2 <<- as.numeric(tclvalue(phi2of2))
        #function for UPPER bounds
        if(function1==3)
        {
          if( !(phi1>0) )
          {
            readyForCalculate<-FALSE
            tkmessageBox(message="Parameter phi in function for UPPER bounds must be >0 !",icon="error",type="ok")
          }
        }
        else if(function1==4)
             {
               if(phi1==0)
               {
                 readyForCalculate<-FALSE
                 tkmessageBox(message="Parameter phi in function for UPPER bounds may NOT be zero!",icon="error",type="ok")
               }
             }

        #same with function for LOWER bounds
        if(function2==3)
        {
          if( !(phi2>0) )
          {
            readyForCalculate<-FALSE
            tkmessageBox(message="Parameter phi in function for LOWER bounds must be >0 !",icon="error",type="ok")
          }
        }
        else if(function2==4)
             {
               if(phi2==0)
               {
                 readyForCalculate<-FALSE
                 tkmessageBox(message="Parameter phi in function for LOWER bounds may NOT be zero!",icon="error",type="ok")
               }
             }
      }#end <--*if(BoundsSymmetry==3)*

      ##################
      # case symmetric #
      ##################
      else #one function cause of symmetric bounds
      {
        alpha1 <<- as.numeric(tclvalue(alpha1of1))
        #check alpha
        if( !(alpha1>=0 & alpha1<=1) )
        {
          readyForCalculate<-FALSE
          tkmessageBox(message="Alpha out of range! Correct it and try again.",icon="error",type="ok")
        }

        #check phi if entered as parameter
        phi1 <<- as.numeric(tclvalue(phi1of1))
        #what function used?
        if(function1==3)
        {
          if( !(phi1>0) )
          {
            readyForCalculate<-FALSE
            tkmessageBox(message="Parameter phi in function must be >0 !",icon="error",type="ok")
          }
        }
        else if(function1==4)
             {
               if(phi1==0)
               {
                 readyForCalculate<-FALSE
                 tkmessageBox(message="Parameter phi in function may NOT be zero!",icon="error",type="ok")
                 }
             }
      }

    }#end <--*else #bounds are computed*
    ##if user typed in unequally spaced times - get them and
    ##check them to be in intervall(0,1] and in right order
    if(!equallySpacedTimesInput)
    {
      tempVal<-0
      interimTimesBad<-FALSE
      for(i in 1:n)
      {
        tempVal[i] <- as.numeric(tclvalue(listOfInputFields.unequalTimes[[i]]))

        if (tempVal[i]<=0) { interimTimesBad<-TRUE }
        if (tempVal[i]>1) { interimTimesBad<-TRUE }
        if (i>1)
        {
          if (tempVal[i]<=tempVal[i-1]) { interimTimesBad<-TRUE }
        }
      }#end <--*for(i in 1:n)*

      ##if times are not good => error and keep old times
        if(interimTimesBad)
        {
          readyForCalculate <- FALSE
          tkmessageBox(message="Bad Interim Times entered - old Times are kept so far! ",icon="error",type="ok")
        }
        ##else take new times
        else
        {
          t<<-tempVal
        }

    }#end <--*if(equallySpacedTimesInput)*
    else
    {
     for(i in 1:n)
     {
       t[i]<<-i/n
     }
    }

    ##if user typed in second time scales - get them and
    ##check them to be in intervall(0,1] and in right order
    if(secondTimeScaleIsUsedInput)
    {
      tempVal<-0
      secondTimesBad<-FALSE
      for(i in 1:n)
      {
        tempVal[i] <- as.numeric(tclvalue(listOfInputFields.secondTimes[[i]]))

        if (tempVal[i]<=0) { secondTimesBad<-TRUE }
        if (tempVal[i]>1) { secondTimesBad<-TRUE }
        if (i>1)
        {
          if (tempVal[i]<=tempVal[i-1]) { secondTimesBad<-TRUE }
        }
      }#end <--*for(i in 1:n)*

      ##if times are not good => error and keep old times
        if(secondTimesBad)
        {
          readyForCalculate <- FALSE
          tkmessageBox(message="Bad Second Time Scale entered - old Times are kept so far! ",icon="error",type="ok")
        }
        ##else take new times
        else
        {
          t2<<-tempVal
        }

      if(readyForCalculate)
      {

        ###########################################################################
        ############### RESCALE SECOND TIME SCALE IF NECESSARY ####################
        ###########################################################################
        # When second time scale is not on (0,1], computation with
        # non-zero drift parameters are incorrect, since the drift
        # always scaled to (0,1].  If the trial is complete (t[n]=1)
        # then t2 can be rescaled as t3 = t2/t2[n].  Otherwise, if
        # the second time scale is to be used for covariances, the
        # user must enter a maximum value.
        #
        # drift[ t[i] ] = drift*t[i]
        #
        # Start with t3 = t2.
        # (t2=t by default e.g. if user did not enter second time scale.)
        t3<<-t2
        t2max<<-0

        ##If t[n]=1, t2[n] is maximum of t2.
        if(t[n]==1)
        {
         tkmessageBox(title="-4- Compute Confidence Interval",message="Second Time scale will be used to determine covariances.",icon="info",type="ok")
         t2max<<-t2[n]
         t3<<-t2/t2max
        }
        else ##Should we try to use 2nd scale?
        {
          response<-tkmessageBox(message="Do you wish to use the 2nd time scale to determine covariances?",
                                 icon="question",type="yesno",default="yes")


          ##--If yes, prompt for maximum of 2nd time scale.--##
          if( as.character(tclvalue(response))=="yes" )
          {
            ##################################################
            # function handles prompting for t2max if needed #
            ##################################################
            t2maxPrompt <- function(title,question,entryInit,entryWidth=4,returnValOnCancel="ID_CANCEL")
            {
          dlg <- tktoplevel(taskWindow)
              tkwm.deiconify(dlg)
              tkgrab.set(dlg)
              tkfocus(dlg)
              tkwm.title(dlg,title)
              textEntryVarTcl <- tclVar(paste(entryInit))
              textEntryWidget <- tkentry(dlg,width=paste(entryWidth),textvariable=textEntryVarTcl)
              tkgrid(tklabel(dlg,text="       "))
              tkgrid(tklabel(dlg,text=question),textEntryWidget)
              tkgrid(tklabel(dlg,text="       "))
              ReturnVal <- returnValOnCancel
              onOK <- function()
              {
                ReturnVal <<- as.numeric(tclvalue(textEntryVarTcl))
                #check whether numeric was entered
                if(is.na(ReturnVal))
                {
                  tkmessageBox(title="ERROR",message="You did not enter a valid numeric value!",icon="error",type="ok")
                }
                #if input numeric check whether the numeric is a valid entry
                else
                {
                  if(ReturnVal<=0)
                  {
                    tkmessageBox(title="ERROR",message="Maximum must be positive! Please try again!",icon="error",type="ok")
                  }
                  else if(ReturnVal<t2[n])
                       {
                         tkmessageBox(title="ERROR",message="The Maximum cannot be smaller than your last seond time scale value!",icon="error",type="ok")
                       }

                       #input ok - go back to main window
                       else
                       {
                         tkgrab.release(dlg)
                         tkdestroy(dlg)
                         tkfocus(task4)
                       }
                }

              }
              onCancel <- function()
              {
                readyForCalculate<<-FALSE
                ReturnVal <<- returnValOnCancel
                tkgrab.release(dlg)
                tkdestroy(dlg)
                tkfocus(task4)
              }
              OK.but     <-tkbutton(dlg,text="   OK   ",command=onOK)
              Cancel.but <-tkbutton(dlg,text=" Cancel ",command=onCancel)
              tkgrid(OK.but,Cancel.but)
              tkgrid(tklabel(dlg,text="    "))

              tkfocus(dlg)
              tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(task4)})
              tkbind(textEntryWidget, "<Return>", onOK)
              tkwait.window(dlg)

              return(ReturnVal)

            }#end <--*t2maxPrompt <- function(...)*


            ReturnVal<-t2maxPrompt("Second Time Scale will be used","Enter the maximum value of the second time scale","" )
            if(ReturnVal!="ID_CANCEL")
            {
              t2max<<-ReturnVal
              #Rescale t2
              t3<<-t2/t2max
            }

          }#end <--*if( as.character(tclvalue(response))=="yes" )*

          else if( as.character(tclvalue(response))=="no" )
          {
            ##Even if the 2nd time scale was on (0,1], if the
            ##maximum information value was not entered, set
            ##t3 to t, which causes t2 to be ignored!#
            t3<<-t

          }

        }#end <--*else ##Should we try to use 2nd scale?*

      }#end <--*if(readyForCalculate)*

    }#end <--*if(secondTimeScaleIsUsedInput)*
    else
    {
      ##set t3 to t, which causes t2 to be ignored!
      t3<<-t
    }

    if(readyForCalculate)
    {
      # second time scale is not used so far --> set t3=t2=t
      t2<-t; t3<-t
      calculateTask4(n,nMax,t,t2,t2max,t3,confidenceLevel,equallySpacedTimesInput,secondTimeScaleIsUsedInput,
                     BoundsSymmetry, c(alpha1,alpha2), c(phi1,phi2), c(function1,function2),TruncateBoundsInput,
                     enterBoundsManually, upperBounds, lowerBounds, Zvalue, taskWindow)
    }
  }#end <--*OnCalculateInputTask1 <- function()*



#######################################################################################################
#################------------------------------------------------------------------####################
##################  FROM HERE ON LABELS AND FRAMES ARE CREATED AND PUT TOGETHER   #####################
#################------------------------------------------------------------------####################
#######################################################################################################

  #Set Toplevel
  task4 <- tktoplevel(taskWindow)
  tkwm.title(task4,"-4- Compute Confidence Interval")

  #Define main Frame
  InputTask4 <- tkframe(task4, relief="groove",borderwidth=2)


  ##--------------------------------------------------------------------------------##
  ##------------------------  number of Interim Times  -----------------------------##
  ##--------------------------------------------------------------------------------##

  #create pull down menu to select interim analyses from n==1 to n==25(=nMax)
  topMenu <- tkmenu(task4)
  tkconfigure(task4,menu=topMenu)
  nMenu <- tkmenu(topMenu,tearoff=FALSE,background="grey",activebackground="red")
  tkadd(nMenu,"command",label="1",command=function() onChangeInterimAnalyses(1))
  tkadd(nMenu,"command",label="2",command=function() onChangeInterimAnalyses(2))
  tkadd(nMenu,"command",label="3",command=function() onChangeInterimAnalyses(3))
  tkadd(nMenu,"command",label="4",command=function() onChangeInterimAnalyses(4))
  tkadd(nMenu,"command",label="5",command=function() onChangeInterimAnalyses(5))
  tkadd(nMenu,"command",label="6",command=function() onChangeInterimAnalyses(6))
  tkadd(nMenu,"command",label="7",command=function() onChangeInterimAnalyses(7))
  tkadd(nMenu,"command",label="8",command=function() onChangeInterimAnalyses(8))
  tkadd(nMenu,"command",label="9",command=function() onChangeInterimAnalyses(9))
  tkadd(nMenu,"command",label="10",command=function() onChangeInterimAnalyses(10))
  tkadd(nMenu,"command",label="11",command=function() onChangeInterimAnalyses(11))
  tkadd(nMenu,"command",label="12",command=function() onChangeInterimAnalyses(12))
  tkadd(nMenu,"command",label="13",command=function() onChangeInterimAnalyses(13))
  tkadd(nMenu,"command",label="14",command=function() onChangeInterimAnalyses(14))
  tkadd(nMenu,"command",label="15",command=function() onChangeInterimAnalyses(15))
  tkadd(nMenu,"command",label="16",command=function() onChangeInterimAnalyses(16))
  tkadd(nMenu,"command",label="17",command=function() onChangeInterimAnalyses(17))
  tkadd(nMenu,"command",label="18",command=function() onChangeInterimAnalyses(18))
  tkadd(nMenu,"command",label="19",command=function() onChangeInterimAnalyses(19))
  tkadd(nMenu,"command",label="20",command=function() onChangeInterimAnalyses(20))
  tkadd(nMenu,"command",label="21",command=function() onChangeInterimAnalyses(21))
  tkadd(nMenu,"command",label="22",command=function() onChangeInterimAnalyses(22))
  tkadd(nMenu,"command",label="23",command=function() onChangeInterimAnalyses(23))
  tkadd(nMenu,"command",label="24",command=function() onChangeInterimAnalyses(24))
  tkadd(nMenu,"command",label="25",command=function() onChangeInterimAnalyses(25))
  tkadd(topMenu,"cascade",label=paste("#Interim Times: K= ",as.character(n)),menu=nMenu)

  tkgrid(tklabel(InputTask4,text="")) # Blank line


  ##--------------------------------------------------------------------------------##
  ##-------------  Interim Times equally or unequally spaced? ----------------------##
  ##-------------       Second Time Scale will be used?       ----------------------##
  ##--------------------------------------------------------------------------------##

  ## prepare Frames
  interimTimesFrame<- tkframe(InputTask4,relief="groove",borderwidth=0)
  equallySpacedTimesFrame <- tkframe(interimTimesFrame,relief="groove",borderwidth=0)
  secondTimesFrame <- tkframe(interimTimesFrame,relief="groove",borderwidth=0)

  #again we need a frame for dynamic working in it to not affect
  #the format of 'equallySpacedTimesFrame' respective 'secondTimesFrame'
  equallySpacedLabelFrame<-tkframe(equallySpacedTimesFrame,relief="groove",borderwidth=0)
  unEquallyDynamicFrame<-tkframe(equallySpacedTimesFrame,relief="groove",borderwidth=0)
  secondTimesLabelFrame<-tkframe(secondTimesFrame,relief="groove",borderwidth=0)
  secondTimesDynamicFrame<-tkframe(secondTimesFrame,relief="groove",borderwidth=0)

  #Default is Equally Spaced Times and no Second Time Scale
  #create Checkboxes
  #equally spaced
  equallySpacedTimesCheckBox<-tkcheckbutton(equallySpacedLabelFrame,command=onCheckBoxEquallySpacedTimes)
  equallySpacedTimesCheckBoxValue <- tclVar(as.character(as.numeric(equallySpacedTimesInput)))
  tkconfigure(equallySpacedTimesCheckBox,variable=equallySpacedTimesCheckBoxValue)
  #second time scale
  secondTimeScaleIsUsedCheckBox<-tkcheckbutton(secondTimesLabelFrame,command=onCheckBoxSecondTimeScale)
  secondTimeScaleIsUsedCheckBoxValue <- tclVar(as.character(as.numeric(secondTimeScaleIsUsedInput)))
  tkconfigure(secondTimeScaleIsUsedCheckBox,variable=secondTimeScaleIsUsedCheckBoxValue)

  #put checkbox and other frames together
  equallyTimesBoxLabel<-tklabel(equallySpacedLabelFrame,text="Equally Spaced Times")
  secondTimesBoxLabel<-tklabel(secondTimesLabelFrame,text="Use Second Time Scale")
  tkgrid(equallyTimesBoxLabel,equallySpacedTimesCheckBox)
  # tkgrid(secondTimesBoxLabel,secondTimeScaleIsUsedCheckBox) - this feature is currently removed
  tkgrid(equallySpacedTimesFrame,secondTimesFrame,sticky="n")
  tkgrid(equallySpacedLabelFrame,sticky="n")
  tkgrid(secondTimesLabelFrame,sticky="n")
  tkgrid(unEquallyDynamicFrame,sticky="nw")
  tkgrid(secondTimesDynamicFrame,sticky="n")
  tkgrid(interimTimesFrame,sticky="w")
  tkgrid(tklabel(InputTask4,text="")) # Blank line


  ##-- Level for Confidence Intervall --##
  ## confidence Intervall replaces last bound with Zvalue
  ZvalueTclVar<-tclVar(as.character(Zvalue))
  ZvalueFrame <- tkframe(InputTask4,relief="groove",borderwidth=0)
  ZvalueLabel<-tklabel(ZvalueFrame,text="Enter the standardized statistic (Z value) at the last analysis:")
  entry.Zvalue <-tkentry(ZvalueFrame,width="6",textvariable=ZvalueTclVar)
  #grid it
  tkgrid(ZvalueLabel,entry.Zvalue, sticky="w")
  tkgrid(ZvalueFrame,sticky="w")


  #Desired Power
  powerTclVar<-tclVar(as.character(confidenceLevel))
  powerFrame <- tkframe(InputTask4,relief="groove",borderwidth=0)
  powerLabel<-tklabel(powerFrame,text="Enter Confidence Level - it must be in (0,1) :")
  entry.power <-tkentry(powerFrame,width="6",textvariable=powerTclVar)
  #grid it
  tkgrid(powerLabel,entry.power, sticky="w")
  tkgrid(tklabel(powerFrame,text="")) # Blank line
  tkgrid(powerFrame,sticky="w")

  ###One- or Two-Sided Bounds or asymmetric Bounds###
  #create frames
  boundsLabelFrame <- tkframe(InputTask4,relief="groove",borderwidth=0)
  boundsRadioButtonFrame <- tkframe(InputTask4,relief="groove",borderwidth=0)

  #create radio buttons
  oneSided <- tkradiobutton(boundsRadioButtonFrame,command=onBoundsChosen)
  twoSided <- tkradiobutton(boundsRadioButtonFrame,command=onBoundsChosen)
  asymmetric <- tkradiobutton(boundsRadioButtonFrame,command=onBoundsChosen)
  SymmetryValue <- tclVar(as.character(BoundsSymmetry))
  tkconfigure(oneSided,variable=SymmetryValue,value="1")
  tkconfigure(twoSided,variable=SymmetryValue,value="2")
  tkconfigure(asymmetric,variable=SymmetryValue,value="3")

  #grid labels and buttons together
  tkgrid(tklabel(boundsLabelFrame,text="One-, Two-sided symmetric or asymmetric bounds?"),sticky="w")
  tkgrid(tklabel(boundsRadioButtonFrame,text="One-Sided "),oneSided)
  tkgrid(tklabel(boundsRadioButtonFrame,text="Two-Sided "),twoSided)
  tkgrid(tklabel(boundsRadioButtonFrame,text="Asymmetric "),asymmetric)

  #put frames
  tkgrid(boundsLabelFrame,sticky="w")
  tkgrid(boundsRadioButtonFrame,sticky="w")
  tkgrid(tklabel(InputTask4,text="")) # Blank line

  ################################################################################
  ### User could enter bounds manually instead of let them be computed (default)##
  ################################################################################
  manualOrComputedBoundsFrame<-tkframe(InputTask4,relief="groove",borderwidth=0)
  LabelManualOrComputedBoundsFrame<-tkframe(manualOrComputedBoundsFrame,relief="groove",borderwidth=0)
  computedBoundsFrame<-tkframe(manualOrComputedBoundsFrame,relief="groove",borderwidth=0)
  manualBoundsFrame<-tkframe(manualOrComputedBoundsFrame,relief="groove",borderwidth=0)
  #upper bounds
  manualBoundsUPPERframe<-tkframe(manualBoundsFrame,relief="groove",borderwidth=0)
  manualBoundsUPPERframe.Label<-tkframe(manualBoundsUPPERframe,relief="groove",borderwidth=0)
  manualBoundsUPPERframe.InputFields<-tkframe(manualBoundsUPPERframe,relief="groove",borderwidth=0)
  #lower bounds
  manualBoundsLOWERframe<-tkframe(manualBoundsFrame,relief="groove",borderwidth=0)
  manualBoundsLOWERframe.Label<-tkframe(manualBoundsLOWERframe,relief="groove",borderwidth=0)
  manualBoundsLOWERframe.InputFields<-tkframe(manualBoundsLOWERframe,relief="groove",borderwidth=0)

  #create checkbox with label
  manualBoundsCheckbox <- tkcheckbutton(LabelManualOrComputedBoundsFrame,command=onManualBounds)
  manualBoundsCheckBoxValue <- tclVar(as.character(as.numeric(enterBoundsManually)))
  tkconfigure(manualBoundsCheckbox,variable=manualBoundsCheckBoxValue)
  manualBoundsCheckboxLabel<-tklabel(LabelManualOrComputedBoundsFrame,text="Enter Bounds Manually")
  tkgrid(manualBoundsCheckboxLabel,manualBoundsCheckbox)
  tkgrid(tklabel(LabelManualOrComputedBoundsFrame,text="")) # Blank line
  tkgrid.configure(manualBoundsCheckboxLabel,sticky="w")
  tkgrid.configure(manualBoundsCheckbox,sticky="w")

  #create labels for manual entering and grid alltogether
  manualBoundsUPPERlabel<-tklabel(manualBoundsUPPERframe.Label,text="Enter UPPER Bounds(standardized)     ")
  tkgrid(manualBoundsUPPERlabel,sticky="w")
  tkgrid(manualBoundsUPPERframe.Label,sticky="w")
  tkgrid(manualBoundsUPPERframe.InputFields,sticky="w")
  manualBoundsLOWERlabel<-tklabel(manualBoundsLOWERframe.Label,text="Enter LOWER Bounds (standardized)")
  tkgrid(manualBoundsLOWERlabel,sticky="w")
  tkgrid(manualBoundsLOWERframe.Label,sticky="w")
  tkgrid(manualBoundsLOWERframe.InputFields,sticky="w")

  ### Significance Level(s) alpha and function(s) to be used to calculate bounds###
  ## if user choses asymmetric bounds two different functions could be used ##
  ##create Frames
  alphaAndFunctionsFrame<-tkframe(computedBoundsFrame,relief="groove",borderwidth=0)
  symmetricBoundsFrame<-tkframe(alphaAndFunctionsFrame,relief="groove",borderwidth=0)
  nonSymmetricBoundsFrame<-tkframe(alphaAndFunctionsFrame,relief="groove",borderwidth=0)

  ##Default alpha1==0.05, alpha2==0
  alpha1of1 <- tclVar(as.character(alphaInput))
  alpha1of2 <- tclVar(as.character("0.025"))
  alpha2of2 <- tclVar(as.character("0.025"))

  #########################################################
  ### case symmetric bounds or one-sided test (default) ###
  #########################################################
  #frames
  alphaFrame1of1 <- tkframe(symmetricBoundsFrame,relief="groove",borderwidth=0)
  functionsFrame1of1 <- tkframe(symmetricBoundsFrame,relief="groove",borderwidth=0)
  additionalParametersFrame1of1<-tkframe(symmetricBoundsFrame,relief="groove",borderwidth=0)

  ##create Labels for alpha
  alphaLabel1of1<-tklabel(alphaFrame1of1,text="Significance Level: alpha=")
  entry.alpha1of1 <-tkentry(alphaFrame1of1,width="6",textvariable=alpha1of1)

  #create Listbox for function choice
  functionLabel1of1<-tklabel(functionsFrame1of1,text="What function should be used?")
  listBoxFunction1of1<-tklistbox(functionsFrame1of1,height=5,width=30,selectmode="single",background="grey")
  functionChoice1of1 <- c("(1) O'Brien-Fleming Type","(2) Pocock Type",
    "(3) Power Family: alpha* t^phi","(4) Hwang-Shih-DeCani Family","(5) Exact Pocock Bounds")
  for (i in (1:5))
  {
    tkinsert(listBoxFunction1of1,"end",functionChoice1of1[i])
  }
  tkselection.set(listBoxFunction1of1, functionInput-1)  # Default function is O'Brien-Fleming Type.  Indexing starts at zero.

  #create and put button to confirm a function because for example in case of 'Power family: alpha* t^phi'
  #user has to enter additional parameter 'phi'
  confirmFun.button1of1 <-tkbutton(functionsFrame1of1,text=" CONFIRM FUNCTION ",command=onConfirmFunction1)

  #create variable for edit box which we will need if additional parameters must be entered
  #edit box is unvisible at beginning since default function O'Brien-Fleming Type does not need any additional parameters
  phi1of1 <- tclVar(as.character(phi1))
  phiLabel1of1<-tklabel(additionalParametersFrame1of1,text="")
  entry.functionParameter1of1 <-tkentry(additionalParametersFrame1of1,width="3",textvariable=phi1of1)

  #grid together
  #alpha
  tkgrid(alphaLabel1of1,entry.alpha1of1)
  tkgrid.configure(alphaLabel1of1,sticky="w")
  tkgrid.configure(entry.alpha1of1,sticky="w")
  tkgrid(tklabel(alphaFrame1of1,text="")) # Blank line
  tkgrid(functionLabel1of1,sticky="w")
  tkgrid(listBoxFunction1of1)

  #put frames and button
  tkgrid(alphaFrame1of1,sticky="w")
  tkgrid(functionsFrame1of1,additionalParametersFrame1of1)
  tkgrid(confirmFun.button1of1)
  tkgrid.configure(functionsFrame1of1,sticky="w")
  tkgrid(tklabel(symmetricBoundsFrame,text="")) # Blank line

  #Finally grid frame for symmetric case as default
  tkgrid(LabelManualOrComputedBoundsFrame,sticky="w")
  tkgrid(symmetricBoundsFrame,sticky="nw")
  tkgrid(alphaAndFunctionsFrame,sticky="w")
  tkgrid(computedBoundsFrame,sticky="w")
  tkgrid(manualOrComputedBoundsFrame,sticky="w")


  ##############################
  ### case Asymmetric bounds ###
  ##############################
  #frames
  alphaFrame1of2 <- tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)
  functionsFrame1of2 <- tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)
  additionalParametersFrame1of2<-tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)
  alphaFrame2of2 <- tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)
  functionsFrame2of2 <- tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)
  additionalParametersFrame2of2<-tkframe(nonSymmetricBoundsFrame,relief="groove",borderwidth=0)

  ##create Labels for alpha
  alphaLabel1of2<-tklabel(alphaFrame1of2,text="UPPER Bounds: alpha=")
  entry.alpha1of2 <-tkentry(alphaFrame1of2,width="6",textvariable=alpha1of2)

  alphaLabel2of2<-tklabel(alphaFrame2of2,text="LOWER Bounds: alpha=")
  entry.alpha2of2 <-tkentry(alphaFrame2of2,width="6",textvariable=alpha2of2)

  #create Listboxes for function choice
  #################
  # List Box 1of2 #
  #################
  functionLabel1of2<-tklabel(functionsFrame1of2,text="Choose Function for UPPER Bounds")
  listBoxFunction1of2<-tklistbox(functionsFrame1of2,height=5,width=30,selectmode="single",background="grey")
  functionChoice1of2 <- c("(1) O'Brien-Fleming Type","(2) Pocock Type",
    "(3) Power Family: alpha* t^phi","(4) Hwang-Shih-DeCani Family","(5) Exact Pocock Bounds")
  for (i in (1:5))
  {
    tkinsert(listBoxFunction1of2,"end",functionChoice1of2[i])
  }

  #create and put first Confirm button which "commands" same function as in symmetric case did
  confirmFun.button1of2 <-tkbutton(functionsFrame1of2,text=" CONFIRM FUNCTION ",command=onConfirmFunction1)
  #edit box for parameter phi[1]
  phi1of2 <- tclVar(as.character(phi1))
  phiLabel1of2<-tklabel(additionalParametersFrame1of2,text="")
  entry.functionParameter1of2 <-tkentry(additionalParametersFrame1of2,width="3",textvariable=phi1of2)

  #################
  # List Box 2of2 #
  #################
  functionLabel2of2<-tklabel(functionsFrame2of2,text="Choose Function for LOWER Bounds")
  listBoxFunction2of2<-tklistbox(functionsFrame2of2,height=5,width=30,selectmode="single",background="grey")
  functionChoice2of2 <- c("(1) O'Brien-Fleming Type","(2) Pocock Type",
    "(3) Power Family: alpha* t^phi","(4) Hwang-Shih-DeCani Family","(5) Exact Pocock Bounds")
  for (i in (1:5))
  {
    tkinsert(listBoxFunction2of2,"end",functionChoice2of2[i])
  }

  #create and put first Confirm button
  confirmFun.button2of2 <-tkbutton(functionsFrame2of2,text=" CONFIRM FUNCTION ",command=onConfirmFunction2)

  #edit box for parameter phi[2]
  phi2of2 <- tclVar(as.character(phi2))
  phiLabel2of2<-tklabel(additionalParametersFrame2of2,text="")
  entry.functionParameter2of2 <-tkentry(additionalParametersFrame2of2,width="3",textvariable=phi2of2)

  #grid together

  #1of2
  tkgrid(alphaLabel1of2,entry.alpha1of2)
  tkgrid.configure(alphaLabel1of2,sticky="w")
  tkgrid.configure(entry.alpha1of2,sticky="w")
  tkgrid(functionLabel1of2,sticky="w")
  tkgrid(listBoxFunction1of2)

  #put frames and button
  tkgrid(alphaFrame1of2,sticky="w")
  tkgrid(functionsFrame1of2,additionalParametersFrame1of2)
  tkgrid(confirmFun.button1of2)
  tkgrid.configure(functionsFrame1of2,sticky="w")
  tkgrid(tklabel(nonSymmetricBoundsFrame,text="")) # Blank line


  #2of2
  tkgrid(alphaLabel2of2,entry.alpha2of2)
  tkgrid.configure(alphaLabel2of2,sticky="w")
  tkgrid.configure(entry.alpha2of2,sticky="w")
  tkgrid(functionLabel2of2,sticky="w")
  tkgrid(listBoxFunction2of2)

  #put frames and button
  tkgrid(alphaFrame2of2,sticky="w")
  tkgrid(functionsFrame2of2,additionalParametersFrame2of2)
  tkgrid(confirmFun.button2of2)
  tkgrid.configure(functionsFrame2of2,sticky="w")
  tkgrid(tklabel(nonSymmetricBoundsFrame,text="")) # Blank line

  #################################################################

  ###Truncate Bounds?###
  PositionTruncateBoundsFrame<-tkframe(InputTask4,relief="groove",borderwidth=0)
  TruncateBoundsFrame <- tkframe(PositionTruncateBoundsFrame,relief="groove",borderwidth=0)
  TruncateLabelFrame <- tkframe(TruncateBoundsFrame,relief="groove",borderwidth=0)
  TruncateDynamicFrame <-tkframe(TruncateBoundsFrame,relief="groove",borderwidth=0)

  #create checkbox
  TruncateBoundsCheckBox<-tkcheckbutton(TruncateLabelFrame,command=onTruncateCheckbox)
  TruncateBoundsCheckBoxValue <- tclVar(as.character(as.numeric(truncateBoundsYesNo)))
  tkconfigure(TruncateBoundsCheckBox,variable=TruncateBoundsCheckBoxValue)

  #create variable for edit box which we will need if user wants truncation of bounds -
  boundsTruncation <- tclVar(as.character(TruncateBoundsInput))
  boundsTruncationLabel<-tklabel(TruncateDynamicFrame,text="Enter Truncation Point:")
  entry.truncationValue <-tkentry(TruncateDynamicFrame,width="3",textvariable=boundsTruncation)

  #put frames
  tkgrid(tklabel(TruncateLabelFrame,text="Truncate standardized Bounds?"),TruncateBoundsCheckBox)
  tkgrid(boundsTruncationLabel,entry.truncationValue,sticky="w")
  tkgrid(TruncateLabelFrame,TruncateDynamicFrame,sticky="w")
  tkgrid(TruncateBoundsFrame,sticky="w")
  tkgrid(PositionTruncateBoundsFrame,sticky="w")
  tkgrid.forget(TruncateDynamicFrame) #default is no Truncating
  tkgrid(tklabel(InputTask4,text="")) # Blank line

  ##put Overall Frame
  tkgrid(InputTask4)


  #frame for the buttons
  buttonFrame<-tkframe(task4,relief="groove",borderwidth=0)

  #create and put button for calculating
  calculate.button <-tkbutton(buttonFrame,text=" CALCULATE ",command=OnCalculateInputTask1)

  # function handles click onto button to Cancel i.e. close current window
  onCancel <- function()
  {
   tkdestroy(task4)
  }
  cancel.button <-tkbutton(buttonFrame,text=" Cancel ",command=onCancel)

  # grid buttons
  tkgrid( tklabel(buttonFrame, text=""))   #blank line
  tkgrid(calculate.button, tklabel(buttonFrame, text="            "),
         cancel.button, sticky="we" )
  tkgrid( tklabel(buttonFrame, text=""))   #blank line
  tkgrid(buttonFrame)

  tkfocus(task4)

}
