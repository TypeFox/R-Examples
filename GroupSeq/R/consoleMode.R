"consoleMode" <-
function()
{
###########################################################################
######################### INITIALIZE VARIABLES ############################
###########################################################################

  #variables for calculation
  t <-0 # the vector of analysis times on (0,1] scale.
  t2<-0 # the second or information time scale, for covariances
  t3<-0 # t2 divided by the maximum information, if needed.
  t2max<-0 #maximum of t2
  n<-0 # number of interim analyses
  alpha<-0 # 'alpha' is the desired overall size.
  lowerBounds <-0 # lowerBounds is the vector of lower standardized boundaries
  upperBounds <-0  # upperBounds is the vector of upper standardized boundaries
  drift<-0 # vector of drift parameters
  Zvalue<-0 # standarized statistic (Z value) at the last analysis
  expectedStoppingTime<-0 #is the expected stopping time (not implemented)
  power <- 0 #power(i) is the prob of rejecting under drift(i)

  probStopping <-0 #is the probaility of reaching ith analysis and stopping.
  probExceedingUpper<-0 # the probaility of reaching ith and exceeding upper.
  probExceedingLower<-0 #is the probaility of reaching ith and exceeding lower.
  probTotal <-0 #is the total probability of rejecting, sum(1,n)(probExceedingUpper+probExceedingLower)
  probDifference<-0 # probDifference(i) = probExit(i)-probExit(i-1) where probExit is a vector of exit probabilities.
  probExit<-0 # probExit is a vector of exit probabilities


  #parameter variables
  OneOrTwoSidedBounds<-1 # One-(OneOrTwoSidedBounds==1) or Two-sided bounds(OneOrTwoSidedBounds==2)
  numberOfDrifts <- 1 #number of drift parameters
  confidenceLevel<-0 #desired power respectively confidence level
  nMax<-25 #Number of interim analyses is limited to 25.
  taskNumber<-0  # number of task the user chooses in interactive dialogue

  #status variables
  t2inUse<-FALSE # second time scale will be used (t2inUse==TRUE) or not used (t2inUse==FALSE(default))
  equallySpacedTimes<-TRUE #variable denoting whether equally spaced times (equallySpacedTimes==TRUE) or different spaced (equallySpacedTimes==FALSE)
  SpendingFunctionIsUsed<-TRUE #denotes whether spending function is used ( SpendingFunctionIsUsed==TRUE ) or not used ( SpendingFunctionIsUsed==FALSE)
  secondTimeScaleIsUsed<-FALSE # second time/information scale (yes=>secondTimeScaleIsUsed==TRUE, no=>secondTimeScaleIsUsed==FALSE)
  symmetricBounds<-FALSE # symmetric bounds (symmetricBounds==TRUE)or asymmetric bounds(symmetricBounds==FALSE)
  driftParametersAreUsed <-FALSE # are drift-parameters used? (refers to Option(3)


###########################################################################
###############__  ___  __  _        _       __       ___##################
############## |_) |_  / _  | |\ |   |  |\ | |_) |  |  |  #################
############## |_) |__ \_/  | | \|   |  | \| |   |__|  |  #################
###############                                          ##################
###########################################################################
  cat("############################################################", "\n")
  cat("#                                                          #", "\n")
  cat("# This function performs computations related to           #", "\n")
  cat("#           group sequential boundaries                    #", "\n")
  cat("#                                                          #", "\n")
  cat("############################################################", "\n")
  cat("\n")
  cat("\n")
  cat("\n")
  cat("\n")

###########################################################################
###########################  CHOOSE TASK ##################################
###########################################################################

  inputCorrect <- FALSE
  while( !inputCorrect )
  {
    cat(">> Enter number for your Task <<", "\n" )
    cat("", "\n" )
    cat("-1- Compute bounds.", "\n")
    cat("-2- Compute drift given power and bounds.", "\n" )
    cat("-3- Compute probabilities given bounds and drift.", "\n" )
    cat("-4- Compute confidence interval.", "\n" )

    # read the input from console
    taskNumber <-scan(file="", n=1, quiet=TRUE)
    cat("\n")
    cat("\n")

    # Execute the chosen option...

    ### Option 1 ###
    if(length(taskNumber)==0)
    {
      taskNumber<-1
    }
    if (taskNumber==1)
    {
      cat("   ----------   ","\n")
      cat("<< Option -1- >>","\n")
      cat("   ----------   ","\n")
      cat("You will be prompted for a spending function.", "\n")
      inputCorrect <- TRUE
    }

    ### Option 2 ###
    else if (taskNumber==2)
         {
           cat("<< Option -2- >>","\n")
           cat("You will be prompted for bounds and a power level", "\n")
   inputCorrect <- TRUE
         }

         ### Option 3 ###
         else if (taskNumber==3)
              {
                cat("   ----------   ","\n")
                cat("<< Option -3- >>","\n")
                cat("   ----------   ","\n")
        cat("You will be prompted for bounds or a spending function", "\n")
        cat("to compute the bounds and for drift parameters, if desired.", "\n")
        inputCorrect <- TRUE
      }

      ### Option 4 ###
      else if (taskNumber==4)
              {
                cat("   ----------   ","\n")
                cat("<< Option -4- >>","\n")
                cat("   ----------   ","\n")
                cat("You will be prompted for bounds and a confidence level.", "\n")
inputCorrect <- TRUE
              }

      ### Invalid Option ###
              else
      {
 cat("   !!!!!!!!!!!!!!", "\n")
                cat("<< INVALID Option >> please try again!", "\n")
                cat("   !!!!!!!!!!!!!!", "\n")
                cat("\n")
   cat("\n")
cat("\n")

              }

  }#end <--*while( !inputCorrect )*




###########################################################################
########################  INTERIM ANALYSIS TIMES  #########################
###########################################################################


  ### How many interim analyses? ###
  numberCorrect<-FALSE
  while (!numberCorrect)
  {
    cat("\n")
    cat("\n")
    cat("###########################", "\n")
    cat("#                         #", "\n")
    cat("# Choosing Analysis Times #", "\n")
    cat("#                         #", "\n")
    cat("###########################", "\n")
    cat("\n")
    cat("\n")
    cat(">> Enter Number of interim analyses! <<", "\n")
    n <-scan(file="", n=1, quiet=TRUE)

    if(length(n)==0)
    {
      n<-1
    }

    ### User typed in a negative digit or 0 ###
    if (n<1)
    {
      cat("   !!!!!!!!!!!!!!!!!!!!!!!!", "\n")
      cat("<< Number must be positive! >> please try again!", "\n")
      cat("   !!!!!!!!!!!!!!!!!!!!!!!!", "\n")
    }

    ### User typed in a digit greater than nMax (default = 25)
    else if (n>nMax)
         {
         cat(" You typed in ",n," but Number must be <= ",nMax, "\n")
         cat("Please try again!", "\n")
         }

         ### User typed in an adequate Number ###
         else
         {
           numberCorrect <- TRUE
         }

  }#end <--*while*


  ### Prompt for equally spaced analyses ###
  inputCorrect<-FALSE
  while (!inputCorrect)
  {
    cat("\n")
    cat("\n")
    cat(">> Do you want equally spaced times between 0 and 1 (1=yes, 0=no)? <<", "\n")
    response <-scan(file="", n=1, quiet=TRUE)

    #handle return without any input
    if(length(response)==0)
    {
      response<-1
    }

    ### User wants equally spaced times ###
    ## Compute n equally spaced times from 1/n to n
    if (response==1)
    {
      equallySpacedTimes<-TRUE
      for(i in 1:n)
      {
        t[i]<- i/n
      }
      cat("   -------------------------------  ","\n")
      cat("<< Analysis times - equally spaced >>", "\n")
      cat("   -------------------------------  ","\n")
      cat("   ",format(t,digits=3),"\n")
      cat("\n")
      inputCorrect<-TRUE
    }

    ###--- User wants to choose his own times ---###
    ## Enter one time for each look
    else if (response==0)
         {
           equallySpacedTimes<-FALSE
           cat("\n")
           cat(">> Enter your times of interim analyses <<", "\n")
           cat("   one for each of the",n,"looks -", "\n")
           cat("   and each must be >0 and <=1 and strictly increasing:", "\n")

   #Get users interim times...
   t<-scan(file="", n=n, quiet=TRUE)

           #...and check them to be in intervall(0,1] and in right order
           interimTimesBad<-FALSE
           for (i in 1:n)
           {
             if (t[i]<=0) { interimTimesBad<-TRUE }
             if (t[i]>1) { interimTimesBad<-TRUE }
             if (i>1)
             {
               if (t[i]<=t[i-1]) { interimTimesBad<-TRUE }
             }
   }#end for

            #The entered Times were OK
            if (!interimTimesBad)
            {
              cat("   --------------------------  ","\n")
              cat("<< Analysis times you entered >>", "\n")
              cat("   --------------------------  ","\n")
              cat("   ",format(t,digits=3),"\n")
              inputCorrect<-TRUE
            }
            else
            {
              #Bad Times were entered - user has to re-enter
              cat("\n")
       cat("\n")
      cat("   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ", "\n")
              cat("<< Times out of range or order is bad >>", "\n")
      cat("   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   ", "\n")
              cat("\n")
              cat("Please re-enter or choose equally spaced times", "\n")
            }
         }#end <--*else if (response==0)*

         ### Invalid Option ###
         else
 {
           cat(" Invalid Option - please try again!", "\n")
           cat("\n")
         }

  }#end <--*while*



  ### Second time scale equal to first by default.###
  t2<-t

  ##user could use a second time/information scale
  if (!taskNumber==1)
  {
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      ###--- ask user for second time/information scale ---###
      cat("\n")
      cat("\n")
      cat("Do you wish to specify a SECOND TIME INFORMATION SCALE?", "\n")
      cat("e.g. the inverse of parameter variance or number of events,","\n")
      cat("as in Lan & DeMets 89?) (1=yes, 0=no)", "\n")
      response <-scan(file="", n=1, quiet=TRUE)

      if(length(response)==0)
      {
        response<-1
      }
      ##user wants second time/information scale - so prompt him for it
      if(response==1)
      {
        secondTimeScaleIsUsed<-TRUE
        cat("\n")
        cat("\n")
        cat("   First time scale will be used in the spending function.", "\n")
        cat("   Second time scale will estimate covariances.", "\n")
        cat(">> Enter Information scale <<", "\n")
        cat("   (Information must be strictly increasing and >0)", "\n")


        #seconds information/time scale is saved in 't2'
        t2 <-scan(file="", n=n, quiet=TRUE)

        ##verify range and ordering of information
        InformationBad<-FALSE

        for (i in 1:n)
        {
          if (t2[i]<=0) { InformationBad<-TRUE }
          if (i>1)
          {
            if (t2[i]<=t2[i-1]) { InformationBad<-TRUE }
          }
}#end for

        #Information entries were OK
        if (!InformationBad)
        {
          cat("Your Information choice:", format(t2,digits=3),"\n")
          cat("\n")
          inputCorrect<-TRUE
        }

        #Information entries were wrong -
        #user could try againg or choose not to use second time scale
        else
        {
          cat("\n")
          cat("\n")
          cat("   !!!!!!!!!!!!!!!!!!!!!!!   ", "\n")
          cat("<< Your Entries were wrong >>", "\n")
          cat("   !!!!!!!!!!!!!!!!!!!!!!!   ", "\n")
          cat("\n")
          cat("Information must be strictly increasing and >0!","\n")
          cat("\n")
          inputCorrect<-FALSE
        }

      }#end <--*if(response==1)*


      ##user does NOT want a second time/information scale
      else if(response==0)
           {
             secondTimeScaleIsUsed<-FALSE
             #Reset Second time scale equal to first as was by default.
             t2<-t
             inputCorrect<-TRUE
           }

           #Invalid option
           else
           {
             cat("!!! Invalid option - please try again !!!","\n")
             cat("\n")
           }

    }#end <--*while(!inputCorrect)*
  }#end <--*(if (!taskNumber==1)*


###########################################################################
###########################  DETERMINE BOUNDS  ############################
###########################################################################
# In Task 1, user has to specify the spending function to compute the bounds
#
# In Tasks 2, 3 and 4 user decides whether to use spending fct.
# or instead of it enter a set of bounds manually.

###########################################################################
########################  TASK 1 WAS CHOSEN  ##############################
############## COMPUTE BOUNDS USING A SPENDING FUNCTION ###################
###########################################################################
  # Task 1 was chosen
  if (taskNumber==1)
  {
    SpendingFunctionIsUsed<-TRUE
  }


  # Task 2,3 or 4 was chosen
  else
  {
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      cat("\n")
      cat("\n")
      cat("    Do you want to use a spending function to determine bounds", "\n")
      cat("    or enter the bounds manually?:", "\n")
      cat(">> (1=yes,use spending fkt, 0=no,enter bounds manually) <<", "\n")
      response <-scan(file="", n=1, quiet=TRUE)

      if(length(response)==0)
      {
        response<-1
      }
      #user wants spending function
      if (response==1)
      {
        SpendingFunctionIsUsed<-TRUE
        inputCorrect<-TRUE
      }
      else if (response==0)
           {
             SpendingFunctionIsUsed<-FALSE
             inputCorrect<-TRUE
            }
           else
           {
             cat("!!! Invalid response - please try again !!!", "\n" )
             cat("\n")
           }
    }#end <--*while*
  }

  ## User wants a spending function to compute the bounds
  ##-----------------------------------------------------------##
  ##--------  DETERMINE BOUNDS USING SPENDING FUNCTION --------##
  ##-----------------------------------------------------------##

  if (SpendingFunctionIsUsed)
  {
    VectorOfBounds <- useSpendingFunction(t, n, t2)
    alpha <- VectorOfBounds[[1]]
    lowerBounds <- VectorOfBounds[[2]]
    upperBounds <- VectorOfBounds[[3]]
    probExit <- VectorOfBounds[[4]]
    probDifference <- VectorOfBounds[[5]]
    symmetricBounds <- VectorOfBounds[[6]]
    spendingFunctionUsed <- VectorOfBounds[[7]]
    OneOrTwoSidedBounds <- VectorOfBounds[[8]]

    ##output bounds
    outputBounds(n,alpha,t,lowerBounds,upperBounds,probDifference,probExit,symmetricBounds,spendingFunctionUsed)
  }


  ##-----------------------------------------------------------##
  ##-----------  USER ENTERS BOUNDS BOUNDS MANUALLY -----------##
  ##-----------------------------------------------------------##
  else
  {

    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      ##User enters choice for one or two sided bounds.##
      cat("\n")
      cat("\n")
      cat(">> Do you want One(1) or two(2)-sided bounds? <<", "\n")
      OneOrTwoSidedBounds <-scan(file="", n=1, quiet=TRUE)

      if(length(OneOrTwoSidedBounds)==0)
      {
        OneOrTwoSidedBounds<-1
      }
      ##One-sided##
      if(OneOrTwoSidedBounds==1)
      {
        cat("\n")
        cat("\n")
        cat("##################", "\n")
        cat("#                #", "\n")
        cat("# One-sided Test #", "\n")
        cat("#                #", "\n")
        cat("##################", "\n")
        cat("\n")
        cat("\n")

        inputCorrect<-TRUE
      }
      ##Two-sided##
      else if(OneOrTwoSidedBounds==2)
           {
             cat("\n")
             cat("\n")
             cat("##################", "\n")
             cat("#                #", "\n")
             cat("# Two-sided Test #", "\n")
             cat("#                #", "\n")
             cat("##################", "\n")
             cat("\n")
             cat("\n")

             inputCorrect<-TRUE

             ##Two-sided bounds may be symmetric or asymmetric - prompt for it
             symmetricInputCorrect<-FALSE
             while (!symmetricInputCorrect)
             {
               cat(">> Do you want Symmetric bounds? (1=yes, 0=no) <<","\n")
               response <-scan(file="", n=1, quiet=TRUE)

               #handle return without any input
               if(length(response)==0)
             {
                 response<-1
               }

               ##Symmetric##
               if(response==1)
               {
                 symmetricBounds<-TRUE
                 cat("\n")
                 cat("\n")
                 cat("<< Two-sided Symmetric Bounds will be used >>", "\n")
 cat("   ----------^----------------------------", "\n")
                 cat("\n")

               symmetricInputCorrect<-TRUE
               }

               ##Asymmetric##
               else if(response==0)
                    {
                      symmetricBounds<-FALSE
                      cat("\n")
                      cat("\n")
                      cat("<< Two-sided Asymmetric Bounds will be used >>", "\n")
                      cat("   ----------^-----------------------------", "\n")
                      cat("\n")

                      symmetricInputCorrect<-TRUE
                    }

                    ##Invalid input
                    else
                    {
                      cat("\n")
                      cat("!!! Invalid option - please try again !!!","\n")
                      cat("\n")
                    }

             }#end <--*while (!symmetricInputCorrect)*

           }#end <--*else if(OneOrTwoSidedBounds==2)*

           ##Invalid input
           else
           {
             cat("\n")
             cat("!!! Invalid option - please try again !!!","\n")
             cat("\n")
           }

    }#end <--*while (!inputCorrect)*





    ### Now user enters the Bounds###
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      cat(">> Enter upper Bounds (standarized) <<", "\n")
      upperBounds <-scan(file="", n=n, quiet=TRUE)

      ## case One-sided ##
      ##the lower bounds are set to minus infinity respectively to -8
      if(OneOrTwoSidedBounds==1)
      {
        lowerBounds <- seq(-8,-8,length=n)
      }

      ## case Two-sided ##
      else if(OneOrTwoSidedBounds==2)
           {
             ##case Symmetric - lower bounds are negative upper bounds##
             if(symmetricBounds)
             {
               lowerBounds <- -upperBounds
             }

             ##case Asymmetric - user also enters lower Bounds##
             else
             {
               cat("\n")
               cat(">> Now Enter lower Bounds (standarized) <<", "\n")
          lowerBounds <-scan(file="", n=n, quiet=TRUE)
             }
      }#end <--*if*

      times <- data.frame(" *Times*"=t,check.names=FALSE)
      bounds <- data.frame(" *Lower Bounds*"=lowerBounds," *Upper Bounds*"=upperBounds,check.names=FALSE)

      cat("\n")
      cat("\n")
      cat(" Your entered Bounds are:","\n")
      cat(" ------------------------","\n")

      print( cbind( format(times,digits=3),format(bounds,digits=5) ) )

      cat(" ------------------------------------------------","\n")
      cat("\n")
      cat(" Use this Bounds, or Re-Enter? (1=yes,0=Re-Enter)","\n")
      yesNo <-scan(file="", n=1, quiet=TRUE)

      #handle a return instead of character input plus return
      if(length(yesNo)==0)
      {
        yesNo<-1
      }

      if(yesNo==1)
      {
        inputCorrect<-TRUE
      }
      else
      {
        ##Re-Enter Bounds##
      }

    }#end <--*while (!inputCorrect)*
  }#end <--  *else 'USER ENTERS BOUNDS MANUALLY'*


###########################################################################
########################  TASK 2 WAS CHOSEN  ##############################
################# COMPUTE DRIFT GIVEN POWER AND BOUNDS#####################
###########################################################################

  ## Default is one drift equal to zero
  numberOfDrifts<-1
  drift[1]<-0

  if (taskNumber==2)
  {
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      ## Enter desired power for Option(2).
      cat("\n")
      cat(">> Enter Desired Power - it must be in (0,1) : <<", "\n")
      confidenceLevel <-scan(file="", n=1, quiet=TRUE)

      ##Range check
      if( !(confidenceLevel>0 & confidenceLevel<1) )
      {
        cat("!!! Invalid response! Confidence Level must be > 0 and < 1 !!!","\n")
        cat(" Please try again.","\n")
      }
      else
      {
        inputCorrect<-TRUE
      }
    }#end <--*while*

  }#end <--*if*





###########################################################################
########################  TASK 3 WAS CHOSEN  ##############################
############ Compute probabilities given bounds and drift.#################
###########################################################################

  if (taskNumber==3)
  {
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      ## For Option(3) are drift parameters to be used?
      cat("\n")
      cat(">> Do you wish to use drift parameters? (1=yes, 0=n0) <<","\n")
      response <-scan(file="", n=1, quiet=TRUE)

      if(length(response)==0)
      {
        response<-1
      }

      ##verify input
      if(response==1)
      {
        driftParametersAreUsed<-TRUE
        inputCorrect<-TRUE
        numberOfDriftParameterCorrect<-FALSE
        while (!numberOfDriftParameterCorrect)
        {
          cat("\n")
          cat("\n")
          cat(">> How many drift parameters do you wish to enter? <<","\n")
          numberOfDrifts <-scan(file="", n=1, quiet=TRUE)

          if(length(numberOfDrifts)==0)
          {
            numberOfDrifts<-1
          }

          ##number of drifts must be >0 and <= nMax (=maximum number of interim analyses)
          if(numberOfDrifts<=0)
          {
            cat("\n")
            cat("!!! Number of drifts must be positive! Try again !!!","\n")
            cat("\n")
          }
          else if(numberOfDrifts>nMax)
               {
                 cat("\n")
                 cat(" Number has to be smaller than maximum number","\n")
                 cat(" of interim analyses, which is by default:",nMax,"\n")
                 cat(" or the maximum number must be increased in function'input.R'. Try again!","\n")
                 cat("\n")
               }
               else
               {
                 numberOfDriftParameterCorrect<-TRUE
               }
        }#end <--*while (!driftParameterInputCorrect)*

        ##-- user now enters the drift parameters --##
        cat("\n")
        cat(">> Enter drift parameters <<","\n")
        drift <-scan(file="", n=numberOfDrifts, quiet=TRUE)

        cat("\n")
        cat(" Drift is equal to the expectation of the Z statistic when time=1","\n")

      }#end <--*if(response==1)*

      ##--User choosed to use no drift parameters--##
      else if(response==0)
           {
             driftParametersAreUsed<-FALSE
             cat("\n")
             cat("   -------------------------------------","\n")
             cat("<< You choose to use no drift parameters >>","\n")
             cat("   -------------------------------------","\n")
             cat("\n")
             inputCorrect<-TRUE
           }
           ##Invalid Option
           else
           {
             cat("\n")
             cat("!!! Invalid Option. Please try again !!!","\n")
             cat("\n")
           }
    }#end <--*while (!inputCorrect)*
  }#end <--*if (taskNumber==3)*


###########################################################################
########################  TASK 4 WAS CHOSEN  ##############################
#################### Compute confidence intervall.#########################
###########################################################################

  if (taskNumber==4)
  {
    ##-- Level for Confidence Intervall --##
    ## confidence Intervall replaces last bound with Zvalue
    cat("\n")
    cat(">> Enter the standarized statistic (Z value) at the last analysis <<","\n")
    Zvalue <-scan(file="", n=1, quiet=TRUE)

    ## prompt for confidence level
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      cat("\n")
      cat(">> Enter confidence Level - it must be > 0 and < 1 <<","\n")
      confidenceLevel <-scan(file="", n=1, quiet=TRUE)

      ##Range check
      if( !(confidenceLevel>0 & confidenceLevel<1) )
      {
        cat("!!! Invalid response! Confidence Level must be > 0 and < 1 !!!","\n")
        cat(" Please try again.","\n")
      }
      else
      {
        inputCorrect<-TRUE
        cat("\n")
        cat("<< Confidence Level in percent: ",confidenceLevel*100,"% >>","\n")
        cat("\n")
      }

    }#end <--*while (!inputCorrect)*

  }#end <--*if (taskNumber==4)*


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
  t3<-t2

# If second time scale has been entered.
  if(secondTimeScaleIsUsed)
  {
    t2max<-0

    ##Options using non-zero drifts:
    if(taskNumber==2 || taskNumber==4 || driftParametersAreUsed)
    {
      ##If t[n]=1, t2[n] is maximum of t2.
      if(t[n]==1)
      {
        cat("\n")
        cat("<< 2nd time scale will be used to determine covariances >>","\n")
        cat("\n")
        t2max<-t2[n]
        t3<-t2/t2max
      }
      else
      {
        ##Should we try to use 2nd scale?
        inputCorrect<-FALSE
        while (!inputCorrect)
        {
          cat("\n")
          cat("\n")
          cat("   Do you wish to use the 2nd time scale","\n")
          cat(">> to determine covariances? (1=yes,0=no) <<","\n")
          response <-scan(file="", n=1, quiet=TRUE)

          if(length(response)==0)
          {
            response<-1
          }

          ##--If yes, prompt for maximum of 2nd time scale.--##
          if(response==1)
          {
            t2inUse<-TRUE
            inputCorrect<-TRUE

            maximumCorrect<-FALSE
            while (!maximumCorrect)
            {
              cat("\n")
              cat(">> Enter the maximum value of the second time scale <<","\n")
              t2max <-scan(file="", n=1, quiet=TRUE)

              ##check input
              if(t2max<=0)
              {
                cat("\n")
                cat("!!! Maximum must be positive! Please try again !!!","\n")
                cat("\n")
              }
              if(t2max<t2[n])
              {
                cat("\n")
                cat("!!! This value has already been exceeded - Please try again !!!","\n")
                cat("\n")
              }
              else
              {
                maximumCorrect<-TRUE

                ##Rescale t2
                t3<-t2/t2max
              }

            }#end <--*while (!maximumCorrect)*
          }#end <--*if(response==1)*

          else if(response==0)
               {
                 t2inUse<-FALSE
                 ##Even if the 2nd time scale was on (0,1], if the
                 ##maximum information value was not entered, set
                 ##t3 to t, which causes t2 to be ignored!

                 inputCorrect<-TRUE
                 t3<-t
               }
               else
               {
                 cat("!!! Invalid response! Pleasy try again !!!","\n")
               }


        }#end <--*while (!inputCorrect)*

      }#end <--*else*

    }#end <--*if(taskNumber==2 || taskNumber==4 || driftParametersAreUsed)*
  }#end <--*if(secondTimeScaleIsUsed)*






###########################################################################
############___      __     __  ___   _       __       ___#################
########### |_  |\ | | \   /  \ |_    |  |\ | |_) |  |  |  ################
########### |__ | \| |_/   \__/ |     |  | \| |   |__|  |  ################
############                                              #################
###########################################################################
  if(taskNumber==1)
  {
    cat("\n")
    cat("\n")
    cat(" Type 'g' to see a graph or anything else to quit", "\n")
    answer <-scan(file="",what='character',n=1, quiet=TRUE)

    if(length(answer)==0)
    {
      answer=0
    }
    if(answer=="g")
    {
      ## if one-Sided-Test we won't see negative Z-Values
      if(OneOrTwoSidedBounds==1)
      {
        xCoordinate<-t
        yCoordinate<-upperBounds

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main="-Option 1-  Group Sequential Boundaries",
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4,
             xlab="Times",ylab="Standarized Z-Value",ylim=c(0,4))

        ##...then add lines between them
        lines(t,upperBounds,col="blue")
      }

      else
      {
        xCoordinate<-c(t,t)
        yCoordinate<-c(lowerBounds,upperBounds)

        ## first plotting bounds as points...
        plot(xCoordinate,yCoordinate,main="-Option 1-  Group Sequential Boundaries",
             pch=21,bg="green",font=4,font.axis=4,font.lab=4,font.main=4,
             xlab="Times",ylab="Standarized Z-Value",ylim=c(-4,4))

        ##...then add lines between them
        lines(t,lowerBounds,col="blue")
        lines(t,upperBounds,col="blue")
      }
    }
    else
    {
      #quit
      cat("quit \n")
    }
  }

  else if(taskNumber==2)
       {
         ##--Return drift for given bounds and power.------------------------------------##
         ##-- Note: the drift return produces exit probability confidenceLevel at t[n].--##
         ##--If t(n) < 1 (or t3(n) < 1) this is not the study's power.-------------------##

         drift <- findDrift(n, t3, t2, lowerBounds, upperBounds, confidenceLevel, drift, nMax)

         ##check whether drift was calculated correctly
         if(is.numeric(drift))
         {
           ##we got result - output drift
           vectorOfResults <- computeAlphaLevel(n, t3, t2, lowerBounds, upperBounds, drift, nMax)
           probStopping <- vectorOfResults[[1]]
           probExceedingUpper <- vectorOfResults[[2]]
           probExceedingLower <- vectorOfResults[[3]]
           expectedStoppingTime <- vectorOfResults[[4]]
           probTotal <- vectorOfResults[[5]]

           ##Print out results from function 'computeAlphaLevel'
           outputDriftWithBounds(n,probTotal,drift,expectedStoppingTime,secondTimeScaleIsUsed,t,t2,t2max,lowerBounds,upperBounds,probStopping,probExceedingUpper,probExceedingLower, confidenceLevel)
         }
         else
         {
           ##something went wrong
           cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!","\n")
           cat("! ERROR WhiLE COMPUTING THE DRIFT !","\n")
           cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!","\n")
         }
       }

       else if(taskNumber==3)
       {
         ##-----------------------------------------------------------##
         ##--Probabilities from bounds, possibly with non zero drift--##
         ##-----------------------------------------------------------##

         cat("\n")
         cat("\n")
         cat("#############################", "\n")
         cat("#                           #", "\n")
         cat("# Output Drifts with Bounds #", "\n")
         cat("#                           #", "\n")
         cat("#############################", "\n")
         cat("\n")
         cat("\n")
         cat("*------------------*","\n")
         cat("*Drift parameter(s):",format(drift,digits=5),"\n")
         cat("*------------------*","\n")
         cat("Drift is equal to the expectation of the Z statistic when time=1.","\n")


         for(i in 1:numberOfDrifts)
         {
           vectorOfResults <- computeAlphaLevel(n,t3,t2,lowerBounds,upperBounds,drift[i],nMax)
           probStopping <- vectorOfResults[[1]]
           probExceedingUpper <- vectorOfResults[[2]]
           probExceedingLower <- vectorOfResults[[3]]
           expectedStoppingTime <- vectorOfResults[[4]]
           probTotal <- vectorOfResults[[5]]

           ##Print out results from function 'computeAlphaLevel'
           cumulativeExitProb <- outputForDifferentDrifts(n,probTotal,drift[i],expectedStoppingTime,secondTimeScaleIsUsed,t,t2,t2max,lowerBounds,upperBounds,probStopping,probExceedingUpper,probExceedingLower)
           power[i] <- cumulativeExitProb
         }

       }#end <--*else if(taskNumber==3)*

       else if(taskNumber==4)
            {
              ##-----------------------------------------------------------##
              ##-Confidence limits from bounds and final statistics value.-##
              ##-----------------------------------------------------------##
              confidenceIntervall <- computeConfidenceIntervall(confidenceLevel,Zvalue,n,t3,t2,lowerBounds,upperBounds,nMax)
              cat("*--------------------*","\n")
              cat("*Confidence Intervall:  <",confidenceIntervall[1]," , ",confidenceIntervall[2],">","\n")
              cat("*--------------------*","\n")
              cat(" Drift is equal to the expectation of the Z statistic when time=1.","\n")
            }


}#end <--*function(x)*
