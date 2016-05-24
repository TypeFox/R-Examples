"useSpendingFunction" <-
function(t, n, t2)
{
##################################################################
#################### INITIALIZE VARIABLES ########################
##################################################################
  probExit<-0 # probExit is a vector of exit probabilities
  probDifference<-0 # probDifference(i) = probExit(i)-probExit(i-1) where probExit is a vector of exit probabilities.
  lowerBounds <-0 # lowerBounds is the vector of lower standardized boundaries
  upperBounds <-0 # upperBounds is the vector of upper standardized boundaries
  negInf<-(-8) # negInf is "negative infinity" for the program.
  boundsTruncation<-abs(negInf) # boundsTruncation is the user selected truncation on integration limits (default is abs(negInf) )
  boundTruncationIsUsed<-FALSE # Truncation point used for bounds?
  symmetricBounds<-FALSE # denotes whether symmetric bounds are used or not
  OneOrTwoSidedBounds<-1 # One-(OneOrTwoSidedBounds==1) or Two-sided bounds(OneOrTwoSidedBounds==2)
  whatSpendingFunctionIsUsed<-0 # denotes the spending function which is used, (e.g. Pocock Type=> whatSpendingFunctionIsUsed<-2, etc...)
  phi<-c(1,1) #parameters for  Power family (spending function) -default is (1,1)
  alpha<-0 # significance level alpha


    cat("", "\n")
    cat("", "\n")
    cat("############################################", "\n")
    cat("#                                          #", "\n")
    cat("# Compute Bounds using a Spending Function #", "\n")
    cat("#                                          #", "\n")
    cat("############################################", "\n")
    cat("", "\n")
    cat("", "\n")
    ##-- Total type I error probability --##
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      cat(">> Enter overall significance level ( >0 and <=1 ) <<", "\n")
      alpha<-scan(file="", n=1, quiet=TRUE)

      if(length(alpha)==0)
      {
        alpha[1]<- 0
      }
      #correct input
      if ( alpha[1]>0 && alpha[1]<=1 )
      {
        cat("", "\n")
        cat("<< You set alpha to",alpha[1],">>","\n")
        cat("", "\n")
        inputCorrect<-TRUE
      }
      #wrong alpha
      else
      {
        cat("", "\n")
        cat("   !!!!!!!!!!!!!!!!", "\n")
        cat("<< INVALID response >>", "\n")
        cat("   !!!!!!!!!!!!!!!!", "\n")
        cat("\n")
        cat("   alpha must be >0 and <=1! Try again!","\n")
        cat("\n")
      }

    }#end <--*while*


    ###---Upper and lower bounds, upper bounds only, or asymmetric? ---###
    inputCorrect<-FALSE
    while (!inputCorrect)
    {
      cat("\n")
      cat("\n")
      cat(">> One(1), two(2)-sided symmetric, or asymmetric(3) bounds? <<","\n")
      boundsOption<-scan(file="", n=1, quiet=TRUE)

      if(length(boundsOption)==0)
      {
        boundsOption<- 1
      }

      ##verify correct input and set status variables for further computation of the bounds
      if (boundsOption==1)
      {
        symmetricBounds<-TRUE
        OneOrTwoSidedBounds<-1
        inputCorrect<-TRUE
      }
      else if (boundsOption==2)
           {
             symmetricBounds<-TRUE
             OneOrTwoSidedBounds<-2
             inputCorrect<-TRUE
           }
           else if (boundsOption==3)
           {
             symmetricBounds<-FALSE
             OneOrTwoSidedBounds<-2
             inputCorrect<-TRUE
           }

           #Invalid response
           else
           {
             cat("", "\n")
             cat("   !!!!!!!!!!!!!!!!", "\n")
             cat("<< INVALID response >>", "\n")
             cat("   !!!!!!!!!!!!!!!!", "\n")
             cat("\n")
             cat("   Enter '1', '2' or '3' - try again!","\n")
             cat("", "\n")
           }

    }#end <--*while*



    ###symmetric bounds, so alpha[2] and whatSpendingFunctionIsUsed[2] are unused - set null values
    if (symmetricBounds)
    {
      alpha[2]<-0
      whatSpendingFunctionIsUsed[2]<-0

      #just print what user choosed
      if (OneOrTwoSidedBounds==1)
      {
        cat("", "\n")
        cat(" ##################", "\n")
        cat(" #                #", "\n")
        cat(" # One-sided Test #", "\n")
        cat(" #                #", "\n")
        cat(" ##################", "\n")
        cat("", "\n")
      }

      else #OneOrTwoSidedBounds==2
      {
        cat("", "\n")
        cat(" ##################", "\n")
        cat(" #                #", "\n")
        cat(" # Two-sided Test #", "\n")
        cat(" #                #", "\n")
        cat(" ##################", "\n")
        cat("", "\n")
      }

    }#end *if*


    ###Asymmetric bounds - prompt user for Type I error for lower bound
    if (!symmetricBounds)
    {
      cat("", "\n")
      cat(" #############################", "\n")
      cat(" #                           #", "\n")
      cat(" # Two-sided asymmetric Test #", "\n")
      cat(" #                           #", "\n")
      cat(" #############################", "\n")
      cat("", "\n")
      cat(" Note: bounds will be approximately only.","\n")
      cat(" Nominal exit probs may be slightly too large.","\n")

      inputCorrect<-FALSE
      while (!inputCorrect)
      {
        ##get alpha for lower bound
        cat("", "\n")
        cat(" >> Enter alpha for LOWER bound <<","\n")
        cat("    consider: 0 < your_alpha_for_lower_bound <",alpha[1],"\n")
        cat("", "\n")
        alpha[2]<-scan(file="", n=1, quiet=TRUE)

        ##input correct? that is if 0 < alpha[2] < alpha[1]
        if (alpha[2]>0 && alpha[2]<alpha[1])
        {
          #compute upper alpha by subtraction
          alpha[1] <- alpha[1] - alpha[2]

          #output upper and lower alpha
          cat("\n")
          cat("<< Upper alpha:",alpha[1],">>","\n")
          cat("<< Lower alpha:",alpha[2],">>","\n")
          cat("\n")
          cat("\n")

          inputCorrect<-TRUE
        }
        else
        {
          cat("", "\n")
          cat("   !!!!!!!!!!!!!", "\n")
          cat("<< INVALID input >>", "\n")
          cat("   !!!!!!!!!!!!!", "\n")
          cat("\n")
          cat("   Lower alpha out of range! - Try again!","\n")
          cat("\n")
        }
      }#end <--*while*

    }#end *if(!symmetricBounds)*



    ###--- Selection of use function from programmed choices ---###
    ## If asymmetric bounds we need TWO functions!

    ####################
    ### Outer For-Loop ##################################
    ## uses the inner Loop twice to select in each case ####
    # a function to get TWO function for asymmetric bounds #

    LoopsFirstTime<-TRUE

    for(i in 1:2)
    {
      if (LoopsFirstTime)
      {
        if (symmetricBounds)
        {
          #symmetric bounds - enter THE function
          cat(">> Enter use function <<", "\n")
        }
        else
        {
          #asymmetric bounds - enter function for lower bound first
          cat(">> Enter use function for lower bound <<","\n")
        }
      }

      ##if this is second loop and we got symmetric bounds - do nothing else
      ##otherwise save the inputs for the functions for lower bounds
      ##and prompt for the functions for upper bounds
      else if(!LoopsFirstTime && symmetricBounds)
           {
             #do nothing but go into while-loop and do nothing else there too
           }

           else ##prompt for 2nd function
           {
             #set spending function2 -> whatSpendingFunctionIsUsed[2] and parameter2 -> phi[2] ...
             whatSpendingFunctionIsUsed[2]<-whatSpendingFunctionIsUsed[1]
             phi[2]<-phi[1]

             #... and go on with the Inner Loop for second (and last) time
             cat("", "\n")
             cat(" Now enter use function for upper bound", "\n")
           }

      ######################
      ### Inner While-Loop #############################
      ## Symmetric bounds - user chooses ONE function: #

      #check whether we got Symmetric Bounds and whether they are computed yet -
      #cause then we do not need a second input
      if (!LoopsFirstTime && symmetricBounds)
      {
        #do nothing
      }

      else #we have to go on
      {
        LoopsFirstTime<-FALSE
        useFunctionIsCorrect<-FALSE
        while (!useFunctionIsCorrect)
        {
          cat("\n")
          cat(" What function should be used? (1-5)", "\n")
          cat("\n")
          cat(" Approximated by spending function:", "\n")
          cat(" ----------------------------------", "\n")
          cat(" (1) O'Brien-Fleming Type.", "\n")
          cat(" (2) Pocock Type.", "\n")
          cat(" (3) Power family: alpha* t^phi.", "\n")
          cat(" (4) Hwang-Shih-DeCani fammily.", "\n")
          cat("\n")
          cat(" Classic Design:", "\n")
          cat(" ---------------", "\n")
          cat(" (5) Pocock - the real Pocock Bounds","\n")
          cat("\n")
          cat("\n")

          response<-scan(file="", n=1, quiet=TRUE)
          if(length(response)==0)
          {
            whatSpendingFunctionIsUsed[1]<- 1
          }
          else
          {
            whatSpendingFunctionIsUsed[1]<- response
          }

          ##(1),(2) or (5) was chosen - input correct
          if (whatSpendingFunctionIsUsed[1]==1 || whatSpendingFunctionIsUsed[1]==2  || whatSpendingFunctionIsUsed[1]==5)
          {
            useFunctionIsCorrect<-TRUE
          }

          ##(3) or (4) was chosen
          ##Enter parameter for spending function family if needed
          else if (whatSpendingFunctionIsUsed[1]==3)
               {
                 powerFamilyParamIsCorrect<-FALSE
                 while (!powerFamilyParamIsCorrect)
                 {
                   cat("\n")
                   cat(">> Enter exponent for Power family (phi > 0) <<","\n")
                   phi<-scan(file="", n=1, quiet=TRUE )
                   if(length(phi)==0)
                   {
                     phi[1]<- 0
                   }

                   #check wheter phi[1] is not zero
                   if (phi[1]<=0)
                   {
                     cat("\n")
                     cat("!!! Parameter must be > 0! Try again !!!","\n")
                     phi[1]<-1 #reset phi[1] to default
                   }
                   else
                   {
                     powerFamilyParamIsCorrect<-TRUE
                     useFunctionIsCorrect<-TRUE
                   }
                 }#end <--*while (!parameterCorrect)*
               }#end <--*else if (whatSpendingFunctionIsUsed[1]==3)*

               else if (whatSpendingFunctionIsUsed[1]==4)
                    {
                      hSDeCaniparamIsCorrect<-FALSE
                      while (!hSDeCaniparamIsCorrect)
                      {
                        cat("\n")
                        cat(" Enter exponent for for Hwang-Shih-DeCani family (phi =/= 0)","\n")
                        cat("\n")
                        phi<-scan(file="", n=1, quiet=TRUE )
                        if(length(phi)==0)
                        {
                          phi[1]<- 0
                        }

                        #check wheter phi[1] is not zero
                       if (phi[1]<=0)
                       {
                         cat("!!! Parameter must be nonzero! Try again !!!","\n")
                         phi[1]<-1 #reset phi[1] to default
                       }
                       else
                       {
                         hSDeCaniparamIsCorrect<-TRUE
                         useFunctionIsCorrect<-TRUE
                       }
                     }#end <--*while (!parameterCorrect)* }
                   }#end <--*if (whatSpendingFunctionIsUsed[1]==4)*

                   ### wrong input
                   else
                        {
          cat("\n")
  cat("!!! The Option you choosed does not exist - Try again !!!","\n")
                          cat("\n")

                        }
        }#end <--*while (!useFunctionIsCorrect)*


      }#end *else #we have to go on*


      ### end <--Inner While-Loop #
      ###########################


    }#end <--*for(i in 1:2)*
    ### end <--Outer For-Loop #
    ###########################
    ###########################


    inputIsCorrect<-FALSE
    while(!inputIsCorrect)
    {
      ###--- Truncation point if needed by user ---###
      cat("\n")
      cat(">> Do you want to truncate the standarized bounds? (1=yes or 'press enter'=no) <<","\n")
      response<-scan(file="", n=1, quiet=TRUE)
      if(length(response)==0)
      {
        response<- 0
      }

      if (response==1)
      {
        boundTruncationIsUsed<-TRUE
        #Prompt for truncation point
        cat("\n")
        cat(">> Enter absolute value of truncation point <<","\n")
        cat("\n")
        boundsTruncation<-scan(file="", n=1, quiet=TRUE)

        #make absolute value of truncation point
        boundsTruncation<-abs(boundsTruncation)
        cat("<< Bounds will be truncated at",boundsTruncation,">>", "\n")
        inputIsCorrect<-TRUE
      }#end <--*if*

       else if (response==0)
            {
              boundTruncationIsUsed<-FALSE
              cat("\n")
              cat("<< Bounds will be not truncated >>", "\n")
              cat("\n")
              inputIsCorrect<-TRUE
            }
            else
            {
              cat("\n")
              cat("!!! Invalid Option - please try again !!!", "\n")
              #inputIsCorrect<-FALSE
            }

    }#end <--*while(!inputIsCorrect)*



###########################################################################
########################## EVALUATE BOUNDS ################################
###########################################################################
#
    ##Symmetric bounds - call function computeBounds(...)

    if (symmetricBounds)
    {
      results<- computeBounds(n, 0, alpha[1], phi[1], t, t2, OneOrTwoSidedBounds, whatSpendingFunctionIsUsed[1], boundsTruncation)
    }

    else
    ## Asymmetric bounds - call function computeBounds(...) twice
    ## first time for upper bounds, second time for lower bounds
    {
      resultsUpperBounds <- computeBounds(n, 0, alpha[1], phi[1], t, t2, 1, whatSpendingFunctionIsUsed[1], boundsTruncation)
      resultsLowerBounds <- computeBounds(n, 0, alpha[2], phi[2], t, t2, 1, whatSpendingFunctionIsUsed[2], boundsTruncation)
    }

    ## get the values depending on one-sided or two-sided test had been made ##
    ##-- symmetric bounds --##
    if (symmetricBounds)
    {
      lowerBounds <- results[[1]]
      upperBounds <- results[[2]]
      probExit <- results[[3]]
      probDifference <- results[[4]]
    }

    ##-- asymmetric bounds --##
    else
    {
      upperBounds <- resultsUpperBounds[[2]]
      lowerBounds <- (-1)*resultsLowerBounds[[2]]
      probExit <- resultsUpperBounds[[3]] + resultsLowerBounds[[3]]
      probDifference <- resultsUpperBounds[[4]] + resultsLowerBounds[[4]]
    }


    ## if -Option 5- '(5) Pocock Type - the real Pocock Bounds' was chosen -
    ## we have to do some extra calculations
    ## The Spending function gives us an approximately Pocock-Design.
    ## To compute the exact Pocock Bounds we will do according to the following pattern:
    ## (1st)we give the bounds with all bounds are equal. As starting value we are using the
    ##      mean of the bounds computed by our Pocock spending function. I figured out that
    ##      in almost every case this is a quite good approximation so far.
    ## (2nd)we compute the probability according to our equal bounds, as we would do, if user
    ##      had chosen Task-3- at the beginning
    ## (3rd)we use Newton Iteration to adjust the bounds in every Iteration until we get the appropriate alpha

    if(whatSpendingFunctionIsUsed[1]==5 || whatSpendingFunctionIsUsed[2]==5)
    {
      ##check for symmetric bounds
      if (symmetricBounds && OneOrTwoSidedBounds==1)
      {
        upperBounds <- calculateEqualBounds(alpha[1],upperBounds,n,t2)
      }

      ##two-sided symmetric
      else if(symmetricBounds && OneOrTwoSidedBounds==2)
           {
             {
               upperBounds <- calculateEqualBounds(alpha[1]/2,upperBounds,n,t2)
               lowerBounds <- -upperBounds
             }
           }

           else
           ## asymmetric bounds -> maybe we have to calculate 2 times
           {
             ##check where -Option 5- was chosen
             if(whatSpendingFunctionIsUsed[1]==5)
             {
               upperBounds <- calculateEqualBounds(alpha[1],upperBounds,n,t2)
             }
             if(whatSpendingFunctionIsUsed[2]==5)
             {
               lowerBounds <- (-1)*calculateEqualBounds(alpha[2],-lowerBounds,n,t2)
             }
           }

      ##-----------------------------------------------------------##
      ##--Probabilities from bounds, possibly with non zero drift--##
      ##-----------------------------------------------------------##

      vectorOfResults <- computeAlphaLevel(n,t2,t2,lowerBounds,upperBounds,0,25)
      probExceedingUpper <- vectorOfResults[[2]]
      probExceedingLower <- vectorOfResults[[3]]

      ## re-compute exit probability and cumulative exit probability
      probDifference<-0
      probExit<-0

      for(i in 1:n)
      {
        probDifference[i] <- probExceedingUpper[i]+probExceedingLower[i]

        if(i==1)
        {
          probExit[i] <- probDifference[i]
        }
        else
        {
          probExit[i] <- probExit[i-1] + probDifference[i]
        }
      }
    }#end <--*if(whatSpendingFunctionIsUsed[1]==5 || whatSpendingFunctionIsUsed[2]==5)*<--


    ##Return alpha,bounds etc...
    toBeReturned<-list(alphaValue=alpha,lowerBounds=lowerBounds, upperBounds=upperBounds, exitProbs=probExit,
                       diffExitProbs=probDifference, symmetricYesNo=symmetricBounds, functionUsed=whatSpendingFunctionIsUsed, oneTwoSided=OneOrTwoSidedBounds)
    return(toBeReturned)
}#end <--*function(x)*
