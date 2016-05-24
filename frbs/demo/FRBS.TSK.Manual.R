 ## This example shows how to use frbs without 
 ## learning process.
 ## Note: some variables might be shared for other examples.

 ## Define shape and parameters of membership functions of input variables.
 ## Please see fuzzifier function to contruct the matrix.
 varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
                       nrow = 5, byrow = FALSE)

 ## Define number of fuzzy terms of input variables.
 ## Suppose, we have 3, 2, 3, and 3 numbers of fuzzy terms 
 ## for first, second, third and fourth variables, respectively.
 num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
 
 ## Give the names of the fuzzy terms of each input variable.
 ## It should be noted that the names of the fuzzy terms must be unique,
 ## so we put a number for making it unique.
 varinput.1 <- c("a1", "a2", "a3")
 varinput.2 <- c("b1", "b2")
 varinput.3 <- c("c1", "c2", "c3")
 varinput.4 <- c("d1", "d2", "d3")
 names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

 ## Set interval of data.
 range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)

 ## Set weighted average method to be used as defuzzification method.
 type.defuz <- "WAM"
 ## We are using standard t-norm and s-norm.
 type.tnorm <- "MIN"
 type.snorm <- "MAX"
 type.implication.func <- "ZADEH"

 ## Give the name of simulation.
 name <- "Sim-0"

 ## Provide new data for testing. 
 newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow= 2, byrow = TRUE)
 ## the names of variables
 colnames.var <- c("input1", "input2", "input3", "input4", "output1")

 #####################################################################
 ## 1b. Using Takagi Sugeno Kang (TSK) Model 
 #####################################################################

 type.model <- "TSK"
 
 ## Define function of TSK 
 func.tsk<-matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), nrow=3, byrow=TRUE)
 
 ## Define the fuzzy IF-THEN rules; 
 ## For TSK model, it isn't necessary to put linguistic term in consequent parts.
 ## Make sure that each rule has a "->" sign. 
 rule <- matrix(c("very a1","and","b1","and","slghtly c1","and","d1","->",
                  "a2","and","extremely b2","and","c2","and","d2", "->",  
                  "a3","and","b2","and","c2","and","dont_care", "->"), 
                  nrow=3, byrow=TRUE) 
				  
 ## Generate a fuzzy model with frbs.gen.
 ## It should be noted that for TSK model, we do not need to input: 
 ## num.fvaloutput, varout.mf, names.varoutput, type.defuz.
 ## 
 object <- frbs.gen(range.data, num.fvalinput, names.varinput, num.fvaloutput = NULL, 
              varout.mf = NULL, names.varoutput = NULL, rule, 
				varinp.mf, type.model, type.defuz = NULL, type.tnorm, type.snorm, func.tsk, colnames.var, 
				type.implication.func, name)
				
## We can plot the membership function
 plotMF(object)

 ## Predicting using new data.
 res <- predict(object, newdata)$predicted.val