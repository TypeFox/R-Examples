 #################################################
 ## The following codes show how to generate a fuzzy model 
 ## using the frbs.gen function for classification tasks using Mamdani model. 
 #################################################
 ## define range of data.
 ## Note. we only define range of input data. 
 range.data.input <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1), nrow=2)
 
 ## Define shape and parameters of membership functions of input variables.
 ## Please see fuzzifier function to construct the matrix.
 ## In this case, we are using TRIANGLE for membership functions.
 varinp.mf <- matrix(c(1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA),
                       nrow = 5, byrow = FALSE)

 ## Define number of fuzzy terms of input variables.
 ## Suppose, we have 3, 3, 3, and 3 numbers of fuzzy terms 
 ## for first up to fourth variables, respectively.
 num.fvalinput <- matrix(c(3, 3, 3, 3), nrow=1)
 
 ## Give the names of the fuzzy terms of each input variable.
 ## It should be noted that the names of the fuzzy terms must be unique,
 ## so we put a number for making it unique.
 varinput.1 <- c("v.1_a.1", "v.1_a.2", "v.1_a.3")
 varinput.2 <- c("v.2_a.1", "v.2_a.2", "v.2_a.3")
 varinput.3 <- c("v.3_a.1", "v.3_a.2", "v.3_a.3")
 varinput.4 <- c("v.4_a.1", "v.4_a.2", "v.4_a.3")
 names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)

 ## Provide inference parameters.
 type.tnorm <- "MIN"
 type.snorm <- "MAX"
 type.implication.func <- "ZADEH"
 type.model <- "FRBCS"

 ## Give the name of simulation.
 name <- "Sim-0"

 ## Provide new data for testing. 
 newdata<- matrix(c(0.45, 0.5, 0.89, 0.44, 0.51, 0.99, 0.1, 0.98, 0.51,
                  0.56, 0.55, 0.5), nrow= 3, byrow = TRUE)

 ## the names of variables
 colnames.var <- c("input1", "input2", "input3", "input4", "output1")
 
 ## Construct rules.
 ## Take into account in consequent part, which expresses classes. 
 rule <- matrix(c("v.1_a.2","and","v.2_a.2","and","v.3_a.3","and","v.4_a.2","->","3",
          "v.1_a.2","and","v.2_a.3","and","v.3_a.1","and","v.4_a.3","->","1",
          "v.1_a.2","and","v.2_a.2","and","v.3_a.2","and","v.4_a.2","->","2"), 
          nrow=3, byrow=TRUE) 
 
 ## Generate frbs object.
 object <- frbs.gen(range.data = range.data.input, num.fvalinput, 
              names.varinput, num.fvaloutput = NULL, varout.mf = NULL, 
              names.varoutput = NULL, rule, varinp.mf, type.model, 
              type.defuz = NULL, type.tnorm, type.snorm, func.tsk = NULL, 
              colnames.var, type.implication.func, name)
				
 ## Plot the shape of membership functions.
 plotMF(object)

 ## Predicting using new data.
res <- predict(object, newdata)