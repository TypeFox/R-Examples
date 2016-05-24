#Simulate a data frame to test fonction tableplot
# by Joris Muller
# 2013

# to be interesting, this data frame must have NA in line, columns and randomly
# and also various kind of class of vector (int, num, factor, date…)

# subjet : int (no)
# initial (char)
# date of birth (date posix)
# sex (factor)
# heigh (num)
# weight (num)
# nsubjects of brother/sister (int)
# pseudo-questionnaire : score from 0 to 4 and more NA at the end 
# Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8, Q9, Q10



simulate_dataframe <- function(
  nsubjects = 100, # number of observations
  nquestions = 10, # number of questions
  includeMatrix = FALSE # Include a matrix inside the data.frame
  ){
  
  ################
  # Helpers functions
  generate_dumb_answer <- function(score_max = 4, nb_obs = 100){
    #generate an random answer for questionnaire
    return(rbinom(nb_obs,score_max,prob=rnorm(1,mean=0.5,sd=0.1) %% 1 ))
  }
  
  generate_initial <- function(nb_obs = 100){
    # generate random initials
    paste(
      letters[rbinom(n=nb_obs, size=25, prob=0.5)],
      letters[rbinom(n=nb_obs, size=25, prob=0.3)],
      sep=""
    )
  }
  
  ################
  # Create a complete data.frame
  
  # Create an empty data.frame and add one per one each column
  dumb_data <- data.frame("subject" = as.integer(1:nsubjects), stringsAsFactors = F)

  # Populate the data.frame with some complete data
  # First letter of the subject
  dumb_data$initial <- as.character(generate_initial())
  # Date of birth
  dumb_data$birth <- as.Date(
    # Use the RANDU dataset with random variables
    x=randu[1:nsubjects,1]*1000, 
    origin=as.Date(x="01/01/2000",format="%d/%m/%Y"))
  
  dumb_data$sex <- factor(
    x=rbinom(n=nsubjects,size=1,prob=0.45),
    labels=c("male","female")
    )
  
  dumb_data$heigh <- round(rnorm(nsubjects,mean=1.70,sd=0.10),digits=2)
  dumb_data$weight <- round(rnorm(nsubjects,mean=70,sd=9),digits=1)
  dumb_data$siblings <- rpois(nsubjects, lambda=3)
  dumb_data$study_level <- factor( x=rbinom(nsubjects,2,0.3), 
                                   levels=0:2, 
                                   labels=c("primary","secondary","superior"),
                                   ordered=T)
  
  # For the questions, should be a matrix (if includeMatrix = T)
  
  # Generate a matrix of 10*nbobs dumb data
  
  questions_matrix <- matrix(
    data=generate_dumb_answer(nb_obs=nquestions*nsubjects),
    ncol = nquestions
  )
  # add names
  colnames(questions_matrix) <- paste("Q",1:nquestions,sep="")



  ################
  # add som NAs
add_NAs <- function(x, proba=0.1){
  
  na_position<-as.logical(rbinom(n=length(x), size=1,prob=proba))
  x[na_position] <- NA
  return(x)
}


dumb_data$birth <- add_NAs(x=dumb_data$birth,0.02)
dumb_data$heigh <- add_NAs(x=dumb_data$heigh,0.04)
dumb_data$weight <- add_NAs(x=dumb_data$weight,0.05)
dumb_data$siblings <- add_NAs(x=dumb_data$siblings,0.30)

# For questions, random NA + one question (last) with a lot of NA

questions_matrix<-add_NAs(x=questions_matrix,0.1)
questions_matrix[,nquestions] <- add_NAs(x=questions_matrix[,nquestions],0.8)


# add NA on entire line (non repondant)
#dumb_data[ as.logical(rbinom(n=nsubjects,size=1,prob=0.07)),4:16] <- NA

######
# Merge demographic data and questions
# if includeMatrix = T then add the matrix as a variable
if(includeMatrix){
  test <- within( data=dumb_data,
                  expr= quest <- questions_matrix)
}else{
  test <- cbind(dumb_data,questions_matrix)
}

return(test)
}
