`make.reg.fnc` <-
function(
  nsubj = 10,           # number of subjects
  nitem = 20,           # number of items
  beta=c(400, 2, 6, 4), # fixed effect coefficients
  learn = FALSE,        # if TRUE, add learning effect
  learnRate = 10,       # default learning rate
  stdevItem = 40,       # st.dev. item random effect
  stdevSubj = 80,       # st.dev. subject random effect
  stdevError = 50)      # st.dev. residual error
{

  # the by-item fixed effect data skeleton 
  # (identical for each subject)
  fixed = data.frame(Intercept = rep(1, nitem), 
    X = 1:nitem, Y = sample(1:nitem), Z = sample(1:nitem),
    Item = as.factor(paste(rep("Item", nitem), 1:nitem, sep="")),
    RanefItem = rnorm(nitem, 0, stdevItem))

  Subject = as.factor(rep(paste("Subj", 1, sep=""), nitem))
  RanefSubj = rnorm(nsubj, 0, stdevSubj)

  # create the data for the first subject 
  Data = fixed
  Data$RanefSubj = rep(RanefSubj[1], nitem)
  Data$Subject = Subject
  # and add further subjects
  for (s in 2:nsubj) {
    Dat = fixed
    Dat$RanefSubj = rep(RanefSubj[s], nitem)
    Dat$Subject = as.factor(rep(paste("Subj", s, sep=""), nitem))
    Data = rbind(Data, Dat)
  }
  # add learning effect if required
  # full counterbalancing: each subject has own list order
  if (learn == TRUE){
    # trial numbers for first subject
    trial = sample(1:nitem)  
    for (s in 2:nsubj){
      # trial numbers for other subjects
      trial = c(trial, sample(1:nitem)) 
    }
  }
  dimnames(Data)[[1]] = 1:nrow(Data)
  # the error vector
  Data$Error = rnorm(nitem*nsubj, 0, stdevError)
  # calculate RTs using matrix multiplication
  Data$RT = as.matrix(Data[,c(1:4,6,7,9)]) %*% c(beta, 1, 1, 1)
  # add effect of learning
  if (learn == TRUE) {
    Data$Trial = trial
    Data$RT = Data$RT + learnRate*Data$Trial
  } 
  return(Data)
}

