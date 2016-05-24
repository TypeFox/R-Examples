stepsTableInitial = as.data.frame(stringsAsFactors = F, matrix(ncol=4, byrow=T, c(
  "1","Not yet done",
"Defining the clinical scenario",
"Who are the patients to help? What are the clinical decision options?"
  ,"2","Not yet done",
"Defining the NNT discomfort zone",
"What region of NNTs (number needed to treat) make both decisions 'treat' and 'wait' uncomfortable?"
  ,"3","Not yet done",
"Choosing target NNTpos and NNTneg outside this zone",
"What NNT's for the BestToAct and BestToWait groups would make the decision clear-cut?"
  ,"4","Not yet done",
"Defining how patients would benefit",
"Specifically how will patients be helped by a test that achieves these NNT s?"
  ,"5","Not yet done",
"Classification performance needed",
"What predictive values do these NNT s correspond to?"
  ,"6","Not yet done",
"Prospective study requirements",
"Given these NNT s, how large should a prospective study be, and how long the follow-up?"
  ,"7","Not yet done",
"Retrospective study requirements",
"Given a prevalence, what sensitivity and specificity do we hope for, and what should the sample sizes be to estimate them sufficiently?"
)) [ , -1])

names(stepsTableInitial) = c("Done?","Stepping stone","Question")

