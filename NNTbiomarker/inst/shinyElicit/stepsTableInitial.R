stepsTableInitial = as.data.frame(stringsAsFactors = F, matrix(ncol=4, byrow=T, c(
  "1","Not yet done",
  "Introduction: the clinical scenario",
  "What is the clinical challenge the biomarker should address?"
  ,"2","Not yet done",
  "Patients and options",
  "Who are the patients to help? What are the clinical decision options?"
  ,"3","Not yet done",
  "Criterion for a useful biomarker: Defining the NNT discomfort zone",
  HTML("What region of NNTs (number needed to treat) make both choices uncomfortable?"
  %&% " <br>What NNT's for positive and negative test groups would make thse decision clear-cut? "
  %&% " <br>What predictive values for positive and negative tests will give NNT's making the decision clear-cut?")
  ,"4","Not yet done",
  "Specific clinical benefit",
  "How, in detail, will patients be helped by a test that achieves the criteria?"
  ,"5","Not yet done",
  "Prospective study requirements",
  "Given the target predictive values how large should a prospective study be, and how long the follow-up?"
  ,"6","Not yet done",
  "Retrospective study requirements",
  "Given a prevalence value, what sensitivity and specificity do we hope for, and what should the sample sizes be to estimate them sufficiently?"
)))

stepsTableInitial = stepsTableInitial[ , -1]
names(stepsTableInitial) = c("Completed?","SteppingStone","Question")

