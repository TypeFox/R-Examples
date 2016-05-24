#
#
#Plots a tableplot of probability of skill mastery for each student
#inputs: skill_probs: an n_student by (number auxilery vars + n_attributes + 1) matrix or dataframe of:
#           auxillery variables use in mplus model, probability of mastery of each attribute, class assignment
#         n_attribtues: a numeric value of number of attributes
PlotSkillMasteryTableplot <- function(dat, ngroups, is.max.class, divide.by = 18){
  indicator.matrix <- matrix(2, nrow(dat), ncol(dat))
  indicator.matrix[ , 1] <- 1
  if (is.max.class){
    indicator.matrix[ , ncol(indicator.matrix)] <- 3
    vparts <- c(1, ngroups, 1)
    cell.specs = list(
      list(5, "black", 0, 0, "red", 1, "white", "white", 4, 1, "black", FALSE, "grey40", 100),
      list(0, "blue", 1, 0, "red", 1, "grey95", "white", 0, 0, "grey50", FALSE, "grey40", 1),
      list(5, "black", 0, 0, "red", 1, "white", "white", 2, 2, "black", FALSE, "grey40", 2))
  } else {
    vparts <- c(1, ngroups)
    cell.specs = list(
      list(5, "black", 0, 0, "red", 1, "white", "white", 4, 1, "black", FALSE, "grey40", 100),
      list(0, "blue", 1, 0, "red", 1, "grey95", "white", 0, 0, "grey50", FALSE, "grey40", 1))
  }
  nparts <- ceiling(nrow(dat)/divide.by)
  dat <- as.matrix(dat)
  for (i in 1:nparts){
    if (i == nparts){
      tableplot(dat[((i-1)*divide.by + 1):nrow(dat),], cell.specs = cell.specs,
                gap = 10, table.label = TRUE, v.parts = vparts,
                , assign.sets = indicator.matrix[((i-1)*divide.by + 1):nrow(dat),])
    }
    else{
      tableplot(dat[((i-1)*divide.by + 1):(i*divide.by),], cell.specs = cell.specs,
                gap = 10, table.label = TRUE, v.parts = vparts,
                , assign.sets = indicator.matrix[((i-1)*divide.by + 1):(i*divide.by),])
    }
  }
}