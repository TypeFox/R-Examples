############################################################################################################################

# Function: CreateTableDesign.
# Argument: data.strucure and label (optional).
# Description: Generate a summary table of design parameters for the report.

CreateTableDesign = function(data.structure, label = NULL) {

  # Number of design
  n.design = length(data.structure$design.parameter.set)

  # Label
  if (is.null(label)) label = paste0("Design ", 1:n.design)
  else label = unlist(label)
  if (length(label) != n.design)
    stop("Summary: Number of the design parameters labels must be equal to the number of design parameters sets.")

  # Summary table
  design.table <- matrix(nrow = n.design, ncol = 9)

  design.parameter.set = data.structure$design.parameter.set

  for (i in 1:n.design) {
    design.table[i, 1] = i
    design.table[i, 2] = label[i]
    design.table[i, 3] = design.parameter.set[[i]]$enroll.period
    if (!is.na(design.parameter.set[[i]]$enroll.dist)){
      if (design.parameter.set[[i]]$enroll.dist=="UniformDist") enroll.dist.par.dummy = list(max = design.parameter.set[[i]]$enroll.period)
      else enroll.dist.par.dummy = design.parameter.set[[i]]$enroll.dist.par
      enroll.dist.desc = do.call(design.parameter.set[[i]]$enroll.dist,list(list("description",enroll.dist.par.dummy)))
      design.table[i, 4] = unlist(enroll.dist.desc[[2]])
      if (!is.na(design.parameter.set[[i]]$enroll.dist.par)){
        enroll.dist.npar = length(enroll.dist.desc[[1]][[1]])
        enroll.dist.par = ""
        for (j in 1: enroll.dist.npar){
          enroll.dist.par = paste0(", ",enroll.dist.desc[[1]][[1]][j],"=",design.parameter.set[[i]]$enroll.dist.par[j])
        }
        enroll.dist.par = sub(", ","",enroll.dist.par)
        design.table[i, 5] = enroll.dist.par
      }
      else design.table[i, 5] = design.parameter.set[[i]]$enroll.dist.par
    }
    else {
      design.table[i, 4] = design.parameter.set[[i]]$enroll.dist
      design.table[i, 5] = design.parameter.set[[i]]$enroll.dist.par
    }
    design.table[i, 6] = design.parameter.set[[i]]$followup.period
    design.table[i, 7] = design.parameter.set[[i]]$study.duration
    if (!is.na(design.parameter.set[[i]]$dropout.dist)){
      dropout.dist.desc = do.call(design.parameter.set[[i]]$dropout.dist,list(list("description",design.parameter.set[[i]]$dropout.dist.par)))
      design.table[i, 8] = unlist(dropout.dist.desc[[2]])
      if (!is.na(design.parameter.set[[i]]$dropout.dist.par)){
        dropout.dist.npar = length(dropout.dist.desc[[1]][[1]])
        dropout.dist.par = ""
        for (j in 1: dropout.dist.npar){
          dropout.dist.par = paste0(", ",dropout.dist.desc[[1]][[1]][j],"=",design.parameter.set[[i]]$dropout.dist.par[j])
        }
        dropout.dist.par = sub(", ","",dropout.dist.par)
        design.table[i, 9] = dropout.dist.par
      }
      else design.table[i, 9] = design.parameter.set[[i]]$enroll.dist.par
    }
    else {
      design.table[i, 8] = design.parameter.set[[i]]$dropout.dist
      design.table[i, 9] = design.parameter.set[[i]]$dropout.dist.par
    }
  }
  design.table = as.data.frame(design.table)
  colnames(design.table) = c("design.parameter", "Design parameter set", "Enrollment period", "Enrollment distribution", "Enrollment distribution parameter", "Follow-up period", "Study duration", "Dropout distribution", "Dropout distribution parameter")
  return(design.table)
}
# End of CreateTableDesign
