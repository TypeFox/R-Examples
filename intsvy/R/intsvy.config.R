intsvy.config <- function(  variables.pvlabelpref,
           variables.pvlabelsuff,
           variables.weight,
           variables.jackknifeZone,
           variables.jackknifeRep,
           parameters.cutoffs,
           parameters.cutoffs2,
           parameters.percentiles,
           parameters.weights,
           parameters.PVreps,
           input.type,
           input.prefixes,
           input.student,
           input.student_colnames1,
           input.student_colnames2,
           input.student_pattern,
           input.homeinput,
           input.home_colnames,
           input.school,
           input.school_colnames,
           input.teacher,
           input.teacher_colnames,
           input.student_ids,
           input.school_ids,
           input.type_part,
           input.cnt_part, base.config = pirls_conf) {
  config <- base.config

  if (!missing(variables.pvlabelpref	)) config$variables$pvlabelpref	 =variables.pvlabelpref
  if (!missing(variables.pvlabelsuff	)) config$variables$pvlabelsuff	 =variables.pvlabelsuff
  if (!missing(variables.weight	)) config$variables$weight	 =variables.weight
  if (!missing(variables.jackknifeZone	)) config$variables$jackknifeZone	 =variables.jackknifeZone
  if (!missing(variables.jackknifeRep	)) config$variables$jackknifeRep	 =variables.jackknifeRep
  
  if (!missing(parameters.cutoffs	)) config$parameters$cutoffs	 =parameters.cutoffs
  if (!missing(parameters.cutoffs2	)) config$parameters$cutoffs2	 =parameters.cutoffs2
  if (!missing(parameters.percentiles	)) config$parameters$percentiles	 =parameters.percentiles
  if (!missing(parameters.weights	)) config$parameters$weights	 =parameters.weights
  if (!missing(parameters.PVreps	)) config$parameters$PVreps	 =parameters.PVreps
  
  if (!missing(input.type	)) config$input$type	 =input.type
  if (!missing(input.prefixes	)) config$input$prefixes	 =input.prefixes
  if (!missing(input.student	)) config$input$student	 =input.student
  if (!missing(input.student_colnames1	)) config$input$student_colnames1	 =input.student_colnames1
  if (!missing(input.student_colnames2	)) config$input$student_colnames2	 =input.student_colnames2
  if (!missing(input.student_pattern	)) config$input$student_pattern	 =input.student_pattern
  if (!missing(input.homeinput	)) config$input$homeinput	 =input.homeinput
  if (!missing(input.home_colnames	)) config$input$home_colnames	 =input.home_colnames
  if (!missing(input.school	)) config$input$school	 =input.school
  if (!missing(input.school_colnames	)) config$input$school_colnames	 =input.school_colnames
  if (!missing(input.teacher	)) config$input$teacher	 =input.teacher
  if (!missing(input.teacher_colnames	)) config$input$teacher_colnames	 =input.teacher_colnames
  if (!missing(input.student_ids	)) config$input$student_ids	 =input.student_ids
  if (!missing(input.school_ids	)) config$input$school_ids	 =input.school_ids
  if (!missing(input.type_part	)) config$input$type_part	 =input.type_part
  if (!missing(input.cnt_part	)) config$input$cnt_part	 =input.cnt_part
  
  config
}

