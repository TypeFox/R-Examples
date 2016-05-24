timssg8.select.merge <-
function(folder=getwd(), countries, student=c(), school, math.teacher, science.teacher, use.labels=TRUE) {
  
  if (!missing(math.teacher) && !missing(science.teacher)) 
    stop("not possible to merge science and math teacher data")
  
  if (!missing(math.teacher)) {
    intsvy.select.merge(folder=folder, countries=countries, student=student, school=school, teacher = math.teacher, use.labels=use.labels,
                        config=timss8_conf)
  } else {
    intsvy.select.merge(folder=folder, countries=countries, student=student, school=school, teacher = science.teacher, use.labels=use.labels,
                        config=timss8_conf)
  }

}
