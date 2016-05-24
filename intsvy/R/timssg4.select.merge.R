timssg4.select.merge <-
function(folder=getwd(), countries, student=c(), home, school, teacher, use.labels=TRUE) {
  intsvy.select.merge(folder=folder, countries=countries, student=student, home=home, school=school, teacher = teacher, use.labels=use.labels,
                      config=timss4_conf)
  
}
