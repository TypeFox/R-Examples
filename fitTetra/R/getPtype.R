getPtype <-
function(model) {
  c("p.free","p.free","p.free","p.free","p.HW","p.HW","p.HW","p.HW")[((model-1)%%8)+1]
}
