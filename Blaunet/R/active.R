active <-
function(blauObj){
  sprintf('These elements are active: %s. Access with object$element.', paste(names(blauObj), collapse = ", "))
}
