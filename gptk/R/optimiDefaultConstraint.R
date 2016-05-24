optimiDefaultConstraint <-
function (constraint) {
  if ( constraint == "positive" ) {
    return (list(func="expTransform", hasArgs=FALSE))
  } else if ( constraint == "zeroone" ) {
    return (list(func="sigmoidTransform", hasArgs=FALSE))
  } else if ( constraint == "bounded" ) {
    return (list(func="boundedTransform", hasArgs=TRUE))
  }
}
