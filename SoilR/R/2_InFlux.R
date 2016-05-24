#
# vim:set ff=unix expandtab ts=2 sw=2:
setClass(
   Class="InFlux",
   contains="VIRTUAL"
)
setMethod(
  f="InFlux",
  signature(object="TimeMap"),
  def=function #create a BoundInFlux from a TimeMap object
  ### The method is used to ensure backward compatibility with the now deprecated
  ### TimeMap class.
  ### The resulting BoundInFlux is created by a call to
  ### the constructor BoundInFlux(object) of that class.
  (object)
  {
    BoundInFlux(object)
  }

)
setMethod(
  f="InFlux",
  signature=signature(object="InFlux"),
  def=function # pass through constructor
  ### This method handles the case that no actual construction is necessary since
  ### the argument is already of a subclass of InFlux 
  ##<<details This is useful to simplify argument handling of functions which rely on 
  ## the presence of some kind of an InFlux 
  ## Due to this method those functions can 
  ## call InFlux(something) without having to check if 
  ## it is necessary.
  (object){
    object
    ### the unchanged argument
  }
)
