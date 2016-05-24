#
# vim:set ff=unix expandtab ts=2 sw=2:

setClass(# decomposition operator 
    Class="DecompOp",
    ,
    contains="VIRTUAL"
)
setMethod(
  f="DecompOp",
  signature=signature(object="DecompOp"),
  def=function # pass through constructor
  ### This method handles the case that no actual construction is necessary since
  ### the argument is already of a subclass of DecompOp 
  ##<<details This is useful to simplify argument handling of functions which rely on 
  ## the presence of a DecompOp. 
  ## Due to this method those functions can always
  ## call DecompOp(something) without having to check if 
  ## it is necessary.
  (object){
    object
    ### the unchanged argument
  }
)
setMethod(
  "DecompOp",
  signature(object="matrix"),
  #valueClass="ConstLinDecompOp",
  def=function # creates a ConstanDecompOp from a matrix
  ### The resulting operator is creted by a call to the constructor of class
  ### ConstLinDecompOp
  (object){
    ConstLinDecompOp(object)
  }
)
setMethod(
  f="DecompOp",
  signature=signature(object="TimeMap"),
  #valueClass="BoundLinDecompOp",
  def=function # creates a BoundLinDecompOp from a TimeMap object
  ### The resulting operator is creted by a call to the constructor of class
  ### BoundLinDecompOp
  ### The method is used to ensure backward compatibility with the now deprecated
  ### TimeMap class
  (object){
    BoundLinDecompOp(object)
  }
)
