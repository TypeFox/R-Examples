if(!isGeneric("df.residual"))
  setGeneric(name = "df.residual", def = function(object, ...) standardGeneric("df.residual"))

### method df.residual for objects of class "glimML" (fitted with betabin or negbin)
setMethod(f = "df.residual", signature = "glimML", definition = function(object, ...) object@df.residual)

### method df.residual for objects of class "gliQML" (fitted with quasibin or quasipois)
setMethod(f = "df.residual", signature = "glimQL", definition = function(object, ...) df.residual(object@fm))
