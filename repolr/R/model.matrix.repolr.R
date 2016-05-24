model.matrix.repolr <-
function(object, ...){
 exdata <- ord.expand(space = object$poly.mod$space, formula = object$orig.formula, 
                times = object$times, poly = object$poly.mod$poly, data = eval(object$data), 
                subjects = object$subjects, categories = object$categories)
 mm <- model.matrix.default(exdata$formula, data = exdata$data)
 mm
}
