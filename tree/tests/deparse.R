## PR#15696

ir.tr <- tree::tree(Species ~., iris)
head(model.frame(ir.tr))
## infinite-looped in 1.0-34

