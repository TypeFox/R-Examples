# generating function for class 'FunSymmList'
FunSymmList <- function(...){ 
    new("FunSymmList", list(...)) 
}

#setAs(from = "FunctionSymmetry", to = "FunSymmList", 
#    def = function(from){ new("FunSymmList", list(from)) })
