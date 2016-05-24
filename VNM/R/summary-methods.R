setMethod("summary",signature("OPT"),
    function(object,...) {
        D<-object@Opt
        eff<-object@Eff
        obj<-object@Par

        if(obj@fid=="ceff1") {
             cat("c-optimal design for ED50:","\n")
             print(D)
             cat("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _","\n")
             cat("c-efficiency for ED50:",eff,"\n")
        }
        if(obj@fid=="ceff2") {
             cat("c-optimal design for MED:","\n")
             print(D)
             cat("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _","\n")
             cat("c-efficiency for MED:",eff,"\n")
        }
        if(obj@fid=="Deff") {
             cat("D-optimal design:","\n")
             print(D)
             cat("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _","\n")
             cat("D-efficiency:",eff,"\n")
        }
        if(obj@fid=="MOPT") {
             cat("Obtained multiple-objective optimal design:","\n")
             cat("X-selected dose levels; W-corresponding proportions of subjects","\n")
             print(D)
        }
    }
)


setMethod("summary",signature("SW"),
    function(object,...) {
        D<-object@Opt.W
        cat("Optimal Proportions of Subjects:","\n")
        print(D)
        cat("_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _","\n")
        f1<-object@First.C
        cat("The first derivative of the criterion:",f1,"\n")
        f2<-object@Second.C
        cat("The second derivative of the criterion:",f2,"\n")
    }
)

        
