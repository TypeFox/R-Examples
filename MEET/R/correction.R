setClass("correction",
        representation(Herror="vector",sderror="vector"),
        prototype(Herror=vector("numeric"),sderror=vector("numeric")))