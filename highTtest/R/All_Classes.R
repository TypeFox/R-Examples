setClassUnion("matrix or NULL", members = c("matrix","NULL"))
setClassUnion("numeric or NULL", members = c("numeric","NULL"))

setClass("highTtest", 
  representation(CK = "matrix or NULL", 
                 pi1 = "numeric or NULL",
                 pvalue = "numeric",
                 ST = "matrix or NULL", 
                 BH = "matrix or NULL",
                 gammas = "numeric"),
  prototype(CK = matrix(NA,1,1), 
            pi1 = NA,
            pvalue = numeric(1),
            ST = matrix(NA,1,1), 
            BH = matrix(NA,1,1),
            gammas = numeric(1))
)

