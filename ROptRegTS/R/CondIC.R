## generating function
CondIC <- function(name, Curve = EuclRandVarList(EuclRandVariable(Map = list(function(x){x[1]*x[2]}),
                                            Domain = EuclideanSpace(dimension = 2),
                                            Range = Reals())),
                    Risks, Infos, CallL2Fam = call("L2RegTypeFamily")){
    if(missing(name))
        name <- "Influence curve for a L_2 differentiable regression type family"
    if(missing(Risks))
        Risks <- list()
    if(missing(Infos))
        Infos <- matrix(c(character(0),character(0)), ncol=2,
                     dimnames=list(character(0), c("method", "message")))
    return(new("CondIC", name = name, Curve = Curve, Risks = Risks,
               Infos = Infos, CallL2Fam = CallL2Fam))
}

## replace methods
setReplaceMethod("CallL2Fam", "CondIC",
    function(object, value){
        object@CallL2Fam <- value
        validObject(object)
        object
    })

setMethod("checkIC", signature(IC = "CondIC", L2Fam = "missing"), 
    function(IC, out = TRUE){ 
        L2Fam <- eval(IC@CallL2Fam)
        K <- L2Fam@RegDistr
        if(is(K, "DiscreteDistribution") || is(K, "DiscreteMVDistribution"))
            cond <- as.matrix(support(K))
        else{
            if(is(K, "AbscontDistribution"))
                cond <- as.matrix(seq(from = q(K)(TruncQuantile), to = q(K)(1-TruncQuantile),
                            length = 100))
            else
                cond <- as.matrix(r(K)(1000))
        }

        trafo <- L2Fam@param@trafo
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        cent <- array(0, c(length(IC1), length(cond), nrow(trafo)))
        for(i in 1:length(IC1)){
            fct <- function(x, cond, f1){ f1(cbind(t(cond),x)) }
            cent[i,,] <- apply(cond, 1, .condE, D1 = L2Fam@distribution, fct = fct, 
                            f1 = IC1@Map[[i]])
        }
        if(out)
            cat("precision of conditional centering:\t", max(abs(cent)), "\n")

        dims <- length(L2Fam@param)
        if(is(L2Fam@distribution, "UnivariateCondDistribution")){
            L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")
            IC.L2 <- IC1 %*% t(L2deriv)
            res <- numeric(length(IC.L2))
            for(i in 1:length(IC.L2)){
                fct <- function(x, cond, f1){ f1(cbind(t(cond),x)) }
                res[i] <- E(K, .condE, D1 = L2Fam@distribution, fct = fct, 
                               f1 = IC.L2@Map[[i]])
            }            
            consist <- matrix(res, nrow = nrow(trafo)) - trafo
            if(out){
                cat("precision of Fisher consistency:\n")
                print(consist)
            }
        }else{
            stop("not yet implemented")
        }
        res <- max(abs(cent), abs(consist))
        names(res) <- "maximum deviation"
        
        return(res)
    })

setMethod("checkIC", signature(IC = "CondIC", L2Fam = "L2RegTypeFamily"), 
    function(IC, L2Fam, out = TRUE){ 
        K <- L2Fam@RegDistr
        if(is(K, "DiscreteDistribution") || is(K, "DiscreteMVDistribution"))
            cond <- as.matrix(support(K))
        else{
            if(is(K, "AbscontDistribution"))
                cond <- as.matrix(seq(from = q(K)(TruncQuantile), to = q(K)(1-TruncQuantile),
                            length = 100))
            else
                cond <- as.matrix(r(K)(1000))
        }

        trafo <- L2Fam@param@trafo
        IC1 <- as(diag(nrow(trafo)) %*% IC@Curve, "EuclRandVariable")
        cent <- array(0, c(length(IC1), length(cond), nrow(trafo)))
        for(i in 1:length(IC1)){
            fct <- function(x, cond, f1){ f1(cbind(t(cond),x)) }
            cent[i,,] <- apply(cond, 1, .condE, D1 = L2Fam@distribution, fct = fct, 
                            f1 = IC1@Map[[i]])
        }
        if(out)
            cat("precision of conditional centering:\t", max(abs(cent)), "\n")

        dims <- length(L2Fam@param)
        if(is(L2Fam@distribution, "UnivariateCondDistribution")){
            L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")
            IC.L2 <- IC1 %*% t(L2deriv)
            res <- numeric(length(IC.L2))
            for(i in 1:length(IC.L2)){
                fct <- function(x, cond, f1){ f1(cbind(t(cond),x)) }
                res[i] <- E(K, .condE, D1 = L2Fam@distribution, fct = fct, 
                               f1 = IC.L2@Map[[i]])                
            }            
            consist <- matrix(res, nrow = nrow(trafo)) - trafo
            if(out){
                cat("precision of Fisher consistency:\n")
                print(consist)
            }
        }else{
            stop("not yet implemented")
        }
        res <- max(abs(cent), abs(consist))
        names(res) <- "maximum deviation"
        
        return(res)
    })
