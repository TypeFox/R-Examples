procformula <- function(object=NULL,value,exo=lava.options()$exogenous,...) {

    ## Split into reponse and covariates by ~ disregarding expressions in parantheses
    ## '(?!...)' Negative lookahead assertion
    regex <- "~(?![^\\(].*\\))"
    yx <- lapply(strsplit(as.character(value),regex,perl=TRUE),function(x) gsub(" ","",x))[-1]
    iscovar <- FALSE
    if (length(yx)==1) {
        lhs <- NULL; xidx <- 1
    } else {
            lhs <- yx[1]; xidx <- 2
            if (yx[[xidx]][1]=="") {
                yx[[xidx]] <- yx[[xidx]][-1]
                iscovar <- TRUE
            }
    }
    ##Check for link function
    invlink <- NULL
    if (xidx==2) {
        if (length(grep("[a-zA-Z0-9_]*\\(.*\\)$",yx[[xidx]]))>0) { ## rhs of the form F(x+y)
            invlink <- strsplit(yx[[xidx]],"\\(.*\\)")[[1]][1]
                if (invlink%in%c("f","v","I","") ||
                    grepl("+",invlink))
                { ## Reserved for setting linear constraints
                    invlink <- NULL
                } else {
                    yx[[xidx]] <- gsub(paste0(invlink,"\\(|\\)$"),"",yx[[xidx]])
                }
            }
    }

    ## Handling constraints with negative coefficients
    ## while not tampering with formulas like y~f(x,-2)
    st <- yx[[xidx]]
    st <- gsub("\\-","\\+\\-",gsub("\\+\\-","\\-",st)) ## Convert - to +- (to allow for splitting on '+')
    ##gsub("[^,]\\-","\\+\\-",st) ## Convert back any - not starting with ','
    st <- gsub(",\\+",",",st) ## Remove + inside 'f' and 'v' constraints
    st <- gsub("^\\+","",st) ## Remove leading plus
    yx[[xidx]] <- st

    ## Match '+' but not when preceeded by ( ... )
    X <- strsplit(yx[[xidx]],"\\+(?![^\\(]*\\))", perl=TRUE)[[1]]
    
    ##regex <- "(?!(\\(*))[\\(\\)]"
    regex <- "[\\(\\)]"
    ## Keep squares brackets and |(...) statements
    ## Extract variables from expressions like
    ## f(x,b) -> x,b  and  2*x -> 2,cx
    ## but avoid to tamper with transformation expressions:
    ## a~(x*b)
    res <- lapply(X,decomp.specials,regex,pattern2="\\*",pattern.ignore="~",reverse=TRUE,perl=TRUE)
    ##OLD:
    ##res <- lapply(X,decomp.specials,pattern2="[*]",reverse=TRUE)
    xx <- unlist(lapply(res, function(x) x[1]))

    xxf <- lapply(as.list(xx),function(x) decomp.specials(x,NULL,pattern2="\\[|~",perl=TRUE))
    xs <- unlist(lapply(xxf,function(x) x[1]))

    ## Alter intercepts?    
    intpos <- which(vapply(xs,function(x) grepl("^[\\+\\-]*[\\.|0-9]+$",x), 0)==1)
    ## Match '(int)'
    intpos0 <- which(vapply(X,function(x) grepl("^\\([\\+\\-]*[\\.|0-9]+\\)$",x),0)==1)
    
    yy <- ys <- NULL
    if (length(lhs)>0) {
        yy <- decomp.specials(lhs)
        yyf <- lapply(yy,function(y) decomp.specials(y,NULL,pattern2="[",fixed=TRUE))
        ys <- unlist(lapply(yyf,function(x) x[1]))
    }
    
    if (!is.null(object)) {
      notexo <- c()
      if (length(lhs)>0) {
        object <- addvar(object,ys,reindex=FALSE,...)
        notexo <- ys
        ## Add link-transformation
        if (!is.null(invlink)) {
            if (invlink=="") {
                object <- transform(object,ys,NULL,post=FALSE)
                covariance(object,ys) <- NA
            } else {
                ff <- function(x) {};  body(ff) <- parse(text=paste0(invlink,"(x)"))
                object <- transform(object,ys,ff,post=FALSE)
                covariance(object,ys) <- 0
            }
        }
      }

      if (length(intpos>0)) {
          xs[intpos[1]] <- gsub("\\+","",xs[intpos[1]])
          if (xs[intpos[1]]==1 && (!length(intpos0)>0) ) {
              xs[intpos[1]] <- NA
          }
          intercept(object,ys) <- as.numeric(xs[intpos[1]])
          xs <- xs[-intpos]
          res[intpos] <- NULL
      }

        object <- addvar(object,xs,reindex=FALSE ,...)
        exolist <- c()

        for (i in seq_len(length(xs))) {

            ## Extract transformation statements: var~(expr)
            xf0 <- strsplit(xx[[i]],"~")[[1]]
            if (length(xf0)>1) {
                myexpr <- xf0[2]
                ftr <- toformula(y="",x=paste0("-1+I(",myexpr,")"))
                xtr <- all.vars(ftr)
                xf0 <- xf0[1]
                transform(object, y=xf0, x=xtr) <- function(x) {
                    structure(model.matrix(ftr,as.data.frame(x)),dimnames=list(NULL,xf0))
                }
            }

            xf <- unlist(strsplit(xf0,"[\\[\\]]",perl=TRUE))
            if (length(xf)>1) {
                xpar <- strsplit(xf[2],":")[[1]]
                if (length(xpar)>1) {
                    val <- ifelse(xpar[2]=="NA",NA,xpar[2])
                    valn <- suppressWarnings(as.numeric(val))
                    covariance(object,xs[i]) <- ifelse(is.na(valn),val,valn)
                }
                val <- ifelse(xpar[1]=="NA",NA,xpar[1])
                valn <- suppressWarnings(as.numeric(val))
                if (is.na(val) || val!=".") {
                    intercept(object,xs[i]) <- ifelse(is.na(valn),val,valn)
                    notexo <- c(notexo,xs[i])
                }
            } else { exolist <- c(exolist,xs[i]) }
        }

        for (i in seq_len(length(ys))) {
            y <- ys[i]
            yf <- unlist(strsplit(yy[i],"[\\[\\]]",perl=TRUE))
            if (length(yf)>1) {
                ypar <- strsplit(yf[2],":")[[1]]
                if (length(ypar)>1) {
                    val <- ifelse(ypar[2]=="NA",NA,ypar[2])
                    valn <- suppressWarnings(as.numeric(val))
                    covariance(object,y) <- ifelse(is.na(valn),val,valn)
                }
                val <- ifelse(ypar[1]=="NA",NA,ypar[1])
                valn <- suppressWarnings(as.numeric(val))
                if (is.na(val) || val!=".") {
                    intercept(object,y) <- ifelse(is.na(valn),val,valn)
                }
            }
        }

        curvar <- index(object)$var
        if (exo) {
            oldexo <- exogenous(object)
            newexo <- setdiff(exolist,c(notexo,curvar,ys))
            exogenous(object) <- union(newexo,setdiff(oldexo,notexo))
        }
    }

    return(list(object=object,
                yx=yx,
                X=X,
                ys=ys,
                xx=xx,
                xs=xs,
                yy=yy,
                ys=ys,
                res=res,
                notexo=notexo,
                intpos=intpos,
                invlink=invlink,
                lhs=lhs,
                iscovar=iscovar))
}


##' Add regression association to latent variable model
##'
##' Define regression association between variables in a \code{lvm}-object and
##' define linear constraints between model equations.
##'
##'
##' The \code{regression} function is used to specify linear associations
##' between variables of a latent variable model, and offers formula syntax
##' resembling the model specification of e.g. \code{lm}.
##'
##' For instance, to add the following linear regression model, to the
##' \code{lvm}-object, \code{m}:
##' \deqn{ E(Y|X_1,X_2) = \beta_1 X_1 + \beta_2 X_2}
##' We can write
##'
##' \code{regression(m) <- y ~ x1 + x2}
##'
##' Multivariate models can be specified by successive calls with
##' \code{regression}, but multivariate formulas are also supported, e.g.
##'
##' \code{regression(m) <- c(y1,y2) ~ x1 + x2}
##'
##' defines
##' \deqn{ E(Y_i|X_1,X_2) = \beta_{1i} X_1 + \beta_{2i} X_2 }
##'
##' The special function, \code{f}, can be used in the model specification to
##' specify linear constraints. E.g. to fix \eqn{\beta_1=\beta_2}
##' , we could write
##'
##' \code{regression(m) <- y ~ f(x1,beta) + f(x2,beta)}
##'
##' The second argument of \code{f} can also be a number (e.g. defining an
##' offset) or be set to \code{NA} in order to clear any previously defined
##' linear constraints.
##'
##' Alternatively, a more straight forward notation can be used:
##'
##' \code{regression(m) <- y ~ beta*x1 + beta*x2}
##'
##' All the parameter values of the linear constraints can be given as the right
##' handside expression of the assigment function \code{regression<-} (or
##' \code{regfix<-}) if the first (and possibly second) argument is defined as
##' well. E.g:
##'
##' \code{regression(m,y1~x1+x2) <- list("a1","b1")}
##'
##' defines \eqn{E(Y_1|X_1,X_2) = a1 X_1 + b1 X_2}. The rhs argument can be a
##' mixture of character and numeric values (and NA's to remove constraints).
##'
##' The function \code{regression} (called without additional arguments) can be
##' used to inspect the linear constraints of a \code{lvm}-object.
##'
##' For backward compatibility the "$"-symbol can be used to fix parameters at
##' a given value. E.g. to add a linear relationship between \code{y} and
##' \code{x} with slope 2 to the model \code{m}, we can write
##' \code{regression(m,"y") <- "x$2"}.  Similarily we can use the "@@"-symbol to
##' name parameters. E.g. in a multiple regression we can force the parameters
##' to be equal: \code{regression(m,"y") <- c("x1@@b","x2@@b")}.  Fixed parameters
##' can be reset by fixing (with \$) them to \code{NA}.
##'
##' @aliases regression regression<- regression<-.lvm regression.lvm regfix
##' regfix regfix<- regfix.lvm regfix<-.lvm
##' @param object \code{lvm}-object.
##' @param value A formula specifying the linear constraints or if
##' \code{to=NULL} a \code{list} of parameter values.
##' @param to Character vector of outcome(s) or formula object.
##' @param from Character vector of predictor(s).
##' @param fn Real function defining the functional form of predictors (for
##' simulation only).
##' @param silent Logical variable which indicates whether messages are turned
##' on/off.
##' @param additive If FALSE and predictor is categorical a non-additive effect is assumed
##' @param y Alias for 'to'
##' @param x Alias for 'from'
##' @param quick Faster implementation without parameter constraints
##' @param \dots Additional arguments to be passed to the low level functions
##' @usage
##' \method{regression}{lvm}(object = lvm(), to, from, fn = NA,
##' silent = lava.options()$silent, additive=TRUE, y, x, value, ...)
##' \method{regression}{lvm}(object, to=NULL, quick=FALSE, ...) <- value
##' @return A \code{lvm}-object
##' @note Variables will be added to the model if not already present.
##' @author Klaus K. Holst
##' @seealso \code{\link{intercept<-}}, \code{\link{covariance<-}},
##' \code{\link{constrain<-}}, \code{\link{parameter<-}},
##' \code{\link{latent<-}}, \code{\link{cancel<-}}, \code{\link{kill<-}}
##' @keywords models regression
##' @examples
##'
##' m <- lvm() ## Initialize empty lvm-object
##' ### E(y1|z,v) = beta1*z + beta2*v
##' regression(m) <- y1 ~ z + v
##' ### E(y2|x,z,v) = beta*x + beta*z + 2*v + beta3*u
##' regression(m) <- y2 ~ f(x,beta) + f(z,beta)  + f(v,2) + u
##' ### Clear restriction on association between y and
##' ### fix slope coefficient of u to beta
##' regression(m, y2 ~ v+u) <- list(NA,"beta")
##'
##' regression(m) ## Examine current linear parameter constraints
##'
##' ## ## A multivariate model, E(yi|x1,x2) = beta[1i]*x1 + beta[2i]*x2:
##' m2 <- lvm(c(y1,y2) ~ x1+x2)
##'
##' @export
"regression<-" <- function(object,...,value) UseMethod("regression<-")


##' @export
"regression<-.lvm" <- function(object, to=NULL, quick=FALSE, ..., value) {
    dots <- list(...)

    if (length(dots$additive)>0 && !dots$additive && !inherits(value,"formula")) {
        regression(object,beta=value,...) <- to
        return(object)
    }
    if (!is.null(to) || !is.null(dots$y)) {
        regfix(object, to=to, ...) <- value
        return(object)
    } else  {
        if (is.list(value)) {
            for (v in value) {
                regression(object,...) <- v
            }
            return(object)
        }

        if (inherits(value,"formula")) {

            fff <- procformula(object,value,...)
            object <- fff$object
            lhs <- fff$lhs
            xs <- fff$xs
            ys <- fff$ys
            res <- fff$res
            X <- fff$X

            
        if (fff$iscovar) {
            ## return(covariance(object,var1=decomp.specials(lhs[[1]]),var2=X))
            covariance(object) <- toformula(decomp.specials(lhs[[1]]),X)
            return(object)
        }
        if (!is.null(lhs) && nchar(lhs[[1]])>2 && substr(lhs[[1]],1,2)=="v(") {
            v <- update(value,paste(decomp.specials(lhs),"~."))
            covariance(object,...) <- v
            return(object)
        }

        if (length(lhs)==0) {
            index(object) <- reindex(object)
            return(object)
        }

        for (i in seq_len(length(ys))) {
        y <- ys[i]
        for (j in seq_len(length(xs))) {
          if (length(res[[j]])>1) {
            regfix(object, to=y[1], from=xs[j],...) <- res[[j]][2]
          } else {
            object <- regression(object,to=y[1],from=xs[j],...)
          }
        }
      }
      object$parpos <- NULL
      return(object)
    }

    if (!is.list(value) | length(value)>2) stop("Value should contain names of outcome (to) and predictors (from)")
    if (all(c("to","from")%in%names(value))) {

      xval <- value$x; yval <- value$y
    } else {
      yval <- value[[1]]; xval <- value[[2]]
    }
    regression(object, to=yval, from=xval,...)
  }
}

##' @export
`regression` <-
  function(object,to,from,...) UseMethod("regression")

##' @export
`regression.lvm` <-
    function(object=lvm(),to,from,fn=NA,silent=lava.options()$silent,
             additive=TRUE, y,x,value,...) {
        if (!missing(y)) {
            if (inherits(y,"formula")) y <- all.vars(y)
            to <- y
        }
        if (!missing(x)) {
            if (inherits(x,"formula")) x <- all.vars(x)
            from <- x
        }
        if (!additive) {
            if (!inherits(to,"formula")) to <- toformula(to,from)
            x <- attributes(getoutcome(to))$x
            K <- object$attributes$nordinal[x]
            if (is.null(K) || is.na(K)) {
                K <- list(...)$K
                if (is.null(K)) stop("Supply number of categories, K (or use method 'categorical' before calling 'regression').")
                object <- categorical(object,x,...)
            }
            dots <- list(...);
            dots$K <- K
            dots$x <- object
            dots$formula <- to
            dots$regr.only <- TRUE
            object <- do.call("categorical",dots)
            return(object)
        }

        if (missing(to)) {
            return(regfix(object))
        }
        if (inherits(to,"formula")) {
            if (!missing(value)) {
                regression(object,to,silent=silent,...) <- value
            } else {
                regression(object,silent=silent,...) <- to
            }
            object$parpos <- NULL
            return(object)
        }
        if (is.list(to)) {
            for (t in to)
                regression(object,silent=silent,...) <- t
            object$parpos <- NULL
            return(object)
        }

        sx <- strsplit(from,"@")
        xx <- sapply(sx, FUN=function(i) i[1])
        ps <- sapply(sx, FUN=function(i) i[2])
        sx <- strsplit(xx,"$",fixed=TRUE)
        xs <- sapply(sx, FUN=function(i) i[1])
        fix <- as.numeric(sapply(sx, FUN=function(i) i[2]))
        allv <- index(object)$vars

        object <- addvar(object, c(to,xs), silent=silent,reindex=FALSE)

        for (i in to)
            for (j in xs) {
                object$M[j,i] <- 1
                if (!is.na(fn))
                    functional(object,j,i) <- fn
            }

        if (lava.options()$exogenous) {
            newexo <- setdiff(xs,c(to,allv))
            exo <- exogenous(object)
            if (length(newexo)>0)
                exo <- unique(c(exo,newexo))
            exogenous(object) <- setdiff(exo,to)
        }

        if (lava.options()$debug) {
            print(object$fix)
        }
        object$fix[xs,to] <- fix
        object$par[xs,to] <- ps
        object$parpos <- NULL

        index(object) <- reindex(object)
        return(object)
    }

