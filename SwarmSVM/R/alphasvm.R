#' Support Vector Machines taking initial alpha values
#' 
#' alphasvm is used to train a support vector machine. It can be used to carry out general regression and classification (of nu and epsilon-type), as well as density-estimation. A formula interface is provided. 
#' 
#' @param formula a symbolic description of the model to be fit.
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which 'svm' is called from.
#' @param x a data matrix, a vector, or a sparse matrix (object of class \code{Matrix} 
#'     provided by the \code{Matrix} package, or of class \code{matrix.csr} provided by the \code{SparseM} package, 
#'     or of class \code{simple_triplet_matrix} provided by the \code{slam package}).
#' @param y a response vector with one label for each row/component of x. Can be either a factor (for classification tasks) or a numeric vector (for regression).
#' @param scale A logical vector indicating the variables to be scaled. 
#'     If scale is of length 1, the value is recycled as many times as needed. 
#'     Per default, data are scaled internally (both x and y variables) to zero mean and unit variance. 
#'     The center and scale values are returned and used for later predictions.
#' @param type svm can be used as a classification machine. 
#'     The default setting for type is C-classification, 
#'     but may be set to nu-classification as well.
#' @param kernel the kernel used in training and predicting. You might consider 
#'     changing some of the following parameters, depending on the kernel type.
#' \describe{
#'     \item{linear:}{\eqn{u'v}{u'*v}}
#'     \item{polynomial:}{\eqn{(\gamma u'v + coef0)^{degree}}{(gamma*u'*v + coef0)^degree}}
#'     \item{radial basis:}{\eqn{e^(-\gamma |u-v|^2)}{exp(-gamma*|u-v|^2)}}
#'     \item{sigmoid:}{\eqn{tanh(\gamma u'v + coef0)}{tanh(gamma*u'*v + coef0)}}
#' }
#' @param degree parameter needed for kernel of type \code{polynomial} (default: 3)
#' @param gamma parameter needed for all kernels except \code{linear} 
#'     (default: 1/(data dimension))
#' @param coef0 parameter needed for kernels of type \code{polynomial} 
#'     and \code{sigmoid} (default: 0)
#' @param cost cost of constraints violation (default: 1)---it is the 
#'     \sQuote{C}-constant of the regularization term in the Lagrange formulation.
#' @param nu parameter needed for \code{nu-classification}
#' @param class.weights a named vector of weights for the different
#'     classes, used for asymmetric class sizes. Not all factor levels have
#'     to be supplied (default weight: 1). All components have to be named.
#' @param cachesize cache memory in MB (default 40)
#' @param tolerance tolerance of termination criterion (default: 0.001)
#' @param epsilon epsilon in the insensitive-loss function (default: 0.1)
#' @param shrinking option whether to use the shrinking-heuristics (default: \code{TRUE})
#' @param cross if a integer value k>0 is specified, a k-fold cross 
#'     validation on the training data is performed to assess the quality 
#'     of the model: the accuracy rate for classification and the Mean 
#'     Squared Error for regression
#' @param fitted logical indicating whether the fitted values should be computed 
#'     and included in the model or not (default: \code{TRUE})
#' @param alpha Initial values for the coefficients (default: \code{NULL}). 
#'     A numerical vector for binary classification or a nx(k-1) matrix for a k-class-classification problem.
#' @param mute a logical value indicating whether to print training information from svm.
#' @param probability logical indicating whether the model should 
#'     allow for probability predictions.
#' @param ... additional parameters for the low level fitting function \code{svm.default}
#' @param subset An index vector specifying the cases to be used in the 
#'     training sample.  (NOTE: If given, this argument must be named.)
#' @param na.action A function to specify the action to be taken if \code{NA}s are 
#'     found. The default action is \code{stats::na.omit}, which leads to rejection of cases  
#'     with missing values on any required variable. An alternative 
#'     is \code{stats::na.fail}, which causes an error if \code{NA} cases 
#'     are found. (NOTE: If given, this argument must be named.)
#' 
#' @details
#'     For multiclass-classification with k levels, k>2, \code{libsvm} uses the
#'     \sQuote{one-against-one}-approach, in which k(k-1)/2 binary classifiers are
#'     trained; the appropriate class is found by a voting scheme.
#' 
#'     \code{libsvm} internally uses a sparse data representation, which is 
#'     also high-level supported by the package \pkg{SparseM}.
#' 
#'     If the predictor variables include factors, the formula interface must be used to get a
#'     correct model matrix.
#' 
#'     \code{plot.svm} allows a simple graphical
#'     visualization of classification models.
#' 
#'     The probability model for classification fits a logistic distribution
#'     using maximum likelihood to the decision values of all binary
#'     classifiers, and computes the a-posteriori class probabilities for the
#'     multi-class problem using quadratic optimization. The probabilistic
#'     regression model assumes (zero-mean) laplace-distributed errors for the
#'     predictions, and estimates the scale parameter using maximum likelihood.
#' 
#' @references
#' 
#' \itemize{
#'     \item
#'         Chang, Chih-Chung and Lin, Chih-Jen:\cr
#'         \emph{LIBSVM: a library for Support Vector Machines}\cr
#'         \url{http://www.csie.ntu.edu.tw/~cjlin/libsvm}
#'     
#'     \item
#'         Exact formulations of models, algorithms, etc. can be found in the
#'         document:\cr
#'         Chang, Chih-Chung and Lin, Chih-Jen:\cr
#'         \emph{LIBSVM: a library for Support Vector Machines}\cr
#'         \url{http://www.csie.ntu.edu.tw/~cjlin/papers/libsvm.ps.gz}
#'     
#'     \item
#'         More implementation details and speed benchmarks can be found on:
#'         Rong-En Fan and Pai-Hsune Chen and Chih-Jen Lin:\cr
#'         \emph{Working Set Selection Using the Second Order Information for Training SVM}\cr
#'         \url{http://www.csie.ntu.edu.tw/~cjlin/papers/quadworkset.pdf}
#' }
#' 
#' @author
#' Tong He (based on package \code{e1071} by David Meyer and C/C++ code by Cho-Jui Hsieh in Divide-and-Conquer kernel SVM (DC-SVM) )
#' 
#' @examples
#' 
#' data(svmguide1)
#' svmguide1.t = svmguide1[[2]]
#' svmguide1 = svmguide1[[1]]
#' 
#' model = alphasvm(x = svmguide1[,-1], y = svmguide1[,1], scale = TRUE)
#' preds = predict(model, svmguide1.t[,-1])
#' table(preds, svmguide1.t[,1])
#' 
#' data(iris)
#' attach(iris)
#' 
#' # default with factor response:
#' model = alphasvm(Species ~ ., data = iris)
#' 
#' # get new alpha
#' new.alpha = matrix(0, nrow(iris),2)
#' new.alpha[model$index,] = model$coefs
#' 
#' model2 = alphasvm(Species ~ ., data = iris, alpha = new.alpha)
#' preds = predict(model2, as.matrix(iris[,-5]))
#' table(preds, iris[,5])
#' 
#' @rdname alphasvm
#' 
#' 
#' @export
#' 
alphasvm <-
function (x, ...)
    UseMethod ("alphasvm")
#'
#' @rdname alphasvm
#' @export
#' 
alphasvm.formula <-
function (formula, data = NULL, ..., subset, na.action = stats::na.omit, scale = FALSE)
{
    call <- match.call()
    if (!inherits(formula, "formula"))
        stop("method is only for formula objects")
    m <- match.call(expand.dots = FALSE)
    if (identical(class(eval.parent(m$data)), "matrix"))
        m$data <- as.data.frame(eval.parent(m$data))
    m$... <- NULL
    m$scale <- NULL
    m[[1]] <- as.name("model.frame")
    m$na.action <- na.action
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 0
    x <- stats::model.matrix(Terms, m)
    y <- stats::model.extract(m, "response")
    attr(x, "na.action") <- attr(y, "na.action") <- attr(m, "na.action")
    if (length(scale) == 1)
        scale <- rep(scale, ncol(x))
    if (any(scale)) {
        remove <- unique(c(which(labels(Terms) %in%
                                 names(attr(x, "contrasts"))),
                           which(!scale)
                           )
                         )
        scale <- !attr(x, "assign") %in% remove
    }
    ret <- alphasvm.default (x, y, scale = scale, ..., na.action = na.action)
    ret$call <- call
    ret$call[[1]] <- as.name("alphasvm")
    ret$terms <- Terms
    if (!is.null(attr(m, "na.action")))
        ret$na.action <- attr(m, "na.action")
    class(ret) <- c("alphasvm.formula", class(ret))
    return (ret)
}

#' 
#' @rdname alphasvm
#' @export
#' 
alphasvm.default <-
function (x,
          y           = NULL,
          scale       = FALSE,
          type        = NULL,
          kernel      = "radial",
          degree      = 3,
          gamma       = if (is.vector(x)) 1 else 1 / ncol(x),
          coef0       = 0,
          cost        = 1,
          nu          = 0.5,
          class.weights = NULL,
          cachesize   = 40,
          tolerance   = 0.001,
          epsilon     = 0.1,
          shrinking   = TRUE,
          cross       = 0,
          probability = FALSE,
          fitted      = TRUE,
          alpha       = NULL,
          mute        = TRUE,
          ...,
          subset,
          na.action = stats::na.omit)
{
    if(inherits(x, "Matrix")) {
        x <- as(x, "matrix.csr")
    }
    if(inherits(x, "simple_triplet_matrix")) {
        ind <- order(x$i, x$j)
        x <- new("matrix.csr",
                 ra = x$v[ind],
                 ja = x$j[ind],
                 ia = as.integer(cumsum(c(1, tabulate(x$i[ind])))),
                 dimension = c(x$nrow, x$ncol))
    }
    sparse <- inherits(x, "matrix.csr")
    
    assertInt(nrow(x), lower = 1)
    assertInt(ncol(x), lower = 1)

    ## NULL parameters?
    if(is.null(degree)) stop(sQuote("degree"), " must not be NULL!")
    if(is.null(gamma)) stop(sQuote("gamma"), " must not be NULL!")
    if(is.null(coef0)) stop(sQuote("coef0"), " must not be NULL!")
    if(is.null(cost)) stop(sQuote("cost"), " must not be NULL!")
    if(is.null(nu)) stop(sQuote("nu"), " must not be NULL!")
    if(is.null(epsilon)) stop(sQuote("epsilon"), " must not be NULL!")
    if(is.null(tolerance)) stop(sQuote("tolerance"), " must not be NULL!")

    xhold   <- if (fitted) x else NA
    x.scale <- y.scale <- NULL
    formula <- inherits(x, "svm.formula")

    if (!is.factor(y)) {
      y = as.factor(y)
      nlvl = length(levels(y))
      if (nlvl>16) 
        stop("Only classification missions with no more than 16 target classes are allowed.")
    }
    assertFactor(y, len = nrow(x), max.levels = 16)

    # determine model type
    if (is.null(type)) type <-
        if (is.null(y)) "one-classification"
        else if (is.factor(y)) "C-classification"
        else "eps-regression"

    #     if (is.factor(y) && length(levels(y))==1)
    #         type = "one-classification"
        
    #type = "C-classification"
    type <- pmatch(type, c("C-classification",
                           "nu-classification"), 99) - 1

    if (type > 10) stop("wrong type specification!")
    

    kernel <- pmatch(kernel, c("linear",
                               "polynomial",
                               "radial",
                               "sigmoid"), 99) - 1

    if (kernel > 10) stop("wrong kernel specification!")

    nac <- attr(x, "na.action")

    ## scaling, subsetting, and NA handling
    if (sparse) {
        scale <- rep(FALSE, ncol(x))
        if(!is.null(y)) stats::na.fail(y)
        x <- SparseM::t(SparseM::t(x)) ## make shure that col-indices are sorted
    } else {
        x <- as.matrix(x)

        ## subsetting and na-handling for matrices
        if (!formula) {
            if (!missing(subset)) x <- x[subset,]
            if (is.null(y))
                x <- na.action(x)
            else {
                df <- na.action(data.frame(y, x))
                y <- df[,1]
                x <- as.matrix(df[,-1])
                nac <-
                    attr(x, "na.action") <-
                        attr(y, "na.action") <-
                            attr(df, "na.action")
            }
        }

        ## scaling
        if (length(scale) == 1)
            scale <- rep(scale, ncol(x))
        if (any(scale)) {
            co <- !apply(x[,scale, drop = FALSE], 2, stats::var)
            if (any(co)) {
                warning(paste("Variable(s)",
                              paste(sQuote(colnames(x[,scale,
                                                      drop = FALSE])[co]),
                                    sep="", collapse=" and "),
                              "constant. Cannot scale data.")
                        )
                scale <- rep(FALSE, ncol(x))
            } else {
                xtmp <- scale(x[,scale])
                x[,scale] <- xtmp
                x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
                if (is.numeric(y) && (type > 2)) {
                    y <- scale(y)
                    y.scale <- attributes(y)[c("scaled:center","scaled:scale")]
                    y <- as.vector(y)
                }
            }
        }
    }

    ## further parameter checks
    nr <- nrow(x)
    if (cross > nr)
        stop(sQuote("cross"), " cannot exceed the number of observations!")

    if (!is.vector(y) && !is.factor (y) && type != 2)
        stop("y must be a vector or a factor.")
    if (type != 2 && length(y) != nr)
        stop("x and y don't match.")

    if (cachesize < 0.1)
        cachesize <- 0.1

    if (type > 2 && !is.numeric(y))
        stop("Need numeric dependent variable for regression.")

    # Type and format check of alpha
    nclass = length(unique(y))
    if (testNull(alpha)) {
      if (nclass>1) {
        alpha = matrix(0, nrow(x), nclass-1)
      } else {
        alpha = matrix(0, nrow(x), 1)
      }
    } 
    if (testMatrix(alpha))
        alpha = as.vector(alpha)
    # Key Point
    alpha = abs(alpha)
    assertVector(alpha)
    if (nclass>1) {
      assertNumeric(alpha, len = nrow(x) * (nclass-1))
    } else {
      assertNumeric(alpha, len = nrow(x))
    }
    alphalen = as.integer(length(alpha))
    assertInt(alphalen)
    
    
    lev <- NULL
    weightlabels <- NULL

    ## in case of classification: transform factors into integers
    if (type == 2) # one class classification --> set dummy
        y <- rep(1, nr)
    else
        if (is.factor(y)) {
            lev <- levels(y)
            y <- as.integer(y)
            if (!is.null(class.weights)) {
                if (is.null(names(class.weights)))
                    stop("Weights have to be specified along with their according level names !")
                weightlabels <- match (names(class.weights), lev)
                if (any(is.na(weightlabels)))
                    stop("At least one level name is missing or misspelled.")
            }
        } else {
            if (type < 3) {
                if(any(as.integer(y) != y))
                    stop("dependent variable has to be of factor or integer type for classification mode.")
                y <- as.factor(y)
                lev <- levels(y)
                y <- as.integer(y)
            } else lev <- unique(y)
        }

    # nclass <- 2
    # if (type < 2) nclass <- length(lev)

    if (type > 1 && length(class.weights) > 0) {
        class.weights <- NULL
        warning(sQuote("class.weights"), " are set to NULL for regression mode. For classification, use a _factor_ for ", sQuote("y"),
", or specify the correct ", sQuote("type"), " argument.")
    }

    err <- empty_string <- paste(rep(" ", 255), collapse = "")

    if (is.null(type)) stop("type argument must not be NULL!")
    if (is.null(kernel)) stop("kernel argument must not be NULL!")
    if (is.null(degree)) stop("degree argument must not be NULL!")
    if (is.null(gamma)) stop("gamma argument must not be NULL!")
    if (is.null(coef0)) stop("coef0 seed argument must not be NULL!")
    if (is.null(cost)) stop("cost argument must not be NULL!")
    if (is.null(nu)) stop("nu argument must not be NULL!")
    if (is.null(cachesize)) stop("cachesize argument must not be NULL!")
    if (is.null(tolerance)) stop("tolerance argument must not be NULL!")
    if (is.null(epsilon)) stop("epsilon argument must not be NULL!")
    if (is.null(shrinking)) stop("shrinking argument must not be NULL!")
    if (is.null(cross)) stop("cross argument must not be NULL!")
    if (is.null(sparse)) stop("sparse argument must not be NULL!")
    if (is.null(probability)) stop("probability argument must not be NULL!")
    
    muteFun({
      cret <- .C ("svmtrain",
                  ## data
                  as.double  (if (sparse) x@ra else t(x)),
                  as.integer (nr), as.integer(ncol(x)),
                  as.double  (y),
                  ## sparse index info
                  as.integer (if (sparse) x@ia else 0),
                  as.integer (if (sparse) x@ja else 0),
  
                  ## parameters
                  as.integer (type),
                  as.integer (kernel),
                  as.integer (degree),
                  as.double  (gamma),
                  as.double  (coef0),
                  as.double  (cost),
                  as.double  (nu),
                  as.integer (weightlabels),
                  as.double  (class.weights),
                  as.integer (length (class.weights)),
                  as.double  (cachesize),
                  as.double  (tolerance),
                  as.double  (epsilon),
                  as.integer (shrinking),
                  as.integer (cross),
                  as.integer (sparse),
                  as.integer (probability),
                  as.double  (alpha),
                  as.integer (alphalen),
  
                  ## results
                  nclasses = integer  (1),
                  nr       = integer  (1), # nr of support vectors
                  index    = integer  (nr),
                  labels   = integer  (nclass),
                  nSV      = integer  (nclass),
                  rho      = double   (nclass * (nclass - 1) / 2),
                  coefs    = double   (nr * (nclass - 1)),
                  sigma    = double   (1),
                  probA    = double   (nclass * (nclass - 1) / 2),
                  probB    = double   (nclass * (nclass - 1) / 2),
  
                  cresults = double   (cross),
                  ctotal1  = double   (1),
                  ctotal2  = double   (1),
                  error    = err,
  
                  PACKAGE = "SwarmSVM")
    }, mute = mute)

    if (cret$error != empty_string)
        stop(paste(cret$error, "!", sep=""))
    if (!is.null(cret$nr) && !is.na(cret$nr)) {
      cret$index <- cret$index[1:cret$nr]
    } else {
      cret$index = NULL
    }
    if (is.null(cret$index) || cret$index == 0) {
      SV = NULL
    } else {
      if (sparse) {
        SV = SparseM::t(SparseM::t(x[cret$index,]))
      } else {
        SV = t(t(x[cret$index,]))
      }
    }

    ret <- list (
                 call     = match.call(),
                 type     = type,
                 kernel   = kernel,
                 cost     = cost,
                 degree   = degree,
                 gamma    = gamma,
                 coef0    = coef0,
                 nu       = nu,
                 epsilon  = epsilon,
                 sparse   = sparse,
                 scaled   = scale,
                 x.scale  = x.scale,
                 y.scale  = y.scale,

                 nclasses = cret$nclasses, #number of classes
                 levels   = lev,
                 tot.nSV  = cret$nr, #total number of sv
                 nSV      = cret$nSV[1:cret$nclasses], #number of SV in diff. classes
                 labels   = cret$label[1:cret$nclasses], #labels of the SVs.
                 SV       = SV, #copy of SV
                 index    = cret$index,  #indexes of sv in x
                 ##constants in decision functions
                 rho      = cret$rho[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 ##probabilites
                 compprob = probability,
                 probA    = if (!probability) NULL else
                 cret$probA[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 probB    = if (!probability) NULL else
                 cret$probB[1:(cret$nclasses * (cret$nclasses - 1) / 2)],
                 sigma    = if (probability) cret$sigma else NULL,
                 ##coefficiants of sv
                 coefs    = if (cret$nr == 0) NULL else
                 t(matrix(cret$coefs[1:((cret$nclasses - 1) * cret$nr)],
                          nrow = cret$nclasses - 1,
                          byrow = TRUE)),
                 na.action = nac
                 )

    ## cross-validation-results
    if (cross > 0)
        if (type > 2) {
            scale.factor     <- if (any(scale)) crossprod(y.scale$"scaled:scale") else 1;
            ret$MSE          <- cret$cresults * scale.factor;
            ret$tot.MSE      <- cret$ctotal1  * scale.factor;
            ret$scorrcoeff   <- cret$ctotal2;
        } else {
            ret$accuracies   <- cret$cresults;
            ret$tot.accuracy <- cret$ctotal1;
        }

    class (ret) <- "alphasvm"

    if (fitted) {
        if (!is.null(ret$coefs)) {
          ret$fitted <- na.action(predict(ret, xhold,
                                          decision.values = TRUE))
          ret$decision.values <- attr(ret$fitted, "decision.values")
          attr(ret$fitted, "decision.values") <- NULL
          if (type > 1) ret$residuals <- y - ret$fitted
        } else {
          ret$fitted = y
        }
    }

    ret
}

#' Prediction function for an alphasvm object
#' 
#' @param object the object trained from \code{alphasvm}
#' @param newdata the test data set
#' @param decision.values a logical variable indicating whether to output the decision values
#' @param probability a logical variable indicating whether to output the classfication probability
#' @param ... currently not used
#' @param na.action A function to specify the action to be taken if 'NA's are found. 
#'     The default action is \code{stats::na.omit}, which leads to rejection of cases with missing values on any required variable. 
#'     An alternative is \code{stats::na.fail}, which causes an error if NA cases are found. (NOTE: If given, this argument must be named.)
#' 
#' @method predict alphasvm
#' 
#' @export
#' 
predict.alphasvm <-
function (object, newdata,
          decision.values = FALSE,
          probability = FALSE,
          ...,
          na.action = stats::na.omit)
{
    if (missing(newdata))
        return(fitted(object))

    if (object$tot.nSV < 1) {
        #stop("Model is empty!")
      assertInt(nrow(newdata), lower = 1)
      n = nrow(newdata)
      res = rep(object$fitted[1], n)
      return(res)
    }


    if(inherits(newdata, "Matrix")) {
        newdata <- as(newdata, "matrix.csr")
    }
    if(inherits(newdata, "simple_triplet_matrix")) {
       ind <- order(newdata$i, newdata$j)
       newdata <- new("matrix.csr",
                      ra = newdata$v[ind],
                      ja = newdata$j[ind],
                      ia = as.integer(cumsum(c(1, tabulate(newdata$i[ind])))),
                      dimension = c(newdata$nrow, newdata$ncol))
   }

    sparse <- inherits(newdata, "matrix.csr")

    act <- NULL
    if ((is.vector(newdata) && is.atomic(newdata)))
        newdata <- t(t(newdata))
    if (sparse)
        newdata <- SparseM::t(SparseM::t(newdata))
    preprocessed <- !is.null(attr(newdata, "na.action"))
    rowns <- if (!is.null(rownames(newdata)))
        rownames(newdata)
    else
        1:nrow(newdata)
    if (!object$sparse) {
        if (inherits(object, "alphasvm.formula")) {
            if(is.null(colnames(newdata)))
                colnames(newdata) <- colnames(object$SV)
            newdata <- na.action(newdata)
            act <- attr(newdata, "na.action")
            newdata <- stats::model.matrix(stats::delete.response(stats::terms(object)),
                                    as.data.frame(newdata))
        } else {
            newdata <- na.action(as.matrix(newdata))
            act <- attr(newdata, "na.action")
        }
    }

    if (!is.null(act) && !preprocessed)
        rowns <- rowns[-act]

    if (any(object$scaled))
        newdata[,object$scaled] <-
            scale(newdata[,object$scaled, drop = FALSE],
                  center = object$x.scale$"scaled:center",
                  scale  = object$x.scale$"scaled:scale"
                  )

    if (ncol(object$SV) != ncol(newdata))
        stop ("test data does not match model !")

    ret <- .C ("svmpredict",
               as.integer (decision.values),
               as.integer (probability),

               ## model
               as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.integer (object$compprob),
               as.double  (if (object$compprob) object$probA else 0),
               as.double  (if (object$compprob) object$probB else 0),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),

               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),

               ## test matrix
               as.double  (if (sparse) newdata@ra else t(newdata)),
               as.integer (nrow(newdata)),
               as.integer (if (sparse) newdata@ia else 0),
               as.integer (if (sparse) newdata@ja else 0),
               as.integer (sparse),

               ## decision-values
               ret = double(nrow(newdata)),
               dec = double(nrow(newdata) * object$nclasses * (object$nclasses - 1) / 2),
               prob = double(nrow(newdata) * object$nclasses),

               PACKAGE = "SwarmSVM"
               )

    ret2 <- if (is.character(object$levels)) # classification: return factors
        factor (object$levels[ret$ret], levels = object$levels)
    else if (object$type == 2) # one-class-classification: return TRUE/FALSE
        ret$ret == 1
    else if (any(object$scaled) && !is.null(object$y.scale)) # return raw values, possibly scaled back
        ret$ret * object$y.scale$"scaled:scale" + object$y.scale$"scaled:center"
    else
        ret$ret

    names(ret2) <- rowns
    ret2 <- stats::napredict(act, ret2)

    if (decision.values) {
        colns = c()
        for (i in 1:(object$nclasses - 1))
            for (j in (i + 1):object$nclasses)
                colns <- c(colns,
                           paste(object$levels[object$labels[i]],
                                 "/", object$levels[object$labels[j]],
                                 sep = ""))
        attr(ret2, "decision.values") <-
            stats::napredict(act,
                      matrix(ret$dec, nrow = nrow(newdata), byrow = TRUE,
                             dimnames = list(rowns, colns)
                             )
                      )
    }

    if (probability && object$type < 2) {
        if (!object$compprob)
            warning("SVM has not been trained using `probability = TRUE`, probabilities not available for predictions.")
        else
            attr(ret2, "probabilities") <-
                stats::napredict(act,
                          matrix(ret$prob, nrow = nrow(newdata), byrow = TRUE,
                                 dimnames = list(rowns, object$levels[object$labels])
                                 )
                          )
    }

    ret2
}

#' Print alphasvm object
#' 
#' @method print alphasvm
#' @rdname alphasvm
#' @export
#' 
print.alphasvm <-
function (x, ...)
{
    cat("\nCall:", deparse(x$call, 0.8 * getOption("width")), "\n", sep="\n")
    cat("Parameters:\n")
    cat("   SVM-Type: ", c("C-classification",
                           "nu-classification",
                           "one-classification",
                           "eps-regression",
                           "nu-regression")[x$type+1], "\n")
    cat(" SVM-Kernel: ", c("linear",
                           "polynomial",
                           "radial",
                           "sigmoid")[x$kernel+1], "\n")
    if (x$type==0 || x$type==3 || x$type==4)
        cat("       cost: ", x$cost, "\n")
    if (x$kernel==1)
        cat("     degree: ", x$degree, "\n")
    cat("      gamma: ", x$gamma, "\n")
    if (x$kernel==1 || x$kernel==3)
        cat("     coef.0: ", x$coef0, "\n")
    if (x$type==1 || x$type==2 || x$type==4)
        cat("         nu: ", x$nu, "\n")
    if (x$type==3) {
        cat("    epsilon: ", x$epsilon, "\n\n")
        if (x$compprob)
            cat("Sigma: ", x$sigma, "\n\n")
    }

    cat("\nNumber of Support Vectors: ", x$tot.nSV)
    cat("\n\n")

}

#' Summary alphasvm object
#' 
#' @method summary alphasvm
#' @rdname alphasvm
#' @export
#' 
summary.alphasvm <-
function(object, ...)
    structure(object, class="summary.alphasvm")

#' Print summary.alphasvm object
#' 
#' @param object An object of class \code{alphasvm}
#' 
#' @method print summary.alphasvm
#' @rdname alphasvm
#' @export
#' 
print.summary.alphasvm <-
function (x, ...)
{
    print.alphasvm(x)
    if (x$type<2) {
        cat(" (", x$nSV, ")\n\n")
        cat("\nNumber of Classes: ", x$nclasses, "\n\n")
        cat("Levels:", if(is.numeric(x$levels)) "(as integer)", "\n", x$levels)
    }
    cat("\n\n")
    if (x$type==2) cat("\nNumber of Classes: 1\n\n\n")

    if ("MSE" %in% names(x)) {
        cat(length (x$MSE), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Mean Squared Error:", x$tot.MSE, "\n")
        cat("Squared Correlation Coefficient:", x$scorrcoef, "\n")
        cat("Mean Squared Errors:\n", x$MSE, "\n\n")
    }
    if ("accuracies" %in% names(x)) {
        cat(length (x$accuracies), "-fold cross-validation on training data:\n\n", sep="")
        cat("Total Accuracy:", x$tot.accuracy, "\n")
        cat("Single Accuracies:\n", x$accuracies, "\n\n")
    }
    cat("\n\n")
}

scale.data.frame <-
function(x, center = TRUE, scale = TRUE)
{
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
        x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], stats::na.omit(center), stats::na.omit(scale))
        if(center || !is.logical(center))
            attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
        if(scale || !is.logical(scale))
            attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
}

#' Plot alphasvm object
#' 
#' @param x An object of class \code{alphasvm}
#' @param data data to visualize. Should be the same used for fitting.
#' @param formula formula selecting the visualized two dimensions. Only needed if more than two input variables are used.
#' @param fill switch indicating whether a contour plot for the class regions should be added.
#' @param grid granularity for the contour plot.
#' @param slice a list of named values for the dimensions held constant (only needed if more than two variables are used). 
#'     The defaults for unspecified dimensions are 0 (for numeric variables) and the first level (for factors). 
#'     Factor levels can either be specified as factors or character vectors of length 1.
#' @param symbolPalette Color palette used for the class the data points and support vectors belong to.
#' @param svSymbol Symbol used for support vectors.
#' @param dataSymbol Symbol used for data points (other than support vectors).
#' @param ... additional graphics parameters passed to \code{graphics::filled.contour} and \code{plot}.
#' 
#' @method plot alphasvm
#' 
#' @export
#' 
plot.alphasvm <-
function(x, data, formula = NULL, fill = TRUE,
         grid = 50, slice = list(), symbolPalette = grDevices::palette(),
         svSymbol = "x", dataSymbol = "o", ...)
{
    if (x$type < 3) {
        if (is.null(formula) && ncol(data) == 3) {
            formula <- formula(stats::delete.response(stats::terms(x)))
            formula[2:3] <- formula[[2]][2:3]
        }
        if (is.null(formula))
            stop("missing formula.")
        if (fill) {
            sub <- stats::model.frame(formula, data)
            xr <- seq(min(sub[, 2]), max(sub[, 2]), length = grid)
            yr <- seq(min(sub[, 1]), max(sub[, 1]), length = grid)
            l <- length(slice)
            if (l < ncol(data) - 3) {
                slnames <- names(slice)
                slice <- c(slice, rep(list(0), ncol(data) - 3 -
                                      l))
                names <- labels(stats::delete.response(stats::terms(x)))
                names(slice) <- c(slnames, names[!names %in%
                                                 c(colnames(sub), slnames)])
            }
            for (i in names(which(sapply(data, is.factor))))
                if (!is.factor(slice[[i]])) {
                    levs <- levels(data[[i]])
                    lev <- if (is.character(slice[[i]])) slice[[i]] else levs[1]
                    fac <- factor(lev, levels = levs)
                    if (is.na(fac))
                        stop(paste("Level", dQuote(lev), "could not be found in factor", sQuote(i)))
                    slice[[i]] <- fac
                }

            lis <- c(list(yr), list(xr), slice)
            names(lis)[1:2] <- colnames(sub)
            new <- expand.grid(lis)[, labels(stats::terms(x))]
            preds <- predict(x, new)
            graphics::filled.contour(xr, yr,
                           matrix(as.numeric(preds),
                                  nrow = length(xr), byrow = TRUE),
                           plot.axes = {
                               graphics::axis(1)
                               graphics::axis(2)
                               colind <- as.numeric(stats::model.response(stats::model.frame(x, data)))
                               dat1 <- data[-x$index,]
                               dat2 <- data[x$index,]
                               coltmp1 <- symbolPalette[colind[-x$index]]
                               coltmp2 <- symbolPalette[colind[x$index]]
                               graphics::points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
                               graphics::points(formula, data = dat2, pch = svSymbol, col = coltmp2)
                           },
                           levels = 1:(length(levels(preds)) + 1),
                           key.axes = graphics::axis(4, 1:(length(levels(preds))) + 0.5,
                           labels = levels(preds),
                           las = 3),
                           plot.title = graphics::title(main = "SVM classification plot",
                           xlab = names(lis)[2], ylab = names(lis)[1]),
                           ...)
        }
        else {
            plot(formula, data = data, type = "n", ...)
            colind <- as.numeric(stats::model.response(stats::model.frame(x,
                                                            data)))
            dat1 <- data[-x$index,]
            dat2 <- data[x$index,]
            coltmp1 <- symbolPalette[colind[-x$index]]
            coltmp2 <- symbolPalette[colind[x$index]]
            graphics::points(formula, data = dat1, pch = dataSymbol, col = coltmp1)
            graphics::points(formula, data = dat2, pch = svSymbol, col = coltmp2)
            invisible()
        }
    }
}

#' Write alphasvm object
#' 
#' @param object Object of class "alphasvm", created by \code{alphasvm}.
#' @param svm.file filename to export the \code{alphasvm} object to.
#' @param scale.file filename to export the scaling data of the explanatory variables to.
#' @param yscale.file filename to export the scaling data of the dependent variable to, if any.
#' 
#' @method write alphasvm
#' 
#' @export
#' 
write.alphasvm <-
function (object, svm.file="Rdata.svm", scale.file = "Rdata.scale",
          yscale.file = "Rdata.yscale")
{

    ret <- .C ("svmwrite",
               ## model
               as.double  (if (object$sparse) object$SV@ra else t(object$SV)),
               as.integer (nrow(object$SV)), as.integer(ncol(object$SV)),
               as.integer (if (object$sparse) object$SV@ia else 0),
               as.integer (if (object$sparse) object$SV@ja else 0),
               as.double  (as.vector(object$coefs)),
               as.double  (object$rho),
               as.integer (object$compprob),
               as.double  (if (object$compprob) object$probA else 0),
               as.double  (if (object$compprob) object$probB else 0),
               as.integer (object$nclasses),
               as.integer (object$tot.nSV),
               as.integer (object$labels),
               as.integer (object$nSV),
               as.integer (object$sparse),

               ## parameter
               as.integer (object$type),
               as.integer (object$kernel),
               as.integer (object$degree),
               as.double  (object$gamma),
               as.double  (object$coef0),

               ## filename
               as.character(svm.file),

               PACKAGE = "SwarmSVM"
               )$ret

    utils::write.table(data.frame(center = object$x.scale$"scaled:center",
                           scale  = object$x.scale$"scaled:scale"),
                file=scale.file, col.names=FALSE, row.names=FALSE)

    if (!is.null(object$y.scale))
        utils::write.table(data.frame(center = object$y.scale$"scaled:center",
                               scale  = object$y.scale$"scaled:scale"),
                    file=yscale.file, col.names=FALSE, row.names=FALSE)
}
