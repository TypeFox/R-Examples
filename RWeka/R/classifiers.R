### Weka classifiers.

## Note that all schemes for numeric or nominal prediction (i.e.,
## classifiers) in Weka extend abstract class "Classifier", and *must*
## provide either distributionForInstance() or classifyInstance().
##
## Note also that class Classifier provides methods
##   getOptions()
##   listOptions()
##   setOptions()
## (in fact, Weka's OptionHandler interface) so that we should be able
## to safely call these methods.

make_Weka_classifier <-
function(name, class = NULL, handlers = list(), init = NULL)
{
    
    ## Return a function interfacing the Weka classification learner
    ## class 'name'.
    
    ## Eventually, add support for more handlers, including:
    ## * a formula handler (e.g., are interactions allowed? etc.)
    ## * a data handler (e.g., are numeric or categorical responses
    ##   allowed? etc.)

    ## Add to registry.
    classes <- c(class, "Weka_classifier")
    kind <- "R_Weka_classifier_interface"
    name <- as_JNI_name(name)
    meta <- make_R_Weka_interface_metadata(name, kind, classes, init)
    Weka_interfaces[[Java_class_base_name(name)]] <- meta

    ## Provide a default data handler.
    if(is.na(match("data", names(handlers))))
        handlers$data <- .default_data_handler_for_classifiers
        
    out <- function(formula, data, subset, na.action,
                    control = Weka_control(), options = NULL)
    {
        ## The "usual" way of creating a model frame from the call.
        mc <- match.call()
        mf <- mc[c(1L, match(c("formula", "data", "subset", "na.action"),
                             names(mc), 0L))]
        ## Need 'stats::' for non-standard evaluation:
        mf[[1L]] <- quote(stats::model.frame)
        mf <- eval(mf, parent.frame())
	
        .structure(c(RWeka_build_classifier(mf, control, name, handlers,
                                            options, init),
                     list(call = mc, handlers = handlers,
                          levels = levels(mf[[1L]]),
                          terms = attr(mf, "terms"))),
                   class = classes)
    }
    make_R_Weka_interface(out, meta)
}

RWeka_build_classifier <-
function(mf, control, name, handlers, options, init)
{
    if(is.function(init)) init()

    out <- list()

    options <- .expand_Weka_classifier_options(options)

    ## If needed save original model frame first, as data handlers might
    ## modify it.
    if(identical(options$model, TRUE))
        out$model <- mf

    mf <- .compose_and_funcall(handlers$data, mf)
    instances <- read_model_frame_into_Weka(mf)

    ## Build the classifier.
    classifier <- .jnew(name)
    control <- as.character(.compose_and_funcall(handlers$control,
                                                 control))
    if(length(control))
        .jcall(classifier, "V", "setOptions", .jarray(control))
    .jcall(classifier, "V", "buildClassifier", instances)

    ## And classify the training instances.
    predictions <- .predictions_for_instances(classifier, instances)
    if(!is.null(levels <- levels(mf[[1L]])))
        predictions <- factor(levels[predictions + 1L], levels = levels)
    
    out <- c(out,
             list(classifier = classifier, predictions = predictions))
    if(identical(options$instances, TRUE))
        out$instances <- instances

    out
}

.expand_Weka_classifier_options <-
function(x)
{
    x <- as.list(x)
    ## List of currently supported options with defaults.
    out <- list(model = FALSE, instances = FALSE)
    ## Now match.
    pos <- charmatch(names(x), names(out), nomatch = 0L)
    out[pos] <- x[pos != 0L]
    out
}

print.Weka_classifier <-
function(x, ...)
{
    writeLines(.jcall(x$classifier, "S", "toString"))
    invisible(x)
}

summary.Weka_classifier <-
function(object, ...)
{
    evaluate_Weka_classifier(object, ...)
}

.predictions_for_instances <-
function(classifier, instances)
{
    ## Get the predictions for a fitted Weka classifier.
   
    ## Weka uses NaN for missing values as we do in RWekaInterfaces.
    ## So we have to map to NA.
    
    if(.has_method(classifier, "classifyInstance")) {
        class <- .jcall("RWekaInterfaces", "[D",
                        "classifyInstances",
                        .jcast(classifier, "weka/classifiers/Classifier"),
                        instances)
        is.na(class) <- is.nan(class)
        class
    }
    else {
        ## If there is no classifyInstance() method, the Weka classifier
        ## must provide a distributionForInstance() method.
        .distribution_for_instances(classifier, instances)
    }
}

.distribution_for_instances <-
function(classifier, instances)
{
    ## Predict the "memberships" for given instances from a fitted Weka
    ## classifier.  (Note that in principle a numeric classifier could
    ## provide just a distributionForInstance() method which in that
    ## case would return the (numeric) predictions.)

    out <- .jcall("RWekaInterfaces", "[D",
                  "distributionForInstances",
                  .jcast(classifier, "weka/classifiers/Classifier"),
                  instances)
    matrix(out, nrow = .jcall(instances, "I", "numInstances"),
           byrow = TRUE)
}

predict.Weka_classifier <-
function(object, newdata = NULL, type = c("class", "probability"), ...)
{
    ## This should work as a general-purpose interface to getting
    ## predictions from Weka.

    type <- match.arg(type)

    if(type == "probability" && is.null(object$levels))
        stop("Can only compute class probabilities for classification problems.")
    instances <- NULL
    if(is.null(newdata)) {
        ## Currently only the class predictions for the training data
        ## are stored:
        if(type == "class") return(object$predictions)
        ## but not the probabilities.  Hence, check whether instances or
        ## the model frame are already known, or try something fancy.
	else {
            instances <- object$instances
            if(is.null(instances)) {
                if(!is.null(mf <- object$model))
                    instances <- read_model_frame_into_Weka(mf)
                else
                    newdata <-
                        eval(object$call$data,
                             environment(formula(object)))
            }
        }
    }
        
    if(is.null(instances)) {
        mf <- model.frame(delete.response(terms(object)), newdata,
                          na.action = object$call$na.action)
        mf <- .compose_and_funcall(object$handlers$data, mf)        

        ## Seems that Weka always needs to have a "class" with its
        ## instances, and even know a factor by its levels ...
        classes <- if(!is.null(object$levels))
            factor(NA, levels = object$levels)
        else {
            ## Use anything *numeric* (could also use NaN or 0).  Just
            ## 'NA' (logical by default) does not always work, as
            ## spotted by Fernando Cela Diaz <fcela@sloan.mit.edu> for
            ## learner weka/classifiers/functions/MultilayerPerceptron.
            NA_real_
        }
        mf <- cbind(CLASS = classes, mf)

        ## Get new instances into Weka.
        instances <- read_model_frame_into_Weka(mf)
    }

    switch(type,
           "class" = {
               ## Get predictions from Weka.
               out <- .predictions_for_instances(object$classifier,
                                                 instances)
               ## Post-process predictions for factors.
               if(!is.null(object$levels))
                   out <- factor(object$levels[out + 1L],
                                 levels = object$levels)    
           },
           "probability" = {
               ## Get predictions from Weka.
               out <- .distribution_for_instances(object$classifier,
                                                  instances)
               ## No Weka_instances class and dimnames method for now.
               rnms <- attr(instances, ".dimnames")[[1L]]
               dimnames(out) <- list(rnms, object$levels)    
           })
    
    out
}

fitted.Weka_classifier <-
function (object, ...) 
{
    predict(object, ...)
}

model.frame.Weka_classifier <-
function(formula, ...)
{
    if(!is.null(mf <- formula$model)) return(mf)
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"),
                        names(dots), 0L)]
    mf <- formula$call
    mf <- mf[c(1L, match(c("formula", "data", "subset", "na.action"),
                         names(mf), 0L))]
    mf$drop.unused.levels <- TRUE
    ## Need 'stats::' for non-standard evaluation:
    mf[[1L]] <- quote(stats::model.frame)
    mf[names(nargs)] <- nargs
    if(is.null(env <- environment(formula$terms)))
        env <- parent.frame()
    eval(mf, env)
}

.default_data_handler_for_classifiers <-
function(mf)
{
    ## A default data handler for classifiers which rejects interaction
    ## terms and drops unused variables, so that e.g.
    ##   J48(Species ~ . - Petal.Width, data = iris)
    ## works "as expected".
    ## Issue raised by David Gleich <dgleich@stanford.edu>.
    
    terms <- attr(mf, "terms")
    ## No interactions.
    if(any(attr(terms, "order") > 1L))
        stop("Interactions are not allowed.")
    factors <- attr(terms, "factors")
    varnms <- rownames(factors)[c(TRUE, rowSums(factors)[-1L] > 0)]
    ## Remove backticks from non-syntactic names.
    mf[, sub("^`(.*)`$", "\\1", varnms), drop = FALSE]
}
