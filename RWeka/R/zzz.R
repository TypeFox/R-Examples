### * Associators
Apriori <-
    make_Weka_associator("weka/associations/Apriori", "Apriori")
Tertius <-
    make_Weka_associator("weka/associations/Tertius", "Tertius",
                         init = make_Weka_package_loader("tertius"))

### * Attribute evaluators

GainRatioAttributeEval <-
    make_Weka_attribute_evaluator("weka/attributeSelection/GainRatioAttributeEval")
InfoGainAttributeEval <-
    make_Weka_attribute_evaluator("weka/attributeSelection/InfoGainAttributeEval")

### * Classifiers

.make_class_name_expander <-
function(packages = NULL)
    function(e) {
        if(is.list(e))
            e[[1L]] <- get_Java_class(e[[1L]], packages)
        else if(is.character(e)) {
            ## Could have seen
            ##   Weka_control(K = "RBFKernel -G 2")
            ## or better/worse ...
            e <- unlist(strsplit(e, "[[:space:]]+"))
            e[1L] <- get_Java_class(e[1L], packages)
            e <- paste(e, collapse = " ")
        }
        e
    }
## Basic class name expander:
.expand_class_name <- .make_class_name_expander()

### ** Functions
LinearRegression <-
    make_Weka_classifier("weka/classifiers/functions/LinearRegression",
                         c("LinearRegression", "Weka_functions"))
Logistic <-
    make_Weka_classifier("weka/classifiers/functions/Logistic",
                         c("Logistic", "Weka_functions"))
## Specialized support vector kernel class name expander:
.expand_kernel_class_name <-
    .make_class_name_expander("weka/classifiers/functions/supportVector")
## Note that weka.classifiers.functions.supportVector.Kernel provides an
## abstract kernel class, so we could interface this (e.g.,
## make_Weka_kernel()) or at least verify that the class name expansion
## gave something extending the Kernel class.
SMO <-
    make_Weka_classifier("weka/classifiers/functions/SMO",
                         c("SMO", "Weka_functions"),
                         .control_handlers("-K", .expand_kernel_class_name))

### ** Lazy
IBk <-
    make_Weka_classifier("weka/classifiers/lazy/IBk",
                         c("IBk", "Weka_lazy"))
LBR <-
    make_Weka_classifier("weka/classifiers/lazy/LBR",
                         c("LBR", "Weka_lazy"),
                         init = make_Weka_package_loader("lazyBayesianRules"))
                            

### ** Rules
JRip <-
    make_Weka_classifier("weka/classifiers/rules/JRip",
                         c("JRip", "Weka_rules"))
M5Rules <-
    make_Weka_classifier("weka/classifiers/rules/M5Rules",
                         c("M5Rules", "Weka_rules"))
OneR <-
    make_Weka_classifier("weka/classifiers/rules/OneR",
                         c("OneR", "Weka_rules"))
PART <-
    make_Weka_classifier("weka/classifiers/rules/PART",
                         c("PART", "Weka_rules"))

### ** Trees
J48 <-
    make_Weka_classifier("weka/classifiers/trees/J48",
                         c("J48", "Weka_tree"))
M5P <-
    make_Weka_classifier("weka/classifiers/trees/M5P",
                         c("M5P", "Weka_tree"))
LMT <-
    make_Weka_classifier("weka/classifiers/trees/LMT",
                         c("LMT", "Weka_tree"))
DecisionStump <-
    make_Weka_classifier("weka/classifiers/trees/DecisionStump",
                         c("DecisionStump", "Weka_tree"))

### ** Meta learners

.expand_base_learner_spec <-
function(e)
{
    ## For some of the meta learners, option '-W' is used for the full
    ## name of the base classifier, and options after '--' are passed to
    ## this.  Let us be nice and make specifications like
    ##   Weka_control(W = list(J48, M = 50))
    ## work in this case, which means having the control handler convert
    ## to character and insert the '--' if needed.
    if(is.list(e)) {
        init <- get_Java_class(e[[1L]])
        if(length(e) > 1L) {
            rest <- as.character(do.call("Weka_control", e[-1L]))
            if(is.na(match("--", rest)))
                rest <- c("--", rest)
            init <- c(init, rest)
        }
        init
    }
    else
        get_Java_class(e)
}
.Weka_meta_classifier_handlers <-
    .control_handlers("-W", .expand_base_learner_spec)
AdaBoostM1 <-
    make_Weka_classifier("weka/classifiers/meta/AdaBoostM1",
                         c("AdaBoostM1", "Weka_meta"),
                         .Weka_meta_classifier_handlers)
Bagging <-
    make_Weka_classifier("weka/classifiers/meta/Bagging",
                         c("Bagging", "Weka_meta"),
                         .Weka_meta_classifier_handlers)
LogitBoost <-
    make_Weka_classifier("weka/classifiers/meta/LogitBoost",
                         c("LogitBoost", "Weka_meta"),
                         .Weka_meta_classifier_handlers)
MultiBoostAB <-
    make_Weka_classifier("weka/classifiers/meta/MultiBoostAB",
                         c("MultiBoostAB", "Weka_meta"),
                         .Weka_meta_classifier_handlers,
                         init = make_Weka_package_loader("multiBoostAB"))
Stacking <-
    make_Weka_classifier("weka/classifiers/meta/Stacking",
                         c("Stacking", "Weka_meta"),
                         .control_handlers(c("-B", "-M"),
                                           .expand_class_name))

.matrix_to_Matlab_text <-
function(M) 
{
   if (is.character(M))
       M
   else
       sprintf("[%s]",
               paste(apply(as.matrix(M), 1, paste, collapse = " "),
                     collapse = ";"))
}
CostSensitiveClassifier <- 
    make_Weka_classifier("weka/classifiers/meta/CostSensitiveClassifier", 
                         c("CostSensitiveClassifier", "Weka_meta"),
                         .control_handlers(list("-cost-matrix", "-W"),
                                           list(.matrix_to_Matlab_text,
                                                .expand_base_learner_spec)))

### * Clusterers

Cobweb <-
    make_Weka_clusterer("weka/clusterers/Cobweb", "Cobweb")
FarthestFirst <-
    make_Weka_clusterer("weka/clusterers/FarthestFirst", "FarthestFirst")
SimpleKMeans <-
    make_Weka_clusterer("weka/clusterers/SimpleKMeans", "SimpleKMeans")
XMeans <-
    make_Weka_clusterer("weka/clusterers/XMeans", "XMeans",
                        init = make_Weka_package_loader("XMeans"))
DBScan <-
    make_Weka_clusterer("weka/clusterers/DBScan", "DBScan",
                        init = make_Weka_package_loader("optics_dbScan"))

### * Converters

### ** Savers

## <FIXME>
## Shouldn't *all* file saver interfaces use the index decrementer as
## control handler?
## If yes, maybe include this in the make_Weka_file_saver() defaults?
.decrement_number_by_one <-
function(o)
{
    ## Weka_control(c = 2)
    if(is.numeric(o))
        o <- o - 1
    else if(is.character(o)
            && (regexpr("^[[:digit:]]*$", o) > -1L))
        o <- as.numeric(o) - 1
    o
}
.Weka_file_saver_handlers <-
    .control_handlers("-c", .decrement_number_by_one)
C45Saver <-
    make_Weka_file_saver("weka/core/converters/C45Saver",
                         .Weka_file_saver_handlers)
XRFFSaver <-
    make_Weka_file_saver("weka/core/converters/XRFFSaver",
                         .Weka_file_saver_handlers)
## </FIXME>
## Could also provide interfaces to ArffSaver, CSVSaver, and
## LibSVMSaver.

### ** Loaders

C45Loader <-
    make_Weka_file_loader("weka/core/converters/C45Loader")
XRFFLoader <-
    make_Weka_file_loader("weka/core/converters/XRFFLoader")

### * Filters

Normalize <-
    make_Weka_filter("weka/filters/unsupervised/attribute/Normalize", 
                     "Normalize")
Discretize <-
    make_Weka_filter("weka/filters/supervised/attribute/Discretize", 
                     "Discretize")

### * Stemmers

IteratedLovinsStemmer <-
    make_Weka_stemmer("weka/core/stemmers/IteratedLovinsStemmer")
LovinsStemmer <-
    make_Weka_stemmer("weka/core/stemmers/LovinsStemmer")

### * Tokenizers

AlphabeticTokenizer <-
    make_Weka_tokenizer("weka/core/tokenizers/AlphabeticTokenizer")
NGramTokenizer <-
    make_Weka_tokenizer("weka/core/tokenizers/NGramTokenizer")
WordTokenizer <-
    make_Weka_tokenizer("weka/core/tokenizers/WordTokenizer")
    

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
