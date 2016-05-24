### <NOTE>
### At least currently, all registries are meant and used for all types
### of clusterings (for the time being, partitions and hierarchies)
### simultaneously.
### </NOTE>

### * Internal stuff.

.make_db_key <- function(name, type) paste(type, name, sep = "_")

### * General-purpose stuff.

### <FIXME>
### This currently insists on a given type: maybe it should simply list
### everything split according to type.  But hey, it's internal stuff
### anyway (at least for the time being ...)
### </FIXME>
get_methods_from_db <- function(db, type) {
    type <- match.arg(type, c("partition", "hierarchy"))
    
    pattern <- sprintf("^%s_", type)
    sub(pattern, "",
        grep(pattern, objects(db), value = TRUE))
}

get_method_from_db <- function(db, type, name, msg) {
    ## <NOTE>
    ## Keep 'msg' here so that gettext()ing could work ...
    ## </NOTE>
    type <- match.arg(type, c("partition", "hierarchy"))

    db_keys <- objects(db)
    ind <- pmatch(.make_db_key(tolower(name), type), tolower(db_keys))
    if(is.na(ind))
        stop(msg, call. = FALSE, domain = NA)
    
    db[[db_keys[ind]]]
}

put_method_into_db <- function(db, type, name, value) {
    type <- match.arg(type, c("partition", "hierarchy"))
    db[[.make_db_key(name, type)]] <- value
}

### * Consensus Method Registration.

cl_consensus_methods_db <- new.env()

get_cl_consensus_methods <- function(type)
    get_methods_from_db(cl_consensus_methods_db, type)
get_cl_consensus_method <- function(name, type) {
    get_method_from_db(cl_consensus_methods_db, type, name,
                       gettextf("Invalid consensus method '%s'.", name))
}
                       
set_cl_consensus_method <-
function(name, type, definition, ...)
{
    ## Register a @code{type} consensus method called @code{name} with
    ## definition @code{definition}.  Provide more information where
    ## appropriate, e.g., @code{dissimilarity} d and @code{exponent} e
    ## for methods minimizing \sum_b d(x_b, x) ^ e.

    put_method_into_db(cl_consensus_methods_db, type, name,
                       .structure(c(list(definition = definition),
                                    list(...)),
                                  class = "cl_consensus_method"))
}

set_cl_consensus_method("DWH", "partition",
                        .cl_consensus_partition_DWH,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("soft/euclidean", "partition",
                        .cl_consensus_partition_soft_euclidean,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("SE", "partition",
                        .cl_consensus_partition_soft_euclidean,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("hard/euclidean", "partition",
                        .cl_consensus_partition_hard_euclidean,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("HE", "partition",
                        .cl_consensus_partition_hard_euclidean,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("soft/manhattan", "partition",
                        .cl_consensus_partition_soft_manhattan,
                        dissimilarity = "manhattan",
                        exponent = 1)
set_cl_consensus_method("SM", "partition",
                        .cl_consensus_partition_soft_manhattan,
                        dissimilarity = "manhattan",
                        exponent = 1)
set_cl_consensus_method("hard/manhattan", "partition",
                        .cl_consensus_partition_hard_manhattan,
                        dissimilarity = "manhattan",
                        exponent = 1)
set_cl_consensus_method("HM", "partition",
                        .cl_consensus_partition_hard_manhattan,
                        dissimilarity = "manhattan",
                        exponent = 1)
set_cl_consensus_method("GV1", "partition",
                        .cl_consensus_partition_GV1,
                        dissimilarity = "GV1",
                        exponent = 2)
set_cl_consensus_method("GV3", "partition",
                        .cl_consensus_partition_GV3,
                        dissimilarity = "comemberships",
                        exponent = 2)
set_cl_consensus_method("soft/symdiff", "partition",
                        .cl_consensus_partition_soft_symdiff,
                        dissimilarity = "symdiff",
                        exponent = 1)
set_cl_consensus_method("hard/symdiff", "partition",
                        .cl_consensus_partition_hard_symdiff,
                        dissimilarity = "symdiff",
                        exponent = 1)

set_cl_consensus_method("cophenetic", "hierarchy",
                        .cl_consensus_hierarchy_cophenetic,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("euclidean", "hierarchy",
                        .cl_consensus_hierarchy_cophenetic,
                        dissimilarity = "euclidean",
                        exponent = 2)
set_cl_consensus_method("manhattan", "hierarchy",
                        .cl_consensus_hierarchy_manhattan,
                        dissimilarity = "manhattan",
                        exponent = 1)
set_cl_consensus_method("majority", "hierarchy",
                        .cl_consensus_hierarchy_majority,
                        dissimilarity = "symdiff",
                        exponent = 1)

### * Dissimilarity Method Registration.

cl_dissimilarity_methods_db <- new.env()

get_cl_dissimilarity_methods <- function(type)
    get_methods_from_db(cl_dissimilarity_methods_db, type)
get_cl_dissimilarity_method <- function(name, type)
    get_method_from_db(cl_dissimilarity_methods_db, type, name,
                       gettextf("Invalid dissimilarity method '%s'.",
                                name))
set_cl_dissimilarity_method <-
    function(name, type, definition, description, ...)
    put_method_into_db(cl_dissimilarity_methods_db, type, name,
                       .structure(c(list(definition = definition,
                                         description = description),
                                    list(...)),
                                  class = "cl_dissimilarity_method"))

set_cl_dissimilarity_method("euclidean", "partition",
                            .cl_dissimilarity_partition_euclidean,
                            "minimal Euclidean membership distance")
set_cl_dissimilarity_method("manhattan", "partition",
                            .cl_dissimilarity_partition_manhattan,
                            "minimal Manhattan membership distance")
set_cl_dissimilarity_method("comemberships", "partition",
                            .cl_dissimilarity_partition_comemberships,
                            "Euclidean comembership distance")
set_cl_dissimilarity_method("symdiff", "partition",
                            .cl_dissimilarity_partition_symdiff,
                            "symmetric difference distance")
set_cl_dissimilarity_method("Rand", "partition",
                            .cl_dissimilarity_partition_Rand,
                            "Rand distance")
set_cl_dissimilarity_method("GV1", "partition",
                            .cl_dissimilarity_partition_GV1,
                            "Gordon-Vichi Delta_1 dissimilarity")
set_cl_dissimilarity_method("BA/A", "partition",
                            .cl_dissimilarity_partition_BA_A,
                            "Boorman/Arabie minimum element moves distance")
set_cl_dissimilarity_method("BA/C", "partition",
                            .cl_dissimilarity_partition_BA_C,
                            "Boorman/Arabie minimum lattice moves distance")
set_cl_dissimilarity_method("BA/D", "partition",
                            .cl_dissimilarity_partition_BA_D,
                            "Boorman/Arabie pair-bonds distance")
set_cl_dissimilarity_method("BA/E", "partition",
                            .cl_dissimilarity_partition_BA_E,
                            "Boorman/Arabie normalized information distance")
set_cl_dissimilarity_method("VI", "partition",
                            .cl_dissimilarity_partition_VI,
                            "Variation of information")
set_cl_dissimilarity_method("Mallows", "partition",
                            .cl_dissimilarity_partition_Mallows,
                            "Mallows dissimilarity")
set_cl_dissimilarity_method("CSSD", "partition",
                            .cl_dissimilarity_partition_CSSD,
                            "Cluster Similarity Sensitive Distance")

set_cl_dissimilarity_method("euclidean", "hierarchy",
                            .cl_dissimilarity_hierarchy_euclidean,
                            "Euclidean ultrametric distance")
set_cl_dissimilarity_method("manhattan", "hierarchy",
                            .cl_dissimilarity_hierarchy_manhattan,
                            "Manhattan ultrametric distance")
set_cl_dissimilarity_method("cophenetic", "hierarchy",
                            .cl_dissimilarity_hierarchy_cophenetic,
                            "cophenetic correlations")
set_cl_dissimilarity_method("gamma", "hierarchy",
                            .cl_dissimilarity_hierarchy_gamma,
                            "rate of inversions")
set_cl_dissimilarity_method("symdiff", "hierarchy",
                            .cl_dissimilarity_hierarchy_symdiff,
                            "symmetric difference distance")
set_cl_dissimilarity_method("Chebyshev", "hierarchy",
                            .cl_dissimilarity_hierarchy_Chebyshev,
                            "Chebyshev distance")
set_cl_dissimilarity_method("Lyapunov", "hierarchy",
                            .cl_dissimilarity_hierarchy_Lyapunov,
                            "Lyapunov distance")
set_cl_dissimilarity_method("BO", "hierarchy",
                            .cl_dissimilarity_hierarchy_BO,
                            "Boorman/Olivier m_delta tree distance")
set_cl_dissimilarity_method("spectral", "hierarchy",
                            .cl_dissimilarity_hierarchy_spectral,
                            "spectral ultrametric distance")

### * Agreement Method Registration.

cl_agreement_methods_db <- new.env()
get_cl_agreement_methods <-
function(type)
    get_methods_from_db(cl_agreement_methods_db, type)
get_cl_agreement_method <-
function(name, type)
    get_method_from_db(cl_agreement_methods_db, type, name,
                       gettextf("Invalid agreement method '%s'.",
                                name))
set_cl_agreement_method <-
function(name, type, definition, description, ...)
    put_method_into_db(cl_agreement_methods_db, type, name,
                       .structure(c(list(definition = definition,
                                         description = description),
                                    list(...)),
                                  class = "cl_agreement_method"))

set_cl_agreement_method("euclidean", "partition",
                        .cl_agreement_partition_euclidean,
                        "minimal euclidean membership distance")
set_cl_agreement_method("manhattan", "partition",
                        .cl_agreement_partition_manhattan,
                        "minimal manhattan membership distance")
set_cl_agreement_method("Rand", "partition",
                        .cl_agreement_partition_Rand,
                        "Rand index")
set_cl_agreement_method("cRand", "partition",
                        .cl_agreement_partition_cRand,
                        "corrected Rand index")
set_cl_agreement_method("NMI", "partition",
                        .cl_agreement_partition_NMI,
                        "normalized mutual information")
set_cl_agreement_method("KP", "partition",
                        .cl_agreement_partition_KP,
                        "Katz-Powell index")
set_cl_agreement_method("angle", "partition",
                        .cl_agreement_partition_angle,
                        "maximal angle between memberships")
set_cl_agreement_method("diag", "partition",
                        .cl_agreement_partition_diag,
                        "maximal co-classification rate")
set_cl_agreement_method("FM", "partition",
                        .cl_agreement_partition_FM,
                        "Fowlkes-Mallows index")
set_cl_agreement_method("Jaccard", "partition",
                        .cl_agreement_partition_Jaccard,
                        "Jaccard index")
set_cl_agreement_method("purity", "partition",
                        .cl_agreement_partition_purity,
                        "purity")
set_cl_agreement_method("PS", "partition",
                        .cl_agreement_partition_PS,
                        "Prediction Strength")

set_cl_agreement_method("euclidean", "hierarchy",
                        .cl_agreement_hierarchy_euclidean,
                        "euclidean ultrametric distance")
set_cl_agreement_method("manhattan", "hierarchy",
                        .cl_agreement_hierarchy_manhattan,
                        "manhattan ultrametric distance")
set_cl_agreement_method("cophenetic", "hierarchy",
                        .cl_agreement_hierarchy_cophenetic,
                        "cophenetic correlations")
set_cl_agreement_method("angle", "hierarchy",
                        .cl_agreement_hierarchy_angle,
                        "angle between ultrametrics")
set_cl_agreement_method("gamma", "hierarchy",
                        .cl_agreement_hierarchy_gamma,
                        "rate of inversions")
