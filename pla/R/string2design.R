.string2design <- function(design) {
    design <- tolower(design)
    if ((design == "lsd") | (design == "latin")
        | (design == "latin square")
        | (design == "latin square design")
        | (design == "latin squares design")
        )
        design <- "LSD"
    if ((design == "rbd") | (design == "blocks")
        | (design == "randomised block")
        | (design == "randomised block design")
        | (design == "randomised blocks design")
        | (design == "randomized block")
        | (design == "randomized block design")
        | (design == "randomized blocks design")
        )
        design <- "RBD"
    if ((design == "crd")
        | (design == "completely randomised")
        | (design == "completely randomised design")
        | (design == "completely randomized")
        | (design == "completely randomized design") | (design == "")
        )
        design <- "CRD"
    design <- tolower(design)
    return(design)
}
