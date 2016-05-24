#------------------------------------
#-- Skinning variables
#------------------------------------


localOptions <- new.env(parent=globalenv())
localOptions$suffixes <- c(".a", ".p", ".s")
localOptions$style <- "behavior"
localOptions$minVar <- 0

role <- list()
role$behavior <- c("Actor", "Partner", "Relationship")
role$perception <- c("Perceiver", "Target", "Relationship")
role$metaperception <- c("Perceiver", "Target", "Relationship")

# Labels actor-partner
unilabels_b <- c("actor variance", "partner variance", "relationship variance", "error variance", "actor-partner covariance", "relationship covariance")
# Labels target-perceiver
unilabels_p <- c("perceiver variance", "target variance", "relationship variance", "error variance", "perceiver-target covariance", "relationship covariance")

unilabels2 <- c("estimate", "standardized", "se", "t.value")

# labels for metaperception
unilabels_b_meta1 <- c("perceiver variance otherperception", "target variance otherperception",  "relationship variance otherperception", "error variance otherperception", "generalized reciprocity otherperception", "dyadic reciprocity otherperception")
unilabels_b_meta2 <- c("perceiver variance metaperception", "target variance metaperception", "relationship variance metaperception", "error variance metaperception", "generalized reciprocity metaperception", "dyadic reciprocity metaperception")

# Labels for bivariate analyses
bilabels_bb <- c("actor-actor covariance","partner-partner covariance",
"actor-partner covariance","partner-actor covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_pp <- c("perceiver-perceiver covariance","target-target covariance",
"perceiver-target covariance","target-perceiver covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_bp <- c("actor-perceiver covariance","partner-target covariance",
"actor-target covariance","partner-perceiver covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")
bilabels_pb <- c("perceiver-actor covariance","target-partner covariance",
"perceiver-partner covariance","target-actor covariance","intrapersonal relationship covariance", "interpersonal relationship covariance")


bilabels_meta <- c("Perceiver assumed reciprocity","Generalized assumed reciprocity",
"Perceiver meta-accuracy", "Generalized meta-accuracy", "Dyadic assumed reciprocity", "Dyadic meta-accuracy")



# set options for printing results etc.
# style = c("behavior", "perception")
#' @export
RR.style <- function(style="behavior", suffixes=NA, minVar=NA) {
	localOptions$style <- style <- match.arg(style, c("behavior", "perception"))
	
	if (is.na(suffixes)) {
		if (style=="behavior") {
			localOptions$suffixes <- c(".a", ".p", ".s")
		} else 
		if (style=="perception") {
			localOptions$suffixes <- c(".p", ".t", ".s")
		}
	} else {
		localOptions$suffixes <- suffixes
	}
	
	if (is.na(minVar)) {
		localOptions$minVar <- 0
	} else {
		localOptions$minVar <- minVar
	}
}

