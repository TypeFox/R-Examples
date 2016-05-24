
gp_globals = new.env()

gp_globals$version =
	tryCatch(
		installed.packages()["gProfileR", "Version"],
		error = function(e) { return("unknown_version") }
	);
gp_globals$rcurl_opts =
	RCurl::curlOptions(useragent = paste("gProfileR/", gp_globals$version, sep=""))
gp_globals$png_magic =
	as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
gp_globals$base_url =
	"http://biit.cs.ut.ee/gprofiler/"

#' Annotate gene list functionally.
#'
#' Interface to the g:Profiler tool for finding enrichments in gene lists. Organism
#' names are constructed by concatenating the first letter of the name and the
#' family name. Example: human - 'hsapiens', mouse - 'mmusculus'. If requesting PNG
#' output, the request is directed to the g:GOSt tool in case 'query' is a vector
#' and the g:Cocoa (compact view of multiple queries) tool in case 'query' is a
#' list. PNG output can fail (return FALSE) in case the input query is too large.
#' In such case, it is advisable to fall back to a non-image request.
#'
#' @param organism organism name.
#' @param query vector of gene IDs or a list of such vectors.
#' @param ordered_query in case output gene lists are ranked this option may be
#'  used to get GSEA style p-values.
#' @param significant whether all or only statistically significant results should
#'  be returned.
#' @param exclude_iea exclude electronic annotations (IEA).
#' @param underrep measure underrepresentation.
#' @param evcodes include GO evidence codes as the final column of output. Note 
#'  that this can decrease performance and make the query slower.
#' @param region_query interpret query as chromosomal ranges.
#' @param max_p_value custom p-value threshold, results with a larger p-value are
#'  excluded.
#' @param min_set_size minimum size of functional category, smaller categories are
#'  excluded.
#' @param max_set_size maximum size of functional category, larger categories are
#'  excluded.
#' @param min_isect_size minimum size of the overlap (intersection) between query 
#'  and functional category, smaller intersections are excluded.
#' @param correction_method the algorithm used for determining the significance
#'  threshold, one of "gSCS", "fdr", "bonferroni".
#' @param hier_filtering hierarchical filtering strength, one of "none",
#' "moderate", "strong".
#' @param domain_size statistical domain size, one of "annotated", "known".
#' @param custom_bg vector of gene names to use as a statistical background.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param png_fn request the result as PNG image and write it to png_fn.
#' @param include_graph request inclusion of network data with the result.
#' @param src_filter a vector of data sources to use. Currently, these include
#'  GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF,
#'  MI, CORUM, HP, HPA, OMIM. Please see the g:GOSt web tool for the comprehensive 
#'  list and details on incorporated data sources.
#' @return A data frame with the enrichment analysis results. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query number'.  When requesting a PNG image, either TRUE or FALSE, depending on
#'  whether a non-empty result was received and a file written or not, respectively.
#'  If 'include_graph' is set, the return value may include the attribute 'networks',
#'  containing a list of all network sources, each in turn containing a list of graph edges.
#'  The edge structure is a list containing the two interacting symbols and two boolean
#'  values (in that order), indicating whether the first or second interactor is part of
#'  the input query (core nodes).
#' @references  J. Reimand, M. Kull, H. Peterson, J. Hansen, J. Vilo: g:Profiler -
#'  a web-based toolset for functional profiling of gene lists from large-scale
#'  experiments (2007) NAR 35 W193-W200
#' @author  Juri Reimand <jyri.reimand@@ut.ee>, Raivo Kolde <rkolde@@gmail.com>,
#'  Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#'  gprofiler(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")
#' @export

gprofiler <- function(
	query,
	organism = "hsapiens",
	ordered_query = F,
	significant = T,
	exclude_iea = F,
	underrep = F,
	evcodes = F,
	region_query = F,
	max_p_value = 1.0,
	min_set_size = 0,
	max_set_size = 0,
	min_isect_size = 0,
	correction_method = "analytical",
	hier_filtering = "none",
	domain_size = "annotated",
	custom_bg = "",
	numeric_ns = "",
	png_fn = NULL,
	include_graph = F,
	src_filter = NULL
) {
	query_url = ""
	my_url = paste(gp_globals$base_url, "gcocoa.cgi", sep="")
	wantpng = ifelse(is.character(png_fn), T, F)
	output_type = ifelse(wantpng, "mini_png", "mini")

	# Query

	if (is.list(query)) {
		qnames = names(query)

		for (i in 1:length(query)) {
			cname = paste("> query", i)
			if (!is.null(qnames) && nchar(qnames[i]) > 0)
				cname = paste(">", qnames[i])
			query_url <- paste(sep="\n", query_url, cname, paste(query[[i]], collapse=" "))
		}
	} else if (is.vector(query)) {
		query_url <- paste(query, collapse=" ")
		if (wantpng)
			my_url = paste(gp_globals$base_url, "index.cgi", sep="")
	} else {
		warning("Missing query")
		return()
	}

	# Significance threshold

	if (correction_method == "gSCS")
		correction_method = "analytical"
	if (!correction_method %in% c("analytical", "fdr", "bonferroni"))
		stop("Multiple testing correction method not recognized")

	# Hierarchical filtering

	if (!hier_filtering %in% c("none", "moderate", "strong"))
		stop("hier_filtering must be one of \"none\", \"moderate\" or \"strong\"")
	if (hier_filtering == "strong") hier_filtering = "compact_ccomp"
	else if (hier_filtering == "moderate") hier_filtering = "compact_rgroups"
	else hier_filtering = ""

	# Domain size

	if (!domain_size %in% c("annotated", "known"))
		stop("domain_size must be one of \"annotated\" or \"known\"")

	# Custom background

	if (is.vector(custom_bg))
		custom_bg = paste(custom_bg, collapse=" ")
	else
		stop("custom_bg must be a vector")

	# Max. set size

	if (max_set_size < 0) max_set_size = 0

	# HTTP request

	query_params = list(
		organism = organism,
		query = query_url,
		output = output_type,
		analytical = "1",
		sort_by_structure = "1",
		ordered_query = ifelse(ordered_query, "1", "0"),
		significant = ifelse(significant, "1", "0"),
		no_iea = ifelse(exclude_iea, "1", "0"),
		underrep = ifelse(underrep, "1", "0"),
		txt_evcodes = ifelse(evcodes, "1", "0"),
		as_ranges = ifelse(region_query, "1", "0"),
		omit_metadata = ifelse(include_graph, "0", "1"),
		user_thr = as.character(max_p_value),
		min_set_size = as.character(min_set_size),
		max_set_size = as.character(max_set_size),
		min_isect_size = as.character(min_isect_size),
		threshold_algo = correction_method,
		hierfiltering = hier_filtering,
		domain_size_type = domain_size,
		custbg_file = "",
		custbg = custom_bg,
		prefix = numeric_ns
	)

	if (!is.null(src_filter)) {
		for (i in as.vector(src_filter))
			query_params[paste("sf_", i, sep="")] = "1"
	}

	raw_query <- RCurl::postForm(
		my_url,
		.opts = gp_globals$rcurl_opts,
		.params = query_params
	)

	# Requested PNG, write to disk and return

	if (wantpng) {
		pngv = as.vector(raw_query)

		if (length(raw_query) == 0)
			return(F)
		if (!isTRUE(all.equal(utils::head(pngv, 8), gp_globals$png_magic))) { # check PNG magic number
			warning(paste(
				"Received non-PNG data, but PNG output requested. Please report this error",
				"with the appropriate details to the package maintainer"
			));
			return(F);
		}

		conn = file(png_fn, "wb")
		writeBin(pngv, conn, useBytes = T)
		close(conn)
		return(T);
	}

	# Requested text

	split_query <- unlist(strsplit(raw_query, split="\n"))
	commented_lines <- grep("^#", split_query)

	# Interactions, later saved in an attribute

	if (include_graph) {
		edges <- grep("BIOGRID INTERACTION", split_query[commented_lines], value=T)
		re <- "\\S+\\s+\\S+$"
		m <- regexpr(re, edges, perl=T);
		edges <- strsplit(substr(edges, m, m+attr(m, "match.length")), '\\s+')

		is_in_core <- function(x) {
			incore <- c(F, F);
			incore[grep("\\*$", x, perl=T)] = T
			x <- sub("\\*$", "", x, perl=T)
			x <- list(x[[1]], x[[2]], incore[[1]], incore[[2]])
			return(x)
		}

		edges <- lapply(edges, is_in_core)
		edges <- list(edges)
		names(edges) <- c("BIOGRID")
	}

	# Parse main result body

	if (length(commented_lines)>0) {
		split_query <- split_query[-commented_lines]
	}

	empty_lines <- which(split_query == "")
	if (length(empty_lines)>0) {
		split_query <- split_query[-empty_lines]
	}

	if(length(split_query) > 0) {
		conn <- textConnection(paste(split_query, collapse = "\n"))
		split_query <- utils::read.table(conn, sep = "\t", quote = "", stringsAsFactors = F)
		close(conn)
	}
	else {
		split_query <- as.data.frame(matrix(NA, 0, 14))
	}
	
	col_names <- c(
		"query.number", "significant", "p.value",
		"term.size", "query.size", "overlap.size",
		"recall", "precision", "term.id",
		"domain", "subgraph.number", "term.name",
		"relative.depth", "intersection"
	)
	if (evcodes)
		col_names <- append(col_names, "evidence.codes")
	
	rownames(split_query) <- NULL
	colnames(split_query) <- col_names
	split_query$term.name <- gsub("^\\s+", "", split_query$term.name)
	split_query$significant <- ifelse(split_query$significant == "!", T, F)

	if(is.list(query) & !is.null(names(query))) {
		split_query$query.number <- names(query)[split_query$query.number]
	}

	if (include_graph && length(edges)) {
		attr(split_query, "networks") <- edges
	}

	return(split_query)
}

#' Convert gene IDs.
#'
#' Interface to the g:Convert tool. Organism names are constructed by
#' concatenating the first letter of the name and the family name. Example: human
#' - 'hsapiens', mouse - 'mmusculus'.
#'
#' @param query list of gene IDs.
#' @param organism organism name.
#' @param target target namespace.
#' @param region_query interpret query as chromosomal ranges.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param mthreshold maximum number of results per initial alias to show.
#' @param filter_na logical indicating whether to filter out results without a
#' corresponding target.
#' @param df logical indicating whether the output will be a data.frame or list.
#' @return The output can be either a list or a data.frame. The list has an entry
#' for every input gene. The data frame is a table closely corresponding to the
#' web interface output.
#' @references  J. Reimand, M. Kull, H. Peterson, J.  Hansen, J. Vilo: g:Profiler
#' - a web-based toolset for functional profiling of gene lists from large-scale
#' experiments (2007) NAR 35 W193-W200
#' @author  Juri Reimand <jyri.reimand@@ut.ee>, Raivo Kolde <rkolde@@gmail.com>,
#' Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#' gconvert(c("POU5F1", "SOX2", "NANOG"), organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")
#' @export

gconvert = function(
	query,
	organism = "hsapiens",
	target = "ENSG",
	region_query = F,
	numeric_ns = "",
	mthreshold = Inf,
	filter_na = T,
	df = T
) {
	url = paste(gp_globals$base_url, "gconvert.cgi", sep="")

	raw_query <- RCurl::postForm(url, .opts = gp_globals$rcurl_opts,
		organism=organism,
		output='mini',
		query=paste(query, collapse="+"),
		target=target,
		region_query = ifelse(region_query, "1", "0"),
		prefix = numeric_ns
	)

	conn <- textConnection(raw_query)
	tab <- utils::read.table(conn, sep = "\t", quote = "")
	close(conn)

	colnames(tab) <- c(
		"alias.number", "alias", "target.number",
		"target", "name", "description", "namespace"
	)

	res = plyr::dlply(tab, "alias", function(x) {
		if (filter_na)
			x = x[x$target != "N/A",]
		if (nrow(x) == 0)
			return(NULL)
		if (nrow(x) > mthreshold)
			x = lapply(x, function(y) { as.character(y)[1:mthreshold] })
		return(x)
	})

	if(df)
		res <- plyr::ldply(res, function(x) data.frame(x, stringsAsFactors = F))
	return(res)
}

#' Find orthologs.
#'
#' Interface to the g:Orth tool. Organism names are constructed by concatenating
#' the first letter of the name and the family name. Example: human - 'hsapiens',
#' mouse - 'mmusculus'.
#'
#' To alleviate the problem of having many orthologs per gene (most of them
#' uninformative) one can set a threshold for the number of results. The program
#' tries to find the most informative by selecting the most popular ones.
#'
#' @param query list of gene IDs to be translated.
#' @param source_organism name of the source organism.
#' @param target_organism name of the target organism.
#' @param region_query interpret query as chromosomal ranges.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param mthreshold maximum number of ortholog names per gene to show.
#' @param filter_na logical indicating whether to filter out results without a
#' corresponding target name.
#' @param df logical indicating whether the output will be a data.frame or list.
#' @return The output can be either a list or a data.frame. The list has an entry
#' for every input gene. The data frame is a table closely corresponding to the
#' web interface output.
#' @references  J. Reimand, M. Kull, H. Peterson, J.  Hansen, J. Vilo: g:Profiler
#' -- a web-based toolset for functional profiling of gene lists from large-scale
#' experiments (2007) NAR 35 W193-W200
#' @author  Raivo Kolde <rkolde@@gmail.com>, Juri Reimand <juri.reimand@@ut.ee>,
#' Tambet Arak <tambet.arak@@gmail.com>
#' @examples
#' gorth(c("Klf4","Pax5","Sox2","Nanog"), source_organism="mmusculus", target_organism="hsapiens")
#' @export

gorth <- function(
	query,
	source_organism = "hsapiens",
	target_organism = "mmusculus",
	region_query = F,
	numeric_ns = "",
	mthreshold = Inf,
	filter_na = T,
	df = T
){
	my_url = paste(gp_globals$base_url, "gorth.cgi", sep="")

	if (length(query) == 0)
		return(NULL)

	raw_query <- RCurl::postForm(my_url, .opts = gp_globals$rcurl_opts,
		output = 'mini',
		query = paste(query, collapse = " "),
		organism = source_organism,
		target = target_organism,
		region_query = ifelse(region_query, "1", "0"),
		prefix = numeric_ns
	)

	conn <- textConnection(raw_query)
	tab <- utils::read.table(conn, sep = "\t", quote = "")
	close(conn)

	colnames(tab) <- c(
		"alias.number", "initial.alias", "initial.ensg",
		"ensg.number", "target.ensg", "target.name",
		"target.description"
	)

	res = plyr::dlply(tab, "initial.alias", function(x) {
		if (filter_na)
			x = x[x$target.name != "N/A",]
		if (nrow(x) == 0)
			return(NULL)
		if (nrow(x) > mthreshold)
			x = lapply(x, function(y) { as.character(y)[1:mthreshold] })
		return(x)
	})

	if(df) {
		res <- plyr::ldply(res, function(x) data.frame(x, stringsAsFactors = F))
		res$initial.alias = as.character(res$initial.alias)
	}

	return(res)
}

#' Get current user agent string.
#'
#' Get the HTTP User-Agent string.
#'
#' @export

get_user_agent = function() {
	gp_globals$rcurl_opts$useragent
}

#' Set custom user agent string.
#'
#' Set the HTTP User-Agent string. Useful for overriding the default user agent for
#' packages that depend on gProfileR functionality.
#'
#' @param ua the user agent string.
#' @param append logical indicating whether to append the passed string to the default user agent string.
#' @export

set_user_agent = function(ua, append = F) {
	rco = gp_globals$rcurl_opts
	rco$useragent = ifelse(
		append, paste(rco$useragent, ua, sep=""), ua)
	assign("rcurl_opts", rco, envir=gp_globals)
}

#' Get the base URL.
#'
#' @export

get_base_url = function() {
	gp_globals$base_url
}

#' Set the base URL.
#'
#' Set the base URL. Useful for overriding the default URL (http://biit.cs.ut.ee/gprofiler)
#' with the bleeding-edge beta or an archived version.
#'
#' @param url the base URL.
#' @export

set_base_url = function(url) {
	url = as.character(url)
	schema = substr(url, 1, 7)
	suffix = substr(url, nchar(url), nchar(url))

	if (schema != "http://")
		stop("The schema portion of the URL must be \"http://\"")
	if (suffix != "/")
		url = paste(url, "/", sep="")

	assign("base_url", url, envir=gp_globals)
}
