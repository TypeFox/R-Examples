.slda.collapsed.gibbs.sampler <-
function (documents, K, vocab, num.iterations, alpha, eta, annotations,
    beta, variance, logistic = FALSE, method = "sLDA", lambda,
    initial = NULL, burnin = NULL, trace = 0L)
{
    retval <- structure(.Call("collapsedGibbsSampler", documents,
        as.integer(K), as.integer(length(vocab)), as.integer(num.iterations),
        as.double(alpha), as.double(eta),if (!logistic) as.double(annotations) else if (method=="sLDA" & logistic) as.integer(annotations) else as.logical(annotations),
        as.double(beta), as.double(variance), pmatch(method,
            c("sLDA", "corrLDA", "prodLDA")), as.double(lambda),
        NULL, NULL, initial, burnin, FALSE, trace, FALSE), names = c("assignments",
        "topics", "topic_sums", "document_sums"))
    colnames(retval$topics) <- vocab
    retval
}

.model.filenames <-
function (data.dir, by.time = TRUE, files = c("elbo", "beta",
    "phi"))
{
    stopifnot(by.time)
    all.files <- list.files(data.dir, "elbo-*", full.names = TRUE)
    golden.file <- all.files[which.max(file.info(all.files)$mtime)]
    iteration <- strsplit(golden.file, "-")
    stopifnot(length(iteration) == 1)
    iteration <- as.numeric(iteration[[1]][length(iteration[[1]])])
    c(structure(paste(data.dir, "/", files, "-", iteration, sep = ""),
        names = files), iteration = iteration)
}
.read.beta <-
function (filename, vocab = NULL, num.topics = NULL, ignore.last.row = TRUE)
{
    stopifnot(is.null(num.topics))
    stopifnot(!is.null(vocab))
    result <- matrix(scan(filename, what = 0), byrow = TRUE,
        nrow = length(vocab) + ifelse(ignore.last.row, 1, 0))
    if (ignore.last.row) {
        result <- result[-dim(result)[1], ]
    }
    if (!is.null(vocab)) {
        rownames(result) <- vocab
    }
    result
}
.pairwise.link.lda.collapsed.gibbs.sampler <-
function (documents, K, vocab, num.iterations, alpha, eta, nbeta,
    net.annotations, initial = NULL, burnin = NULL, trace = 0L)
{
    retval <- structure(.Call("collapsedGibbsSampler", documents,
        as.integer(K), as.integer(length(vocab)), as.integer(num.iterations),
        as.double(alpha), as.double(eta), NULL, NULL, NULL, NULL,
        NULL, nbeta, as.logical(net.annotations), initial, burnin,
        FALSE, trace, FALSE), names = c("assignments", "topics",
        "topic_sums", "document_sums", if (is.null(burnin)) NA else "document_expects",
        "net.assignments.left", "net.assignments.right", "blocks.neg",
        "blocks.pos"))
    colnames(retval$topics) <- vocab
    retval
}

.documents.as.Matrix <-
function (docs, vocab)
{
    ii <- rep(1:length(docs), times = sapply(docs, function(x) dim(x)[2]))
    both <- do.call(cbind, docs)
    condensed <- xtabs(both[2, ] ~ ii + both[1, ], sparse = TRUE)
    stopifnot(ncol(condensed) == length(vocab))
    colnames(condensed) <- vocab
    condensed
}
