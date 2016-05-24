# binary formats cannot be compared line-by-line
binary.formats <- c("hdf5", "rds");

# use print-output-based object comparison tests for noncomparable objects
noncomparable.formats <- c("xml");

# exempt some of the externally implemented formats from connection tests
conn.exempt.formats <- union(c("xml", "json"), binary.formats);

suggested.packages <- list(xml="XML", hdf5="rhdf5", yaml="yaml", lst="yaml", json="jsonlite");
suggested.package.versions <- list(
	xml = list(name = "XML", op = ">=", version = package_version("3.98-1.1")),
	hdf5 = list(name = "rhdf5", op = ">=", version = package_version("2.10.0")),
	yaml = list(name = "yaml", op = ">=", version = package_version("2.1.13")),
	jsonlite = list(name = "jsonlite", op = ">=", version = package_version("0.9.14"))
);

.requireNamespace <- function(x, versionCheck) {
	x <- try(loadNamespace(x, versionCheck=versionCheck), silent=TRUE);
	if (class(x) == "try-error") {
		FALSE
	} else {
		TRUE
	}
}

# read an input file, write it to a temporary file
# and test that the files are the same;
# write a data object to file, read it back in
# and test that the data objects are the same
test_read_write_read <- function(infile) {

	ext <- tolower(as.filename(infile)$ext);

	if (ext %in% names(suggested.packages)) {
		if (!.requireNamespace(
			suggested.packages[[ext]], suggested.package.versions[[ext]]
		)) {
			# optional depedency is missing: skip test
			return(invisible());
		}
	}

	outfile <- tempfile("test-out", fileext=paste(".", ext, sep=""));

	# read-write-read cycle
	x <- qread(infile);
	qwrite(x, outfile);
	y <- qread(outfile);

	if (! ext %in% binary.formats) {

		inlines <- readLines(infile);

		test_that("read-and-write does not change file", {
			outlines <- readLines(outfile);
			expect_equal(inlines, outlines);
		});

	}

	if (ext %in% noncomparable.formats) {

		test_that("write-and-read does not change data (printed)", {
			expect_equal(capture.output(print(x)), capture.output(print(y)));
		});

	} else {

		test_that("write-and-read does not change data", {
			expect_equal(x, y);
		});

	}

	if (! ext %in% conn.exempt.formats) {

		test_that("read from text connection does not change content", {
			conn <- textConnection(inlines, "rt");
			z <- qread(conn, type=ext);
			close(conn);
			expect_equal(x, z);
		});

		test_that("write to text connection does not change content", {
			textlines <- character();
			conn <- textConnection("textlines", "wt", local=TRUE);
			# type has to be specified since c2 is a connection (no file extension!)
			qwrite(x, conn, type=ext);
			close(conn);
			expect_equal(inlines, textlines);
		});

	}

	file.remove(outfile);
}
