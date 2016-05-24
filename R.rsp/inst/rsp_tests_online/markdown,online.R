library("R.rsp")

if (Sys.getenv("_R_CHECK_FULL_") != "" && isCapableOf(R.rsp, "markdown")) {
  type <- "application/x-markdown";
  urlPath <- "https://daringfireball.net/projects/markdown"
  filenames <- c(
    site="basics.text"
  );
  urls <- file.path(urlPath, filenames);
  names(urls) <- names(filenames);

  outPath <- file.path("demos", "markdown");
  for (kk in seq_along(urls)) {
    url <- urls[kk];
    workdir <- file.path(outPath, names(url));
    output <- rfile(url, workdir=workdir, type=type, verbose=TRUE);
    print(output)
  }
}
