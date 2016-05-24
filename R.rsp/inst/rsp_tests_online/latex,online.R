library("R.rsp")

if (Sys.getenv("_R_CHECK_FULL_") != "" && isCapableOf(R.rsp, "latex")) {
  urlPath <- "http://latex-project.org/guides"
  filenames <- c(
    usrguide="usrguide.tex"
  );
  urls <- file.path(urlPath, filenames);
  names(urls) <- names(filenames);

  outPath <- file.path("demos", "tex");
  for (kk in seq_along(urls)) {
    url <- urls[kk];
    workdir <- file.path(outPath, names(url));
    output <- rfile(url, workdir=workdir, verbose=TRUE);
    print(output)
  }
}
