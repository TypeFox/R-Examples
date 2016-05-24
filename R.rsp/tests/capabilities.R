library("R.rsp")

cat("Tools supported by the package:\n")
print(capabilitiesOf(R.rsp))

# Check capabilities one by one
for (what in c("asciidoc", "knitr", "markdown", "latex", "pandoc", "sweave", "unknown")) {
  cat(sprintf(" - %s: %s\n", what, isCapableOf(R.rsp, what)))
}
