library("R.rsp")
library("R.utils") # Arguments

path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.tex")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  outPath <- file.path("LoremIpsum", "tex");
  pathnameR <- compileLaTeX(pathname, outPath=outPath, verbose=-10)
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}


# Fake LaTeX compilation
if (Sys.getenv("_R_CHECK_FULL_") != "") {
  outPath <- file.path("LoremIpsum", "tex-fake");
  local({
    oenv <- Sys.getenv("R_RSP_COMPILELATEX_FALLBACK")
    on.exit(Sys.setenv("R_RSP_COMPILELATEX_FALLBACK"=oenv))
    Sys.setenv("R_RSP_COMPILELATEX_FALLBACK"="copy-force")
    pathnameR <- compileLaTeX(pathname, outPath=outPath, verbose=-10)
  })
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}

