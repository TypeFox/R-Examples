path <- dirname(parent.frame(2)$ofile)
for (file in Sys.glob(file.path(path, 'R', '*.R')))
  source(file)
