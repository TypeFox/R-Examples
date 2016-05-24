# filename does not begin with `test` so not run by `testthat::test_dir()`

# This test assumes the working directory contains all the XML files provided 
# here: https://github.com/rvosa/supertreebase/tree/master/data/treebase

files <- system("ls *.xml", intern=TRUE)

print("testing parsing only")

parses <- sapply(files, function(x){
                out <- try(xmlParse(x))
                if(is(out, "try-error"))
                  out <- x
                else {
                  free(out)
                  out = "success"
                }
                out
                })

fails <- parses[parses!="success"]
works <- files[parses == "success"]

writeLines(fails, "unparseable.txt")
print("testing parsing only")


treebase <- sapply(works, 
                function(x){
                  print(x)
                  tree <- try(nexml_read(x, "nexml"))
                  if(is(tree, "try-error"))
                    out = "read failed:"
                  else {
                    tree <- try(as(tree, "phylo"))
                    if(is(tree, "try-error"))
                      out = "conversion failed:"
                    else 
                      out = "success"
                  }
                  rm(tree)
                  out
                })
save(list=ls(), file = "RNeXML_test_results.rda")

table(treebase)

