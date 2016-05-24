test.for.each <- function() {
  checkEquals(capture.output(for.each(function(a, b)
                                      cat(a + b, ',', sep=''),
                                      list(1, 2),
                                      list(3, 4))),
              "4,6,",
              'list-oriented for.each')
}

test.pair.fold.right <- function() {
  checkEquals(pair.fold.right(function(a, b, accum)
                              append(a, append(b, accum)) ,
                              NULL,
                              list(1, 2),
                              list(3, 4)),
              list(1, 2, 3, 4, 2, 4),
              'pair.fold.right on multiple lists')
}

test.zip.with.names <- function() {
  checkEquals(zip.with.names(j=list(a=1, x=4, h=3),
                             k=list(b=2, y=4, i=10),
                             l=list(c=3, z=5, 11)),
              list(j=list(a=1, b=2, c=3),
                   k=list(x=4, y=4, z=5),
                   l=list(h=3, i=10, 11)),
              'zip.with.names with meta-names, too')
}

test.last <- function() {
  checkEquals(last(list(1, 2, 3)),
              3,
              'last element of trinary list')
}
