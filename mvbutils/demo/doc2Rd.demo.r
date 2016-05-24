# Show example of flat-format documentation, and conversion to Rd format
# Work out where this file lives!

textfilename <- attr( as.environment( 'package:mvbutils'), 'path')
textfilename <- file.path( textfilename, 'demostuff', 'sample.fun.rrr')

sample.fun <- source.mvb( textfilename)

cat( "'sample.fun' has been read in from \"sample.fun.rrr\" by 'source.mvb'.")
cat( "Press <ENTER> to see its informal help: ")
readline()
help( sample.fun)

cat( "Press <ENTER> to see 'cat( doc2Rd( sample.fun), sep='\n')': ")
readline()
cat( doc2Rd( sample.fun), sep='\n')

