# Reading in combined code & documentation, and showing help
# Work out where this file lives!

textfilename <- attr( as.environment( 'package:mvbutils'), 'path')
textfilename <- file.path( textfilename, 'demostuff', 'original.dochelp.rrr')

original.dochelp <- source.mvb( textfilename)

cat( "'original.dochelp' has been read in from \"original.dochelp.rrr\" by 'source.mvb'.")
cat( "Press <ENTER> to see 'print( original.dochelp)': ")
readline()
print( original.dochelp)

cat( "Press <ENTER> to see 'write.sourceable.function( original.dochelp, stdout())': ")
readline()
write.sourceable.function( original.dochelp, stdout())

cat( "Press <ENTER> to see the original flat-format help for 'dochelp': ")
readline()

help( original.dochelp)

