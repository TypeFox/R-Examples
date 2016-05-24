`.onLoad` <-
function(libname, pkgname)
{
  .afmpath <<- file.path(libname, pkgname, 'fonts', 'afm')
  .pfbpath <<- file.path(libname, pkgname, 'fonts', 'type1')

  pdfFonts(CMRoman = Type1Font('CMRoman',
             file.path(.afmpath,
                       c('sfrm1000.afm', 'sfbx1000.afm',
                         'sfti1000.afm', 'sfbi1000.afm',
                         'cmsyase.afm'))))
  pdfFonts(CMSans  = Type1Font('CMSans',
             file.path(.afmpath,
                       c('sfss1000.afm', 'sfsx1000.afm',
                         'sfsi1000.afm', 'sfso1000.afm',
                         'cmsyase.afm'))))
}
