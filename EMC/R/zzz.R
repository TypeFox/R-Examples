
### $Id: zzz.R,v 1.10 2008/02/05 20:21:23 goswami Exp $

.onLoad <-
    function (libname, pkgname)
{
    this.year <- substr(as.character(Sys.Date( )), 1, 4)
    packageStartupMessage('##\n',
                          '## Evolutionary Monte Carlo Package (EMC)\n',
                          '##\n',
                          '## Functionality: random walk Metropolis, general Metropolis-Hastings\n',
                          '## parallel tempering, target oriented EMC (TOEMC), temperature ladder\n',
                          '## construction and placement\n',
                          '##\n',
                          '## Use: "help(package = EMC)" at the R prompt for more info\n',
                          '##\n',
                          '## Copyright (C) 2006-', this.year, ' Gopi Goswami\n',
                          '##\n',
                          '##    Created by: Gopi Goswami <goswami@stat.harvard.edu>\n',
                          '## Maintained by: Gopi Goswami <grgoswami@gmail.com>\n',
                          '##\n')        

    library.dynam(pkgname, pkgname, lib.loc = libname)
}






