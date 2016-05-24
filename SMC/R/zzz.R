
.onLoad <-
    function (libname, pkgname)
{
    this.year <- substr(as.character(Sys.Date( )), 1, 4)    
    packageStartupMessage('##\n',
                          '## Sequential Monte Carlo Package (SMC)\n',
                          '##\n',
                          '## Functionality: sequential Monte Carlo (SMC) or sequential importance\n',
                          '## sampling (SIS) or hidden Markov models (HMM), particle filter (PF)\n',
                          '## and auxiliary particle filter (APF)\n',
                          '##\n',
                          '## Use: "help(package = SMC)" at the R prompt for more info\n',
                          '##\n',
                          '## Copyright (C) 2006-', this.year, ' Gopi Goswami\n',
                          '##\n',
                          '##    Created by: Gopi Goswami <goswami@stat.harvard.edu>\n',
                          '## Maintained by: Gopi Goswami <grgoswami@gmail.com>\n',
                          '##\n')
    
    library.dynam(pkgname, pkgname, lib.loc=libname)
}






