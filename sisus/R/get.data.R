get.data <-
function# Get data from Excel spreadsheet
### internal function for sisus
(filename
### internal variable
)
{
  ##details<<
  ## interal function for sisus.run()

#
#   When making changes to the workbook, need to:
#     1. get.data.R define sheet names
#     2. assign.variables.R define variables
#     3. sisus.R assign variables
#     4. write.input.R define sheet names
#
  # read data file and convert to a matrix
  # attached via DESCRIPTION # library("gdata"); # to read Excel files

  #Perl.path = "C:\\Perl\\bin\\perl.exe";
  Perl.path = "perl";

  n.sheets = 10;  # version v0.09 has 10 sheets
      p.o = paste("   reading up to ", n.sheets, "sheets, as requested: "); write.out(p.o);

  DATA = new.env();  # create a list to return with all data

  # for each sheet
  for (i.sheet in seq(1, n.sheets)) {
      write.out(paste(" ", i.sheet));

         #DATA         = read.csv(filename, header=FALSE, sep=",", quote="\"", dec=".", fill = TRUE, comment.char="");

    # assign each worksheet to a matrix in the DATA list

    if (i.sheet ==  1){
                                                                         write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );                                                                              DATA$SISUS.PARAMETERS                = DATA.TEMP;
      SW = assign.variables(DATA, 1);  # only read the sheets we have to based on the switches in the first worksheet
    };
    if (i.sheet ==  2){                                                  write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );                                                                              DATA$MIXTURES.ISOTOPE.RATIOS         = DATA.TEMP;};
    if (i.sheet ==  3){                                                  write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );                                                                              DATA$SOURCES.ISOTOPE.RATIOS          = DATA.TEMP;};
    if (i.sheet ==  4){ if (SW$discrimination.include.sw         == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.DISCRIMINATION         not read");};  DATA$SOURCES.DISCRIMINATION          = DATA.TEMP;};
    if (i.sheet ==  5){ if (SW$concentration.include.sw          == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.CONCENTRATION          not read");};  DATA$SOURCES.CONCENTRATION           = DATA.TEMP;};
    if (i.sheet ==  6){ if (SW$assimeffic.include.sw             == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.EFFICIENCY             not read");};  DATA$SOURCES.EFFICIENCY              = DATA.TEMP;};
    if (i.sheet ==  7){ if (SW$lc.include.sw                     == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("LINEAR.CONSTRAINTS             not read");};  DATA$LINEAR.CONSTRAINTS              = DATA.TEMP;};
    # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
    #if (i.sheet ==  8){ if (SW$prior.include.sw                  == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("PRIOR.DIRICHLET                not read");};  DATA$PRIOR.DIRICHLET                 = DATA.TEMP;};
    if (i.sheet ==  8){ if (SW$isotope.stddev.include.sw         == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("MIXTURES.ISOTOPE.STDDEV        not read");};  DATA$MIXTURES.ISOTOPE.STDDEV         = DATA.TEMP;};
    if (i.sheet ==  9){ if (SW$isotope.stddev.include.sw         == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.ISOTOPE.STDDEV         not read");};  DATA$SOURCES.ISOTOPE.STDDEV          = DATA.TEMP;};
    if (i.sheet == 10){ if (SW$biomass.per.individual.include.sw == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.BIOMASS.PER.INDIVIDUAL not read");};  DATA$SOURCES.BIOMASS.PER.INDIVIDUAL  = DATA.TEMP;};
    ##if (i.sheet == 12){ if (SW$number.of.individuals.include.sw  == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.NUMBER.OF.INDIVIDUALS  not read");};  DATA$SOURCES.NUMBER.OF.INDIVIDUALS   = DATA.TEMP;};
    #if (i.sheet == 11){ if (SW$concentration.stddev.include.sw  == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("MIXTURES.CONCENTRATION.STDDEV  not read");};  DATA$MIXTURES.CONCENTRATION.STDDEV   = DATA.TEMP;};
    #if (i.sheet == 12){ if (SW$concentration.stddev.include.sw  == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.CONCENTRATION.STDDEV   not read");};  DATA$SOURCES.CONCENTRATION.STDDEV    = DATA.TEMP;};
    #if (i.sheet == 13){ if (SW$efficiency.stddev.include.sw     == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("MIXTURES.EFFICIENCY.STDDEV     not read");};  DATA$MIXTURES.EFFICIENCY.STDDEV      = DATA.TEMP;};
    #if (i.sheet == 14){ if (SW$efficiency.stddev.include.sw     == 1) { write.out("r"); DATA.TEMP = as.matrix( gdata::read.xls(filename, sheet=i.sheet, verbose=FALSE, header=FALSE, perl=Perl.path) );} else { DATA.TEMP = as.matrix("SOURCES.EFFICIENCY.STDDEV      not read");};  DATA$SOURCES.EFFICIENCY.STDDEV       = DATA.TEMP;};

  } # for i.sheet
      write.out(paste("\n"));

  return( as.list(DATA) );
  ### internal variable
}
