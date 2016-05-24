write.input <-
function# Write input values out to a file
### internal function for sisus
(DATA
### internal variable
, filename
### internal variable
, output.inputs.filename
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
  #fn = paste(".\\output\\", output.inputs.filename, sep="");
  fn = output.inputs.filename;
  h = "########"; s = ":: sheet";
  unlink(fn);   # delete old file
  i.s =  1; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SISUS.PARAMETERS              ), file = fn, ncolumns = dim(DATA$SISUS.PARAMETERS              )[2], append = TRUE, sep = ",");
  i.s =  2; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$MIXTURES.ISOTOPE.RATIOS       ), file = fn, ncolumns = dim(DATA$MIXTURES.ISOTOPE.RATIOS       )[2], append = TRUE, sep = ",");
  i.s =  3; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.ISOTOPE.RATIOS        ), file = fn, ncolumns = dim(DATA$SOURCES.ISOTOPE.RATIOS        )[2], append = TRUE, sep = ",");
  i.s =  4; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.DISCRIMINATION        ), file = fn, ncolumns = dim(DATA$SOURCES.DISCRIMINATION        )[2], append = TRUE, sep = ",");
  i.s =  5; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.CONCENTRATION         ), file = fn, ncolumns = dim(DATA$SOURCES.CONCENTRATION         )[2], append = TRUE, sep = ",");
  i.s =  6; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.EFFICIENCY            ), file = fn, ncolumns = dim(DATA$SOURCES.EFFICIENCY            )[2], append = TRUE, sep = ",");
  i.s =  7; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$LINEAR.CONSTRAINTS            ), file = fn, ncolumns = dim(DATA$LINEAR.CONSTRAINTS            )[2], append = TRUE, sep = ",");
    # prior removed at v0.09 until these ideas are developed 5/22/2008 3:42PM
  #i.s =  8; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$PRIOR.DIRICHLET               ), file = fn, ncolumns = dim(DATA$PRIOR.DIRICHLET               )[2], append = TRUE, sep = ",");
  i.s =  8; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$MIXTURES.ISOTOPE.STDDEV       ), file = fn, ncolumns = dim(DATA$MIXTURES.ISOTOPE.STDDEV       )[2], append = TRUE, sep = ",");
  i.s =  9; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.ISOTOPE.STDDEV        ), file = fn, ncolumns = dim(DATA$SOURCES.ISOTOPE.STDDEV        )[2], append = TRUE, sep = ",");
  i.s = 10; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.BIOMASS.PER.INDIVIDUAL), file = fn, ncolumns = dim(DATA$SOURCES.BIOMASS.PER.INDIVIDUAL)[2], append = TRUE, sep = ",");
  ##i.s = 12; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.NUMBER.OF.INDIVIDUALS ), file = fn, ncolumns = dim(DATA$SOURCES.NUMBER.OF.INDIVIDUALS )[2], append = TRUE, sep = ",");
 # i.s = 11; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$MIXTURES.CONCENTRATION.STDDEV ), file = fn, ncolumns = dim(DATA$MIXTURES.CONCENTRATION.STDDEV )[2], append = TRUE, sep = ",");
 # i.s = 12; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.CONCENTRATION.STDDEV  ), file = fn, ncolumns = dim(DATA$SOURCES.CONCENTRATION.STDDEV  )[2], append = TRUE, sep = ",");
 # i.s = 13; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$MIXTURES.EFFICIENCY.STDDEV    ), file = fn, ncolumns = dim(DATA$MIXTURES.EFFICIENCY.STDDEV    )[2], append = TRUE, sep = ",");
 # i.s = 14; s.n = paste(h, filename, s, i.s, h); write(s.n, file = fn, append = TRUE, sep = ","); write(t(DATA$SOURCES.EFFICIENCY.STDDEV     ), file = fn, ncolumns = dim(DATA$SOURCES.EFFICIENCY.STDDEV     )[2], append = TRUE, sep = ",");

  ### internal variable
}
