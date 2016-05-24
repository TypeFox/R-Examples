# TESTING ERRORS ------------------------------------------------------------------------------

data( hosp )

test_that( 'a warning is printed when t_cens is the only var used',
           expect_warning( augment( hosp, subj, adm_number, label_3,
                                  t_start = dateIN, t_end = dateOUT, t_cens = dateCENS,
                                  verbose = F ) ) )

test_that( 'n_events must be an integer',
           expect_error( augment( hosp, subj, !as.integer( adm_number ), label_3,
                                  t_start = dateIN, t_end = dateOUT, t_cens = dateCENS,
                                  t_death = dateCENS,
                                  verbose = F ) ) )
