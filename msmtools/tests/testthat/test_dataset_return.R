
# TESTING 2 OR 3 VALUES FOR PATTERN -----------------------------------------------------------

data( hosp )
hosp_aug_2 = augment( hosp, subj, adm_number, label_2,
                      t_start = dateIN, t_end = dateOUT, t_cens = dateCENS, verbose = F )

hosp_aug_3 = augment( hosp, subj, adm_number, label_3,
                      t_start = dateIN, t_end = dateOUT, t_cens = dateCENS, verbose = F )

test_that( "passing to pattern a var with 2 or 3 values is identical",
           expect_identical( hosp_aug_2, hosp_aug_3 ) )
#
#
#
# # TESTING IF N_EVENTS IS MISSING --------------------------------------------------------------
#
data( hosp )
hosp_aug = augment( hosp, subj, adm_number, label_3,
                    t_start = dateIN, t_end = dateOUT, t_cens = dateCENS, verbose = F )
hosp_aug_no_events = augment( hosp, subj, pattern = label_3,
                              t_start = dateIN, t_end = dateOUT, t_cens = dateCENS, verbose = F )

test_that( "adm_number is identical to what augment created",
           expect_identical( hosp_aug$adm_number, hosp_aug_no_events$n_events ) )

