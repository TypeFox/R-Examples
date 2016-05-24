## Set global variables to avoid CRAN check errors - see:
## http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
utils::globalVariables(c('elapsed', 'value', 'ma2', 'sampleid', 'conc', 'samptype', 'concunits', 'group_idx', 'first_col',
                         'gc_model', 'wellid', 'cut_val', 'drc_model', 'y', 'row_number', 'test_drc', 'celltype', 'growthcondition',
                         'platename', 'first', 'n', 'first_row', 'idx', 'col_id', '.', 'Elapsed', 'Date.Time'))


