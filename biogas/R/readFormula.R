# Modified: 4 April 2016 SDH

readFormula <-
function(
  form,
  elements = NULL,        # Set of elements returned, all others ignored, e.g., c('C', 'H', 'N', 'O')
  min.elements = NULL,    # Minimum set of elements, will return error if these at least are not included
  cdigits = 6,
  value = 'numeric'       # Type of output, 'numeric' for named vector, 'shortform' for shortened formula
  ) {
  #fc <- as.list(rep(0, length(elements)))

  # For the rare case with something like C0.00001, to avoid C1e-5 which will result in an error
  oldscipen <- options(scipen = 999)
  on.exit(options(scipen = oldscipen))

  form.orig <- form
  #form <- toupper(form)

  # Remove spaces
  form <- gsub(' ', '', form)

  # Add implied coefficients of 1 (also after ")")
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form)
  form <- gsub('([a-zA-Z\\)])([A-Z\\)\\(])', '\\11\\2', form) # Repeated for e.g., COOH
  form <- gsub('([a-zA-Z\\)])$', '\\11', form)

  # Find parentheses and remove them, multipying coefficients inside by coefficient at end
  # So (CH2)2 ---> C2H4
  # First add ( after N), e.g., (CH2)2CH3 ---> (CH2)2(CH3 for separation below
  form <- gsub('(\\)[0-9\\.]+)', '\\1(', form)
  # Drop extra (
  form <- gsub('^\\(', '', form)
  form <- gsub('\\($', '', form)
  form <- gsub('\\(\\(', '(', form)
  s1  <- strsplit(form, '\\(')[[1]]

  # Build up elementwise formula piecewise
  formpw <- NULL
  for(i in 1:length(s1)) {
    xx <- s1[i]
    if(grepl('\\)', xx)) {
      nn <- as.numeric(gsub('.+\\)','',xx))
      ff <- gsub('\\).+', '', xx)
      cc <- nn*as.numeric(strsplit(ff, '[A-Za-z]+')[[1]][-1])
      ee <- strsplit(ff, '[0-9.]+')[[1]]
      formpw <- paste0(formpw, paste0(ee, cc, collapse = ''))
    } else {
      formpw <- paste0(formpw, xx)
    }
  }

  form <- formpw

  # Extract integer coefficients 
  cc <- as.numeric(strsplit(form, '[A-Za-z]+')[[1]][-1])
  names(cc) <- strsplit(form, '[0-9.]+')[[1]]

  # Sort out elements to return
  if(is.null(elements)) elements <- unique(names(cc))
  fc <- numeric(length(elements))
  names(fc) <- elements

  # Fill in fc, summing elements of cc if required (if elements are repeated)
  for(i in elements) {
    for(j in 1:length(cc)) {
      if(names(cc)[j]==i) fc[i] <- fc[i] + cc[j]
    }
  }

  # Simplify form based on fc (for output only)
  form <- paste0(names(fc), signif(fc/min(fc), cdigits), collapse = '')
  # And drop coefficients of 1
  form <- gsub('([a-zA-Z])1([a-zA-Z])', '\\1\\2', form)
  form <- gsub('([a-zA-Z])1$', '\\1\\2', form)
  
  # Check for minimum set of elements
  if(!is.null(min.elements)) if(any(!min.elements %in% names(fc)) | any(fc[min.elements] == 0)) stop('Minimum elements required are ', min.elements, ' (from min.elements argument), but form is ', form.orig, ', interpreted as ', form)

  if(value == 'numeric') return(fc)
  if(value == 'shortform') as.vector(form)

}
