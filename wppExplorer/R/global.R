utils::globalVariables("wpp.data.env")
get.indicator.choices <- function() {
	ind.names <- c('Total Fertility Rate', 'Female Life Expectancy', 'Male Life Expectancy', 
					'Total Population', 'Female Population', 'Male Population', 
					'Net Migration', 'Net Migration Rate', 
					'Sex Ratio at Birth', 'Median Age', 'Mean Age at Childbearing', 'Mean Age of Women in Childbearing Ages',
					'Total Dependency Ratio', 'Child Dependency Ratio', 'Old-age Dependency Ratio','Potential Support Ratio',
					'Mean Annual Population Growth',
					'Population by sex and age', 'Mortality Rate by sex and age', 'Age-specific Fertility Rate', 'Percent Age-specific Fertility')
	ind.def <- c('', '', '',
				'Total population in thousands', 'Female population in thousands', 'Male population in thousands', 
				'Net migration counts in thousands per 5 years', 'Annual net migration rate (per thousand population; denominator is approx. average population)',
				'Ratio of male to female', '', 'Mean age of mothers at the birth of their children', 'Mean age of women that are in childbearing ages',
				'Ratio of population age 0-14 and 65+ to population age 15-64',
				'Ratio of population age 0-14 to population age 15-64', 
				'Ratio of population age 65+ to population age 15-64', 
				'Ratio of population age 15-64 to population age 65+', 
				'log(P_t/P_{t-1})/5', 'Population in thousands', '', '', '')
	funcs <- c('fert', 'leF', 'leM', 'tpop', 'tpopF', 'tpopM', 'mig', 'migrate',
				'sexratio', 'medage', 'meanagechbear', 'meanageinchbearage',
				'tdratio', 'chdratio', 'oadratio', 'psratio',
				'popgrowth',
				'popagesex', 'mortagesex', 'fertage', 'pfertage')
	# if a new indicator is added, change also the condition in ui.R for displaying age-specific stuff
	l <- length(ind.names)
	ini <- rep(FALSE, l)
	ind.df <- data.frame(by.age=ini, no.age.sum=ini, sum.in.table=ini, low.high=ini, prob.ci=ini, mid.years=ini,
							digits=rep(2, l)) 
	rownames(ind.df) <- funcs
	ind.df[c('popagesex', 'mortagesex', 'fertage', 'pfertage'), 'by.age'] <- TRUE  # display sex and age menu
	ind.df[c('mortagesex','fertage'), 'no.age.sum'] <- TRUE                        # don't allow multiple age- and sex-selection
	ind.df[c('tpop', 'tpopF', 'tpopM', 'mig','popagesex'), 'sum.in.table'] <- TRUE # show sum in the trend table
	ind.df[c('fert', 'leF', 'leM', 'tpop', 'popagesex'), 'low.high'] <- TRUE       # has uncertainty
	ind.df[c('fert', 'leF', 'leM', 'mig', 'sexratio', 'mortagesex', 'fertage', 'pfertage'), 'mid.years'] <- TRUE # use mid years in slider (not implemented)
	ind.df[c('tpop', 'tpopF', 'tpopM','popagesex'), 'digits'] <- 0                 # number of digits in trend table (not implemented)
	ind.df['mortagesex', 'digits'] <- 4
	ind.df['fertage', 'digits'] <- 3
	ind.df['pfertage', 'digits'] <- 1
	
	structure(
		as.character(1:length(ind.names)),
		names = ind.names,
		definition = ind.def,
		settings = ind.df
	)
}


assign("wpp.data.env", new.env(), envir=parent.env(environment())
	#envir = .GlobalEnv
	)
data('iso3166', envir=wpp.data.env)
wpp.data.env$indicators <- get.indicator.choices()
wpp.data.env$package <- "wpp2015"
# Filter out non-used countries
do.call('data', list("popM", package=wpp.data.env$package, envir=wpp.data.env))
wpp.data.env$iso3166 <- wpp.data.env$iso3166[is.element(wpp.data.env$iso3166$uncode, wpp.data.env$popM$country_code),]
