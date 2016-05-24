# edit definitions
# Derived from Basic questionnaire list of Dutch Labour Force Survey (2009)
# The edits follow routing of a phone interview, consisting of a 
# number of questions on the household ('householdbox') and a number
# of questions on payed work ('workbox')
#
# Names of variables and categories are loosely translated from Dutch.
# mvdl, 12.12.2011


BOOL <- c(TRUE, FALSE)
AGE <- c('0-11','12-15','16+')

# note: this wil be treated as character
NUMHOUSEHOLD <- as.character(1:4)

HOUSEHOLDTYPE <- c(
    'partners, no children',
    'partners with children',
    'partners with children and others',
    'partners with others',
    'single parent with children',
    'single parent with children and others',
    'other')

GENDER <- c('male','female')

MARITALSTATUS <- c(
    'married',
    'divorced',
    'widowed',
    'never married')

RELATION <- c(
    'child',
    'parent',
    'parent in law',
    'sibling',
    'sibling in law',
    'grandchild',
    'other, family',
    'other, not family')

PARTNERRELATION <- c(
    'spouse',
    'partner',
    'not partner')

WORKENVIRONMENTS <- c('none','single','multiple')
HOURSPERWEEK <- c('0-4','5-11','12-29','30+')

COMPANYOWNER <- c('self','partner','in-laws')

COMPANYCAT <- c(
    'agro',
    'energy/water',
    'construction',
    'education',
    'healthcare',
    'government',
    'culture, sports, tourism',
    'horeca',
    'trade/retail',
    'finacial',
    'transport and communication',
    'service',
    'industry',
    'other')

PROFESSIONS <- c(
    'teacher/instructor',
    'agro',
    'science, math',
    'technical',
    'transport',
    'medical',
    'business administration',
    'law, management',
    'horeca, wellness',
    'language, cultural',
    'social',
    'other')

WANTWORK <-c('no','yes','found','cant')


# Establishing interview
# ap = answering person
age_ap %in% AGE
interview %in% BOOL

if ( age_ap == '0-11' ) proxyavailable %in% BOOL
if ( age_ap == '12-15') permission %in% BOOL

if ( proxyavailable ) interview
if ( permission ) interview 
if ( age_ap == "16+" ) interview 




# Establish household type
if (interview) numhousehold %in% NUMHOUSEHOLD
if (interview) householdtype %in% HOUSEHOLDTYPE
if (interview) respondentinhhcore %in% BOOL

if (numhousehold == '1') householdtype == 'other'

if ( numhousehold == '2' ) 
    householdtype %in% c( 
        'partners, no children', 
        'partners with others', 
        'single parent with children','other')

if ( numhousehold %in% '3' ) 
    householdtype %in% c(
       'partners with children',
       'partners with others',
       'single parent with children',
       'single parent with children and others',
       'other')



# data on first person (head of household)
if (interview) gender_1 %in% GENDER
if (interview) age_1 %in% AGE
  # if respondent is in hh core, age_1 == age.
if (respondentinhhcore) age_1 == '16+'
if (interview) maritalstatus_1 %in% MARITALSTATUS
if ( age_1 %in% c('0-11','12-15') ) maritalstatus_1 == 'never married'


# 2nd person, if applicable
if ( numhousehold %in% as.character(2:4) ) gender_2 %in% GENDER
if ( numhousehold %in% as.character(2:4) ) 
    relation_2 %in% c(PARTNERRELATION,RELATION)

if (relation_2 == 'spouse') maritalstatus_1 == 'married'

if (householdtype %in% c(
    'single parent with children',
    'single parent with children and others')) relation_2 %in% RELATION

if (householdtype %in% grep('partners',HOUSEHOLDTYPE,value=TRUE)) 
    relation_2 %in% PARTNERRELATION

if ( numhousehold %in% as.character(2:4) & householdtype == 'single parent with children')
    relation_2 == 'child'

if (numhousehold %in% as.character(3:4) & householdtype == 'single parent with children' )
   relation_3 == 'child' 

if (numhousehold == '4' & householdtype == 'single parent with children' )
   relation_4 == 'child' 


if (numhousehold %in% as.character(3:4) & householdtype == 'partners with children')
    relation_3 == 'child'

if (numhousehold == '4' & householdtype == 'partners with children')
    relation_4 == 'child'


# 3rd person, if applicable
if ( numhousehold %in% as.character(3:4) ) gender_3 %in% GENDER
if ( numhousehold %in% as.character(3:4) ) age_3 %in% AGE
if ( numhousehold %in% as.character(3:4) ) relation_3 %in% RELATION


# 3rd person, if applicable
if ( numhousehold == '4' ) gender_4 %in% GENDER
if ( numhousehold == '4' ) age_4 %in% AGE
if ( numhousehold == '4' ) relation_4 %in% RELATION


# can the workbox be filled in?
workbox %in% BOOL
if ( !interview ) !workbox
if ( interview & age_ap == '16+' ) workbox

# we assume we are only interested in work of answering person.
if (workbox) haspayedjob %in% BOOL
if ( !haspayedjob ) owncompany %in% BOOL
if ( !owncompany ) familycompany %in% BOOL

# haswork 
haswork %in% BOOL
if ( haspayedjob ) haswork
if (!haspayedjob & owncompany ) haswork
if (!haspayedjob & familycompany ) haswork

if ( !haspayedjob & !owncompany & !familycompany ) !haswork

if ( haspayedjob ) workhours %in% HOURSPERWEEK

if ( haswork ) workenvironments %in% WORKENVIRONMENTS
if ( familycompany ) workenvironments %in% WORKENVIRONMENTS

if ( haswork ) employee %in% BOOL
if ( familycompany ) employee %in% BOOL

if ( !employee & haswork ) companyowner %in% COMPANYOWNER
if ( companyowner == 'self' ) companycat %in%  COMPANYCAT
if ( owncompany ) companycat %in% COMPANYCAT


if (employee) profession %in% PROFESSIONS
if (!employee & familycompany ) profession %in% PROFESSIONS
if ( companyowner %in% c('partner','in-laws')) profession %in% PROFESSIONS


# is someone looking for work?
unemploymentbox %in% BOOL

if (!isemployee & age_ap == '16+' ) unemploymentbox
if (workhours %in% c('0-4','5-11')) unemploymentbox

if (unemploymentbox) wantwork %in% WANTWORK

if (wantwork == 'found') willwork %in% HOURSPERWEEK
if (wantwork == 'cant' ) diffwork %in% BOOL
if (wantwork == 'yes' & !haswork) wantworkhrs %in% HOURSPERWEEK


if (wantwork=='yes'  & workhours %in% c('0-4','5-11') ) want12hrowncompany %in% BOOL
if (wantwork=='yes'  & wantworkhrs %in% c('12-29','30+')) want12hrowncompany %in% BOOL










