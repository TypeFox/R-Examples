#This R script setups up the PISA data for North America included in this package.
pisa.colnames = as.data.frame(matrix(c(
	"CNT",      'CNT',     "Country",
	"SCHOOLID", 'SchoolId', "SchoolID",
	"StIDStd", 'StudentId', "Student ID",
	"ST01Q01", 'Grade',   "Grade",
	"ST04Q01", 'Sex',     "Sex",
	"ST05Q01", 'Attend',  "Attend",
	"ST06Q01", 'Age',     "Age",
	"ST07Q01", 'Repeat',  "Repeat",
	"ST08Q01", 'Mother',  "At home mother",
	"ST08Q02", 'Father',  "At home father",
	"ST08Q03", 'Brother', "At home brothers",
	"ST08Q04", 'Sister',  "At home sisters",
	"ST08Q05", 'GrandPa', "At home grandparents",
	"ST08Q06", 'Other',   "At home others",
	"ST10Q01", 'MomEd',   "Mother highest schooling",
	"ST12Q01", 'MomJob',  "Mother current job status",
	"ST14Q01", 'DadEd',   "Father highest schooling",
	"ST16Q01", 'DadJob',  "Father current job status",
	"ST19Q01", 'Lang',    "Language at home",
	"ST20Q01", 'Desk',    "Desk",
	"ST20Q02", 'OwnRoom', "Own room",
	"ST20Q03", 'StudyPl', "Study place",
	"ST20Q04", 'Computer',"Computer",
	"ST20Q05", 'Software',"Software",
	"ST20Q06", 'Internet',"Internet",
	"ST20Q07", 'Lit',     "Literature",
	"ST20Q08", 'Poetry',  "Poetry",
	"ST20Q09", 'Art',     "Art",
	"ST20Q10", 'TxtBooks',"Textbooks",
	"ST20Q12", 'Dict',    "Dictionary",
	"ST20Q13", 'DishW',   "Dishwasher",
	"ST20Q14", 'DVD',     "DVD",
	"ST21Q01", 'CellPh',  "How many cellphones",
	"ST21Q02", 'TVs',     "How many TVs",
	"ST21Q03", 'nComp',   "How many computers",
	"ST21Q04", 'nCars',   "How many cars",
	"ST21Q05", 'nBaths',  "How many rooms bath or shower",
	"ST22Q01", 'nBooks',  "How many books",
	"ST23Q01", 'Reading',      "Reading enjoyment time",
	"ST31Q01", 'EnrichLang',   "Enrich in test language",
	"ST31Q02", 'EnrichMath',   "Enrich in mathematics",
	"ST31Q03", 'EnrichScie',   "Enrich in science",
	"ST31Q05", 'RemedialLang', "Remedial in test language",
	"ST31Q06", 'RemedialMath', "Remedial in mathematics",
	"ST31Q07", 'RemedialScie', "Remedial in science",
	"ST32Q01", 'LangLessons',  "Out of school lessons in test language",
	"ST32Q02", 'MathLessons',  "Out of school lessons maths",
	"ST32Q03", 'ScieLessons',  "Out of school lessons in science",
	"PUBPRIV", "PublicPrivate", "Public or private school",
	"STRATIO", "StudentTeacherRatio", "Student to teacher ratio in school"
), ncol=3, byrow=TRUE))
names(pisa.colnames) = c('Variable', 'ShortDesc', 'Desc')

recode.pisa <- function(x) {
	x[,"ST04Q01"] = factor(x[,'ST04Q01'], levels=c(1,2), labels=c('Female', 'Male')) #Sex
	x[,"ST05Q01"] = factor(x[,'ST05Q01'], levels=c(1,2,3), labels=c('No', 'Yes <1 yr', 'Yes >1 yr')) #Attend
	x[,"ST07Q01"] = factor(x[,'ST07Q01'], levels=c(1,2,3), labels=c('No', 'Once', '2 or more')) #Repeat
	x[,"ST08Q01"] = factor(x[,'ST08Q01'], levels=c(1,2), labels=c('Yes', 'No')) #At home mother 1=Yes; 2=No
	x[,"ST08Q02"] = factor(x[,'ST08Q02'], levels=c(1,2), labels=c('Yes', 'No')) #At home father 1=Yes; 2=No
	x[,"ST08Q03"] = factor(x[,'ST08Q03'], levels=c(1,2), labels=c('Yes', 'No')) #At home brothers 1=Yes; 2=No
	x[,"ST08Q04"] = factor(x[,'ST08Q04'], levels=c(1,2), labels=c('Yes', 'No')) #At home sisters
	x[,"ST08Q05"] = factor(x[,'ST08Q05'], levels=c(1,2), labels=c('Yes', 'No')) #At home grandparents
	x[,"ST08Q06"] = factor(x[,'ST08Q06'], levels=c(1,2), labels=c('Yes', 'No')) #At home others
	x[,"ST10Q01"] = as.integer(x[,'ST10Q01']) #Mother highest schooling
	x[,"ST12Q01"] = factor(x[,'ST12Q01'], levels=c(1,2,3,4), labels=c('Full-time','Part-Time','Looking','Other')) #Mother current job status
	x[,"ST14Q01"] = as.integer(x[,'ST14Q01']) #Father highest schooling
	x[,"ST16Q01"] = factor(x[,'ST16Q01'], levels=c(1,2,3,4), labels=c('Full-time','Part-Time','Looking','Other')) #Father current job status
	x[,"ST19Q01"] = factor(x[,'ST19Q01'], levels=c(1,2), labels=c('Test Lang', 'Other Lang')) #Language at home
	x[,"ST20Q01"] = factor(x[,'ST20Q01'], levels=c(1,2), labels=c('Yes', 'No')) #Desk
	x[,"ST20Q02"] = factor(x[,'ST20Q02'], levels=c(1,2), labels=c('Yes', 'No')) #Own room
	x[,"ST20Q03"] = factor(x[,'ST20Q03'], levels=c(1,2), labels=c('Yes', 'No')) #Study place
	x[,"ST20Q04"] = factor(x[,'ST20Q04'], levels=c(1,2), labels=c('Yes', 'No')) #Computer
	x[,"ST20Q05"] = factor(x[,'ST20Q05'], levels=c(1,2), labels=c('Yes', 'No')) #Software
	x[,"ST20Q06"] = factor(x[,'ST20Q06'], levels=c(1,2), labels=c('Yes', 'No')) #Internet
	x[,"ST20Q07"] = factor(x[,'ST20Q07'], levels=c(1,2), labels=c('Yes', 'No')) #Literature
	x[,"ST20Q08"] = factor(x[,'ST20Q08'], levels=c(1,2), labels=c('Yes', 'No')) #Poetry
	x[,"ST20Q09"] = factor(x[,'ST20Q09'], levels=c(1,2), labels=c('Yes', 'No')) #Art
	x[,"ST20Q10"] = factor(x[,'ST20Q10'], levels=c(1,2), labels=c('Yes', 'No')) #Textbooks
	x[,"ST20Q12"] = factor(x[,'ST20Q12'], levels=c(1,2), labels=c('Yes', 'No')) #Dictionary
	x[,"ST20Q13"] = factor(x[,'ST20Q13'], levels=c(1,2), labels=c('Yes', 'No')) #Dishwasher
	x[,"ST20Q14"] = factor(x[,'ST20Q14'], levels=c(1,2), labels=c('Yes', 'No')) #DVD
	x[,"ST21Q01"] = as.integer(x[,'ST21Q01']) #How many cellphones
	x[,"ST21Q02"] = as.integer(x[,"ST21Q02"]) #How may TVs
	x[,"ST21Q03"] = as.integer(x[,"ST21Q03"]) #How many computers
	x[,"ST21Q04"] = as.integer(x[,"ST21Q04"]) #How many cars
	x[,"ST23Q01"] = as.integer(x[,"ST23Q01"]) #How many rooms bath or shower
	x[,"ST23Q01"] = as.integer(x[,"ST23Q01"]) #How many books 1=0-10; 2=11-25; 3=26-100; 4=101-200; 5=201-500; 6= >500
	x[,"ST23Q01"] = as.integer(x[,"ST23Q01"]) #Reading enjoyment time 1=Don't read; 2= <30 minutes/day; 3=1-2 hours/day; 4= >2hours/day
	x[,"ST31Q01"] = factor(x[,'ST31Q01'], levels=c(1,2), labels=c('Yes', 'No')) #Enrich in test language
	x[,"ST31Q02"] = factor(x[,'ST31Q02'], levels=c(1,2), labels=c('Yes', 'No')) #Enrich in mathematics
	x[,"ST31Q03"] = factor(x[,'ST31Q03'], levels=c(1,2), labels=c('Yes', 'No')) #Enrich in science
	x[,"ST31Q05"] = factor(x[,'ST31Q05'], levels=c(1,2), labels=c('Yes', 'No')) #Remedial in test language
	x[,"ST31Q06"] = factor(x[,'ST31Q06'], levels=c(1,2), labels=c('Yes', 'No')) #Remedial in mathematics
	x[,"ST31Q07"] = factor(x[,'ST31Q07'], levels=c(1,2), labels=c('Yes', 'No')) #Remedial in science
	x[,"ST32Q01"] = factor(x[,'ST32Q01'], levels=c(1,2,3,4,5), labels=c('Do not attend', '<2 hrs/wk', '2-4 hrs/wk', '4-6 hrs/wk', '>6 hrs/wk'))  #Out of school lessons in test language
	x[,"ST32Q02"] = factor(x[,'ST32Q02'], levels=c(1,2,3,4,5), labels=c('Do not attend', '<2 hrs/wk', '2-4 hrs/wk', '4-6 hrs/wk', '>6 hrs/wk')) #Out of school lessons <maths>
	x[,"ST32Q03"] = factor(x[,'ST32Q03'], levels=c(1,2,3,4,5), labels=c('Do not attend', '<2 hrs/wk', '2-4 hrs/wk', '4-6 hrs/wk', '>6 hrs/wk')) #Out of school lessons <science>
	
	x$mathscore = apply(x[,c("PV1MATH", "PV2MATH", "PV3MATH", "PV4MATH", "PV5MATH")], 1, sum) / 5 #Math scores
	x$readscore = apply(x[,c("PV1READ", "PV2READ", "PV3READ", "PV4READ", "PV5READ")], 1, sum) / 5 #Reading scores
	x$sciescore = apply(x[,c("PV1SCIE", "PV2SCIE", "PV3SCIE", "PV4SCIE", "PV5SCIE")], 1, sum) / 5 #Science scores
	return(x)
}

pisa.countries = as.data.frame(matrix(c(
	"ALB","Albania",
	"ARG","Argentina",
	"AUS","Australia",
	"AUT","Austria",
	"AZE","Azerbaijan",
	"BEL","Belgium",
	"BRA","Brazil",
	"BGR","Bulgaria",
	"CAN","Canada",
	"CHL","Chile",
	"TAP","Chinese Taipei",
	"COL","Colombia",
	"HRV","Croatia",
	"CZE","Czech Republic",
	"DNK","Denmark",
	"QAR","Dubai (UAE)",
	"EST","Estonia",
	"FIN","Finland",
	"FRA","France",
	"DEU","Germany",
	"GRC","Greece",
	"HKG","Hong Kong-China",
	"HUN","Hungary",
	"ISL","Iceland",
	"IDN","Indonesia",
	"IRL","Ireland",
	"ISR","Israel",
	"ITA","Italy",
	"JPN","Japan",
	"JOR","Jordan",
	"KAZ","Kazakhstan",
	"KOR","Korea",
	"KGZ","Kyrgyzstan",
	"LVA","Latvia",
	"LIE","Liechtenstein",
	"LTU","Lithuania",
	"LUX","Luxembourg",
	"MAC","Macao-China",
	"MEX","Mexico",
	"MNE","Montenegro",
	"NLD","Netherlands",
	"NZL","New Zealand",
	"NOR","Norway",
	"PAN","Panama",
	"PER","Peru",
	"POL","Poland",
	"PRT","Portugal",
	"QAT","Qatar",
	"ROU","Romania",
	"RUS","Russian Federation",
	"SRB","Serbia",
	"QCN","Shanghai-China",
	"SGP","Singapore",
	"SVK","Slovak Republic",
	"SVN","Slovenia",
	"ESP","Spain",
	"SWE","Sweden",
	"CHE","Switzerland",
	"THA","Thailand",
	"TTO","Trinidad and Tobago",
	"TUN","Tunisia",
	"TUR","Turkey",
	"GBR","United Kingdom",
	"USA","United States",
	"URY","Uruguay"
), ncol=2, byrow=TRUE))
names(pisa.countries) = c('CNT3', 'Country')
pisa.countries$CNT3 = as.character(pisa.countries$CNT3)
pisa.countries$Country = as.character(pisa.countries$Country)

