#' Madeyski15SQJ.NDC data
#'
#' If you use this data set please cite: Lech Madeyski and Marian Jureczko, "Which Process Metrics Can Significantly Improve Defect Prediction Models? An Empirical Study," Software Quality Journal, 2015. DOI: 10.1007/s11219-014-9241-7
#'
#' "This paper presents an empirical evaluation in which several process metrics were investigated in order to identify the ones which significantly improve the defect prediction models based on product metrics. Data from a wide range of software projects (both, industrial and open source) were collected. The predictions of the models that use only product metrics (simple models) were compared with the predictions of the models which used product metrics, as well as one of the process metrics under scrutiny (advanced models). To decide whether the improvements were significant or not, statistical tests were performed and effect sizes were calculated. The advanced defect prediction models trained on a data set containing product metrics and additionally Number of Distinct Committers (NDC) were significantly better than the simple models without NDC, while the effect size was medium and the probability of superiority (PS) of the advanced models over simple ones was high (p=.016, r=-.29, PS=.76), which is a substantial finding useful in defect prediction. A similar result with slightly smaller PS was achieved by the advanced models trained on a data set containing product metrics and additionally all of the investigated process metrics (p=.038, r=-.29, PS=.68). The advanced models trained on a data set containing product metrics and additionally Number of Modified Lines (NML) were significantly better than the simple models without NML, but the effect size was small (p=.038, r=.06). Hence, it is reasonable to recommend the NDC process metric in building the defect prediction models." [http://dx.doi.org/10.1007/s11219-014-9241-7]
#'
#'
#' @format A data frame with variables:
#' \describe{
#' \item{Project}{In case of open source projects this field includes the name of the project as well as its version. In case of industrial projects this field includes the string "properietary" (we were not allowed to disclose the names of the analyzed industrial software projects developed by Capgemini Polska).}
#' \item{simple}{The percentage of classes that must be tested in order to find 80\% of defects in case of simple defect prediction models, i.e., using only software product metrics as predictors.}
#' \item{advanced}{The percentage of classes that must be tested in order to find 80\% of defects in case of advanced defect prediction models, using not only software product metrics but also the NDC (Number of distinct committers) process metric.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15SQJ.NDC
#'
"Madeyski15SQJ.NDC"

#' Madeyski15EISEJ.OpenProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' This paper presents an analysis of 84 versions of industrial, open-source and academic projects. We have empirically evaluated whether those project types constitute separate classes of projects with regard to defect prediction. The predictions obtained from the models trained on the data from the open source projects were compared with the predictions from the other models (built on proprietary, i.e. industrial, student, open source, and not open source projects).
#'
#'
#' @format A data frame with variables:
#' \describe{
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{NOTOPEN}{The percentage of classes of projects which are not open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on open source projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.OpenProjects
#'
"Madeyski15EISEJ.OpenProjects"

#' Madeyski15EISEJ.PropProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' @format A data frame with variables:
#' \describe{
#' \item{NOTPROP}{The percentage of classes of non-proprietary (i.e., non-industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on proprietary (i.e., industrial) projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.PropProjects
#'
"Madeyski15EISEJ.PropProjects"


#' Madeyski15EISEJ.StudProjects data
#'
#' If you use this data set please cite: Marian Jureczko and Lech Madeyski, "Cross-Project Defect Prediciton: An Empirical Study," (under review), 2015.
#'
#' @format A data frame with variables:
#' \describe{
#' \item{PROP}{The percentage of classes of proprietary (i.e., industrial) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{NOTSTUD}{The percentage of classes of projects which are not student projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{STUD}{The percentage of classes of student (i.e., academic) projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' \item{OPEN}{The percentage of classes of open source projects that must be tested in order to find 80\% of defects in case of software defect prediction models built on student (i.e., academic) projects.}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Madeyski15EISEJ.StudProjects
#'
"Madeyski15EISEJ.StudProjects"

#' Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR data form a set of primary studies on reading methods for software inspections. They were reported and analysed by M. Ciolkowski ("What do we know about perspective-based reading? an approach for quantitative aggregation in software engineering", in Proceedings of the 3rd International Symposium on Empirical Software Engineering and Measurement, ESEM'09, pp. 133-144, IEEE Computer Society, 2009), corrected and re-analysed by Madeyski and Kitchenham ("How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted)).
#'
#' If you use this data set please cite: Lech Madeyski and Barbara Kitchenham, "How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted).
#'
#' @format A data frame with 21 rows and 7 variables:
#' \describe{
#' \item{Study}{Name of empirical study}
#' \item{Ref.}{Reference to the paper reporting primary study or experimental run where data were originally reported}
#' \item{Control}{Control treatment: Check-Based Reading (CBR) or Ad-hoc Reading (AR)}
#' \item{Within-subjects}{Yes - if the primary study used the within-subjects experimental design, No - if the primary study did not use the within-subjects experimental design}
#' \item{Cross-over}{Yes - if the primary study used the cross-over experimental design, No - if the primary study did not use the cross-over experimental design}
#' \item{d_ByCiolkowski}{d effect size calculated by Ciolkowski}
#' \item{d_ByOriginalAuthors}{d effect size as reported by the original authors}
#' }
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR
#'
"Ciolkowski09ESEM.MetaAnalysis.PBRvsCBRorAR"


#' MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR data form a set of primary studies on reading methods for software inspections. They were analysed by L. Madeyski and B. Kitchenham ("How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted)).
#'
#' If you use this data set please cite: Lech Madeyski and Barbara Kitchenham, "How variations in experimental designs impact the construction of comparable effect sizes for meta-analysis", 2015 (to be submitted).
#'
#' @format A data frame with 17 rows and 26 variables:
#' \describe{
#' \item{Study}{Name of empirical study}
#' \item{Ref.}{Reference to the paper reporting primary study or experimental run where data were originally reported}
#' \item{Teams}{The number of teams including both, PBR and Control teams}
#' \item{DesignDesc}{Experimental design description: Before-after, Between-groups, Cross-over}
#' \item{ExpDesign}{Experimental design: between-groups (BG), within-subjects cross-over (WSCO), within-subjects before-after (WSBA)}
#' \item{M_PBR}{The average proportion of defects found by teams using PBR}
#' \item{M_C}{The average proportion of defects found by teams using Control treatment: Check-Based Reading (CBR) or Ad-Hoc Reading (AR)}
#' \item{Diff}{The difference between M_PBR and M_C, i.e. Diff = M_PBR - M_C}
#' \item{Inc}{The percentage increase in defect rate detection, i.e. Inc=100*[(M_PBR-M_C)/M_C]}
#' \item{SD_C_ByAuthors}{The standard deviation of the control group values reproted by the original Authors, i.e., obtained from the papers/raw data}
#' \item{SD_C}{The standard deviation of the control group values equals SD_C_ByAuthors for studies for which the data was available OR the weighted average of SD_C_ByAuthors (i.e., 0.169) for studies where SD_C_ByAuthors is missing.}
#' \item{V_C}{The variance of the Control group observations, i.e., the variance obtained from the teams using the Control method V_C=SD_C^2}
#' \item{V_D}{The variance of the unstandardized mean difference D (between the mean value for the treatment group and the mean value for the Control group)}
#' \item{SD_C_Alt}{This is the equivalent of SD_C (the standard deviation of the control group) based on a different variance for the student studies or the practitioner studies depending on the subject type of the study with the missing value.}
#' \item{V_Alt}{The variance of the mean difference in the meta-analysis based on SD_C_Alt}
#' \item{SS_C}{The sum of squares of the Control group values. For within subjects studies SS=V_C*(n-1). For between subjects studies SS=V_C*(n_C-1)}
#' \item{n_PBR}{The number of PBR teams}
#' \item{n_C}{The number of Control (CBR or AR) teams}
#' \item{ControlType}{Type of Control treatment: CRB or AR}
#' \item{ParticipantsType}{Type of participants: Engineers or Students}
#' \item{TeamType}{Type of team: Nominal or Real}
#' \item{TwoPersonTeamVsLargerTeam}{Reflects size of the teams: 2-PersonTeam or LargerTeam}
#' \item{ArtefactType}{The type of artefact: Requirements or Other}
#' \item{AssociatedWithBasili}{Whether study is associated with Basili (the forerunner): Yes or No}
#' \item{ControlType_Basili}{Combined ControlType and AssociatedWithBasili: AH_AssociatedWithBasili, CBR_AssociatedWithBasili, CBR_NotAssociatedWithBasili}
#' }
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR
#'
"MadeyskiKitchenham.MetaAnalysis.PBRvsCBRorAR"


#' KitchenhamMadeyskiBudgen16.FINNISH data
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#'Data set collected from 9 Finish companies by Mr Hanna M\"aki from the TIEKE organisation see Barbara Kitchneham and Kari K\"{a}ns\"{a}l\"{a}, Inter-item correlations among function points, Proceedings ICSE 15, 1983, pp 477-480
#'
#' @format A data frame with variables:
#' \describe{
#' \item{Project}{Project ID}
#' \item{DevEffort}{Development Effort measured in hours}
#' \item{UserEffort}{Effort provided by the customer/user organisation measured in hours}
#' \item{Duration}{Project duration measurted in months}
#' \item{HWType}{A catagorical variable defining the hardware type}
#' \item{AppType}{A categorical variable definiting the application type}
#' \item{FP}{Function Points measured using the TIEKE organisation method}
#' \item{Co}{A categorical variable defining the company}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.FINNISH
#'
"KitchenhamMadeyskiBudgen16.FINNISH"

#' KitchenhamMadeyskiBudgen16.PolishSubjects data
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#' Data set collected at Wroclaw University of Technology (POLAND) by Lech Madeyski includes separate entries for each abstract assessed by a judge, that is 4 entries for each judge. Data collected from 16 subjects recruited from Wroclaw  University of Technology who were each asked to assess 4 abstracts.
#'
#' Note Only completeness question 2 was expected to be context dependent and have a NA (not applicable)  answer, if other completeness answers were left blank, BAK coded the answer as NA
#'
#' polishsubjects.txt
#' @format A data frame with variables:
#' \describe{
#' \item{Judge}{The identifier for each subject}
#' \item{Abstract}{The identifier for each abstract - the code starts with a three alphanumeric string that defines the source of the abstract}
#' \item{OrderViewed}{Each judge assessed 4 abstracts in sequence, this data item identifies the order in which the subject viewed the specified abstract}
#' \item{Completness1}{Assessment by judge of question 1:Is the reason for the project clear? Can take values:  Yes/No/Partly}
#' \item{Completness2}{Assessment by judge of question 2: Is the specific aim/purpose of the study clear? Can take values:  Yes/No/Partly}
#' \item{Completness3}{Assessment by judge of question 3: If the aim is to describe a new or enhanced software technology (e.g. method, tool, procedure or process) is the method used to develop this technology defined? Can take values:  Yes/No/Partly/NA}
#' \item{Completness4}{Assessment by judge of question 4: Is the form (e.g. experiment, general empirical study, data mining, case study, survey, simulation etc.) that was used to evaluate the technology made clear? Can take values:  Yes/No/Partly}
#' \item{Completness5}{Assessment by judge of question 5: Is there a description of how the evaluation process was organised? Can take values:  Yes/No/Partly}
#' \item{Completness6}{Assessment by judge of question 6: Are the results of the evaluation clearly described? Can take values:  Yes/No/Partly}
#' \item{Completness7}{Assessment by judge of question 7: Are any limitations of the study reported?:  Yes/No/Partly}
#' \item{Completness8}{Assessment by judge of question 8: Are any ideas for future research presented?:  Yes/No/Partly}
#' \item{Clarity}{Assessment by judge of question regarding the overall understandability of the abstract: Please give an assessment of the clarity of this abstract by circling a number on the scale of 1-10 below, where a value of 1 represents Very Obscure and 10 represents Extremely Clearly Written.}
#' \item{Completness1NumValue}{A numerical value for completeness question 1 where 0=No, Partly=0.5, yes =1}
#' \item{Completness2NumValue}{A numerical value for completeness question 2 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Completness3NumValue}{A numerical value for completeness question 3 where 0=No, Partly=0.5, yes =1, NA means not applicable or not answered}
#' \item{Completness4NumValue}{A numerical value for completeness question 4 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Completness5NumValue}{A numerical value for completeness question 5 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Completness6NumValue}{A numerical value for completeness question 6 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Completness7NumValue}{A numerical value for completeness question 7 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Completness8NumValue}{A numerical value for completeness question 8 where 0=No, Partly=0.5, yes =1, NA means not applicable}
#' \item{Sum}{The sum of the numerical completeness questions excluding those labelled NA}
#' \item{TotalQuestions}{The count of the number of question related to completeness excluding questions considered not applicable }
#' \item{Completeness}{Sum/TotalQuestions}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.PolishSubjects
#'
"KitchenhamMadeyskiBudgen16.PolishSubjects"


#' KitchenhamMadeyskiBudgen16.SubjectData
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#'Data set collected from 16 judges assessing 4 abstracts at 6 sites: Lincoln University NZ=1, Hong Kong Polytechnic University=2, PSu Thailand=3, Durham=4, Keele=5, Hong Kong City University=6
#'
#' subjectdata.txt: Judge	Institution  JudgeID age eng1st years.study abs.read Absid Treat TreatID Order Com.1 Com.2 Com.3 Com.4 Com.5 Com.6
#' Com.7 Com.8 Clarity num.questions total.score av.score Site
#' @format A data frame with variables:
#' \describe{
#' \item{Judge}{Alphanumeric identifier for each judge}
#' \item{Institution}{Numerical value identifying each site from which data was collected}
#' \item{JudgeID}{Numerical value odetifying eacjh judge}
#' \item{Age}{Age of the judge in years}
#' \item{Eng1st}{Whether the judge's first langauage was Enlish: Yes/No}
#' \item{YearsStudy}{The number of years have student been studying computing at University: 1, 2, 3, 4}
#' \item{AbstractsRead}{Number of abstracts the judge had read prior to the study" 0, 1 to 10, 10+}
#' \item{AbstractsWritten}{Whether the judge had ever written an abstract for a scientific report/article}
#' \item{AbstractID}{Alphanumeric identifier for an abstract. The first character identies the journal, I=IST, J=JSS, the third digit identifies the time period as 1 or 2, the remaining digits identify the abstract number within the set of abstracts found for the specified journal and time period}
#' \item{Treat}{The initial 3 characters of AbstractID}
#' \item{TreatID}{A numeric identifier for the journal and time period, 1=IB1, 2=IB2, 3=JB1, 4=JB2}
#' \item{Order}{The order in which the judge should have viewed the specified abstract}
#' \item{Completness1NumValue}{The numeric answer to completeness question 1}
#' \item{Completness2NumValue}{The numeric answer to completeness question 2}
#' \item{Completness3NumValue}{The numeric answer to completeness question 3}
#' \item{Completness4NumValue}{The numeric answer to completeness question 4}
#' \item{Completness5NumValue}{The numeric answer to completeness question 5}
#' \item{Completness6NumValue}{The numeric answer to completeness question 6}
#' \item{Completness7NumValue}{The numeric answer to completeness question 7}
#' \item{Completness8NumValue}{The numeric answer to completeness question 8}
#' \item{Clarity}{The response to the clarity question or NA if not answered}
#' \item{NumberOfAnsweredCompletnessQuestions}{The number of completeness questions excluding those with NA}
#' \item{TotalScore}{Sum of the numeric values of the 8 completeness questions}
#' \item{MeanScore}{Sum of the completeness questions 1 to 8 divided by TotalScore}
#' \item{Site}{The name of the site which provided the data. HongKong refers to the Polytechnic University, HongKong.2 refers to the City University}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.SubjectData
#'
"KitchenhamMadeyskiBudgen16.SubjectData"


#' KitchenhamMadeyskiBudgen16.PolishData data
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#' Data set derived from PolishSubjects data set collected at Wroclaw University. It summarizes the completeness and clarity data collected from 4 judges about the same abstract.
#'
#' PolishData.txt
#' @format A data frame with variables:
#' \describe{
#' \item{Abstract}{The abstract identifier}
#' \item{Site}{Numeric identifier for the site}
#' \item{Treatment}{The first three characters of the Abstract field which identies the jounral and time period of the abstract}
#' \item{Journal}{An acronym for the journal from which the abstarct was obtained: IST or JSS}
#' \item{Timeperiod}{The Time period in which the abstarct was found: 1 or 2}
#' \item{J1}{The identifer for the judge who made the next 2 assessments}
#' \item{J1Completeness}{The average completeness made by judge J1 based on the 8 completeness questions}
#' \item{J1Clarity}{The clarity assessment made by judge J1}
#' \item{J2}{The identifer for the judge who made the next 2 assessments}
#' \item{J2Completeness}{The average completeness made by judge J2 based on the 8 completeness questions}
#' \item{J2Clarity}{The clarity assessment made by judge J2}
#' \item{J3}{The identifer for the judge who made the next 2 assessments}
#' \item{J3Completeness}{The average completeness made by judge J3 based on the 8 completeness questions}
#' \item{J3Clarity}{The clarity assessment made by judge J3}
#' \item{J4}{The identifer for the judge who made the next 2 assessments}
#' \item{J4Completeness}{The average completeness made by judge J4 based on the 8 completeness questions}
#' \item{J4Clarity}{The clarity assessment made by judge J4}
#' \item{MedianCompleteness}{The median of J1Completeness, J2Completeness, J3Completeness, J4Completeness}
#' \item{MedianClarity}{The median of J1Clarity, J2Clarity, J3Clarity, J4Clarity}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.PolishData
#'
"KitchenhamMadeyskiBudgen16.PolishData"


#' KitchenhamMadeyskiBudgen16.DiffInDiffData data
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#'Data set was derived from the data reported in the SubjectData data set (subjectdata.txt). It contains the summary completeness and clarity data from 4 judges who assessed the same abstract. Only the initial 5 sites are included.
#'
#' dinddata.txt
#' @format A data frame with variables:
#' \describe{
#' \item{Abstract}{The abstract identifier}
#' \item{Site}{A numeric identifier of the site}
#' \item{Treatment}{A three character alphanumeric identifying the jounral and time period of the abstract}
#' \item{Journal}{The journal in which the abstract was published: IST or JSS}
#' \item{Timeperiod}{The time period in which the abstract: 1 or 2}
#' \item{J1}{The identifer for the judge who made the next 2 assessments}
#' \item{J1Completeness}{The average completeness made by judge J1 based on the 8 completeness questions}
#' \item{J1Clarity}{The clarity assessment made by judge J1}
#' \item{J2}{The identifer for the judge who made the next 2 assessments}
#' \item{J2Completeness}{The average completeness made by judge J2 based on the 8 completeness questions}
#' \item{J2Clarity}{The clarity assessment made by judge J2}
#' \item{J3}{The identifer for the judge who made the next 2 assessments}
#' \item{J3Completeness}{The average completeness made by judge J3 based on the 8 completeness questions}
#' \item{J3Clarity}{The clarity assessment made by judge J3}
#' \item{J4}{The identifer for the judge who made the next 2 assessments}
#' \item{J4Completeness}{The average completeness made by judge J4 based on the 8 completeness questions}
#' \item{J4Clarity}{The clarity assessment made by judge J4}
#' \item{MeanCompleteness}{The mean of J1Completeness, J2Completeness, J3Completeness, J4Completeness}
#' \item{MedianCompleteness}{The median of J1Completeness, J2Completeness, J3Completeness, J4Completeness}
#' \item{MedianClarity}{The median clarity of J1Clarity, J2Clarity, J3Clarity, J4Clarity}
#' \item{MeanClarity}{The mean clarity of J1Clarity, J2Clarity, J3Clarity, J4Clarity}
#' \item{VarCompleteness}{The variance of J1Completeness, J2Completeness, J3Completeness, J4Completeness}
#' \item{VarClarity}{The variance clarity of J1Clarity, J2Clarity, J3Clarity, J4Clarity}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.DiffInDiffData
#'
"KitchenhamMadeyskiBudgen16.DiffInDiffData"


#' KitchenhamMadeyskiBudgen16.COCOMO data
#'
#' If you use this data set please cite this R package and the following paper when accepted: Barbara Kitchenham, Lech Madeyski, David Budgen, Jacky Keung et al. "Robust Statistical Methods for Empirical Software Engineering".
#'
#'Data set collected at TRW by Barry Boehm see: B.W. Boehm. 1981.  Software Engineering Economics. Prentice-Hall.
#'
#'Explanations by Barbara Kitchehnam / https://terapromise.csc.ncsu.edu:8443/!/#repo/view/head/effort/cocomo/cocomo1/nasa93/nasa93.arff
#'
#' COCOMO.txt: pro	type	year	Lang	Rely	Data	CPLX	aaf	time	store	virt	turn	type2	acap	aexp	pcap	vexp	lexp	cont	modp	TOOL	TOOLcat	SCED	RVOL	Select	rvolcat	Modecat	Mode1	Mode2	Mode3	KDSI	AKDSI	Effort	Dur	Productivity
#' @format A data frame with variables:
#' \describe{
#' \item{Project}{Project ID}
#' \item{Type}{A categorical variable describing the type of the project}
#' \item{Year}{The year the project was completed}
#' \item{Lang}{A categorical variable describing the development language used}
#' \item{Rely}{Ordinal value defining the required software reliability}
#' \item{Data}{Ordinal value defining the data complexity / Data base size}
#' \item{Cplx}{Ordinal value defining the complexity of the software / Process complexity}
#' \item{Aaf}{??}
#' \item{Time}{Ordinal value defining the stringency of timing constraints / Time constraint for cpu}
#' \item{Stor}{Ordinal value defining the stringency of the data storage requirements / Main memory constraint}
#' \item{Virt}{Virtual Machine volatility}
#' \item{Turn}{Turnaround time}
#' \item{Type2}{A categorical variable defining the hardware type: mini, max=mainframe, midi}
#' \item{Acap}{Ordinal value defining the analyst capability}
#' \item{Aexp}{Ordinal value defining the analyst experience / application experience}
#' \item{Pcap}{Ordinal value defining the programming capability of the team / Programmers capability}
#' \item{Vexp}{Ordinal value defining the virtual machine experience of the team}
#' \item{Lexp}{Ordinal value defining the programming language experience of the team}
#' \item{Cont}{??}
#' \item{Modp}{ / Modern programing practices}
#' \item{Tool}{Ordinal value defining the extent of tool use / Use of software tools}
#' \item{ToolCat}{Recoding of Tool to labelled ordinal scale}
#' \item{Sced}{Ordinal value defining the stringency of the schedule requirements / Schedule constraint}
#' \item{Rvol}{Ordinal value defining the requriements volatility of the project}
#' \item{Select}{Categorical value calculated by BAK for an analysis example}
#' \item{Rvolcat}{Recoding of Rvol to a labelled ordinal scale}
#' \item{Modecat}{Mode of the projects: O=Organic, E=Embedded, SD-Semi-Detached}
#' \item{Mode1}{Dummy variable calculated by BAK: 1 if the project is Organic, 0 otherwise}
#' \item{Mode2}{Dummy variable calculated by BAK: 1 if the project is Semi-detached, 0 otherwise}
#' \item{Mode3}{Dummy variable calculated by BAK: 1 if the project is Embedded, 0 otherwise}
#' \item{KDSI}{Product Size Thousand of Source Instructions}
#' \item{AKDSI}{Adjusted Product Size for Project in Thousand Source Instructions - differs from KDSI for enhancement projects}
#' \item{Effort}{Project Effort in Man months}
#' \item{Duration}{Duration in months}
#' \item{Productivity}{Productivity of project calculated by BAK as AKDSI/Effort, so the the larger the value the better the productivity}
#' }
#'
#' @source \url{http://madeyski.e-informatyka.pl/reproducible-research/}
#' @examples
#' KitchenhamMadeyskiBudgen16.COCOMO
#'
"KitchenhamMadeyskiBudgen16.COCOMO"
