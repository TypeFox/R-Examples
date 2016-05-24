
cancerCodes = list(
		site = c('All Sites'=0,'Oral Cavity and Pharynx'=20000,'Lip'=20010,'Tongue'=20020, 
				'Salivary Gland'=20030,'Floor of Mouth'=20040,'Gum and Other Mouth'=20050,
				'Nasopharynx'=20060,'Tonsil'=20070,'Oropharynx'=20080,'Hypopharynx'=20090,
				'Other Oral Cavity and Pharynx'=20100,'Digestive System'=21000,
				'Esophagus'=21010,'Stomach'=21020,'Small Intestine'=21030,'Colon and Rectum'=1,
				'Early Stage Colon and Rectum'=6,'Late Stage Colon and Rectum'=7,
				'Colon excluding Rectum'=21040,'Cecum'=21041,'Appendix'=21042,'Ascending Colon'=21043,
				'Hepatic Flexure'=21044,'Transverse Colon'=21045,'Splenic Flexure'=21046,
				'Descending Colon'=21047,'Sigmoid Colon'=21048,'Large Intestine, NOS'=21049,
				'Rectum and Rectosigmoid Junction'=21050,'Rectosigmoid Junction'=21051,
				'Rectum'=21052,'Anus, Anal Canal, and Anorectum'=21060,
				'Liver and Intrahepatic Bile Duct'=2,'Liver'=21071,'Intrahepatic Bile Duct'=21072,
				'Gallbladder'=21080,'Other Biliary'=21090,'Pancreas'=21100,'Retroperitoneum'=21110,
				'Peritoneum, Omentum, and Mesentery'=21120,'Other Digestive Organs'=21130,
				'Respiratory System'=22000,'Nose, Nasal Cavity, and Middle Ear'=22010,
				'Larynx'=22020,
				'Lung and Bronchus'=22030,'Pleura'=22050,
				'Trachea, Mediastinum and Other Respiratory Organs'=22060,'Bones and Joints'=23000,
				'Soft Tissue including Heart'=24000,'Skin excluding Basal and Squamous'=25000,
				'Melanoma of the Skin'=25010,'Other NonEpithelial Skin'=25020,'Breast'=26000,
				'In Situ Breast'=3,'Early Stage Breast'=4,'Late Stage Breast'=5,
				'Female Genital System'=27000,'Cervix Uteri'=27010,'Early Stage Cervix Uteri'=8,
				'Late Stage Cervix Uteri'=9,'Corpus Uteri'=27020,'Uterus, NOS'=27030,
				'Ovary'=27040,'Vagina'=27050,'Vulva'=27060,'Other Female Genital Organs'=27070,
				'Male Genital System'=28000,'Prostate'=28010,'Testis'=28020,'Penis'=28030,
				'Other Male Genital Organs'=28040,'Urinary System'=29000,'Urinary Bladder'=29010,
				'Kidney and Renal pelvis'=29020,'Ureter'=29030,'Other Urinary Organs'=29040,
				'Eye and Orbit'=30000,'Brain and Other Nervous System'=31000,'Brain'=31010,
				'Cranial Nerves and Other Nervous System'=31040,'Endocrine System'=32000,
				'Thyroid'=32010,'Other Endocrine including Thymus'=32020,'Lymphoma'=33000,
				'Hodgkin Lymphoma'=33010,'Hodgkin Lymphoma  Nodal'=33011,
				'Hodgkin Lymphoma  Extranodal'=33012,'NonHodgkin Lymphoma'=33040,
				'NonHodgkin Lymphoma  Nodal'=33041,'NonHodgkin Lymphoma  Extranodal'=33042,
				'Myeloma'=34000,'Leukemia'=35000,'Lymphocytic Leukemia'=35010,
				'Acute Lymphocytic Leukemia'=35011,'Chronic Lymphocytic Leukemia'=35012,
				'Other Lymphocytic Leukemia'=35013,'Myeloid and Monocytic Leukemia'=35020,
				'Acute Myeloid Leukemia'=35021,'Acute Monocytic Leukemia'=35031,
				'Chronic Myeloid Leukemia'=35022,'Other Myeloid/Monocytic Leukemia'=35023,
				'Other Leukemia'=35040,'Other Acute Leukemia'=35041,
				'Aleukemic, Subleukemic and NOS Leukemia'=35043,
				'Mesothelioma'=36010,'Kaposi Sarcoma'=36020,'Miscellaneous'=37000),
		state = c(Texas='tx', Georgia='ga', Kentucky = 'ky',
				Michigan='mi', Arkansas='ar', Mississippi='ms',
				Wisconsin='wi', Iowa='ia', 'New Mexico' ='nm',
				Utah='ut', California='ca', Seattle='se',
				Conneticut = 'ct', 'New Jersey' = 'nj'),
		sex = c(M=1, F=2, both=0)

)	


# minnisota
# https://apps.health.state.mn.us/mndata/cancer_query?p_auth=Wy7tICix&p_p_id=springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w&p_p_lifecycle=1&p_p_state=normal&p_p_mode=view&p_p_col_id=column-1&p_p_col_count=2&p_p_col_pos=1&_springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w__facesViewIdRender=%2Fpages%2Findex.xhtml
# https://apps.health.state.mn.us/mndata/cancer_query?p_auth=Wy7tICix&p_p_id=springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w&p_p_lifecycle=1&p_p_state=normal&p_p_mode=view&p_p_col_id=column-1&p_p_col_count=2&p_p_col_pos=1&_springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w__facesViewIdRender=%2Fpages%2Findex.xhtml
# https://apps.health.state.mn.us/mndata/cancer_query?p_auth=Wy7tICix&p_p_id=springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w&p_p_lifecycle=1&p_p_state=normal&p_p_mode=view&p_p_col_id=column-1&p_p_col_count=2&p_p_col_pos=1&_springQueryPortlet_WAR_mndataspringQueryportlet_INSTANCE_TQg4MoZMNm0w__facesViewIdRender=%2Fpages%2Findex.xhtml

usCancer = function(
		state='Kentucky',
		site='Lung',
		year = c(2004,2008),
		sex='both'
		) {
			
			

			
			Ssite = cancerCodes$site[grep(paste(site, collapse='|'), names(cancerCodes$site), ignore.case=TRUE)]
			Ssex = cancerCodes$sex[grep(paste(sex, collapse='|'), names(cancerCodes$sex), ignore.case=TRUE)]
			Sstate = cancerCodes$state[grep(paste(state, collapse='|'), names(cancerCodes$state), ignore.case=TRUE)]
			
			if(is.matrix(year)) {
				SstartYear = year[,1]
				SendYear = year[,2]
			} else {
				SstartYear = min(year)
				SendYear = max(year)
			}

			forUrl = expand.grid(siteCode=Ssite, sexCode=Ssex, stateCode=Sstate, startYear = SstartYear, endYear=SendYear)
			for(D in c('sex','state','site')){
				forUrl[,D] = names(cancerCodes[[D]])[match(forUrl[,paste(D, 'Code', sep='')], cancerCodes[[D]])]
			}

			
			
			allCases = data.frame()
			
			for(D in 1:nrow(forUrl))	{
				
				if(any(is.na(forUrl[D,]))){
					warning('site, sex or state not found')
				}
				
				if(forUrl[D, 'stateCode'] %in% c('ar')) {
					middleUrl = 'beta/common/v1'
          depth=5
				} else {
					middleUrl = 'common'
					depth=5
				}
				
				kUrl = paste(
						'http://www.cancer-rates.info/',middleUrl, '/index.php?',
						'std=us2000m&geography=1&syear=', forUrl[D, 'startYear'], 
						'&eyear=', forUrl[D, 'endYear'], '&site=', forUrl[D, 'siteCode'], 
						'&race=0&sex=', forUrl[D, 'sexCode'],	'&dataset=I&database=', forUrl[D, 'stateCode'], 
						'&datasource=inv&m_color=1&c_intrvls=0&title=Stuff&r=370,*&c=538', 
						sep='')
				
				namesHere = forUrl[D, c('startYear','endYear','sex','state','site')]
				
				myDir = file.path(tempdir(), gsub("[[:space:]]+", "_", paste(c("cfiles",as.character(namesHere)), collapse='_')))
				
				myCommand = paste("httrack --depth=",depth, " --priority=1 -N1 -O ", myDir, " \'",  kUrl,  "\'", sep='')
				cat('\ndownloading', paste(namesHere, collapse=" "))
				sRes = try(system(myCommand))
				cat(' done\n')
				if(class(sRes)=='try-error') {
					stop("install httrack from www.httrack.com")
				}
				

				fname = system(paste("ls ", myDir, "/web/newalldetails*.html", sep=''), TRUE)
				
				
				fname = grep('newalldetails[[:alnum:]]+\\.html$', fname, value=TRUE)

				datHeader = XML::readHTMLTable(fname[1], isUrl=FALSE, which=2)[1,]
				
				datText = scan(fname[1], what=character(), quiet=TRUE)
				
				startTable = grep("<TABLE", datText)
				startTable = startTable[length(startTable)]

				endTable = grep("</TABLE", datText)
				endTable = endTable[length(endTable)]
				
				datText = datText[startTable:endTable]
				
				startTr = grep("<TR", datText)
				if(min(startTr)>1)
					datText = datText[-seq(1,min(startTr)-1)]
				startTr = grep("<TR", datText)
				startTd = grep("<TD", datText)
				
				earlyTr = max(which(startTr <= startTd[1]))
				startTr = startTr[seq(earlyTr, length(startTr))]
				
				datText = paste(datText, collapse='')
				datText = unlist(strsplit(datText, "<TR>"))
				datText = grep("^[[:space:]]?$|^<TH", datText, invert=TRUE, value=TRUE)

				datSplit = strsplit(datText, "<TD")
				datLen = unlist(lapply(datSplit, length))
				datSplit = unlist(datSplit[which(datLen ==6)])
				datSplit = gsub("</A>$", "", datSplit)
				datSplit = gsub("^[[:print:]]+>", "", datSplit)
				datSplit = gsub("~", "1", datSplit)
				
				dat  = as.data.frame(matrix(datSplit, ncol=6, byrow=TRUE),
						stringsAsFactors=FALSE)[,-1]
				colnames(dat) = as.character(unlist(datHeader))
				dat = dat[grep("^$", dat$County, invert=TRUE),]
				for(Dcol in grep("County", colnames(dat), invert=TRUE)) {
					dat[,Dcol] = as.numeric(as.character(dat[,Dcol]))
				}
				
				cases = dat[,c('County','Cases')]
				
				if(nrow(cases)){
					cases = cbind(cases, namesHere, stateCode=forUrl[D,'stateCode'])
					allCases = rbind(allCases, cases)
				}
			}		
			
			getRid = grep("^[[:space:]]?(Note|Data|~|The population estimates):?[[:space:]]", allCases$County)
			getRidState = grep("^[[:space:]]?(STATE|\\*\\*\\*Counts|Unknown)[[:space:]]?$", allCases$County)
			getRid = c(getRid, getRidState)
			if(length(getRid))
				allCases = allCases[-getRid,]
			
			allCases$Cases = as.numeric(gsub("~", "1", allCases$Cases))
	
			allCases
			
		}