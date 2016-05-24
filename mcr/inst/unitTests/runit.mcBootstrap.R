cat("\n\nmcBootstrap.R method comparison test cases\n\n")


test.mcBootstrap.call <- function() 
{
    data(creatinine)
    set.seed(19061978)                          
    index <- sample(1:nrow(creatinine), size=30)
    crea.data <- creatinine[index,]                  # no NAs for this choice of seed
    crea.data <- as.matrix(crea.data)
    
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife="non-logical") )           # check parameter jacknife
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=as.logical(NA)) )
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=logical(0)) )
    
    checkException( mcr:::mc.bootstrap(method.reg="something", X=crea.data[,1], Y=crea.data[,2]) )                                    # wrong method
    
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[1:10,1], Y=crea.data[1:20,2]) )                               # check data
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=numeric(0), Y=numeric(0))  ) 
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=LETTERS, Y=LETTERS)  )
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=as.numeric(NA), Y=as.numeric(NA))  )
    checkException( mcr:::mc.bootstrap(method.reg="LinReg", X=NA, Y=NA)  )
    
    checkException( mcr:::mc.bootstrap(method.reg="Deming", X=crea.data[,1], Y=crea.data[,2], error.ratio=numeric(0)) )               # check parameter 'error.ratio'
    checkException( mcr:::mc.bootstrap(method.reg="Deming", X=crea.data[,1], Y=crea.data[,2], error.ratio=as.numeric(NA)) ) 
    checkException( mcr:::mc.bootstrap(method.reg="Deming", X=crea.data[,1], Y=crea.data[,2], error.ratio=-1) )
    
    checkException( mcr:::mc.bootstrap(method.reg="WDeming", X=crea.data[,1], Y=crea.data[,2], error.ratio=1, iter.max=NA) )          # check parameter 'iter.max'
    checkException( mcr:::mc.bootstrap(method.reg="WDeming", X=crea.data[,1], Y=crea.data[,2], error.ratio=1, iter.max=numeric(0)) )
    checkException( mcr:::mc.bootstrap(method.reg="WDeming", X=crea.data[,1], Y=crea.data[,2], error.ratio=1, iter.max=0) )
    
    obj.LR.jk <- list(  glob.coef=c(-0.1141941, 1.0796718), glob.sigma=c(0.06192228, 0.04633909), xmean=1.236667,
                        B0jack=c(-0.1094335, -0.1153060, -0.1076864, -0.1126950, -0.1036905, -0.1073347, -0.1175696, -0.1122742, -0.1107537, -0.1155753, -0.1076899, -0.1156347, -0.1089353, -0.1154491, -0.1089383,
                                 -0.1197399, -0.1163754, -0.1238692, -0.1270153, -0.1064683, -0.1220275, -0.1289707, -0.1400572, -0.1153418, -0.1139551, -0.1140733, -0.1132589, -0.1057711, -0.1157885, -0.1018171),
                        B1jack=c( 1.077704, 1.083021, 1.077581, 1.083118, 1.073545, 1.066255, 1.083320, 1.078590, 1.078081, 1.076496, 1.070464, 1.080338, 1.078418, 1.080233, 1.072720, 1.082063, 1.080776, 1.084568,
                                  1.099984, 1.076408, 1.082189, 1.093007, 1.093287, 1.080246, 1.080221, 1.079607, 1.076855, 1.075238, 1.080052, 1.072741),
                        npoints=30, cimeth="jackknife", weight=rep(1,dim(crea.data)[1])) 
    
    checkEquals( mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=TRUE ), obj.LR.jk, tolerance=10e-7 )
    
    set.seed(42) 
    smpls <- sample(1:100, 10)
    
    # Linear Regression
   
    set.seed(42)
    obj.mc.bs.lr <- mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="bootstrap", nsamples=100)
    checkEquals(obj.mc.bs.lr$glob.coef, c(-0.114194075011703, 1.079671758769571))
    checkEquals(obj.mc.bs.lr$B0[smpls], c(-0.149442172003223, -0.10009824901689, -0.182070795157603, -0.00963083414263788, -0.0417436651764671, -0.0898870402490386, -0.0732871321216504, -0.116277937500499, -0.0956027811456019, -0.174617127982557))
    checkEquals(obj.mc.bs.lr$sigmaB1[smpls], c(0.046706203428574, 0.0532798572886871, 0.0563029333649874, 0.0780892665640541, 0.090346807359435, 0.0480789433811427, 0.0648781950690985, 0.0498024672564715, 0.0696705739759542, 0.0443377393718819))
    checkEquals(obj.mc.bs.lr$cimeth, "bootstrap")
    

    set.seed(42)
    obj.mc.nbs.lr <- mcr:::mc.bootstrap(method.reg="LinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="nestedbootstrap", nsamples=100, nnested=10)
    checkEquals(obj.mc.nbs.lr$glob.coef, c(-0.114194075011703, 1.079671758769571), tolerance=10e-7)
    checkEquals(obj.mc.nbs.lr$B0[smpls], c(0.01837780, -0.12765121, -0.13129154, -0.09558644, -0.14156246, -0.10937373, -0.19590455, -0.07164572, -0.16956044, -0.08179280), tolerance=10e-7)
    checkEquals(obj.mc.nbs.lr$sigmaB1[smpls], c(0.13387268, 0.05228267, 0.05913064, 0.02716473, 0.02139479, 0.03588489, 0.08150116, 0.01431682, 0.01891546, 0.06971891), tolerance=10e-7)
    checkEquals(obj.mc.nbs.lr$cimeth, "nestedbootstrap")
    
    # Weighted Linear Regression
    
    set.seed(42)
    obj.mc.bs.wlr <- mcr:::mc.bootstrap(method.reg="WLinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="bootstrap", nsamples=100)
    checkEquals(obj.mc.bs.wlr$glob.coef, c(-0.13798868997394, 1.10026004389133))
    checkEquals(obj.mc.bs.wlr$B0[smpls], c(-0.209591364393151, -0.0224109937189423, -0.214064995126739, -0.021265410135924, -0.0786128865381872, -0.139363303793797, -0.061884485304752, -0.0770094089388964, -0.109362300760702, -0.228145123255878))
    checkEquals(obj.mc.bs.wlr$sigmaB1[smpls], c(0.0598680111823546, 0.0607445012399577, 0.0676523402815724, 0.0847176875735508, 0.096281244472021, 0.0724625952830802, 0.0819794075846059, 0.0770717568618718, 0.0675026489278241, 0.0509761319533487))
    checkEquals(obj.mc.bs.wlr$cimeth, "bootstrap")
    
    set.seed(42)
    obj.mc.nbs.wlr <- mcr:::mc.bootstrap(method.reg="WLinReg", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="nestedbootstrap", nsamples=100, nnested=10)
    checkEquals(obj.mc.nbs.wlr$glob.coef, c(-0.137988689973938, 1.100260043891328))
    checkEquals(obj.mc.nbs.wlr$B0[smpls], c(-0.0220899939152511, -0.225168625246526, -0.175640840045823, -0.12664436178391, -0.128315905320232, -0.106298039742549, -0.185812712172651, -0.0974580947884642, -0.16749998224116, -0.0731344061899082))
    checkEquals(obj.mc.nbs.wlr$sigmaB1[smpls], c(0.113046909404017, 0.0545211760156853, 0.053243797863501, 0.0361216906153157, 0.0269697448867467, 0.0547747779061429, 0.0671336434784934, 0.0312002802965443, 0.0300078319457522, 0.0915308544221075))
    checkEquals(obj.mc.nbs.wlr$cimeth, "nestedbootstrap")
    
    # Deming
    
    set.seed(42)
    obj.mc.bs.dem <- mcr:::mc.bootstrap(method.reg="Deming", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="bootstrap", nsamples=100, error.ratio=1)
    checkEquals(obj.mc.bs.dem$glob.coef, c(-0.151736318860087, 1.110029368350475))
    checkEquals(obj.mc.bs.dem$B0[smpls], c(-0.187352238567723, -0.155342047907383, -0.23478551494312, -0.114947501264896, -0.179151541150708, -0.128415240132832, -0.149168586537764, -0.160808205527603, -0.181989271732016, -0.209758348912587))
    checkEquals(obj.mc.bs.dem$sigmaB1[smpls], c(0.0480029520696056, 0.0553851325960843, 0.0584478441138261, 0.0851996608660044, 0.10094459923718, 0.0495465121033404, 0.068539489782525, 0.0514384209129926, 0.0742177058410856, 0.0454368247765725))
    checkEquals(obj.mc.bs.dem$cimeth, "bootstrap")
    
    set.seed(42)
    obj.mc.nbs.dem <- mcr:::mc.bootstrap(method.reg="Deming", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="nestedbootstrap", nsamples=100, nnested=10, error.ratio=1)
    checkEquals(obj.mc.nbs.dem$glob.coef, c(-0.151736318860087, 1.110029368350475))
    checkEquals(obj.mc.nbs.dem$B0[smpls], c(-0.100353687441834, -0.164165155770712, -0.14850823154574, -0.112283018542142, -0.173858501673744, -0.147970291542156, -0.262749133692674, -0.0922917602426494, -0.192178619749573, -0.181815862124773))
    checkEquals(obj.mc.nbs.dem$sigmaB1[smpls], c(0.155663347996796, 0.0560480385829639, 0.0873466220333589, 0.0357497458369248, 0.0115681703850787, 0.0447257168543137, 0.0852965111592678, 0.0183325120557489, 0.017983210520142, 0.0793436978380041))
    checkEquals(obj.mc.nbs.dem$cimeth, "nestedbootstrap")
    
    # Weighted Deming
    
    set.seed(42)
    obj.mc.bs.wdem <- mcr:::mc.bootstrap(method.reg="WDeming", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="bootstrap", nsamples=100, error.ratio=1)
    checkEquals(obj.mc.bs.wdem$glob.coef, c(-0.20833788085213, 1.15884355708208))
    checkEquals(obj.mc.bs.wdem$B0[smpls], c(-0.253732510159959, -0.0787643881048334, -0.288529356609159, -0.121138866165996, -0.209391448955941, -0.208048795145741, -0.178919028906314, -0.176708853713483, -0.188818656463374, -0.257541374176516))
    checkEquals(obj.mc.bs.wdem$MX[smpls], c(1.01611323598916, 1.09023823698536, 1.00308992400947, 1.03453640915801, 0.988881922604329, 0.970960828197418, 1.0309963242759, 1.06268271958762, 1.07317934400164, 0.989284432999055))
    checkEquals(obj.mc.bs.wdem$cimeth, "bootstrap")
    
    set.seed(42)
    obj.mc.nbs.wdem <- mcr:::mc.bootstrap(method.reg="WDeming", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="nestedbootstrap", nsamples=100, nnested=10, error.ratio=1)
    checkEquals(obj.mc.nbs.wdem$glob.coef, c( -0.20833788085213, 1.15884355708208))
    checkEquals(obj.mc.nbs.wdem$B0[smpls], c(-0.147843778652592, -0.279903535435754, -0.242473261007011, -0.19386806149172, -0.169351800694116, -0.191662942949458, -0.238951513557179, -0.149145296900057, -0.223964100998818, -0.184997804213853))
    checkEquals(obj.mc.nbs.wdem$sigmaB1[smpls], c(0.11183633771502, 0.0494769836857074, 0.0752383830743171, 0.040083857908686, 0.0236556930999204, 0.065690363014448, 0.0825827303361371, 0.0276343974159533, 0.0301389475800413, 0.0790736204920856))
    checkEquals(obj.mc.nbs.wdem$cimeth, "nestedbootstrap")
    
    # Passing-Bablok
    
    set.seed(42)
    obj.mc.bs.pb <- mcr:::mc.bootstrap(method.reg="PaBa", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="bootstrap", nsamples=100, error.ratio=1)
    checkEquals(obj.mc.bs.pb$glob.coef, c(-0.184623847841471, 1.157751706128444))
    checkEquals(obj.mc.bs.pb$B0[smpls], c(-0.19373417721519, -0.127139984000515, -0.246613402933132, -0.163472457210465, -0.22, -0.198636363636364, -0.145384615384615, -0.167272727272727, -0.161269841269841, -0.251594202898551))
    checkEquals(obj.mc.bs.pb$MX[smpls], c(1.23266666666667, 1.36466666666667, 1.18366666666667, 1.18133333333333, 1.14133333333333, 1.177, 1.25833333333333, 1.26633333333333, 1.23766666666667, 1.27066666666667))
    checkEquals(obj.mc.bs.pb$cimeth, "bootstrap")
    
    set.seed(42)
    obj.mc.nbs.pb <- mcr:::mc.bootstrap(method.reg="PaBa", X=crea.data[,1], Y=crea.data[,2], jackknife=FALSE, bootstrap="nestedbootstrap", nsamples=100, nnested=10, error.ratio=1)
    checkEquals(obj.mc.nbs.pb$glob.coef, c(-0.184623847841471, 1.157751706128444))
    checkEquals(obj.mc.nbs.pb$B0[smpls], c(-0.140882352941176, -0.193951612903226, -0.209460217580462, -0.134941520467837, -0.145833333333333, -0.172857142857143, -0.271683168316832, -0.134941520467837, -0.220769230769231, -0.171578947368421))
    checkEquals(obj.mc.nbs.pb$sigmaB1[smpls], c(0.0962574390955512, 0.0713521702605754, 0.054047154373797, 0.034149540914331, 0.0194428549936299, 0.0621111889553757, 0.0853939447603508, 0.051585238550563, 0.038646889179085, 0.0692536020998787))
    checkEquals(obj.mc.nbs.pb$cimeth, "nestedbootstrap")
}
    
    