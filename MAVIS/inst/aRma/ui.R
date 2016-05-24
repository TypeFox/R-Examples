library(shiny)
library(shinyAce)
library(shinyBS)
library(meta)
library(metafor)
#library(metamisc)
library(MAd)
library(MAc)
library(quantreg)
library(ggplot2)
library(compute.es)
library(SCMA)
library(SCRT)
shinyUI(navbarPage("aRma:Altyapısı R'a dayalı Meta Analiz  ",
                   
                   tabPanel("Giriş", icon = icon("floppy-save", lib = "glyphicon"),
                            wellPanel(width = 5,
                              
                              radioButtons("type", strong("Analiz ve Veri Giriş Seçenekleri:"),
                                           list("Ham Ortalamaların Farkı (n, M, SD)" = "mdms",
                                                "Standartlaştırılmış Ortalamaların Farkı (n, Etki büyüklüğü d)" = "mdes",
                                                "Korelasyon (n, r)" = "cor",
                                                "İkili veriler"="or"
                                           ),
                              ),
                              bsTooltip("type", "Analizlerini koşmadan önce araç çubuğunda yer alan seçenekleri kontrol ediniz.",
                                        "right", trigger = "click", options = list(container = "body")),
                              helpText("İkili veriler seçeneğini kontrol ediniz, varsayılan seçenek logaritmik olasılık oranıdır (log odds ratio)"),
                              checkboxInput("moderator", label = ("Altgrup verisi mevcut ise işaretleyiniz."), value = T),
                              br(),
                              submitButton("Yenile"),
                              helpText("Eğer datayı, modeli yada seçenekleri değiştirdiyseniz sonuçları yenileyiniz"),
                              br(),
                              actionButton("quit", "Çıkış", icon("sign-out")),
                              helpText("Uygulamadan çıkış")
                              
                              
                            ),
                            br(),
                            
                            p('Not: Sütünlar birbirinden tab ile ayrılmalıdır.Daha kolay kullanım için lütfen verileri excel dosyası olarak hazırlayın kopyalayın ve yapıştırın.'),
                            p("Excel veri dosyasının ilk satırı (başlıkları) örnek verilerin başlıkları ile birebir aynı olmalıdır."),
                            p("Örnek veri girişi için lütfen Örnek Veriler sekmesine gidiniz"),
                            
                            aceEditor("text", value="Veri\tModerator\tN1\tM1\tSD1\tN2\tM2\tSD2\nVeri 01\tJH\t30\t51.57\t16.5\t30\t72.97\t13.23\nVeri 02\tUNI\t23\t75.09\t23.01\t24\t81.63\t14.42\nVeri 03\tSH\t83\t30.08\t14.29\t81\t35.38\t16.13\nVeri 04\tSH\t21\t2.95\t1.28\t21\t3.48\t0.68\nVeri 05\tSH\t15\t53.8\t17.4\t15\t60.47\t17.37\nVeri 06\tSH\t7\t15.7\t4.1\t7\t27.3\t4.1\nVeri 07\tSH\t28\t27.9\t9.57\t28\t33.2\t15.65\nVeri 08\tUNI\t40\t17.53\t8.87\t40\t19.23\t9.55\nVeri 09\tUNI\t18\t11.86\t13.24\t17\t29.92\t16.67\nVeri 10\tUNI\t21\t29.76\t16\t25\t27.98\t16.52\nVeri 11\tUNI\t26\t8.23\t3.59\t26\t9.65\t2.99\nVeri 12\tUNI\t49\t13.71\t4.07\t48\t16\t3.47\nVeri 13\tUNI\t27\t2.8\t1.7\t27\t5.9\t1.4\nVeri 14\tSH\t41\t10.05\t2.52\t34\t11.03\t1.78\nVeri 15\tUNI\t58\t3.62\t1.79\t57\t4.26\t1.61\nVeri 16\tSH\t60\t7.36\t2.8\t63\t8.82\t2.5\nVeri 17\tUNI\t15\t5.93\t3.55\t15\t12.27\t4.95\nVeri 18\tJH\t37\t13.68\t3.68\t142\t17.53\t4.34\nVeri 19\tJH\t27\t3.3\t2.3\t54\t12.98\t7.67\nVeri 20\tJH\t35\t5.49\t3.88\t39\t12.36\t7.68\nVeri 21\tJH\t32\t5.81\t3.14\t34\t12.44\t5.66\nVeri 22\tJH\t62\t17.84\t4.09\t60\t18.18\t4.09\nVeri 23\tSH\t39\t8.77\t5\t39\t13.72\t5.32\nVeri 24\tSH\t213\t59.8\t15.3\t39\t79.8\t9.5\nVeri 25\tUNI\t34\t14.32\t2.79\t42\t16\t2.05\nVeri 26\tUNI\t77\t70.85\t11.74\t56\t78.17\t9.94\nVeri 27\tUNI\t28\t80.83\t22.47\t28\t85.06\t23.52\nVeri 28\tUNI\t33\t25.38\t4.71\t36\t25.02\t3.36\nVeri 29\tUNI\t66\t0.45\t0.29\t66\t0.93\t0.59",
                                      mode="r", theme="mono_industrial"),
                            
                            br(),
                            
                            h3("Etki Büyüklüğü ve Örneklem Varyansı (iç varyans)"),
                            
                            verbatimTextOutput("data.out"),
                            
                            br(),
                            
                            h3("Sabit Etki Modeli (FE Model)"),
                            verbatimTextOutput("fe.out"),
                            
                            br(),
                            
                            h3("Rastgele Etkiler Modeli"),
                            verbatimTextOutput("re.out"),
                            
                            p('[Heterojenliği Tanımlama Çabaları]',
                              br(),
                              br(),
                              'I^2 (Etki büyüklüğü çalışmalar arasında ne kadar değişiyor?)', br(),
                              "25-50: Az Değisken", br(),
                              '50-75: Değisken', br(),
                              '75-100: Çok Değişken', br(),
                              br(),
                              'Heterojenlik testi: p-val < .05 (homojen degil)', br(),
                              br(),
                              'H > 1: Açıklanmamış heterojenlik var.'
                            ),
                            
                            br(),
                            br(),
                            
                            
                            h3("Diyagram (Sabit Etki Modeli)"),
                            downloadButton('downloadfePlot', 'Diyagrami pdf olarak indir'),
                            plotOutput("fePlot", height = "550px"),
                            
                            
                            
                            br(),
                            
                            h3("Diyagram (Rastgele Etkiler Modeli)"),
                            downloadButton('downloadrePlot', 'Diyagrami pdf olarak indir'),
                            plotOutput("rePlot", height = "550px"),
                            
                            br(),
                            
                            h3("Huni grafiği (Sabit Etki Modeli)"),
                            downloadButton('downloadFunFixPlot', 'Grafiği pdf olarak indir'),
                            plotOutput("FunFixPlot"),
                            p('Eğer içi beyaz olan daireler varsa bunlar grafiğe kırpma ve doldurma metodu ile eklenmiştir.'),
                            
                            br(),
                            
                            h3("Huni grafigi (Rastgele Etkiler Modeli)"),
                            downloadButton('downloadFunRandPlot', 'Grafiği pdf olarak indir'),
                            plotOutput("FunRandPlot"),
                            p('Eğer içi beyaz olan daireler varsa bunlar grafiğe kırpma ve doldurma metodu ile eklenmiştir.'),
                            br(),
                            br(),
                            h3("Yayın Yanlılığı"),
                            verbatimTextOutput("asy.out"), # regression tests for funnel plot asymmetry
                            p('Korumalı N: Bir meta-analiz çalışması yapılsın ve bulunan p değeri istatistiksel olarak anlamlı olsun. Korumalı N , çıkan anlamlı sonucun (p<.05), istatistiksel olarak anlamsızlaşması için (p>.05) kaç adet etkisiz çalışmanın eklenmesini gerektiğini verir. Eğer bu sayı büyük ise , bulunan anlamlı farkın yayın yanlılığına karşı dirençli olduğu söylenebilir. Daha detali yorumlar için bknz: ',
                              a('(Oswald & Plonsky, 2010, p. 92)', href='http://dx.doi.org/10.1017/S0267190510000115', target="_blank"), '.'),
                            br(),
                            
                            br(),
                            
                            # Display this only if "moderator" is checked
                            conditionalPanel(condition = "input.moderator == true",
                                             h3("Alt grup analizleri"),
                                             verbatimTextOutput("modAnalysis.out")
                            ),
                            
                            br(),
                            
                            # Display this only if "moderator" is checked
                            conditionalPanel(condition = "input.moderator == true",
                                             h4("Moderator grafiği (Sabit Etki Modeli)"),
                                             plotOutput("ModFixGraph")
                            ),
                            
                            br(),
                            
                            # Display this only if "moderator" is checked
                            conditionalPanel(condition = "input.moderator == true",
                                             h4("Moderator Grafiği(Rasal Etki Modeli)"),
                                             plotOutput("ModRandGraph")
                            ),
                            
                            br(),
                            br(),
                            
                            strong('Oturum Bilginiz'),
                            verbatimTextOutput("info.out")
                   ),
                   
                   
                   
                   tabPanel("Örnek Veriler", icon = icon("table", lib = "font-awesome"),
                            
                            p('Not: Sutunlar birbirinden tab ile ayrılmalıdır . Daha kolay kullanım için lütfen verileri excel dosyası olarak hazırlayın kopyalayın ve yapıştırın'),
                            
                            p(HTML("<b><div style='background-color:#FADDF2;border:1px solid black;'>Excel veri dosyasının ilk satırı (başlıkları) ornek verilerin başlıkları ile birebir aynı olmalıdır.</div></b>")),
                                                        
                            br(),
                            
                            p(strong("Ham Ortalamalar (n, M, SD)")),
                            aceEditor("text1", value="Veri\tModerator\tN1\tM1\tSD1\tN2\tM2\tSD2\nVeri 01\tJH\t30\t51.57\t16.5\t30\t72.97\t13.23\nVeri 02\tUNI\t23\t75.09\t23.01\t24\t81.63\t14.42\nVeri 03\tSH\t83\t30.08\t14.29\t81\t35.38\t16.13\nVeri 04\tSH\t21\t2.95\t1.28\t21\t3.48\t0.68\nVeri 05\tSH\t15\t53.8\t17.4\t15\t60.47\t17.37\nVeri 06\tSH\t7\t15.7\t4.1\t7\t27.3\t4.1\nVeri 07\tSH\t28\t27.9\t9.57\t28\t33.2\t15.65\nVeri 08\tUNI\t40\t17.53\t8.87\t40\t19.23\t9.55\nVeri 09\tUNI\t18\t11.86\t13.24\t17\t29.92\t16.67\nVeri 10\tUNI\t21\t29.76\t16\t25\t27.98\t16.52\nVeri 11\tUNI\t26\t8.23\t3.59\t26\t9.65\t2.99\nVeri 12\tUNI\t49\t13.71\t4.07\t48\t16\t3.47\nVeri 13\tUNI\t27\t2.8\t1.7\t27\t5.9\t1.4\nVeri 14\tSH\t41\t10.05\t2.52\t34\t11.03\t1.78\nVeri 15\tUNI\t58\t3.62\t1.79\t57\t4.26\t1.61\nVeri 16\tSH\t60\t7.36\t2.8\t63\t8.82\t2.5\nVeri 17\tUNI\t15\t5.93\t3.55\t15\t12.27\t4.95\nVeri 18\tJH\t37\t13.68\t3.68\t142\t17.53\t4.34\nVeri 19\tJH\t27\t3.3\t2.3\t54\t12.98\t7.67\nVeri 20\tJH\t35\t5.49\t3.88\t39\t12.36\t7.68\nVeri 21\tJH\t32\t5.81\t3.14\t34\t12.44\t5.66\nVeri 22\tJH\t62\t17.84\t4.09\t60\t18.18\t4.09\nVeri 23\tSH\t39\t8.77\t5\t39\t13.72\t5.32\nVeri 24\tSH\t213\t59.8\t15.3\t39\t79.8\t9.5\nVeri 25\tUNI\t34\t14.32\t2.79\t42\t16\t2.05\nVeri 26\tUNI\t77\t70.85\t11.74\t56\t78.17\t9.94\nVeri 27\tUNI\t28\t80.83\t22.47\t28\t85.06\t23.52\nVeri 28\tUNI\t33\t25.38\t4.71\t36\t25.02\t3.36\nVeri 29\tUNI\t66\t0.45\t0.29\t66\t0.93\t0.59", mode="r", theme="monokai"),
                            
                            
                            br(),
                            p(strong("Standardlaştırılmış Ortalamalar  (n, Effect size d)")),
                            aceEditor("text2", value="Veri\tModerator\tN1\tN2\td\nVeri 01\tJH\t30\t30\t-1.431\nVeri 02\tUNI\t23\t24\t-0.3423\nVeri 03\tSH\t83\t81\t-0.3481\nVeri 04\tSH\t21\t21\t-0.5171\nVeri 05\tSH\t15\t15\t-0.3837\nVeri 06\tSH\t7\t7\t-2.8293\nVeri 07\tSH\t28\t28\t-0.4086\nVeri 08\tUNI\t40\t40\t-0.1845\nVeri 09\tUNI\t18\t17\t-1.2039\nVeri 10\tUNI\t21\t25\t0.1093\nVeri 11\tUNI\t26\t26\t-0.4298\nVeri 12\tUNI\t49\t48\t-0.605\nVeri 13\tUNI\t27\t27\t-1.9907\nVeri 14\tSH\t41\t34\t-0.4422\nVeri 15\tUNI\t58\t57\t-0.3758\nVeri 16\tSH\t60\t63\t-0.5508\nVeri 17\tUNI\t15\t15\t-1.4719\nVeri 18\tJH\t37\t142\t-0.9136\nVeri 19\tJH\t27\t54\t-1.5079\nVeri 20\tJH\t35\t39\t-1.111\nVeri 21\tJH\t32\t34\t-1.4368\nVeri 22\tJH\t62\t60\t-0.0831\nVeri 23\tSH\t39\t39\t-0.9588\nVeri 24\tSH\t213\t39\t-1.3729\nVeri 25\tUNI\t34\t42\t-0.6976\nVeri 26\tUNI\t77\t56\t-0.6642\nVeri 27\tUNI\t28\t28\t-0.1839\nVeri 28\tUNI\t33\t36\t0.0886\nVeri 29\tUNI\t66\t66\t-1.0326", mode="r", theme="monokai"),
                            
                            
                            br(),
                            p(strong("Korelasyonlar (n, r)")),
                            aceEditor("text3", value="Veri\tN\tr\tModerator\nIzumi (2000)\t175\t0.78\tcollege\nYu (2009)\t53\t0.38\tJS high\nThuy (1996)\t250\t0.69\tcollege\nOckey (2002)\t90\t0.89\tcollege\nAraru (2005)\t86\t0.52\tJS high\nWee (1997)\t182\t0.98\tcollege\nOzoda (2007)\t591\t0.91\tcollege\nHala (2004)\t30\t0.95\tcollege\nTapio (2008)\t37\t0.47\tJS high\nAndarani (2008)\t107\t0.84\tcollege\nDavis (1999)\t74\t0.99\tcollege\nPlonsky (2002)\t217\t0.86\tcollege\nGassel (1993)\t203\t0.99\tcollege",mode="r", theme="monokai"),
                            
                            br(),
                            p(strong("İkili veriler (olay görülme, olay görülmeme, n)")),
                            aceEditor("text4", value="Veri\tModerator\tupoz\tuneg\tNU\tkpoz\tkneg\tNK\nVeri 01\tra\t4\t119\t123\t11\t128\t139\nVeri 02\tra\t6\t300\t306\t29\t274\t303\nVeri 03\tra\t3\t228\t331\t11\t209\t220\nVeri 04\tsys\t17\t1699\t1716\t65\t1600\t1665\nVeri 05\tsys\t5\t2493\t2498\t3\t2338\t2341\nVeri 06\tra\t29\t7470\t7499\t45\t7232\t7277", mode="r", theme="monokai"),
                            
                            br()
                            
                   ),
                   
                   navbarMenu("Model seçenekleri ve ayarları ", icon = icon("cog", lib = "font-awesome"),
                              #                       tabPanel("Bayesian Model Options", icon = icon("tasks", lib = "font-awesome"),
                              #                                
                              #                      strong('Bayesian Analysis Options'),
                              #                      selectInput("bayoption1", label = h3("Run Bayesian Analysis"), 
                              #                                  choices = list("No" = FALSE, "Yes" = TRUE)
                              #                                
                              #                       )),
                              tabPanel("Korelasyon modeli seçenekleri", icon = icon("line-chart", lib = "font-awesome"),
                                       
                                       radioButtons("cormeasures", strong("Korelasyon modeli seçenekleri"),
                                                    c("Ham korelasyon katsayısı " = "COR",
                                                      "Düzeltilmiş ham korelasyon, Olkin & Pratt, 1958" = "UCOR",
                                                      "Fisher Z dönüşümlü korelasyon" = "ZCOR"
                                                    ), selected = "ZCOR"),
                                       p(h6('Aksi belirtilmedikçe, sisteme sunulan korelasyon verisi Fisher Z dönüşümlü korelasyon olarak işlem görür.')),
                                       
                                       verbatimTextOutput('cormeasures.out')
                                       
                              ),
                              tabPanel("İkili veri seçenekleri ", icon = icon("ellipsis-v", lib = "font-awesome"),
                                       
                                       radioButtons("dichotomousoptions", strong("İndeks seçimi"),
                                                    c("log risk oranı" = "RR",
                                                      "log olasılık oranı" = "OR",
                                                      "Risk farkı" = "RD",
                                                      "arcsine karekök transorfmlu risk farkı (Rücker et al., 2009)." = "AS",
                                                      "Peto metodu ile tahminlenmiş olasılık oranı (Yusuf et al., 1985)." = "PETO",
                                                      "Standardlaştırılmış ortalama farkı tahmini olarak probit transformlu risk farkı." = "PBIT",
                                                      "Standardlaştırılmış ortalama farkı tahmini olarak transforme olasılık oranı (normal dağılım)." = "OR2DN",
                                                      "Standardlaştırılmış ortalama farkı tahmini olarak transforme olasılık oranı  (logistic dağılım)." = "OR2DL"
                                                    ), selected = "OR"),
                                       p(h6('Aksi belirtilmediğinde, analizlerin yapıldığı R paketi log olasılık oranı kullanır')),
                                       
                                       verbatimTextOutput('dichotomousoptions.out')
                                       
                              ),
                              tabPanel("Rassal etki modeli tahminleyici seçenekleri", icon = icon("random", lib = "glyphicon"),
                                       
                                       radioButtons("model", strong("Tahminleyici seçenekleri"),
                                                    c("DerSimonian-Laird" = "DL",
                                                      "Hedges" = "HE",
                                                      "Hunter-Schmidt" = "HS",
                                                      "Sidik-Jonkman" = "SJ",
                                                      "Maximum-likelihood" = "ML",
                                                      "Restricted maximum-likelihood" = "REML",
                                                      "Empirical Bayes (Paule-Mandel)" = "EB",
                                                      "Generalized Q-statistic" = "GENQ"
                                                    ), selected = "REML"),
                                       p(h6('Aksi belirtilmediğinde, analizlerin yapıldığı R paketi REML kullanır')),
                                       
                                       
                                       checkboxInput("khadjust", label = "Knapp & Hartung Düzeltmesi", value = FALSE),
                                       p(h6('Aksi belirtilmedikçe Knapp & Hartung düzeltmesi yapılmaz')),
                                       p("The Knapp and Hartung (2003) metodunun standard hataları düzeltmesi beklenir."),
                                       h3("References"),
                                       p("Knapp, G., & Hartung, J. (2003). Improved tests for a random effects meta-regression with a single covariate. Statistics in Medicine, 22, 2693–2710.")
                                       
                                       )),
                   
                   navbarMenu("Yayın Yanlılığı ", icon = icon("book", lib = "font-awesome"),
                              tabPanel("Kırpma ve Doldurma", icon = icon("chevron-right", lib = "font-awesome"),
                                       
                                       radioButtons("trimfillopt", strong("Kırpma ve Doldurma Tahminleyiceleri"),
                                                    c("L0" = "L0",
                                                      "R0" = "R0",
                                                      "Q0" = "Q0"
                                                    ), selected = "L0"),
                                       p(h6('Bu üç farklı tahminleyici için bak. Duval ve Tweedie (2000a, 2000b). Aksi belirtilmediyse, L0. ')),
                                       
                                       verbatimTextOutput('trimfillopt.out'),
                                       h3("References"),
                                       p("Duval, S. J., & Tweedie, R. L. (2000a). Trim and fill: A simple funnel-plot-based method of testing and adjusting for publication bias in meta-analysis. Biometrics, 56, 455–463."),
                                       p("Duval, S. J., & Tweedie, R. L. (2000b). A nonparametric trim and fill method of accounting for publication bias in meta-analysis. Journal of the American Statistical Association, 95, 89–98."),
                                       p("Duval, S. J. (2005). The trim and fill method. In H. R. Rothstein, A. J. Sutton, & M. Borenstein (Eds.) Publication bias in meta-analysis: Prevention, assessment, and adjustments (pp. 127–144). Chichester, England: Wiley.")
                                       
                                       
                              ),
                              
                              tabPanel("Huni grafiğinin asitmeriliği için regresyon testi", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(4,
                                                    p(strong("Regresyon testi seçenekleri")),
                                                    radioButtons("regtestpredictor", strong("Yordayıcı"),
                                                                 c("standard hata" = "sei",
                                                                   "örneklem varyansı" = "vi",
                                                                   "örneklem sayısı" = "ni",
                                                                   "örneklem sayısının tersi (inverse)" = "ninv"
                                                                 ), selected = "sei"),
                                                    radioButtons("regtestmodeltype", strong("Model Seçimi "),
                                                                 c("Weighted Regression with a Multiplicative Dispersion Term" = "lm",
                                                                   "Meta-analytic Models" = "rma"
                                                                 ), selected = "lm"),
                                                    
                                                    
                                                    p(br())
                                             ),
                                             column(4,
                                                    strong('Huni grafiği seçenekleri'),
                                                    checkboxInput("contourenhancedbox", "Contour enhanced (Dış çizgisi geliştirilmiş) grafik", FALSE),
                                                    helpText("Eğer huni grafiğini dış çizgisi geliştirilmiş (Peters et al., 2008) olarak rapor etmek istiyorsanız işaretleyiniz."),
                                                    checkboxInput("regtestfullmodel", "Regresyon modeli sonuçları", FALSE),
                                                    helpText("Eğer regresyon modelinin bütün sonuçlarını görmek istiyorsanız işaretleyiniz.")
                                                    
                                             )
                                             
                                           )),
                                         p("Yayın yanlılığını tespit etmek amaçlı diğer metodlar için bknz: (Jin, Zhou, & He, 2015)"),
                                         h3("Referaslar"),
                                         p("Egger, M., Davey Smith, G., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. British Medical Journal, 315, 629--634."),
                                         p("Jin, Zhi-Chao, Zhou, Xiao-Hua & He, Jia (2015). Statistical methods for dealing with publication bias in meta-analysis. Statistics in Medicine, 34, 343-360."),
                                         p("Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2008). Contour-enhanced meta-analysis funnel plots help distinguish publication bias from other causes of asymmetry. Journal of Clinical Epidemiology, 61(10), 991–-996."),
                                         p("Sterne, J. A. C., & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: Guidelines on choice of axis. Journal of Clinical Epidemiology, 54(10), 1046--1055."),
                                         br()
                                         
                                       )
                                       
                              ),
                              
                              tabPanel("Eksik çalışma (File Drawer) Analizleri ", icon = icon("chevron-right", lib = "font-awesome"),
                                       
                                       radioButtons("filedraweranalysis", strong("Eksik Çalışma Analizleri"),
                                                    c("Rosenthal" = "Rosenthal",
                                                      "Orwin" = "Orwin",
                                                      "Rosenberg" = "Rosenberg"
                                                    ), selected = "Rosenthal"),
                                       p(h6('Eksik çalışma analizleri için metod seçimi. Aksi belirtilmedikçe Rosental yöntemi kullanılır.')),
                                       p('Rosenthal metodu için hesaplamalar Stouffer yönetimini kullanır (Rosenthal, 1979)'),
                                       p('Orwin (1983).'),
                                       p('Rosenberg (2005)'),
                                       
                                       verbatimTextOutput('filedraweranalysis.out'),
                                       h3("Referans"),
                                       p("Rosenthal, R. (1979). The file drawer problem and tolerance for null results. Psychological Bulletin, 86, 638--641."),
                                       p("Orwin, R. G. (1983). A fail-safe N for effect size in meta-analysis. Journal of Educational Statistics, 8, 157--159."),
                                       p("Rosenberg, M. S. (2005). The file-drawer problem revisited: A general weighted method for calculating fail-safe numbers in meta-analysis. Evolution, 59, 464--468.")
                                       
                                       
                                       )),
                   navbarMenu("Etki Büyüklüğü Hesaplaması", icon = icon("calculator", lib = "font-awesome"),
                              tabPanel("Ortalama, standard sapma(SS) ve örneklem sayısı", icon = icon("chevron-right", lib = "font-awesome"),
                                       
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("Grup 1:")),
                                                    
                                                    numericInput("nx", " Örneklem Sayısı (n)", 21),
                                                    
                                                    numericInput("mx", " Ortalama", 61.33),
                                                    
                                                    numericInput("sdx", " SS", 16.43),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    p(strong("Grup 2:")),
                                                    
                                                    numericInput("ny", "Örneklem Sayısı (n)", 24),
                                                    
                                                    numericInput("my", " Ortalama", 59.79),
                                                    
                                                    numericInput("sdy", " SS", 18.50),
                                                    
                                                    p(br())
                                             ),
                                             column(4,
                                                    strong('Option:'),
                                                    
                                                    
                                                    checkboxInput("varequal", "Varyans eşitliği varsayımı ile t test", FALSE),
                                                    
                                                    
                                                    checkboxInput("vartest", "Varyans eşitliği testi sonuçları", FALSE),
                                                    helpText("Yenilemek için tıklayınız."),
                                                    submitButton("Yenile")
                                             )
                                             
                                           )),
                                         
                                         h3("Verileri kontrol etme"),
                                         tableOutput("values"),
                                         
                                         br(),
                                         
                                         h3("Ortalamaların farkı ve %95 güven aralığı"),
                                         verbatimTextOutput("difference.out"),
                                         
                                         br(),
                                         
                                         h3("t-test"),
                                         verbatimTextOutput("ttest.out"),
                                         h3(""),
                                         verbatimTextOutput("vartest.out"),
                                         
                                         br(),
                                         
                                         h3("Etki büyüklüğü"),
                                         verbatimTextOutput("es.out"),                   
                                         br()
                                         
                                       )
                                       
                              ),
                              tabPanel("ANCOVA F istatistiğinden etki büyüklüğü hesaplama", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("ANCOVA F istatistiğinden etki büyüklüğü")),
                                                    
                                                    numericInput("ancovaf", " ANCOVA nın F değeri", 21),
                                                    
                                                    numericInput("ancovafn1", " Uygulama grubu örneklem sayısı", 50),
                                                    
                                                    numericInput("ancovafn2", " Kontrol grubu örneklem sayısı", 50),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    
                                                    numericInput("anovafcovar", " Yordayıcı-yordanan yada multiple korelasyon", 0.24),
                                                    
                                                    numericInput("anovafcovarnum", " Yordayıcı Sayısı", 3),
                                                    
                                                    numericInput("sdy", " SS", 18.50),
                                                    helpText("Yenilemek için tıklayınız."),
                                                    submitButton("Yenile"),
                                                    p(br())
                                             )
                                             
                                             
                                           )),
                                         
                                         
                                         h3("Etki Büyüklüğü"),
                                         verbatimTextOutput("ancovaf.out"),
                                         p(br())
                                         
                                       )
                                       
                              ),
                              tabPanel("ANCOVA düzeltilmiş (adjusted) ortalamalardan Etki büyüklüğü", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("ANCOVA düzeltilmiş ortalamalardan Etki büyüklüğü")),
                                                    
                                                    numericInput("ancovamean1", " ANCOVA ile düzeltilmiş uygulama grubu ortalaması", 21.7),
                                                    
                                                    numericInput("ancovamean2", " ANCOVA ile düzeltilmiş kontrol grubu ortalaması", 33.5),
                                                    
                                                    numericInput("ancovameansd", " Düzeltilmiş SS", 50),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    
                                                    numericInput("ancovameann1", " Uygulama grubu örneklem sayısı", 142),
                                                    
                                                    numericInput("ancovameann2", " Kontrol grubu örneklem sayısı", 133),
                                                    
                                                    numericInput("ancovameancovar", " Yordayıcı-yordanan yada multiple korelasyon", 0.24),
                                                    
                                                    numericInput("ancovameancovarnumber", " Yordayıcı sayısı", 3),
                                                    
                                                    
                                                    helpText("Yenilemek için tıklayınız"),
                                                    submitButton("Yenile"),
                                                    p(br())
                                             )
                                             
                                             
                                           )),
                                         
                                         
                                         h3("Etki Büyüklüğü"),
                                         verbatimTextOutput("ancovamean.out"),
                                         p(br())
                                         
                                       )
                                       
                              ),
                              tabPanel("Chi-Squared istatistiğinden etki büyüklüğü", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("Chi-Squared istatistiğinden etki büyüklüğü")),
                                                    
                                                    numericInput("chisquaredstat", " Chi-Squared istatistiği.", 5.3),
                                                    
                                                    numericInput("chisquaredn1", " Örneklem Sayısı.", 50),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    helpText("Yenilemek için tıklayınız"),
                                                    submitButton("Yenile"),
                                                    p(br())
                                             )
                                             
                                             
                                           )),
                                         
                                         
                                         h3("Etki Büyüklüğü"),
                                         verbatimTextOutput("chisquared.out"),
                                         p(br())
                                         
                                       )
                                       
                              ),
                              tabPanel("İki grup kıyaslama", icon = icon("chevron-right", lib = "font-awesome"),
                                       
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("Grup 1:")),
                                                    
                                                    numericInput("bune", "Veri 1", 100),
                                                    
                                                    numericInput("buneb", "Veri 2", 21),
                                                    
                                                    #numericInput("n1i", "Total", 121),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    p(strong("Grup 2:")),
                                                    
                                                    numericInput("bunec", "Veri 1", 120),
                                                    
                                                    numericInput("buned", "Veri 2", 67),
                                                    
                                                    #numericInput("n2i", "Total", 187),
                                                    
                                                    p(br())
                                             ),
                                             column(4,
                                                    radioButtons("twoxtwovalue", strong("Indeks seçimi"),
                                                                 c("log risk oranı" = "RR",
                                                                   "log olasılık oranı" = "OR",
                                                                   "Risk farkı" = "RD",
                                                                   "arcsine karekok transformlu risk farkı(Rücker et al., 2009)." = "AS",
                                                                   "Peto metodu ile tahminlenmis olasılık oranı (Yusuf et al., 1985)." = "PETO"
                                                                 ), selected = "OR"),
                                                    submitButton("Yenile")
                                             )
                                             
                                           )),
                                         
                                         
                                         h3("Etki Büyüklüğü and Örneklem Varyansı"),
                                         verbatimTextOutput("twobytwogroups.out"),
                                         
                                         br()
                                         
                                       )
                                       
                              ),
                              tabPanel("p değerinden etki büyüklüğü hesaplama", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("p değerinden etki büyüklüğü hesaplama")),
                                                    
                                                    numericInput("pvaluenum", " p-değeri", 0.01),
                                                    numericInput("pvaluen1", " Uygulama grubu örneklem sayısı", 50),                      
                                                    numericInput("pvaluen2", " Kontrol grubu uygulama sayısı", 50),
                                                    
                                                    radioButtons("pvaluetail", strong("Tek yada çift taraflı p değeri"),
                                                                 c("tek" = "one",
                                                                   "Çift" = "two"
                                                                 ), selected = "two"),
                                                    
                                                    p(br())
                                             ),
                                             column(4, offset = 1,
                                                    helpText("Yenilemek için tıklayınız"),
                                                    submitButton("Yenile"),
                                                    p(br())
                                             )
                                             
                                             
                                           )),
                                         
                                         
                                         h3("Etki Büyüklüğü"),
                                         verbatimTextOutput("pvaluees.out"),
                                         p(br()),
                                         
                                         br()
                                         
                                       )
                                       
                              ),
                              tabPanel("Tek denekli dizayn", icon = icon("chevron-right", lib = "font-awesome"),
                                       verticalLayout(
                                         
                                         wellPanel(
                                           fluidRow(
                                             column(3,
                                                    p(strong("Tek denekli dizayn")),
                                                    
                                                    radioButtons("SCDtype", strong("Tek denekli dizayn tipi"),
                                                                 c("AB" = "AB",
                                                                   "ABA" = "ABA",
                                                                   "ABAB" = "ABAB",
                                                                   "Tamamen rassal desen" = "CRD",
                                                                   "Rassal Bloklu desen" = "RBD"
                                                                   #"Alternating Treatments Design" = "ATD",
                                                                   #"Multiple-baseline AB design" = "MBD"
                                                                 ), selected = "AB"),
                                                    radioButtons("SCDes", strong("Effect Size"),
                                                                 c("Standardize ortalama farkı" = "SMD",
                                                                   "bileşkeli (pooled)  Standardize ortalama farkı" = "SMDpool"
                                                                   #"Percentage of Nonoverlapping Data (Positive)" = "PND+",
                                                                   #"Percentage of Nonoverlapping Data (Negative)" = "PND-",
                                                                   #"Percentage of Data Points Exceeding the Median (Positive)" = "PND+",
                                                                   #"Percentage of Data Points Exceeding the Median (Negative)" = "PND-"
                                                                 ), selected = "SMD"),
                                                    helpText("Yenilemek için tıklayınız"),
                                                    bsAlert("alert"),
                                                    submitButton("Yenile"),
                                                    
                                                    p(br())
                                             ),
                                             p(strong("Tek denekli desen veri girişi")),
                                             p("Sol sütun koşulu, sağ sütun puanları içermelidir. "),
                                             aceEditor("SCDdata", value="A, 9.523465\nA, 12.371462\nA, 13.265618\nA, 10.182837\nA, 10.987079\nA, 8.161392\nA, 10.655287\nA, 9.563863\nA, 9.381336\nA, 8.822936\nA, 10.227932\nA, 11.961484\nA, 9.425201\nA, 12.199128\nB, 16.212489\nB, 17.657583\nB, 18.45166\nB, 16.645105\nB, 14.618445\nB, 15.769643\nB, 16.017145\nB, 14.000921\nB, 17.081538\nB, 14.06722\nB, 20.423526\nB, 14.123096\nB, 16.728538", mode="r", theme="terminal"),
                                             p("Hesaplanmış etki büyüklüğünüz"),
                                             verbatimTextOutput('SCDES.out'),
                                             p(br()),
                                             h3("Referans"),
                                             p("Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467--478."),
                                             p("Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477--485.")
                                             
                                             
                                             
                                            )),
                                          
                                          p(br())
                                          
                                        )
                                     
                              )
                              ),
                   navbarMenu("Hakkında", icon = icon("tag", lib = "font-awesome"),
                              tabPanel("Hakkında", icon = icon("bookmark", lib = "font-awesome"),
                                       
                                       HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2015/04/mavis3.png" alt="" style="float: top; margin-right:5px" /></div>'),
                                       br(),
                                       strong('aRma Hakkında'),
                                       p("aRma, Meta Analyses via Shiny (MAVIS) yazılımının Türkçeleştirilmiş halidir. 
                                         MAVIS'in, dolayısıyla aRma'nın geliştirilme amacı, meta-analizlerin mümkün olduğunca kolay yapılabilmesidir. 
                                         Bu amaçla R programlama dilini Shiny arayüzü ile birleştirmeye çalışmıştır. aRma şu an test aşamasındadır. 
                                         Herhangi bir finansal destek almamaktadır ve almayacaktır, işe yaraması umuduyla küçük bir aRmağandır."),
                                     
                                       br(),
                                       strong("MAVIS Version 1.1"),
                                       p("MAVIS'in aylık indirilme sayısı "),
                                       img(src = "http://cranlogs.r-pkg.org/badges/MAVIS", seamless=NA),
                                       
                                       
                                       br()
                                       
                                       ),
                              tabPanel("Yazarlar ve Katkıda bulunanlar", icon = icon("users", lib = "font-awesome"),
                                       
                                       
                                       strong('Yazarlar'),
                                       
                                       HTML('<div style="clear: left;"><img src="http://oi59.tinypic.com/2mnrcci.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                                       p(a("Burak Aydın, PhD ", href="http://aydinburak.wix.com/test", target="_blank"),br(),
                                         p("İngilizce yazılımın geliştirilmesine yardım etmiş ve uygulamayı Türkçeye çevirmiştir.")
                                       ),
                                       
                                       br(),
                                                                           
                                       HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/kyle80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                                       p(a("William Kyle Hamilton - University of California, Merced", href="http://www.kylehamilton.com", target="_blank")),
                                       p("William Kyle Hamilton bu uygulamanın Ingilizce versiyonunun sahibi ve sorumlusudur. Kyle'ın websitesinden MAVIS'e ve kodlara ulaşabilirsiniz."),
                                       
                                       br(),
                                       HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2014/11/atsushi80.jpg" alt="" style="float: left; margin-right:5px" /></div>'),
                                       p(a("Atsushi Mizumoto, PhD - Kansai University", href="http://mizumot.com", target="_blank"),br(),
                                         p("Atsushi Mizumoto bu uygulamanın temelini atmıştır; uygulamanın ham hali için :", a("bkn", href="https://github.com/mizumot/meta", target="_blank"))
                                         
                                         
                                       ),
                                       
                                       br(),
                                       
                                       strong('Katkıda bulunanlar '),
                                       

                                       HTML('<div style="clear: left;"><img src="http://kylehamilton.com/wp-content/uploads/2015/04/katie80.png" alt="" style="float: left; margin-right:5px" /></div>'),
                                       p(a("Kathleen Coburn - University of California, Merced", href="http://psychology.ucmerced.edu/content/kathleen-coburn", target="_blank"),br(),
                                         p("Kathleen teknik destek vermiştir.")
                                       ),
                                       br()
                              ),
#                               tabPanel("Bug Reports", icon = icon("bug", lib = "font-awesome"),
#                                        
#                                        strong('Bug Reports'),
#                                        
#                                        p("If you discover a problem with MAVIS please submit it to the project GitHub page", 
#                                          a("https://github.com/kylehamilton/MAVIS/issues", href="https://github.com/kylehamilton/MAVIS/issues", target="_blank"),br()),
#                                        
#                                        p("MAVIS is an Open Source project, you are more than welcome to submit patches or features and help the project grow."),
#                                        
#                                        
#                                        br()
#                                        
#                               ),
                              tabPanel("Geri Dönüt", icon = icon("comments", lib = "font-awesome"),
                                       
                                       strong('aRma için geri dönüt'),
                                  
                                       p("Görüş, öneri ve geliştirme istekleri için burakaydin@ufl.edu"),
                                    
                                       
                                       br()
                                       
                                       ),
                              
                              tabPanel("License", icon = icon("legal", lib = "font-awesome"),
                                       
                                       strong('License'),
                                       
                                       p("MAVIS: Meta Analysis via Shiny"),
                                       p(" Copyright 2015  William Kyle Hamilton and Atsushi Mizumoto"),
                                       
                                       p(" This program is free software you can redistribute it and or modify
                                         it under the terms of the GNU General Public License as published by
                                         the Free Software Foundation either version 3 of the License or
                                         at your option any later version."),
                                       
                                       p("This program is distributed in the hope that it will be useful,
                                         but WITHOUT ANY WARRANTY; without even the implied warranty of
                                         MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
                                         GNU General Public License for more details."),
                                       
                                       p("You should have received a copy of the GNU General Public License
                                         along with this program.  If not, see", a("http://www.gnu.org/licenses/gpl.html", href="http://www.gnu.org/licenses/gpl.html", target="_blank"),br()),
                                       img(src = "http://www.gnu.org/graphics/gplv3-127x51.png", seamless=NA),
                                       
                                       
                                       br(),
                                       
                                       strong('Futher Infomation'),
                                       p("If you would like to learn more about the GNU General Public License and what it means tl'dr legal has a simple explaination which can be found here", a("https://www.tldrlegal.com/l/gpl-3.0", href="https://www.tldrlegal.com/l/gpl-3.0", target="_blank"),br()),
                                       
                                       
                                       
                                       br()
                                       
                                       ),
                              
                              tabPanel("Support", icon = icon("chevron-right", lib = "font-awesome"),
                                       
                                       
                                       strong('Support'),
                                       
                                       p("If you're having problems with MAVIS feel free to refer to our GitHub wiki or the documentation available on CRAN."),
                                       a("CRAN page for MAVIS", href="http://cran.r-project.org/web/packages/MAVIS/index.html", target="_blank"),
                                       br(),
                                       a("GitHub Wiki page for MAVIS", href="https://github.com/kylehamilton/MAVIS/wiki", target="_blank"),
                                       br(),
                                       p("As always you are more than welcome to contact the project maintainer at kyle.hamilton@gmail.com"),
                                       br()
                                       
                              )
),
                   
                   #This is just so I can get ui.R to run, I'll fix this later
                   tabPanel(" ",
                            h5(" ")
                   ),
                   
                   p(br())
                   
                              )
                              )