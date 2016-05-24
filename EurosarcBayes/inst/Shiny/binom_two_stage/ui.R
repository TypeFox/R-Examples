library(shiny)


numberInput<-function (inputId,label, value = "",...){
  div(class="form-group",style="display:block",tags$label(label, `for` = inputId),tags$input(id = inputId, type = "number", value = value,class="input-small",...))
}

shinyUI(
  fluidPage(
    theme="bootstrap_cerulean.css",

    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap_cerulean.css")
    ),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")
    ),

    # Application title
    div(style="text-align:center;",titlePanel("Single arm two stage phase II trial with a single binomial endpoint")),

    # Sidebar
    fluidRow(
      column(
        style="height:750px; text-align:right",
        wellPanel(
          tabsetPanel(
            tabPanel(
              "Main parameters",
              h4("Frequentist Constraints"),
              numberInput(
                "alpha",label=HTML("One sided type I error (&alpha;)"),value = 0.1,min=0,max=1,step=0.01),
              numberInput("power", label=HTML("Power (1-&beta;)"),value = 0.8,min=0,max=1,step=0.01),
              h4("Hypotheses"),
              numberInput("p0",label=HTML("Unacceptable probability (p<sub>0</sub>) of success"),value = 0.1,min=0,max=1,step=0.01),
              numberInput("p1",label=HTML("Acceptable probability (p<sub>1</sub>) of success"),value = 0.3,min=0,max=1,step=0.01)
            ),
            tabPanel(
              "Extra parameters",
              h4("Bayesian parameters"),
              numberInput("prior.a",label="Prior alpha",value = 0,min=0,step=1),
              numberInput("prior.b", label="Prior beta",value = 0,min=0,step=1),
              h4("Specific design"),
              numberInput("s1",label="Cut point futility at interim",value = 0,min=-1,step=1),
              numberInput("r1",label="Cut point efficacy at interim",value = 3,min=0,step=1),
              numberInput("n1",label="Number of patients at interim",value = 7,min=0,step=1),
              numberInput("r",label="Cut point efficacy at final analysis",value = 4,min=0,step=1),
              numberInput("n",label="Number of patients at final analysis",value = 18,min=0,step=1)
            ))),
        width = 4),

      # Main panels
      column(
        width=8,
        tabsetPanel(
          tabPanel(
            "Introduction",
            p(),
            p("This app calculates the required sample sizes cut points and properties for Simon's two stage design (Simon 1989). Other frequentist and Bayesian two stage designs for single endpoint trials are discussed and a calculator is provided for computing the properties of specific two stage designs in this setting. A Bayesian prior can be provided, and Bayesian parameters are also calculated. The default for this is a Beta(0,0) prior."),
            h4("Design space"),
            p(HTML("Designs discussed here are for single arm trials with a single binary endpoint, typically efficacy (response to treatment vs no response to treatment). There are advantages and disadvantages to running a multi-stage trial in a frequentist setting. The main advantage to a multi-stage trial is the ability to stop early if the treatment is not providing sufficient benefit to patients thus reducing the number of patients who received a treatment which does not work. The trial can also stop early if showing sufficient benefit to patients speeding up the process of movng to a confirmatory trial and thus getting drug approved and available quicker. The main drawback of frequentist approaches is that analysis must occur at a specific number of patients and the frequentist properties of these trials are often inflexible to changes to these numbers. As a consequence it is recommended that the trial stop recruiting once the required number of patients has been reached for interim analysis and only reopened after the DSMC has met and made their recomendations. Some designs do address this concern though.")),
            h4("Simon's two stage design"),
            p(HTML("Simon's design is the exact (analytic) version of Fleming's two stage design for a binomial endpoint. Fleming's design relies on Gaussian asymptotics which are poor for small sample sizes. The widespread use of computers means that the exact approach is easily computed and thus preferable. Both designs allow stopping at the single interim analysis for futility only.")),
            p(HTML("Simon proposed two designs optimising for two different properties of the design. He proposed an optimal design where the expected sample size under H<sub>0</sub> is minimsed. This design often required a higher maximum number of patients than a single stage design. To combat this he also proposed a min-max design (mini-max) This design aims to minimise the maximum sample size and then minimised the expected sample size under H<sub>0</sub>. This design can require a lower maximum number of patients than a single stage design but often it will require the same number of patients. This is due to the discrete properties of binomial distributions.")),
            h4("Design properties"),
            p(HTML("In the frequentist setting the type I error (&alpha;) is the probability that the data observed is at least as extreme as the critical value (r) if the true probability of a given event is p<sub>0</sub>. The power is the probability that the data observe is at least as extreme as the critical value (r) if the true probability of a given event is p<sub>1</sub>. These properties are often conisered as the long run proportion of trial that will satisfy the critical value given p<sub>0</sub> and p<sub>1</sub> respectively.")),
            p(HTML("In the Bayesian setting a single distribution for the unknown probability parameter p is computed using the prior information and data collected on trial. We require that with probability &eta; that the unknown p is greater than p<sub>0</sub> to stop the trial and conclude efficacy, and similarly with probability &zeta; that the unknown p is less than p<sub>1</sub>.")),
            p(HTML("The frequentist type I error, &alpha;, is approximately equivalent to the Bayesian parameter 1-&zeta; whilst the frequenitst power, &beta;, is approximately equivalent to the Bayesian parameter &eta;.")),
            h4("References"),
            p(HTML("Simon R. (1989). Optimal Two-Stage Designs for Phase II Clinical Trials. <em>Controlled Clinical Trials</em> ,10, 1-10."))
          ),
          tabPanel(
            "Simon's Two stage design",
            h3("Optimal design"),
            htmlOutput("optimal"),
            h3("Minimax design"),
            htmlOutput("minimax"),
            h4("Programs"),
            p(HTML("The frequentist design was computed using the function <strong>ph2simon</strong> which is part of the R package <strong>clinfun</strong> built and maintained by Venkatraman E. Seshan. These calculations can be reproduced using the function <strong>freq.binom.one.simons.twostage</strong> which is a wrapper for <strong>ph2simon</strong> and the design properties. See the help page for details."))
          ),
          tabPanel(
            "Other two stage designs",
            br(),
            h4("Stopping early for efficacy as well as futility"),
            p("A simple extention of Simon's two stage design is to allow stopping for efficacy only or both efficacy and futility at the interim analysis (Mander, 2010) This has the added benefit of speeding up very effective drugs research cycle to approval. Whilst this approach still requires stopping recruitment for interim anaylsis there are no more drawbacks than to Simon's two stage except for the lack of information about the drug"),
            p(HTML("Mander AP, Thompson SG. Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials. <em>Contemp Clin Trials</em> 2010; 31: 572-578.")),
            h4("Optimisation under another criteria"),
            p(HTML("Simon originally proposed two methods of optimisation, optimisation under H<sub>0</sub> and minimax under H<sub>0</sub>. There are many other proposals for the choice of optimisation strategy. These include: optimisation under H<sub>1</sub> and minimax under H<sub>1</sub> (Mander, 2010) and balanced designs (Ye, 2007)")),
            p(HTML("Ye F, Shyr Y. Balanced two-stage designs for phase II clinical trials. <em>Clin Trials</em> 2007; 4: 514-524.")),
            h4("Flexible designs"),
            p(HTML("Flexible designs were a response to the acceptance that practical considerations make it difficult to arrive at the planned sample size exactly (Chen 1998).A flexible design is defined as a collection of two-stage designs where the first stage size is in a set of consecutive values (n<sub>1</sub>,n<sub>2</sub>,... n<sub>k</sub>) and the second stage size is also in another set ofconsecutive values (N<sub>1</sub>,N<sub>2</sub>,...N<sub>k</sub>), and each of k<sup>2</sup> possible designs has the same probability of occurrence. Optimisation is dones over each collection of two-stage designs which will propose a design collection which is flexible to the sample size at interim and final analysis.")),
            p(HTML("Chen TT, Ng TH. Optimal flexible designs in phase II clinical trials. <em>Stat Med</em> 1998; 17: 2301-2312.")),
            h4("More than two stages"),
            p(HTML("The expected sample size saving for a three stage trial is around 10%. This is considerably less than the expected sample size saving for a 2 stage design. In may situations having to stop the trial an additional time will delay the trial more than having to recruit fewer patients on average. Computing optimal exact trial designs with more stages is exponentially computationally expensive.0")),
            p(HTML("Ensign LG, Gehan EA, Kamen DS, Thall PF. An optimal three-stage design for phase II clinical trials. <em>Stat Med</em> 1994; 13: 1727-1736.")),
            p(HTML("Chen TT. Optimal three-stage designs for phase II cancer clinical trials. <em>Stat Med</em> 1997; 16: 2701-2711.")),
            p(HTML("Chen K, Shan M. Optimal and minimax three-stage designs for phase II oncology clinical trials. <em>Contemp Clin Trials</em> 2008; 29: 32-41.")),
            h4("More endpoints"),
            p(HTML("It is possible to add an additional primary endpoint. The most common choice is toxicity. This endpoint can be made binary by defining  time period and a toxic event. There are a number of papers on this topic.")),
            p(HTML("Bryant J, Day R. Incorporating toxicity considerations into the design of two-stage phase II clinical trials. <em>Biometrics</em> 1995; 51: 1372-1383.")),
            p(HTML("Conaway MR, Petroni GR. Bivariate Sequential Designs for Phase II Trials. <em>Biometrics</em> 1995; 51: 656-664.")),
            p(HTML("Conaway MR, Petroni GR. Designs for Phase II Trials Allowing for a Trade-Off between Response and Toxicity. <em>Biometrics</em> 1996; 52: 1375-1386.")),
            h4("Bayesian designs"),
            p("All designs above are based on a frequentist framework. There are also a number of approaches which are based on a Bayesian frame work. In some cases these designs are identical, but typically Bayesian designs have smaller sample sizes. There is no free lunch however since designs with smaller sample sizes will have poorer frequentist properties."),
            p(HTML("Jung SH, Lee T, Kim K, George SL. Admissible two-stage designs for phase II cancer clinical trials. <em>Stat Med</em> 2004; 23: 561-569.")),
            p(HTML("Zhao L, Taylor JM, Schuetze SM. Bayesian decision theoretic two-stage design in phase II clinical trials with survival endpoint. <em>Stat Med</em> 2012; 31: 1804-1820."))

          ),
          tabPanel(
            "A note on analysis",
            h3("Analysis of a multistage trial is not the same as analysis of a single stage trial"),
            p(HTML("A common mistake in multi stage trials is to report the p-values the way you would report them for a single stage trial. The porblem with doing this is that if the trial goes to final analysis the standard p-value will not account for the fact that you could have stopped at interim analysis. An analytic adjustment must be made for this if p-values are to be calculated")),
            p(HTML("Koyama T, Chen H. Proper inference from Simon's two-stage designs. <em>Stat Med</em> 2008; 27: 3145-3154."))
          ),
          tabPanel(
            "Specific design",
            br(),
            p("The design here is calculated from the specific trial design. Stopping early for efficacy is also allowed. If you do not wish to stop for futility enter -1 for the interim futility cut point. If you do not wish to stop for efficacy enter a value greater than the number of patients at interim into the interim efficacy cut point."),
            h3("Specified design"),
            htmlOutput("custom"),
            h4("Programs"),
            p(HTML("The properties here were computed using the function <strong>binom.single.trial.properties</strong>."))
          )
        )

      )
    )
  )
  )
