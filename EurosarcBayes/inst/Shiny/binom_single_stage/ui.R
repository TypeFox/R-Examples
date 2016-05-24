library(shiny)


numberInput<-function (inputId,label, value = "",...){
  div(class="form-group",style="display:block",tags$label(label, `for` = inputId),tags$input(id = inputId, type = "number", value = value,class="input-small",...))
}

shinyUI(
  fluidPage(
    theme="bootstrap_cerulean.css",

    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap_cerulean.css")),
    tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css")),


    # Application title
    div(style="text-align:center;",titlePanel("Single arm phase II trial for a single binomial endpoint")),

    # Sidebar
    sidebarLayout(
      sidebarPanel(
        style="text-align:right",
        numberInput("prior.a",label="Prior alpha",value = 1,min=0,step=1),
        numberInput("prior.b", label="Prior beta",value = 1,min=0,step=1),
        numberInput("p0",label=HTML("Unacceptable probability (p<sub>0</sub>) of success"),value = 0.2,min=0,max=1,step=0.01),
        numberInput("p1",label=HTML("Acceptable probability (p<sub>1</sub>) of success"),value = 0.35,min=0,max=1,step=0.01),
        numberInput("eta",label=HTML("Probability of being greater than p<sub>0</sub> (&eta;)"),value = 0.9,min=0,max=1,step=0.01),
        numberInput("zeta",label=HTML("Probability of being less than p<sub>1</sub> (&zeta;)"),value = 0.9,min=0,max=1,step=0.01),
        p(HTML("1-&zeta; is the Bayesian equivalent of frequentist type I error (&alpha;).")),
        p(HTML("&eta; is the Bayesian equivalent of frequentist power (&beta;).")),
        width = 4),


      mainPanel(
        width=8,
        tabsetPanel(
          tabPanel(
            "Introduction",
            p(),
            p("This app calculates the required sample size for a single arm binomial trial in both the Bayesian and frequentist methodologies."),
            p(HTML("For both methodologies the data is binary and so modelled from a binomial distribution with some probability parameter p. The hypotheses are: H<sub>0</sub>:p=p<sub>0</sub> and H<sub>0</sub>:p=p<sub>1</sub>. In the bayesian setting a beta prior is used to include prior information about p.")),
            h4("Frequentist"),
            p(HTML("In the frequentist setting the type I error (&alpha;) is the probability that the data observed is at least as extreme as the critical value (r) if the true probability of a given event is p<sub>0</sub>. The power is the probability that the data observe is at least as extreme as the critical value (r) if the true probability of a given event is p<sub>1</sub>. These properties are often conisered as the long run proportion of trial that will satisfy the critical value given p<sub>0</sub> and p<sub>1</sub> respectively.")),
            h4("Bayesian"),
            p(HTML("In the Bayesian setting a single distribution for the unknown probability parameter p is computed using the prior information and data collected on trial. We require that with probability &eta; that the unknown p is greater than p<sub>0</sub> to stop the trial and conclude efficacy, and similarly with probability &zeta; that the unknown p is less than p<sub>1</sub>.")),
            h4("Minimising sample size"),
            p("In both setting the smallest sample size n and a corresponding cut-off can be computed to satisfy both conditions."),
            p(HTML("The frequentist type I error, &alpha;, is approximately equivalent to the Bayesian parameter 1-&zeta; whilst the frequenitst power, &beta;, is approximately equivalent to the Bayesian parameter &eta;."))


          ),
          tabPanel(
            "Results",
            h3("Bayesian single arm design"),
            htmlOutput("eff.sample.size"),
            htmlOutput("sample.size"),
            htmlOutput("sample.success"),
            p(strong("Bayesian Properties")),
            htmlOutput("sample.eta"),
            htmlOutput("sample.zeta"),
            p(),
            p(strong("Frequentist Properties")),
            htmlOutput("sample.alpha"),
            htmlOutput("sample.power"),
            p(),
            h3("Frequentist single arm design"),
            htmlOutput("f.sample.size"),
            htmlOutput("f.sample.success"),
            p(),
            p(strong("Frequentist Properties")),
            htmlOutput("f.sample.alpha"),
            htmlOutput("f.sample.power"),
            p(),
            p(strong("Bayesian Properties")),
            htmlOutput("f.sample.eta"),
            htmlOutput("f.sample.zeta")
          ),
          tabPanel(
            "Visualise",
            div(style="text-align:center;",h3("Prior distribution with required probability for a successful trial")),
            p("This plot allows you to compare the chosen prior distribution (blue) from the smallest positive result (vertical lines). Changing the prior distribution will affect the Bayesian probability threshold but not the frequntist one. It indicates that a conservative prior distribution negatively affect the required observed probabily to report a successful trial (increases the required probability above that of the frequentist probability), whilst an optimistic does the converse."),
            plotOutput("priorplot")
          )
        )
      )
    )
  ))
