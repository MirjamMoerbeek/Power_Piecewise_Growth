library(shiny)
library(shinydashboard)
library(shinythemes)
library(plotly)
library(DT)

ui <- dashboardPage(skin="yellow",   
  dashboardHeader(title = "Power analysis of longitudinal studies with piecewise linear growth and attrition.",disable = TRUE),
  dashboardSidebar(disable = TRUE),
  dashboardBody(
    titlePanel("Power analysis of longitudinal studies with piecewise linear growth and attrition"),
       fluidRow(
         
      column(width=3,  height = 450,
        box(title = "Design", width = NULL, 
          numericInput("D", "Duration of the study", 12, min = 5, max = 12),
          textInput('distr', 'Enter distribution of turning points: for each turning point T=0,..., D enter proportion (comma delimited).', "0.0,0.0,0.0,0.0,0.1,0.2,0.4,0.2,0.1,0.0,0.0,0.0,0.0"),
          plotOutput("distributionplot")
        )
      ),
      
      column(width=3,height = 450,
             box(title = "Correlation parameters", width = NULL, 
               box(width=6,title="Variances",
                   numericInput("var.u0", "Variance of random intercept", 0.2),
                   numericInput("var.u1", "Variance of random slope 1", 0.1),
                   numericInput("var.u2", "Variance of random slope 2", 0.16),
                   numericInput("var.e", "Residual variance", 0.2)
               ),
             
               box(width=6,title="Covariances",
                   numericInput("covar.u01", "Covariance of random intercept and random slope 1", 0.1),
                   numericInput("covar.u02", "Covariance of random intercept and random slope 2", 0.1),
                   numericInput("covar.u12", "Covariance of random slope 1 and random slope 2", 0.1))
               ),
             box(title = "Create graphs",  width=NULL, background = "yellow",
                 submitButton("Submit")
             )
      ),
      
     column(width=3,height = 450,
            box(title = "Growth rates and test choices", width = NULL, 
              numericInput("beta1", "Growth rate in phase 1", 0.16, min = 0, max = 0.8),
              numericInput("beta2", "Growth rate in phase 2", 0.11, min = 0, max = 0.8),
              numericInput("alpha", "Type I error rate (alpha)", 0.05, min = 0, max = 0.8),
              selectInput("test", label = "Type of test", 
                          choices = list("One-sided" = 1, "Two-sided" = 2), 
                          selected = 2))
     ),
     
     column(width=3,height = 450,
            box(title = "Attrition function", width = NULL, 
                 numericInput("omega", "Omega", 0.5, min = 0, max = 1),
                 numericInput("gamma", "Gamma", 1, min = 0, max = 10),
                 plotOutput("survivalplots")
        )
      )
    ),

    fluidRow(  

   column(width=12,height = 450,
         box( title = "Power versus sample size for mean growth rate in phase 1",width=4,height=540,
              plotlyOutput("Resultsplot1", height=480)),
         box( title = "Power versus sample size for mean growth rate in phase 2",width=4,height=540,
              plotlyOutput("Resultsplot2",height=480)),
         box( title = "Power versus sample size for difference in mean growth rates",width=4,height=540,
              plotlyOutput("Resultsplot3",height=480))
    ))
  )
)
