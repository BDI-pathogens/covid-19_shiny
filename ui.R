ui <- navbarPage("Digital contact tracing for SARS-COV-2", 
                 
                 tags$head(includeHTML(("google-analytics.html"))), # google analytics token
                 
  theme = shinytheme("spacelab"), # change the theme here; use "themeSelector()" below for an easy way to pick
  #shinythemes::themeSelector(), # use this to dynamically select a theme
  
  tabsetPanel(
    tabPanel(
      "Infectiousness",
      HTML('<meta name="viewport" content="width=1024">'), # forces "desktop view" on mobile
      
      # Sidebar with a slider inputs etc
      sidebarLayout(
        sidebarPanel(
          id = "sidePanel",
          style = "overflow-y: auto; max-height: 90vh; position:relative;", # make it scrollable
          useShinyjs(), # for the reset button
          withMathJax(), # for LaTeXing
          
          h5("Incubation period (days):"),
          sliderInput("incperMedian",
                      h6("Median:"),
                      min = 1,
                      max = 12, # https://github.com/HopkinsIDD/ncov_incubation#parameter-estimates
                      step = 0.01,
                      value = exp(1.64)),
          
          sliderInput("incperSdlog",
                      h6("sdlog parameter:"),
                      min = 0.01,
                      max = 0.99,
                      step = 0.01,
                      value = 0.36),
          
          # mini plot to show the chosen incubation period distribution
          h6("Chosen incubation period lognormal distribution"),
          withSpinner(plotOutput("IncperDistribPlot", height="130px"), type=7),
          
          hr(),
          
          h5("Generation time (days):"),
          sliderInput("serIntShape",
                      h6("Shape:"),
                      min = 1,
                      max = 7,
                      step = 0.01,
                      value = 2.83),
          
          sliderInput("serIntScale",
                      h6("Scale:"),
                      min = 1,
                      max = 7,
                      step = 0.01,
                      value = 5.67),
         
          # mini plot to show the chosen generation time distribution
          h6("Chosen generation time Weibull distribution"),
          withSpinner(plotOutput("SerintDistribPlot", height="130px"), type=7),
          
          hr(),
          
          sliderInput("doublingTime",
                      h5("Epidemic doubling time (days):"),
                      min = 1,
                      max = 10,
                      step = 0.01,
                      value = 5),
          
          hr(),
          
          sliderInput("xa",
                      h5("Relative infectiousness of asymptomatic compared with symptomatic individuals:"),
                      min = 0,
                      max = 2,
                      step = 0.01,
                      value = 0.1),
          
          hr(),
          
          sliderInput("P.a",
                      h5("Fraction of infected individuals who are asymptomatic:"),
                      min = 0,
                      max = 0.99, # 1 gives an error
                      step = 0.01,
                      value = 0.4),
          
          sliderInput("frac.Re",
                      h5("Fraction of transmissions that are environmentally mediated:"),
                      min = 0,
                      max = 0.99, # 1 gives an error
                      step = 0.01,
                      value = 0.1),
          
          hr(),
          
          h5("Environmental contribution:"),
          radioButtons("env.type",
                       label = "",
                       choices = list("Environmental infectiousness is constant" = "constant",
                                      "Environmental infectiousness decays exponentially" = "exp.decay"),
                       selected = "constant"),
          
          # dependent on the selection of env inf type above, give an appropriate slider to determine the rate
          uiOutput("environment_sliders"),
          
          hr()
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          id = "mainPanel",
          style = "overflow-y: auto; max-height: 90vh; position:relative;", # make it scrollable
          
          # make a "reset" button
          fluidRow(column=3, 
                   align="right", 
                   actionButton("reset_input", "Reset values")
          ),
          
          h3("Average rate of infecting others as a function of the time since infection"),
          
          # I would love to control the height of the plot dynamically so that it's, say, 50% of the window height,
          # but the documentation for "plotOutput" says it might behave weirdly, and it does!
          # Tried using shinyjqui resizeable but it's a bit buggy. Other suggestions welcome!
          withSpinner(plotOutput("mainPlot", height="400px"), type=7),
          
          withSpinner(plotOutput("decompositionPlot", height="200px"), type=7),
          
          withMathJax(uiOutput("parameterSummary")), # withMathJax enables LaTeX to be rendered
          
          hr()
        ) # end "main" panel (the bit with the graphs and the explanations)
      ) # end "sidebarLayout"
   ),# end "Infectiousness" tab
   
   
   tabPanel(
     "Interventions",
     HTML('<meta name="viewport" content="width=1024">'), # forces "desktop view" on mobile
     
     # Sidebar with a slider inputs etc
     sidebarLayout(
       sidebarPanel(
         id = "sidePanel2",
         style = "overflow-y: auto; max-height: 90vh; position:relative;", # make it scrollable
         useShinyjs(), # for the reset button
         withMathJax(), # for LaTeXing
         
         h5("Delay:"),
         sliderInput("delay",
                     h6("Hours from symptoms to isolation and quarantining of contacts:"),
                     min = 0,
                     max = 72,
                     step = 1,
                     value = 0),
         
         hr(),
         
         h5("Incubation period (days):"),
         sliderInput("incperMedianContour",
                     h6("Median:"),
                     min = 1,
                     max = 12, # https://github.com/HopkinsIDD/ncov_incubation#parameter-estimates
                     step = 0.01,
                     value = exp(1.64)),
         
         sliderInput("incperSdlogContour",
                     h6("sdlog parameter:"),
                     min = 0.01,
                     max = 0.99,
                     step = 0.01,
                     value = 0.36),
         
         # mini plot to show the chosen incubation period distribution
         h6("Chosen incubation period lognormal distribution"),
         withSpinner(plotOutput("IncperDistribPlotContour", height="130px"), type=7),
         
         hr(),
         
         h5("Generation time (days):"),
         sliderInput("serIntShapeContour",
                     h6("Shape:"),
                     min = 1,
                     max = 7,
                     step = 0.01,
                     value = 2.83),
         
         sliderInput("serIntScaleContour",
                     h6("Scale:"),
                     min = 1,
                     max = 7,
                     step = 0.01,
                     value = 5.67),
         
         # mini plot to show the chosen generation time distribution
         h6("Chosen generation time Weibull distribution"),
         withSpinner(plotOutput("SerintDistribPlotContour", height="130px"), type=7),
         
         hr(),
         
         sliderInput("doublingTimeContour",
                     h5("Epidemic doubling time before interventions (days):"),
                     min = 1,
                     max = 10,
                     step = 0.01,
                     value = 5),
         
         hr(),
         
         sliderInput("xaContour",
                     h5("Relative infectiousness of asymptomatic compared with symptomatic individuals:"),
                     min = 0,
                     max = 2,
                     step = 0.01,
                     value = 0.1),
         
         hr(),
         
         sliderInput("P.aContour",
                     h5("Fraction of infected individuals who are asymptomatic:"),
                     min = 0,
                     max = 0.99, # 1 gives an error
                     step = 0.01,
                     value = 0.4),
         
         sliderInput("frac.ReContour",
                     h5("Fraction of transmissions that are environmentally mediated:"),
                     min = 0,
                     max = 0.99, # 1 gives an error
                     step = 0.01,
                     value = 0.1),
         
         hr(),
         
         h5("Environmental contribution:"),
         radioButtons("env.typeContour",
                      label = "",
                      choices = list("Environmental infectiousness is constant" = "constant",
                                     "Environmental infectiousness decays exponentially" = "exp.decay"),
                      selected = "constant"),
         
         # dependent on the selection of env inf type above, give an appropriate slider to determine the rate
         uiOutput("environment_slidersContour"),
         
         hr()
       ),
       
       # Show a plot of the generated distribution
       mainPanel(
         id = "mainPanel2",
         style = "overflow-y: auto; max-height: 90vh; position:relative;", # make it scrollable
         
         # make a "reset" button
         fluidRow(column=3, 
                  align="right", 
                  actionButton("reset_inputContour", "Reset values")
         ),
         
         h3("Daily epidemic growth rate r, with interventions of varying efficiencies"),
         
         withMathJax(uiOutput("contourDescribeDelay")),
         
         withSpinner(plotOutput("contourPlot", height="600px", width="775px"), type=7),
         
         withMathJax(uiOutput("contourDescribeR0")),
         
         h5("Please note that our published results use a higher resolution and a tighter tolerance for convergence;
                  for ease of use we perform a faster approximation here."),
         
         hr()
       ) # end "main" panel (the bit with the graphs and the explanations)
     ) # end "sidebarLayout"
   ), # end "Interventions" tab
  
  tabPanel("Details", 
           withMathJax(includeMarkdown("markdown/details.md"))
  ), # end "Details" tab
  
  tabPanel("About", 
           withMathJax(includeMarkdown("markdown/about.md")),
           verbatimTextOutput("systeminfo") # server info
  ) # end "Details" tab
  
  ) # end tabsetPanel
  
                 
                 
) # end ui

