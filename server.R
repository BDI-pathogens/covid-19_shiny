# loading libraries, global values and initial starting parameters is handled in "global.R"

# Define server logic
server <- function(input, output, session) {

  # define the reset button
  observeEvent(input$reset_input, { 
    shinyjs::reset("sidePanel")
  })
  
  ################### 
  # Get parameters 
  ###################
  
  # Parameters start using the values from "initial_params.RData" and then can be altered by the user
  getDoublingTime <- reactive({
    input$doublingTime
  })
  
  getEnvInfType <- reactive({
    input$env.type
  })
  
  getEnvDecay <- reactive({
    input$envDecayRate
  })
  
  getEnvConst <- reactive({
    input$envConstant
  })
  
  getIncperMeanlog <- reactive({
    log(input$incperMedian)
  })
  
  getIncperSdlog <- reactive({
    input$incperSdlog
  })

  getSerIntShape <- reactive({
    input$serIntShape
  })  
  
  getSerIntScale <- reactive({
    input$serIntScale
  })
  
  # xp used to be reactive, now we fix it to equal 1
  # getXp <- reactive({
  #   input$xp # the relative infectiousness of pre-symptomatic c.f. symptomatic individuals
  # })
  getXp <- reactive({
    1
  })
  
  getXa <- reactive({
    input$xa # the relative infectiousness of asymptomatic c.f. symptomatic individuals
  })
  
  getP.a <- reactive({
    input$P.a # the fraction of all infections that are asymptomatic
  })
  
  getFrac.Re <- reactive({
    input$frac.Re # the fraction of all transmissions that are environmentally mediated
  })
  
  ################
  # Calculations
  ###############
  
  # following Chris's method, create "dummy" which uses model.gen.solve and other functions to create a list of the values for
  # env.scale.constant
  # R0 
  # RSorP
  # RA
  # RS 
  # RP
  # RE
  # theta.obs.predicted 
  getDummy <- reactive({
    # get parameters
    frac.Re <- getFrac.Re()
    P.a <- getP.a()
    doubling.time <- getDoublingTime()
    xp <- getXp()
    xa <- getXa()
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    serint.shape <- getSerIntShape()
    serint.scale <- getSerIntScale()
    theta.obs <- 0.83
    env.decay.rate <- getEnvDecay()
    env.constant.duration <- getEnvConst()
    env.infectiousness.type <- getEnvInfType() 
    
    # get reactive functions 
    model.gen.beta.env.div.by.E.RSorP <- getModelGenBetaEnvDivByERSorP()
    model.gen.full.beta.div.by.RSorP <- getModelGenFullBetaDivByRSorP()

    # use "model.gen.solve" in a reactive way:
    r <- log(2) / doubling.time # units of per day
    
    integral.of.model.gen.beta.s.div.by.RSorP <-
      integrate(model.gen.beta.s.div.by.RSorP, 
                lower = 0, 
                upper = Inf,
                incper.meanlog = incper.meanlog, 
                incper.sdlog = incper.sdlog,
                serint.scale = serint.scale, 
                serint.shape = serint.shape,
                P.a = P.a, 
                xp = xp)$value
    
    env.scale.constant <- (1 + P.a * xa * integral.of.model.gen.beta.s.div.by.RSorP) / 
      (((1 / frac.Re) - 1) * integrate(model.gen.beta.env.div.by.E.RSorP,
                                       lower = 0, 
                                       upper = Inf,
                                       serint.scale = serint.scale,
                                       serint.shape = serint.shape,
                                       incper.meanlog = incper.meanlog,
                                       incper.sdlog = incper.sdlog,
                                       P.a = P.a,
                                       xp = xp,
                                       env.decay.rate = env.decay.rate,
                                       env.constant.duration = env.constant.duration,
                                       env.infectiousness.type = env.infectiousness.type)$value)
    
    
    RSorP <- 1 / integrate(function(tau) {
      model.gen.full.beta.div.by.RSorP(tau = tau, 
                                       serint.scale = serint.scale,
                                       serint.shape = serint.shape,
                                       P.a = P.a, 
                                       xa = xa, 
                                       incper.meanlog = incper.meanlog,
                                       incper.sdlog = incper.sdlog, 
                                       xp = xp,
                                       env.scale.constant = env.scale.constant,
                                       env.decay.rate = env.decay.rate,
                                       env.constant.duration = env.constant.duration,
                                       env.infectiousness.type = env.infectiousness.type) *
        exp(-r * tau)}, lower = 0, upper = Inf)$value
    
    R0 <- RSorP * integrate(model.gen.full.beta.div.by.RSorP, 
                            lower = 0, 
                            upper = Inf,
                            serint.scale = serint.scale,
                            serint.shape = serint.shape,
                            P.a = P.a, 
                            xa = xa, 
                            incper.meanlog = incper.meanlog,
                            incper.sdlog = incper.sdlog, 
                            xp = xp,
                            env.scale.constant = env.scale.constant,
                            env.decay.rate = env.decay.rate,
                            env.constant.duration = env.constant.duration,
                            env.infectiousness.type = env.infectiousness.type)$value
    
    should.equal.one <- integrate(function(tau) {exp(-r * tau) * RSorP *
        model.gen.full.beta.div.by.RSorP(tau = tau, 
                                         serint.scale = serint.scale,
                                         serint.shape = serint.shape,
                                         P.a = P.a, 
                                         xa = xa, 
                                         incper.meanlog = incper.meanlog,
                                         incper.sdlog = incper.sdlog, 
                                         xp = xp,
                                         env.scale.constant = env.scale.constant,
                                         env.decay.rate = env.decay.rate,
                                         env.constant.duration = env.constant.duration,
                                         env.infectiousness.type = env.infectiousness.type)},
        lower = 0, upper = Inf)$value
    
    RS <- integrate(model.gen.beta.sym.tot, 
                    lower = 0, 
                    upper = Inf,
                    incper.meanlog = incper.meanlog,
                    incper.sdlog = incper.sdlog,
                    serint.scale = serint.scale,
                    serint.shape = serint.shape,
                    P.a = P.a, 
                    xp = xp, 
                    RSorP = RSorP)$value
    
    RP <- integrate(model.gen.beta.presym.tot, 
                    lower = 0, 
                    upper = Inf,
                    incper.meanlog = incper.meanlog,
                    incper.sdlog = incper.sdlog,
                    serint.scale = serint.scale,
                    serint.shape = serint.shape,
                    P.a = P.a, 
                    xp = xp, 
                    RSorP = RSorP)$value
    
    RA <- RSorP * P.a * xa * integral.of.model.gen.beta.s.div.by.RSorP
    
    RE <- env.scale.constant * 
      RSorP * 
      integrate(function(tau) {
        model.gen.beta.env.div.by.E.RSorP(tau = tau, 
                                          serint.scale = serint.scale,
                                          serint.shape = serint.shape, 
                                          incper.meanlog = incper.meanlog,
                                          incper.sdlog = incper.sdlog, 
                                          P.a = P.a, 
                                          xp = xp, 
                                          env.decay.rate = env.decay.rate,
                                          env.constant.duration = env.constant.duration,
                                          env.infectiousness.type = env.infectiousness.type)},
        lower = 0, upper = Inf)$value

      theta.obs.predicted <- 1 - integrate(function(tau) {
        model.gen.beta.sym.tot(tau = tau,
                               incper.meanlog = incper.meanlog,
                               incper.sdlog = incper.sdlog,
                               serint.scale = serint.scale,
                               serint.shape = serint.shape,
                               P.a = P.a, 
                               xp = xp, 
                               RSorP = RSorP) * 
          exp(-r * tau)},
        lower = 0, upper = Inf)$value
                             
    # return
    list(env.scale.constant = env.scale.constant, 
         R0 = R0, 
         RSorP = RSorP, 
         RA = RA,
         RS = RS, 
         RP = RP, 
         RE = RE, 
         theta.obs.predicted = theta.obs.predicted)
  })
  
  # extract each component from the list "dummy"
  getR0 <- reactive({
    dummy <- getDummy()
    dummy$R0
  })
  
  getRSorP <- reactive({
    dummy <- getDummy()
    dummy$RSorP
  })
  
  getRe <- reactive({
    dummy <- getDummy()
    dummy$RE
  })
  
  getRa <- reactive({
    dummy <- getDummy()
    dummy$RA
  })
  
  getRs <- reactive({
    dummy <- getDummy()
    dummy$RS
  })
  
  getRp <- reactive({
    dummy <- getDummy()
    dummy$RP
  })
  
  getTheta <- reactive({
    rp <- getRp()
    R0 <- getR0()
    
    1 - rp / R0
  })
  
  getEnvScaleConst <- reactive({
    dummy <- getDummy()
    dummy$env.scale.constant
  })


  # following "model.gen.beta.env.div.by.E.RSorP"
  getModelGenBetaEnvDivByERSorP <- reactive({
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    P.a <- getP.a()
    xp <- getXp()
    env.decay.rate <- getEnvDecay()
    env.constant.duration <- getEnvConst()
    env.infectiousness.type <- getEnvInfType()
    
    # because of the way that "env.decay.rate" and "env.constant.duration" depend on
    # "env.infectiousness.type", there is a delay when shiny is launched, when these are undefined
    # so we use "req" to wait until they are defined before continuing the calculation
    if (env.infectiousness.type == "constant") {
      req(env.constant.duration)
    } else { 
      req(env.decay.rate)
    }
    
    Vectorize(function(tau, 
                       incper.meanlog, 
                       incper.sdlog, 
                       serint.scale,
                       serint.shape, 
                       P.a, 
                       xp, 
                       env.decay.rate,
                       env.constant.duration, 
                       env.infectiousness.type) {
      integrate(function(l) {
        model.gen.beta.s.div.by.RSorP(tau = tau - l,
                                      incper.meanlog = incper.meanlog,
                                      incper.sdlog = incper.sdlog,
                                      serint.scale = serint.scale,
                                      serint.shape = serint.shape,
                                      P.a = P.a,
                                      xp = xp) *
          p.el(l = l, 
               env.decay.rate = env.decay.rate,
               env.constant.duration = env.constant.duration,
               env.infectiousness.type = env.infectiousness.type)
        },
        lower = 0, upper = tau)$value
    }, vectorize.args = "tau")
      
    
  })
  
  # following model.gen.full.beta.div.by.RSorP
  getModelGenFullBetaDivByRSorP <- reactive({
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    P.a <- getP.a()
    xa <- getXa()
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    xp <- getXp()
    env.decay.rate <- getEnvDecay()
    env.constant.duration <- getEnvConst()
    env.infectiousness.type <- getEnvInfType()
    
    model.gen.beta.env.div.by.E.RSorP <- getModelGenBetaEnvDivByERSorP()
    
    function(tau, serint.scale, serint.shape, P.a, xa, incper.meanlog,
             incper.sdlog, xp, env.scale.constant, env.decay.rate,
             env.constant.duration, env.infectiousness.type) {
      serint(x = tau, 
             serint.scale = serint.scale,
             serint.shape = serint.shape) *
        (1 + P.a * xa / model.gen.f(tau = tau, 
                                    incper.meanlog = incper.meanlog,
                                    incper.sdlog = incper.sdlog, 
                                    P.a = P.a, 
                                    xp = xp)) +
        env.scale.constant * # this is a free variable of this function
        model.gen.beta.env.div.by.E.RSorP(tau = tau, 
                                          serint.scale = serint.scale,
                                          serint.shape = serint.shape, 
                                          incper.meanlog = incper.meanlog,
                                          incper.sdlog = incper.sdlog, 
                                          P.a = P.a, 
                                          xp = xp, 
                                          env.decay.rate = env.decay.rate,
                                          env.constant.duration = env.constant.duration,
                                          env.infectiousness.type = env.infectiousness.type)
    }
  })
  

  ################################
  # Build data frames of results
  ################################

  # a reactive version of df.beta.p
  getDFBetaP <- reactive({
    # get parameters
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    P.a <- getP.a()
    xp <- getXp()
    RSorP <- getRSorP()
    
    data.frame(tau = tau.test, 
               label = "pre-symptomatic", 
               beta = vapply(tau.test, 
                             model.gen.beta.presym.tot, 
                             numeric(1), 
                             incper.meanlog = incper.meanlog,
                             incper.sdlog = incper.sdlog, 
                             serint.scale = serint.scale,
                             serint.shape = serint.shape, 
                             P.a = P.a, 
                             xp = xp, 
                             RSorP = RSorP))
  })


  # a reactive version of df.beta.s
  getDFBetaS <- reactive({
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    P.a <- getP.a()
    xp <- getXp()
    RSorP <- getRSorP()
    
    data.frame(tau = tau.test, 
               label = "symptomatic", 
               beta = vapply(tau.test, 
                             model.gen.beta.sym.tot, 
                             numeric(1), 
                             incper.meanlog = incper.meanlog,
                             incper.sdlog = incper.sdlog, 
                             serint.scale = serint.scale,
                             serint.shape = serint.shape, 
                             P.a = P.a, 
                             xp = xp, 
                             RSorP = RSorP))
  })
  
  # a reactive version of df.beta.a
  getDFBetaA <- reactive({
    xp <- getXp()
    xa <- getXa()
    P.a <- getP.a()
    RSorP <- getRSorP()
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    
    data.frame(tau = tau.test, 
               label = "asymptomatic", 
               beta = vapply(tau.test, function(tau) {
                 RSorP * P.a * xa * model.gen.beta.s.div.by.RSorP(tau = tau,
                                                                  incper.meanlog = incper.meanlog, 
                                                                  incper.sdlog = incper.sdlog,
                                                                  serint.scale = serint.scale,
                                                                  serint.shape = serint.shape,
                                                                  P.a = P.a, xp = xp)} , numeric(1)))
  })
  
  # a reactive version of df.beta.e
  getDFBetaE <- reactive({
    model.gen.beta.env.div.by.E.RSorP <- getModelGenBetaEnvDivByERSorP()
    xp <- getXp()
    P.a <- getP.a()
    RSorP <- getRSorP()
    serint.scale <- getSerIntScale()
    serint.shape <- getSerIntShape()
    incper.meanlog <- getIncperMeanlog()
    incper.sdlog <- getIncperSdlog()
    env.decay.rate <- getEnvDecay()
    env.constant.duration <- getEnvConst()
    env.infectiousness.type <- getEnvInfType()
    env.scale.constant <- getEnvScaleConst()
    
    df.beta.e <- data.frame(tau = tau.test, 
                            label = "environmental", 
                            beta = vapply(tau.test, 
                                          model.gen.beta.env.div.by.E.RSorP, 
                                          numeric(1),
                                          serint.scale = serint.scale,
                                          serint.shape = serint.shape, 
                                          incper.meanlog = incper.meanlog,
                                          incper.sdlog = incper.sdlog, 
                                          P.a = P.a,
                                          xp = xp, 
                                          env.decay.rate = env.decay.rate,
                                          env.constant.duration = env.constant.duration,
                                          env.infectiousness.type = env.infectiousness.type)
                            )
    
    df.beta.e$beta <- df.beta.e$beta * env.scale.constant * RSorP
    df.beta.e
  })

  # collect each df.beta into a data frame, for plotting
  getDF <- reactive({
    validate(
      need(try(df.beta.p <- getDFBetaP()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.s <- getDFBetaS()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.e <- getDFBetaE()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.a <- getDFBetaA()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    
    rbind(df.beta.p, df.beta.s, df.beta.e, df.beta.a) # if you change the order of these you also need to change the order of the legend labels!
  })


  ################
  # For plotting
  ################
  
  # get ymax for plotting
  getYmax <- reactive({
    df.plot <- getDF()
    df.plot.wide <- spread(df.plot, label, beta)
    df.plot.wide$max <- df.plot.wide$`pre-symptomatic` + df.plot.wide$symptomatic +
      df.plot.wide$environmental + df.plot.wide$asymptomatic
    ymax <- max(df.plot.wide$max * 1.05)
  })
  
  # an example of how to set the colours; this can be a user input if desired
  getCols <- reactive({
    brewer.pal(8, "Paired")[c(1,3,5,7,2)] # palette which is colour-blind friendly but doesn't risk implying any pairing. Four colours for main plot, final colour for mini plots
  })
  
  # get the main plot using the latest data frame
  getPlot <- reactive({
    df.plot <- getDF()
    ymax <- getYmax()
    cols <- getCols()
    RS <- getRs()
    RP <- getRp()
    RE <- getRe()
    RA <- getRa()
    R0 <- getR0()
    
    mainPlot <- ggplot(df.plot, aes(x=tau, y=beta)) + #, color=label)) +
      theme_bw(base_size = 18) +
      geom_area(aes(fill=label)) +
      labs(x = expression(paste(tau, " (days)")),
           y = expression(paste(beta, "(", tau,
                                ")  (new infections per day)")),
           fill = bquote(paste('R'['0'] * ' = ' * .(format(round(R0, 1), nsmall = 1)) * ":" ))
      ) +
      coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
    scale_fill_manual(values = cols,
                       labels = c(
                         bquote(paste('R'['p'] * ' = ' * .(format(round(RP, 1), nsmall = 1)) * " from pre-symptomatic")),
                         bquote(paste('R'['s'] * ' = ' * .(format(round(RS, 1), nsmall = 1)) * " from symptomatic")),
                         bquote(paste('R'['e'] * ' = ' * .(format(round(RE, 1), nsmall = 1)) * " from environmental")),
                         bquote(paste('R'['a'] * ' = ' * .(format(round(RA, 1), nsmall = 1)) * " from asymptomatic"))))
    #scale_x_continuous(minor_breaks = minor_breaks, breaks = breaks)
    
    grid.arrange(mainPlot, ncol=1)
  })
  
  # create four mini plots, decomposing the main plot into each category
  getDecompositionPlot <- reactive({
    validate(
      need(try(df.beta.p <- getDFBetaP()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.s <- getDFBetaS()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.e <- getDFBetaE()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    validate(
      need(try(df.beta.a <- getDFBetaA()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    
    ymax <- getYmax()
    cols <- getCols()
    RS <- getRs()
    RP <- getRp()
    RE <- getRe()
    RA <- getRa()
    R0 <- getR0()
    
    p <- ggplot(df.beta.p, aes(x=tau, y=beta)) +
      theme_bw(base_size = 18) +
      theme(legend.position="none") +
      geom_area(aes(fill=label)) +
      labs(x = expression(paste(tau, " (days)")),
           y = expression(paste("contribution to ", beta, "(", tau,")"))
      ) +
      coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
      scale_fill_manual(values = cols[[1]])
    
    s <- ggplot(df.beta.s, aes(x=tau, y=beta)) +
      theme_bw(base_size = 18) +
      theme(legend.position="none") +
      geom_area(aes(fill=label)) +
      labs(x = expression(paste(tau, " (days)")),
           y = ""
      ) +
      coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
      scale_fill_manual(values = cols[[2]])
    
    e <- ggplot(df.beta.e, aes(x=tau, y=beta)) +
      theme_bw(base_size = 18) +
      theme(legend.position="none") +
      geom_area(aes(fill=label)) +
      labs(x = expression(paste(tau, " (days)")),
           y = ""
      ) +
      coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
      scale_fill_manual(values = cols[[3]])
    
    a <- ggplot(df.beta.a, aes(x=tau, y=beta)) +
      theme_bw(base_size = 18) +
      theme(legend.position="none") +
      geom_area(aes(fill=label)) +
      labs(x = expression(paste(tau, " (days)")),
           y = ""
      ) +
      coord_cartesian(xlim = c(0, max(tau.test)), ylim = c(0, ymax), expand = F) +
      scale_fill_manual(values = cols[[4]])
                          
                          
    grid.arrange(p,s,e,a, ncol=4)
  })
  
  # get the parameter summary to be written below the plot
  getParameterSummary <- reactive({
    validate(
      need(try(theta <- getTheta()), "Parameters are too extreme for the integral to converge, please adjust one or more of the sliders")
    )
    
    HTML(paste0("<h4>
                <br/>
                <br/> Proportion of infections which are not from direct contact with symptomatic individuals &Theta; = ", round(theta, 2),
                "</h4>")) # using html "&theta;" to get symbol theta
  })
  
  ### Make mini plot to demonstrate the user-defined incubation period distribution with
  # the mean and sd they have specified. 
  getIncperDistributionPlot <- reactive({
    req(getIncperMeanlog()) # can't plot it until these values have been calculated
    req(getIncperSdlog())
    
    m <- getIncperMeanlog()
    s <- getIncperSdlog()
    cols <- getCols()
    
    # Make sure we capture at least 99% of the distribution in the plot
    xmax <- 20 #max(qlnorm(0.99, meanlog = m, sdlog = s), 20)
    x <- seq(0,xmax,by=0.01)
    y <- sapply(x, function(x) plnorm(x,m,s))
    
    ggplot(cbind.data.frame(x,y), aes(x=x, y=y)) + 
      geom_line(col=cols[[5]]) +
      theme_classic() +
      coord_cartesian(expand = F, xlim = c(0, xmax), ylim = c(0, 1.05 * max(y))) + 
      labs(x = "Incubation period (days)", y = "Probability density")
  })
  
  # mini plot for generation time distribution
  getSerintDistributionPlot <- reactive({
    req(getSerIntShape()) # can't plot it until these values have been calculated
    req(getSerIntScale())
      
    m <- getSerIntShape()
    s <- getSerIntScale()
    cols <- getCols()
    
    # Make sure we capture at least 99% of the distribution in the plot.
    xmax <- 20 # max(qlnorm(0.99, meanlog = m, sdlog = s), 20)
    x <- seq(0,xmax,by=0.01)
    y <- sapply(x, function(x) dweibull(x,m,s))
      
    ggplot(cbind.data.frame(x,y), aes(x=x, y=y)) + 
      geom_line(col=cols[[5]]) +
      theme_classic() + 
      coord_cartesian(expand = F, xlim = c(0, xmax), ylim = c(0, 1.05 * max(y))) + 
      labs(x = "Generation time (days)", y = "Probability density")
  })
  
  ###############################################
  # Reactive sidebar input / outputs
  ###############################################
  
  getEnvSliders <- reactive({
    type <- input$env.type
    if (type=="exp.decay") {
      sliderInput("envDecayRate",
                  h6("Environmental infectiousness exponential decay rate (per day):"),
                  min = 0.05,
                  max = 5,
                  step = 0.01,
                  value = 2.3)
    } else {
      sliderInput("envConstant",
                  h6("Environmental infectiousness duration (days):"),
                  min = 1,
                  max = 10,
                  step = 0.01,
                  value = 3)
    }
    
  })


  
  ###########
  # OUTPUTS
  ##########
  
  # main plot
  output$mainPlot <- renderPlot({
    validate(need(getPlot(), ""))
  })
  
  # decomposed version of plot
  output$decompositionPlot <- renderPlot({
    validate(need(getDecompositionPlot(), ""))
  })
  
  # parameter summary below plot
  output$parameterSummary  <- renderUI({
    validate(need(getParameterSummary(), ""))
    getParameterSummary()
  })
  
  # SIDEBAR:
  
  # mini plot to demonstrate incper distribution as chosen by user
  output$IncperDistribPlot <- renderPlot({
    getIncperDistributionPlot()
  })
  
  # mini plot to demonstrate serint distribution as chosen by user
  output$SerintDistribPlot <- renderPlot({
    getSerintDistributionPlot()
  })
  
  # environment sliders, depending on type
  output$environment_sliders <- renderUI({
    getEnvSliders()
  })
  
  #####################################################################################################
  
  ###################
  # CONTOUR PLOT TAB
  ###################
  
  # define the reset button
  observeEvent(input$reset_inputContour, { 
    shinyjs::reset("sidePanel2")
    shinyjs::reset("mainPanel2")
  })
  
  # Largely a repeat of everything above, with "Contour" appended! Plus delay slider(s)
  getDelay <- reactive({
    input$delay
  })
  
  getIncperMedianlogContour <- reactive({
    log(input$incperMedianContour)
  })
  
  getIncperSdlogContour <- reactive({
    input$incperSdlogContour
  })
  
  getSerIntShapeContour <- reactive({
    input$serIntShapeContour
  })  
  
  getSerIntScaleContour <- reactive({
    input$serIntScaleContour
  })
  
  getDoublingTimeContour <- reactive({
    input$doublingTimeContour
  })
  
  getXpContour <- reactive({
    1
  })
  
  getXaContour <- reactive({
    input$xaContour # the relative infectiousness of asymptomatic c.f. symptomatic individuals
  })
  
  getP.aContour <- reactive({
    input$P.aContour # the fraction of all infections that are asymptomatic
  })
  
  getFrac.ReContour <- reactive({
    input$frac.ReContour # the fraction of all transmissions that are environmentally mediated
  })
  
  
  ### Make mini plot to demonstrate the user-defined incubation period distribution with
  # the mean and sd they have specified. 
  getIncperDistributionPlotContour <- reactive({
    req(getIncperMedianlogContour()) # can't plot it until these values have been calculated
    req(getIncperSdlogContour())
    
    m <- getIncperMedianlogContour()
    s <- getIncperSdlogContour()
    cols <- getCols()
    
    # Make sure we capture at least 99% of the distribution in the plot
    xmax <- 20 #max(qlnorm(0.99, meanlog = m, sdlog = s), 20)
    x <- seq(0,xmax,by=0.01)
    y <- sapply(x, function(x) plnorm(x,m,s))
    
    ggplot(cbind.data.frame(x,y), aes(x=x, y=y)) + 
      geom_line(col=cols[[5]]) +
      theme_classic() +
      coord_cartesian(expand = F, xlim = c(0, xmax), ylim = c(0, 1.05 * max(y))) + 
      labs(x = "Incubation period (days)", y = "Probability density")
  })
  
  # mini plot for generation time distribution
  getSerintDistributionPlotContour <- reactive({
    req(getSerIntShapeContour()) # can't plot it until these values have been calculated
    req(getSerIntScaleContour())
    
    m <- getSerIntShapeContour()
    s <- getSerIntScaleContour()
    cols <- getCols()
    
    # Make sure we capture at least 99% of the distribution in the plot.
    xmax <- 20 # max(qlnorm(0.99, meanlog = m, sdlog = s), 20)
    x <- seq(0,xmax,by=0.01)
    y <- sapply(x, function(x) dweibull(x,m,s))
    
    ggplot(cbind.data.frame(x,y), aes(x=x, y=y)) + 
      geom_line(col=cols[[5]]) +
      theme_classic() + 
      coord_cartesian(expand = F, xlim = c(0, xmax), ylim = c(0, 1.05 * max(y))) + 
      labs(x = "Generation time (days)", y = "Probability density")
  })
  
  getFa <- reactive({
    P.a <- getP.aContour()
    x.a <- getXaContour()
    
    fa = P.a*x.a / (P.a*x.a + (1 - P.a) )
  })
  
  
  # make interactive sliders for environmental decay type
  getEnvSlidersContour <- reactive({
    type <- input$env.typeContour
    if (type=="exp.decay") {
      sliderInput("envDecayRateContour",
                  h6("Environmental infectiousness exponential decay rate (per day):"),
                  min = 0.05,
                  max = 5,
                  step = 0.01,
                  value = 2.3)
    } else {
      sliderInput("envConstantContour",
                  h6("Environmental infectiousness duration (days):"),
                  min = 1,
                  max = 10,
                  step = 0.01,
                  value = 3)
    }
    
  })
  
  
  getEnvInfTypeContour <- reactive({
    input$env.typeContour
  })
  
  getEnvDecayContour <- reactive({
    input$envDecayRateContour
  })
  
  getEnvConstContour <- reactive({
    input$envConstantContour
  })
  
  # following "model.gen.beta.env.div.by.E.RSorP"
  getModelGenBetaEnvDivByERSorPContour <- reactive({
    incper.meanlog <- getIncperMedianlogContour()
    incper.sdlog <- getIncperSdlogContour()
    serint.scale <- getSerIntScaleContour()
    serint.shape <- getSerIntShapeContour()
    P.a <- getP.aContour()
    xp <- getXpContour()
    env.decay.rate <- getEnvDecayContour()
    env.constant.duration <- getEnvConstContour()
    env.infectiousness.type <- getEnvInfTypeContour()
    
    # because of the way that "env.decay.rate" and "env.constant.duration" depend on
    # "env.infectiousness.type", there is a delay when shiny is launched, when these are undefined
    # so we use "req" to wait until they are defined before continuing the calculation
    if (env.infectiousness.type == "constant") {
      req(env.constant.duration)
    } else { 
      req(env.decay.rate)
    }
    
    Vectorize(function(tau, 
                       incper.meanlog, 
                       incper.sdlog, 
                       serint.scale,
                       serint.shape, 
                       P.a, 
                       xp, 
                       env.decay.rate,
                       env.constant.duration, 
                       env.infectiousness.type) {
      integrate(function(l) {
        model.gen.beta.s.div.by.RSorP(tau = tau - l,
                                      incper.meanlog = incper.meanlog,
                                      incper.sdlog = incper.sdlog,
                                      serint.scale = serint.scale,
                                      serint.shape = serint.shape,
                                      P.a = P.a,
                                      xp = xp) *
          p.el(l = l, 
               env.decay.rate = env.decay.rate,
               env.constant.duration = env.constant.duration,
               env.infectiousness.type = env.infectiousness.type)
      },
      lower = 0, upper = tau)$value
    }, vectorize.args = "tau")
    
    
  })
  
  # following model.gen.full.beta.div.by.RSorP
  getModelGenFullBetaDivByRSorPContour <- reactive({
    serint.scale <- getSerIntScaleContour()
    serint.shape <- getSerIntShapeContour()
    P.a <- getP.aContour()
    xa <- getXaContour()
    incper.meanlog <- getIncperMedianlogContour()
    incper.sdlog <- getIncperSdlogContour()
    xp <- getXpContour()
    env.decay.rate <- getEnvDecayContour()
    env.constant.duration <- getEnvConstContour()
    env.infectiousness.type <- getEnvInfTypeContour()
    
    model.gen.beta.env.div.by.E.RSorP <- getModelGenBetaEnvDivByERSorPContour()
    
    function(tau, serint.scale, serint.shape, P.a, xa, incper.meanlog,
             incper.sdlog, xp, env.scale.constant, env.decay.rate,
             env.constant.duration, env.infectiousness.type) {
      serint(x = tau, 
             serint.scale = serint.scale,
             serint.shape = serint.shape) *
        (1 + P.a * xa / model.gen.f(tau = tau, 
                                    incper.meanlog = incper.meanlog,
                                    incper.sdlog = incper.sdlog, 
                                    P.a = P.a, 
                                    xp = xp)) +
        env.scale.constant * # this is a free variable of this function
        model.gen.beta.env.div.by.E.RSorP(tau = tau, 
                                          serint.scale = serint.scale,
                                          serint.shape = serint.shape, 
                                          incper.meanlog = incper.meanlog,
                                          incper.sdlog = incper.sdlog, 
                                          P.a = P.a, 
                                          xp = xp, 
                                          env.decay.rate = env.decay.rate,
                                          env.constant.duration = env.constant.duration,
                                          env.infectiousness.type = env.infectiousness.type)
    }
  })
  
  
  # following "getDummy()" from the infectiousness tab, but slimmed down because we only need to get R0 out
  getR0Contour <- reactive({
    # get parameters
    frac.Re <- getFrac.ReContour()
    P.a <- getP.aContour()
    doubling.time <- getDoublingTimeContour()
    xp <- getXpContour()
    xa <- getXaContour()
    incper.meanlog <- getIncperMedianlogContour()
    incper.sdlog <- getIncperSdlogContour()
    serint.shape <- getSerIntShapeContour()
    serint.scale <- getSerIntScaleContour()
    theta.obs <- 0.83
    env.decay.rate <- getEnvDecayContour()
    env.constant.duration <- getEnvConstContour()
    env.infectiousness.type <- getEnvInfTypeContour() 
    
    # get reactive functions 
    model.gen.beta.env.div.by.E.RSorP <- getModelGenBetaEnvDivByERSorPContour()
    model.gen.full.beta.div.by.RSorP <- getModelGenFullBetaDivByRSorPContour()
    
    # use "model.gen.solve" in a reactive way:
    r <- log(2) / doubling.time # units of per day. This is r before interventions
    
    integral.of.model.gen.beta.s.div.by.RSorP <-
      integrate(model.gen.beta.s.div.by.RSorP, 
                lower = 0, 
                upper = Inf,
                incper.meanlog = incper.meanlog, 
                incper.sdlog = incper.sdlog,
                serint.scale = serint.scale, 
                serint.shape = serint.shape,
                P.a = P.a, 
                xp = xp)$value
    
    env.scale.constant <- (1 + P.a * xa * integral.of.model.gen.beta.s.div.by.RSorP) / 
      (((1 / frac.Re) - 1) * integrate(model.gen.beta.env.div.by.E.RSorP,
                                       lower = 0, 
                                       upper = Inf,
                                       serint.scale = serint.scale,
                                       serint.shape = serint.shape,
                                       incper.meanlog = incper.meanlog,
                                       incper.sdlog = incper.sdlog,
                                       P.a = P.a,
                                       xp = xp,
                                       env.decay.rate = env.decay.rate,
                                       env.constant.duration = env.constant.duration,
                                       env.infectiousness.type = env.infectiousness.type)$value)
    
    
    RSorP <- 1 / integrate(function(tau) {
      model.gen.full.beta.div.by.RSorP(tau = tau, 
                                       serint.scale = serint.scale,
                                       serint.shape = serint.shape,
                                       P.a = P.a, 
                                       xa = xa, 
                                       incper.meanlog = incper.meanlog,
                                       incper.sdlog = incper.sdlog, 
                                       xp = xp,
                                       env.scale.constant = env.scale.constant,
                                       env.decay.rate = env.decay.rate,
                                       env.constant.duration = env.constant.duration,
                                       env.infectiousness.type = env.infectiousness.type) *
        exp(-r * tau)}, lower = 0, upper = Inf)$value
    
    R0 <- RSorP * integrate(model.gen.full.beta.div.by.RSorP, 
                            lower = 0, 
                            upper = Inf,
                            serint.scale = serint.scale,
                            serint.shape = serint.shape,
                            P.a = P.a, 
                            xa = xa, 
                            incper.meanlog = incper.meanlog,
                            incper.sdlog = incper.sdlog, 
                            xp = xp,
                            env.scale.constant = env.scale.constant,
                            env.decay.rate = env.decay.rate,
                            env.constant.duration = env.constant.duration,
                            env.infectiousness.type = env.infectiousness.type)$value
    
    R0
  })
  
  getM <- reactive({
    fa <- getFa()
    fe <- getFrac.ReContour()
    validate(need(try(R0 <- getR0Contour()), "")) # this is not immediately available because of the reactive "environment" sliders
    delay <- getDelay()
    
    # incubation from Lauer et al
    log_incubation_sd <- getIncperSdlogContour()
    log_incubation_median <- getIncperMedianlogContour()
    
    S<-function(x){1-(1-fa)*plnorm(x-delay/24,meanlog = log_incubation_median, sdlog = log_incubation_sd)}
    
    # Weibull generation time
    beta<-function(x){
      dweibull(x, shape = getSerIntShapeContour(), scale = getSerIntScaleContour())
    }
    # auxiliary functions 
    M1<-function(x){R0*beta(x)*(1-ei+ei*S(x))}
    M2<-function(x,y){(1-et+et*S(x+y)/S(y))}
    
    # nmax, n, ndiscr are defined in global.R
   
    # initialization
    Y<-rep(1,nmax)/nmax
    r<-0
    
    m<-matrix(0,nrow=n-1,ncol=n-1)
    for(ci in 1:(n-1)){
      for(cj in 1:(n-1)){
        ei<-ci/n ; et<-cj/n*(1-fe)
        eigen<-function(my_r){
          r<-my_r
          Y<-rep(1,nmax)/nmax;
          Yold<-Y
          # for(i in 1:maxiter){ # removed in favour of the three lines below, suggested by Chris for speed-up
          #   Y<-M1(v)*exp(-v*r)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})
          temp.vec <- M1(v)*exp(-v*r)
          for(i in 1:maxiter){
            Y<-temp.vec * sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})
            Y<-Y/sum(Y)
            if(sum(abs(Y-Yold))<Yerror & i>=miniter){break}
            Yold<-Y
          }
          return(lm(I(M1(v)*exp(-v*r)*sapply(v,function(z){sum(M2(z,v)*Y)/ndiscr})) ~ Y)$coeff[2]-1)
        }
        m[ci,cj]<-tryCatch(uniroot(eigen,interval = c(-2,2))$root, error = function(e) {return(NA)} ) }
    }
    colnames(m)<-c(1:(n-1))/n
    rownames(m)<-c(1:(n-1))/n
    
    df <- tibble(x=numeric((n-1)^2), y=numeric((n-1)^2), z=numeric((n-1)^2))
    count <- 0
    for (i in 1:(n-1)) {
      for (j in 1:(n-1)) {
        count <- count + 1
        df[count,] <- list(x = rownames(m)[[i]],
                           x = colnames(m)[[j]],
                           z = m[i,j])
      }
    }
    
    # remove the "X" at the start of the y values, and convert to percentages
    df$y <- vapply(df$y, function(value) substring(value, 2), character(1))
    df$x <- 100 * as.numeric(df$x)
    df$y <- 100 * as.numeric(df$y)
    df
  })
  
  getContourPlot <- reactive({
    validate(need(try(df <- getM()), ""))
    
    ggplot(df, aes(x, y, fill = z)) + 
      theme_bw(base_size = 18) +
      coord_cartesian(expand = F) +
      geom_tile() +
      scale_fill_gradient2(low = "#66C2A5", 
                           mid = "white", 
                           high = "#FC8D62", 
                           midpoint = 0) +
      labs(x = "% success in isolating cases",
           y = "% success in quarantining contacts",
           fill = "growth rate r \nafter interventions") +
      stat_contour(aes(z=z), breaks=0, color = "black", size = 1) 
  })
  
  # describe the delay assumption above the plot
  getContourDescribeDelay <- reactive({
    validate(need(d <- getDelay(), ""))
    
    if (d==0) { HTML(paste0("<h4>With instant interventions:</h4>")) }
    else if (d==1) {HTML(paste0("<h4>With interventions after a delay of ",d," hour:</h4>")) } # ("hour" singular)
    else {HTML(paste0("<h4>With interventions after a delay of ",d," hours:</h4>")) }
  })
  
  # describe the delay assumption above the plot
  getContourDescribeR0 <- reactive({
    validate(need(R0 <- getR0Contour(), ""))
    r <- log(2) / getDoublingTimeContour()
    
    HTML(paste0("<h4>The current choice of input parameters describes an epidemic 
                where the reproduction number R0 = ",round(R0,1)," and the 
                growth rate r = ",round(r,2)," <b>before</b> interventions are applied.
                The colours in the plot indicate what the growth rate would be if interventions were applied
                with a range of success rates: red corresponds to a growing epidemic, green to a declining epidemic.
                The black line shows the combinations of success rates of interventions which are needed to achieve r = 0, the threshold for epidemic control.
                </h4>"))
  })

  ##########
  # OUTPUTS
  ##########
  
  # mini plot to demonstrate incper distribution as chosen by user
  output$IncperDistribPlotContour <- renderPlot({
    getIncperDistributionPlotContour()
  })
  
  # mini plot to demonstrate serint distribution as chosen by user
  output$SerintDistribPlotContour <- renderPlot({
    getSerintDistributionPlotContour()
  })
  
  # environment sliders, depending on type
  output$environment_slidersContour <- renderUI({
    getEnvSlidersContour()
  })
  
  output$contourPlot <- renderPlot({
    getContourPlot()
  })
  
  output$contourDescribeDelay <- renderUI({
    getContourDescribeDelay()
  })
  
  output$contourDescribeR0 <- renderUI({
    getContourDescribeR0()
  })
  
}

