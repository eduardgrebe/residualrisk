library(shiny)
library(shinyWidgets)
library(shinyjs)
library(shinydashboard)
library(tidyverse)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Residual risk of viral transfusion transmission (alpha version)"),
  fluidRow(
    tabsetPanel(id = "tabset", type = "tabs",
                tabPanel("Estimate infectious window period",
                         fluidRow(
                           column(2,
                                  wellPanel(
                                    h4("Viral parameters"),
                                    numericInput("copies_per_virion",
                                                 label = "Nucleic acid copies per virion",
                                                 value = 2, min = 1, max = 2, step = 1),
                                    numericInput("doubling_time_hrs",
                                                 label = "Viral doubling time (hours)",
                                                 value = 20.5,
                                                 min = 1.0, max = 120, step = 0.1),
                                    #conditionalPanel(
                                    #  condition = "input.iwp_uncertainty == true",
                                      numericInput("doubling_time_hrs_sd",
                                                   label = "Viral doubling time (hours) standard deviation",
                                                   value = 1.326531, min = 0.1, max = 12, step = 0.1)
                                    #)
                                  ),
                                  wellPanel(
                                    h4("Transfusion parameters"),
                                    sliderInput("transfused_volume", 
                                                label = "Transfused plasma volume (mL)", 
                                                value = 20, 
                                                min = 10, max = 300, step = 5),
                                    #conditionalPanel(
                                    #  condition = "input.iwp_uncertainty == true",
                                      sliderInput("transfused_volume_range", "Transfused plasma volume range:",
                                                  min = 10, max = 300, step = 5,
                                                  value = c(15,50))
                                    #)
                                  )
                           ),
                           column(2,
                                  wellPanel(
                                    h4("Screening parameters"),
                                    radioButtons("lod_units", label = "Assay LoD units",
                                                 c("copies/mL" = "copies",
                                                   "IU/mL" = "iu"
                                                 ),
                                                 selected = "iu"),
                                    conditionalPanel(
                                      condition = "input.lod_units == 'copies'",
                                      numericInput("lod95cpml",
                                                   label = "95% LoD (copies/mL)",
                                                   value = 12.32558, min = 1.0, max = 100, step = 0.1),
                                      numericInput("lod50cpml",
                                                   label = "50% LoD (copies/mL)",
                                                   value = 2.732558, min = 0.1, max = 50, step = 0.1),
                                      #conditionalPanel(
                                      #  condition = "input.iwp_uncertainty == true",
                                        numericInput("lod50cpml_sd",
                                                     label = "50% LoD standard deviation",
                                                     value = 0.1928097, min = 0.01, max = 10, step = 0.01)
                                      #)
                                    ),
                                    conditionalPanel(
                                      condition = "input.lod_units == 'iu'",
                                      numericInput("lod95iupml",
                                                   label = "95% LoD (IU/mL)",
                                                   value = 21.2, min = 1.0, max = 150, step = 0.1),
                                      numericInput("lod50iupml",
                                                   label = "50% LoD (IU/mL)",
                                                   value = 4.7, min = 0.1, max = 75, step = 0.1),
                                      #conditionalPanel(
                                      #  condition = "input.iwp_uncertainty == true",
                                        numericInput("lod50iupml_sd",
                                                     label = "50% LoD standard deviation",
                                                     value = 0.3316327, min = 0.01, max = 10, step = 0.01),
                                      #),
                                      numericInput("iu_cp_conversion_factor",
                                                   label = "Conversion factor (IUs/copy)",
                                                   value = 1.72, min = 0.10, max = 5.00, step = 0.01)
                                    ),
                                    sliderInput("pool_size",
                                                label = "Minipool size (1 for ID-NAT)",
                                                value = 16, min = 1, max = 64, step = 1),
                                    sliderInput("retests",
                                                label = "Number of retests (that could result in release)",
                                                value = 1, min = 0, max = 5, step = 1)
                                  )
                                  
                           ),
                           column(2,
                                  wellPanel(
                                    h4("Uncertainty analysis"),
                                    #checkboxInput("iwp_uncertainty", "Compute credible intervals"),
                                    #conditionalPanel(
                                    #  condition = "input.iwp_uncertainty == true",
                                      sliderInput("confidence_level", 
                                                  label = "Confidence level (%)", 
                                                  value = 95, 
                                                  min = 80, max = 99, step = 1),
                                    sliderInput("n_bs", 
                                                label = "Bootstrapping iterations", 
                                                value = 5000, 
                                                min = 5000, max = 25000, step = 5000)
                                    #),
                                    #conditionalPanel(
                                    #  condition = "input.iwp_uncertainty == true",
                                    #  actionButton("bs_iwp", "Perform bootstrapping", width = "100%")
                                    #)
                                  )
                           ),
                           column(6,
                                  wellPanel(
                                    h2("Infectious Window Period estimate", align = "center"),
                                    br(), br(),
                                    strong(textOutput("iwp_pe"), align = "center"),
                                    strong(textOutput("iwp_ci"), align = "center"),
                                    br(), br(),
                                    plotOutput("iwp_plot"),
                                    br(),br(),
                                    textOutput("written")
                                  )
                           )
                         )
                ),
                tabPanel("Estimate residual risk",
                         fluidRow(
                           column(4,
                                  wellPanel(
                                    h4("Incidence"),
                                    #conditionalPanel(
                                    #  condition = "typeof output.iwp_ci == 'undefined'",
                                    #  em("To perform uncertainty analysis, compute CIs on IWP",
                                    #     style = "color:orange;")
                                    #),
                                    numericInput("incidence",
                                                 label = "Donor incidence (cases/100,000PY)",
                                                 value = 3.123, min = 0, max = 5000, step = 0.001),
                                    numericInput("incidence_sd",
                                                 label = "Donor incidence standard deviation",
                                                 value = 0.912, min = 0, max = 5000, step = 0.001)
                                  )
                                  #conditionalPanel(
                                  #  condition = "typeof output.iwp_ci !== 'undefined'",
                                    
                                  #)
                           ),
                           column(8,
                                  wellPanel(
                                    h2("Residual risk estimate"),
                                    strong(textOutput("rr_pe"), align = "center"),
                                    strong(textOutput("rr_ci"), align = "center")
                                  )
                           )
                         )
                ),
                tabPanel("About",
                         fluidRow(
                           column(12,
                                  includeMarkdown("about.md")
                                  )
                           
                         )
                )
    )
  ),
  fluidRow(
    hr()
  ),
  fluidRow(
    hidden(p("Starting parallel computing engines...",
             id = "starting",
             style = "font-weight:bold;color:red;",
             align = "center")),
    hidden(p("Warning: fewer than 10,000 bootstrap iterations are unlikely to render meaningful CIs.",
             id = "warn_bs",
             style = "font-weight:bold;color:orange;",
             align = "center")),
    em(textOutput("n_cpu"), align = "center")
  ),
  fluidRow(
    hr()
  ),
  fluidRow(
    div(img(src="vri_logo.png", height = "100px"), style="text-align: center;")
  )
)

server <- function(input, output, session) {
  #browser()
  
  # -------- Set up Python environment (compatible with shinyapps.io) -------- #
  PYTHON_DEPENDENCIES = c('pip', 'numpy','scipy')
  
  
  # ALTON replace "alton" with your local login username, then this should run on 
  # your machine and be deployable to shinyapps.io, which requires virtualenv 
  # instead of conda environments and the installation of packages upon deploy. 
  # I've configured my machine up using .Rprofile to accept the shinyapps setup
  if (Sys.info()["user"] == "alton") {
    VIRTUALENV_NAME = "residualrisk"
    options(shiny.port = 7450)
    reticulate::use_condaenv(condaenv = VIRTUALENV_NAME, conda = "auto", required = T)
  } else if (Sys.info()["user"] == "eduard") {
    VIRTUALENV_NAME = "residualrisk"
    options(shiny.port = 7450)
    python_path = "python3"
    virtualenv_dir = paste0('~/.virtualenvs/', VIRTUALENV_NAME, '/')
    reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
    reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES)
    Sys.setenv(RETICULATE_PYTHON = paste0('~/.virtualenvs/', VIRTUALENV_NAME, '/bin/python'))
    reticulate::use_virtualenv(virtualenv_dir, required = T)
    
  } else if (Sys.info()[['user']] == 'shiny') {
    virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
    python_path = Sys.getenv('PYTHON_PATH')
    reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
    reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES)
    reticulate::use_virtualenv(virtualenv_dir, required = T)
  }
  
  
  reticulate::source_python("residualrisk.py")
  
  # virtualenv_create(envname = "residualriskapp", python= "python3") # uncomment deployment
  # virtualenv_install("residualrisk", packages = c('numpy','scipy')) # uncomment for deployment
  # use_virtualenv("residualrisk", required = TRUE)
  # source_python("residualrisk.py")
  
  load_default_sims <- function() {
    iwp_bs_default <- read_rds("iwp_bs_default.rds")
    return(iwp_bs_default)
  }
  
  check_all_defaults <- function() {
    #browser()
    if(
      input$copies_per_virion == 2 &
      input$doubling_time_hrs == 20.5 &
      input$doubling_time_hrs_sd == 1.326531 &
      input$transfused_volume == 20 &
      input$transfused_volume_range[1] == 15 &
      input$transfused_volume_range[2] == 50 &
      ( (input$lod_units == "iu" & 
         input$lod95iupml == 21.2 & 
         input$lod50iupml == 4.7 & 
         input$lod50iupml_sd == 0.3316327 & 
         input$iu_cp_conversion_factor == 1.72) | 
        (input$lod_units == "copies" &
         input$lod95cpml == 12.32558 &
         input$lod50cpml == 2.732558 &
         input$lod50cpml_sd == 0.1928097) ) &
      input$pool_size == 16 &
      input$retests == 1 &
      input$n_bs == 5000
    ) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  calculate_iwp_pe <- reactive({
    #browser()
    if (input$lod_units == "copies") {
      lod50 <- input$lod50cpml
      lod95_lod50_ratio <- input$lod95cpml / input$lod50cpml
    } else if (input$lod_units == "iu") {
      lod50 <- input$lod50iupml / input$iu_cp_conversion_factor
      lod95_lod50_ratio <- input$lod95iupml / input$lod50iupml
    }
    
    C0 <- 0.00025
    doubling_time <- input$doubling_time_hrs / 24
    k <- 0.02433568
    
    iwp_pe <- risk_days(input$copies_per_virion,
                        C0,
                        doubling_time,
                        input$transfused_volume,
                        k,
                        as.integer(input$pool_size),
                        lod50,
                        lod95_lod50_ratio,
                        as.integer(input$retests))
    return(iwp_pe)
    
  })
  
  calculate_iwp_ci <- reactive({
    if (check_all_defaults()) {

      iwp_vec <- load_default_sims()[[3]]

    } else {
      
    shinyjs::show("running")
    if (input$n_bs < 10000) {
      shinyjs::show("warn_bs")
    } else {
      shinyjs::hide("warn_bs")
    }
    
    #if (input$iwp_uncertainty) {
      if (input$lod_units == "copies") {
        lod50 <- input$lod50cpml
        lod95_lod50_ratio <- input$lod95cpml / input$lod50cpml
        lod50_sd <- input$lod50cpml_sd
      } else if (input$lod_units == "iu") {
        lod50 <- input$lod50iupml / input$iu_cp_conversion_factor
        lod95_lod50_ratio <- input$lod95iupml / input$lod50iupml
        lod50_sd <- input$lod50iupml_sd / input$iu_cp_conversion_factor
      }
      
      C0 <- 0.00025
      doubling_time <- input$doubling_time_hrs / 24
      doubling_time_sd <- input$doubling_time_hrs_sd / 24
      
      k <- 0.02433568
      k_gamma_shape <- 3.985468
      k_gamma_scale <- 0.006784984
      
      #alpha <- 1 - input$confidence_level/100
      withProgress(message = 'Performing bootstrapping', value = 0, {
        n <- 100
        seed <- 126887
        iwp_vec <- vector()
        for (i in 1:100) {
          iwp_vec <- append(iwp_vec, iwp_bs_par(k, 
                                 k_gamma_shape, 
                                 k_gamma_scale, 
                                 doubling_time, 
                                 doubling_time_sd, 
                                 lod50, 
                                 lod50_sd, 
                                 lod95_lod50_ratio, 
                                 input$transfused_volume,
                                 input$transfused_volume_range,
                                 as.integer(input$pool_size), 
                                 as.integer(input$retests), 
                                 alpha = 0.05, # Not relevant since we get quantiles outside
                                 n_bs = as.integer(input$n_bs/n),
                                 seed = as.integer(seed))[[3]])
          incProgress(1/n)
          seed <- seed+1+input$n_bs/n
        }
        
      })
    }
      return(iwp_vec)
  })
  
  output$iwp_plot <- renderPlot({
    n_bins <- floor(input$n_bs / 10)
    tibble::tibble(val = calculate_iwp_ci()) %>%
      ggplot(aes(x = val)) + 
      geom_histogram(aes(y=..density..), colour="#0072B2", fill="white", bins = n_bins) + 
      geom_density(alpha=.2, fill="#FF6666") + theme_bw() -> plot
    return(plot)
  })
  
  calculate_rr_pe <- reactive({
    rr_pe <- input$incidence / 1e5 * calculate_iwp_pe()/365.25
    return(rr_pe)
  })
  
  calculate_rr_ci <- reactive({
    if (!is.null(iwp_bs)) {
      alpha <- 1 - input$confidence_level/100
      iwp <- calculate_iwp_pe()
      iwpbs <- calculate_iwp_ci()
      rr_list <- residual_risk_iwp(iwp,
                                   iwpbs,
                                   input$incidence / 1e5,
                                   input$incidence_sd /1e5,
                                   per = 1)
      return(rr_list[[2]])
    }
  })
  
  output$iwp_pe <- renderText({
    iwp_pe <- paste0("Point estimate: ", format(round(calculate_iwp_pe(),2), nsmall=2), " days")
  })
  
  output$iwp_ci <- renderText({
    #ci <- calculate_iwp_ci()[[2]]
    iwpbs <- calculate_iwp_ci()
    alpha <- 1 - input$confidence_level/100
    ci <- quantile(iwpbs, probs = c(alpha/2, 1 - alpha/2))
    ci_text <- paste0(input$confidence_level, "% CI: ", format(round(ci[1],2), nsmall=2), "-", format(round(ci[2],2), nsmall=2))
    return(ci_text)
  })
  
  output$rr_pe <- renderText({
    rr <- calculate_rr_pe()
    rr_pe_text <- paste0("Point estimate: ", 
                         format(round(rr*1e6, 2), nsmall=2), 
                         " transmissions/million transfusions (1:", 
                         format(round(1/rr), big.mark = ","), ").")
    return(rr_pe_text)
  })
  
  output$rr_ci <- renderText({
    #browser()
    ci <- calculate_rr_ci()
    rr_ci_text <- paste0(
      input$confidence_level, "% CI: ", 
      format(round(ci[1]*1e6,2), nsmall=2), "-", format(round(ci[2]*1e6,2), nsmall=2),
      " transmission/million transfusions (1:",
      format(round(1/ci[1]), big.mark=","),
      "-1:",
      format(round(1/ci[2]), big.mark=","), ").")
    return(rr_ci_text)
  })
  
  output$n_cpu <- renderText({
    return(paste0("Number of CPUs: ", count_cores()))
  })
  
  # write_iwp_bs <- reactive({
  #   iwp_bootstraps <- calculate_iwp_ci()
  #   write_rds(iwp_bootstraps, "iwp_bs_default.rds")
  # })
  # 
  # output$written <- renderText({
  #   write_iwp_bs()
  #   return("Written.")
  # })
  
}

shinyApp(ui, server)

