library(deSolve)
library(ggplot2)
library(shiny)
library(plotly)
# library(gtools)
library(reshape2)
pres_ready <- theme(axis.title = element_text(size=16), axis.text = element_text(size=16), legend.text = element_text(size=16), title = element_text(size=20),legend.key.height = unit(1.0, 'cm'))  

library(reshape2)



### parameters for the models



parameter_genpop <- c("beta_A", "beta_I", "sigma", "gamma_A", "gamma_I", "nu_H", "nu_A", "nu_I")
parameter_def_genpop <- c(2, 2, 0.1, .1, .1, 0, 0, 0.0)


parameter_hcw <- c("beta_AH", "beta_IH","beta_AG", "beta_IG",  "sigma_H", "gamma_AH", "gamma_IH", "nu_H", "nu_AH", "nu_IH")
parameter_def_hcw <- c(2, 2, 2, 2, 0.1, 0.1, 0.1, 0, 0, 0.0)

initial_val <- c("S_g" = 14.57e6,  A_g = 1e4, I_g = 1e4, R_g = 0, S_h = 5e3, A_h = 50, I_h = 50, R_h = 0)
init_N <- sum(initial_val[grep(names(initial_val), pat="_g")])
init_N_h <- sum(initial_val[grep(names(initial_val), pat="_h")])
initial_val <- c(initial_val, c("N" = init_N, "N_h" = init_N_h))



ui <- fluidPage(
    titlePanel("NOT READY FOR USE, WORK IN PROGRESS! Exploring local COVID Stats"),
    fluidRow(column(3, 
                    fluidRow(h4("General Population Parameters")),
                    lapply(seq_along(parameter_genpop), function(i) {
                        fluidRow(numericInput(inputId = paste0(parameter_genpop[i]), label = paste(parameter_genpop[i]),
                                              value = parameter_def_genpop[i]))
                    })
    ), 
    column(3, 
           fluidRow(h4("HCW Parameters")),
           lapply(seq_along(parameter_hcw), function(i) {
               fluidRow(numericInput(inputId = paste0(parameter_hcw[i]), label = paste(parameter_hcw[i]),
                                     value = parameter_def_hcw[i]))
           })
    ),
    column(6, fluidRow(plotlyOutput("HCWPlot")), fluidRow(plotlyOutput("GenPopPlot")))
    ),
    fluidRow(column(6, sliderInput("days", "Number of Days to project", 7, 365, 24, round = TRUE)))
    # uiOutput("parameter_input")
    # uiOutput("projectionPlot")
)

server <- function(input, output) {
    
    
    
    
    model_def <- function(time, state, parameters) {
        par <- as.list(c(state, parameters))
        with(par, {
            ## General Pop Equations
            dS_g <- -beta_A*(I_g+A_g)*S_g/N - beta_I*(I_g+A_g)*S_g/N - nu_H*S_g
            dA_g <- beta_A*(I_g+A_g)*S_g/N - sigma*A_g - gamma_A*A_g - nu_A*A_g
            dI_g <- beta_I*(I_g+A_g)*S_g/N + sigma*A_g - gamma_I*I_g - nu_I*I_g
            dR_g <- gamma_I*I_g + gamma_A*A_g - nu_H*R_g
            dN <-  - nu_H*S_g - nu_A*A_g - nu_I*I_g - nu_H*R_g
            ## HCW equations
            dS_h <- -beta_AH*(A_h)*S_h/N_h - beta_IH*(A_h)*S_h/N_h - beta_AG*(I_g+A_g)*S_h/(N) - beta_IG*(I_g+A_g)*S_h/(N) - nu_H*S_h ## TODO:: What is the right normalization factor for these equations for beta?
            dA_h <- beta_AH*(A_h)*S_h/N_h + beta_AG*(I_g+A_g)*S_h/(N) - sigma_H*A_h - gamma_AH*A_h - nu_AH*A_h
            dI_h <- beta_IH*(A_h)*S_h/N_h + beta_IG*(I_g+A_g)*S_h/(N) + sigma_H*A_h - gamma_IH*I_h - nu_IH*I_h
            dR_h <- gamma_IH*I_h + gamma_AH*A_h - nu_H*R_h  
            dN_h <- - nu_AH*A_h - nu_IH*I_h - nu_H*R_h - nu_H*S_h
            list(c(dS_g, dA_g, dI_g, dR_g, dS_h, dA_h, dI_h, dR_h, dN, dN_h)) # Order matches order of initial values
        })
    }
    
    
    #     RSS <- function(parameters) {
    #         names(parameters) <- c("beta", "gamma")
    #         out <- ode(y = init, times = Day, func = SIR, parms = parameters)
    #         fit <- out[, 3]
    #         sum((Infected - fit)^2)
    #     }
    #     
    
    input_parameters <- reactive({
        vals <- lapply(c(parameter_genpop, parameter_hcw), function(nm) {return(input[[nm]])})
        names(vals) <- c(parameter_genpop, parameter_hcw)
        vals
    })
    duration <- reactive({
        input[["days"]]
    })
    cumInc <- reactive({
        pars <- input_parameters()
        days <- duration()
        t <- seq_len(days)
        cumulative_incidence <- data.frame(ode(y = initial_val, times = t, func = model_def, parms = pars))
        cumulative_incidence
    })
    output$HCWPlot <- renderPlotly({
        days <- duration()
        t <- seq_len(days)
        cumulative_incidence <- cumInc()
        # browser()
        t <- t[seq_along(cumulative_incidence[,"I_h"])]
        # toPlot <- data.frame("Date" = t, "Predicted Symptomatic HCW" = cumulative_incidence[,"I_h"], check.names=FALSE)
        # 
        # ggplot(toPlot, aes(x=Date, y=`Predicted Symptomatic HCW`)) + geom_line() + pres_ready
        toPlot <- cbind("Date" = t, as.data.frame(cumulative_incidence))
        
        # ggplot(toPlot) + geom_line(aes(t, I_g), color="red") + geom_line(aes(t, A_g), color="cyan") + geom_line(aes(t, S_g), color="black") + geom_line(aes(t, R_g), color="green") + pres_ready
        colors <- c("Infected HCW Symptomatic" = "red", "Infected HCW Asymptomatic" = "cyan", "Recovered HCW" = "green", "Susceptible HCW" = "black")
        
        ggplot(toPlot) + geom_line(aes(t, I_h, color="Infected HCW Symptomatic")) + geom_line(aes(t, A_h, color="Infected HCW Asymptomatic")) + 
                         geom_line(aes(t, S_h, color="Susceptible HCW")) + 
                         geom_line(aes(t, R_h, color="Recovered HCW")) + pres_ready + labs(x = "Days", y = "People", color = "Legend") +
                         scale_color_manual(values = colors)
        
    })
    output$GenPopPlot <- renderPlotly({
        days <- duration()
        t <- seq_len(days)
        cumulative_incidence <- cumInc()
        # browser()
        t <- t[seq_along(cumulative_incidence[,"I_h"])]
        # toPlot <- data.frame("Date" = t, "Predicted Symptomatic HCW" = cumulative_incidence[,"I_h"], check.names=FALSE)
        # 
        # ggplot(toPlot, aes(x=Date, y=`Predicted Symptomatic HCW`)) + geom_line() + pres_ready
        toPlot <- cbind("Date" = t, as.data.frame(cumulative_incidence))
        
        # ggplot(toPlot) + geom_line(aes(t, I_g), color="red") + geom_line(aes(t, A_g), color="cyan") + geom_line(aes(t, S_g), color="black") + geom_line(aes(t, R_g), color="green") + pres_ready
        colors <- c("Infected GP Symptomatic" = "red", "Infected GP Asymptomatic" = "cyan", "Recovered GP" = "green", "Susceptible GP" = "black")
        
        ggplot(toPlot) + geom_line(aes(t, I_g, color="Infected GP Symptomatic")) + geom_line(aes(t, A_g, color="Infected GP Asymptomatic")) + 
            geom_line(aes(t, S_g, color="Susceptible GP")) + 
            geom_line(aes(t, R_g, color="Recovered GP")) + pres_ready + labs(x = "Days", y = "People", color = "Legend") +
            scale_color_manual(values = colors)
        
    })

    # 
    # output$parameter_input <- renderUI(
    #     )

    
    # 
    # output$mainPage <- renderUI(
    #     fluidRow(fluidRow(column(6,selectInput("whichRegion", "Region", c("Toronto", "Ontario", "Canada")))),
    #              tabsetPanel(
    #                 tabPanel("SIR Model for local cases only", 
    #                      fluidRow(column(12,plotOutput("currentPlot"))), fluidRow(column(3, h3("Reproduction Rate: ")), column(6,h3(textOutput("reproductionRateSIR")))),
    #                      fluidRow(column(12,plotOutput("projectionPlot"))),
    #                      fluidRow(column(12, paste0("Data from: https://art-bd.shinyapps.io/covid19canada/, work still in progress. Modelling with simple SIR model now. Dates cut to after:",format(cutOff.date))))
    #                 ),
    #             tabPanel("EarlyR model using SI estimates", 
    #                     fluidRow(column(12, plotOutput("plotEarlyR"))),
    #                     fluidRow(column(12, plotOutput("plotProjectionEarlyR"))),
    #                     fluidRow(column(12, paste0("Using estimates of SI from Li et al. 2020 of: µ=", mu_si," sigma=", sd_si)))))
    # ))
    #     
    # 
    
}

shinyApp(ui = ui, server = server)

#column(6, sliderInput("daysProject", "Number of Days to project", 7, 365, 24, round = TRUE))
