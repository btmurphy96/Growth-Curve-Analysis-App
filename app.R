library(shiny)
library(DT)
library(tidyverse)
library(scales)
library(ggsci)
library(ggpubr)
library(shinyjqui)

# Define UI
ui <- fluidPage(
  
  #Title of application
  titlePanel("Growth Curve Analysis"),
  
  #Layout of application 
  sidebarLayout(
    sidebarPanel(
      fileInput('target_upload', 'Choose CSV file to upload',
                accept = '.csv'),
      tags$hr(),
      h3("Plot customization"),
      h5(strong("Plot dimensions")), 
      p("Adjust plot dimensions by clicking the bottom right tab and dragging to the desired width or height."),
      selectInput("color", "Plot color selection", choices = c("plasma", "viridis", "cividis", "turbo"), selected = "viridis"),
      sliderInput("labels", "Label text size", min = 12, max = 30, step = 2, value = 20),
      sliderInput("text", "Descriptor text size", min = 12, max = 30, step = 2, value = 16),
      numericInput("point", "Graph point size", min = 1, max = 5, step = 1, value = 3),
      numericInput("line", "Graph line size", min = 0.5, max = 2.5, step = 0.5, value = 1),
      numericInput("error", "95% CI error bar line size", min = 0.2, max = 1.4, step = 0.2, value = 0.6),
      tags$hr(),
      h3("Statistics"),
      uiOutput("selectComp")
      
    ),
    mainPanel(
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", jqui_resizable(plotOutput('plot', width = '750px', height = '650px'))),
                  tabPanel("T-test", DT::dataTableOutput("stats")),
                  tabPanel("Tables", h3("Input Table"), p("Your CSV file converted to an R data table."), DT::dataTableOutput("sample_table"), 
                           tags$hr(), h3("Data Summary"), p("A summary of the CSV input data. The interaction between sample and condition is grouped by day to compute the average concentraction, standard deviation, number of replicates, standard error, and 95% confidence interval for each sample under its specific condition on each day."), DT::dataTableOutput("tidy_table"))
      )
    )
  )
)

# Define server logic
server <- shinyServer(function(input, output, session) {
  
  #Select T-test reference interaction from tidy CSV data below
  output$selectComp <- renderUI(
    selectInput("comp_sel","Reference (control) group for T-test", choices= 
                  as.character(unique(unlist(tidy_data()$interaction))))
  )
  
  # Capture CSV data in reactive data frame
  df_products_upload <- reactive({
    inFile <- input$target_upload
    if (is.null(inFile))
      return(NULL)
    df <- read_csv(inFile$datapath)
    return(df)
  })
  
  # Print CSV data to confirm proper upload
  output$sample_table<- DT::renderDataTable({
    DT::datatable(df_products_upload(), options = list(searching = FALSE),  rownames= FALSE)
  })

  # Tidy the CSV data 
  tidy_data <- reactive({
    if (is.null(df_products_upload()))
      return(NULL)
    tidy_data <- df_products_upload()%>% 
      mutate(sample = as_factor(sample), condition = as_factor(condition), interaction = factor(interaction(sample, condition, sep = " "))) %>% 
      pivot_longer(-c("sample", "condition", "interaction"), names_to = "day", names_transform = list(day = as.integer), values_to = "conc")%>% 
      group_by(day, interaction) %>% 
      summarise(conc = conc, avg = mean(conc), sd = sd(conc), n = length(conc), 
                se = sd/sqrt(n), MoE = 1.96*se)
    return(tidy_data)
  })
  
  # Print the tidy CSV data summary (replicate data required for plotting is removed to allow for more intuitive analysis)
  output$tidy_table <- DT::renderDataTable({
    if (is.null(df_products_upload()))
      return(NULL)
    tidy_table <- df_products_upload() %>% 
      mutate(sample = as_factor(sample), condition = as_factor(condition), interaction = factor(interaction(sample, condition, sep = " "))) %>% 
      pivot_longer(-c("sample", "condition", "interaction"), names_to = "day", names_transform = list(day = as.integer), values_to = "conc")%>% 
      group_by(day, interaction) %>% 
      summarise(average = mean(conc), st.dev = sd(conc), n = length(conc), 
                se = st.dev/sqrt(n), `95% CI` = 1.96*se)
    DT::datatable(tidy_table, options = list(searching = FALSE),  rownames= FALSE)
    })
  
  # Print a statistical analysis of the data 
  output$stats <- DT::renderDataTable({
    if (is.null(df_products_upload()))
      return(NULL)
    t.test <- compare_means(conc ~ interaction, tidy_data(), method = "t.test", 
                            group.by = "day", ref.group = input$comp_sel) 
    DT::datatable(t.test, options = list(searching = FALSE, lengthChange = FALSE, paging = FALSE),  rownames= FALSE)
    })
  
  # Plot the tidy CSV data
  output$plot <- renderPlot({
    if (is.null(df_products_upload()))
      return(NULL)
    
    breaks <- 10^(-10:10)
    minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))
    labels <- input$labels
    text <- input$text
    
    (ggplot(data = tidy_data()) +
      geom_line(aes(x = day, y = avg, color = interaction), size = input$line, alpha = 1/2) + 
      geom_errorbar(aes(x = day, ymin = avg - MoE, ymax = avg + MoE, color = interaction), 
                    size = input$error, width = 0.5, position=position_dodge(width=0.5)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.5), 
                 aes(x = day, y = conc, color = interaction), alpha = 1/3, 
                 size = input$point) +
      scale_y_log10(breaks = breaks, minor_breaks = minor_breaks, 
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(sides = "l") +
      labs(x = "Day", y = "Cell Concentration (cells/mL)") +
      theme_bw() + #black and white theme
      theme(panel.grid.minor.x = element_blank(), #remove minor x-axis gridlines
            panel.border = element_blank(), #remove all border
            axis.line = element_line(), #add border on both axes 
            axis.text = element_text(size = text),
            axis.title = element_text(size = labels),
            axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
            axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
      theme(legend.title = element_text(size= labels), #change legend title font size
            legend.text = element_text(size= text)) + #change legend text font size
      scale_color_viridis_d(option = input$color, "Sample"))
    })
  })

# Run the application 
shinyApp(ui = ui, server = server)

