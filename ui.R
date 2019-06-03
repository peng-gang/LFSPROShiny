library(shiny)

shinyUI(
  fluidPage(
    includeCSS("style.css"),
    # Application title
    titlePanel("Predict future cancer risk for families with Li-Fraumeni syndrome"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        p(em("TP53"), "germline mutations are the main cause of Li-Fraumeni Syndrome. 
          This application is designed to estimate probabilities that: 
          1) the counselee is a ", em("TP53"), "germline mutation carrier, 
          2) the counselee develops any cancer in future, 
          3) the counselee develops breast cancer, sarcoma or any other cancers in future, 
          4) the counselee develops a first or second primary cancer in future, on the basis of his/her family cancer history. 
          The package also provides functions for using the LFS classic and Chompret criteria."),
        br(),
        
        p("The online tool requires two files as input: Family Data File and Cancer Data File.
  Please check the two links below to download the example files:"),
        p(a("Family Data File", href="fam.data.csv", download="fam.data.csv")),
        p(a("Cancer Data File", href="cancer.data.csv", download="cancer.data.csv")),
        
        br(),
        p("After running the LFSPRO, there will be a figure and a table shown on the right. 
  The figure shows the pedigree structure of the input family. 
  The family member with cancer is denoted by solid circle/sqare. 
  The table shows the results of LFRPRO, Classic and Chompret(2015) criteria.",
          span("The cancer history will be shown by clicking the sample in the table.", style = "color:blue")),
        
        br(),
        fileInput("file1", "Choose Family Data File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        fileInput("file2", "Choose Cancer Data File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        tags$hr(),
        uiOutput('ui.action'), # instead of conditionalPanel
        uiOutput("ui.cutoff")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        fluidRow(plotOutput("distPlot")),
        fluidRow(verbatimTextOutput("sInfo")),
        fluidRow(DT::dataTableOutput("table")),
        fluidRow(column(12, align="center",
                        uiOutput('ui.download')))
      )
    )
  )
)