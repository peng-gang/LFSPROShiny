library(shiny)
library(kinship2)
library(LFSPRO)
library(DT)

source("functions.R")

#Style
sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, 'ID'),
      th(rowspan = 2, 'ProbLFSPRO'),
      th(rowspan = 2, 'LFSPRO-carrier'),
      th(rowspan = 2, 'Chompret criteria'),
      th(rowspan = 2, 'Classic criteria'),
      th(colspan = 3, 'Cancer Risk')
    ),
    tr(
      lapply(c("5 Years", "10 Years", "15 Years"), th)
    )
  )
))


LFSPRO.rlt <- NULL
cutoff <- 0.2

shinyServer(function(input, output) {
  famdata <- reactive({
    infile <- input$file1
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = TRUE,
             sep = ",",
             stringsAsFactors = FALSE)
  })
  
  cancerdata <- reactive({
    infile <- input$file2
    if (is.null(infile)){
      return(NULL)      
    }
    read.csv(infile$datapath,
             header = TRUE,
             sep = ",",
             stringsAsFactors = FALSE)
  })
  
  # cutoff <- reactive({
  #   cutoff <- NULL
  #   if(is.null(input$cutoff)){
  #     cutoff <- 0.2
  #   } else {
  #     cutoff <- input$cutoff
  #   }
  #   cutoff
  # })
  
  output$ui.action <- renderUI({
    if (is.null(famdata())) return()
    if (is.null(cancerdata())) return()
    actionButton("action", "Run LFSPRO")
  })
  
  output$ui.cutoff <- renderUI({
    if (is.null(input$action)) return()
    if (input$action==0) return()
    sliderInput("cutoff", "Cutoff for probability",
                min=0, max=1, value = 0.2, step = 0.05)
  })
  
  output$download <- downloadHandler('LFSPRO.csv', content = function(file) {
    write.csv(LFSPRO.rlt, file, row.names = FALSE)
  })
  
  output$ui.download <- renderUI({
    if (is.null(input$action)) return()
    if (input$action==0) return()
    downloadButton('download', "Download Results")
  })
  
  output$sInfo <- renderPrint({
    rlt <- NULL
    s = input$table_rows_selected
    if(length(s)){
      s <- sort(s)
      fam.data <- famdata()
      cancer.data <- cancerdata()
      rlt <- 'Selected Sample Information:\n\n'
      for(idx in s){
        id.sel <- fam.data$id[idx]
        idx.sel <- which(cancer.data$id==id.sel)
        if(length(idx.sel)){
          rlt <- paste0(rlt, "Sample ", fam.data$id[idx], " has the following cancer(s):\n")
          for(i in idx.sel){
            rlt <- paste0(rlt, cancer.data$cancer.type[i], " at age ", cancer.data$diag.age[i], "\n")
          }
          rlt <- paste0(rlt, "\n")
        } else {
          rlt <- paste0(rlt, "Sample ", fam.data$id[idx], " has no cancer at age ", fam.data$age[idx], ".\n\n")
        }
      }
      cat(rlt)
    }
  })
  
  output$distPlot <- renderPlot({
    if (is.null(input$action)) return()
    if (input$action==0) return()
    isolate({
      fam.data <- famdata()
      #fam.data <- LFSPRO::fam.data
      if (is.null(fam.data)) return(NULL)
      cancer.data <- cancerdata()
      #cancer.data <- LFSPRO::cancer.data
      if (is.null(cancer.data)) return(NULL)
      
      aff <- fam.data$id %in% cancer.data$id
      
      ped <- pedigree(id =  fam.data$id, 
                      dadid = fam.data$fid,
                      momid = fam.data$mid,
                      sex = ifelse(fam.data$gender==0, 2, 1),
                      famid = rep("fam", nrow(fam.data)),
                      affected = aff)
      plot(ped['fam'])
    })
  })
  
  output$table <- DT::renderDataTable({
    if (is.null(input$action)) return(DT::datatable(NULL))
    if (input$action==0) return(DT::datatable(NULL))
    
    DT::datatable({
      isolate({
        fam.data <- famdata()
        #print(famdata)
        #fam.data <- LFSPRO::fam.data
        if (is.null(fam.data)) return(NULL)
        cancer.data <- cancerdata()
        #print(cancer.data)
        #cancer.data <- LFSPRO::cancer.data
        if (is.null(cancer.data)) return(NULL)
        
        fam.data$fam.id <- "fam"
        cancer.data$fam.id <- "fam"
        
        allef.g <- list(c(0.9999,0.0001))
        mRate.g <- 5e-4
        
        counselee.id <- data.frame(fam.id=fam.data$fam.id, id = fam.data$id)
        rlt <- lfspro.mode(fam.data, cancer.data, counselee.id, "1st.all")
        rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
        rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
        rlt.mpc <- lfspro.mode(fam.data, cancer.data, counselee.id, "mpc")
        prob <- rlt[,3]
        risk.mpc <- rlt.mpc[[2]][,3:5]
        rlt <- data.frame(id = factor(rlt[,2], levels =  rlt[,2]),
                          ProbLFSPRO = as.numeric(prob),
                          LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                          levels = c("Yes", "No")),
                          Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                            levels = c("Yes", "No")),
                          Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                           levels = c("Yes", "No")),
                          Cancer5y = risk.mpc[,1],
                          Cancer10y = risk.mpc[,2],
                          Cancer15y = risk.mpc[,3])
        colnames(rlt) <- c("ID", "ProbLFSPRO", "LFSPRO-carrier", "Chompret criteria", "Classic criteria",
                           "5 Years", "10 Years", "15 Years")
        LFSPRO.rlt <<- rlt
        rlt
      })
    }, rownames= FALSE, container = sketch, filter = 'top') %>% DT::formatRound('ProbLFSPRO', 3) %>%
      DT::formatRound(c("5 Years", "10 Years", "15 Years"), 3) %>%
      formatStyle(c('LFSPRO-carrier', 'Chompret criteria', 'Classic criteria'),
                  color = styleEqual("Yes", 'red'))
  })
  
  observeEvent(eventExpr = input$cutoff, handlerExpr = {
    cutoff <<- input$cutoff
    
    output$distPlot <- renderPlot({
      if (is.null(input$action)) return()
      if (input$action==0) return()
      isolate({
        fam.data <- famdata()
        #fam.data <- LFSPRO::fam.data
        if (is.null(fam.data)) return(NULL)
        cancer.data <- cancerdata()
        #cancer.data <- LFSPRO::cancer.data
        if (is.null(cancer.data)) return(NULL)
        
        aff <- fam.data$id %in% cancer.data$id
        
        ped <- pedigree(id =  fam.data$id, 
                        dadid = fam.data$fid,
                        momid = fam.data$mid,
                        sex = ifelse(fam.data$gender==0, 2, 1),
                        famid = rep("fam", nrow(fam.data)),
                        affected = aff)
        plot(ped['fam'])
        #pedigree.legend( ped, location="topright",radius=.2)
      })
    })
    
    
    output$table <- DT::renderDataTable({
      if (is.null(input$action)) return(DT::datatable(NULL))
      if (input$action==0) return(DT::datatable(NULL))
      
      DT::datatable({
        isolate({
          fam.data <- famdata()
          #print(famdata)
          #fam.data <- LFSPRO::fam.data
          if (is.null(fam.data)) return(NULL)
          cancer.data <- cancerdata()
          #print(cancer.data)
          #cancer.data <- LFSPRO::cancer.data
          if (is.null(cancer.data)) return(NULL)
          
          fam.data$fam.id <- "fam"
          cancer.data$fam.id <- "fam"
          
          allef.g <- list(c(0.9999,0.0001))
          mRate.g <- 5e-4
          
          counselee.id <- data.frame(fam.id=fam.data$fam.id, id = fam.data$id)
          rlt <- lfspro.mode(fam.data, cancer.data, counselee.id, "1st.all")
          rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
          rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
          rlt.mpc <- lfspro.mode(fam.data, cancer.data, counselee.id, "mpc")
          prob <- rlt[,3]
          risk.mpc <- rlt.mpc[[2]][,3:5]
          rlt <- data.frame(id = factor(rlt[,2], levels =  rlt[,2]),
                            ProbLFSPRO = as.numeric(prob),
                            LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                            levels = c("Yes", "No")),
                            Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                              levels = c("Yes", "No")),
                            Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                             levels = c("Yes", "No")),
                            Cancer5y = risk.mpc[,1],
                            Cancer10y = risk.mpc[,2],
                            Cancer15y = risk.mpc[,3])
          colnames(rlt) <- c("ID", "ProbLFSPRO", "LFSPRO-carrier", "Chompret criteria", "Classic criteria",
                             "5 Years", "10 Years", "15 Years")
          LFSPRO.rlt <<- rlt
          rlt
        })
      }, rownames= FALSE, container = sketch, filter = 'top') %>% DT::formatRound('ProbLFSPRO', 3) %>%
        DT::formatRound(c("5 Years", "10 Years", "15 Years"), 3) %>%
        formatStyle(c('LFSPRO-carrier', 'Chompret criteria', 'Classic criteria'),
                    color = styleEqual("Yes", 'red'))
    })
  })
})


