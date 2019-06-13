library(shiny)
library(kinship2)
library(LFSPRO)
library(DT)
library(data.table)

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
      th(rowspan = 2, 'Cancer Type'),
      th(colspan = 3, 'Cancer Risk')
    ),
    tr(
      lapply(c("5 Years", "10 Years", "15 Years"), th)
    )
  )
))


LFSPRO.rlt <- NULL
cutoff <- 0.2
lfs.mode <- NULL

shinyServer(function(input, output) {
  
  buttonInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  
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
  
  output$ui.action <- renderUI({
    if (is.null(famdata())) return()
    if (is.null(cancerdata())) return()
    actionButton("action", "Run LFSPRO (1st)")
  })
  
  output$ui.action2 <- renderUI({
    if (is.null(famdata())) return()
    if (is.null(cancerdata())) return()
    actionButton("action2", "Run LFSPRO (mpc)")
  })
  
  observeEvent(
    eventExpr = input$action,
    handlerExpr = {
      lfs.mode <<- "1st.cs"
      #lfs.mode <<- "mpc"
      
      output$distPlot <- renderPlot({
        if (is.null(input$action) || is.null(input$action2)) return()
        if (input$action==0 && input$action2==0) return()
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
        if (is.null(input$action) || is.null(input$action2)) return(DT::datatable(NULL))
        if (input$action==0 && input$action2==0) return(DT::datatable(NULL))
        
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
            rlt.lfs <- lfspro.mode(fam.data, cancer.data, counselee.id, lfs.mode)
            prob <- rlt.lfs[[1]][,3]
            if(lfs.mode=="1st.cs"){
              risk.lfs <- rlt.lfs[[2]][[1]][,4:6]
            } else if(lfs.mode=="mpc"){
              risk.lfs <- rlt.lfs[[2]][,3:5]
            }
            
            rlt <- data.frame(id = factor(rlt[,2], levels =  rlt[,2]),
                              ProbLFSPRO = as.numeric(prob),
                              LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                              levels = c("Yes", "No")),
                              Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                                levels = c("Yes", "No")),
                              Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                               levels = c("Yes", "No")),
                              Cancer5y = risk.lfs[,1],
                              Cancer10y = risk.lfs[,2],
                              Cancer15y = risk.lfs[,3])
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
    }
    
  )
  
  observeEvent(
    eventExpr = input$action2,
    handlerExpr = {
      lfs.mode <<- "mpc"
      
      output$distPlot <- renderPlot({
        if (is.null(input$action) || is.null(input$action2)) return()
        if (input$action==0 && input$action2==0) return()
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
        if (is.null(input$action) || is.null(input$action2)) return(DT::datatable(NULL))
        if (input$action==0 && input$action2==0) return(DT::datatable(NULL))
        
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
            
            #allef.g <- list(c(0.9999,0.0001))
            #mRate.g <- 5e-4
            
            counselee.id <- data.frame(fam.id=fam.data$fam.id, id = fam.data$id)
            rlt <- lfspro.mode(fam.data, cancer.data, counselee.id, "1st.all")
            rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
            rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
            rlt.lfs <- lfspro.mode(fam.data, cancer.data, counselee.id, lfs.mode)
            prob <- rlt.lfs[[1]][,3]
            if(lfs.mode=="1st.cs"){
              risk.lfs <- rlt.lfs[[2]][[1]][,4:6]
            } else if(lfs.mode=="mpc"){
              risk.lfs <- rlt.lfs[[2]][,3:5]
            }
            
            rlt <- data.table(id = factor(rlt[,2], levels =  rlt[,2]),
                              ProbLFSPRO = as.numeric(prob),
                              LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                              levels = c("Yes", "No")),
                              Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                                levels = c("Yes", "No")),
                              Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                               levels = c("Yes", "No")),
                              CancerType = buttonInput(
                                FUN = actionButton,
                                len = length(prob),
                                id = 'button_',
                                label = "All Cancers",
                                onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
                              ),
                              Cancer5y = risk.lfs[,1],
                              Cancer10y = risk.lfs[,2],
                              Cancer15y = risk.lfs[,3])
            colnames(rlt) <- c("ID", "ProbLFSPRO", "LFSPRO-carrier", "Chompret criteria", "Classic criteria",
                               "Cancer Type", "5 Years", "10 Years", "15 Years")
            LFSPRO.rlt <<- rlt
            rlt
          })
        }, rownames= FALSE, container = sketch, filter = 'top', escape = F) %>% DT::formatRound('ProbLFSPRO', 3) %>%
          DT::formatRound(c("5 Years", "10 Years", "15 Years"), 3) %>%
          formatStyle(c('LFSPRO-carrier', 'Chompret criteria', 'Classic criteria'),
                      color = styleEqual("Yes", 'red'))
      })
      
    }
  )
  
  output$ui.cutoff <- renderUI({
    if (is.null(input$action) || is.null(input$action2)) return()
    if (input$action==0 && input$action2==0) return()
    #print(input$action)
    #print(input$action2)
    sliderInput("cutoff", "Cutoff for probability",
                min=0, max=1, value = 0.2, step = 0.05)
  })
  
  output$download <- downloadHandler('LFSPRO.csv', content = function(file) {
    write.csv(LFSPRO.rlt, file, row.names = FALSE)
  })
  
  output$ui.download <- renderUI({
    if (is.null(input$action) || is.null(input$action2)) return()
    if (input$action==0 && input$action2==0) return()
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
  
  
  observeEvent(eventExpr = input$cutoff, handlerExpr = {
    cutoff <<- input$cutoff
    
    output$distPlot <- renderPlot({
      if (is.null(input$action) || is.null(input$action2)) return()
      if (input$action==0 && input$action2==0) return()
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
      if (is.null(input$action) || is.null(input$action2)) return(DT::datatable(NULL))
      if (input$action==0 && input$action2==0) return(DT::datatable(NULL))
      
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
          
          counselee.id <- data.frame(fam.id=fam.data$fam.id, id = fam.data$id)
          rlt <- lfspro.mode(fam.data, cancer.data, counselee.id, "1st.all")
          rlt.chompret <- lfsChompret2015(fam.data, cancer.data, counselee.id)
          rlt.classic <- lfsClassic(fam.data, cancer.data, counselee.id)
          rlt.lfs <- lfspro.mode(fam.data, cancer.data, counselee.id, lfs.mode)
          prob <- rlt.lfs[[1]][,3]
          if(lfs.mode=="1st.cs"){
            risk.lfs <- rlt.lfs[[2]][[1]][,4:6]
            rlt <- data.frame(id = factor(rlt[,2], levels =  rlt[,2]),
                              ProbLFSPRO = as.numeric(prob),
                              LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                              levels = c("Yes", "No")),
                              Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                                levels = c("Yes", "No")),
                              Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                               levels = c("Yes", "No")),
                              Cancer5y = risk.lfs[,1],
                              Cancer10y = risk.lfs[,2],
                              Cancer15y = risk.lfs[,3])
            colnames(rlt) <- c("ID", "ProbLFSPRO", "LFSPRO-carrier", "Chompret criteria", "Classic criteria",
                               "5 Years", "10 Years", "15 Years")
          } else if(lfs.mode=="mpc"){
            risk.lfs <- rlt.lfs[[2]][,3:5]
            
            rlt <- data.table(id = factor(rlt[,2], levels =  rlt[,2]),
                              ProbLFSPRO = as.numeric(prob),
                              LFSPRO = factor(ifelse(prob>cutoff, "Yes", "No"), 
                                              levels = c("Yes", "No")),
                              Chompret = factor(ifelse(rlt.chompret$result, "Yes", "No"),
                                                levels = c("Yes", "No")),
                              Classic = factor(ifelse(rlt.classic$result, "Yes", "No"),
                                               levels = c("Yes", "No")),
                              CancerType = buttonInput(
                                FUN = actionButton,
                                len = length(prob),
                                id = 'button_',
                                label = "All Cancers",
                                onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
                              ),
                              Cancer5y = risk.lfs[,1],
                              Cancer10y = risk.lfs[,2],
                              Cancer15y = risk.lfs[,3])
            colnames(rlt) <- c("ID", "ProbLFSPRO", "LFSPRO-carrier", "Chompret criteria", "Classic criteria",
                               "Cancer Type", "5 Years", "10 Years", "15 Years")
          }
          
          LFSPRO.rlt <<- rlt
          rlt
        })
      }, rownames= FALSE, container = sketch, filter = 'top', escape = F) %>% DT::formatRound('ProbLFSPRO', 3) %>%
        DT::formatRound(c("5 Years", "10 Years", "15 Years"), 3) %>%
        formatStyle(c('LFSPRO-carrier', 'Chompret criteria', 'Classic criteria'),
                    color = styleEqual("Yes", 'red'))
    })
  })
})


