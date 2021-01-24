library(shiny)
library(kinship2)
library(LFSPRO)
library(DT)
library(data.table)
library(shinyjs)
library(shinyalert)
library(reshape2)
library(ggplot2)
library(ggsci)

source("functions.R")

#Style
sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 3, ' '),
      th(rowspan = 3, 'Details'),
      th(rowspan = 3, 'ID'),
      th(rowspan = 3, "Info"),
      th(rowspan = 3, 'ProbLFSPRO'),
      th(rowspan = 3, 'LFSPRO-carrier'),
      th(rowspan = 3, 'Chompret criteria'),
      th(rowspan = 3, 'Classic criteria'),
      th(colspan = 13, 'Cancer Risk')
    ),
    tr(
      th(colspan = 3, "Breast Cancer"),
      th(colspan = 3, "Sarcoma"),
      th(colspan = 3, "Other Cancers"),
      th(colspan = 3, "Second Primary Cancer"),
      th(rowspan = 2, "Figure")
    ),
    tr(
      lapply(rep(c('5 years', '10 years', '15 years'), 4), th)
    )
  )
))


cutoff <- 0.2

shinyServer(function(input, output) {
  LFSPRO.rlt <- NULL
  myValue <- reactiveValues(idx.button = '')
  
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
  
  cid <- reactive({
    cid.tmp <- input$txtSampleId
    cid <- NULL
    if(length(cid.tmp)>0){
      cid <- strsplit(cid.tmp, ";")[[1]]
    }
    cid
  })
  
  output$ui.action <- renderUI({
    if (is.null(famdata())) return()
    if (is.null(cancerdata())) return()
    actionButton("action", "Run LFSPRO")
  })
  
  observeEvent(
    eventExpr = input$action,
    handlerExpr = {
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
          
          aff <- data.frame(Cancer = fam.data$id %in% cancer.data$id,
                            stringsAsFactors = FALSE)
          
          ped <- pedigree(id =  fam.data$id, 
                          dadid = fam.data$fid,
                          momid = fam.data$mid,
                          sex = ifelse(fam.data$gender==0, 2, 1),
                          famid = rep("fam", nrow(fam.data)),
                          status = ifelse(fam.data$vital == "A", 0, 1),
                          affected = as.matrix(aff))
          plot(ped['fam'], col = ifelse(fam.data$proband == "Y", "red", "black"))
          pedigree.legend(ped['fam'], location="topright",radius=.15)
        })
      })
      
      output$table <- DT::renderDataTable({
        if (is.null(input$action)) return(DT::datatable(NULL))
        if (input$action==0) return(DT::datatable(NULL))
        
        DT::datatable({
          isolate({
            fam.data <- famdata()
            if (is.null(fam.data)) return(NULL)
            cancer.data <- cancerdata()
            if (is.null(cancer.data)) return(NULL)
            
            fam.data$fam.id <- "fam"
            cancer.data$fam.id <- "fam"
            
            cid <- cid()
            #print(cid)
            if(is.null(cid) || length(cid)==0){
              counselee.id <- data.frame(fam.id=fam.data$fam.id, id = fam.data$id)
              #print(counselee.id)
            } else {
              idx.cid <- cid %in% fam.data$id
              if(sum(!idx.cid) > 0){
                cid.notfound <- cid[!idx.cid]
                cid <- cid[idx.cid]
                msg <- paste0("The following id(s) are not found in the input data: ",
                              paste(cid.notfound, collapse = ", "), ".")
                if(sum(idx.cid)==0){
                  msg <- paste0(msg, " All samples in the input file are selected for the calculation. ")
                } else {
                  msg <- paste0(msg, " Only the following samples are selected for the calculation: ",
                                paste(cid, collapse = ", "), ".")
                }
                shinyalert("ID Not Found", msg, type = "warning")
              }
              counselee.id <- data.frame(fam.id = "fam", id = cid)
            }
            
            ## patients are removed with age >= 80 or dead
            id.rm <- fam.data$id[fam.data$age >= 80 | fam.data$vital == "D"]
            counselee.id <- counselee.id[!(counselee.id$id %in% id.rm),]
            
            info <- NULL
            for(id in counselee.id$id){
              idx.sel <- which(cancer.data$id == id)
              if(length(idx.sel)){
                info.tmp <- NULL
                for(i in idx.sel){
                  info.tmp <- paste0(info.tmp, cancer.data$cancer.type[i], " at age ", cancer.data$diag.age[i], "; ")
                }
                info <- c(info, info.tmp)
              } else {
                info <- c(info, "No Cancer")
              }
            }
            
            rltTmp <- runLFSPRO(fam.data, cancer.data, counselee.id)
            rltTmp$info <- info
            LFSPRO.rlt <<- rltTmp
            
            rlt <- data.frame(
              id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
              info = info, 
              ProbLFSPRO = LFSPRO.rlt$carrier,
              LFSPRO = factor(ifelse(LFSPRO.rlt$carrier>cutoff, "Yes", "No"), 
                              levels = c("Yes", "No")),
              Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                                levels = c("Yes", "No")),
              Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                               levels = c("Yes", "No")),
              breast.5 = LFSPRO.rlt$breast.5,
              breast.10 = LFSPRO.rlt$breast.10,
              breast.15 = LFSPRO.rlt$breast.15,
              sarcoma.5 = LFSPRO.rlt$sarcoma.5,
              sarcoma.10 = LFSPRO.rlt$sarcoma.10,
              sarcoma.15 = LFSPRO.rlt$sarcoma.15,
              other.5 = LFSPRO.rlt$other.5,
              other.10 = LFSPRO.rlt$other.10,
              other.15 = LFSPRO.rlt$other.15,
              second.5 = LFSPRO.rlt$second.5,
              second.10 = LFSPRO.rlt$second.10,
              second.15 = LFSPRO.rlt$second.15,
              #nc.5 = LFSPRO.rlt$nc.5,
              #nc.10 = LFSPRO.rlt$nc.10,
              #nc.15 = LFSPRO.rlt$nc.15,
              figure = buttonInput(
                FUN = actionButton,
                len = nrow(LFSPRO.rlt),
                id = "button_",
                label = "Risk Trend",
                onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
              )
              )
            
            rlt <- cbind('Details' = '&oplus;', rlt)
            rlt
          })
        }, container = sketch, filter = 'top', escape = FALSE,
        options = list(
          columnDefs = list(
            list(visible = FALSE, targets = c(0, 3)),
            list(orderable = FALSE, className = 'details-control', targets = 1)
          )
        ),
        callback = JS(
        "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
        )) %>% 
          DT::formatRound('ProbLFSPRO', 2) %>%
          DT::formatRound(
            c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
              "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
          formatStyle(c('LFSPRO', 'Chompret', 'Classic'),
                      color = styleEqual("Yes", 'red'))
      })
    }
  )
  
  output$ui.cutoff <- renderUI({
    if (is.null(input$action) ) return()
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
  
  observeEvent(eventExpr = input$cutoff, handlerExpr = {
    cutoff <<- input$cutoff

    output$table <- DT::renderDataTable({
      if (is.null(input$action)) return(DT::datatable(NULL))
      if (input$action==0) return(DT::datatable(NULL))
      
      DT::datatable({
        isolate({
          rlt <- data.frame(
            id = factor(LFSPRO.rlt$id, levels =  LFSPRO.rlt$id),
            info = LFSPRO.rlt$info, 
            ProbLFSPRO = LFSPRO.rlt$carrier,
            LFSPRO = factor(ifelse(LFSPRO.rlt$carrier>cutoff, "Yes", "No"), 
                            levels = c("Yes", "No")),
            Chompret = factor(ifelse(LFSPRO.rlt$chompret, "Yes", "No"),
                              levels = c("Yes", "No")),
            Classic = factor(ifelse(LFSPRO.rlt$classic, "Yes", "No"),
                             levels = c("Yes", "No")),
            breast.5 = LFSPRO.rlt$breast.5,
            breast.10 = LFSPRO.rlt$breast.10,
            breast.15 = LFSPRO.rlt$breast.15,
            sarcoma.5 = LFSPRO.rlt$sarcoma.5,
            sarcoma.10 = LFSPRO.rlt$sarcoma.10,
            sarcoma.15 = LFSPRO.rlt$sarcoma.15,
            other.5 = LFSPRO.rlt$other.5,
            other.10 = LFSPRO.rlt$other.10,
            other.15 = LFSPRO.rlt$other.15,
            second.5 = LFSPRO.rlt$second.5,
            second.10 = LFSPRO.rlt$second.10,
            second.15 = LFSPRO.rlt$second.15,
            #nc.5 = LFSPRO.rlt$nc.5,
            #nc.10 = LFSPRO.rlt$nc.10,
            #nc.15 = LFSPRO.rlt$nc.15,
            figure = buttonInput(
              FUN = actionButton,
              len = nrow(LFSPRO.rlt),
              id = "button_",
              label = "Risk Trend",
              onclick = 'Shiny.onInputChange(\"lastClick\",  this.id)'
            )
          )
          
          rlt <- cbind('Details' = '&oplus;', rlt)
          rlt
        })
      }, container = sketch, filter = 'top', escape = FALSE,
      options = list(
        columnDefs = list(
          list(visible = FALSE, targets = c(0, 3)),
          list(orderable = FALSE, className = 'details-control', targets = 1)
        )
      ),
      callback = JS(
        "
        table.column(1).nodes().to$().css({cursor: 'pointer'});
        
        var format = function(d) {
          return '<div style=\"background-color:#eee; padding: .5em;\"> ' +
            d[3] + '</div>';
        };
        
        table.on('click', 'td.details-control', function() {
            var td = $(this), row = table.row(td.closest('tr'));
            if (row.child.isShown()) {
            row.child.hide();
            td.html('&oplus;');
            } else {
            row.child(format(row.data())).show();
            td.html('&CircleMinus;');
            }
        });
        "
      )) %>% 
        DT::formatRound('ProbLFSPRO', 2) %>%
        DT::formatRound(
          c("breast.5", "breast.10", "breast.15", "sarcoma.5", "sarcoma.10", "sarcoma.15", 
            "other.5", "other.10", "other.15", "second.5", "second.10", "second.15"), 2) %>%
        formatStyle(c('LFSPRO', 'Chompret', 'Classic'),
                    color = styleEqual("Yes", 'red'))
    })
  })
  
  
  observeEvent(eventExpr = input$lastClick, handlerExpr = {
    myValue$idx.button <<- as.numeric(strsplit(input$lastClick, "_")[[1]][2])
    
    showModal(modalDialog(
      title =  "Cancer Risk",
      plotOutput("cancerrisk"),
      fluidRow(
        column(
          6,
          verbatimTextOutput("sampleInfo")
        ),
        column(
          6,
          verbatimTextOutput("cancerInfo")
        )
      )
    ))
  })
  
  output$sampleInfo <- renderText({
    idx.button <- myValue$idx.button
    id.sel <- LFSPRO.rlt$id[idx.button]
    fam.data <- famdata()
    idx.sel <- which(fam.data$id == id.sel)
    rlt <- paste0(
      "Sample id: ", id.sel,  "\n",
      ifelse(
        fam.data$vital[idx.sel] == "A",
        paste0("Age: ", fam.data$age[idx.sel], "\n"),
        paste0("Died at age of: ", fam.data$age[idx.sel], "\n")
      ),
      "Sex: ", ifelse(fam.data$gender[idx.sel] == 0, "Female", "Male"), "\n")
    rlt
  })
  
  output$cancerInfo <- renderText({
    idx.button <- myValue$idx.button
    id.sel <- LFSPRO.rlt$id[idx.button]
    cancer.data <- cancerdata()
    idx.sel <- which(cancer.data$id == id.sel)
    rlt <- paste0("Sample id: ", id.sel, "\n")
    if(length(idx.sel)==0){
      rlt <- paste0(rlt, "No Cancer")
    } else {
      for(i in 1:length(idx.sel)){
        rlt <- paste0(rlt, cancer.data$cancer.type[idx.sel[i]], "\t", cancer.data$diag.age[idx.sel[i]], "\n")
      }
    }
    rlt
  })
  
  output$cancerrisk <- renderPlot({
    idx.button <- myValue$idx.button
    dplot <- data.frame(
      year = c(5, 10, 15),
      Breast = c(LFSPRO.rlt$breast.5[idx.button], LFSPRO.rlt$breast.10[idx.button], LFSPRO.rlt$breast.15[idx.button]),
      Sarcoma = c(LFSPRO.rlt$sarcoma.5[idx.button], LFSPRO.rlt$sarcoma.10[idx.button], LFSPRO.rlt$sarcoma.15[idx.button]),
      Other = c(LFSPRO.rlt$other.5[idx.button], LFSPRO.rlt$other.10[idx.button], LFSPRO.rlt$other.15[idx.button]),
      Second = c(LFSPRO.rlt$second.5[idx.button], LFSPRO.rlt$second.10[idx.button], LFSPRO.rlt$second.15[idx.button]),
      NonCarrier = c(LFSPRO.rlt$nc.5[idx.button], LFSPRO.rlt$nc.10[idx.button], LFSPRO.rlt$nc.15[idx.button]),
      
      stringsAsFactors = FALSE
    )
    idx.rm.col <- colSums(is.na(dplot))==0
    dplot <- dplot[,idx.rm.col]
    if(is.null(ncol(dplot))){
      return(NULL)
    }
    
    if(ncol(dplot)==1){
      return(NULL)
    }
    
    dplot.2 <- reshape2::melt(dplot, id.vars = "year", variable.name = "type", value.name = "risk")
    
    gp <- ggplot(dplot.2, aes(x=year, y = risk, fill = type)) + 
      geom_bar(stat="identity", position=position_dodge(), width = 2) + 
      theme_light() + 
      labs(x="Year", y="Cancer Risk") + 
      scale_x_continuous(breaks = c(5, 10, 15)) + 
      scale_fill_npg() + 
      theme(legend.title = element_blank(), legend.position = "bottom") + 
      theme(text = element_text(size=18))
    
    gp
    
  })
  
})


