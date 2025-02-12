Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

server <- function(input, output, session){
  options(shiny.maxRequestSize=1000000*1024^2)

  #Hide tabs
  hideTab("navbar", target = "output_panel")
  
  #################################################################################
  
  # No Example
  
  ################################################################################
  
  #################################################################################
  
  # Read data
  
  ################################################################################

  observeEvent(input$startAnalysis,{
    
    # Data path
    path <- reactive({
      input$uploadcsv$datapath
    })
    
    # Load data
    data_all <- reactive({
      
      # Require data input
      req(path())
      
      # Show loading message
      showModal(modalDialog(title = h4(strong("Reading data..."),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      # Message #1
      message("Start data reading")
      
      # Read data
      data_all <- fread(path())
      
      # Format data
      data_all <- as.data.frame(data_all)
      rownames(data_all) <- data_all[,1]
      data_all <- data_all[,-1]
      
      # Message #2
      message("End data reading")
      
      # Return output
      return(data_all)
      
    })
    
    #################################################################################
    
    # Make predictions
    
    ################################################################################
    
    # Run MPS models
    predictedScore_factors <- reactive({
      
      # data_all is required
      req(data_all())
      
      showModal(modalDialog(title = h4(strong("Calculation of methylation profile scores..."),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      # Message #3
      message("Start MPS calculation")
      
      # EXTEND models
      extend_out <- extendModel(data = data_all(), 
                                finalModels = finalModels,
                                all_cpgs =  all_cpgs)
      
      # Marioni models
      marioni_out <- marioniModel(data = data_all(), cpgs = cpgs)
      
      
      # Combine MRSs in a dataframe
      predictedScore_factors <- cbind.data.frame(extend_out,
                                                 marioni_out[,c("EpiAge",
                                                                "Alcohol",
                                                                "BMI", 
                                                                "HDL",
                                                                "Smoking")])
      rownames(predictedScore_factors) <- rownames(extend_out)
      
      # Message #4
      message("End MPS calculation")
      
      # Return output
      return(predictedScore_factors)
    })
    
    # Run MMRS-MCI model
    predictDF <- reactive({
      
      # Require predictedScore_factors
      req(predictedScore_factors())
      
      
      showModal(modalDialog(title = h4(strong("Calculation of multivariate methylation risk scores..."),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      # Message #6
      message("Start MMRS calculation")
      
      # Predict
      pred_RF <- predict(fit, predictedScore_factors(), type = "prob")
      
      # Combine into data frame
      predictDF <- data.frame(SampleID = rownames(predictedScore_factors()), 
                              MMRS = log(pred_RF$MCI/(1-pred_RF$MCI)))
      
      # Message #6
      message("End MMRS calculation")
      
      # Return output
      return(predictDF)
    })
    
    # Go the next step if previous steps are completed
    observeEvent(if (length(predictDF()) > 0){input$startAnalysis}, {
      
      removeModal()
      showTab("navbar", target = "output_panel")
      updateTabsetPanel(session, "navbar",
                        selected = "output_panel")
    })
    
    #################################################################################
    
    # Make output
    
    ################################################################################
    
    observe({
      # Output table predictedScore_factors
      output$predictedScore_factors <- DT::renderDataTable({
        
        req(predictedScore_factors())
        
        # Format table
        factors <- colnames(predictedScore_factors())
        output_table <- predictedScore_factors()
        output_table$SampleID <- rownames(output_table)
        output_table <- output_table[,c("SampleID", factors)]
        
        # Return table
        return(output_table)
      }, server=TRUE,
      options = list(pageLength = 10), rownames= FALSE)
      
      
      # Download table predictedScore_factors
      output$download_predictedScore_factors <- downloadHandler(
        filename = "predictedScore_MPSs.csv",
        content = function(file){
          factors <- colnames(predictedScore_factors())
          output_table <- predictedScore_factors()
          output_table$SampleID <- rownames(output_table)
          output_table <- output_table[,c("SampleID", factors)]
          write.table(output_table,file,row.names = FALSE,sep = ",",quote = FALSE)
        }
      )
      
      # Output table MMRS
      output$predictDF <- DT::renderDataTable({
        
        req(predictDF())
        
        # Format table
        output_table <- predictDF()
        output_table$` ` <- " "
        
        # Return table
        return(output_table)
        
      }, server=TRUE,
      options = list(pageLength = 10), rownames= FALSE)
      
      # Download table epi-MCI
      output$download_predictDF <- downloadHandler(
        filename = "predictedScore_MMRS.csv",
        content = function(file){
          write.table(predictDF(), file, row.names = FALSE, sep = ",",quote = FALSE)
        }
      )
      
    })
  })
  
  #################################################################################
  
  # Example
  
  ################################################################################
  
  
  #################################################################################
  
  # Read data
  
  ################################################################################
  
  observeEvent(input$example, {
    
    
    # Load data
   
      data_all <- reactive({
        
        # Show loading message
        showModal(modalDialog(title = h4(strong("Reading data..."),
                                         align = "center"), 
                              footer = NULL,
                              h5("This might take a while. Please be patient.", 
                                 align = "center")))
        
        # Message #1
        message("Start data reading")
        
        # Read data
        data_all <- fread("example.csv")
        
        # Format data
        data_all <- as.data.frame(data_all)
        rownames(data_all) <- data_all[,1]
        data_all <- data_all[,-1]
        
        # Message #2
        message("End data reading")
        
        # Return output
        return(data_all)
        
      })
    
    
    
    
    #################################################################################
    
    # Make predictions
    
    ################################################################################
    
    # Run MPS models
    predictedScore_factors <- reactive({
      
      # data_all is required
      req(data_all())
      
      showModal(modalDialog(title = h4(strong("Calculation of methylation profile scores..."),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      # Message #3
      message("Start MPS calculation")
      
      # EXTEND models
      extend_out <- extendModel(data = data_all(), 
                                finalModels = finalModels,
                                all_cpgs =  all_cpgs)
      
      # Marioni models
      marioni_out <- marioniModel(data = data_all(), cpgs = cpgs)
      
      
      # Combine MRSs in a dataframe
      predictedScore_factors <- cbind.data.frame(extend_out,
                                                 marioni_out[,c("EpiAge",
                                                                "Alcohol",
                                                                "BMI", 
                                                                "HDL",
                                                                "Smoking")])
      rownames(predictedScore_factors) <- rownames(extend_out)
      
      # Message #4
      message("End MPS calculation")
      
      # Return output
      return(predictedScore_factors)
    })
    
    # Run MMRS-MCI model
    predictDF <- reactive({
      
      # Require predictedScore_factors
      req(predictedScore_factors())
      
      showModal(modalDialog(title = h4(strong("Calculation of multivariate methylation risk scores..."),
                                       align = "center"), 
                            footer = NULL,
                            h5("This might take a while. Please be patient.", 
                               align = "center")))
      
      # Message #6
      message("Start MMRS calculation")
      
      # Predict
      pred_RF <- predict(fit, predictedScore_factors(), type = "prob")
      
      # Combine into data frame
      predictDF <- data.frame(SampleID = rownames(predictedScore_factors()), 
                              MMRS = log(pred_RF$MCI/(1-pred_RF$MCI)))
      
      # Message #6
      message("End MMRS calculation")

      # Return output
      return(predictDF)
    })
  
    
    # Go the next step if previous steps are completed
    observeEvent(if (length(predictDF()) > 0){input$example}, {
      
      removeModal()
      showTab("navbar", target = "output_panel")
      updateTabsetPanel(session, "navbar",
                        selected = "output_panel")
      
    })
    
    #################################################################################
    
    # Make output
    
    ################################################################################
    
    observe({
      # Output table predictedScore_factors
      output$predictedScore_factors <- DT::renderDataTable({
        
        req(predictedScore_factors())
        
        # Format table
        factors <- colnames(predictedScore_factors())
        output_table <- predictedScore_factors()
        output_table$SampleID <- rownames(output_table)
        output_table <- output_table[,c("SampleID", factors)]
        
        # Return table
        return(output_table)
      }, server=TRUE,
      options = list(pageLength = 10), rownames= FALSE)
      
      
      # Download table predictedScore_factors
      output$download_predictedScore_factors <- downloadHandler(
        filename = "predictedScore_MPSs.csv",
        content = function(file){
          write.table(predictedScore_factors(),file,row.names = FALSE,sep = ",",quote = FALSE)
        }
      )
      
      # Output table MMRS
      output$predictDF <- DT::renderDataTable({
        
        req(predictDF())
        
        # Format table
        output_table <- predictDF()
        output_table$` ` <- " "
        
        # Return table
        return(output_table)
        
      }, server=TRUE,
      options = list(pageLength = 10), rownames= FALSE)
      
      # Download table epi-MCI
      output$download_predictDF <- downloadHandler(
        filename = "predictedScore_MMRS.csv",
        content = function(file){
          write.table(predictDF(), file, row.names = FALSE, sep = ",",quote = FALSE)
        }
      )
      
    })
  })
  
  
}

