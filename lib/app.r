#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(tidyr)
library(DT)
library(colourpicker)
source("nurigFunctions.R")
# Define UI for application

finaliza <- function(dir) {
  unlink(dir, recursive = TRUE, force = TRUE)
}

ui <- fluidPage(
  titlePanel(windowTitle = "NuRIG: Nucmer Ring Image Generator",{imageOutput("Logo", width = 30, height = 20)}),
  fluidRow(
    column(
      width = 4 ,
      wellPanel(
        wellPanel(
          h3("Input Data"),
          p(""),
          fileInput(
            inputId = "reference",
            multiple = FALSE,
            label = "Select Reference Genome"
          ),
          fileInput(
            inputId = "query",
            multiple = TRUE,
            label = "Select Query Genomes"
          ),
          
          actionButton(
            inputId = "button",
            label = "Apply"
          )
        ),
        uiOutput("imagePanel"),
        uiOutput("downloadPanel")
      )
    ),
    column(
      width = 8,tabsetPanel(
          tabPanel("Image",conditionalPanel(condition = "!is.null(output.figure)",imageOutput("figure", width = "auto"))),
          tabPanel("Annotation",
                   wellPanel(
                      fileInput(inputId = "annotation",
                        multiple = FALSE,
                        label = "Annotation file"
                      ),
                      uiOutput("filterPanel")
                   ),
                   dataTableOutput(outputId = "annotationTable")
          )
      )
    )
  )
)


server <- function(input, output, session) {
  datos <-
    reactiveValues(tmpPath = tempfile(),
                   currentPath = getwd(),
                   condicion = 0,
                   separador = " ")
  
  output$Logo = renderImage({
    return(list(
      src = "./www/Logo.png",
      width = 100,
      height = 100
    ))
  })
  
  
  
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  observeEvent(input$width, {
    return(NULL)
  })
  observeEvent(input$height, {
    return(NULL)
  })
  observeEvent(input$palette, {
    return(NULL)
  })
  observeEvent(input$reference, {
    return(NULL)
  })
  observeEvent(input$query, {
    return(NULL)
  })
  
  
  observeEvent(input$button, {
    dir.create(datos$tmpPath)
    setwd(datos$tmpPath)
    
    req(input$reference)
    req(input$query)
    
    withProgress(message = "Performing Nucmer", {
      datos$FinalTable <<- nurigCalc(reference = input$reference,
                                     query = as.data.frame(input$query))
    })
    
    output$imagePanel <- renderUI({
      wellPanel(
        h3("Figure Options"),
        
        selectInput(
          inputId = "palette",
          label = "Select palette color",
          choices = c("rainbow", "topo", "terrain"),
          selected = TRUE
        ),
        textInput(inputId = "height",
          label = "Image Height",
          value = 800
        ),
        textInput(inputId = "width",
          label = "Image Width",
          value = 800
        ),
        actionButton(inputId = "redraw", label = "Draw")
        
      )
    })
  })
  
  observeEvent(input$redraw, once = TRUE, {
    output$downloadPanel <- renderUI({
      wellPanel(
        h3("Dowload Options"),
        selectInput(
          inputId = "format",
          label = "Format Image",
          choices = c("SVG", "PNG", "JPG", "SVGZ"),
          selected = "SVG"
        ),
        downloadButton(outputId = "download", label = "Download Image")
      )
    })
  })
  observeEvent(input$redraw, {
    withProgress(message = "Rendering Image", {
      nurigRender(
        finalTable = datos$FinalTable,
        reference = input$reference,
        pal = input$palette,
        annotation =datos$view,
        width = input$width,
        height = input$height,
        binPath = datos$currentPath
      )
    })
    
   
    output$figure = renderImage({
      return(list(
        src = "fig.svg",
        alt = "Image",
        width = input$width,
        height = input$height
      ))
    })
  })
  output$download <- downloadHandler(
    filename = function() {
      paste("FigureDownload.", tolower(input$format), sep = "")
    },
    content <- function(file) {
      system(
        paste(
          "java -jar ",
          datos$currentPath,
          "/bin/cgview/cgview.jar -i out.xml -o figure.",
          input$format,
          " -f ",
          input$format,
          " -H ",
          input$height,
          " -W ",
          input$width,
          collapse = "",
          sep = ""
        )
      )
      
      file.copy(paste(
        "figure.",
        input$format,
        collapse = "",
        sep = ""
      ), file)
      
      
    }
  )
  
########################### Annotation TAB ######################################################################
  observeEvent(input$annotation,{
    withProgress(message = "Parsing Annotation File",{datos$table = gffImporter(input$annotation$datapath)})
    datos$table$color = "black"
    datos$view = datos$table
   
    
    output$filterPanel  <- renderUI({
      flowLayout(
        selectInput(inputId = "selectFeatureField",label ="Feature Type",choices = c(levels(as.factor(datos$table$feature))), multiple = TRUE),
        selectInput(inputId = "selectLabelField", label = "Label Field", choices = c(colnames(datos$table)[8:length(colnames(datos$table))-1]), multiple = TRUE),
        selectInput(inputId = "separador", label = "Field Separator",choices = c("Space","Tab",";",":","|")),
        checkboxInput(inputId = "removeEmpty", label = "Remove Empty cells?"),
        colourInput(inputId = "colorSelector", label ="Set Color", value = "black")
        
      )
    })
    output$annotationTable <- renderDataTable(datos$view)
  })
  
  
  observeEvent(c(input$selectFeatureField,input$selectLabelField,input$removeEmpty),{
    
    if(length(input$selectFeatureField)>0){
      datos$view = datos$table %>%filter(feature == input$selectFeatureField)
    }else{
      datos$view = datos$table
      
    }
    
    if(length(input$selectLabelField) > 1)
    {
      datos$view[is.na(datos$view)] =
      datos$view = datos$view %>% select(one_of(c(colnames(datos$table)[1:7],input$selectLabelField))) %>% unite_("label",c(input$selectLabelField), sep = datos$separador)
    }else if(length(input$selectLabelField) == 1) {
      datos$view = datos$view %>% select(one_of(c(colnames(datos$table)[1:7])), label = one_of(input$selectLabelField), color)
    }else{
      datos$view = datos$view %>% select(everything())
    }
    if(input$removeEmpty && length(input$selectLabelField) > 0)
    {
      datos$view = datos$view %>% filter(label != "" & !is.na(filter ))
      
    }
    
    
    output$annotationTable <- renderDataTable(datos$view, options = list(searching = TRUE))
    renderPrint({output$nAnnot = nrow(datos$view)})  
  })
  
  # observeEvent(input$tableId_cell_clicked,{
  #   
  #   
  # })
  
  observeEvent(input$separador,{###c("Space","Tab",";",":","|")
    if(input$separador == "Space")
    {
       datos$separador = " "
    }else if(input$separador == "Tab")
    {
      datos$separador = '\t'
    }else if(input$separador == ":")
    {
      datos$separador = ":"
    }else if(input$separador == ";")
    {
      satos$separador = ";"
    }else if (input$separador == "|")
    {
      datos$separador = "|"
    }
            
  })
 
  
 
  #observeEvent(input$applyFilters,{output$annotationTable <- renderDataTable(datos$view)})
  
  
  #session$onSessionEnded(finaliza(datos$tmpPath))
}

# Run the application

shinyApp(ui = ui, server = server)
