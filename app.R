library(shiny)
library(shinythemes)
library(DT)
library(shinyWidgets)
library(plotly)
source("cross-species phenotypes-disease.R")
Do_disease <- read.delim("data/Do_disease.txt",header = TRUE,sep=",")
ui <- fluidPage(
  navbarPage(theme = shinytheme("flatly"),title = h2("Cell-Phenotype"),
             tabPanel(h3("Cross-species Phenotype-Disease"),fluidRow(
               column(1),
               column(11,h2("Explore the cell-phenotype-disease association ")),
               column(1),
               column(4, wellPanel(
                 selectizeInput('disease', label = h4("Enter DO or MONDO disease:"), choices = c(Do_disease$x),
                selected = "Charcot-Marie-Tooth disease type 1A",options = list(create = TRUE)),
                h5("Note: Only diseases connected with both mammalian phenotypes and human phenotypes are on the choice list") )
                 ),
               column(7,
                      br(),
                      h4("The connection between mammalian phenotypes and human diseases comes from Human - Mouse: Disease Connection (HMDC)"),
                  h4("The connection between human phenotypes and human diseases comes from Mondo Disease Ontology (Mondo)")
               )),
               fluidRow(
                 column(2),
                 column(10,h4("Cell-Phenotype-Disease Network:"),
                        h6("the lightblue boxes are denoted as MP, the lightgreen boxes are denoted as HP")),
                 visNetworkOutput("network_cross_do", height = "690px"),
                 plotlyOutput("pheno_cell_plot",height='100%'),
                fluidPage(DTOutput('table_do1')),
                fluidPage(DTOutput('table_do2'))
               ))
  )
    
)

server <- function(input, output) {

  output$network_cross_do <- renderVisNetwork({
    net_cross_do(disease = input$disease)
  })
  output$pheno_cell_plot <- plotly::renderPlotly(plot_pheno_cell(disease = input$disease))
  output$table_do1 <- renderDT(table_do_m(disease = input$disease), 
                              rownames = FALSE,extensions = 'Buttons', 
                              options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
  output$table_do2 <- renderDT(table_do_h(disease = input$disease), 
                               rownames = FALSE,extensions = 'Buttons', 
                               options = list(dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))
}

# Run the application 
shinyApp(ui = ui, server = server)
