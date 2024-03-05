library(shiny)
library(shinydashboard)
library(ComplexHeatmap)
library(ggplot2)
library(DT)

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(
    title = "Cell Line Selector (working title)",
    titleWidth = 200,
    tags$li(class = "dropdown", tags$a(href = "https://github.com/CaioMussatto", icon("github"), "My github"), target = "_blank"),
    tags$li(class = "dropdown", tags$a(href = "https://www.linkedin.com/in/caiomussatto", icon("linkedin"), "My Profile"), target = "_blank"),
    tags$li(class = "dropdown", tags$a(href = "https://orcid.org/0000-0003-4850-0293", icon("orcid"), "My Profile"), target = "_blank")
  ),
  dashboardSidebar(
    sidebarMenu(
      id = "sidebar",
      menuItem("Correlation Analisys", tabName = "data", icon = icon("database")),
      selectizeInput(
        inputId = "human",
        label = "Select Genes Associated with Human Cancer:",
        choices = names(lineplot),
        selected = "Lung"
      ),
      selectizeInput(
        inputId = "cell",
        label = "Choose Cell Line Dataset:",
        choices = names(lineplot[[1]]),
        selected = "LUNG"
      ),
      selectizeInput(
        inputId = "gene_set",
        label = "Select the gene set to plot the Heatmap",
        choices = names(Genes),
        selected = "All"
      )
    ),
    sidebarMenu(
      menuItem("Enrichr Analisys", tabName = "enrichr", icon = icon("chart-gantt")),
      selectInput(
        inputId = "enrichr",
        label = "Please, select the DEGS to use",
        choices = names(vias)
      ),
      selectInput(
        inputId = "enrichr2",
        label = "Please, select the pathways to use:",
        choices = dbs,
        selected = "KEGG_2019_Human"
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        "data",
        fluidRow(
          box(
            title = "Heatmap", width = 5, solidHeader = TRUE, status = "primary",
            plotOutput("heatmap", height = "900px", width = "575px")
          ),
          box(
            title = "Select the Cell lines to display", width = 5, solidHeader = TRUE, status = "primary",
            selectizeInput(
              inputId = "filter",
              label = "Select Lines to Display",
              choices = names(lineplot),
              selected = names(lineplot),
              multiple = TRUE,
              options = list(plugins = list("remove_button"))
            ),
            actionButton("filter_button", "Filter Table")
          ),
          box(
            title = "Table", width = 5, solidHeader = TRUE, status = "primary",
            column(
              width = 6,
              downloadButton("download_table", "Download Table")
            ),
            dataTableOutput("tabela")
          ),
          tags$style("
            .content-wrapper, .right-side {
                overflow-x: auto;
            }
            .content {
                min-width:1500px;
            }
        ")
        )
      ),
      tabItem(
        tabName = "enrichr",
        fluidRow(
          box(
            title = "Enrichr plot", width = 6, solidHeader = TRUE, status = "primary",
            plotOutput("enrichr_plot", height = "550px", width = "600px")
          ),
          box(
            title = "Dataframe used", width = 5, solidHeader = T, status = "primary",
            column(
              width = 6,
              downloadButton("download_enrichr", "Download Table")
            ),
            dataTableOutput("tabela_enrichr")
          ),
          box(
            title = "Genes used", width = 5, solidHeader = TRUE, status = "primary",
            dataTableOutput("genes")
          )
        )
      )
    )
  )
)
