library(shiny)
library(heatmaply)
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)

lineplot <- readRDS("lineplot.rds")

cell_annotation <- read_csv("Model.csv")
cell_annotation <- cell_annotation %>%
  select(OncotreeCode, CCLEName)

Morpheus <- function(y) {
  malignos <- y %>% select(ends_with("Malignant"))
  linhagens2 <- y %>% select(-ends_with("Malignant"))
  
  if (nrow(linhagens2) == nrow(malignos)) {
    correlation_2 <- as.data.frame(cor(linhagens2, malignos, method = "pearson"))
    correlation_2 <- cbind("CCLEName" = rownames(correlation_2), correlation_2)
    correlation_2$CCLEName <- rownames(correlation_2)
    cell_annotation <- cell_annotation[which(cell_annotation$CCLEName %in% rownames(correlation_2)), ]
    correlation_2 <- merge(correlation_2, cell_annotation, by = "CCLEName", all.x = TRUE)
    rownames(correlation_2) <- correlation_2$CCLEName
    correlation_2 <- correlation_2 %>% select(-"CCLEName")
    correlation_2 <- correlation_2[, c("OncotreeCode", setdiff(names(correlation_2), "OncotreeCode"))]
    colnames(correlation_2)[1] <- "Subtype"
    return(correlation_2)
  } else {
    print("linhagens2 e malignos têm comprimentos diferentes. A análise não pode ser realizada.")
    return(NULL)
  }
}

cor_p_value <- function(dataframe){
  resultados <- data.frame(
    Human_Cancer = character(),
    Cell_lines = character(),
    r_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  malignos <- dataframe %>% select(ends_with("Malignant"))
  linhagens2 <- dataframe %>% select(-ends_with("Malignant"))
  malignos <- mutate_all(malignos, function(x) as.numeric(as.character(x)))
  linhagens2 <- mutate_all(linhagens2, function(x) as.numeric(as.character(x)))
  
  for (col_malignos in colnames(malignos)) {
    for (col_linhagens2 in colnames(linhagens2)) {
      cor_test <- cor.test(malignos[[col_malignos]], linhagens2[[col_linhagens2]])
      
      resultados <- rbind(
        resultados,
        data.frame(
          Human_Cancer = col_malignos,
          Cell_lines = col_linhagens2,
          r_value = cor_test$estimate,
          p_value = cor_test$p.value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  resultados <- resultados %>% arrange(desc(r_value))
  return(resultados)
}

ui <- fluidPage(
  titlePanel("OncolinesExpressionHub"),
  sidebarLayout(
    sidebarPanel(
      tags$style(".selectize-input, .selectize-dropdown { font-size: 12px; }"),
      selectInput("selected_list", "Please, select the Human Dataset", choices = c(names(lineplot))),
      uiOutput("select_df"),
      br(),
      checkboxInput("show_violin", "Show Violin Plot", value = FALSE) # Checkbox para mostrar o violin plot
    ),
    mainPanel(
      fluidRow(
        column(width = 8,
               plotlyOutput("heatmap", height = "900px"),
               br(),
               conditionalPanel(
                 condition = "input.show_violin == true",
                 plotlyOutput("violin_plot")  # Condicional para mostrar o violin plot
               )
        ),
        column(width = 4,
               br(),
               downloadButton("download_table", "Download Table"),
               dataTableOutput("table_output")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$select_df <- renderUI({
    req(input$selected_list)
    df_choices <- names(lineplot[[input$selected_list]])
    selectInput("selected_df", "Please, select the Cell Line Dataset", choices = df_choices)
  })
  
  output$heatmap <- renderPlotly({
    req(input$selected_list, input$selected_df)
    
    selected_df <- lineplot[[input$selected_list]][[input$selected_df]]
    
    morpheus_result <- Morpheus(selected_df)
    
    if (!is.null(morpheus_result)) {
      title <- paste("Human primary", input$selected_list, "cancer vs", str_to_title(input$selected_df), "cell lines")
      heatmaply(
        morpheus_result[, 2:length(morpheus_result)], colors = viridis(256, option = "inferno"),  
        limits = c(-1,1), row_side_colors = data.frame("Subtypes" = morpheus_result[1]), row_text_angle = 20, 
        Colv = TRUE, grid_gap = 0.25, showticklabels = c(TRUE, TRUE), main = title
      )
    } else {
      print("Erro ao processar os dados para o heatmap.")
    }
  })
  
  output$violin_plot <- renderPlotly({
    req(input$selected_list, input$selected_df)
    
    if(input$show_violin) {  # Verifica se o checkbox está marcado para exibir o violin plot
      selected_df <- lineplot[[input$selected_list]][[input$selected_df]]
      violin <- Morpheus(selected_df)
      violin <- violin[, -1] 
      violin <- as.data.frame(t(violin)) 
      df_long <- gather(data = violin, key = "LinhasCelulares", value = "Correlacao")
      if (!is.null(df_long)) {
        ggplot(df_long, aes(x = LinhasCelulares, y = Correlacao)) +
          geom_violin() +
          labs(x = "Cell Lines", y = "Correlation", 
               title = "Violin Plot of Correlation by Cell Lines and Tumor Types") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
      } else {
        print("Error")
      }
    }
  })
  
  output$table_output <- renderDataTable({
    req(input$selected_list, input$selected_df)
    
    selected_df <- lineplot[[input$selected_list]][[input$selected_df]]
    
    morpheus_result <- Morpheus(selected_df)
    cor_p_table <- cor_p_value(selected_df) 
    
    if (!is.null(cor_p_table)) {
      return(cor_p_table)
    } else {
      print("Erro ao processar os dados para a tabela.")
      return(NULL)
    }
  })
  output$download_table <- downloadHandler(
    filename = function() {
      paste("correlation_table", ".csv", sep = "")
    },
    content = function(file) {
      selected_df <- lineplot[[input$selected_list]][[input$selected_df]]
      write.csv(cor_p_value(selected_df), file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)
