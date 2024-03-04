server = function(input, output, session) {
  
  selected_df <- reactive ({
    req(input$human, input$cell)
    data <- lineplot[[input$human]][[input$cell]]
    if (input$gene_set != "All") {
      selected_genes <- Genes[[input$gene_set]]$GENES
      data <- data[which(rownames(data) %in% selected_genes), ]
    } else {
      data <- data
    }
  })
  
  morpheus_result <- reactive({
    req(input$human, input$cell)
    morpheus_result <- Morpheus(selected_df())
    return(morpheus_result)
  })
  
  
  pathways_result <- reactive({
    req(input$enrichr, input$enrichr2)
    pathways_result <- pathways(vias[[input$enrichr]], dbs[[input$enrichr2]])
  })
  
  genes_used <- reactive({
    req(input$enrichr)
    genes_used <- dataframe(vias[[input$enrichr]])
  })
  
  output$heatmap <- renderPlot({
    morpheus_data <- morpheus_result()
    title <- paste("Human primary", input$human, "cancer vs", str_to_title(input$cell), "cell lines")
    nomes <- HeatmapAnnotation(df = as.data.frame(morpheus_data["Subtype"]), which = "row", annotation_width = unit(c(1, 3), 'cm'),
                               gap = unit(1, 'mm'), annotation_name_side = "top", annotation_name_rot = 45)
    morpheus_data <- as.matrix(morpheus_data[,2:length(morpheus_data)])
    Heatmap(morpheus_data, right_annotation = nomes, rect_gp = gpar(col = "white", lwd = 0.8),
            column_title = title, column_title_gp = gpar(fontsize = 14, fontface = "bold"),
            column_names_side = "top", column_names_gp = gpar(fontsize = 12), column_names_rot = 80,
            row_names_gp = gpar(fontsize = 10),cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.3f", morpheus_data[i, j]), x, y, gp = gpar(fontsize = 10))},
            heatmap_legend_param = list(title = "R value", at = seq(-1, 1, 0.5)), width = 500)
  })
  
  tabela_result_filtered <- reactiveVal(NULL)  
  
  tabela_result <- reactive({
    req(input$human, input$cell)
    tabela_result <- cor_p_value(selected_df())
  })
  
  observeEvent(input$filter_button, { 
    tabela_filtered <- tabela_result() %>% 
      filter(str_detect(Human_Cancer, paste0("^", input$filter, collapse = "|")))
    tabela_result_filtered(tabela_filtered)  
  })
  
  output$tabela <- DT::renderDT({
    if (!is.null(tabela_result_filtered())) {
      DT::datatable(tabela_result_filtered(), rownames = F, options = list(pageLength = 15))
    } else {
      DT::datatable(tabela_result(), rownames = F, options = list(pageLength = 15))
    }
  })
  output$download_table <- downloadHandler(
    filename = function() {
      paste(input$human, "_cancer_markers", ".csv", sep = "")
    },
    content = function(file) {
      if (!is.null(tabela_result_filtered())) {
        write.csv(tabela_result_filtered(), file, row.names = FALSE)
      } else {
        write.csv(tabela_result, file, row.names = FALSE)
      }
    })
  
  output$enrichr_plot <- renderPlot({
    grafico <- pathways_result()
    ggplot(grafico[1:10,], aes(x=`Combined.Score`, y=`Term`, label = round(`Combined.Score`))) +
      geom_segment(aes(x = 0, y = `Term`, xend = `Combined.Score`, yend = `Term`), linetype = "longdash") +
      geom_point(size = 8, color = "#f2545b") + theme(text = element_text(size = 14)) + 
      geom_text(color = 'black', size = 3, position = "identity", fontface = "bold.italic") +
      theme_bw() +
      xlab("") + ylab("")
  }, res = 100)
  
  output$genes <- DT::renderDT({
    DT::datatable(genes_used(), rownames = F)
  })
  
  output$tabela_enrichr <- DT::renderDT({
    DT::datatable(pathways_result()[,c(1,4,8)], rownames = F, options = list(pageLength = 5))
  })
}