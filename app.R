library(shiny)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(readr)
library(shinyWidgets)
library(colourpicker)

ui <- fluidPage(
  titlePanel("Volcano Plot Generator"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload CSV File", accept = ".csv"),
      
      uiOutput("columnSelectors"),
      uiOutput("baitSelector"),
      
      tags$hr(),
      numericInput("fc_thresh", "Fold Change Threshold", value = 1.5, step = 0.1),
      numericInput("pval_thresh", "Adjusted P-value Threshold", value = 0.05, step = 0.01),
      
      numericInput("top_up", "Number of Upregulated Genes/Sites to Label", value = 15, min = 0),
      numericInput("top_down", "Number of Downregulated Genes/Sites to Label", value = 15, min = 0),
      
      tags$hr(),
      checkboxInput("custom_colors", "Use Custom Colors", FALSE),
      
      conditionalPanel(
        condition = "input.custom_colors == true",
        colourInput("up_color", "Upregulated Color", value = "#D26E01"),
        colourInput("down_color", "Downregulated Color", value = "darkorchid4"),
        colourInput("nonsig_color", "Non-significant Color", value = "grey70")
      ),
      
      tags$hr(),
      downloadButton("downloadPlot", "Download Plot as PDF")
    ),
    
    mainPanel(
      plotOutput("volcanoPlot", height = "800px")
    )
  )
)

server <- function(input, output, session) {
  
  # Load data
  dataset <- reactive({
    req(input$file)
    read_csv(input$file$datapath)
  })
  
  # Dynamic selectors for columns
  output$columnSelectors <- renderUI({
    req(dataset())
    cols <- names(dataset())
    
    tagList(
      selectInput("logfc_col", "Select log2 Fold Change Column", choices = cols),
      selectInput("pval_col", "Select Adjusted P-value Column", choices = cols),
      selectInput("gene_col", "Select Gene or Site Name Column", choices = cols)
    )
  })
  
  # Processed data for plotting
  plotData <- reactive({
    req(dataset(), input$logfc_col, input$pval_col, input$gene_col)
    
    df <- dataset()
    
    df <- df %>%
      mutate(
        logFC = .[[input$logfc_col]],
        adjP = .[[input$pval_col]],
        Gene = .[[input$gene_col]],
        neglog10p = -log10(adjP),
        group = case_when(
          logFC >= input$fc_thresh & adjP <= input$pval_thresh ~ "Upregulated",
          logFC <= -input$fc_thresh & adjP <= input$pval_thresh ~ "Downregulated",
          TRUE ~ "Non-significant"
        )
      )
    
    df
  })
  
  # Top labels
  topLabels <- reactive({
    df <- plotData()
    up <- df %>% filter(group == "Upregulated") %>% arrange(adjP) %>% slice_head(n = input$top_up)
    down <- df %>% filter(group == "Downregulated") %>% arrange(adjP) %>% slice_head(n = input$top_down)
    bind_rows(up, down)
  })
  
  # Bait gene selector
  output$baitSelector <- renderUI({
    req(plotData())
    genes <- unique(plotData()$Gene)
    selectizeInput("bait_genes", "Select Highlight Gene(s)/Site(s)", 
                   choices = genes, multiple = TRUE)
  })
  
  # Subset of bait genes
  baitData <- reactive({
    req(plotData())
    if (is.null(input$bait_genes) || length(input$bait_genes) == 0) {
      return(NULL)
    }
    plotData() %>% filter(Gene %in% input$bait_genes)
  })
  
  
  # Volcano plot
  volcanoPlot <- reactive({
    df <- plotData()
    top_labels <- topLabels()
    bait_labels <- baitData()
    
    color_values <- if (input$custom_colors) {
      c("Upregulated" = input$up_color,
        "Downregulated" = input$down_color,
        "Non-significant" = input$nonsig_color)
    } else {
      c("Upregulated" = "#D26E01",
        "Downregulated" = "darkorchid4",
        "Non-significant" = "grey70")
    }
    
    p <- ggplot(df, aes(x = logFC, y = neglog10p, color = group)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = color_values) +
      geom_hline(yintercept = -log10(input$pval_thresh), linetype = "dashed", color = "grey30") +
      geom_vline(xintercept = c(-input$fc_thresh, input$fc_thresh), linetype = "dashed", color = "grey30") +
      geom_label_repel(
        data = top_labels, aes(label = Gene),
        show.legend = FALSE, size = 4,
        box.padding = 0.5, point.padding = 0.5,
        force = 6, max.overlaps = 50
      ) +
      labs(
        x = bquote(Log[2] ~ "Fold Change"),
        y = bquote(-Log[10] ~ "Adjusted P-value"),
        color = NULL
      ) +
      theme_minimal(base_size = 15) +
      theme(
        legend.position = "top",
        axis.line = element_line(size = 0.8, color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()
      )
    
    # Only add bait highlighting if some genes are selected
    if (!is.null(bait_labels) && nrow(bait_labels) > 0) {
      p <- p +
        geom_point(data = bait_labels, shape = 21, fill = "firebrick", color = "black", size = 4, stroke = 1, show.legend = FALSE) +
        geom_label_repel(
          data = bait_labels, aes(label = Gene),
          color = "firebrick", fill = "white",
          size = 4, box.padding = 0.5, point.padding = 0.5,
          force = 6, max.overlaps = 100, show.legend = FALSE
        )
    }
    
    p
  })
  
  
  output$volcanoPlot <- renderPlot({
    req(volcanoPlot())
    volcanoPlot()
  })
  
  # PDF download
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = volcanoPlot(), width = 10, height = 8, dpi = 300)
    }
  )
}

shinyApp(ui, server)
