library(DT)
library(shinyjs)
library(dplyr)
library(tibble)
library(data.table)
library(magrittr)


server <- function(input, output, session) {
  values <- reactiveValues(
    pathways = defaultPathways(),
    pathwayMetadataFields = NULL,
    signatures = defaultSignatures(),
    signatureMetadataFields = NULL,
    showHumanReference = TRUE
  )

  output$filterPathwaysOutput <- renderUI({
    selectInput("filterPathways",
      "Pathways",
      getPathways(ist.result),
      selected = values$pathways,
      multiple = TRUE,
      width = "100%"
    )
  })

  output$filterSignaturesOutput <- renderUI({
    selectInput("filterSignatures",
      "Signatures",
      getSignatures(ist.result),
      selected = values$signatures,
      multiple = TRUE,
      width = "100%"
    )
  })

  observeEvent(input$filterUpdateButton, {
    if (!is.null(input$filterPathways)) {
      values$pathways <- input$filterPathways
    }
    if (!is.null(input$filterPathwayMetadataFields)) {
      values$pathwayMetadataFields <- input$filterPathwayMetadataFields
    }
    if (!is.null(input$filterSignatures)) {
      values$signatures <- input$filterSignatures
    }
    if (!is.null(input$filterSignatureMetadataFields)) {
      values$signatureMetadataFields <- input$filterSignatureMetadataFields
    }
    values$showHumanReference <- input$filterShowHumanReference
  })

  output$boxplots <- renderPlot({
    if (values$showHumanReference) {
      signatures <- c(values$signatures, getGroupLevels(ist.result))
    } else {
      signatures <- values$signatures
    }
    plot.ist.boxplots(
      ist.result,
      y = values$pathways,
      sig.ids = signatures,
      mapping = aes(fill = study.id),
      facet_rows = vars(sig.org)
    )
  })

  observe({
    metaheatmap <- plot.ist.pathwaymaps(
      ist.result,
      y = values$pathways,
      sig.ids = values$signatures,
      type = "pheatmap",
      vars.meta.sig = values$signatureMetadataFields,
      vars.meta.path = values$pathwayMetadataFields,
      args.pheatmap = append(default.args.pheatmap("pathwaymap"), list(silent = TRUE))
    )
    output$metaheatmap <- renderPlot(metaheatmap$plot.obj)

    genemaps <- lapply(values$pathways, function(pathway) {
      tryCatch(
        {
          plot <- plot.ist.genemaps(
            ist.result,
            y = pathway,
            type = "pheatmap",
            sig.ids = values$signatures,
            max.genes = 100,
            vars.meta.sig = values$signatureMetadataFields,
            args.pheatmap = list(cluster_rows = TRUE, cluster_cols = TRUE, fontsize = 8, silent = TRUE)
          )
          if (!is.null(plot$plot.obj)) {
            list(pathway, plot$plot.obj)
          } else {
            NULL
          }
        },
        error = function(err) {
          NULL
        }
      )
    })
    output$genemaps <- renderUI({
      lapply(genemaps, function(plot) {
        if (!is.null(plot)) {
          box(
            title = plot[1],
            width = 12,
            renderPlot(plot[2], height = 400),
            tags$p(tags$small("NOTE: For readability this plot shows only up to 100 genes."))
          )
        }
      })
    })
  })

  observeEvent(input$compareButton, {
    pathwayDelta <- getTabDelta(ist.result, "pathway")
    data <- pathwayDelta[sig.id %in% c(input$compareSignature1, input$compareSignature2)] %>%
      data.table::dcast(pathway ~ sig.id, value.var = "total.delta.percent", fill = 0) %>%
      dplyr::mutate(
        path.type = gsub("(\\[.*\\])(.*)", "\\1", pathway),
        diff = .[[3]] - .[[2]],
        largediff = abs(diff) >= input$compareTreshold,
        n.genes = as.numeric(gsub("(.+\\()(\\d+)(\\))", "\\2", pathway))
      )
    data <- subset(data, largediff)
    session$userData$compareList <- data
    output$compareList <- DT::renderDataTable(
      server = FALSE,
      DT::datatable(
        as.data.frame(data),
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "csv", filename = "compare"),
            list(extend = "excel", filename = "compare")
          ),
          pageLength = 10
        )
      ) %>% formatRound(c(2, 3, 5), 1)
    )
    output$comparePlot <- renderPlot({
      ggplot(
        data,
        aes_string(
          x = isolate(input$compareSignature1),
          y = isolate(input$compareSignature2),
          label = "pathway"
        )
      ) +
        geom_hline(yintercept = 0, colour = "gray50") +
        geom_vline(xintercept = 0, colour = "gray50") +
        geom_abline(slope = 1, intercept = 0) +
        geom_point() +
        ggrepel::geom_text_repel(aes(label = pathway), cex = 2.75) +
        scale_colour_brewer(palette = "Set1") +
        theme_bw()
    })
  })

  observeEvent(input$compareUpdateButton, {
    if (!is.null(input$compareList_rows_selected)) {
      values$pathways <- session$userData$compareList$pathway[input$compareList_rows_selected]
    }
    values$signatures <- c(input$compareSignature1, input$compareSignature2)
  })

  output$prioritisePathwaysOutput <- renderUI({
    selectInput("prioritisePathways",
      "Pathways",
      getPathways(ist.result),
      selected = values$pathways,
      multiple = TRUE,
      width = "100%"
    )
  })

  output$prioritiseSignaturesOutput <- renderUI({
    selectInput("prioritiseSignatures",
      "Signatures",
      getSignatures(ist.result),
      selected = values$signatures,
      multiple = TRUE,
      width = "100%"
    )
  })

  observeEvent(input$prioritiseButton, {
    if (input$prioritiseWhat == "signatures") {
      if (is.null(input$prioritisePathways)) {
        return()
      }
      rankedList <- rankPathwaymaps(ist.result,
        what = "signatures",
        max.out = 50,
        return.table = TRUE,
        id.path = input$prioritisePathways
      )
      rankedList <- subset(rankedList,
        select = c(
          "sig.id",
          "perc_mean",
          "perc_mae100",
          "perc_mae100trim",
          "z_mean",
          "z_mae100",
          "z_mae100trim",
          "rank_aggr",
          "z_aggr"
        )
      )
    } else {
      if (is.null(input$prioritiseSignatures)) {
        return()
      }
      rankedList <- rankPathwaymaps(ist.result,
        what = "pathways",
        max.out = 50,
        return.table = TRUE,
        sig.ids = input$prioritiseSignatures
      )
      rankedList <- subset(rankedList,
        select = c(
          "pathway",
          "perc_mean",
          "perc_mae100",
          "perc_mae100trim",
          "z_mean",
          "z_mae100",
          "z_mae100trim",
          "rank_aggr",
          "z_aggr"
        )
      )
    }
    session$userData$prioritisedList <- rankedList
    output$prioritisedList <- DT::renderDataTable(
      server = FALSE,
      DT::datatable(
        rankedList,
        extensions = "Buttons",
        options = list(
          dom = "Bfrtip",
          buttons = list(
            list(extend = "csv", filename = "prioritise"),
            list(extend = "excel", filename = "prioritise")
          ),
          pageLength = 10
        )
      ) %>% formatRound(2:9, 2)
    )
  })

  updatePrioritiseWhat <- function() {
    if (input$prioritiseWhat == "signatures") {
      enable("prioritisePathways")
      disable("prioritiseSignatures")
    } else {
      disable("prioritisePathways")
      enable("prioritiseSignatures")
    }
  }

  observeEvent(input$prioritiseWhat, {
    updatePrioritiseWhat()
  })

  observeEvent(input$prioritiseUpdateButton, {
    if (input$prioritiseWhat == "signatures") {
      if (!is.null(input$prioritiseSignatures)) {
        values$pathways <- input$prioritisePathways
        values$signatures <- unique(as.character(session$userData$prioritisedList[input$prioritisedList_rows_selected]$sig.id))
      }
    } else {
      if (!is.null(input$prioritisePathways)) {
        values$pathways <- unique(as.character(session$userData$prioritisedList[input$prioritisedList_rows_selected]$pathway))
        values$signatures <- input$prioritiseSignatures
      }
    }
  })

  output$filterPathwayMetadata <- DT::renderDataTable(
    DT::datatable(
      getMetaPathways(ist.result),
      selection = list(
        target = "row",
        selected = which(getMetaPathways(ist.result)$path.id %in% isolate(input$filterPathways))
      ),
      options = list(
        pageLength = 10
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    ),
    server = FALSE
  )
  filterPathwayMetadataProxy <- DT::dataTableProxy("filterPathwayMetadata")

  observeEvent(input$filterPathwaysBox, {
    if (input$filterPathwaysBox) {
      disable("filterPathways")
    } else {
      enable("filterPathways")
    }
  })

  observeEvent(input$filterPathways, {
    if (!is.null(input$filterPathwaysBox)) {
      if (!input$filterPathwaysBox) {
        selectRows(
          filterPathwayMetadataProxy,
          which(getMetaPathways(ist.result)$path.id %in% input$filterPathways)
        )
      }
    }
  })

  observeEvent(input$filterPathwayMetadata_rows_selected, {
    if (input$filterPathwaysBox) {
      updateSelectInput(session, "filterPathways",
        choices = getPathways(ist.result),
        selected = getMetaPathways(ist.result)[input$filterPathwayMetadata_rows_selected]$path.id
      )
    }
  })

  output$filterSignatureMetadata <- DT::renderDataTable(
    DT::datatable(
      subset(getMetaSig(ist.result), select = defaultSignatureMetadata()),
      selection = list(
        target = "row",
        selected = which(getMetaSig(ist.result)$sig.id %in% isolate(input$filterSignatures))
      ),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    ),
    server = FALSE
  )
  filterSignatureMetadataProxy <- DT::dataTableProxy("filterSignatureMetadata")

  observeEvent(input$filterSignaturesBox, {
    if (input$filterSignaturesBox) {
      disable("filterSignatures")
    } else {
      enable("filterSignatures")
    }
  })

  observeEvent(input$filterSignatures, {
    if (!is.null(input$filterSignaturesBox)) {
      if (!input$filterSignaturesBox) {
        selectRows(
          filterSignatureMetadataProxy,
          which(getMetaSig(ist.result)$sig.id %in% input$filterSignatures)
        )
      }
    }
  })

  observeEvent(input$filterSignatureMetadata_rows_selected, {
    if (input$filterSignaturesBox) {
      updateSelectInput(session, "filterSignatures",
        choices = getSignatures(ist.result),
        selected = getMetaSig(ist.result)[input$filterSignatureMetadata_rows_selected]$sig.id
      )
    }
  })

  observeEvent(input$tabBox, {
    if (input$tabBox == "Prioritise") {
      updatePrioritiseWhat()
    }
  })

  observeEvent(input$filterResetButton, {
    js$collapse("filterPathwaysBox")
    js$collapse("filterSignaturesBox")
    values$pathways <- defaultPathways()
    values$pathwayMetadataFields <- NULL
    values$signatures <- defaultSignatures()
    values$signatureMetadataFields <- NULL
  })

  output$dataGenes <- DT::renderDataTable({
    DT::datatable(
      genesByPathway[, c("label", "path.id")],
      options = list(
        pageLength = 10,
        dom = "ltp"
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    )
  })

  output$dataGenesVsSignatures <- DT::renderDataTable({
    DT::datatable(
      genesVsSignatures,
      options = list(
        pageLength = 10,
        dom = "ltp"
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    ) %>% formatRound(3)
  })

  output$dataGenesVsSignaturesAndPathways <- DT::renderDataTable({
    DT::datatable(
      genesVsSignaturesAndPathways,
      options = list(
        pageLength = 10,
        dom = "ltp"
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    ) %>% formatRound(4)
  })

  output$dataExpressionGenes <- DT::renderDataTable(
    DT::datatable(
      genesByPathway[path.id == input$dataExpressionPathway, c("label", "path.id")],
      selection = list(
        target = "row"
      ),
      options = list(
        pageLength = 10,
        dom = "ltp"
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    )
  )

  observe({
    genes <- genesByPathway[path.id == input$dataExpressionPathway][input$dataExpressionGenes_rows_selected]$gene.id
    data <- IST:::getDataMatrix(ist.result)
    genes <- genes[genes %in% colnames(data)]
    data <- data[, genes, drop = FALSE] %>% set_colnames(genes)
    if (dim(data)[2] >= 2) {
      output$dataExpressionOutput <- renderUI(plotOutput("dataExpression", height = 600))
      plot <- data %>% pheatmap(
        color = grDevices::colorRampPalette(c("indianred1", "white", "turquoise2"))(15),
        breaks = sym.breaks(as.matrix(.), n = 15),
        scale = "column",
        clustering_method = "ward.D2",
        annotation_row = data.frame(Group = IST::getGroups(ist.result)) %>% set_rownames(rownames(data)),
        cluster_rows = TRUE,
        show_rownames = FALSE,
        main = "Original (scaled) expression values in human reference",
        silent = TRUE
      )
      output$dataExpression <- renderPlot(plot)
    } else {
      output$dataExpressionOutput <- renderUI(tags$p("Please select more genes."))
    }
  })

  observeEvent(input$dataFCPathway, {
    output$dataFCSignaturesOutput <- renderUI({
      selectInput(
        "dataFCSignatures",
        "Signatures",
        genesVsSignaturesAndPathways[path.id == input$dataFCPathway]$sig.id,
        multiple = TRUE
      )
    })
  })

  output$dataFCGenes <- DT::renderDataTable({
    genes <- genesVsSignaturesAndPathways[path.id == input$dataFCPathway & sig.id %in% input$dataFCSignatures]
    genes <- unique(genes[order(genes$label), "label"])
    DT::datatable(
      genes,
      selection = list(target = "row"),
      options = list(
        pageLength = 10,
        dom = "ltp"
      ),
      filter = list(
        position = "top",
        clear = FALSE,
        plain = TRUE
      )
    )
  })

  observe({
    data <- genesVsSignaturesAndPathways[path.id == input$dataFCPathway & sig.id %in% input$dataFCSignatures]
    genes <- unique(data[order(data$label), "label"])[input$dataFCGenes_rows_selected]
    if (!is.null(dim(genes)) && dim(genes)[1] >= 2) {
      data <- data[label %in% genes$label]
      if (!is.null(dim(data)) && dim(data)[1] >= 2) {
        output$dataFCOutput <- renderUI({
          tags$div(
            plotOutput("dataFCBarplot", height = 800),
            plotOutput("dataFCHeatmap", height = 800)
          )
        })

        output$dataFCBarplot <- renderPlot({
          data %>%
            ggplot(aes(x = label, y = logFC, fill = logFC)) +
            geom_hline(yintercept = 0, lty = 1, lwd = .25, colour = "gray30") +
            geom_bar(stat = "identity", colour = "gray30") +
            scale_fill_gradient2(low = "indianred1", mid = "white", high = "turquoise2") +
            facet_wrap(~sig.id, nrow = 1) +
            coord_flip() +
            ylab("logFC (significant genes only)") +
            xlab("Custom gene selection") +
            theme_bw() +
            theme(panel.grid.major.y = element_blank(), legend.position = "bottom") +
            ggtitle("Original fold changes")
        })

        plot <- data %>%
          data.table::dcast(sig.id ~ label, value.var = "logFC", fill = 0) %>%
          column_to_rownames("sig.id") %>%
          pheatmap(
            color = grDevices::colorRampPalette(c("indianred1", "white", "turquoise2"))(15),
            breaks = sym.breaks(as.matrix(.), n = 15),
            clustering_method = "ward.D2", display_numbers = TRUE,
            cluster_rows = FALSE,
            number_format = "%.1f", fontsize_number = 8,
            main = paste0(input$dataFCPathway, "\nsignificant log2 fold changes (non-significant collapsed to 0)"),
            silent = TRUE
          )
        output$dataFCHeatmap <- renderPlot(plot)
      }
    } else {
      output$dataFCOutput <- renderUI(tags$p("Please select some genes."))
    }
  })
}
