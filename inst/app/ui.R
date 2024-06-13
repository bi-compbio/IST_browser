library(shinydashboard)
library(DT)
library(shinyWidgets)
library(shinyjs)
library(shinyBS)
library(ggplot2)

metadataBox <- function(..., id = "box", title = NULL) {
  boxClass <- "box box-solid collapsed-box"
  collapseTag <- NULL
  collapseTag <- div(
    class = "box-tools pull-right",
    tags$button(
      class = paste0("btn btn-box-tool"),
      `data-widget` = "collapse",
      shiny::icon("plus")
    )
  )
  headerTag <- NULL
  if (!is.null(title) || !is.null(collapseTag)) {
    headerTag <- div(class = "box-header", tags$label(title), collapseTag)
  }
  div(
    class = "col-sm-100%",
    div(
      id = id,
      class = boxClass,
      headerTag,
      div(class = "box-body", ...)
    ),
    tags$script(
      sprintf(
        "$('#%s').on('hidden.bs.collapse', function () { Shiny.setInputValue('%s', false, {priority: 'event'}); })",
        id, id
      ),
      sprintf(
        "$('#%s').on('shown.bs.collapse', function () { Shiny.setInputValue('%s', true, {priority: 'event'}); })",
        id, id
      ),
      sprintf(
        "$(document).on('shiny:connected', function(event) { Shiny.setInputValue('%s', false, {priority: 'event'}); });",
        id
      )
    )
  )
}

exploreTab <- fluidPage(
  fluidRow(
    div(
      class = "col-sm-12",
      div(
        class = "box box-solid",
        div(
          class = "box-body",
          tags$h1(getShinyOption("ist.browser.title", basename(getShinyOption("ist.result.path", "IST browser")))),
          tags$p(getShinyOption("ist.browser.description", ""))
        )
      )
    ),
    tabBox(
      id = "tabBox",
      tabPanel(
        "Filter",
        uiOutput("filterPathwaysOutput"),
        metadataBox(
          id = "filterPathwaysBox",
          title = "Pathway metadata",
          DT::dataTableOutput("filterPathwayMetadata")
        ),
        selectInput("filterPathwayMetadataFields",
          "Pathway metadata fields",
          names(getMetaPathways(ist.result)),
          multiple = TRUE,
          width = "100%"
        ),
        uiOutput("filterSignaturesOutput"),
        metadataBox(
          id = "filterSignaturesBox",
          title = "Signature metadata",
          DT::dataTableOutput("filterSignatureMetadata")
        ),
        selectInput("filterSignatureMetadataFields",
          "Signature metadata fields",
          names(getMetaSig(ist.result)),
          multiple = TRUE,
          width = "100%"
        ),
        checkboxInput("filterShowHumanReference", "Show human reference in box plot", TRUE),
        bsButton("filterUpdateButton", "Update plots", style = "primary"),
        actionButton("filterResetButton", "Reset to defaults")
      ),
      tabPanel(
        "Compare",
        selectInput("compareSignature1",
          "Signature 1",
          getSignatures(ist.result),
          selected = getSignatures(ist.result)[1]
        ),
        selectInput("compareSignature2",
          "Signature 2",
          getSignatures(ist.result),
          selected = getSignatures(ist.result)[2]
        ),
        numericInput("compareTreshold", "Treshold", 10, min = 1),
        tags$p(tags$small("Only differences above percent threshold will be displayed.")),
        bsButton("compareButton", "Compare", style = "primary"),
        hr(),
        plotOutput("comparePlot", height = 400),
        box(
          title = span(icon("question"), "Captions"),
          width = "100%",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          tags$div(
            tags$dl(
              tags$dt("<contrast name>"),
              tags$dd("Respective contrast recapitulation percentage."),
              tags$dt("path.type"),
              tags$dd("Pathway metadata."),
              tags$dt("diff"),
              tags$dd("Percent difference in recapitulations."),
              tags$dt("largediff"),
              tags$dd("Is the difference large, according to the user-defined threshold? Value: true/false."),
              tags$dt("n.genes"),
              tags$dd("Number of genes in original gene set.")
            )
          )
        ),
        DT::dataTableOutput("compareList"),
        bsButton("compareUpdateButton", "Update plots", style = "primary")
      ),
      tabPanel(
        "Prioritise",
        selectInput("prioritiseWhat",
          "I would like to prioritise",
          c("signatures", "pathways"),
          selected = "signatures"
        ),
        uiOutput("prioritisePathwaysOutput"),
        uiOutput("prioritiseSignaturesOutput"),
        bsButton("prioritiseButton", "Prioritise", style = "primary"),
        hr(),
        box(
          title = span(icon("question"), "Captions"),
          width = "100%",
          status = "info",
          solidHeader = TRUE,
          collapsible = TRUE,
          collapsed = TRUE,
          tags$div(
            tags$dl(
              tags$dt("perc_"),
              tags$dd("Prefix indicates that the column displays percentages. Higher is better."),
              tags$dt("rank_"),
              tags$dd("Prefix indicates that the column displays rankings. Lower is better."),
              tags$dt("z_"),
              tags$dd("Prefix indicates that the column displays z-scores. Higher is better."),
              tags$dt("_mean/_mae/_mae100"),
              tags$dd("Metrics to prioritise the contrasts."),
              tags$dt("perc_mean"),
              tags$dd(
                "Mean recapitulation of gene sets.",
                "Higher is better, if we accept overshoot.",
                "Ideally close to 100%."
              ),
              tags$dt("_mae100/_mae100trim"),
              tags$dd(
                "Mean absolute error from 100%, i.e. how far recapitulations are, in average, from 100%.",
                "Lower is better.",
                "Mae100 penalises overshoot, Mae100trim does not, i.e. a 120% would have an error of 0%."
              ),
              tags$dt("_aggr"),
              tags$dd(
                "Aggregated mean ranking, comes from averaging the rankings from respective column.",
                "Lower is better."
              )
            )
          )
        ),
        DT::dataTableOutput("prioritisedList"),
        bsButton("prioritiseUpdateButton", "Update plots", style = "primary"),
      ),
      width = 12
    )
  ),
  fluidRow(),
  fluidRow(
    box(
      title = "Boxplots",
      width = 12,
      plotOutput("boxplots", height = 800)
    ),
    box(
      title = "Metaheatmap",
      width = 12,
      plotOutput("metaheatmap", width = 1000, height = 800)
    ),
    box(
      title = "Genemaps",
      width = 12,
      uiOutput("genemaps")
    )
  )
)

dataTab <- fluidPage(
  div(
    class = "col-sm-12",
    div(
      class = "box box-solid",
      div(
        class = "box-body",
        tags$h1(getShinyOption("ist.browser.title", basename(getShinyOption("ist.result.path", "IST browser")))),
        tags$p(getShinyOption("ist.browser.description", ""))
      )
    )
  ),
  box(
    title = "Gene sets",
    width = 12,
    DT::dataTableOutput("dataGenes")
  ),
  box(
    title = "Genes vs signatures",
    width = 12,
    DT::dataTableOutput("dataGenesVsSignatures")
  ),
  box(
    title = "Genes vs signatures and pathways",
    width = 12,
    DT::dataTableOutput("dataGenesVsSignaturesAndPathways")
  ),
  box(
    title = "Gene expression",
    width = 12,
    selectInput(
      "dataExpressionPathway",
      "Pathway",
      getPathways(ist.result),
      selected = defaultPathways()[1]
    ),
    DT::dataTableOutput("dataExpressionGenes"),
    hr(),
    uiOutput("dataExpressionOutput")
  ),
  box(
    title = "Fold changes",
    width = 12,
    selectInput(
      "dataFCPathway",
      "Pathway",
      getPathways(ist.result),
      selected = defaultPathways()[1]
    ),
    uiOutput("dataFCSignaturesOutput"),
    DT::dataTableOutput("dataFCGenes"),
    hr(),
    uiOutput("dataFCOutput")
  )
)

body <- dashboardBody(
  useShinyjs(),
  extendShinyjs(
    text = "shinyjs.collapse = function(boxid) {
                    $('#' + boxid + ':not(.collapsed-box)').closest('.box').find('[data-widget=collapse]').click();
                }",
    functions = c("collapse")
  ),
  tabItems(
    tabItem(tabName = "explore", exploreTab),
    tabItem(tabName = "data", dataTab)
  )
)

sidebar <- dashboardSidebar(
  collapsed = TRUE,
  sidebarMenu(
    id = "tabs",
    menuItem("Explore", tabName = "explore", icon = icon("dashboard")),
    menuItem("Data overview", tabName = "data", icon = icon("table"))
  )
)

ui <- dashboardPage(
  dashboardHeader(title = "IST Browser"),
  sidebar,
  body
)
