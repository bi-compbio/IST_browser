#' @title Run ISTbrowser app
#'
#' @description Runs \code{shiny::runApp()} on the ISTbrowser shiny app,
#' as stored in the inst folder within this package.
#' This function has no arguments.
#'
#' @examples
#' \dontrun{run()}
#'
#' @importFrom shiny runApp
#' @export
run <- function() {
  shiny::runApp(system.file("app", package = "ISTBrowser"), display.mode = "normal")
}
