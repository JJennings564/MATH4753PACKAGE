#' SHINY MLE
#'
#' @returns shiny app
#' @export
#'
#' @examples
#' \dontrun{shinymle()}
shinymle <- function(){
  shiny::runApp(system.file("SHINY", package="JesseLabPkg"),launch.browser = TRUE)
}
