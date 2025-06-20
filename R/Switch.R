#' Função de converter arquivos ARP
#'
#' Esta função executa a conversão de planilhas para ARP sem complicações.
#' @export
Convert <- function() {
check_and_install_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

required_packages <- c("shiny", "openxlsx")
check_and_install_packages(required_packages)

ui <- fluidPage(
  titlePanel("Conversor de Excel para ARP (Arlequin)"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Selecione o arquivo Excel", accept = c(".xlsx", ".xls")),
      uiOutput("column_selector"),
      textInput("output_name", "Nome do arquivo de saída:", value = "dados_output.arp"),
      actionButton("convert", "Converter Arquivo"),
      tags$hr(),
      helpText("Selecione a coluna para filtrar (será renomeada para 'pops' no arquivo ARP).")
    ),
    mainPanel(
      verbatimTextOutput("conversion_status"),
      tags$hr(),
      uiOutput("download_ui")
    )
  )
)

server <- function(input, output) {

  output$column_selector <- renderUI({
    req(input$file)
    df <- read.xlsx(input$file$datapath)
    locus_cols <- grep("A[0-9]|B[0-9]|DRB1_[0-9]", names(df), value = TRUE)
    available_cols <- setdiff(names(df), locus_cols)
    selectInput("filter_column", "Coluna para filtrar (será 'pops'):",
                choices = available_cols)
  })

  observeEvent(input$convert, {
    req(input$file, input$filter_column)

    tryCatch({
      df <- read.xlsx(input$file$datapath)

      names(df)[names(df) == input$filter_column] <- "pops"

      allele_cols <- grep("A[0-9]|B[0-9]|DRB1_[0-9]", names(df), value = TRUE)
      populations <- unique(df$pops)
      npops <- length(populations)
      nloci <- length(allele_cols)/2

      header <- c(
        paste0("[Profile]"),
        paste0("  Title=\"Arquivo convertido de ", input$file$name, "\""),
        paste0("  NbSamples=", npops),
        paste0("  GenotypicData=1"),
        paste0("  LocusSeparator=WHITESPACE"),
        paste0("  DataType=STANDARD"),
        paste0("  MissingData=\"?\""),
        paste0("  GameticPhase=0"),
        paste0(""),
        paste0("[Data]"),
        paste0("  [[Samples]]")
      )

      content <- header

      for (pop in populations) {
        pop_data <- df[df$pops == pop, ]
        content <- c(content, paste0("    SampleName=\"", pop, "\""))
        content <- c(content, paste0("    SampleSize=", nrow(pop_data)))
        content <- c(content, "    SampleData={")

        for (i in 1:nrow(pop_data)) {
          ind_id <- ifelse(grepl("case", pop, ignore.case = TRUE), paste0("P", i), paste0("C", i))

          alleles <- sapply(allele_cols, function(col) as.character(pop_data[i, col]))

          locus_groups <- split(alleles, rep(1:nloci, each = 2))
          formatted_alleles <- sapply(locus_groups, function(x) paste(x, collapse = " "))
          ind_line <- paste0(ind_id, "\t", paste(formatted_alleles, collapse = "\t"))

          content <- c(content, paste0("      ", ind_line))
        }

        content <- c(content, "    }")
      }

      output_file <- ifelse(
        grepl("\\.arp$", input$output_name, ignore.case = TRUE),
        input$output_name,
        paste0(input$output_name, ".arp")
      )

      writeLines(content, output_file)

      output$conversion_status <- renderText({
        paste("Conversão concluída!\nArquivo salvo como:", output_file,
              "\nColuna original usada como 'pops':", input$filter_column)
      })

      output$download_ui <- renderUI({
        downloadButton("download", "Download do Arquivo ARP")
      })

      output$download <- downloadHandler(
        filename = function() { output_file },
        content = function(file) { file.copy(output_file, file) }
      )

    }, error = function(e) {
      output$conversion_status <- renderText({
        paste("ERRO NA CONVERSÃO:\n", e$message)
      })
    })
  })
}

shinyApp(ui = ui, server = server)
}
