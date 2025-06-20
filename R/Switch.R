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

required_packages <- c("shiny", "readxl", "dplyr", "tidyr")
check_and_install_packages(required_packages)

ui <- fluidPage(
  titlePanel("Conversor de Excel para ARP"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Carregar arquivo Excel (.xlsx)", accept = ".xlsx"),
      uiOutput("column_selector"),
      actionButton("convert", "Converter"),
      uiOutput("download_ui")
    ),
    mainPanel(
      verbatimTextOutput("messages")
    )
  )
)

server <- function(input, output, session) {
  # Armazena os dados carregados
  data_loaded <- reactiveVal(NULL)

  # Atualiza quando um arquivo é carregado
  observeEvent(input$file, {
    req(input$file)

    tryCatch({
      # Ler o arquivo Excel
      data <- read_excel(input$file$datapath)
      data_loaded(data)

      # Atualizar seletor de colunas
      output$column_selector <- renderUI({
        selectInput("group_column", "Selecione a coluna para agrupamento",
                    choices = colnames(data))
      })

      output$messages <- renderText("Arquivo carregado com sucesso!")
    }, error = function(e) {
      output$messages <- renderText(paste("Erro ao carregar arquivo:", e$message))
    })
  })

  # Gera o arquivo ARP quando o botão é clicado
  observeEvent(input$convert, {
    req(input$file, input$group_column, data_loaded())

    tryCatch({
      data <- data_loaded()

      # Verificar se as colunas necessárias existem
      required_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2")
      if(!all(required_cols %in% colnames(data))) {
        stop("O arquivo Excel não contém todas as colunas necessárias (A1, A2, B1, B2, DRB1_1, DRB1_2)")
      }

      # Renomear a coluna de grupo para "disease" para compatibilidade com a função original
      data <- data %>% rename(disease = !!sym(input$group_column))

      # Processar metadados e estrutura
      profile <- c(
        "[Profile]",
        '  Title = "Dados_convertidos"',
        paste("  NbSamples =", length(unique(data$disease))),
        "  DataType = MICROSAT",
        "  GenotypicData = 1",
        "  LocusSeparator = WHITESPACE",
        '  MissingData = "?"',
        "  GameticPhase = 0",
        "  RecessiveData = 0"
      )

      # Processar dados genéticos
      diseases <- unique(data$disease)
      data_sections <- list()

      for(disease in diseases) {
        sample_data <- data %>% filter(disease == disease)

        # Gerar IDs aleatórios no formato C + número sequencial
        sample_data$IndividualID <- paste0("C", 1:nrow(sample_data))

        # Formatar os dados genéticos
        formatted_data <- sample_data %>%
          mutate(
            allele1 = paste(
              paste0("A*", A1),
              paste0("B*", B1),
              paste0("DRB1*", DRB1_1),
              sep = " "),
            allele2 = paste(
              paste0("A*", A2),
              paste0("B*", B2),
              paste0("DRB1*", DRB1_2),
              sep = " ")
          ) %>%
          select(IndividualID, allele1, allele2)

        # Criar strings no formato correto
        genetic_lines <- character()
        for(i in 1:nrow(formatted_data)) {
          genetic_lines <- c(
            genetic_lines,
            sprintf("%s\t1\t%s", formatted_data$IndividualID[i], formatted_data$allele1[i]),
            sprintf("\t\t%s", formatted_data$allele2[i])
          )
        }

        # Criar seção da amostra
        data_sections[[disease]] <- c(
          "[[Samples]]",
          sprintf('  SampleName = "%s"', disease),
          sprintf("  SampleSize = %d", nrow(sample_data)),
          "  SampleData = {",
          genetic_lines,
          "  }"
        )
      }

      # Criar conteúdo do arquivo
      arp_content <- c(profile, "", "[Data]", "", unlist(data_sections))

      # Armazenar o conteúdo para download
      output$download_ui <- renderUI({
        downloadButton("download_arp", "Baixar arquivo .arp")
      })

      output$download_arp <- downloadHandler(
        filename = function() {
          paste0(tools::file_path_sans_ext(input$file$name), ".arp")
        },
        content = function(file) {
          writeLines(arp_content, file)
        }
      )

      output$messages <- renderText("Conversão concluída! Clique em 'Baixar arquivo .arp'")

    }, error = function(e) {
      output$messages <- renderText(paste("Erro na conversão:", e$message))
    })
  })
}

shinyApp(ui, server)
}
