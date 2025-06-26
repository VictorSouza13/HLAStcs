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
      verbatimTextOutput("messages"),
      tableOutput("conversion_info")
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

  # Função para converter numérico em categórico usando quantis
  auto_convert_to_category <- function(x) {
    if(is.numeric(x)) {
      # Usar quartis como padrão, mas ajusta o número de grupos conforme o tamanho dos dados
      n_groups <- min(4, length(unique(na.omit(x))))
      if(n_groups > 1) {
        breaks <- quantile(x, probs = seq(0, 1, length.out = n_groups + 1), na.rm = TRUE)
        labels <- paste0("Q", 1:n_groups)
        return(cut(x, breaks = unique(breaks), labels = labels, include.lowest = TRUE))
      }
    }
    return(as.character(x)) # Mantém como está se não for numérico ou não puder ser agrupado
  }

  # Função para limpar os prefixos dos alelos
  clean_allele <- function(allele) {
    str_replace(allele, ".*\\*", "")
  }

  # Gera o arquivo ARP quando o botão é clicado
  observeEvent(input$convert, {
    req(input$file, input$group_column, data_loaded())

    tryCatch({
      data <- data_loaded()

      # Verificar colunas necessárias
      required_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2", input$group_column)
      if(!all(required_cols %in% colnames(data))) {
        stop("O arquivo Excel não contém todas as colunas necessárias")
      }

      # Converter coluna de agrupamento se for numérica
      original_values <- data[[input$group_column]]
      converted_values <- auto_convert_to_category(original_values)

      # Mostrar informações da conversão (APENAS SE FOR NUMÉRICO)
      if(is.numeric(original_values)) {
        output$conversion_info <- renderTable({
          # Criar um resumo da conversão em vez de mapear todos os valores
          conversion_summary <- data.frame(
            Categoria = levels(converted_values),
            Minimo = tapply(original_values, converted_values, min, na.rm = TRUE),
            Maximo = tapply(original_values, converted_values, max, na.rm = TRUE),
            Observacoes = as.numeric(table(converted_values))
          )
          names(conversion_summary) <- c("Categoria", "Valor Mínimo", "Valor Máximo", "Nº Observações")
          conversion_summary
        }, caption = "Resumo da Conversão Automática por Quartis", caption.placement = "top")
      } else {
        output$conversion_info <- NULL
      }

      data[[input$group_column]] <- converted_values

      # Criar coluna temporária para o grupo
      data <- data %>% mutate(group_column = !!sym(input$group_column))

      # Processar metadados e estrutura
      profile <- c(
        "[Profile]",
        '  Title = "Dados_convertidos"',
        paste("  NbSamples =", length(unique(data$group_column))),
        "  DataType = MICROSAT",
        "  GenotypicData = 1",
        "  LocusSeparator = WHITESPACE",
        '  MissingData = "?"',
        "  GameticPhase = 0",
        "  RecessiveData = 0"
      )

      # Processar dados genéticos
      groups <- unique(data$group_column)
      data_sections <- list()

      for(group in groups) {
        sample_data <- data %>% filter(group_column == group)

        # Gerar IDs
        sample_data$IndividualID <- paste0("C", 1:nrow(sample_data))

        # Limpar prefixos dos alelos
        sample_data <- sample_data %>%
          mutate(across(c(A1, A2, B1, B2, DRB1_1, DRB1_2), clean_allele))

        # Formatar dados genéticos
        formatted_data <- sample_data %>%
          mutate(
            allele1 = paste(paste0("A*", A1), paste0("B*", B1), paste0("DRB1*", DRB1_1), sep = " "),
            allele2 = paste(paste0("A*", A2), paste0("B*", B2), paste0("DRB1*", DRB1_2), sep = " ")
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
        data_sections[[as.character(group)]] <- c(
          "[[Samples]]",
          sprintf('  SampleName = "%s"', group),
          sprintf("  SampleSize = %d", nrow(sample_data)),
          "  SampleData = {",
          genetic_lines,
          "  }"
        )
      }

      # Criar conteúdo do arquivo
      arp_content <- c(profile, "", "[Data]", "", unlist(data_sections))

      # Preparar download
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

      output$messages <- renderText({
        if(is.numeric(original_values)) {
          paste("Conversão concluída! Valores numéricos foram automaticamente agrupados por quartis. Clique em 'Baixar arquivo .arp'")
        } else {
          "Conversão concluída! Clique em 'Baixar arquivo .arp'"
        }
      })

    }, error = function(e) {
      output$messages <- renderText(paste("Erro na conversão:", e$message))
    })
  })
}

shinyApp(ui, server)
}
