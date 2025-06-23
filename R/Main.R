#' Função principal do pacote HLAStcs
#'
#' Esta função executa a análise principal do pacote.
#' @export
HLAStcs <- function() {
  required_packages <- c(
    "openxlsx", "shiny", "DT"
  )

  install_missing_packages <- function(packages) {
    for(pkg in packages) {
      if(!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
      }
    }
  }

  install_missing_packages(required_packages)

  library(openxlsx)
  library(shiny)
  library(DT)

  read_hla_data <- function(file_path) {
    tryCatch({
      if(!file.exists(file_path)) {
        stop("Arquivo não encontrado. Verifique o caminho.")
      }

      hla_data <- read.xlsx(file_path)

      if(is.null(hla_data) || nrow(hla_data) == 0) {
        stop("Arquivo vazio ou formato inválido.")
      }

      drb1_cols <- grep("^DRB1", colnames(hla_data), value = TRUE, ignore.case = TRUE)

      if(length(drb1_cols) < 2) {
        stop("São necessárias pelo menos duas colunas DRB1 (DRB1_1 e DRB1_2)")
      }

      if(!all(c("DRB1_1", "DRB1_2") %in% colnames(hla_data))) {
        colnames(hla_data)[match(drb1_cols[1:2], colnames(hla_data))] <- c("DRB1_1", "DRB1_2")
      }

      required_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2")
      missing_cols <- setdiff(required_cols, colnames(hla_data))

      if(length(missing_cols) > 0) {
        stop(paste("Colunas faltando:", paste(missing_cols, collapse = ", ")))
      }

      loci <- c("A", "B", "DRB1")
      hla_alleles <- list()

      for(locus in loci) {
        col1 <- ifelse(locus == "DRB1", paste0(locus, "_1"), paste0(locus, "1"))
        col2 <- ifelse(locus == "DRB1", paste0(locus, "_2"), paste0(locus, "2"))

        alleles <- c(trimws(as.character(hla_data[[col1]])), trimws(as.character(hla_data[[col2]])))
        alleles <- alleles[!is.na(alleles) & alleles != "" & alleles != "NA"]

        if(locus == "DRB1") {
          alleles <- alleles[grepl("^DRB1\\*\\d+", alleles)]

          if(length(alleles) == 0) {
            stop("Nenhum alelo DRB1 válido encontrado. Verifique se os dados seguem o formato 'DRB1*NN'")
          }
        }

        hla_alleles[[locus]] <- alleles
      }

      return(list(data = hla_data, alleles = hla_alleles))

    }, error = function(e) {
      stop(paste("Erro na leitura:", e$message))
    })
  }

  calculate_allele_frequencies <- function(allele_list) {
    if(is.null(allele_list)) stop("Dados de alelos ausentes")

    freq_list <- list()

    for(locus in names(allele_list)) {
      alleles <- allele_list[[locus]]
      if(length(alleles) == 0) next

      freq_table <- as.data.frame(table(alleles))
      if(nrow(freq_table) == 0) next

      freq_table$Frequency <- freq_table$Freq / sum(freq_table$Freq)
      freq_table$Percentage <- round(freq_table$Frequency * 100, 2)
      freq_table$Count <- freq_table$Freq
      freq_table <- freq_table[order(-freq_table$Frequency), c("alleles", "Count", "Percentage")]

      freq_list[[locus]] <- freq_table
    }

    if(length(freq_list) == 0) stop("Não foi possível calcular frequências")
    return(freq_list)
  }

  test_hwe <- function(hla_data, loci = c("A", "B", "DRB1")) {
    hwe_results <- list()

    for(locus in loci) {
      col1 <- paste0(locus, ifelse(locus == "DRB1", "_1", "1"))
      col2 <- paste0(locus, ifelse(locus == "DRB1", "_2", "2"))

      if(!all(c(col1, col2) %in% colnames(hla_data))) {
        warning(paste("Colunas", col1, "ou", col2, "não encontradas"))
        next
      }

      allele1 <- trimws(as.character(hla_data[[col1]]))
      allele2 <- trimws(as.character(hla_data[[col2]]))

      valid <- !is.na(allele1) & !is.na(allele2) &
        allele1 != "" & allele2 != "" &
        !grepl("^NA$", allele1) & !grepl("^NA$", allele2)
      allele1 <- allele1[valid]
      allele2 <- allele2[valid]

      if(length(allele1) < 5) {
        warning(paste("Amostra muito pequena para", locus, "(n =", length(allele1), ")"))
        next
      }

      allele1_simple <- gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", allele1)
      allele2_simple <- gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", allele2)

      geno_data <- paste(allele1_simple, allele2_simple, sep = "/")

      geno_counts <- table(geno_data)
      alleles <- c(allele1_simple, allele2_simple)
      allele_freq <- table(alleles)/(length(alleles))

      expected_counts <- numeric(length(geno_counts))
      names(expected_counts) <- names(geno_counts)

      for(geno in names(geno_counts)) {
        alleles_in_geno <- unlist(strsplit(geno, "/"))
        if(length(alleles_in_geno) != 2) next

        if(alleles_in_geno[1] == alleles_in_geno[2]) {
          expected <- allele_freq[alleles_in_geno[1]]^2
        } else {
          expected <- 2 * allele_freq[alleles_in_geno[1]] * allele_freq[alleles_in_geno[2]]
        }
        expected_counts[geno] <- expected * length(allele1)
      }

      observed <- as.numeric(geno_counts)
      expected <- as.numeric(expected_counts)

      valid_categories <- expected > 0
      observed <- observed[valid_categories]
      expected <- expected[valid_categories]

      if(length(observed) < 2) {
        warning(paste("Não há diversidade suficiente para testar HWE em", locus))
        next
      }

      df <- length(observed) - length(allele_freq)
      if(df < 1) df <- 1

      chisq <- sum((observed - expected)^2 / expected)
      p.value <- pchisq(chisq, df, lower.tail = FALSE)

      hwe_results[[locus]] <- list(
        chisq = round(chisq, 3),
        p.value = format.pval(p.value, digits = 4),
        df = df,
        n_alleles = length(allele_freq),
        n_individuals = length(allele1),
        observed = geno_counts,
        expected = round(expected_counts, 2)
      )
    }

    return(hwe_results)
  }

  calculate_all_ld <- function(hla_data) {
    loci <- c("A", "B", "DRB1")
    ld_results <- list()

    locus_pairs <- combn(loci, 2, simplify = FALSE)

    for(pair in locus_pairs) {
      locus1 <- pair[1]
      locus2 <- pair[2]

      tryCatch({
        col1_1 <- paste0(locus1, ifelse(locus1 == "DRB1", "_1", "1"))
        col1_2 <- paste0(locus1, ifelse(locus1 == "DRB1", "_2", "2"))
        col2_1 <- paste0(locus2, ifelse(locus2 == "DRB1", "_1", "1"))
        col2_2 <- paste0(locus2, ifelse(locus2 == "DRB1", "_2", "2"))

        geno1 <- paste(
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[col1_1]]),
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[col1_2]]),
          sep = "/"
        )

        geno2 <- paste(
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[col2_1]]),
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[col2_2]]),
          sep = "/"
        )

        valid_rows <- !is.na(geno1) & !is.na(geno2) &
          geno1 != "NA/NA" & geno2 != "NA/NA" &
          !grepl("^NA/", geno1) & !grepl("/NA$", geno1) &
          !grepl("^NA/", geno2) & !grepl("/NA$", geno2)

        geno1 <- geno1[valid_rows]
        geno2 <- geno2[valid_rows]

        if(length(geno1) < 10) next

        alleles1 <- unique(unlist(strsplit(geno1, "/")))
        alleles2 <- unique(unlist(strsplit(geno2, "/")))

        if(length(alleles1) < 2 || length(alleles2) < 2) next

        n <- length(geno1)
        cont_table <- matrix(0, nrow = length(alleles1), ncol = length(alleles2))
        rownames(cont_table) <- alleles1
        colnames(cont_table) <- alleles2

        for(i in 1:n) {
          a1 <- unlist(strsplit(geno1[i], "/"))
          a2 <- unlist(strsplit(geno2[i], "/"))

          for(x in a1) {
            for(y in a2) {
              cont_table[x, y] <- cont_table[x, y] + 0.5
            }
          }
        }

        p_A <- rowSums(cont_table)/n
        p_B <- colSums(cont_table)/n

        D <- cont_table/n - outer(p_A, p_B)

        Dmax <- pmin(outer(p_A, 1-p_B), outer(1-p_A, p_B))
        Dmin <- pmax(outer(-p_A, p_B), outer(p_A, -p_B))

        D_prime <- matrix(0, nrow = nrow(D), ncol = ncol(D))
        for(i in 1:nrow(D)) {
          for(j in 1:ncol(D)) {
            if(D[i,j] >= 0) {
              D_prime[i,j] <- D[i,j]/Dmax[i,j]
            } else {
              D_prime[i,j] <- D[i,j]/(-Dmin[i,j])
            }
          }
        }

        r_squared <- (D^2) / outer(p_A*(1-p_A), p_B*(1-p_B))

        chi_sq <- n * r_squared
        p_values <- 1 - pchisq(chi_sq, df = 1)

        allele_pairs <- expand.grid(allele1 = alleles1, allele2 = alleles2)
        ld_details <- data.frame(
          allele1 = allele_pairs$allele1,
          allele2 = allele_pairs$allele2,
          D = as.numeric(D),
          D_prime = as.numeric(D_prime),
          r_squared = as.numeric(r_squared),
          chi_sq = as.numeric(chi_sq),
          p_value = as.numeric(p_values)
        )

        ld_details <- ld_details[order(-ld_details$r_squared), ]

        mean_D_prime <- sum(abs(D_prime) * outer(p_A, p_B)) / sum(outer(p_A, p_B))
        mean_r_squared <- mean(r_squared, na.rm = TRUE)

        ld_results[[paste0(locus1, "-", locus2)]] <- list(
          locus1 = locus1,
          locus2 = locus2,
          contingency_table = cont_table,
          D_matrix = D,
          D_prime_matrix = D_prime,
          r_squared_matrix = r_squared,
          chi_sq_matrix = chi_sq,
          p_value_matrix = p_values,
          ld_details = ld_details,
          mean_D_prime = mean_D_prime,
          mean_r_squared = mean_r_squared,
          n_individuals = n,
          alleles_locus1 = alleles1,
          alleles_locus2 = alleles2
        )

      }, error = function(e) {
        message(paste("Erro no cálculo de LD para", locus1, "-", locus2, ":", e$message))
      })
    }

    return(ld_results)
  }

  calculate_complete_haplotypes <- function(hla_data) {
    tryCatch({
      loci <- c("A", "B", "DRB1")

      allele_cols <- list()
      for(locus in loci) {
        col1 <- paste0(locus, ifelse(locus == "DRB1", "_1", "1"))
        col2 <- paste0(locus, ifelse(locus == "DRB1", "_2", "2"))
        allele_cols[[locus]] <- c(col1, col2)
      }

      missing_cols <- setdiff(unlist(allele_cols), colnames(hla_data))
      if(length(missing_cols) > 0) {
        stop(paste("Colunas faltando:", paste(missing_cols, collapse = ", ")))
      }

      alleles <- list()
      for(locus in loci) {
        cols <- allele_cols[[locus]]
        allele1 <- gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[cols[1]]])
        allele2 <- gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[cols[2]]])

        valid <- !is.na(allele1) & !is.na(allele2) &
          allele1 != "" & allele2 != ""

        alleles[[paste0(locus, "1")]] <- allele1[valid]
        alleles[[paste0(locus, "2")]] <- allele2[valid]
      }

      n_individuals <- unique(sapply(alleles, length))
      if(length(n_individuals) != 1) {
        stop("Número inconsistente de indivíduos entre os loci")
      }
      n_individuals <- n_individuals[1]

      if(n_individuals < 10) {
        stop("Amostra muito pequena para cálculo de haplótipos (n < 10)")
      }

      haplotypes <- character(n_individuals * 8)
      index <- 1

      for(i in 1:n_individuals) {
        for(a in c(alleles[["A1"]][i], alleles[["A2"]][i])) {
          for(b in c(alleles[["B1"]][i], alleles[["B2"]][i])) {
            for(drb in c(alleles[["DRB11"]][i], alleles[["DRB12"]][i])) {
              haplotypes[index] <- paste(a, b, drb, sep = "-")
              index <- index + 1
            }
          }
        }
      }

      haplo_counts <- as.data.frame(table(haplotypes))
      haplo_counts$Frequency <- haplo_counts$Freq / sum(haplo_counts$Freq)
      haplo_counts$Percentage <- round(haplo_counts$Frequency * 100, 2)
      haplo_counts <- haplo_counts[order(-haplo_counts$Frequency), ]
      colnames(haplo_counts) <- c("Haplótipo", "Contagem", "Frequência", "Porcentagem")

      haplo_alleles <- strsplit(as.character(haplo_counts$Haplótipo), "-")
      haplo_counts$A <- sapply(haplo_alleles, `[`, 1)
      haplo_counts$B <- sapply(haplo_alleles, `[`, 2)
      haplo_counts$DRB1 <- sapply(haplo_alleles, `[`, 3)

      return(list(
        haplotype_counts = haplo_counts,
        n_individuals = n_individuals,
        n_haplotypes = nrow(haplo_counts)
      ))

    }, error = function(e) {
      stop(paste("Erro no cálculo de haplótipos:", e$message))
    })
  }

  run_logistic_regression <- function(data, outcome_var, group1, group2, allele_var, covariates = NULL) {
    tryCatch({
      df <- data.frame(data)

      count_alleles <- function(locus) {
        cols <- if (locus == "DRB1") c("DRB1_1", "DRB1_2") else c(paste0(locus, "1"), paste0(locus, "2"))
        alleles <- unique(na.omit(unlist(df[cols])))
        alleles <- alleles[alleles != "" & !is.na(alleles)]
        length(alleles)
      }

      current_locus <- gsub("\\*.*", "", allele_var)
      n_alleles <- count_alleles(current_locus)

      if (is.numeric(df[[outcome_var]])) {
        median_val <- median(df[[outcome_var]], na.rm = TRUE)
        df$outcome_binary <- as.integer(df[[outcome_var]] > median_val)
        group_labels <- c(paste("≤", round(median_val, 2)), paste(">", round(median_val, 2)))
      } else {
        if (!all(c(group1, group2) %in% unique(df[[outcome_var]]))) {
          stop("Grupos selecionados não existem na variável de desfecho")
        }
        df$outcome_binary <- ifelse(df[[outcome_var]] == group1, 0,
                                    ifelse(df[[outcome_var]] == group2, 1, NA))
        group_labels <- c(group1, group2)
      }

      df <- df[!is.na(df$outcome_binary), ]

      allele_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2")
      df$allele_present <- as.integer(
        df$A1 == allele_var | df$A2 == allele_var |
          df$B1 == allele_var | df$B2 == allele_var |
          df$DRB1_1 == allele_var | df$DRB1_2 == allele_var
      )

      if (!is.null(covariates)) {
        clean_covars <- list()

        for (covar in covariates) {
          if (!covar %in% colnames(df)) {
            stop(paste("Covariável não encontrada:", covar))
          }

          if (is.numeric(df[[covar]])) {
            clean_covars[[covar]] <- df[[covar]]
          } else {
            levels <- unique(na.omit(df[[covar]]))
            if (length(levels) != 2) {
              stop(paste("Covariável categórica deve ter exatamente 2 níveis:", covar))
            }
            clean_covars[[covar]] <- as.integer(df[[covar]] == levels[2])
          }
        }

        df <- cbind(df, as.data.frame(clean_covars))
      }

      if (is.null(covariates)) {
        formula <- as.formula("outcome_binary ~ allele_present")
      } else {
        formula <- as.formula(paste("outcome_binary ~ allele_present +",
                                    paste(covariates, collapse = " + ")))
      }

      model <- glm(formula, data = df, family = binomial())

      coefs <- summary(model)$coefficients
      or <- exp(coef(model))
      ci <- exp(confint.default(model))

      p_values <- coefs[-1, 4]
      p_corrected <- rep(NA, length(p_values))

      is_allele <- grepl("allele_present", names(p_values))
      p_corrected[is_allele] <- p_values[is_allele] * n_alleles
      p_corrected[!is_allele] <- p_values[!is_allele]

      p_corrected <- ifelse(p_corrected > 1, ">0.99",
                            ifelse(p_corrected < 0.0001, "<0.0001",
                                   round(p_corrected, 4)))

      results <- list(
        model = model,
        summary = data.frame(
          Variável = c(allele_var, covariates),
          OR = round(or[-1], 4),
          IC_95 = sapply(2:nrow(ci), function(i) {
            paste0("(", round(ci[i, 1], 3), "-", round(ci[i, 2], 3), ")")
          }),
          p_valor = ifelse(p_values < 0.0001, "<0.0001", round(p_values, 4)),
          Pc = p_corrected,
          stringsAsFactors = FALSE,
          row.names = NULL
        ),
        groups = group_labels,
        allele = allele_var,
        n_alleles = n_alleles,
        n = nrow(df),
        covariates = covariates
      )

      return(results)

    }, error = function(e) {
      stop(paste("Erro na regressão:", e$message))
    })
  }

  ui <- fluidPage(
    titlePanel("Análise HLA"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Carregar arquivo Excel", accept = c(".xlsx", ".xls")),
        actionButton("analyze", "Analisar Dados", class = "btn-primary"),
        downloadButton("download", "Exportar Resultados", class = "btn-success"),
        helpText("Analise frequências alélicas, equilíbrio de Hardy-Weinberg, desequilíbrio de ligação e haplótipos completos.")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Resumo",
                   verbatimTextOutput("summary"),
                   uiOutput("errorMessage")),
          tabPanel("Frequências",
                   h3("Frequências HLA-A"),
                   textOutput("freqSummaryA"),
                   DTOutput("freqTableA"),
                   h3("Frequências HLA-B"),
                   textOutput("freqSummaryB"),
                   DTOutput("freqTableB"),
                   h3("Frequências HLA-DRB1"),
                   textOutput("freqSummaryDRB1"),
                   DTOutput("freqTableDRB1")),
          tabPanel("Hardy-Weinberg",
                   tableOutput("hweTable"),
                   verbatimTextOutput("hweDetails")),
          tabPanel("Desequilíbrio",
                   h4("Resultados de Desequilíbrio de Ligação"),
                   uiOutput("ldTabs")),
          tabPanel("Haplótipos",
                   h3("Haplótipos Completos (A-B-DRB1)"),
                   verbatimTextOutput("haplotypeSummary"),
                   DTOutput("haplotypeTable")),
          tabPanel("Regressão Logística",
                   h3("Análise de Regressão Logística"),
                   selectInput("outcome_var", "Variável de Desfecho:", choices = NULL),
                   uiOutput("group_selection_ui"),
                   selectInput("allele_var", "Alelo para análise:", choices = NULL),
                   uiOutput("covariate_selection_ui"),
                   actionButton("run_regression", "Calcular", class = "btn-primary"),
                   tags$hr(),
                   h4("Resultados:"),
                   verbatimTextOutput("regression_summary"),
                   tableOutput("regression_table"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    analysis_results <- reactiveValues()
    last_error <- reactiveVal(NULL)

    output$summary <- renderPrint({
      req(analysis_results$data)
      cat("### Resumo da Análise HLA ###\n")
      cat("Indivíduos analisados:", nrow(analysis_results$data), "\n")
      cat("Loci analisados: A, B, DRB1\n")
      cat("\nUse as outras abas para ver os resultados detalhados.")
    })

    output$errorMessage <- renderUI({
      if(!is.null(last_error())) {
        tagList(
          h4("Erro na análise:", style = "color:red;"),
          p(last_error()),
          p("Verifique o formato do arquivo e tente novamente.")
        )
      }
    })

    output$freqSummaryA <- renderText({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$A
      if(is.null(freq_table)) return(NULL)
      paste("Total de alelos únicos HLA-A:", nrow(freq_table))
    })

    output$freqSummaryB <- renderText({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$B
      if(is.null(freq_table)) return(NULL)
      paste("Total de alelos únicos HLA-B:", nrow(freq_table))
    })

    output$freqSummaryDRB1 <- renderText({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$DRB1
      if(is.null(freq_table)) return(NULL)
      paste("Total de alelos únicos HLA-DRB1:", nrow(freq_table))
    })

    output$freqTableA <- renderDT({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$A
      if(is.null(freq_table)) return(NULL)

      datatable(
        freq_table,
        colnames = c("Alelo", "Contagem", "Frequência (%)"),
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100, "All"),
          dom = 'Blfrtip',
          buttons = c('copy', 'csv', 'excel'),
          ordering = TRUE,
          order = list(2, 'desc')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        selection = 'none'
      ) %>%
        formatRound("Percentage", 2) %>%
        formatStyle(
          'Percentage',
          background = styleColorBar(freq_table$Percentage, 'lightblue'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })

    output$freqTableB <- renderDT({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$B
      if(is.null(freq_table)) return(NULL)

      datatable(
        freq_table,
        colnames = c("Alelo", "Contagem", "Frequência (%)"),
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100, "All"),
          dom = 'Blfrtip',
          buttons = c('copy', 'csv', 'excel'),
          ordering = TRUE,
          order = list(2, 'desc')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        selection = 'none'
      ) %>%
        formatRound("Percentage", 2) %>%
        formatStyle(
          'Percentage',
          background = styleColorBar(freq_table$Percentage, 'lightblue'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })

    output$freqTableDRB1 <- renderDT({
      req(analysis_results$allele_freqs)
      freq_table <- analysis_results$allele_freqs$DRB1
      if(is.null(freq_table)) return(NULL)

      datatable(
        freq_table,
        colnames = c("Alelo", "Contagem", "Frequência (%)"),
        options = list(
          pageLength = 10,
          lengthMenu = c(10, 25, 50, 100, "All"),
          dom = 'Blfrtip',
          buttons = c('copy', 'csv', 'excel'),
          ordering = TRUE,
          order = list(2, 'desc')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        selection = 'none'
      ) %>%
        formatRound("Percentage", 2) %>%
        formatStyle(
          'Percentage',
          background = styleColorBar(freq_table$Percentage, 'lightblue'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })

    output$hweTable <- renderTable({
      req(analysis_results$hwe_results)

      results <- lapply(names(analysis_results$hwe_results), function(locus) {
        x <- analysis_results$hwe_results[[locus]]
        data.frame(
          Locus = locus,
          "Chi-Sq" = x$chisq,
          "p-value" = x$p.value,
          "df" = x$df,
          "Alelos" = x$n_alleles,
          "Indivíduos" = x$n_individuals,
          check.names = FALSE
        )
      })

      do.call(rbind, results)
    }, rownames = FALSE)

    output$hweDetails <- renderPrint({
      req(analysis_results$hwe_results)

      for(locus in names(analysis_results$hwe_results)) {
        cat("\n=== Locus", locus, "===\n")
        x <- analysis_results$hwe_results[[locus]]

        cat("\nTeste de Hardy-Weinberg:\n")
        cat("Chi-Sq =", x$chisq, ", df =", x$df, ", p-value =", x$p.value, "\n")
      }
    })

    output$ldTabs <- renderUI({
      req(analysis_results$ld_result)

      ld_tabs <- lapply(names(analysis_results$ld_result), function(pair_name) {
        tabPanel(
          pair_name,
          verbatimTextOutput(paste0("ldResults_", pair_name)),
          h4("Detalhes por Par de Alelos"),
          DTOutput(paste0("ldDetailsTable_", pair_name))
        )
      })

      do.call(tabsetPanel, ld_tabs)
    })

    observe({
      req(analysis_results$ld_result)

      for(pair_name in names(analysis_results$ld_result)) {
        local({
          local_pair_name <- pair_name
          ld_data <- analysis_results$ld_result[[local_pair_name]]

          output[[paste0("ldResults_", local_pair_name)]] <- renderPrint({
            cat("### Desequilíbrio de Ligação ###\n")
            cat("Entre", ld_data$locus1, "e", ld_data$locus2, "\n\n")
            cat("Indivíduos analisados:", ld_data$n_individuals, "\n")
            cat("Média de |D'|:", round(ld_data$mean_D_prime, 4), "\n")
            cat("Média de r²:", round(ld_data$mean_r_squared, 4), "\n\n")
            cat("Interpretação:\n")
            cat("- D' varia de -1 a 1, com |D'| > 0.8 indicando desequilíbrio forte\n")
            cat("- r² varia de 0 a 1, com r² > 0.33 indicando desequilíbrio forte\n")
          })

          output[[paste0("ldDetailsTable_", local_pair_name)]] <- renderDT({
            datatable(
              ld_data$ld_details,
              options = list(
                pageLength = 10,
                scrollX = TRUE
              ),
              rownames = FALSE
            ) %>%
              formatRound(columns = c("D", "D_prime", "r_squared", "chi_sq"), digits = 4) %>%
              formatRound(columns = "p_value", digits = 6) %>%
              formatStyle(
                'p_value',
                backgroundColor = styleInterval(0.05, c('lightgreen', 'white')))
          })
        })
      }
    })

    output$haplotypeSummary <- renderPrint({
      req(analysis_results$haplotype_result)
      cat("### Resumo de Haplótipos Completos ###\n")
      cat("Loci analisados: A-B-DRB1\n")
      cat("Indivíduos analisados:", analysis_results$haplotype_result$n_individuals, "\n")
      cat("Haplótipos únicos encontrados:", analysis_results$haplotype_result$n_haplotypes, "\n")
      cat("\nOs haplótipos são estimados considerando fase desconhecida (todas as combinações possíveis).\n")
    })

    output$haplotypeTable <- renderDT({
      req(analysis_results$haplotype_result)
      datatable(
        analysis_results$haplotype_result$haplotype_counts,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'Blfrtip',
          buttons = c('copy', 'csv', 'excel'),
          ordering = TRUE,
          order = list(3, 'desc')
        ),
        rownames = FALSE,
        extensions = 'Buttons',
        selection = 'none'
      ) %>%
        formatRound(columns = c("Frequência", "Porcentagem"), digits = 4) %>%
        formatStyle(
          'Porcentagem',
          background = styleColorBar(analysis_results$haplotype_result$haplotype_counts$Porcentagem, 'lightgreen'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })

    observeEvent(input$analyze, {
      req(input$file)

      tryCatch({
        last_error(NULL)
        showModal(modalDialog("Processando dados...", footer = NULL))

        hla_data <- read_hla_data(input$file$datapath)

        allele_freqs <- calculate_allele_frequencies(hla_data$alleles)

        hwe_results <- test_hwe(hla_data$data)

        ld_result <- calculate_all_ld(hla_data$data)

        haplotype_result <- calculate_complete_haplotypes(hla_data$data)

        analysis_results$data <- hla_data$data
        analysis_results$allele_freqs <- allele_freqs
        analysis_results$hwe_results <- hwe_results
        analysis_results$ld_result <- ld_result
        analysis_results$haplotype_result <- haplotype_result

        removeModal()
        showNotification("Análise concluída!", type = "message")

      }, error = function(e) {
        removeModal()
        last_error(e$message)
        showNotification(paste("Erro:", e$message), type = "error", duration = 10)
      })
    })

    output$download <- downloadHandler(
      filename = function() {
        paste("hla_results_", format(Sys.time(), "%Y%m%d_%H%M"), ".xlsx", sep = "")
      },
      content = function(file) {
        wb <- createWorkbook()

        addWorksheet(wb, "Frequencias")
        writeData(wb, "Frequencias", do.call(rbind, analysis_results$allele_freqs), rowNames = TRUE)

        hwe_df <- data.frame(
          Locus = names(analysis_results$hwe_results),
          ChiSq = sapply(analysis_results$hwe_results, function(x) x$chisq),
          p.value = sapply(analysis_results$hwe_results, function(x) x$p.value),
          df = sapply(analysis_results$hwe_results, function(x) x$df),
          n_alleles = sapply(analysis_results$hwe_results, function(x) x$n_alleles),
          n_individuals = sapply(analysis_results$hwe_results, function(x) x$n_individuals)
        )
        addWorksheet(wb, "Hardy-Weinberg")
        writeData(wb, "Hardy-Weinberg", hwe_df)

        hwe_details <- list()
        for(locus in names(analysis_results$hwe_results)) {
          x <- analysis_results$hwe_results[[locus]]
          hwe_details[[paste0(locus, "_observed")]] <- as.data.frame(x$observed)
          hwe_details[[paste0(locus, "_expected")]] <- as.data.frame(x$expected)
        }
        addWorksheet(wb, "HWE_Detalhes")
        writeData(wb, "HWE_Detalhes", hwe_details, rowNames = TRUE)

        addWorksheet(wb, "Desequilíbrio")
        ld_summary <- data.frame(
          Par = names(analysis_results$ld_result),
          Locus1 = sapply(analysis_results$ld_result, function(x) x$locus1),
          Locus2 = sapply(analysis_results$ld_result, function(x) x$locus2),
          Mean_D_prime = sapply(analysis_results$ld_result, function(x) x$mean_D_prime),
          Mean_r_squared = sapply(analysis_results$ld_result, function(x) x$mean_r_squared),
          N_Indivíduos = sapply(analysis_results$ld_result, function(x) x$n_individuals)
        )
        writeData(wb, "Desequilíbrio", ld_summary)

        addWorksheet(wb, "Haplótipos")
        writeData(wb, "Haplótipos", data.frame(
          N_Indivíduos = analysis_results$haplotype_result$n_individuals,
          N_Haplótipos = analysis_results$haplotype_result$n_haplotypes
        ), startRow = 1)

        writeData(wb, "Haplótipos", analysis_results$haplotype_result$haplotype_counts, startRow = 5)

        addWorksheet(wb, "Dados_Originais")
        writeData(wb, "Dados_Originais", analysis_results$data)

        saveWorkbook(wb, file)
      }
    )
    observe({
      req(analysis_results$data)
      data <- analysis_results$data

      hla_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2")
      outcome_options <- setdiff(colnames(data), hla_cols)
      updateSelectInput(session, "outcome_var", choices = outcome_options)

      allele_options <- unique(unlist(data[, hla_cols]))
      allele_options <- allele_options[!is.na(allele_options) & allele_options != ""]
      updateSelectInput(session, "allele_var", choices = sort(allele_options))

      covar_options <- setdiff(colnames(data), c(hla_cols, outcome_options))
      updateSelectInput(session, "covariates", choices = covar_options)
    })

    output$group_selection_ui <- renderUI({
      req(input$outcome_var, analysis_results$data)
      data <- analysis_results$data

      if (!input$outcome_var %in% colnames(data)) return(NULL)

      if (is.numeric(data[[input$outcome_var]])) {
        median_val <- median(data[[input$outcome_var]], na.rm = TRUE)
        tagList(
          helpText(paste("Variável numérica - divisão pela mediana (", round(median_val, 2), "):")),
          textInput("group1_name", "Grupo 1 (Referência):",
                    value = paste("≤", round(median_val, 2))),
          textInput("group2_name", "Grupo 2 (Comparação):",
                    value = paste(">", round(median_val, 2)))
        )
      } else {
        values <- unique(data[[input$outcome_var]])
        values <- values[!is.na(values)]
        tagList(
          selectInput("group1", "Grupo 1 (Referência):", choices = values),
          selectInput("group2", "Grupo 2 (Comparação):", choices = values,
                      selected = ifelse(length(values) >= 2, values[2], values[1]))
        )
      }
    })

    regression_results <- eventReactive(input$run_regression, {
      req(analysis_results$data, input$outcome_var, input$allele_var)

      if (is.numeric(analysis_results$data[[input$outcome_var]])) {
        group1 <- NA
        group2 <- NA
      } else {
        group1 <- input$group1
        group2 <- input$group2
      }

      run_logistic_regression(
        data = analysis_results$data,
        outcome_var = input$outcome_var,
        group1 = group1,
        group2 = group2,
        allele_var = input$allele_var,
        covariates = input$covariates
      )
    })

    output$regression_summary <- renderPrint({
      req(regression_results())
      res <- regression_results()

      cat("=== REGRESSÃO LOGÍSTICA ===\n")
      cat("Alelo analisado:", res$allele, "\n")
      cat("Grupo de referência (0):", res$groups[1], "\n")
      cat("Grupo de comparação (1):", res$groups[2], "\n")
      cat("Amostra utilizada:", res$n, "indivíduos\n\n")
      cat("Método: Modelo linear generalizado (família binomial)\n")
    })

    output$regression_table <- renderTable({
      req(regression_results())
      res <- regression_results()$summary

      res$Variável <- gsub("_", " ", res$Variável)
      res$Variável <- tools::toTitleCase(res$Variável)

      res <- res[, c("Variável", "OR", "IC_95", "p_valor", "Pc")]

      output$regression_note <- renderText({
        paste("Correção de Bonferroni aplicada multiplicando pelo número de alelos no locus (n =",
              regression_results()$n_alleles, ")")
      })

      res
    }, striped = TRUE, digits = 4, align = 'lrrrr')


    output$covariate_selection_ui <- renderUI({
      req(analysis_results$data, input$outcome_var)

      data <- analysis_results$data
      hla_cols <- c("A1", "A2", "B1", "B2", "DRB1_1", "DRB1_2")

      covar_options <- setdiff(
        colnames(data),
        c(hla_cols, input$outcome_var)
      )

      covar_options <- covar_options[sapply(data[, covar_options], function(x) {
        is.numeric(x) || length(unique(na.omit(x))) == 2
      })]

      if (length(covar_options) == 0) {
        return(helpText("Nenhuma covariável adequada encontrada (precisa ser numérica ou binária)"))
      }

      var_types <- sapply(data[, covar_options], function(x) {
        if (is.numeric(x)) " (numérico)" else " (binário)"
      })

      choices <- setNames(covar_options, paste0(covar_options, var_types))

      checkboxGroupInput(
        "covariates",
        "Selecione as covariáveis:",
        choices = choices,
        inline = FALSE
      )
    })
  }

  shinyApp(ui = ui, server = server)
}
