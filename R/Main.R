#' Main function of the HLAStcs package
#'
#' This function performs the main analysis of the package.
#' @export
HLAStcs <- function() {
  required_packages <- c(
    "openxlsx", "shiny", "DT", "shinyjs"
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
  library(shinyjs)

  read_hla_data <- function(file_path) {
    tryCatch({
      if(!file.exists(file_path)) {
        stop("File not found. Please check the path.")
      }

      hla_data <- read.xlsx(file_path)

      if(is.null(hla_data) || nrow(hla_data) == 0) {
        stop("Empty file or invalid format.")
      }

      return(list(data = hla_data))

    }, error = function(e) {
      stop(paste("Error reading file:", e$message))
    })
  }

  calculate_allele_frequencies <- function(hla_data, selected_loci) {
    if(is.null(hla_data)) stop("Missing allele data")

    freq_list <- list()

    for(locus in selected_loci) {
      col1 <- paste0(locus, "_1")
      col2 <- paste0(locus, "_2")

      if(!all(c(col1, col2) %in% colnames(hla_data))) {
        stop(paste("Columns for locus", locus, "not found (expected:", col1, "and", col2, ")"))
      }

      alleles <- c(trimws(as.character(hla_data[[col1]])), trimws(as.character(hla_data[[col2]])))
      alleles <- alleles[!is.na(alleles) & alleles != "" & alleles != "NA"]

      if(length(alleles) == 0) next

      freq_table <- as.data.frame(table(alleles))
      if(nrow(freq_table) == 0) next

      freq_table$Frequency <- freq_table$Freq / sum(freq_table$Freq)
      freq_table$Percentage <- round(freq_table$Frequency * 100, 2)
      freq_table$Count <- freq_table$Freq
      freq_table <- freq_table[order(-freq_table$Frequency), c("alleles", "Count", "Percentage")]

      freq_list[[locus]] <- freq_table
    }

    if(length(freq_list) == 0) stop("Could not calculate frequencies")
    return(freq_list)
  }

  test_hwe <- function(hla_data, selected_loci) {
    hwe_results <- list()

    for(locus in selected_loci) {
      col1 <- paste0(locus, "_1")
      col2 <- paste0(locus, "_2")

      if(!all(c(col1, col2) %in% colnames(hla_data))) {
        warning(paste("Columns", col1, "or", col2, "not found"))
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
        warning(paste("Sample size too small for", locus, "(n =", length(allele1), ")"))
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
        warning(paste("Not enough diversity to test HWE in", locus))
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

  calculate_all_ld <- function(hla_data, selected_loci) {
    if(length(selected_loci) < 2) return(NULL)

    ld_results <- list()
    locus_pairs <- combn(selected_loci, 2, simplify = FALSE)

    for(pair in locus_pairs) {
      locus1 <- pair[1]
      locus2 <- pair[2]

      tryCatch({
        col1_1 <- paste0(locus1, "_1")
        col1_2 <- paste0(locus1, "_2")
        col2_1 <- paste0(locus2, "_1")
        col2_2 <- paste0(locus2, "_2")

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
        message(paste("Error calculating LD for", locus1, "-", locus2, ":", e$message))
      })
    }

    return(ld_results)
  }

  calculate_complete_haplotypes <- function(hla_data, selected_loci) {
    tryCatch({
      if(length(selected_loci) < 2) return(NULL)

      allele_cols <- list()
      for(locus in selected_loci) {
        col1 <- paste0(locus, "_1")
        col2 <- paste0(locus, "_2")

        if(!all(c(col1, col2) %in% colnames(hla_data))) {
          stop(paste("Missing columns for locus", locus, ":",
                     paste(setdiff(c(col1, col2), colnames(hla_data)), collapse = ", ")))
        }

        allele_cols[[locus]] <- c(col1, col2)
      }

      haplotypes <- character()

      for(i in 1:nrow(hla_data)) {
        haplo1 <- sapply(selected_loci, function(locus) {
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[allele_cols[[locus]][1]]][i])
        })

        haplo2 <- sapply(selected_loci, function(locus) {
          gsub("([A-Z]+\\*\\d+):?\\d*", "\\1", hla_data[[allele_cols[[locus]][2]]][i])
        })

        if(!any(is.na(haplo1))) {
          haplotypes <- c(haplotypes, paste(haplo1, collapse = "-"))
        }

        if(!any(is.na(haplo2))) {
          haplotypes <- c(haplotypes, paste(haplo2, collapse = "-"))
        }
      }

      haplotypes <- haplotypes[!is.na(haplotypes) & haplotypes != "" & !grepl("NA", haplotypes)]

      if(length(haplotypes) == 0) {
        stop("No valid haplotypes found after processing")
      }

      haplo_counts <- as.data.frame(table(haplotypes))
      haplo_counts$Frequency <- haplo_counts$Freq / sum(haplo_counts$Freq)
      haplo_counts$Percentage <- round(haplo_counts$Frequency * 100, 2)
      haplo_counts <- haplo_counts[order(-haplo_counts$Frequency), ]
      colnames(haplo_counts) <- c("Haplotype", "Count", "Frequency", "Percentage")

      haplo_alleles <- strsplit(as.character(haplo_counts$Haplotype), "-")
      for(i in seq_along(selected_loci)) {
        haplo_counts[[selected_loci[i]]] <- sapply(haplo_alleles, `[`, i)
      }

      return(list(
        haplotype_counts = haplo_counts,
        n_individuals = nrow(hla_data),
        n_haplotypes = nrow(haplo_counts),
        loci = selected_loci
      ))

    }, error = function(e) {
      stop(paste("Error calculating haplotypes:", e$message))
    })
  }

  run_logistic_regression <- function(data, outcome_var, group1, group2, allele_var, covariates = NULL) {
    tryCatch({
      df <- data.frame(data)

      count_alleles <- function(locus) {
        cols <- paste0(locus, c("_1", "_2"))
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
          stop("Selected groups don't exist in the outcome variable")
        }
        df$outcome_binary <- ifelse(df[[outcome_var]] == group1, 0,
                                    ifelse(df[[outcome_var]] == group2, 1, NA))
        group_labels <- c(group1, group2)
      }

      df <- df[!is.na(df$outcome_binary), ]

      allele_cols <- paste0(current_locus, c("_1", "_2"))
      df$allele_present <- as.integer(
        df[[allele_cols[1]]] == allele_var | df[[allele_cols[2]]] == allele_var
      )

      if (!is.null(covariates)) {
        clean_covars <- list()

        for (covar in covariates) {
          if (!covar %in% colnames(df)) {
            stop(paste("Covariate not found:", covar))
          }

          if (is.numeric(df[[covar]])) {
            clean_covars[[covar]] <- df[[covar]]
          } else {
            levels <- unique(na.omit(df[[covar]]))
            if (length(levels) != 2) {
              stop(paste("Categorical covariate must have exactly 2 levels:", covar))
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
      suppressMessages(ci <- exp(confint(model)))
      aic <- AIC(model)

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
          Variable = c(allele_var, covariates),
          OR = round(or[-1], 4),
          CI_95 = sapply(2:nrow(ci), function(i) {
            paste0("(", round(ci[i, 1], 3), "-", round(ci[i, 2], 3), ")")
          }),
          p_value = ifelse(p_values < 0.0001, "<0.0001", round(p_values, 4)),
          Pc = p_corrected,
          stringsAsFactors = FALSE,
          row.names = NULL
        ),
        groups = group_labels,
        allele = allele_var,
        n_alleles = n_alleles,
        n = nrow(df),
        covariates = covariates,
        aic = aic
      )

      return(results)

    }, error = function(e) {
      stop(paste("Regression error:", e$message))
    })
  }

  ui <- fluidPage(
    useShinyjs(),
    titlePanel("HLA Analysis"),
    sidebarLayout(
      sidebarPanel(
        fileInput("file", "Load Excel file", accept = c(".xlsx", ".xls")),

        h4("HLA Column Mapping"),
        div(id = "loci_mapping",
            uiOutput("loci_mapping_ui")
        ),
        fluidRow(
          column(6, actionButton("add_locus", "Add Locus", icon = icon("plus"), width = "100%")),
          column(6, actionButton("remove_locus", "Remove Locus", icon = icon("minus"), width = "100%"))
        ),

        actionButton("analyze", "Analyze Data", class = "btn-primary",
                     style = "width: 100%; height: 50px; font-size: 16px;"),
        helpText("Analyze allele frequencies, Hardy-Weinberg equilibrium, linkage disequilibrium and complete haplotypes.")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Frequencies",
                   uiOutput("freq_tables_ui")),
          tabPanel("Hardy-Weinberg",
                   tableOutput("hweTable"),
                   verbatimTextOutput("hweDetails")),
          tabPanel("Linkage Disequilibrium",
                   h4("Linkage Disequilibrium Results"),
                   uiOutput("ldTabs")),
          tabPanel("Haplotypes",
                   h3("Complete Haplotypes"),
                   verbatimTextOutput("haplotypeSummary"),
                   DTOutput("haplotypeTable")),
          tabPanel("Logistic Regression",
                   h3("Logistic Regression Analysis"),
                   selectInput("outcome_var", "Outcome Variable:", choices = NULL),
                   uiOutput("group_selection_ui"),
                   selectInput("allele_var", "Allele for analysis:", choices = NULL),
                   uiOutput("covariate_selection_ui"),
                   actionButton("run_regression", "Calculate", class = "btn-primary"),
                   tags$hr(),
                   h4("Results:"),
                   verbatimTextOutput("regression_summary"),
                   tableOutput("regression_table"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    analysis_results <- reactiveValues()
    last_error <- reactiveVal(NULL)

    loci_mapping <- reactiveVal(list())

    observeEvent(input$add_locus, {
      current <- loci_mapping()
      new_id <- paste0("locus_", length(current) + 1)
      current[[new_id]] <- list(locus = "", col1 = "", col2 = "")
      loci_mapping(current)
    })

    observeEvent(input$remove_locus, {
      current <- loci_mapping()
      if(length(current) > 0) {
        current[[length(current)]] <- NULL
        loci_mapping(current)
      }
    })

    output$loci_mapping_ui <- renderUI({
      loci_list <- loci_mapping()
      if(length(loci_list) == 0) return(NULL)

      cols_available <- if(!is.null(analysis_results$raw_data)) {
        colnames(analysis_results$raw_data)
      } else {
        character(0)
      }

      tagList(
        lapply(names(loci_list), function(id) {
          fluidRow(
            column(4, textInput(paste0(id, "_name"), "Locus Name:",
                                value = loci_list[[id]]$locus)),
            column(4, selectInput(paste0(id, "_col1"), "Allele 1 Column:",
                                  choices = cols_available,
                                  selected = loci_list[[id]]$col1)),
            column(4, selectInput(paste0(id, "_col2"), "Allele 2 Column:",
                                  choices = cols_available,
                                  selected = loci_list[[id]]$col2))
          )
        })
      )
    })

    observeEvent(input$file, {
      tryCatch({
        hla_data <- read_hla_data(input$file$datapath)
        analysis_results$raw_data <- hla_data$data

        loci_mapping(list())

      }, error = function(e) {
        last_error(e$message)
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

    observeEvent(input$analyze, {
      req(analysis_results$raw_data)

      tryCatch({
        showModal(modalDialog("Processing data...", footer = NULL))

        selected_loci <- c()
        processed_data <- analysis_results$raw_data

        loci_inputs <- reactiveValuesToList(input)
        loci_inputs <- loci_inputs[grepl("^locus_", names(loci_inputs))]

        for(id in names(loci_mapping())) {
          locus_name <- input[[paste0(id, "_name")]]
          col1 <- input[[paste0(id, "_col1")]]
          col2 <- input[[paste0(id, "_col2")]]

          if(!is.null(locus_name) && !is.null(col1) && !is.null(col2) &&
             locus_name != "" && col1 != "" && col2 != "") {
            if(!all(c(col1, col2) %in% colnames(processed_data))) {
              stop(paste("Columns", col1, "or", col2, "not found in data"))
            }

            colnames(processed_data)[colnames(processed_data) == col1] <- paste0(locus_name, "_1")
            colnames(processed_data)[colnames(processed_data) == col2] <- paste0(locus_name, "_2")
            selected_loci <- c(selected_loci, locus_name)
          }
        }

        if(length(selected_loci) == 0) {
          stop("No valid loci were mapped. Please add at least one locus.")
        }

        allele_freqs <- calculate_allele_frequencies(processed_data, selected_loci)
        hwe_results <- test_hwe(processed_data, selected_loci)
        ld_result <- calculate_all_ld(processed_data, selected_loci)
        haplotype_result <- calculate_complete_haplotypes(processed_data, selected_loci)

        analysis_results$data <- processed_data
        analysis_results$allele_freqs <- allele_freqs
        analysis_results$hwe_results <- hwe_results
        analysis_results$ld_result <- ld_result
        analysis_results$haplotype_result <- haplotype_result
        analysis_results$selected_loci <- selected_loci

        removeModal()
        showNotification("Analysis complete!", type = "message")

      }, error = function(e) {
        removeModal()
        last_error(e$message)
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
      })
    })

    output$freq_tables_ui <- renderUI({
      req(analysis_results$allele_freqs, analysis_results$selected_loci)

      tagList(
        lapply(analysis_results$selected_loci, function(locus) {
          tagList(
            h3(paste("Frequencies", locus)),
            textOutput(paste0("freqSummary_", locus)),
            DTOutput(paste0("freqTable_", locus))
          )
        })
      )
    })

    observe({
      req(analysis_results$allele_freqs, analysis_results$selected_loci)

      for(locus in analysis_results$selected_loci) {
        local({
          local_locus <- locus
          freq_table <- analysis_results$allele_freqs[[local_locus]]

          output[[paste0("freqSummary_", local_locus)]] <- renderText({
            paste("Total unique alleles", local_locus, ":", nrow(freq_table))
          })

          output[[paste0("freqTable_", local_locus)]] <- renderDT({
            datatable(
              freq_table,
              colnames = c("Allele", "Count", "Frequency (%)"),
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
        })
      }
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
          "Alleles" = x$n_alleles,
          "Individuals" = x$n_individuals,
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

        cat("\nHardy-Weinberg Test:\n")
        cat("Chi-Sq =", x$chisq, ", df =", x$df, ", p-value =", x$p.value, "\n")
      }
    })

    output$ldTabs <- renderUI({
      req(analysis_results$ld_result)

      ld_tabs <- lapply(names(analysis_results$ld_result), function(pair_name) {
        tabPanel(
          pair_name,
          verbatimTextOutput(paste0("ldResults_", pair_name)),
          h4("Allele Pair Details"),
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
            cat("### Linkage Disequilibrium ###\n")
            cat("Between", ld_data$locus1, "and", ld_data$locus2, "\n\n")
            cat("Individuals analyzed:", ld_data$n_individuals, "\n")
            cat("Mean |D'|:", round(ld_data$mean_D_prime, 4), "\n")
            cat("Mean r²:", round(ld_data$mean_r_squared, 4), "\n\n")
            cat("Interpretation:\n")
            cat("- D' ranges from -1 to 1, with |D'| > 0.8 indicating strong disequilibrium\n")
            cat("- r² ranges from 0 to 1, with r² > 0.33 indicating strong disequilibrium\n")
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
      cat("### Complete Haplotypes Summary ###\n")
      cat("Analyzed loci:", paste(analysis_results$haplotype_result$loci, collapse = "-"), "\n")
      cat("Individuals analyzed:", analysis_results$haplotype_result$n_individuals, "\n")
      cat("Unique haplotypes found:", analysis_results$haplotype_result$n_haplotypes, "\n")
      cat("\nHaplotypes are estimated considering unknown phase (all possible combinations).\n")
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
        formatRound(columns = c("Frequency", "Percentage"), digits = 4) %>%
        formatStyle(
          'Percentage',
          background = styleColorBar(analysis_results$haplotype_result$haplotype_counts$Percentage, 'lightgreen'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })

    observe({
      req(analysis_results$data)
      data <- analysis_results$data

      hla_cols <- unlist(lapply(analysis_results$selected_loci, function(l) paste0(l, c("_1", "_2"))))
      outcome_options <- setdiff(colnames(data), hla_cols)
      updateSelectInput(session, "outcome_var", choices = outcome_options)

      allele_options <- unique(unlist(data[, hla_cols]))
      allele_options <- allele_options[!is.na(allele_options) & allele_options != ""]
      updateSelectInput(session, "allele_var", choices = sort(allele_options))

      covar_options <- setdiff(colnames(data), c(hla_cols, input$outcome_var))
      updateSelectInput(session, "covariates", choices = covar_options)
    })

    output$group_selection_ui <- renderUI({
      req(input$outcome_var, analysis_results$data)
      data <- analysis_results$data

      if (!input$outcome_var %in% colnames(data)) return(NULL)

      if (is.numeric(data[[input$outcome_var]])) {
        median_val <- median(data[[input$outcome_var]], na.rm = TRUE)
        tagList(
          helpText(paste("Numeric variable - split by median (", round(median_val, 2), "):")),
          textInput("group1_name", "Group 1 (Reference):",
                    value = paste("≤", round(median_val, 2))),
          textInput("group2_name", "Group 2 (Comparison):",
                    value = paste(">", round(median_val, 2)))
        )
      } else {
        values <- unique(data[[input$outcome_var]])
        values <- values[!is.na(values)]
        tagList(
          selectInput("group1", "Group 1 (Reference):", choices = values),
          selectInput("group2", "Group 2 (Comparison):", choices = values,
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

      cat("=== LOGISTIC REGRESSION ===\n")
      cat("Analyzed allele:", res$allele, "\n")
      cat("Reference group (0):", res$groups[1], "\n")
      cat("Comparison group (1):", res$groups[2], "\n")
      cat("Sample used:", res$n, "individuals\n")
      cat("AIC (Akaike Information Criterion):", round(res$aic, 2), "\n\n")
      cat("Method: Generalized linear model (binomial family)\n")
      cat("AIC interpretation:\n")
      cat("- Lower AIC indicates better model\n")
      cat("- Differences >2 between models are considered significant\n")
    })

    output$regression_table <- renderTable({
      req(regression_results())
      res <- regression_results()

      model_summary <- res$summary

      model_summary$Variable <- gsub("_", " ", model_summary$Variable)
      model_summary$Variable <- tools::toTitleCase(model_summary$Variable)

      model_summary <- model_summary[, c("Variable", "OR", "CI_95", "p_value", "Pc")]

      aic_row <- data.frame(
        Variable = "Model AIC",
        OR = round(res$aic, 2),
        CI_95 = "",
        p_value = "",
        Pc = ""
      )

      final_table <- rbind(model_summary, aic_row)

      output$regression_note <- renderText({
        paste("Bonferroni correction applied by multiplying by number of alleles in locus (n =",
              res$n_alleles, ")")
      })

      final_table
    }, striped = TRUE, digits = 4, align = 'lrrrr')

    output$covariate_selection_ui <- renderUI({
      req(analysis_results$data, input$outcome_var)

      data <- analysis_results$data
      hla_cols <- unlist(lapply(analysis_results$selected_loci, function(l) paste0(l, c("_1", "_2"))))

      covar_options <- setdiff(
        colnames(data),
        c(hla_cols, input$outcome_var)
      )

      covar_options <- covar_options[sapply(data[, covar_options], function(x) {
        is.numeric(x) || length(unique(na.omit(x))) == 2
      })]

      if (length(covar_options) == 0) {
        return(helpText("No suitable covariates found (must be numeric or binary)"))
      }

      var_types <- sapply(data[, covar_options], function(x) {
        if (is.numeric(x)) " (numeric)" else " (binary)"
      })

      choices <- setNames(covar_options, paste0(covar_options, var_types))

      checkboxGroupInput(
        "covariates",
        "Select covariates:",
        choices = choices,
        inline = FALSE
      )
    })
  }

  shinyApp(ui = ui, server = server)
}
