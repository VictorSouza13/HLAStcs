#' Main function of the HLAStcs package
#'
#' This function performs the main analysis of the package.
#' @export
HLAStcs <- function() {
  required_packages <- c(
    "openxlsx", "shiny", "DT", "shinyjs", "shinythemes"
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
  library(shinythemes)

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

  power_analysis <- function(A, B, C, D, alpha = 0.05, desired_power = 0.8) {
    OR <- (A * D) / (B * C)
    log_OR <- log(OR)
    se_log_OR <- sqrt(1/A + 1/B + 1/C + 1/D)
    z <- abs(log_OR) / se_log_OR
    current_power <- pnorm(z - qnorm(1 - alpha / 2))

    if (current_power >= desired_power) {
      return(paste0(
        "The estimated statistical power is ", round(current_power, 4),
        ", which meets or exceeds the desired threshold of ", desired_power, "."
      ))
    }

    p0 <- C / (C + D)
    r <- (A + B) / (C + D)

    simulate_power <- function(n1, n0) {
      p1 <- OR * p0 / (1 + p0 * (OR - 1))
      A1 <- round(n1 * p1); B1 <- n1 - A1
      C1 <- round(n0 * p0); D1 <- n0 - C1
      if (min(A1, B1, C1, D1) == 0) return(0)
      se_log <- sqrt(1/A1 + 1/B1 + 1/C1 + 1/D1)
      z <- abs(log(OR)) / se_log
      pnorm(z - qnorm(1 - alpha / 2))
    }

    for (n1 in 10:1000000) {
      n0 <- round(n1 / r)
      estimated_power <- simulate_power(n1, n0)
      if (estimated_power >= desired_power) {
        return(paste0(
          "The current power of this analysis is ", round(current_power, 4),
          ", which is below the desired threshold of ", desired_power,
          ". To achieve at least ", desired_power,
          " power, you would need approximately ", n1,
          " cases and ", n0, " controls. The estimated power with this sample size would be ",
          round(estimated_power, 4), "."
        ))
      }
    }

    return(paste0(
      "The current power is ", round(current_power, 4),
      ", and a sufficient sample size to achieve ", desired_power,
      " power could not be found within the tested range (up to 1,000,000 cases)."
    ))
  }

  run_logistic_regression <- function(data, outcome_var, reference_groups, comparison_groups, allele_var, covariates = NULL) {
    tryCatch({
      if(is.null(data)) stop("Dados não fornecidos")
      if(!outcome_var %in% names(data)) stop("Variável de resultado não encontrada")

      df <- as.data.frame(data, stringsAsFactors = FALSE)

      safe_count_alleles <- function(locus) {
        cols <- paste0(locus, c("_1", "_2"))
        if(!all(cols %in% names(df))) return(0)

        alleles <- unlist(df[cols])
        alleles <- alleles[!is.na(alleles) & alleles != ""]
        if(length(alleles) == 0) return(0)

        length(unique(alleles))
      }

      current_locus <- sub("\\*.*", "", allele_var)
      n_alleles <- safe_count_alleles(current_locus)

      ref_name <- paste(reference_groups, collapse = " + ")
      comp_name <- paste(setdiff(comparison_groups, reference_groups), collapse = " + ")

      df$analysis_group <- ifelse(
        df[[outcome_var]] %in% reference_groups, ref_name,
        ifelse(
          df[[outcome_var]] %in% comparison_groups, comp_name,
          NA_character_
        )
      )

      df <- df[!is.na(df$analysis_group), ]

      if(length(unique(df$analysis_group)) < 2) {
        stop("Não há grupos suficientes para comparação")
      }

      df$outcome_binary <- as.integer(df$analysis_group == comp_name)

      allele_cols <- paste0(current_locus, c("_1", "_2"))
      if(!all(allele_cols %in% names(df))) {
        stop("Colunas de alelo não encontradas nos dados")
      }

      df$allele_present <- as.integer(
        df[[allele_cols[1]]] == allele_var | df[[allele_cols[2]]] == allele_var
      )

      original_var_names <- list()
      if(!is.null(covariates)) {
        missing_covars <- setdiff(covariates, names(df))
        if(length(missing_covars) > 0) {
          stop(paste("Covariáveis não encontradas:", paste(missing_covars, collapse = ", ")))
        }

        for(covar in covariates) {
          if(!is.numeric(df[[covar]])) {
            original_var_names[[covar]] <- covar
            df[[covar]] <- factor(df[[covar]])
          }
        }
      }

      formula_terms <- "allele_present"
      if(!is.null(covariates)) {
        formula_terms <- c(formula_terms, covariates)
      }
      formula <- reformulate(termlabels = formula_terms, response = "outcome_binary")

      model <- tryCatch({
        glm(formula, data = df, family = binomial())
      }, error = function(e) {
        stop(paste("Falha ao ajustar modelo:", e$message))
      })

      safe_extract_results <- function(model, allele_name, n_alleles, original_names) {
        coefs <- summary(model)$coefficients
        or <- exp(coef(model))
        ci <- suppressMessages(exp(confint(model)))

        result_table <- data.frame(
          Term = rownames(coefs)[-1],
          OR = round(or[-1], 4),
          CI_lower = round(ci[-1, 1], 4),
          CI_upper = round(ci[-1, 2], 4),
          p_value = coefs[-1, 4],
          stringsAsFactors = FALSE
        )

        is_allele_term <- result_table$Term == "allele_present"
        result_table$Pc <- result_table$p_value
        result_table$Pc[is_allele_term] <- pmin(1, result_table$p_value[is_allele_term] * n_alleles)

        restore_original_names <- function(term) {
          if(term == "allele_present") return(allele_name)

          for(var in names(original_names)) {
            if(grepl(paste0("^", var), term)) {
              return(original_names[[var]])
            }
          }

          term
        }

        result_table$Term <- sapply(result_table$Term, restore_original_names)
        result_table <- result_table[!duplicated(result_table$Term), ]

        result_table
      }

      safe_contingency_table <- function(df) {
        tbl <- table(
          Group = df$analysis_group,
          Allele = factor(df$allele_present, levels = 0:1, labels = c("Absent", allele_var))
        )
        as.data.frame.matrix(tbl)
      }

      power_analysis <- function(A, B, C, D, comp_name, ref_name, alpha = 0.05, desired_power = 0.8) {
        if(A == 0 || B == 0 || C == 0 || D == 0) {
          return("Could not calculate power - insufficient data")
        }

        OR <- (A * D) / (B * C)
        log_OR <- log(OR)
        se_log_OR <- sqrt(1/A + 1/B + 1/C + 1/D)

        z_power <- qnorm(1 - alpha/2)
        current_power <- pnorm(abs(log_OR)/se_log_OR - z_power)

        if(current_power >= desired_power) {
          return(paste0(
            "Current power: ", round(current_power*100, 1), "% (adequate)"
          ))
        }

        p0 <- C / (C + D)
        p1 <- (OR * p0) / (1 + p0 * (OR - 1))

        n <- (qnorm(1-alpha/2) + qnorm(desired_power))^2 *
          (1/(p1*(1-p1)) + 1/(p0*(1-p0))) /
          (log(OR)^2)

        n <- ceiling(n)

        paste0(
          "Current power: ", round(current_power*100, 1), "%. ",
          "To achieve 80% power, approximately ",
          n, " total samples would be needed (",
          round(n * (A+B)/(A+B+C+D)), " ", comp_name,
          " and ", round(n * (C+D)/(A+B+C+D)), " ", ref_name, ")."
        )
      }

      contingency <- safe_contingency_table(df)
      power_result <- tryCatch({
        power_analysis(
          A = contingency[comp_name, allele_var],
          B = contingency[comp_name, "Absent"],
          C = contingency[ref_name, allele_var],
          D = contingency[ref_name, "Absent"],
          comp_name = comp_name,
          ref_name = ref_name
        )
      }, error = function(e) {
        paste("Não foi possível calcular o poder estatístico:", e$message)
      })

      list(
        results = safe_extract_results(model, allele_var, n_alleles, original_var_names),
        contingency = contingency,
        model_info = list(
          reference = ref_name,
          comparison = comp_name,
          n_samples = nrow(df),
          aic = round(AIC(model), 2),
          n_alleles = n_alleles,
          power = power_result
        )
      )

    }, error = function(e) {
      stop(paste("Erro na análise:", e$message))
    })
  }

  ui <- fluidPage(
    tags$head(
      tags$style(HTML("
    /* Increase font size for select boxes */
    select {
      font-size: 16px !important;
      height: auto !important;
      padding: 8px !important;
    }

    /* Increase font size for dropdown items */
    select option {
      font-size: 12px !important;
      padding: 8px !important;
    }

    /* Adjust dropdown size */
    .selectize-dropdown {
      font-size: 10px !important;
    }

    /* Increase font size for text inputs */
    input[type='text'] {
      font-size: 16px !important;
      height: auto !important;
      padding: 8px !important;
    }
  "))
    ),
    tags$head(
      tags$style(HTML("
      /* Main container adjustment */
      .container-fluid {
        max-width: 95%;
        width: 95%;
        margin-left: auto;
        margin-right: auto;
      }

      /* Adjustment for main area */
      .main-container {
        width: 100% !important;
      }

      /* Adjustment for main panel */
      .col-sm-9 {
        width: 65% !important;
      }

      /* Adjustment for side panel */
      .col-sm-3 {
        width: 35% !important;
      }

      /* Adjustment for tables */
      .dataTables_wrapper {
        width: 100% !important;
        margin-left: 0 !important;
        margin-right: 0 !important;
      }
    "))
    ),
    theme = shinytheme("flatly"),
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .well {
          background-color: #f8f9fa;
          border-radius: 5px;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        .btn-primary {
          background-color: #2c3e50;
          border-color: #2c3e50;
        }
        .btn-primary:hover {
          background-color: #1a252f;
          border-color: #1a252f;
        }
        .nav-tabs>li.active>a, .nav-tabs>li.active>a:focus, .nav-tabs>li.active>a:hover {
          background-color: #f8f9fa;
          border-bottom-color: transparent;
          color: #2c3e50;
          font-weight: bold;
        }
        h3, h4 {
          color: #2c3e50;
          margin-top: 20px;
        }
        .tab-content {
          background-color: white;
          padding: 15px;
          border-radius: 0 0 5px 5px;
          border: 1px solid #ddd;
          border-top: none;
        }
        .dataTables_wrapper .dataTables_info {
          color: #7b8a8b;
        }
        .shiny-notification {
          font-size: 16px;
          padding: 15px 20px;
        }
        .shiny-notification-error {
          background-color: #e74c3c;
          color: white;
        }
        .shiny-notification-success {
          background-color: #2ecc71;
          color: white;
        }
      "))
    ),

    titlePanel(
      div(
        icon("dna", class = "fa-2x", style = "margin-right: 15px; color: #2c3e50;"),
        "HLA/SNP Statistical Analysis Tool",
        style = "display: flex; align-items: center; color: #2c3e50;"
      )
    ),

    sidebarLayout(
      sidebarPanel(
        width = 3,
        wellPanel(
          style = "background-color: #f8f9fa;",
          h4("Data Input", icon("database")),
          fileInput("file", "",
                    accept = c(".xlsx", ".xls"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected")
        ),

        wellPanel(
          style = "background-color: #f8f9fa;",
          h4("Variant Mapping", icon("map-marked-alt")),
          div(id = "loci_mapping",
              uiOutput("loci_mapping_ui")
          ),
          fluidRow(
            column(6, actionButton("add_locus", "Add",
                                   icon = icon("plus-circle"),
                                   class = "btn-success",
                                   width = "100%")),
            column(6, actionButton("remove_locus", "Del",
                                   icon = icon("minus-circle"),
                                   class = "btn-danger",
                                   width = "100%"))
          )
        ),

        actionButton("analyze", "Run analysis",
                     icon = icon("play-circle"),
                     class = "btn-primary",
                     style = "background-color: #9b59b6; color: white; border: none;width: 100%; height: 50px; font-size: 24px; margin-top: 20px;"),
        helpText(icon("info-circle"),
                 "Analyze allele frequencies, Hardy-Weinberg equilibrium, linkage disequilibrium and complete haplotypes.")
      ),

      mainPanel(
        width = 9,
        tabsetPanel(
          type = "tabs",
          tabPanel(icon("chart-bar"), "Frequencies",
                   uiOutput("freq_tables_ui")),
          tabPanel(icon("balance-scale"), "Hardy-Weinberg",
                   div(class = "well",
                       h4("Hardy-Weinberg Equilibrium Results"),
                       tableOutput("hweTable"),
                       verbatimTextOutput("hweDetails"))),
          tabPanel(icon("link"), "Linkage Disequilibrium",
                   div(class = "well",
                       h4("Linkage Disequilibrium Results"),
                       uiOutput("ldTabs"))),
          tabPanel(icon("dna"), "Haplotypes",
                   div(class = "well",
                       h3(icon("dna"), "Complete Haplotypes"),
                       verbatimTextOutput("haplotypeSummary"),
                       DTOutput("haplotypeTable"))),
          tabPanel(icon("line-chart"), "Logistic Regression",
                   div(class = "well",
                       h3(icon("line-chart"), "Logistic Regression Analysis"),
                       fluidRow(
                         column(6, selectInput("outcome_var", "Outcome Variable:",
                                               choices = NULL, width = "100%")),
                         column(6, uiOutput("group_selection_ui"))
                       ),
                       fluidRow(
                         column(4, selectInput("locus_var", "Select Locus:",
                                               choices = NULL, width = "100%")),
                         column(4, uiOutput("allele_selection_ui")),
                         column(4, uiOutput("covariate_selection_ui"))
                       ),
                       actionButton("run_regression", "Run Regression",
                                    icon = icon("calculator"),
                                    class = "btn-primary"),
                       tags$hr(),
                       h4(icon("table"), "Results:"),
                       verbatimTextOutput("regression_summary"),
                       tableOutput("regression_table"),
                       uiOutput("regression_footer")
                   )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    analysis_results <- reactiveValues()
    last_error <- reactiveVal(NULL)

    loci_mapping <- reactiveVal(list())

    power_analysis <- function(A, B, C, D, comp_name, ref_name, alpha = 0.05, desired_power = 0.8) {
      if(A == 0 || B == 0 || C == 0 || D == 0) {
        return("Could not calculate power - insufficient data")
      }

      OR <- (A * D) / (B * C)
      log_OR <- log(OR)
      se_log_OR <- sqrt(1/A + 1/B + 1/C + 1/D)

      z_power <- qnorm(1 - alpha/2)
      current_power <- pnorm(abs(log_OR)/se_log_OR - z_power)

      if(current_power >= desired_power) {
        return(paste0(
          "Current power: ", round(current_power*100, 1), "% (adequate)"
        ))
      }

      p0 <- C / (C + D)
      p1 <- (OR * p0) / (1 + p0 * (OR - 1))

      n <- (qnorm(1-alpha/2) + qnorm(desired_power))^2 *
        (1/(p1*(1-p1)) + 1/(p0*(1-p0))) /
        (log(OR)^2)

      n <- ceiling(n)

      paste0(
        "Current power: ", round(current_power*100, 1), "%. ",
        "To achieve 80% power, approximately ",
        n, " total samples would be needed (",
        round(n * (A+B)/(A+B+C+D)), " ", comp_name,
        " and ", round(n * (C+D)/(A+B+C+D)), " ", ref_name, ")."
      )
    }

    observeEvent(input$add_locus, {
      current <- loci_mapping()
      new_id <- paste0("locus_", length(current) + 1)
      current[[new_id]] <- list(locus = "", col = "")
      loci_mapping(current)
    })

    observeEvent(input$remove_locus, {
      current <- loci_mapping()
      if(length(current) > 0) {
        current[[length(current)]] <- NULL
        loci_mapping(current)
      }
    })

    format_p_value <- function(p) {
      if(is.na(p)) return(NA)
      if(p == "<0.0001") return("<0.001")
      if(p == ">0.99") return(">0.99")
      if(grepl("<|>", p)) return(p)
      num_p <- suppressWarnings(as.numeric(p))
      if(is.na(num_p)) return(p)
      if(num_p < 0.001) return("<0.001")
      return(format(round(num_p, 3), nsmall = 3))
    }

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
            column(4,
                   div(style = "font-size: 12px;",
                       textInput(
                         inputId = paste0(id, "_name"),
                         label = "Locus Name:",
                         value = loci_list[[id]]$locus
                       )
                   )
            ),
            column(8,
                   div(style = "font-size: 12px;",
                       selectInput(
                         inputId = paste0(id, "_col"),
                         label = "Genotype Column:",
                         choices = cols_available,
                         selected = loci_list[[id]]$col
                       )
                   )
            )
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
          genotype_col <- input[[paste0(id, "_col")]]

          if(!is.null(locus_name) && !is.null(genotype_col) &&
             locus_name != "" && genotype_col != "") {

            if(!genotype_col %in% colnames(processed_data)) {
              stop(paste("Column", genotype_col, "not found in data"))
            }

            split_alleles <- strsplit(as.character(processed_data[[genotype_col]]), "/")

            processed_data[[paste0(locus_name, "_1")]] <- sapply(split_alleles, `[`, 1)
            processed_data[[paste0(locus_name, "_2")]] <- sapply(split_alleles, `[`, 2)

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

      if(!is.null(analysis_results$selected_loci)) {
        updateSelectInput(session, "locus_var", choices = analysis_results$selected_loci)
      }

      covar_options <- setdiff(colnames(data), c(hla_cols, input$outcome_var))
      updateSelectInput(session, "covariates", choices = covar_options)
    })

    output$allele_selection_ui <- renderUI({
      req(input$locus_var, analysis_results$data)

      locus <- input$locus_var
      col1 <- paste0(locus, "_1")
      col2 <- paste0(locus, "_2")

      if(all(c(col1, col2) %in% colnames(analysis_results$data))) {
        alleles <- unique(c(analysis_results$data[[col1]], analysis_results$data[[col2]]))
        alleles <- alleles[!is.na(alleles) & alleles != ""]

        selectInput("allele_var", "Select Allele:",
                    choices = sort(alleles),
                    width = "100%")
      }
    })

    output$group_selection_ui <- renderUI({
      req(input$outcome_var, analysis_results$data)
      data <- analysis_results$data

      if (!input$outcome_var %in% colnames(data)) return(NULL)

      if (is.numeric(data[[input$outcome_var]])) {
        quantiles <- quantile(data[[input$outcome_var]],
                              probs = seq(0, 1, 0.25),
                              na.rm = TRUE)
        groups <- c("Q1 (0-25%)", "Q2 (25-50%)", "Q3 (50-75%)", "Q4 (75-100%)")
        tagList(
          selectInput("reference_groups", "Reference Group(s):",
                      choices = groups, selected = groups[1],
                      multiple = TRUE),
          selectInput("comparison_groups", "Comparison Group(s):",
                      choices = groups, selected = groups[-1],
                      multiple = TRUE),
          helpText("Select one or more groups for each category")
        )
      } else {
        values <- unique(data[[input$outcome_var]])
        values <- values[!is.na(values)]
        tagList(
          selectInput("reference_groups", "Reference Group(s):",
                      choices = values, selected = values[1],
                      multiple = TRUE),
          selectInput("comparison_groups", "Comparison Group(s):",
                      choices = values, selected = values[-1],
                      multiple = TRUE),
          helpText("Select one or more groups for each category")
        )
      }
    })

    regression_results <- eventReactive(input$run_regression, {
      req(analysis_results$data, input$outcome_var, input$allele_var,
          input$reference_groups, input$comparison_groups)

      comparison_groups <- setdiff(input$comparison_groups, input$reference_groups)

      if(length(input$reference_groups) == 0) {
        stop("Please select at least one reference group")
      }
      if(length(comparison_groups) == 0) {
        stop("Please select at least one comparison group that's different from reference")
      }

      run_logistic_regression(
        data = analysis_results$data,
        outcome_var = input$outcome_var,
        reference_groups = input$reference_groups,
        comparison_groups = comparison_groups,
        allele_var = input$allele_var,
        covariates = input$covariates
      )
    })

    output$regression_summary <- renderPrint({
      req(regression_results())
      res <- regression_results()

      if(is.null(res$model_info) || is.null(res$contingency)) {
        stop("Invalid results structure")
      }

      cat("=== LOGISTIC REGRESSION ANALYSIS ===\n\n")
      cat("Reference Groups:", paste(input$reference_groups, collapse = " + "), "\n")
      cat("Comparison Groups:", paste(setdiff(input$comparison_groups, input$reference_groups), collapse = " + "), "\n")
      cat("Number of Samples:", res$model_info$n_samples, "\n")
      cat("Number of Unique Alleles:", res$model_info$n_alleles, "\n")
      cat("\nAllele Distribution:\n")
      print(res$contingency)
    })

    output$regression_table <- renderTable({
      req(regression_results())
      res <- regression_results()$results

      data.frame(
        Variable = res$Term,
        OR = sprintf("%.2f", res$OR),
        `95% CI` = sprintf("(%.2f-%.2f)", res$CI_lower, res$CI_upper),
        `p-value` = ifelse(res$p_value < 0.001, "<0.001", sprintf("%.3f", res$p_value)),
        `Pc` = ifelse(res$Pc < 0.001, "<0.001", sprintf("%.3f", res$Pc)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }, rownames = FALSE, striped = TRUE, hover = TRUE)

    observe({
      req(regression_results())
      res <- regression_results()

      for(comparison in res) {
        local({
          comp <- comparison
          output_name <- paste0("regression_table_", comp$comparison_group)

          output[[output_name]] <- renderTable({
            model_summary <- comp$summary

            model_summary$Variable <- gsub("_", " ", model_summary$Variable)
            model_summary$Variable <- tools::toTitleCase(model_summary$Variable)

            model_summary$`p-value` <- sapply(model_summary$p_value, format_p_value)
            model_summary$`Pc` <- sapply(model_summary$Pc, format_p_value)

            model_summary$`OR (95% CI)` <- paste0(
              format(round(model_summary$OR, 2), nsmall = 2),
              " (",
              gsub("[\\(\\)]", "", model_summary$CI_95),
              ")"
            )

            final_table <- model_summary[, c("Variable", "OR (95% CI)", "p-value", "Pc")]
            colnames(final_table) <- c("Variable", "OR (95% CI)", "p-value", "Pc*")

            final_table
          },
          align = "llrr",
          rownames = FALSE,
          striped = TRUE,
          width = "100%",
          hover = TRUE,
          bordered = TRUE)
        })
      }
    })

    output$regression_footer <- renderUI({
      req(regression_results())
      res <- regression_results()

      tagList(
        tags$div(
          style = "margin-top: 20px; background: #f8f9fa; padding: 10px; border-radius: 5px;",
          tags$p(tags$strong("Model AIC:"), tags$code(res$model_info$aic)),
          tags$p(tags$strong("Statistical Power:"), res$model_info$power),
          tags$p(tags$strong("Bonferroni Correction:"),
                 "Pc corrected for", res$model_info$n_alleles, "alleles"),
          tags$hr(),
          tags$p(tags$em("Interpretation:"),
                 paste("OR > 1 indicates association with",
                       res$model_info$comparison,
                       "group(s) compared to",
                       res$model_info$reference))
        )
      )
    })

    output$covariate_selection_ui <- renderUI({
      req(analysis_results$data, input$outcome_var)

      data <- analysis_results$data
      hla_cols <- unlist(lapply(analysis_results$selected_loci, function(l) paste0(l, c("_1", "_2"))))

      covar_options <- setdiff(
        colnames(data),
        c(hla_cols, input$outcome_var)
      )

      covar_options <- covar_options[sapply(data[, covar_options, drop = FALSE], function(x) {
        is.numeric(x) || length(unique(na.omit(x))) == 2
      })]

      if (length(covar_options) == 0) {
        return(helpText("No suitable covariates found (must be numeric or binary)"))
      }

      var_types <- sapply(data[, covar_options, drop = FALSE], function(x) {
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
