#' Uso de predição de HLA
#'
#' Esta função ativa o software.
#' @examples
#' InSilico()
#' @export
InSilico <- function() {
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }

  required_packages <- c(
    "readr", "httr", "bio3d", "tidyr", "dplyr", "NGLVieweR",
    "shiny", "ggplot2", "DT", "stringr", "gtsummary", "shinyjs",
    "reshape", "openxlsx"
  )

  invisible(sapply(required_packages, install_if_missing))

  alelosDisponiveis <- "HLA-C*12:03"

  if (interactive()) {
    suppressPackageStartupMessages({
      library(readr)
      library(httr)
      library(bio3d)
      library(tidyr)
      library(dplyr)
      library(NGLVieweR)
      library(shiny)
      library(ggplot2)
      library(DT)
      library(stringr)
      library(gtsummary)
      library(shinyjs)
      library(reshape)
      library(openxlsx)
    })

    ui <- fluidPage(
      titlePanel(""),

      tags$script(HTML('
        $(document).ready(function() {
          Shiny.addCustomMessageHandler("updateButtonColor", function(color) {
            $(".btn-analisar").css("background-color", color);
          });
        });
      ')),

      sidebarLayout(
        sidebarPanel(
          splitLayout(cellWidths = c("50%","50%"),
                      textInput(inputId = "alelo",
                                label = "Allele:",
                                value = "HLA-C*12:03"),
                      textInput("pdb.code", "PDB RCSB code (protein):", "6VYB")
          ),
          splitLayout(cellWidths = c("23%","23%","4%","50%"),
                      actionButton("remove", "Reset"),
                      actionButton("execute", "Run", class = "btn-analisar"),
                      uiOutput("circuloEstilizado"),
                      sliderInput("mer","Peptide size (mer):",min = 8,max = 30,value = 9)
          ),
        ),

        mainPanel(
          splitLayout(cellWidths = c("30%","50%","20%"),
                      dataTableOutput("table"),
                      NGLVieweROutput("structure"),
                      plotOutput("plot"))
        )
      )
    )

    server <- function(input, output, session) {
      output$analizarButtonColor <- renderUI({
        alelo <- input$alelo
        cor <- aleloCor(alelo)
        tags$style(HTML(sprintf(".btn-analisar {background-color: %s;}", cor)))
      })

      observeEvent(input$alelo, {
        cor <- aleloCor(input$alelo)
        session$sendCustomMessage(type = 'updateButtonColor', message = cor)
      })

      aleloCor <- function(alelo) {
        if (alelo %in% alelosDisponiveis) {
          return("green")
        } else {
          return("red")
        }
      }

      tipo <- "cartoon"

      output$structure <- renderNGLVieweR({
        NGLVieweR(input$pdb.code) %>%
          stageParameters(backgroundColor = "white") %>%
          setQuality("high") %>%
          setFocus(0) %>%
          setSpin(FALSE) %>%
          addRepresentation(tipo,
                            param = list(name = tipo,
                                         colorValue = "black",
                                         colorScheme = "uniform")
          )
      })

      observeEvent(input$execute, {
        pdb.code <- input$pdb.code

        is_nucleotide <- function(sequence) {
          if (grepl("U", sequence)) {
            return(TRUE)
          }
          valid_nucleotides <- c("A", "T", "C", "G")
          all(strsplit(sequence, "")[[1]] %in% valid_nucleotides)
        }

        is_amino_acid <- function(sequence) {
          amino_acid_letters <- c("B", "J", "O", "U", "X", "Z")
          !any(strsplit(sequence, "")[[1]] %in% amino_acid_letters)
        }

        fasta_file <- readLines(paste0("https://www.rcsb.org/fasta/entry/", pdb.code, "/display"))
        all_sequences <- c()
        amino_acid_sequences <- c()
        nucleotide_sequences <- c()

        for (line in fasta_file) {
          if (grepl("^>", line)) {
            next
          } else {
            all_sequences <- c(all_sequences, line)
            sequence <- gsub("-", "", line)
            sequence <- gsub("[0-9]", "", sequence)
            sequence <- gsub("[()d]", "", sequence)
            if (is_nucleotide(sequence)) {
              nucleotide_sequences <- c(nucleotide_sequences, sequence)
            } else if (is_amino_acid(sequence)) {
              amino_acid_sequences <- c(amino_acid_sequences, sequence)
            }
          }
        }

        Sprotein <- paste(amino_acid_sequences, collapse = "")
        ATGC <- paste(nucleotide_sequences, collapse = "")
        Total <- paste(all_sequences, collapse = "")

        positions <- gregexpr(ATGC, Total)
        start_positions <- unlist(positions)
        end_positions <- start_positions + nchar(ATGC) - 1

        sele5 <- paste(start_positions, end_positions, sep = "-")

        alelo <- input$alelo

        InSilico <- function(Epitope, alelo) {
          if (substr(alelo, 1, 5) == "HLA-D") {
            POST(url = "http://tools-cluster-interface.iedb.org/tools_api/mhcii/",
                 body = list(method = "nn_align",
                             sequence_text = Epitope,
                             allele = alelo,
                             length = input$mer)) -> rest
            ress <- content(rest, as = "text") %>% read_delim(delim = "\t")
          } else {
            POST(url = "http://tools-cluster-interface.iedb.org/tools_api/mhci/",
                 body = list(method = "smm",
                             sequence_text = Epitope,
                             allele = alelo,
                             length = input$mer)) -> rest
            ress <- content(rest, as = "text") %>% read_delim(delim = "\t")
          }
          return(ress)
        }

        temp <- InSilico(Sprotein, alelo)

        file_name <- paste0(pdb.code, "_ic50_data.xlsx")
        write.xlsx(temp, file = file_name)

        num_colunas <- temp$length[1] - 1
        for (i in 1:num_colunas) {
          nome_coluna <- paste0("p", i)
          temp[nome_coluna] <- temp$start + i
        }
        colnames(temp)[colnames(temp) == "start"] <- "p0"

        if (substr(alelo, 1, 5) == "HLA-D") {
          temp2 <- subset(temp,
                          select = c(-allele, -seq_num, -end, -length, -core_peptide, -peptide, -rank, -adjusted_rank))
        } else {
          temp2 <- subset(temp,
                          select = c(-allele, -seq_num, -end, -length, -peptide, -percentile_rank))
        }

        temp3 <- temp2 %>%
          gather(key = "peptide", value = "pos", -ic50) %>%
          group_by(pos) %>%
          summarize(avg_ic50 = mean(ic50))

        if (substr(temp$allele[1], 1, 5) == "HLA-D") {
          temp3$color <- ifelse(temp3$avg_ic50 >= 0 & temp3$avg_ic50 < 50, "red",
                                ifelse(temp3$avg_ic50 >= 50 & temp3$avg_ic50 < 1000, "yellow",
                                       ifelse(temp3$avg_ic50 >= 1000 & temp3$avg_ic50 < 5000,
                                              "blue",
                                              "black")))
        } else {
          temp3$color <- ifelse(temp3$avg_ic50 >= 0 & temp3$avg_ic50 < 50, "red",
                                ifelse(temp3$avg_ic50 >= 50 & temp3$avg_ic50 < 500, "yellow",
                                       ifelse(temp3$avg_ic50 >= 500 & temp3$avg_ic50 < 5000,
                                              "blue",
                                              "black")))
        }

        temp3 <- subset(temp3, select = c(-avg_ic50))

        collapse_seq <- function(x) {
          seqs <- split(x, cumsum(c(1, diff(x) != 1)))
          seq_strings <- sapply(seqs, function(seq) {
            if (length(seq) > 1) {
              paste0(min(seq), "-", max(seq))
            } else {
              as.character(seq)
            }
          })
          paste(seq_strings, collapse = " or ")
        }

        temp4 <- aggregate(pos ~ color, temp3, collapse_seq)
        sele1 <- temp4[temp4$color == "red", "pos"]
        sele2 <- temp4[temp4$color == "yellow", "pos"]
        sele3 <- temp4[temp4$color == "blue", "pos"]
        sele4 <- temp4[temp4$color == "black", "pos"]

        tabela <- function(alelo) {
          if (substr(alelo, 1, 5) == "HLA-D") {
            temp1 <- data.frame(allele = alelo,
                                nonbin = round((length(which(temp$ic50 >= 5000))), digits = 2),
                                weakbin = length(which(temp$ic50 >= 1000 & temp$ic50 < 5000)),
                                regbin = length(which(temp$ic50 >= 50 & temp$ic50 < 1000)),
                                strbin = length(which(temp$ic50 >= 0 & temp$ic50 < 50)))
          } else {
            temp1 <- data.frame(allele = alelo,
                                nonbin = round((length(which(temp$ic50 >= 5000))), digits = 2),
                                weakbin = length(which(temp$ic50 >= 500 & temp$ic50 < 5000)),
                                regbin = length(which(temp$ic50 >= 50 & temp$ic50 < 500)),
                                strbin = length(which(temp$ic50 >= 0 & temp$ic50 < 50)))
          }
          return(temp1)
        }

        my_colors <- c("Non Binders" = "black",
                       "Weak Binders" = "dark blue",
                       "Regular Binders" = "yellow",
                       "Strong Binders" = "red")

        tempS <- tabela(alelo)
        tempS <- reshape::melt(tempS, id.vars = "allele")
        tempS$value <- as.numeric(tempS$value)

        tempS2 <- data.frame(allele = alelo,
                             nonbin = (tempS$value[1]) / (sum(tempS$value)),
                             weakbin = (tempS$value[2]) / (sum(tempS$value)),
                             regbin = (tempS$value[3]) / (sum(tempS$value)),
                             strbin = (tempS$value[4]) / (sum(tempS$value)))

        tempS2 <- reshape::melt(tempS2, id.vars = "allele")
        tempS2$value <- as.numeric(tempS2$value)

        tempS2 %>%
          filter(variable %in% c("nonbin", "weakbin", "regbin", "strbin"), value != 0) %>%
          mutate(variable = recode(variable,
                                   "regbin" = "Regular Binders",
                                   "strbin" = "Strong Binders",
                                   "weakbin" = "Weak Binders",
                                   "nonbin" = "Non Binders")) %>%
          ggplot(aes(x = allele, y = value, fill = variable)) +
          geom_col(width = 0.5) +
          scale_y_continuous(limits = NULL, breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                             labels = c("0%", "25%", "50%", "75%", "100%")) +
          scale_fill_manual(values = my_colors) +
          theme(axis.ticks.x = element_blank(), aspect.ratio = 5/1,
                axis.text.x = element_blank(), axis.title = element_blank(),
                legend.title = element_blank(), legend.position = "bottom",
                legend.direction = "vertical",
                legend.text = element_text(size = 11),
                plot.title = element_text(hjust = 0.5, size = 8, face = "bold")) +
          labs(title = "", x = NULL, y = NULL) -> barra

        lista <- subset(temp, select = c(ic50, peptide))

        normalizar_mp <- function(mp) {
          mp_normalizada <- (mp - 0.25) / 0.75
          return(round(mp_normalizada, 2))
        }

        pesos <- c(0.25, 0.5, 0.75, 1.0)
        MP <- sum(tempS2$value[1:4] * pesos)
        MP_normalizada <- normalizar_mp(MP)

        output$plot <- renderPlot({
          barra +
            annotate("text", x = 1, y = 1.1, label = paste("MP:", MP_normalizada),
                     size = 5, hjust = 0.5, color = "blue")
        })

        output$table <- renderDataTable({
          datatable(lista,
                    options = list(
                      pageLength = 5,
                      lengthMenu = c(1, 3, 5),
                      dom = 'tp',
                      pagingType = 'simple',
                      order = list(0, 'asc')
                    ),
                    rownames = FALSE) %>%
            formatStyle(1:ncol(lista), fontSize = "6")
        })

        NGLVieweR_proxy("structure") %>%
          addSelection(tipo,
                       param = list(
                         name = "sel1",
                         sele = isolate(sele1),
                         colorValue = "red",
                         colorScheme = "uniform"
                       )
          ) %>%
          addSelection(tipo,
                       param = list(
                         name = "sel2",
                         sele = isolate(sele2),
                         colorValue = "yellow",
                         colorScheme = "uniform"
                       )
          ) %>%
          addSelection(tipo,
                       param = list(
                         name = "sel3",
                         sele = isolate(sele3),
                         colorValue = "blue",
                         colorScheme = "uniform"
                       )
          ) %>%
          addSelection(tipo,
                       param = list(
                         name = "sel4",
                         sele = isolate(sele4),
                         colorValue = "black",
                         colorScheme = "uniform"
                       )
          ) %>%
          addSelection(tipo,
                       param = list(
                         name = "sel5",
                         sele = isolate(sele5),
                         colorValue = "black",
                         colorScheme = "uniform"
                       )
          )
      })

      observeEvent(input$remove, {
        NGLVieweR_proxy("structure") %>% removeSelection("sel1") %>% removeSelection("sel2") %>%
          removeSelection("sel3") %>% removeSelection("sel4") %>% removeSelection("sel5")
      })
    }

    shinyApp(ui = ui, server = server)
  }
}
