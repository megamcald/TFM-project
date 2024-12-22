# Packages
# Shiny
library(shiny)
library(shinydashboard)
library(shinyjs)
library(rsconnect)
# Path
library(fs)
# Visualization: table, plot, cloud
library(DT)
library(wordcloud)
library(ggplot2)
# PubMed
library(pubmed.mineR)
library(easyPubMed)
# Text mining
library(lsa)
library(tokenizers)
# Others
library(viridis)

# Global config
# Directorios de datos
raw_data_dir <- file.path("data", "raw", "crohns_disease")
processed_data_dir <- file.path("data", "processed")
word_freq_dir <- file.path("data", "results", "word_frequency")
gene_analysis_dir <- file.path("data", "results", "gene_analysis")

# Crear directorios si no existen
dirs_to_create <- c(raw_data_dir, processed_data_dir, word_freq_dir, gene_analysis_dir)
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Functions
# Function 1. Generate PubMed query
generate_pubmed_query <- function(keywords, start_date, end_date) {
  # Ensure dates are in the correct format
  start_date <- format(as.Date(start_date), "%Y-%m-%d")
  end_date <- format(as.Date(end_date), "%Y-%m-%d")
  query <- paste(c(keywords, "[MH]", " AND ", start_date, ":", end_date, "[dp]"), collapse = "")
  return(query)
}

# Function 2. Download PubMed data
download_pubmed_data <- function(query, dest_dir, file_prefix) {
  # Ensure destination directory exists
  if (!dir.exists(dest_dir)) dir.create(dir, recursive = TRUE)
  batch_pubmed_download(pubmed_query_string = query, dest_dir, file_prefix, "abstract", batch_size = 5000)
}

# Function 3. Concatenate text files
concatenate_text_files <- function(files_list, output_file) {
  out_file <- file(description = output_file, open = "w")
  for (i in files_list) {
    x <- readLines(i)
    writeLines(x, out_file)
  }
  close(out_file)
}

# Function 4. Create and save barplots
create_bar_plot <- function(data, category, value, output_file, title) {
  # Create barplot
  plot <- ggplot(data, aes_string(x = category, y = value)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    labs(title = title, x = category, y = value) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Save barplot
  ggsave(output_file, plot = plot, width = 6.67, height = 6.67, dpi = 300)
}

# Function 5. Generate the HTML of the abstract
generate_abstract_html <- function(selected_pmid, abstracts_reactive, doi_reactive) {
  idx_in_corpus <- which(abstracts_reactive()@PMID == selected_pmid)

  selected_abstract <- abstracts_reactive()@Abstract[idx_in_corpus]
  selected_journal  <- abstracts_reactive()@Journal[idx_in_corpus]
  selected_doi      <- doi_reactive()[idx_in_corpus]

  abstract_sentences <- tokenize_sentences(selected_abstract, simplify = TRUE)

  to_print <- paste('<div class="abstract-card">', '\n')

  to_print <- paste(to_print,
                    '<div class="journal-title">', selected_journal, '</div>', '\n')

  for (i in seq_along(abstract_sentences)) {
    if (i <= 2) {
      to_print <- paste(to_print,
                        '<p class="abstract-sentence highlight">', abstract_sentences[i], '</p>', '\n')
    } else {
      to_print <- paste(to_print,
                        '<p class="abstract-sentence">', abstract_sentences[i], '</p>', '\n')
    }
  }

  pubmed_link <- paste0(
    '<p class="links"><a href="https://www.ncbi.nlm.nih.gov/pubmed/',
    selected_pmid,
    '" target="_blank"><i class="fas fa-pager"></i> Link PubMed</a></p>\n'
  )

  doi_link <- ""
  if (!is.na(selected_doi) && selected_doi != "") {
    doi_link <- paste0(
      '<p class="links"><a href="https://doi.org/',
      selected_doi,
      '" target="_blank"><i class="fas fa-link"></i> Link DOI</a></p>\n'
    )
  }

  to_print <- paste(to_print, pubmed_link, doi_link, '</div>')
  HTML(to_print)
}

# Function 6. Generate QuickGo y UniProt links
generate_gene_links <- function(gene_symbol) {
  quickgo_link <- paste0(
    '<p class="links"><a href="https://www.ebi.ac.uk/QuickGO/searchproducts/',
    gene_symbol,
    '?taxonId=9606" target="_blank"><i class="fas fa-link"></i> Link QuickGO</a></p>'
  )

  uniprot_link <- paste0(
    '<p class="links"><a href="https://www.uniprot.org/uniprotkb?query=',
    gene_symbol,
    '&facets=model_organism%3A9606" target="_blank"><i class="fas fa-link"></i> Link UniProt</a></p>'
  )

  HTML(
    paste0("<div class='abstract-card'>",
           "<div class='journal-title'>Información Gene_symbol: ", gene_symbol, "</div>",
           quickgo_link,
           uniprot_link,
           "</div>")
  )
}

### Script body: app ###
# UI
ui <- fluidPage(
  tags$head(
    tags$link(rel = "icon", type = "image/png", href = "logos/lupa.png"),
    # Font Awesome for icons in links
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css"),
    # CSS style file
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  useShinyjs(),
  tags$div(
    id = "loading_spinner",
    style = "display: none;
            position: fixed;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            text-align: center;",  # spinner loading (hide)
    tags$img(src = "spinnerv2.gif", alt = "Cargando...", width = 50, height = 50)
  ),
  # Header
  tags$div(
    class = "header",
    tags$img(src = "logos/logo_appv2.png", height = "175px", width = "auto", style = "margin-right: 15px;")
  ),

  titlePanel("PubExplorer", windowTitle = "PubExplorer: Text Mining Analysis"),

  navlistPanel(
    widths = c(2, 10),

    # PubMed Search Tab
    tabPanel(title = tagList(icon("search"), "Búsqueda en PubMed"),
             h1("Búsqueda en PubMed"),
             p("Puedes ingresar palabras clave y un rango de fechas para buscar artículos relevantes"),
             sidebarLayout(
               sidebarPanel(
                 h3("Parámetros de búsqueda"),
                 # Enter keywords
                 textInput("keywords", label = "Palabras clave", value = "Crohn's disease", placeholder = "E.g., Crohn's disease"),

                 # Enter dates
                 dateRangeInput("start_date", label = "Rango de fechas", start = Sys.Date() - (5 * 365), end = Sys.Date(), format = "dd-mm-yyyy", startview = "year", weekstart = 1, language = "es", separator = "hasta"),
                 p(),
                 p(strong("Consulta a PubMed")),
                 verbatimTextOutput("keyw"),

                 # Buttons for searching (2) search and new search
                 fluidRow(
                   column(6, align = "center",
                          actionButton("search", "Buscar en PubMed", class = "btn-primary")
                   ),
                   column(6, align = "center",
                          actionButton("new_search", "Nueva Búsqueda", class = "btn-secondary")
                   )
                 )
               ),

               column(8,
                      # Show message "Búsqueda completada"
                      textOutput("message"),

                      # Abstracts Table
                      DT::dataTableOutput("articles_table"),

                      # Show abstract
                      uiOutput("abstractText")
               )
             )
    ),

    # Word Frequency Tab
    # 3 display options: table, graph and word cloud
    tabPanel(title = tagList(icon("chart-bar"), "Palabras"),
             h1("Tabla Frecuencia Palabras"),
             p("Aquí se muestra una tabla con las palabras más frecuentes"),
             sidebarLayout(
               sidebarPanel(
                 selectInput("view_words", "Ver como:", choices = c("Tabla", "Gráfico", "Nube de palabras")),
                 numericInput("min_freq_words", "Frecuencia mínima:", value = 2),
                 numericInput("max_words", "Máximo número de palabras:", value = 50),
                 actionButton("draw_wordcloud", "Dibujar Nube")
               ),
               mainPanel(
                 conditionalPanel(
                   condition = "input.view_words == 'Tabla'",
                   DT::dataTableOutput("word_table")
                 ),
                 conditionalPanel(
                   condition = "input.view_words == 'Gráfico'",
                   plotOutput("word_plot")
                 ),
                 conditionalPanel(
                   condition = "input.view_words == 'Nube de palabras'",
                   plotOutput("wordcloud", height = "600px")
                 )
               )
             )
    ),

    # Gene Frequency Tab
    # 3 display options: table, graph and gene cloud
    tabPanel(title = tagList(icon("dna"), "Genes"),
             h1("Tabla Frecuencia de Genes"),
             p("Aquí se muestra una tabla con los genes más frecuentes"),
             sidebarLayout(
               sidebarPanel(
                 selectInput("view_genes", "Ver como:", choices = c("Tabla", "Gráfico", "Nube de genes")),
                 numericInput("min_freq_genes", "Frecuencia mínima:", value = 2),
                 numericInput("max_genes", "Máximo número de genes:", value = 50),
                 actionButton("draw_genecloud", "Dibujar Nube")
               ),
               mainPanel(
                 conditionalPanel(
                   condition = "input.view_genes == 'Tabla'",
                   DT::dataTableOutput("gene_table"),
                   uiOutput("gene_links")
                 ),
                 conditionalPanel(
                   condition = "input.view_genes == 'Gráfico'",
                   plotOutput("gene_plot")
                 ),
                 conditionalPanel(
                   condition = "input.view_genes == 'Nube de genes'",
                   plotOutput("genewordcloud", height = "600px")
                 )
               )
             )
    ),

    # Genes Navigation (genes)
    tabPanel(title = tagList(icon("search-plus"), "Navegación Genes"),
             h1("Navegación Genes"),
             p("Aquí se puede navegar por aquellas publicaciones que contienen determinados términos (Gene_symbol) en el abstract de los artículos relacionados con la enfermedad de Crohn."),
             htmlOutput("gene_info"),
             DT::dataTableOutput("filtered_articles"),
             htmlOutput("abstractNavText")
    ),

    # Term Navigation (palabras)
    tabPanel(
      title = tagList(icon("book"), "Navegación palabras"),
      h1("Navegación términos"),
      p("Aquí se puede navegar por aquellas publicaciones que contienen determinados términos (word) en el abstract de los artículos relacionados con la enfermedad de Crohn."),
      htmlOutput("word_info"),
      DT::dataTableOutput("filtered_term_articles"),
      htmlOutput("abstractTermText")
    ),

    # Co-ocurrencia of terms Tab
    tabPanel(
      title = tagList(icon("link"), "Co-ocurrencia de términos"),
      h1("Co-ocurrencia de términos"),
      sidebarLayout(
        sidebarPanel(
          # Selection Terms1 (words)
          uiOutput("cooc_term1_ui"),

          # Selection Terms2(genes)
          uiOutput("cooc_term2_ui"),

          # Parameter selection n = 0, 1, 2
          selectInput("cooc_n", "Proximidad de Co-ocurrencia (n):",
                      choices = list(
                        "Mismo párrafo" = 0,
                        "Párrafos consecutivos" = 1,
                        "Párrafos separados por uno" = 2
                      ),
                      selected = 0),

          # Button co-occurrence search
          actionButton("coocurrence_search", "Buscar co-ocurrencias")
        ),
        mainPanel(
          DT::dataTableOutput("coocurrence_table")
        )
      )
    ),

    # About Tab
    tabPanel(title = tagList(icon("info-circle"), "Acerca de"),
             h1("Acerca de PubExplorer"),
             p("Aplicación web para el estudio y análisis de la enfermedad de Crohn mediante minería de textos de PubMed."),
             h3("Objetivos"),
             p("El objetivo principal es el desarrollo de una aplicación web interactiva, que use el sistema de búsqueda de PubMed, y permita a los usuarios extraer la información contenida en los abstracts de las publicaciones, aplicando técnicas de minería de textos."),
             tags$ul(
               tags$li("Facilitar la búsqueda de artículos en PubMed sobre la enfermedad de Crohn."),
               tags$li("Proporcionar estadísticas sobre la frecuencia de palabras clave y genes asociados a la enfermedad."),
               tags$li("Análisis de co-ocurremcia entre términos: palabras clave - genes."),
               tags$li("Generar visualizaciones como nubes de palabras, gráficos y tablas para facilitar el análisis de datos."),
             ),
             h3("Paquetes y librerías utilizadas"),
             p("Para el desarrollo de la aplicación PubExplorer, se utilizó el lenguaje de programación R y diversos paquetes para las distintas funcionalidades."),
             tags$ul(
               tags$li("Shiny: Creación de la aplicación web interactiva."),
               tags$li("DT, wordcloud, ggplot2: Visualización de datos."),
               tags$li("pubmed.mineR, easyPubMed: Acceso y análisis de PubMed."),
               tags$li("lsa, tokenizers: Minería de textos"),
               tags$li("fs: Manejo de archivos y rutas")
             ),
             h3("Créditos"),
             p("Este proyecto fue desarrollado como parte del Trabajo Final del ",
               a("Máster en Bioinformática y Bioestadística en la UOC-UB (Universitat Oberta de Catalunya - Universitat de Barcelona)", href = "https://www.uoc.edu/es/estudios/masters/master-universitario-bioinformatica-bioestadistica", target = "_blank"),
             ),
             p("Desarrollado por: ",
               a(icon("linkedin"),"Megam Calderón", href = "https://www.linkedin.com/in/megam-calder%C3%B3n/", target = "_blank")
             ),

             h3("Referencias"),
             p("Repositorio del proyecto: ",
               a("PubExplorer", href = "https://github.com/megamcald", target = "_blank")
             ),
    )
  ),
  # Footer
  tags$footer(
    style = "background-color: #f8f9fa; padding: 10px; text-align: center; border-top: 1px solid #dee2e6;",
    tags$img(src = "logos/logo_appv2.png", height = "75px", style = "margin-right: 10px;"),
    tags$img(src = "logos/logo_uoc.png", height = "75px", style = "margin-left: 10px;"),
    tags$p("2024 | TFM Máster Universitario en Bioinformática y Bioestadística (UOC-UB)", style = "margin: 0; font-size: 14px; color: #555;")
  )
)


# SERVER
server <- function(input, output, session) {


  # reactiveVal
  doi_reactive <- reactiveVal(NULL)
  abstracts_reactive <- reactiveVal(NULL)
  word_tokens_reactive <- reactiveVal(NULL)
  top_words_reactive <- reactiveVal(NULL)
  top_words_fixed_reactive <- reactiveVal(NULL)
  genes_df_reactive <- reactiveVal(NULL)
  top_genes_reactive <- reactiveVal(NULL)
  selected_gene <- reactiveVal(NULL)
  selected_word <- reactiveVal(NULL)
  tfidf_reactive <- reactiveVal(NULL)

  # PubMed Query
  generate_query <- reactive({
    validate(
      need(input$keywords != "", "Es necesario ingresar una palabra clave."),
      need(input$start_date[1] <= input$start_date[2], "La fecha de inicio no puede ser posterior a la fecha de fin.")
    )

    query <- generate_pubmed_query(input$keywords, input$start_date[1], input$start_date[2])
    return(query)
  })

  # Show PubMed Query
  output$keyw <- renderText({
    generate_query()
  })

  # Logic for searching in PubMed
  observeEvent(input$search, {
    shinyjs::show("loading_spinner") # Show loading spinner

    query <- generate_query()
    dest_dir <- file.path(raw_data_dir)
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

    download_pubmed_data(query, dest_dir, file_prefix = "crohn_disease")
    shinyjs::hide("loading_spinner") # Hide loading spinner
    output$message <- renderText({
      "Búsqueda completada."
    })

    # Concatenate text files
    raw_files <- dir_ls(file.path(raw_data_dir),
                        regexp = "crohn_disease")
    concatenate_text_files(raw_files, output_file = file.path(processed_data_dir, "all_abstracts.txt"))

    # Read the abstracts
    # readabs
    abstracts <- readabs(file.path(processed_data_dir, "all_abstracts.txt"))
    corpus_primario <- abstracts
    abstracts_reactive(corpus_primario)
    # Delete unnecessary txt files
    raw_files_to_delete <- dir_ls(path = dest_dir, pattern = "\\.txt$")
    file.remove(raw_files_to_delete)

    # Extract DOI using regex
    # DOI Link
    doi_pattern <- "doi:\\s?(10[.][0-9]+/[A-Za-z0-9._%+-]+)"
    doisc <- sapply(corpus_primario@Abstract, function(abstract) {
      match <- regexpr(doi_pattern, abstract, perl = TRUE, ignore.case = TRUE)
      if (match != -1) {
        doi <- regmatches(abstract, match)
        # removing ‘doi:’ and spaces
        doi_clean <- sub("(?i)doi:\\s*", "", doi, perl = TRUE)
        doi_clean
      } else {
        NA
      }
    })

    doi_reactive(doisc)

    # df PMID, Publicación
    articles_table <- data.frame(
      PMID = corpus_primario@PMID,
      Publicación = corpus_primario@Journal,
      stringsAsFactors = FALSE
    )

    # note: Spanish.json to display content in Spanish
    output$articles_table <- DT::renderDataTable({
      datatable(articles_table,
                selection = list(mode = 'single', selected = 1),
                options = list(language = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")))
    })

    # word_tokens
    # word_atomizations function
    word_tokens <- word_atomizations(abstracts)
    word_tokens_reactive(word_tokens)
    save(word_tokens, file = file.path(processed_data_dir, "word_tokens.RData"))

    # word_table.txt
    write.table(word_tokens, file = file.path(word_freq_dir, "word_tokens.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # top 10 words
    top_words <- head(word_tokens, 10)
    top_words_fixed_reactive(top_words) # co-ocurrence
    top_words_reactive(top_words)       # wordcloud, barplot
    create_bar_plot(top_words, "words", "Freq", file.path(word_freq_dir, "word_tokens_barplot.png"), "Palabras más frecuentes")

    # word table
    output$word_table <- DT::renderDataTable({
      req(word_tokens_reactive())
      datatable(word_tokens_reactive(),
                selection = list(mode = 'single', selected = 1),
                options = list(pageLength = 10, language = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")))
    })

    # word barplot
    output$word_plot <- renderPlot({
      ggplot(top_words_fixed_reactive(), aes(x = reorder(words, Freq), y = Freq)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        theme_classic() +
        labs(x = "Palabras", y = "Frecuencia") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })

    # genes (Gene_symbol)
    # genes gene_atomization function
    genes <- gene_atomization(abstracts)
    genes_df <- as.data.frame(genes)
    colnames(genes_df) <- c("Gene_symbol", "Genes", "Freq")
    genes_df$Freq <- as.numeric(genes_df$Freq)

    # genes_df.txt
    genes_df_reactive(genes_df)
    save(genes_df, file = file.path(processed_data_dir, "genes_df.RData"))
    write.table(genes_df,
                file = file.path(gene_analysis_dir, "genes_df.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # top 10 genes
    top_genes <- head(genes_df, 10)
    top_genes_reactive(top_genes)
    create_bar_plot(top_genes, "Gene_symbol", "Freq",
                    file.path(gene_analysis_dir, "genes_barplot.png"),
                    "Genes más frecuentes")

    output$gene_table <- DT::renderDataTable({
      req(genes_df_reactive())
      datatable(genes_df_reactive(),
                selection = list(mode = 'single', selected = 1),
                options = list(pageLength = 10, language = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")))
    })

    output$gene_plot <- renderPlot({
      req(top_genes_reactive())
      ggplot(top_genes_reactive(), aes(x = reorder(Gene_symbol, Freq), y = Freq)) +
        geom_bar(stat = "identity", fill = "skyblue") +
        theme_classic() +
        labs(x = "Genes", y = "Frecuencia") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
  })

  # Display the selected abstract with PubMed and DOI links
  output$abstractText <- renderUI({
    row_selected <- input$articles_table_rows_selected
    if (length(row_selected) == 0) {
      return(NULL)
    }
    selected_pmid <- abstracts_reactive()@PMID[row_selected]

    generate_abstract_html(selected_pmid, abstracts_reactive, doi_reactive)
  })


  # Disable the search button if entries are invalid
  observe({
    # Check if dates are correctly defined
    if (input$keywords == "" || is.null(input$start_date) || is.null(input$start_date[1]) || is.null(input$start_date[2]) || input$start_date[1] > input$start_date[2]) {
      shinyjs::disable("search")
    } else {
      shinyjs::enable("search")
    }
  })

  # Logic for a new PubMed search
  observeEvent(input$new_search, {
    # Reset to default values
    updateTextInput(session, "keywords", value = "Crohn's disease")
    updateDateRangeInput(session, "start_date",
                         start = Sys.Date() - (2 * 365),
                         end = Sys.Date())

    # Remove results if it is a new search
    output$message <- renderText({})
    output$articles_table <- DT::renderDataTable({})
    output$abstractText <- renderUI({ NULL })
    output$word_table <- DT::renderDataTable({})
    output$word_plot <- renderPlot({})
    output$wordcloud <- renderPlot({})
    output$gene_table <- DT::renderDataTable({})
    output$gene_plot <- renderPlot({})
    output$genewordcloud <- renderPlot({})
    output$gene_links <- renderUI({ NULL })
    output$filtered_articles <- DT::renderDataTable({})
    output$abstractNavText <- renderUI({ NULL })
    output$word_info <- renderUI({ NULL })
    output$filtered_term_articles <- DT::renderDataTable({})
    output$abstractTermText <- renderUI({ NULL })
    output$coocurrence_table <- DT::renderDataTable({})
  })

  # Logic for enabling/disabling controls on the Word Frequency tab
  # are only activated if the wordcloud option is selected.
  observe({
    # Check the value of ‘view_words’
    if (input$view_words == "Nube de palabras") {
      shinyjs::enable("min_freq_words")
      shinyjs::enable("max_words")
      shinyjs::enable("draw_wordcloud")
    } else {
      shinyjs::disable("min_freq_words")
      shinyjs::disable("max_words")
      shinyjs::disable("draw_wordcloud")
    }
  })

  # Logic for enabling/disabling controls on the Gene Frequency tab
  # are only activated if the wordcloud option is selected.
  observe({
    # Check the value of ‘view_genes’
    if (input$view_genes == "Nube de genes") {
      shinyjs::enable("min_freq_genes")
      shinyjs::enable("max_genes")
      shinyjs::enable("draw_genecloud")
    } else {
      shinyjs::disable("min_freq_genes")
      shinyjs::disable("max_genes")
      shinyjs::disable("draw_genecloud")
    }
  })

  # Logic to draw the wordcloud when pressing the button
  observeEvent(input$draw_wordcloud, {
    req(word_tokens_reactive()) # Make sure word_tokens_reactive is not empty

    # Filter words by minimum frequency
    filtered_words <- word_tokens_reactive()[word_tokens_reactive()$Freq >= input$min_freq_words, ]

    # Select maximum number of words
    if (nrow(filtered_words) > input$max_words) {
      filtered_words <- head(filtered_words, input$max_words)
    }


    top_words_reactive(filtered_words)


    output$wordcloud <- renderPlot({
      if (nrow(filtered_words) == 0) {
        plot.new()
        text(0.5, 0.5, "No hay palabras que cumplan con los criterios.", cex = 1.2)
      } else {
        wordcloud(words = filtered_words$words,
                  freq = filtered_words$Freq,
                  min.freq = input$min_freq_words,
                  max.words = input$max_words,
                  scale = c(3, 2), # word size
                  random.order = FALSE,
                  rot.per = 0.35,
                  colors = viridis(n = input$max_words)) # viridis colour palette
      }
    })
  })

  # Logic to draw the gencloud when pressing the button (= wordcloud)
  observeEvent(input$draw_genecloud, {
    req(genes_df_reactive())
    filtered_genes <- genes_df_reactive()[genes_df_reactive()$Freq >= input$min_freq_genes, ]
    if (nrow(filtered_genes) > input$max_genes) {
      filtered_genes <- head(filtered_genes, input$max_genes)
    }

    output$genewordcloud <- renderPlot({
      if (nrow(filtered_genes) == 0) {
        plot.new()
        text(0.5, 0.5, "No hay genes que cumplan con los criterios.", cex = 1.2)
      } else {
        wordcloud(words = filtered_genes$Gene_symbol,
                  freq = filtered_genes$Freq,
                  min.freq = input$min_freq_genes,
                  max.words = input$max_genes,
                  scale = c(3, 2),
                  random.order = FALSE,
                  rot.per = 0.35,
                  colors = viridis(n = input$max_genes))
      }
    })
  })

  # Selection of a gene in the table
  # show QuickGO and UniProt Links
  observeEvent(input$gene_table_rows_selected, {
    req(genes_df_reactive())
    sel <- input$gene_table_rows_selected
    if (length(sel) > 0) {

      selected_gene(genes_df_reactive()[sel, "Gene_symbol"])

      output$gene_links <- renderUI({
        generate_gene_links(selected_gene())
      })
    } else {
      output$gene_links <- renderUI({ NULL })
    }
  })

  # Filter abstracts containing the selected gene
  filtered_abstracts <- reactive({
    req(selected_gene())
    searchabsL(abstracts_reactive(), include = selected_gene())
  })

  tabla_filtered <- reactive({

    req(filtered_abstracts())

    data.frame(
      PMID        = filtered_abstracts()@PMID,
      Publicación = filtered_abstracts()@Journal,
      stringsAsFactors = FALSE
    )
  })

  output$filtered_articles <- DT::renderDataTable({
    datatable(
      tabla_filtered(),
      selection = list(mode = 'single', selected = 1),
      options = list(
        pageLength = 10,
        language   = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")
      )
    )
  })

  output$gene_info <- renderUI({
    req(selected_gene())

    gene_symbol <- selected_gene()

    quickgo_link <- paste0(
      "<a href='https://www.ebi.ac.uk/QuickGO/searchproducts/",
      gene_symbol,
      "?taxonId=9606' target='_blank'>Link QuickGO</a>"
    )

    uniprot_link <- paste0(
      "<a href='https://www.uniprot.org/uniprotkb?query=",
      gene_symbol,
      "&facets=model_organism%3A9606' target='_blank'>Link UniProt</a>"
    )

    HTML(
      paste0(
        "<p><strong>Artículos Gene_symbol: ", gene_symbol, "</strong></p>",
        "<p>", quickgo_link, " | ", uniprot_link, "</p>"
      )
    )
  })

  output$abstractNavText <- renderUI({
    row_selected <- input$filtered_articles_rows_selected
    if (length(row_selected) == 0) {
      return(NULL)
    }

    selected_pmid <- filtered_abstracts()@PMID[row_selected]

    generate_abstract_html(selected_pmid, abstracts_reactive, doi_reactive)
  })

  # Observer to select a word in the word table
  observeEvent(input$word_table_rows_selected, {
    req(word_tokens_reactive())
    sel <- input$word_table_rows_selected
    if (length(sel) > 0) {
      the_word <- word_tokens_reactive()[sel, "words"]
      selected_word(the_word)
    } else {
      selected_word(NULL)
    }
  })

  filtered_abstracts_by_word <- reactive({
    req(selected_word())
    # searchabsL to filter
    searchabsL(abstracts_reactive(), include = selected_word())
  })

  tabla_filtered_words <- reactive({
    req(filtered_abstracts_by_word())
    data.frame(
      PMID = filtered_abstracts_by_word()@PMID,
      Publicación = filtered_abstracts_by_word()@Journal,
      stringsAsFactors = FALSE
    )
  })

  output$filtered_term_articles <- DT::renderDataTable({
    datatable(
      tabla_filtered_words(),
      selection = list(mode = 'single', selected = 1),
      options = list(
        pageLength = 10,
        language   = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")
      )
    )
  })

  output$word_info <- renderUI({
    req(selected_word())
    HTML(paste0("<p><strong>Artículos que contienen la palabra: ",
                selected_word(), "</strong></p>"))
  })

  output$abstractTermText <- renderUI({
    row_selected <- input$filtered_term_articles_rows_selected
    if (length(row_selected) == 0) {
      return(NULL)
    }

    selected_pmid <- filtered_abstracts_by_word()@PMID[row_selected]

    generate_abstract_html(selected_pmid, abstracts_reactive, doi_reactive)
  })

  # Co-occurrence of terms
  # Generate the dynamic selectInput for Co-occurrence of terms
  output$cooc_term1_ui <- renderUI({
    req(top_words_fixed_reactive())
    selectInput("cooc_term1", "Término 1",
                choices = top_words_fixed_reactive()$words, # top 10 words
                selected = NULL,
                multiple = TRUE)
  })

  output$cooc_term2_ui <- renderUI({
    req(top_genes_reactive())
    selectInput("cooc_term2", "Término 2",
                choices = top_genes_reactive()$Gene_symbol, # top 10 genes
                selected = NULL,
                multiple = TRUE)
  })

  # co-occurrence search logic
  observeEvent(input$coocurrence_search, {
    req(input$cooc_term1, input$cooc_term2, input$cooc_n)

    # Validate that at least one term is selected in each set.
    if(length(input$cooc_term1) == 0 || length(input$cooc_term2) == 0){
      showNotification("Por favor, selecciona al menos un término en ambos conjuntos.", type = "error")
      return(NULL)
    }

    shinyjs::show("loading_spinner") # Show loading spinner

    term1 <- input$cooc_term1
    term2 <- input$cooc_term2
    n <- as.numeric(input$cooc_n)

    # co_occurrence_advance function
    cooc_data <- co_occurrence_advance(abstracts_reactive(), term1, term2, n)

    # Check for results
    if(nrow(cooc_data) == 0){
      output$coocurrence_table <- DT::renderDataTable({
        datatable(
          data.frame(Mensaje = "No se encontraron co-ocurrencias con los criterios seleccionados."),
          options = list(
            paging = FALSE,
            searching = FALSE,
            info = FALSE
          )
        )
      })
      shinyjs::hide("loading_spinner") # Hide loading spinner
      return(NULL)
    }

    output$coocurrence_table <- DT::renderDataTable({
      datatable(
        cooc_data,
        options = list(
          pageLength = 10,
          autoWidth = TRUE,
          language = list(url = "https://cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")
        ),
        rownames = FALSE
      )
    })

    shinyjs::hide("loading_spinner") # Hide loading spinner
  })
}

# Run app
shinyApp(ui = ui, server = server)
