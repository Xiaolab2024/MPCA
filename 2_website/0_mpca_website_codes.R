# Load required libraries
library(shiny)
library(tidyverse)
library(ggpubr)
library(visNetwork)
library(png)
library(dqshiny)
library(DT)
library(gsubfn)
library(shinyjs)
library(glue)
library(shinydashboard)
library(readr)
library(readxl)
library(plotly)
library(igraph)
library(base64enc)
library(htmltools)
library(dplyr)


# Set working directory
#setwd('/Users/shelleywei/Partners HealthCare Dropbox/Shelley Wei/MPA_revision')

# Load data
bat_pearson_data <- read_tsv('bat_pearson.tsv')
liver_pearson_data <- read_tsv('liver_pearson_new.tsv')
liver_pearson_data[1, "uniprot.name"] <- "-"


prot_abun_bat <- read_tsv("bat_prot_abun_updated_jan23.tsv")
met_abun_bat <- read_tsv("bat_met_abunt.tsv", col_types = cols())

bpa <- read_tsv("bat_prot_abunt.tsv")
prot_met_abun_bat <- cbind(bpa, met_abun_bat)
prot_met_abun_bat <- prot_met_abun_bat[ , -10231]
colnames(prot_met_abun_bat) <- toupper(colnames(prot_met_abun_bat))
colnames(prot_met_abun_bat) <- gsub("[^a-zA-Z0-9 ]", "", colnames(prot_met_abun_bat))


prot_abun_liver <- read_tsv("liver_prot_abun_updated_jan23.tsv")
met_abun_liver <- read_tsv("liver_met_abunt.tsv")
colnames(met_abun_liver)[1] <- "analyte"

lpa <- read_tsv("liver_prot_abunt.tsv")
prot_met_abun_liver <- cbind(lpa, met_abun_liver)
prot_met_abun_liver <- prot_met_abun_liver[ , -10231]
colnames(prot_met_abun_liver) <- toupper(colnames(prot_met_abun_liver))
colnames(prot_met_abun_liver) <- gsub("[^a-zA-Z0-9 ]", "", colnames(prot_met_abun_liver))

map_bat <- read_excel("bat_to_chebi_mapping.xlsx")
map_liver <- read_excel("liver_to_chebi_mapping.xlsx")

reactome_bat1 <- read_tsv("reactome_bat.tsv")
reactome_bat <- left_join(reactome_bat1, unique(map_bat[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))
rhea_bat1 <- read_tsv("rhea_bat.tsv")
rhea_bat <- left_join(rhea_bat1, unique(map_bat[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))
tcdb_bat1 <- read_tsv("tcdb_bat.tsv")
tcdb_bat <- left_join(tcdb_bat1, unique(map_bat[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))


#transporter_proteins_bat <- as.data.frame(unique(tcdb_bat$gene.symbol.x))
#transporter_proteins_bat$Transporter <- "transporter"
#colnames(transporter_proteins_bat)[1] <- "gene.symbol"
#transporter_proteins_bat$gene.symbol <- toupper(transporter_proteins_bat$gene.symbol)

rhea_proteins_bat <- as.data.frame(rhea_bat$gene.symbol.x)
colnames(rhea_proteins_bat)[1] <- "gene.symbol"

reactome_proteins_bat <- as.data.frame(reactome_bat$gene.symbol.x)
colnames(reactome_proteins_bat)[1] <- "gene.symbol"

#metabolic_proteins_bat <- as.data.frame(rbind(rhea_proteins_bat, reactome_proteins_bat))
#metabolic_proteins_bat$gene.symbol <- toupper(metabolic_proteins_bat$gene.symbol)
#metabolic_proteins_bat <- unique(metabolic_proteins_bat)
#metabolic_proteins_bat$Metabolic <- "metabolic"

reactome_liver1 <- read_tsv("reactome_liver.tsv")
reactome_liver <- left_join(reactome_liver1, unique(map_liver[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))
rhea_liver1 <- read_tsv("rhea_liver.tsv")
rhea_liver <- left_join(rhea_liver1, unique(map_liver[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))
tcdb_liver1 <- read_tsv("tcdb_liver.tsv")
tcdb_liver <- left_join(tcdb_liver1, unique(map_liver[ , c(1, 3)]), by = c("chebi.identifier" = "ChEBI.identifier"))


#transporter_proteins_liver <- as.data.frame(unique(tcdb_liver$gene.symbol.x))
#transporter_proteins_liver$Transporter <- "transporter"
#colnames(transporter_proteins_liver)[1] <- "gene.symbol"
#transporter_proteins_liver$gene.symbol <- toupper(transporter_proteins_liver$gene.symbol)

rhea_proteins_liver <- as.data.frame(rhea_liver$gene.symbol.x)
colnames(rhea_proteins_liver)[1] <- "gene.symbol"

reactome_proteins_liver <- as.data.frame(reactome_liver$gene.symbol.x)
colnames(reactome_proteins_liver)[1] <- "gene.symbol"

#metabolic_proteins_liver <- as.data.frame(rbind(rhea_proteins_liver, reactome_proteins_liver))
#metabolic_proteins_liver$gene.symbol <- toupper(metabolic_proteins_liver$gene.symbol)
#metabolic_proteins_liver <- unique(metabolic_proteins_liver)
#metabolic_proteins_liver$Metabolic <- "metabolic"

validation_scores <- read_csv("0_validation_scores.csv")
vs_bat <- validation_scores[ , c(2:5)]
vs_liver <- validation_scores[ , c(2, 6:8)]

bat_lasso <- read_tsv("bat_lasso_list.tsv")
liver_lasso <- read_tsv("liver_lasso_list.tsv")

prot_descript_all <- read_tsv("uniprotkb_mouse_AND_model_organism_1009_2024_01_22.tsv")
prot_descript_all$`Gene Names (primary)` <- toupper(prot_descript_all$`Gene Names (primary)`)

#mito_proteins <- read_csv("mouse_mito_proteins.csv")
#mito_proteins <- as.data.frame(mito_proteins$Symbol)
#mito_proteins$Mito <- "mitochondrial"
#colnames(mito_proteins)[1] <- "gene.symbol"
#mito_proteins$gene.symbol <- toupper(mito_proteins$gene.symbol)

lasso_annotations_bat <- read_csv("lasso_annotations_bat.csv")
lasso_annotations_liver <- read_csv("lasso_annotations_liver.csv")

# Replace this path with the correct path to your image file
img_path <- "MPCA_logo_clean.png"
img_base64 <- base64encode(img_path)

# Define UI
ui <- dashboardPage(
  dashboardHeader(
    title = "MPCA"
  ),
  dashboardSidebar(
    width = 280,
    tags$style(HTML(paste0("
      /* Header and Logo Styling */
      .main-header .logo {
        margin: 0;
        padding: 0;
        text-align: center;
        font-weight: bold;
        font-size: 24px; /* Modify font size */
        background: url('data:image/png;base64,", img_base64, "') no-repeat 47px 5px; /* Shift up by reducing the vertical value */
        background-size: 48px 40px; /* Logo size: width 50px, height 40px */
        padding-left: 30px; /* Adjust for logo spacing */
      }
    
      /* Sidebar Menu Styling */
      .sidebar-menu > li > a {
        font-size: 15px;  /* Tab font size */
        padding: 30px;    /* Increase padding for larger tabs */
      }
    
      /* Optional: Sidebar Icon Size */
      .sidebar-menu > li > a > .fa {
        font-size: 25px;  /* Icon size */
      }
    "))),
    sidebarMenu(
      menuItem("Home", tabName = "home"),
      #menuItem("Tour Guide", tabName = "tour"),
      menuItem("Abundance", tabName = "prot_met_abun"),
      menuItem("Correlations", tabName = "correlations"),
      menuItem("Reaction/Pathway Recapitulation", tabName = "recapitulation"),
      menuItem("LASSO Prediction", tabName = "predictions"),
      menuItem("Download", tabName = "download"),
      menuItem("Acknowledgements", tabName = "acknowledgements")
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        /* Set the font family to Helvetica */
        body, h1, h2, h3, h4, h5, h6, p, .box, .sidebar-menu, .navbar, .main-header .logo {
          font-family: 'Helvetica', sans-serif;
        }

        /* Optional: Ensure it falls back to Arial and sans-serif if Helvetica is unavailable */
        body {
          font-family: 'Helvetica', 'Arial', sans-serif;
          font-size: 12px;
        }
      "))
    ),
    tags$head(
      tags$style(HTML('
        body {
          overflow-y: scroll !important;
          height: 100%;')
      )
    ),
    tabItems(
      # tab 0
      tabItem(
        tabName = "home",
        tags$h1("Welcome to MPCA!"),
        #HTML("insert updated abstract."),
        imageOutput("logo"), # LOGO MADE WITH BIORENDER IDK WHAT CITATIONS ARE NECESSARY
        tags$h1("Citation"),
        HTML("<p>If this website is useful to you, please consider citing <a href='blanks'> space holder</a>!</p>"),
        #imageOutput("image_tab1")
      ),
      #tab 0.5
      #tabItem(
      #  tabName = "tour",
      #  tags$h1("Tour Guide of MPCA Application"),
      #  tags$h3("Pearson Correlations"),
      #  tags$h4("Select a protein and metabolite to explore their correlation"),
      #  tags$video(id="video1", width="400px",height="400px",type = "video/mp4",src = "correlation_tab.mp4", controls = "controls"),
      
      #  tags$h3("Database Recapitulations"),
      #  tags$h4("Select a protein, metabolite, or database to explore pathway and reaction networks"),
      #  tags$video(id="video2", width="400px",height="400px",type = "video/mp4",src = "recapitulation_tab.mp4", controls = "controls"),
      
      #  tags$h3("Machine Learning Predictions"),
      #  tags$h4("Select a protein or metabolite to explore positive and negative protein predictors"),
      #  tags$video(id="video3", width="400px",height="400px",type = "video/mp4",src = "predictions_tab.mp4", controls = "controls"),
      
      #  tags$h3("Protein/Metabolite Abundance Search"),
      #  tags$h4("Select a metabolite to explore its abundance"),
      #  tags$video(id="video4", width="400px",height="400px",type = "video/mp4",src = "single_metabolite_tab.mp4", controls = "controls"),
      #),
      # tab 1
      tabItem(
        tabName = "correlations",
        h2("Protein-metabolite Pearson Correlations:"),
        selectInput("tissue_type", "Choose a Tissue:", c("BAT", "Liver"), selected = "BAT", multiple = FALSE),
        tabsetPanel(
          id = "tabs1",
          tabPanel("Edge Query",
                   tags$h4("Search a protein-metabolite edge for correlation statistics"),
                   sidebarLayout(
                     sidebarPanel(
                       fluidRow(
                         column(6, selectizeInput("search_protein", "Search for a protein:", choices = NULL, multiple = FALSE)),
                         column(6, selectInput("search_type", "Search for:", c("Gene Symbol", "Uniprot ID"), selected = "Gene Symbol")),
                         column(6, selectizeInput("search_metabolite", "Search for a metabolite:", choices = NULL, multiple = FALSE))
                       ),
                       actionButton("search_button", "Search")
                     ),
                     mainPanel(
                       conditionalPanel(
                         condition = "input.search_button > 0",
                         tags$h3("Correlation Statistics:"),
                         paste(" "),
                         uiOutput("corr"),
                         paste(" "),
                         tags$h3("Pearson Plot: "),
                         plotOutput("pearson_plot", height = "450px")
                       )
                     )
                   )
                   
          ),
          tabPanel("Node Query",
                   tags$h4("Search a protein or metabolite for top correlators"),
                   sidebarLayout(
                     sidebarPanel(
                       fluidRow(
                         column(6, selectizeInput("search_nodes", "Enter your search query:", choices = NULL, multiple = FALSE)),
                         column(6, selectInput("search_type101", "Search for:", c("Metabolite", "Gene Symbol", "Uniprot ID"), selected = "Metabolite")),
                       ),
                       actionButton("search_button101", "Search")
                     ),
                     mainPanel(
                       conditionalPanel(
                         condition = "input.search_button101 > 0 && (input.search_type101 == 'Gene Symbol' || input.search_type101 == 'Uniprot ID') && input.search_nodes != ''",
                         tags$h3("All Correlates"),
                         HTML(paste("<font color='red'>Positive edges in red.</font>")),
                         HTML(paste("<font color='royalblue'>Negative edges in blue.</font>")),
                         visNetworkOutput("correlations", width = "100%", height = "500px")
                       ),
                       conditionalPanel(
                         condition = "input.search_button101 > 0 && input.search_nodes != ''",
                         tags$h3("Postive correlators (Ranked):"),
                         DTOutput("pos"),
                         tags$h3("Negative correlators (Ranked):"),
                         DTOutput("neg")
                       )
                     )
                   )
                   
          )
        )
        
      ),
      # tab 2
      tabItem(
        tabName = "recapitulation",
        style = "overflow-y: scroll; height: 1100px;",  # Set the desired height,
        h2("Pathway/Reaction Recapitulations"),
        p("Search for significant protein-metabolite correlations recapitulating established biochemical reactions, metabolite-transporter relationships, and pathways from Rhea, TCDB, and Reactome."),
        selectInput("tissue_type2", "Choose a Tissue:", c("BAT", "Liver"), selected = "BAT", multiple = FALSE),
        sidebarPanel(
          fluidRow(
            column(6, selectizeInput("search_recapitulation", "Enter your query:", choices = NULL, multiple = FALSE)),
            column(6, selectInput("search_type2", "Search for:", c("Metabolite", "Protein (Gene Symbol)", "Protein (Uniprot ID)"), selected = "Metabolite"))
          ),
          actionButton("search_button2", "Search")
        ),
        mainPanel(
          conditionalPanel(
            condition = "input.search_button2 > 0",
            tabsetPanel(
              id = "tabs",
              tabPanel("Table",
                       tags$h3("All pathways and reactions"),
                       tags$h4("TCDB"),
                       p("The database details a comprehensive IUBMB approved classification system for membrane transport proteins known as the Transporter Classification (TC) system. The TC system is analogous to the Enzyme Commission (EC) system for classification of enzymes, except that it incorporates both functional and phylogenetic information. Curated annotations, TC numbers, and external references for 1922 families of transport proteins are provided. Transport systems are classified on the basis of five criteria, and each of these criteria corresponds to one of the five numbers or letters within the TC accession for a particular type of transporter."),
                       tags$h4("Reactome"),
                       p("REACTOME is an open-source, open access, manually curated and peer-reviewed pathway database. Our goal is to provide intuitive bioinformatics tools for the visualization, interpretation and analysis of pathway knowledge to support basic and clinical research, genome analysis, modeling, systems biology and education. Founded in 2003, the Reactome project is led by Lincoln Stein of OICR, Peter D’Eustachio of NYU Langone Health, Henning Hermjakob of EMBL-EBI, and Guanming Wu of OHSU."),
                       tags$h4("Rhea"),
                       p("Rhea is an expert-curated knowledgebase of chemical and transport reactions of biological interest - and the standard for enzyme and transporter annotation in UniProtKB. Rhea uses the chemical dictionary ChEBI (Chemical Entities of Biological Interest) to describe reaction participants."),
                       p("note: rows with empty statistical values indicates it was not recapitulated in the database."),
                       p("Swipe for all columns"),
                       DTOutput("recapitulation_dt"),
                       
              ),
              tabPanel("Network",
                       tags$h3("Pathway/Reaction Interactive"),
                       conditionalPanel(
                         condition = "['TCDB ID', 'Reactome ID', 'Rhea ID'].indexOf(input.search_type2) !== -1",
                         HTML(paste("Full reaction network:")),
                         HTML(paste("<font color='purple'>Metabolite nodes in purple.</font>")),
                         HTML(paste("<font color='green'>Protein nodes in green.</font>")),
                         HTML(paste("<font color='red'>Positive recapitulated edges in red.</font>")),
                         HTML(paste("<font color='blue'>Negative recapitulated edges in blue.</font>")),
                         HTML(paste("Hover over protein nodes for more info.")),
                         visNetworkOutput("interactive", width = "100%", height = "500px")
                       ),
                       conditionalPanel(
                         condition = "['Metabolite', 'Protein (Gene Symbol)', 'Protein (Uniprot ID)'].indexOf(input.search_type2) !== -1",
                         paste("Select a row in the table to see reaction/pathway."),
                         conditionalPanel(
                           HTML(paste("Full reaction network:")),
                           HTML(paste("<font color='purple'>Metabolite nodes in purple.</font>")),
                           HTML(paste("<font color='green'>Protein nodes in green.</font>")),
                           HTML(paste("<font color='red'>Positive recapitulated edges in red.</font>")),
                           HTML(paste("<font color='blue'>Negative recapitulated edges in blue.</font>")),
                           HTML(paste("Hover over protein nodes for more info.")),
                           condition = "input.recapitulation_dt_rows_selected > 0",
                           visNetworkOutput("selected_row", width = "100%", height = "500px")
                         )
                       )
              )
            )
          )
        )
      ),
      # tab 3
      tabItem(
        tabName = "predictions",
        h2("Protein predictors of metabolite abundance"),
        selectInput("tissue_type3", "Choose a Tissue:", c("BAT", "Liver"), selected = "BAT", multiple = FALSE),
        tabsetPanel(
          id = "lasso",
          tabPanel(
            "Prediction",
            tags$h4("Search a protein or metabolite:"),
            sidebarPanel(
              fluidRow(
                column(6, selectizeInput("search_predictions", "Enter your query:", choices = NULL, multiple = FALSE)),
                column(6, selectInput("search_type3", "Search for:", c("Metabolite", "Protein (Gene Symbol)", "Protein (Uniprot ID)"), selected = "Metabolite"))
              ),
              actionButton("search_button3", "Search")
            ),
            mainPanel(
              #style = "overflow-y: scroll; overflow-x: scroll; height: 900px; width: 1200px; position: relative; margin-top: 20px;",  # Set the desired height,
              conditionalPanel(
                condition = "input.search_button3 > 0",
                conditionalPanel(
                  condition = "input.search_type3 == 'Metabolite'",
                  HTML(paste("<font color='#E75480'>Mitochondrial proteins are colored in pink.</font>")),
                  HTML(paste("<font color='royalblue'>Transporter proteins are colored in blue.</font>")),
                  HTML(paste("<font color='#FFCC00'>Metabolic proteins are colored in yellow.<br></font>")),
                  HTML(paste("<font color='red'>Extreme outliers are colored in red</font>")),
                  HTML(paste("<font color='grey'>Unspecified proteins are colored in grey.<br></font>")),
                  HTML(paste("<font color='grey'><br></font>")),
                  HTML(paste("<font color='black'>\nProteins that are classified more than once are colored as the overlap between the categories.</font>")),
                  HTML(paste("<font color='grey'>    <br></font>")),            ), 
                plotlyOutput("ml_plot")
              )
            )
          ),
          tabPanel(
            "Validation Score",
            tags$h4("Search for a metabolite:"),
            sidebarPanel(
              fluidRow(
                column(6, selectizeInput("search_validations", "Enter your query:", choices = NULL, multiple = FALSE)),
              ),
              actionButton("search_button4", "Search")
            ),
            mainPanel(
              #style = "overflow-y: scroll; overflow-x: scroll; height: 900px; width: 1200px; position: relative; margin-top: 20px;",  # Set the desired height,
              conditionalPanel(
                condition = "input.search_button4 > 0",
                p("Zoom out to see all metabolites."),
                plotlyOutput("vs_plot")
              )
            )
          )
        )
      ),
      # tab 3.5
      tabItem(
        tabName = "prot_met_abun",
        h2("Single Analyte Search"),
        selectInput("tissue_type35", "Choose a Tissue:", c("BAT", "Liver"), selected = "BAT", multiple = FALSE),
        sidebarPanel(
          fluidRow(
            column(6, selectizeInput("search_abun", "Enter your query:", choices = NULL, multiple = FALSE)),
            column(6, selectInput("search_type35", "Search for:", c("Metabolite", "Protein (Uniprot ID)", "Protein (Gene Symbol)"), selected = "Metabolite"))
          ),
          actionButton("search_button35", "Search")
        ),
        mainPanel(
          conditionalPanel(
            condition = "input.search_button35 > 0",
            tabsetPanel(
              id = "abundance",
              tabPanel("Abundance",
                       h3("Protein/Metabolite Abundance"),
                       textOutput("abun_text"),
                       plotlyOutput("abun_plot")
              ),
              tabPanel("CV",   # Corrected here
                       h3("Coefficient of Variation Density Plot"),
                       textOutput("dist_text"),
                       plotOutput("dist_plot")
                       # Add additional UI elements for the Distribution tab if needed
              )
            )
          )
        )
      ),
      # tab 4
      tabItem(
        tabName = "download",
        h2("Bulk download data files"),
        p("Pearson Correlation Files:"),
        downloadButton("download_btn_1", "BAT (csv)"),
        downloadButton("download_btn_2", "Liver (csv)"),
        p("Database Recapitulation Files:"),
        downloadButton("download_btn_4", "Rhea BAT (csv)"),
        downloadButton("download_btn_5", "Reactome BAT (csv)"),
        downloadButton("download_btn_6", "TCDB BAT (csv)"),
        downloadButton("download_btn_7", "Rhea Liver (csv)"),
        downloadButton("download_btn_8", "Reactome Liver (csv)"),
        downloadButton("download_btn_9", "TCDB Liver (csv)"),
        p("Machine Learning Metabolite Lists:"),
        downloadButton("download_btn_13", "BAT (csv)"),
        downloadButton("download_btn_14", "Liver (csv)"),
      ),
      # tab 5
      tabItem(
        tabName = "acknowledgements",
        h2("Acknowledgements"),
        p("this is the content of the acknowledgements tab")
      )
    )
  ),
  skin = "purple"
)

server <- function(input, output, session) {
  
  #tab0
  output$logo <- renderImage({
    list(src="MPCA_logo_clean.png", width=300, height=270)
  },
  deleteFile = FALSE)
  #output$image_tab1 <- renderImage({
  #  list(src="data/final/OPABAT.png",width=240,height=180)
  #},
  #deleteFile = FALSE)
  
  #tab1
  filtered_pearson <- reactiveVal(NULL)
  filtered_pearson2 <- reactiveVal(NULL)
  sp <- reactiveVal(NULL)
  sm <- reactiveVal(NULL)
  gn <- reactiveVal(NULL)
  gs <- reactiveVal()
  joined1 <- reactiveVal(NULL)
  
  pearson_data <- reactiveVal()
  prot_met_abun <- reactiveVal()
  met_abun <- reactiveVal()
  prot_abun <- reactiveVal()
  abun_all <- reactiveVal()
  abun_all2 <- reactiveVal()
  
  observeEvent(input$tissue_type, {
    if (input$tissue_type == "Liver") {
      pearson_data(liver_pearson_data)
      prot_met_abun(prot_met_abun_liver)
      met_abun(met_abun_liver)
      prot_abun(prot_abun_liver)
    } else {
      pearson_data(bat_pearson_data)
      prot_met_abun(prot_met_abun_bat)
      met_abun(met_abun_bat)
      prot_abun(prot_abun_bat)
      
    }
    
    #gene_entry <- prot_descript_all[ , c(1, 7)]
    #prots <- as.data.frame(colnames(prot_abun()))
    
    prots <- data.frame(
      uniprot.id = colnames(prot_abun()),  # Column names
      gene.symbol = as.character(prot_abun()[1, ])  # First row values
    )
    prots <- as.data.frame(prots[-1, ])
    #colnames(prots)[1] <- "Entry"
    #abun_all <- left_join(prots, gene_entry, by = "Entry")
    #abun_all <- abun_all[complete.cases(abun_all), ]
    new <- data.frame(uniprot.id = "-", gene.symbol = "-")
    #colnames(new)[2] <- "Gene Names (primary)"
    #abun_all <- rbind(new, abun_all)
    abun_all <- rbind(new, prots)
    abun_all <- unique(abun_all)
    abun_all(abun_all)
  })
  
  
  
  # Process protein searches
  autocomplete_choices_protein <- reactive({
    column_to_search <- switch(input$search_type,
                               "Uniprot ID" = "uniprot.id",
                               "Gene Symbol" = "gene.symbol"
    )
    unique(abun_all()[[column_to_search]])
  })
  
  observe({
    updateSelectizeInput(session, "search_protein", choices = autocomplete_choices_protein(), server = TRUE)
  })
  
  # Process metabolite searches
  autocomplete_choices_metabolite <- reactive({
    list <- as.data.frame(unique(colnames(met_abun())))
    list[1, 1] <- "-"
    colnames(list) <- " "
    return(list)
  })
  
  observe({
    updateSelectizeInput(session, "search_metabolite", choices = autocomplete_choices_metabolite())
  })
  
  
  observeEvent(input$search_button, {
    search_protein <- toupper(input$search_protein)
    search_metabolite <- toupper(input$search_metabolite)
    
    sm(search_metabolite)
    
    search_type <- input$search_type
    
    if (search_type == "Uniprot ID") {
      column_to_search <- "uniprot.id"
    } else if (search_type == "Gene Symbol") {
      column_to_search <- "gene.symbol"
      joined1(left_join(pearson_data(), abun_all(), by = "uniprot.id") %>% mutate(!!names(.)[3] := .[[9]]) %>% select(-9) %>% rename("gene.symbol" = colnames(.)[3]))
    } else {
      column_to_search <- NULL
    }
    
    if (!is.null(column_to_search)) {
      if (search_type != "Gene Symbol") {
        filtered_pearson(pearson_data()[toupper(pearson_data()[[column_to_search]]) %in% search_protein, ])
      } else if (search_type == "Gene Symbol") {
        filtered_pearson(joined1()[toupper(joined1()[[column_to_search]]) %in% search_protein, ])
      }
      if (nrow(filtered_pearson()) == 0) {
        gn(search_metabolite)
        if(column_to_search == "uniprot.id"){
          sp(search_protein)
          gs(abun_all()[abun_all()[, 1] == search_protein, 2])
        } else {
          sp(abun_all()[abun_all()[, 2] == search_protein, 1])
          gs(search_protein)
        }
      } else {
        sp(filtered_pearson()[2, 1])
        gn(filtered_pearson()[1, 4])
        if(!is.na(filtered_pearson()[1, 3])){
          gs(filtered_pearson()[1, 3])
        } else {
          gs(filtered_pearson()[1, 1])
        }
      }
    }
    
    filtered_pearson2(filtered_pearson()[toupper(filtered_pearson()[["column"]]) == search_metabolite, ])
    
    
    # output
    output$corr <- renderUI({
      corr_stats <- as.data.frame(filtered_pearson2()[1, ])
      if (!is.na(corr_stats[1, 5])) {
        HTML(paste(
          "Correlation Coefficient: ", corr_stats$cor, "<br>",
          "p Value: ", corr_stats$p, "<br>",
          "p adj Value: ", corr_stats$p_adj, "<br>",
          "n Number: ", corr_stats$n
        ))
      }else{
        paste("Edge not included in MPCA network because n<50")
      }
    })
    corr_stats <- as.data.frame(filtered_pearson2()[1, ])
    p1 <- as.data.frame(sp())
    col <- as.data.frame(gn())
    gene <- as.data.frame(gs())
    p2 <- as.character(p1[1, 1])
    p <- as.character(gsub("[^a-zA-Z0-9 ]", "", p2))
    m1 <- sm()
    m <- gsub("[^a-zA-Z0-9 ]", "", m1)
    m2 <- as.character(sub("^([0-9])", " \\1", m))
    prot_met_abun <- prot_met_abun()
    colnames(prot_met_abun) <- gsub("-.*", "", colnames(prot_met_abun))
    colnames(prot_met_abun) <- sub("^([0-9])", " \\1", colnames(prot_met_abun))
    
    if(p %in% colnames(prot_met_abun) && m2 %in% colnames(prot_met_abun)) {
      output$pearson_plot <- renderPlot({
        ggscatter(prot_met_abun, x = p, y = m2,
                  add = "reg.line", conf.int = TRUE, 
                  conf.int.level = 0.95,
                  cor.coef = TRUE, cor.method = "pearson",
                  color = "black", 
                  alpha=0.65,
                  add.params = list(color = "darkred",
                                    fill = "gray"),
                  xlab = paste0("rel. abundance of ",gene), ylab = paste0("rel. abundance of ",m1))+
          theme(
            axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10)
          )
      }, res=96)
    } else {
      output$pearson_plot <- renderPlot({
        ggscatter(data.frame(x = numeric(0), y = numeric(0)),
                  x = "x", y = "y") +
          ggtitle("Invalid search query") +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank()
          )
      })
    }
  })
  
  
  #tab 1 pt 2
  filtered_data <- reactiveVal(NULL)
  joined <- reactiveVal(NULL)
  
  autocomplete_choices_nodes <- reactive({
    # Choose the column or file based on input$search_type101
    search_col <- switch(input$search_type101,
                         "Uniprot ID" = "uniprot.id",  # Reference 'uniprot.id' from abun_all2
                         "Gene Symbol" = "gene.symbol",  # Reference 'gene.symbol' from abun_all2
                         "Metabolite" = "column"  # Reference 'metabolite' from a different file (e.g., metabolite_data)
    )
    if (input$search_type101 == "Metabolite") {
      unique(pearson_data()[[search_col]])  
    } else {
      unique(abun_all()[[search_col]])
    }
  })
  
  observe({
    updateSelectizeInput(session, "search_nodes", choices = autocomplete_choices_nodes())
  })
  
  
  observeEvent(input$search_button101, {
    search_nodes <- toupper(input$search_nodes)
    search_type101 <- input$search_type101
    
    if (search_type101 == "Uniprot ID") {
      col_to_search <- "uniprot.id"
    } else if (search_type101 == "Gene Symbol") {
      col_to_search <- "gene.symbol"
      joined(left_join(pearson_data(), abun_all(), by = "uniprot.id") %>% mutate(!!names(.)[3] := .[[9]]) %>% select(-9) %>% rename("gene.symbol" = colnames(.)[3]))
    } else if (search_type101 == "Metabolite") {
      col_to_search <- "column"
    } else {
      col_to_search <- NULL
    }
    
    
    if (!is.null(col_to_search) && col_to_search != "gene.symbol") {
      filtered_data(pearson_data()[toupper(pearson_data()[[col_to_search]]) %in% search_nodes, ])
    } else if(!is.null(col_to_search) && col_to_search == "gene.symbol") {
      filtered_data(joined()[toupper(joined()[[col_to_search]]) %in% search_nodes, ])
    }

    #des <- prot_descript_all[ , c(7, 4)]
    #colnames(des) <- c("gene.symbol", "description")
    #des <- des[!duplicated(des$gene.symbol), ]
    
    if (input$search_type101 != "Metabolite"){
      output$correlations <- renderVisNetwork({
        edges_pearson <- filtered_data()
        edges_pearson <- edges_pearson[abs(as.numeric(edges_pearson$p_adj)) < 0.05 | is.na(edges_pearson$p_adj), ]
        edges_pearson <- edges_pearson[complete.cases(edges_pearson$gene.symbol, edges_pearson$column), ]
        
        edges_pearson$stats <- paste("cor:", edges_pearson$cor, "| n:", edges_pearson$n, "| p adj:", edges_pearson$p_adj)
        
        #metabolite_names <- unique(edges_pearson[, c("column")])
        #metabolite_names$description <- metabolite_names$column
        #colnames(metabolite_names)[1] <- "gene.symbol"
        #des_met <- rbind(des, metabolite_names)
        
        edges <- data.frame(from = edges_pearson$gene.symbol, to = edges_pearson$column, title = edges_pearson$stats, cor = edges_pearson$cor)
        nodes <- unique(c(edges$from, edges$to))
        nodes <- as.data.frame(nodes) %>% dplyr::rename("label" = "nodes")
        #nodes <- left_join(nodes, des_met, by = c("label" = "gene.symbol"))
        nodes$id <- nodes$label
        nodes$title <- nodes$label
        #nodes<-nodes%>%dplyr::rename("title"="description")
        nodes<-nodes%>%mutate(color=case_when(
          label%in%unique(edges$from)~"darkgreen",
          TRUE~"#b19cd7"))
        
        edges$color <- ifelse(!is.na(edges$cor), ifelse(edges$cor >= 0, "red", "blue"), "grey")
        edges <- edges %>% select(-cor)
        
        nodes <- nodes %>% mutate(
          title = title
        )
        
        visNetwork(nodes, edges) %>%
          visIgraphLayout() %>%
          visNodes(
            size = nodes$value, 
            label = nodes$label,
            title = nodes$title
          ) %>%
          visEdges(color = list(color = edges$color), width = 5) %>%
          visOptions(highlightNearest = list(enabled = T, hover = T), 
                     nodesIdSelection = T)
      })
    }

    
    if(input$search_type101 %in% "Metabolite"){
      prot_descript_all2 <- prot_descript_all[ , c(1, 4)]
      filtered_data_descript <- left_join(filtered_data(), prot_descript_all2, by = c("uniprot.id" = "Entry"))
      filtered_negative <- filtered_data_descript[!is.na(filtered_data_descript$cor) & filtered_data_descript$cor < 0, c(3, 1, 2, 9, 5, 6, 7, 8)]
      filtered_negative <- filtered_negative[abs(as.numeric(filtered_negative$p_adj)) < 0.05, ]
      filtered_negative$cor <- as.numeric(filtered_negative$cor)
      filtered_negative <- filtered_negative[order(filtered_negative$cor), ]
      top_neg <- filtered_negative
      top_neg$rank <- seq_len(nrow(top_neg))
      top_neg <- top_neg[ , c(9, 1:8) ]
      colnames(top_neg) <- c("Rank", "Gene Symbol", "Uniprot ID", "Uniprot Name", "Description", "Correlation Coefficient", "p", "n", "p_adj")
      
      filtered_positive <- filtered_data_descript[!is.na(filtered_data_descript$cor) & filtered_data_descript$cor > 0, c(3, 1, 2, 9, 5, 6, 7, 8)]
      filtered_positive <- filtered_positive[abs(as.numeric(filtered_positive$p_adj)) < 0.05, ]
      filtered_positive$cor <- as.numeric(filtered_positive$cor)
      filtered_positive <- filtered_positive[order(-filtered_positive$cor), ]
      top_pos <- filtered_positive
      top_pos$rank <- seq_len(nrow(top_pos))
      top_pos <- top_pos[ , c(9, 1:8) ]
      colnames(top_pos) <- c("Rank", "Gene Symbol", "Uniprot ID", "Uniprot Name", "Description", "Correlation Coefficient", "p", "n", "p_adj")
      
    } else if(input$search_type101 %in% c("Uniprot ID", "Gene Symbol")){
      filtered_negative <- filtered_data()[!is.na(filtered_data()$cor) & filtered_data()$cor < 0, c(4:8)]
      filtered_negative <- filtered_negative[abs(as.numeric(filtered_negative$p_adj)) < 0.05, ]
      filtered_negative$cor <- as.numeric(filtered_negative$cor)
      filtered_negative <- filtered_negative[order(filtered_negative$cor), ]
      top_neg <- filtered_negative
      top_neg$rank <- seq_len(nrow(top_neg))
      top_neg <- top_neg[ , c(6, 1:5) ]
      colnames(top_neg) <- c("Rank", "Metabolite", "Correlation Coefficient", "p", "n", "p_adj")
      
      filtered_positive <- filtered_data()[!is.na(filtered_data()$cor) & filtered_data()$cor > 0, c(4:8)]
      filtered_positive <- filtered_positive[abs(as.numeric(filtered_positive$p_adj)) < 0.05, ]
      filtered_positive$cor <- as.numeric(filtered_positive$cor)
      filtered_positive <- filtered_positive[order(-filtered_positive$cor), ]
      top_pos <- filtered_positive
      top_pos$rank <- seq_len(nrow(top_pos))
      top_pos <- top_pos[ , c(6, 1:5) ]
      colnames(top_pos) <- c("Rank", "Metabolite", "Correlation Coefficient", "p", "n", "p_adj")
    }
    
    output$pos <- renderDT({
      datatable(
        top_pos, 
        options = list(
          pageLength = 10,            # Default number of entries per page
          lengthMenu = c(10, 20, 50, 100), # Options for entries per page
          dom = 'lfrtip',            # 'l' for length changing input, 'f' for filtering input (search bar), 'r' for processing display element, 't' for table, 'i' for table information summary, 'p' for pagination
          autoWidth = TRUE, 
          scrollX = TRUE
        ), 
        rownames = FALSE
      )
    })
    
    output$neg <- renderDT({
      datatable(
        top_neg, 
        options = list(
          pageLength = 10,            # Default number of entries per page
          lengthMenu = c(10, 20, 50, 100), # Options for entries per page
          dom = 'lfrtip',            # 'l' for length changing input, 'f' for filtering input (search bar), 'r' for processing display element, 't' for table, 'i' for table information summary, 'p' for pagination
          autoWidth = TRUE, 
          scrollX = TRUE
        ), 
        rownames = FALSE
      )
    })
    
  })
  #tab2
  reactome <- reactiveVal()
  rhea <- reactiveVal()
  tcdb <- reactiveVal()
  map <- reactiveVal()
  pd <- reactiveVal()
  
  c2s <- reactiveVal(NULL)
  db <- reactiveVal(NULL)
  data_table <- reactiveVal(NULL)
  
  observeEvent(input$tissue_type2, {
    if (input$tissue_type2 == "Liver") {
      reactome(reactome_liver)
      rhea(rhea_liver)
      tcdb(tcdb_liver)
      pd(liver_pearson_data)
      map(map_liver)
    } else {
      reactome(reactome_bat)
      rhea(rhea_bat)
      tcdb(tcdb_bat)
      pd(bat_pearson_data)
      map(map_bat)
    }
  })
  
  # Reactive expression to process searches
  autocomplete_choices_recapitulation <- reactive({
    pearson_data_value <- pd()
    recapitulated_rhea <- rhea()
    sorted <- rhea_bat %>% arrange(desc(row))
    sorted <- sorted[c(which(sorted[[1]] == "CDO1"), which(sorted[[1]] != "CDO1")), ]    
    unique_values <- unique(sorted[[1]])
    pearson_data_value <- pearson_data_value[order(match(pearson_data_value[[3]], unique_values, incomparables = NA)), ]
    pearson_data_value <- pearson_data_value[c(which(pearson_data_value[[1]] == "-"), which(pearson_data_value[[1]] != "-")), ]
    
    if (input$search_type2 == "Protein (Gene Symbol)") {
      c2s("gene.symbol")
      unique(pearson_data_value$gene.symbol)
    } else if (input$search_type2 == "Protein (Uniprot ID)") {
      c2s("uniprot.id")
      unique(pearson_data_value$uniprot.id)
      #} else if (input$search_type2 == "Protein (Entry Name)") {
      #  c2s("uniprot.name")
      #  unique(pearson_data_value$uniprot.name)
    } else if (input$search_type2 == "Metabolite") {
      c2s("column")
      unique(pearson_data_value$column)
      #} else if (input$search_type2 == "TCDB ID") {
      #  db(tcdb())
      #  unique(tcdb()$ID)
      #} else if (input$search_type2 == "Reactome ID") {
      #  db(reactome())
      #  unique(reactome()$ID)
      #} else if (input$search_type2 == "Rhea ID") {
      #  db(rhea())
      #  unique(rhea()$ID)
    }
  })
  
  # Update select input choices using an observer
  observe({
    updateSelectizeInput(session, "search_recapitulation", choices = autocomplete_choices_recapitulation(), server = TRUE)
  })
  
  observeEvent(input$search_button2, {
    search_recapitulation <- input$search_recapitulation
    des <- prot_descript_all[ , c(7, 4)]
    colnames(des) <- c("gene.symbol", "description")
    des <- des[!duplicated(des$gene.symbol), ]
    
    #if(input$search_type2 %in% c("TCDB ID", "Reactome ID", "Rhea ID")){
    #  data <- db()[db()$ID %in% search_recapitulation, c("ID", "gene.symbol.x", "chebi.identifier", "chebi.name", "edge_db", "cor", "p", "n", "p_adj")]
    #  
    #  data2 <- data
    #  colnames(data2) <- c("ID", "Gene Symbol", "Chebi Identifier", "Chebi Name", "Edge", "Correlation Coefficient", "p", "n", "p_adj")
    #  
    #  output$recapitulation_dt <- renderDT({
    #    if (!is.null(data2)) {
    #      datatable(
    #        data2,
    #        options = list(
    #          pageLength = 10,
    #          dom = 'Bfrtip',  # Adjusts the position of the search box
    #          autoWidth = TRUE,
    #          scrollX = TRUE,
    #          search = list(search = "", regex = TRUE)  # Enables the search box
    #        ),
    #        rownames = FALSE,
    #        selection = "single"  # Optional: allows for row selection
    #      )
    #    }
    #  })
    
    
    #  output$interactive <- renderVisNetwork({
    #    if (input$search_type2 == "TCDB ID") {
    #      edges_all <- tcdb()
    #    } else if (input$search_type2 == "Reactome ID") {
    #      edges_all <- reactome()
    #    } else if (input$search_type2 == "Rhea ID") {
    #      edges_all <- rhea()
    #    }
    #    
    #    edges_all$stats <- paste("cor:", edges_all$cor, "| n:", edges_all$n, "| p adj:", edges_all$p_adj)
    #    
    #    edges_reaction <- edges_all[edges_all$ID %in% search_recapitulation, ]
    #    metabolite_names <- unique(edges_reaction[, c("chebi.name", "chebi.identifier")])
    #    colnames(metabolite_names)[1] <- "gene.symbol"
    #    colnames(metabolite_names)[2] <- "description"
    #    des_met <- rbind(des, metabolite_names)
    #    edges <- data.frame(from = edges_reaction$gene.symbol.x, to = edges_reaction$chebi.name, row = edges_reaction$row, title = edges_reaction$stats)
    #    nodes <- unique(c(edges$from, edges$to))
    #    nodes <- as.data.frame(nodes) %>% dplyr::rename("label" = "nodes")
    #    nodes <- left_join(nodes, des_met, by = c("label" = "gene.symbol"))
    #    nodes$id <- nodes$label
    #    nodes<-nodes%>%dplyr::rename("title"="description")
    #    nodes<-nodes%>%mutate(color=case_when(
    #      label%in%unique(edges$from)~"darkgreen",
    #      TRUE~"#b19cd7"))
    #    
    #    edges$color <- ifelse(!is.na(edges_reaction$cor), ifelse(edges_reaction$cor >= 0, "red", "blue"), "grey")
    #    
    #    nodes <- nodes %>% mutate(
    #      title = title
    #    )
    #    
    #    visNetwork(nodes, edges) %>%
    #      visIgraphLayout() %>%
    #      visNodes(
    #        size = nodes$value, 
    #        label = nodes$label,
    #        title = nodes$title
    #      ) %>%
    #      visEdges(color = list(color = edges$color), width = 5) %>%
    #      visOptions(highlightNearest = list(enabled = T, hover = T), 
    #                 nodesIdSelection = T)
    #  })
    #  
    #}else{
    c2s2 <- c2s()
    if(c2s2 == "column") {
      data <- map()[map()$Metabolite %in% search_recapitulation, c(1,3)]
      data_table(rbind(tcdb()[tcdb()$chebi.identifier %in% data$ChEBI.identifier, c(1, 2, 4, 18, 5, 26, 6, 12, 13, 14, 15, 25)], rhea()[rhea()$chebi.identifier %in% data$ChEBI.identifier, c(2, 1, 4, 17, 3, 25, 5, 11, 12, 13, 14, 24)], reactome()[reactome()$chebi.identifier %in% data$ChEBI.identifier, c(1:3, 17, 4, 25, 5, 11, 12, 13, 14, 24)]))
      data_table(data_table()[!is.na(data_table()$cor), ])
      data_table(data_table() %>%
                   mutate(
                     cor = formatC(cor, format = "e", digits = 2),
                     p = formatC(p, format = "e", digits = 2),
                     p_adj = formatC(p_adj, format = "e", digits = 2)))
      data_table(data_table()%>%
                   mutate(Metabolite = ifelse(is.na(Metabolite), chebi.name, Metabolite)) %>%
                   select(-chebi.name))
      #data_table()$Metabolite <- search_recapitulation
    }else{
      data <- pd()[pd()[[c2s2]] %in% search_recapitulation, ]
      #gs <- data[1, 4] #this is only for proteins, need to figure out what to do about metabolites
      data_table(rbind(tcdb()[tcdb()$gene.symbol.x %in% data$gene.symbol, c(1, 2, 4, 18, 5, 26, 6, 12, 13, 14, 15, 25)], rhea()[rhea()$gene.symbol.x %in% data$gene.symbol, c(2, 1, 4, 17, 3, 25, 5, 11, 12, 13, 14, 24)], reactome()[reactome()$gene.symbol.x %in% data$gene.symbol, c(1:3, 17, 4, 25, 5, 11, 12, 13, 14, 24)]))
      data_table(data_table()[!is.na(data_table()$cor), ])
      data_table(data_table() %>%
                   mutate(
                     cor = formatC(cor, format = "e", digits = 2),
                     p = formatC(p, format = "e", digits = 2),
                     p_adj = formatC(p_adj, format = "e", digits = 2)))
      data_table(data_table()%>%
                   mutate(Metabolite = ifelse(is.na(Metabolite), chebi.name, Metabolite)) %>%
                   select(-chebi.name))
      #data_table()$Metabolite <- search_recapitulation
    }
    
    output$recapitulation_dt <- renderDT({
      data_dt <- data_table()
      colnames(data_dt) <- c("ID", "Gene Symbol", "Chebi Identifier (Edge)", "Chebi Identifier (Lowest)", "Metabolite Name", "Edge", "Correlation Coefficient", "p", "n", "p_adj", "Database")
      
      if (!is.null(data_dt)) {
        if(nrow(data_dt) > 0) {
          datatable(
            data_dt,
            options = list(
              pageLength = 10,
              dom = 'Bfrtip',  # Adjusts the position of the search box
              autoWidth = TRUE,
              scrollX = TRUE,
              search = list(search = "", regex = TRUE),  # Enables the search box
              initComplete = JS("function(settings, json) { 
                           this.api().columns().every(function() {
                             var column = this;
                             var input = document.createElement('input');
                             $(input).appendTo($(column.footer()).empty())
                             .on('change', function() {
                               column.search($(this).val(), false, false, true).draw();
                             });
                           });
                         }")
            ),
            rownames = FALSE,
            selection = "single"
          )
        } else {
          no_data_df <- data.frame(. = "No MVPA edge involving the search entry was recapitulated.")
          
          # Render the dataframe with the message in a datatable
          datatable(
            no_data_df,
            options = list(
              searching = FALSE, # Disable searching because it's a single message
              pageLength = 1,     # Set page length to 1 since it's a single row
              dom = 't'           # Only display the table (no other controls)
            ),
            rownames = FALSE,
            selection = "none"
          )        }
      }
    })
    
    
    observeEvent(input$recapitulation_dt_rows_selected, {
      selected <- input$recapitulation_dt_rows_selected
      selected_row <- data_table()[selected, ]
      if (!all(is.na(selected_row))){
        output$selected_row <- renderVisNetwork({
          if(selected_row$database == "Rhea"){
            edges_all <- rhea()
          }else if(selected_row$database == "Reactome"){
            edges_all <- reactome()
          }else if(selected_row$database == "TCDB"){
            edges_all <- tcdb()
          }else{
            edges_all <- NULL
          }
          
          
          edges_all <- edges_all %>%
            mutate(
              cor = formatC(cor, format = "e", digits = 2),
              p = formatC(p, format = "e", digits = 2),
              p_adj = formatC(p_adj, format = "e", digits = 2))
          edges_all$stats <- paste("cor:", edges_all$cor, "|n:", edges_all$n, "|p adj:", edges_all$p_adj)
          
          edges_reaction <- edges_all[edges_all$ID %in% selected_row$ID, ]
          
          metabolite_names <- unique(edges_reaction[, c("Metabolite", "chebi.identifier")])
          colnames(metabolite_names)[1] <- "gene.symbol"
          colnames(metabolite_names)[2] <- "description"
          des_met <- rbind(des, metabolite_names)
          des_met <- des_met %>% distinct(gene.symbol, .keep_all = TRUE)
          edges <- data.frame(
            from = edges_reaction$gene.symbol.x,
            to = ifelse(edges_reaction$Metabolite == "" | is.na(edges_reaction$Metabolite), 
                        edges_reaction$chebi.name, 
                        edges_reaction$Metabolite),
            row = edges_reaction$row,
            title = edges_reaction$stats,
            cor = edges_reaction$cor
          )
          # Check the resulting data frame
          head(edges)
          nodes<-unique(c(edges$from,edges$to))
          nodes<-as.data.frame(nodes)%>%dplyr::rename("label"="nodes")
          nodes <- left_join(nodes, des_met, by = c("label" = "gene.symbol"))
          #nodes$title <- nodes$label
          nodes$id<-nodes$label
          #nodes<-nodes%>%dplyr::rename("title"="description")
          #nodes$title <- paste(nodes[[1]], nodes[[2]], sep = ": ")
          nodes$title <- ifelse(is.na(nodes[[2]]) | nodes[[2]] == "", nodes[[1]], paste(nodes[[1]], nodes[[2]], sep = ": "))
          
          nodes<-nodes%>%mutate(color=case_when(
            label%in%unique(edges$from)~"darkgreen",
            TRUE~"#b19cd7"))
          
          edges$color <- ifelse(edges$cor != " NA", ifelse(edges$cor >= 0, "red", "blue"), "grey")
          
          edges <- edges %>% select(-cor)
          
          #edges$group <- ifelse(!is.na(edges_reaction$cor), ifelse(edges_reaction$cor >= 0, "Positive Correlation", "Negative Correlation"), "No Correlation")
          
          visNetwork(nodes, edges) %>%
            visIgraphLayout(layout = "layout_in_circle") %>%
            visNodes(
              size = nodes$value, 
              label = nodes$label,
              title = nodes$title
            ) %>%
            visEdges(color = list(color = edges$color), width = 5) %>%
            visOptions(highlightNearest = list(enabled = T, hover = T), 
                       nodesIdSelection = T)
          
        })
      }
    })
    DT::selectRows(proxy = dataTableProxy("recapitulation_dt"), NULL)
    
  })
  
  #tab 3
  lasso <- reactiveVal()
  lasso_annotations <- reactiveVal()
  val_score <- reactiveVal(NULL)
  
  
  observeEvent(input$tissue_type3, {
    if (input$tissue_type3 == "Liver") {
      lasso(liver_lasso)
      lasso_annotations(lasso_annotations_liver)
      val_score(vs_liver)
    } else {
      lasso(bat_lasso)
      lasso_annotations(lasso_annotations_bat)
      val_score(vs_bat)
    }
  })
  
  # Reactive expression to process lasso data
  lasso_processed <- reactive({
    lasso_data <- lasso()
    
    # Create new column gene.symbol.conc
    lasso_data$gene.symbol.conc <- paste(lasso_data$gene.symbol, lasso_data$uniprot.id, sep = "_")
    
    # Add a new row with all "-"
    new_row <- data.frame(matrix("-", ncol = ncol(lasso_data)))
    colnames(new_row) <- colnames(lasso_data)
    lasso_data <- rbind(new_row, lasso_data)
    
    return(lasso_data)
  })
  
  # Process searches
  autocomplete_choices_prediction <- reactive({
    lasso_data <- lasso_processed()
    if (input$search_type3 == "Protein (Gene Symbol)") {
      unique(lasso_data$gene.symbol.conc)
    } else if (input$search_type3 == "Protein (Uniprot ID)") {
      unique(lasso_data$uniprot.id)
    } else if (input$search_type3 == "Metabolite") {
      metab <- colnames(lasso_data)
      metab <- metab[-c(1:2, length(metab))]
      metab <- c("-", metab)
      unique(metab)
    }
  })
  
  # Update select input choices using an observer
  observe({
    updateSelectizeInput(session, "search_predictions", choices = autocomplete_choices_prediction(), server = TRUE)
  })
  
  
  observeEvent(input$search_button3, {
    search_predictions <- input$search_predictions
    
    if(input$search_type3 == "Metabolite"){
      
      if(search_predictions %in% colnames(lasso())) {
        met_lasso <- subset(lasso(), get(search_predictions) != 0, select = c(uniprot.id, gene.symbol, get(search_predictions)))
        met_lasso3 <- left_join(met_lasso, prot_descript_all, by = c("uniprot.id" = "Entry"))
        colnames(met_lasso3)[8] <- "Description"
        
        #met_lasso3_mito <- left_join(met_lasso3, mito_proteins)
        #met_lasso3_mito_trans <- left_join(met_lasso3_mito, transporter_proteins())
        #met_lasso3_mito_trans_meta <- left_join(met_lasso3_mito_trans, metabolic_proteins())
        
        met_lasso3_annotated <- left_join(
          met_lasso3, 
          lasso_annotations()[, c(2, 5, 7:9, 15)], 
          by = c("gene.symbol" = "Gene.symbol")
        )
        
        met_lasso3_annotated <- met_lasso3_annotated %>%
          filter(Variable == search_predictions)
        
        met_lasso3_annotated$Color <- NA
        
        # Apply the color blending logic to the dataframe
        met_lasso3_annotated <- met_lasso3_annotated %>%
          mutate(Color = case_when(
            Outlier == "Extreme Outlier" ~ "red",
            mito == "yes" & transporter == "yes" & meta == "yes" ~ "saddlebrown",
            mito == "yes" & transporter == "yes" ~ "plum",
            mito == "yes" & meta == "yes" ~ "orange",
            transporter == "yes" & meta == "yes" ~ "darkseagreen",
            mito == "yes" ~ "pink",
            transporter == "yes" ~ "lightblue",
            meta == "yes" ~ "lemonchiffon",
            TRUE ~ "grey"
          ))
        
        met_lasso3_annotated$Annotation <- NA
        met_lasso3_annotated <- met_lasso3_annotated %>%
          mutate(Annotation = case_when(
            mito == "yes" & transporter == "yes" & meta == "yes" & Outlier == "Extreme Outlier" ~ "Mitochondrial, Transport, Metabolic, Extreme Outlier",
            mito == "yes" & transporter == "yes" & meta == "yes" ~ "Mitochondrial, Transport, Metabolic",
            
            mito == "yes" & transporter == "yes" & Outlier == "Extreme Outlier" ~ "Mitochondria, Transport, Extreme Outlier",
            mito == "yes" & transporter == "yes" ~ "Mitochondria, Transport",
            
            mito == "yes" & meta == "yes" & Outlier == "Extreme Outlier" ~ "Mitochondrial, Metabolic, Extreme Outlier",
            mito == "yes" & meta == "yes" ~ "Mitochondrial, Metabolic",
            
            transporter == "yes" & meta == "yes" & Outlier == "Extreme Outlier" ~ "Transport, Metabolic, Extreme Outlier",
            transporter == "yes" & meta == "yes" ~ "Transport, Metabolic",
            
            mito == "yes" & Outlier == "Extreme Outlier" ~ "Mitochondrial, Extreme Outlier",
            mito == "yes" ~ "Mitochondrial",
            
            transporter == "yes" & Outlier == "Extreme Outlier" ~ "Transport, Extreme Outlier",
            transporter == "yes" ~ "Transport",
            
            meta == "yes" & Outlier == "Extreme Outlier" ~ "Metabolic, Extreme Outlier",
            meta == "yes" ~ "Metabolic",
            
            Outlier == "Extreme Outlier" ~ "Extreme Outlier",
            
            TRUE ~ ""
          ))
        met_lasso3_annotated <- unique(met_lasso3_annotated)
        met_lasso3_annotated2 <- met_lasso3_annotated[ , -c(4:9)]
        met_lasso3_annotated2 <- met_lasso3_annotated2[complete.cases(met_lasso3_annotated2), ]
        
        output$ml_plot <- renderPlotly({
          plot_ly(data = met_lasso3_annotated2, x = ~get(search_predictions), y = ~reorder(gene.symbol, get(search_predictions)),
                  type = 'bar', marker = list(color = ~Color), hoverinfo = 'x+y+text', text = met_lasso3_annotated2$Annotation) %>%
            layout(title = paste(search_predictions, "Protein Predictors"),
                   xaxis = list(title = "Coefficient"), yaxis = list(title = "Protein Predictors", automargin = FALSE, tickmode = "linear", tickfont = list(size = 8)),
                   height = 480, width = 755)
        })
        
      } else {
        output$ml_plot <- renderPlotly({
          dummy_plot <- plot_ly(data = data.frame(x = numeric(0), y = numeric(0)), type = 'scatter', mode = 'markers') %>%
            layout(
              title = "Invalid search query.",
              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
            )
        })
      }
      
      
    }else{
      lassot <- t(lasso())
      row_names <- rownames(lassot)
      lassot <- cbind(metabolite = row_names, lassot)
      lassot <- lassot[-2, ]
      
      if(input$search_type3 == "Protein (Gene Symbol)") {
        colnames(lassot) <- lassot[1, ]
        lassot <- lassot[-c(1,nrow(lassot)), ]
        colnames(lassot)[1] <- "metabolite"
        prediction_split <- strsplit(search_predictions, "_")[[1]]
        prediction_search_value <- prediction_split[2]
        gene_uni <- search_predictions
      } else if(input$search_type3 == "Protein (Uniprot ID)") {
        colnames(lassot) <- lassot[1, ]
        lassot <- lassot[-c(1,nrow(lassot)), ]
        colnames(lassot)[1] <- "metabolite"
        prediction_search_value <- search_predictions
        gene_uni <- lasso_processed()[lasso_processed()$uniprot.id == search_predictions, "gene.symbol.conc"]
      }
      
      if(prediction_search_value %in% colnames(lassot)) {
        prot_lasso <- as.data.frame(lassot[ , c("metabolite", prediction_search_value)])
        prot_lasso <- prot_lasso[as.numeric(prot_lasso[, 2]) != 0, ]
        if(nrow(prot_lasso) == 0) {
          output$ml_plot <- renderPlotly({
            # Create a dummy plotly plot
            dummy_plot <- plot_ly(data = data.frame(x = numeric(0), y = numeric(0)), type = 'scatter', mode = 'markers') %>%
              layout(
                title = paste0(gene_uni, " does not have a nonzero coefficient in predicting any metabolite"),
                xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
              )
          })
        } else {
          output$ml_plot <- renderPlotly({
            # Convert prot_lasso to a dataframe if it's not already
            prot_lasso <- as.data.frame(prot_lasso)
            # Ensure that search_predictions is converted to numeric
            prot_lasso[[prediction_search_value]] <- as.numeric(prot_lasso[[prediction_search_value]])
            
            prot_lasso_ordered <- prot_lasso %>%
              arrange(!!sym(prediction_search_value))  # Arrange in ascending order
            
            # Convert metabolites to a factor with the desired order
            prot_lasso_ordered$metabolite <- factor(prot_lasso_ordered$metabolite,
                                                    levels = prot_lasso_ordered$metabolite)
            # Plot the ordered dataframe
            plot_ly(data = prot_lasso_ordered, x = ~get(prediction_search_value), y = ~metabolite,
                    type = 'bar', orientation = 'h', marker = list(color = 'lemonchiffon'), hoverinfo = 'x+y') %>%
              layout(title = paste(gene_uni, "Metabolite Predictors"),
                     yaxis = list(title = "Metabolite Predictors"), xaxis = list(title = "Coefficient"),
                     height = 480, width = 755)
          })
          
        }
      } else{
        output$ml_plot <- renderPlotly({
          dummy_plot <- plot_ly(data = data.frame(x = numeric(0), y = numeric(0)), type = 'scatter', mode = 'markers') %>%
            layout(
              title = "Invalid search query.",
              xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
              yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
            )
        })
      }
      
      
    }
    
  })
  
  # Process searches
  autocomplete_choices_validations <- reactive({
    vs <- val_score()
    vs_met <- vs[ , 1]
    blank <- data.frame(Variable = "-")
    vs_met <- rbind(blank, vs_met)
    colnames(vs_met) <- " "
    unique(vs_met)
  })
  
  # Update select input choices using an observer
  observe({
    updateSelectizeInput(session, "search_validations", choices = autocomplete_choices_validations())
  })
  
  
  observeEvent(input$search_button4, {
    search_validations <- input$search_validations
    vs <- val_score()[!rowSums(is.na(val_score()[, (ncol(val_score())-2):ncol(val_score())])) == 3, ]
    colnames(vs) <- c("Variable", "Newfound", "Established", "Score")
    if(search_validations %in% vs$Variable) {
      vs$Color <- ifelse(vs$Variable == search_validations, "orchid", "grey")
      vs <- vs[order(vs$Newfound, decreasing = TRUE), ]
      vs <- vs[order(vs$Score, decreasing = TRUE), ]
      vs$Rank <- 1:nrow(vs)
      search_hit_row <- which(vs$Variable == input$search_validations)
      
      # Calculate the start and end rows for the range (20 before and after)
      start_row <- max(0, search_hit_row - 20)  # Ensure the start row is at least 1
      end_row <- min(nrow(vs), search_hit_row + 20)
      output$vs_plot <- renderPlotly({
        plot_ly(data = vs, 
                x = ~Rank, 
                y = ~Score, 
                type = 'bar', 
                text = ~paste("Metabolite:", vs$Variable, 
                              "<br>Validation Score:", vs$Score,
                              "<br>Established:", vs$Established, 
                              "<br>Newfound:", vs$Newfound,
                              "<br>Rank:", vs$Rank),
                textposition = 'none',
                hoverinfo = 'text',  # Display x (Rank), y (Score), and text (Variable) on hover
                marker = list(color = ~Color)  # Color bars based on the Color column
        ) %>%
          layout(
            title = paste(search_validations, "Validation Score"),
            xaxis = list(title = "Rank", tickmode = "linear", range = c(start_row, end_row)),
            yaxis = list(title = "Validation Score"),
            bargap = 0.2  # Space between bars
          )
      })        
    }else{
      output$vs_plot <- renderPlotly({
        dummy_plot <- plot_ly(data = data.frame(x = numeric(0), y = numeric(0)), type = 'scatter', mode = 'markers') %>%
          layout(
            title = "This metabolite was not found in this tissue.",
            xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
            yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE)
          )
      })
    }
    
  })
  
  #tab 3.5
  met_abun <- reactiveVal()
  prot_abun <- reactiveVal()
  abun <- reactiveVal()
  axis <- reactiveVal()
  total <- reactiveVal()
  total_pre <- reactiveVal()
  met_total_pre <- reactiveVal()
  prot_total_pre <- reactiveVal()
  met_total <- reactiveVal()
  prot_total <- reactiveVal()
  #cv_matrix_total <- reactiveVal()
  
  
  observeEvent(input$tissue_type35, {
    if (input$tissue_type35 == "Liver") {
      met_abun(met_abun_liver)
      prot_abun(prot_abun_liver)
      
      hold <- prot_abun_liver
      firstrow <- as.data.frame(t(as.data.frame(sapply(hold[1, ], function(x) strsplit(as.character(x), "_")[[1]][1]))))
      hold <- rbind(firstrow, hold)
      hold <- hold[-2, ]
      rownames(hold) <- NULL
      prot_abun(hold)
      
      calculate_cv <- function(row) {
        non_na_values <- na.omit(as.numeric(row))
        if (length(non_na_values) < 2) {
          return(NA)  # Return NA if there are not enough non-NA values for CV calculation
        }
        non_na_values2 <- 2^non_na_values
        cv <- sd(non_na_values2) / mean(non_na_values2) * 100
        return(cv)
      }
      
      met_total_t <- as.data.frame(t(met_abun_liver))
      colnames(met_total_t) <- met_total_t[1, ]
      met_total_t <- met_total_t[-1, ]
      cv_matrix_met <- apply(met_total_t, 1, calculate_cv)
      cv_matrix_met <- cv_matrix_met[!is.na(cv_matrix_met)]
      met_total(cv_matrix_met)
      met_total_pre(met_total_t)
      
      prot_total_t <- as.data.frame(t(prot_abun_liver))
      colnames(prot_total_t) <- prot_total_t[1, ]
      prot_total_t <- prot_total_t[-1, ]
      cv_matrix_prot <- apply(prot_total_t, 1, calculate_cv)
      cv_matrix_prot <- cv_matrix_prot[!is.na(cv_matrix_prot)]
      prot_total(cv_matrix_prot)
      prot_total_pre(prot_total_t)
      
      
    } else {
      met_abun(met_abun_bat)
      prot_abun(prot_abun_bat)
      
      hold <- prot_abun_bat
      firstrow <- as.data.frame(t(as.data.frame(sapply(hold[1, ], function(x) strsplit(as.character(x), "_")[[1]][1]))))
      hold <- rbind(firstrow, hold)
      hold <- hold[-2, ]
      rownames(hold) <- NULL
      prot_abun(hold)
      
      calculate_cv <- function(row) {
        non_na_values <- na.omit(as.numeric(row))
        if (length(non_na_values) < 2) {
          return(NA)  # Return NA if there are not enough non-NA values for CV calculation
        }
        non_na_values2 <- 2^non_na_values
        cv <- sd(non_na_values2) / mean(non_na_values2) * 100
        return(cv)
      }
      
      met_total_t <- as.data.frame(t(met_abun_bat))
      colnames(met_total_t) <- met_total_t[1, ]
      met_total_t <- met_total_t[-1, ]
      cv_matrix_met <- apply(met_total_t, 1, calculate_cv)
      cv_matrix_met <- cv_matrix_met[!is.na(cv_matrix_met)]
      met_total(cv_matrix_met)
      met_total_pre(met_total_t)
      
      prot_total_t <- as.data.frame(t(prot_abun_bat))
      colnames(prot_total_t) <- prot_total_t[1, ]
      prot_total_t <- prot_total_t[-1, ]
      cv_matrix_prot <- apply(prot_total_t[ , -1], 1, calculate_cv)
      cv_matrix_prot <- cv_matrix_prot[!is.na(cv_matrix_prot)]
      prot_total(cv_matrix_prot)
      prot_total_pre(prot_total_t)
    }
    
    updateSelectizeInput(session, "search_abun", selected = "")
    
  })
  
  
  metabolites_abun <- reactive({
    df <- data.frame(list = colnames(met_abun()))
    df[1, 1] <- "-"
    return(df)
  })
  proteins_abun_uniprot <- reactive({
    df <- data.frame(list = colnames(prot_abun()))
    df[1, 1] <- "-"
    return(df)
  })
  proteins_abun_gene <- reactive({
    col_names <- colnames(prot_abun())
    #first_row <- sapply(prot_abun()[1, ], function(x) strsplit(as.character(x), "_")[[1]][1])
    first_row <- prot_abun()[1, ]
    new_names <- paste(first_row, col_names, sep = "_")
    df <- data.frame(list = new_names)
    df[1, 1] <- "-"
    rownames(df) <- NULL
    df <- na.omit(df)
    return(df)
  })
  
  # Update select input choices using an observer
  observe({
    updateSelectizeInput(session, "search_abun", choices = autocomplete_choices_abun(), server = TRUE)
  })
  
  # Reactive expression to get unique metabolites for autocomplete
  autocomplete_choices_abun <- reactive({
    if(input$search_type35 == "Metabolite"){
      unique(metabolites_abun()$list)
    }else if(input$search_type35 == "Protein (Gene Symbol)"){
      unique(proteins_abun_gene()$list)
    }else if(input$search_type35 == "Protein (Uniprot ID)"){
      unique(proteins_abun_uniprot()$list)
    }
  })
  
  observeEvent(input$search_button35, {
    search_abun <- input$search_abun
    
    # Set the appropriate abundance data based on the selected type
    if (input$search_type35 == "Metabolite") {
      abun(met_abun())
      total(met_total())
      total_pre(met_total_pre())
      axis("Abundance (log2 sample-to-pool ratio)")
    } else if(input$search_type35 == "Protein (Uniprot ID)") {
      abun(prot_abun()[-1, ])
      total(prot_total())
      total_pre(prot_total_pre())
      axis("Abundance (log2 sample-to-bridge ratio)")
    } else if(input$search_type35 == "Protein (Gene Symbol)") {
      gs <- prot_abun()[1, ]
      us <- colnames(prot_abun())
      gs_us <- paste(gs, us, sep = "_")
      
      pa <- prot_abun()
      colnames(pa) <- gs_us
      colnames(pa)[1] <- "analyte"
      abun(pa[-1, ])
      
      total(prot_total())
      
      gs2 <- prot_total_pre()[ , 1]
      us2 <- rownames(prot_total_pre())
      conc <- paste(gs2, us2, sep = "_")
      
      symbol_switch <- as.data.frame(prot_total_pre())
      rownames(symbol_switch) <- conc
      symbol_switch <- symbol_switch[ , -1]
      total_pre(symbol_switch)
      
      axis("Abundance (log2 sample-to-bridge ratio)")
    }
    
    
    output$abun_text <- renderText({
      # Check if the selected entry exists in the dataset
      if (search_abun %in% colnames(abun())) {
        " "
      } else {
        "Entry not found in dataset."
      }
    })
    
    output$dist_text <- renderText({
      # Check if the selected entry exists in the dataset
      if (search_abun %in% colnames(abun())) {
        " "
      } else {
        "Entry not found in dataset."
      }
    })
    
    output$abun_plot <- renderPlotly({
      if (search_abun %in% colnames(abun())) {
        # Create the bar chart with hoverinfo
        abun1 <- abun()
        abun1[ , -1] <- lapply(abun1[ , -1], as.numeric)
        plot_ly(data = abun1, x = ~analyte, y = ~.data[[search_abun]], type = 'bar', 
                marker = list(color = 'navy'), 
                hoverinfo = 'text',
                text = ~paste(analyte, "<br>", .data[[search_abun]])) %>%
          layout(title = paste("Abundance of", search_abun),
                 xaxis = list(title = "Mouse"), yaxis = list(title = axis(), zeroline = TRUE),
                 height = 480, width = 740)
      } else{
        NULL
      }
    })
    
    
    ## distribution tab ##
    
    calculate_cv <- function(row) {
      non_na_values <- na.omit(as.numeric(row))
      if (length(non_na_values) < 2) {
        return(NA)  # Return NA if there are not enough non-NA values for CV calculation
      }
      non_na_values2 <- 2^non_na_values
      cv <- sd(non_na_values2) / mean(non_na_values2) * 100
      return(cv)
    }
    
    total_t <- as.data.frame(total_pre())
    analyte <- as.data.frame(total_t[search_abun, ])
    cv_matrix2 <- apply(analyte, 1, calculate_cv)
    cv_matrix2 <- cv_matrix2[!is.na(cv_matrix2)]
    
    cv_matrix1 <- total()
    
    if (input$search_type35 == "Metabolite") {
      # Combine the data into a data frame
      data <- data.frame(matrix = rep(c("All Metabolites", search_abun), times = c(length(cv_matrix1), length(cv_matrix2))),
                         cv = c(cv_matrix1, cv_matrix2))
      data$cv_graph <- ifelse(data$matrix == "All Metabolites", data$cv, ifelse(data$cv > 200, 200, data$cv))
      data$color <- ifelse(data$matrix == search_abun, "grey", ifelse(data$matrix == "All Metabolites", "firebrick3", "default_color") )
      data$matrix <- factor(data$matrix, levels = unique(data$matrix))
      
      output$dist_plot <- renderPlot({
        if (search_abun %in% colnames(abun())) {
          ggplot(data, aes(x = cv_graph, fill = matrix)) +
            geom_density(alpha = 0.8, color = NA) +
            labs(
              x = "Coefficient of Variation (%)",
              y = "Density"
            ) +
            theme_minimal() +
            scale_fill_manual(values = unique(data$color)) +
            geom_vline(
              xintercept = mean(cv_matrix1),
              linetype = "solid",
              color = "firebrick",
              size = 1
            ) +
            geom_vline(
              xintercept = mean(as.numeric(data[data[ , "matrix"] == search_abun, "cv_graph"])),
              linetype = "dashed",
              color = "grey",
              size = 1
            ) +
            annotate(
              "text",
              x = mean(cv_matrix1),
              y = 0.02,
              label = paste("Avg CV All Metabolites:", round(mean(cv_matrix1), 2)),
              color = "firebrick4",
              vjust = -0.5,
              size = 4
            ) +
            annotate(
              "text",
              x = mean(data[data[ , "matrix"] == search_abun, "cv_graph"] - 15),
              y = 0.015,
              label = paste(search_abun, "CV:", round(mean(cv_matrix2), 2)),
              color = "black",
              vjust = -0.5,
              size = 4
            ) +
            theme(
              axis.title = element_text(color = "black", size = 15),
              axis.text = element_text(color = "black"),
              axis.line = element_line(color = "black"),  # Set axis line color
              panel.grid = element_blank(),  # Remove inner gridlines
              plot.background = element_rect(fill = "white"),  # Set plot background color
              legend.position = "bottom",  # Move legend below the graph
              legend.title = element_text(size = 12, face = "bold"),
              legend.text = element_text(size = 10)
            ) +
            scale_x_continuous(limits = c(-5, 207))  # Set x-axis limits
          
        } else {
          NULL
        }
      }) 
    } else {
      data <- data.frame(matrix = rep(c("All Proteins", search_abun), times = c(length(cv_matrix1), length(cv_matrix2))),
                         cv = c(cv_matrix1, cv_matrix2))
      data$cv_graph <- ifelse(data$matrix == "All Proteins", data$cv, ifelse(data$cv > 130, 130, data$cv))
      data$color <- ifelse(data$matrix == search_abun, "grey", ifelse(data$matrix == "All Proteins", "firebrick3", "default_color") )
      data$matrix <- factor(data$matrix, levels = unique(data$matrix))
      
      output$dist_plot <- renderPlot({
        if (search_abun %in% colnames(abun())) {
          ggplot(data, aes(x = cv_graph, fill = matrix)) +
            geom_density(alpha = 0.8,color=NA) +
            labs(x = "Coefficient of Variation",
                 y = "Density") +
            theme_minimal() +
            scale_fill_manual(values = unique(data$color)) +  
            geom_vline(xintercept = mean(cv_matrix1), linetype = "solid", color = "firebrick", size = 1) +
            geom_vline(xintercept = mean(as.numeric(data[data[ , "matrix"] == search_abun, "cv_graph"])), linetype = "dashed", color = "grey", size = 1) +
            annotate("text", x = mean(cv_matrix1+15), y = 0.02, label = paste("Avg CV All Proteins:", round(mean(cv_matrix1), 2)), color = "firebrick4", vjust = -0.5, size = 4) +
            annotate("text", x = mean(data[data[ , "matrix"] == search_abun, "cv_graph"]-10), y = 0.015, label = paste(search_abun, "CV:", round(mean(cv_matrix2), 2)), color = "black", vjust = -0.5, size = 4) +
            theme(axis.title = element_text(color = "black", size = 15),
                  axis.text = element_text(color = "black"),
                  axis.line = element_line(color = "black"),  # Set axis line color
                  panel.grid = element_blank(),  # Remove inner gridlines
                  plot.background = element_rect(fill = "white")) +  # Set plot background color
            scale_x_continuous(limits = c(-5, 133))  # Set x-axis limits
        } else{
          NULL
        }
      }) 
    }
    
  })
  
  
  #tab4
  output$download_btn_1 <- downloadHandler(
    filename = function() {
      paste("bat_pearson.csv", sep = "")
    },
    content = function(file) {
      write.csv(bat_pearson_data, file)
    }
  )
  output$download_btn_2 <- downloadHandler(
    filename = function() {
      paste("liver_pearson.csv", sep = "")
    },
    content = function(file) {
      write.csv(liver_pearson_data, file)
    }
  )
  output$download_btn_4 <- downloadHandler(
    filename = function() {
      paste("rhea_bat.csv", sep = "")
    },
    content = function(file) {
      write.csv(rhea_bat, file)
    }
  )
  output$download_btn_5 <- downloadHandler(
    filename = function() {
      paste("reactome_bat.csv", sep = "")
    },
    content = function(file) {
      write.csv(reactome_bat, file)
    }
  )
  output$download_btn_6 <- downloadHandler(
    filename = function() {
      paste("tcdb_bat.csv", sep = "")
    },
    content = function(file) {
      write.csv(tcdb_bat, file)
    }
  )
  output$download_btn_7 <- downloadHandler(
    filename = function() {
      paste("rhea_liver.csv", sep = "")
    },
    content = function(file) {
      write.csv(rhea_liver, file)
    }
  )
  output$download_btn_8 <- downloadHandler(
    filename = function() {
      paste("reactome_liver.csv", sep = "")
    },
    content = function(file) {
      write.csv(reactome_liver, file)
    }
  )
  output$download_btn_9 <- downloadHandler(
    filename = function() {
      paste("tcdb_liver.csv", sep = "")
    },
    content = function(file) {
      write.csv(tcdb_liver, file)
    }
  )
  output$download_btn_13 <- downloadHandler(
    filename = function() {
      paste("bat_lasso.csv", sep = "")
    },
    content = function(file) {
      write.csv(bat_lasso, file)
    }
  )
  output$download_btn_14 <- downloadHandler(
    filename = function() {
      paste("liver_lasso.csv", sep = "")
    },
    content = function(file) {
      write.csv(liver_lasso, file)
    }
  )
  
}

shinyApp(ui = ui, server = server)