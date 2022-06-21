library(shiny)
library(tinytex)
library(magrittr)
library(DT) 
library(ggplot2)
library(purrr)
library(shinythemes)
library(plotly)
library(dplyr)
library(tidyverse)
library(knitr)

##Author: ELENAOR SEABY##

#load in GATK tables
Metrics_Collection <- read.table(file = 'MetricsCollection.sp.tsv_table', sep = '\t', header = TRUE)
CompOverlap <- read.table(file = 'CompOverlap.sp.tsv_table', sep = '\t', header = TRUE)
Count_Variants <- read.table(file = 'CountVariants.sp.tsv_table', sep = '\t', header = TRUE)
TiTv <- read.table(file = 'TiTvVariantEvaluator.sp.tsv_table', sep = '\t', header = TRUE)
Multiallelic_Summary <- read.table(file = 'MultiallelicSummary.sp.tsv_table', sep = '\t', 
                                   header = TRUE)
Indel_Summary <- read.table(file = 'IndelSummary.sp.tsv_table', sep = '\t', header = TRUE)
All_tables <- read.table(file = 'joint_called.QC', sep = '\t', header = TRUE) 

#load in picard tables
Detail_metrics <- read.table(file = 'joint_called.variant_calling_detail_metrics.table', 
                             sep = '\t', header = TRUE)
Summary_metrics <- read.table(file = 'joint_called.variant_calling_summary_metrics.table', 
                              sep = '\t', header = TRUE)

#load indel histogram
Indel_histogram <- read.table(file = 'histogram.tsv', 
                              sep = '\t', header = TRUE)

#modify Detail_metrics
Detail_metrics_short <- Detail_metrics %>% select(-c("FILTERED_SNPS", "FILTERED_INDELS", "NUM_IN_DB_SNP_INDELS",
                                                     "DBSNP_INS_DEL_RATIO", "NOVEL_INS_DEL_RATIO",
                                                     "PCT_GQ0_VARIANTS", "PCT_DBSNP_INDELS",
                                                     "NUM_IN_DB_SNP_MULTIALLELIC", "NUM_IN_DB_SNP_COMPLEX_INDELS",
                                                     "SNP_REFERENCE_BIAS", "NUM_SINGLETONS"))

colnames(Detail_metrics_short) <- c("SAMPLE",	"HET:HOM",	"No_GQ0_VAR",	"No_SNPs", "No_in_DBSNP",
                                    "NOVEL_SNPs",	"%_DBSNP",	"DBSNP_TITV",	"NOVEL_TITV",	"No_INDELS","NOVEL_INDELS",
                                    "No_MULTIALLELICs",	"No_COMPLEX_INDELS")

#load peddy
tmp <- list.files(pattern="*het_check.csv")
Ethnicity = lapply(tmp, read.csv)
Ethnicity <- Ethnicity[[1]]

#load predicted ethnicity into detailed metrics
Detail_metrics$ancestry_prediction <- Ethnicity$ancestry.prediction[match(Ethnicity$sample_id, 
                                                                          Detail_metrics$SAMPLE_ALIAS)]

Detail_metrics$ancestry_probability <- Ethnicity$ancestry.prob[match(Ethnicity$sample_id, 
                                                                     Detail_metrics$SAMPLE_ALIAS)]

#load predicted into detailed metrics
tmp <- list.files(pattern="*sex_check.csv")
Sex = lapply(tmp, read.csv)
Sex <- Sex[[1]]
Sex %<>% rename("sample"=sample_id)
Sex %<>% relocate(sample)

#load pca plot
pca <- list.files(pattern="*pca_check.png")
#move file to www folder for it to render correctly
dir.create("www")
file.copy(pca[[1]], "www")

#load pca html page
html_peddy <- list.files(pattern="*.html")

#load predicted ethnicity into detailed metrics
Detail_metrics$sex_prediction <- Sex$predicted_sex[match(Sex$sample, 
                                                         Detail_metrics$SAMPLE_ALIAS)]
###########################
# make relatedness matrix #
###########################
tmp <- list.files(pattern="*ped_check.csv")
Relatedness = lapply(tmp, read.csv)
Relatedness <- Relatedness[[1]]
related_matrix <- Relatedness %>% select(sample_a, sample_b, rel)

un2 <- sort(unique(unlist(related_matrix[1:2])))
out2_new <- related_matrix %>% 
  complete(sample_a = un2, sample_b = un2) %>% 
  pivot_wider(names_from = sample_b, values_from = rel)

tmp <- map2_dfc(data.table::transpose(out2_new, make.names = 'sample_a'), 
                out2_new[-1], coalesce) %>% 
  bind_cols(out2_new %>%
              select(sample_a), .)

tmp2 <- column_to_rownames(tmp, var = "sample_a")
tmp2 %<>% as.matrix(tmp2, rownames.force=TRUE)

#make plotly heatmap
relatedness_heatmap <- plot_ly(x=colnames(tmp2), y=rownames(tmp2), z = tmp2, type = "heatmap")

#make ggplot heatmap
rel_mat2 <- related_matrix %>% select(sample_b, sample_a, rel)
colnames(rel_mat2) <- c("X", "Y", "rel")
tmp <- related_matrix
colnames(tmp) <- c("X", "Y", "rel")

tmp2 <- bind_rows(tmp, rel_mat2)

relatedness_heatmap2 <- ggplot(tmp2, aes(X, Y, fill= rel)) + 
  geom_tile()  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=12))


#generate indel histogram
Indel_Histogram <- ggplot(Indel_histogram, aes(Length, Freq)) + geom_col()

######################
# deal with coverage #
######################

#load in coverage.hist files
coverage_hist = list.files(pattern="*_coverage.hist")
for (i in 1:length(coverage_hist)) assign(coverage_hist[i], read.table(coverage_hist[i]))
hist_files = lapply(coverage_hist, read.table)

#load in coverage.stats files
coverage_stats = list.files(pattern="*_coverage.stats")
#for (i in 1:length(coverage_stats)) assign(coverage_stats[i], read.table(coverage_stats[i]))
cov_stats_files = lapply(coverage_stats, read.table, fill=T)

#plot histogram function
plot_histogram <- function(x) {
  
  x <- ggplot(x, aes(V2, V5)) + geom_col() + xlim(0,500) + ylim(0,0.01)+
    ylab("Proportion") + xlab("Read Depth")
  
}

#lapply(hist_files, plot_histogram)

#handle coverage histogram for facet wrap
histo_prep <- function(x) {
  
  x %<>% select(V2, V5, V6) %>% filter(V2 <= 500)
  x %<>% rename("sample"="V6")
  x %<>% column_to_rownames("V2")
  x %<>% rownames_to_column(var = "depth")
}

hist_tables <- lapply(hist_files, histo_prep)

#replace with sample IDs
coverage_hist <- str_replace_all(coverage_hist, "_coverage.hist", "")
names(hist_files) = coverage_hist

#now merge tables
tmp <- Reduce(rbind, hist_tables)
#make depth numeric
tmp$depth <- as.numeric(tmp$depth)

#Generate coverage histograms
coverage_histograms <- tmp %>% group_by(sample) %>% ggplot(aes(depth, V5)) + geom_col() + xlim(0,500) + ylim(0,0.03)+
  ylab("Proportion") + xlab("Read Depth") + facet_wrap(~ sample) + theme(axis.text=element_text(size=12))

#Deal with cov stats table
#First transform tables
transform_tables <- function(x) {
  x %<>% t()
  x <- as.data.frame(x)
  x %<>% setNames(as.character(x[1,]))
  x <- x[-1,] 
}

cov_stats_files <- lapply(cov_stats_files, transform_tables)

#now merge tables
cov_tables <- Reduce(rbind, cov_stats_files)
cov_tables %<>% remove_rownames()

colnames(cov_tables) <- c("sample",	"no_reads",	"mapped_to_target",	"%",	"mapped_+150bp",
                          "%",	"mean_cov", "accessible_target_bases",	"accessible_1x",	"%",
                          "accessible_5x","%", "accessible_10x", "%", "target_20x",	"%")


##################
# Make Shiny App #
##################

ui <- fluidPage(theme = shinytheme("united"),
                titlePanel("QC output"),
                navbarPage("Menu",
                           tabPanel("QC tables",
                                    sidebarLayout(
                                      sidebarPanel(
                                        selectInput("QC", "QC tables",
                                                    choices = c("Summary metrics" ="Summary_metrics", 
                                                                "Detailed metrics" = "Detail_metrics", 
                                                                "Sex" ="Sex", 
                                                                "Ancestry" ="Ethnicity", 
                                                                "Relatedness" ="Relatedness", 
                                                                "All tables" ="All_tables", 
                                                                "CompOverlap" ="CompOverlap", 
                                                                "Metrics collection" ="Metrics_Collection", 
                                                                "Count variants" ="Count_Variants", 
                                                                "TiTv" ="TiTv",
                                                                "Multiallelic summary" ="Multiallelic_Summary", 
                                                                "Indel summary" ="Indel_Summary", 
                                                                "Coverage" ="cov_tables")),
                                        selectInput("plotInput", "Detailed metrics variable",
                                                    choices = c("HET_HOMVAR_RATIO", "PCT_GQ0_VARIANTS", "TOTAL_GQ0_VARIANTS",
                                                                "TOTAL_HET_DEPTH", "TOTAL_SNPS", "NUM_IN_DB_SNP", "NOVEL_SNPS",
                                                                "FILTERED_SNPS", "PCT_DBSNP", "DBSNP_TITV", "NOVEL_TITV",
                                                                "TOTAL_INDELS", "NOVEL_INDELS", "FILTERED_INDELS","PCT_DBSNP_INDELS",
                                                                "NUM_IN_DB_SNP_INDELS", "DBSNP_INS_DEL_RATIO", "NOVEL_INS_DEL_RATIO",
                                                                "TOTAL_MULTIALLELIC_SNPS", "NUM_IN_DB_SNP_MULTIALLELIC", 
                                                                "TOTAL_COMPLEX_INDELS", "NUM_IN_DB_SNP_COMPLEX_INDELS",
                                                                "SNP_REFERENCE_BIAS", "NUM_SINGLETONS")),
                                        selectInput("fill", "Plot fill",
                                                    choices = c("Sex prediction" ="sex_prediction", "Ancestry prediction" ="ancestry_prediction")), width=3,
                                        a(img(src=pca[[1]], height = 500, width = 300, slign="center", 
                                              target="_blank"), 
                                          href=pca[[1]])),
                                      mainPanel(
                                        dataTableOutput("results"),
                                        br(),
                                        br(),
                                        "Detailed Metrics Summary Plot",
                                        plotlyOutput("scatter", height=400),
                                        br(),
                                        verbatimTextOutput("info"), width=9)
                                    )),
                           tabPanel("Plots",
                                    sidebarLayout(
                                      sidebarPanel(
                                        selectInput("more_plots", "Select Plot",
                                                    choices = c("Relatedness matrix" ="relatedness_heatmap", 
                                                                "Indel histogram" ="Indel_Histogram", 
                                                                "Coverage histogram" ="coverage_histogram"))
                                        , width=4),
                                      mainPanel(uiOutput("my_plot"), height="100%", width=8))
                           ),
                           tabPanel("Coverage",
                                    sidebarLayout(
                                      sidebarPanel(selectInput("histo", "Select sample",
                                                               choices = coverage_hist),
                                                   br(),
                                                   plotOutput("histograms")),
                                      mainPanel(dataTableOutput("cov_table"))))
                           ,
                           tabPanel("Report",
                                    sidebarLayout(
                                      sidebarPanel(downloadButton("report", "Generate report"), width=0
                                      ),
                                      mainPanel(tags$h2("Detailed metrics"),
                                        tableOutput("detail"),
                                                br(),
                                        tags$h2("Sex inference"),
                                                tableOutput("sex"),
                                                br(),
                                        tags$h2("Ancestry prediction Peddy"),
                                                tableOutput("ancestry"),
                                                a(img(src=pca[[1]], height = 500, width = 300, slign="center",
                                                      target="_blank"),
                                                  href=pca[[1]]),
                                                br(),
                                        tags$h2("Coverage metrics"),
                                                tableOutput("cove"),
                                        tags$h2("Coverage plot"),
                                                plotOutput("cov", height="700px"),
                                                br(),
                                        tags$h2("Relatedness matrix"),
                                                plotOutput("heat2"))))
                           
                           
                ))

server <- function(input, output) {
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.pdf",
    content = function(file) {
      tmp_dir <- tempdir()
      tempReport <- file.path(tmp_dir, "report.Rmd")
      tmp_pic <- file.path(tmp_dir, pca[[1]])
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      file.copy(pca[[1]], tmp_pic, overwrite = TRUE)
      # Set up parameters to pass to Rmd document
      
      # Set up parameters to pass to Rmd document
      params <- list(sex = Sex,
                     detail = Detail_metrics_short,
                     ancestry = Ethnicity,
                     coverage = cov_tables,
                     cov_plot = coverage_histograms,
                     relatedness_mat = relatedness_heatmap2,
                     pca = pca[[1]])
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  
  output$scatter <- renderPlotly({
    plot_ly(data = Detail_metrics, x = ~SAMPLE_ALIAS, 
            y = ~get(input$plotInput), color = ~get(input$fill), type="scatter",
            marker = list(size = 12), colors = "Dark2")
  })
  
  #tableOutput facilitates renderTable
  output$detail <- renderTable({Detail_metrics_short}, striped=TRUE, spacing="xs")
  output$sex <- renderTable({Sex}, striped=TRUE, spacing="xs")
  output$ancestry <- renderTable({Ethnicity}, striped=TRUE, spacing="xs")
  output$related <- renderTable({Relatedness}, striped=TRUE, spacing="xs")
  output$cove <- renderTable({cov_tables}, striped=TRUE, spacing="xs")
  
  output$hist <- renderPlot({
    Indel_Histogram  ##  assuming you already did this histogram 
  }, height=600)
  
  output$heat <- renderPlotly({
    relatedness_heatmap  ## assuming you already have this heatmap
  })
  
  output$heat2 <- renderPlot({relatedness_heatmap2}, width=1000, height=800)
  
  output$cov <- renderPlot({
    coverage_histograms  ## assuming you already have this heatmap
  })
  
  output$my_plot <- renderUI({
    if (input$more_plots=="Indel_Histogram"){
      plot <- plotOutput("hist")}
    else if (input$more_plots=="coverage_histogram"){
      plot <- plotOutput("cov", height="600px")}
    else plot <- plotlyOutput("heat",height="600px")
    
  })
  
  output$results <- renderDataTable(
    {get(input$QC)},
    options = list(scrollX = TRUE))
  
  output$cov_table <- renderDataTable(
    {cov_tables},
    options = list(scrollX = TRUE))
  
  output$plots <- renderPlot({ 
    ggplot(Detail_metrics, aes(SAMPLE_ALIAS, get(input$plotInput), color=get(input$fill))) + 
      geom_point(size = 3) +ylab("variable")+
      theme(axis.title.x = element_text(vjust=-2))
  })
  
  output$histograms <- renderPlot({ 
    ggplot(hist_files[[input$histo]], aes(V2, V5)) + geom_col() + xlim(0,500) + ylim(0,0.03)+
      ylab("Proportion") + xlab("Read Depth") 
  })
  
}

shinyApp(ui = ui, server = server)
