library(shiny)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(data.table)
library(DT)
library(htmltools)

# LOAD DATA
#https://www.illumina.com/systems/sequencing-platforms.html
sequencers <- read.csv("data/sequencers.csv")
sequencers2 <- sequencers %>% gather("min_max","output_Gbp",c("min_output_Gbp","max_output_Gbp")) %>%
  mutate(min_max = gsub("_output_Gbp", "", min_max)) %>%
  mutate(output_bp=output_Gbp*(10^9)) %>% 
  mutate(output_reads=ifelse(paired==1,
                             ((output_bp/length/2)*(10^-6)),
                             ((output_bp/length)*(10^-6)))) %>%
  mutate(cost_per_lane=cost/num_lanes) %>%
  mutate(output_Gbp_per_lane=output_Gbp/num_lanes) %>%
  mutate(output_bp_per_lane=output_bp/num_lanes) %>%
  mutate(output_reads_per_lane=output_reads/num_lanes)

{sketch <- withTags(
  table(
    class = "display",
    thead(
      tr(
        th(colspan = 1, "Sequencer", style = "border-right: solid 2px;"),
        th(colspan = 1, "No. Lanes", style = "border-right: solid 2px;"),
        th(colspan = 2, "Flow Cell Gbp", style = "border-right: solid 2px;"),
        th(colspan = 2, "Lane Gbp", style = "border-right: solid 2px;"),
        th(colspan = 1, "Flow Cell Cost", style = "border-right: solid 2px;"),
        th(colspan = 1, "Lane Cost", style = "border-right: solid 2px;"),
        th(colspan = 2, "% Flow Cell", style = "border-right: solid 2px;"),
        th(colspan = 2, "% Lane", style = "border-right: solid 2px;"),
        th(colspan = 2, "% Cost", style = "border-right: solid 2px;"),
        th(colspan = 1, "Overfilled?", style = "border-right: solid 2px;"),
      ),
      tr(
        th(colspan = 1, "", style = "border-right: solid 2px;"),
        th(colspan = 1, "", style = "border-right: solid 2px;"),
        th(colspan = 1, "min"),
        th(colspan = 1, "max", style = "border-right: solid 2px;"),
        th(colspan = 1, "min"),
        th(colspan = 1, "max", style = "border-right: solid 2px;"),
        th(colspan = 1, "", style = "border-right: solid 2px;"),
        th(colspan = 1, "", style = "border-right: solid 2px;"),
        th(colspan = 1, "min"),
        th(colspan = 1, "max", style = "border-right: solid 2px;"),
        th(colspan = 1, "min"),
        th(colspan = 1, "max", style = "border-right: solid 2px;"),
        th(colspan = 1, "min"),
        th(colspan = 1, "max", style = "border-right: solid 2px;"),
        th(colspan = 1, "", style = "border-right: solid 2px;"),
      )
    )
  )
)
  }


# Define UI for dataset viewer app ----
ui <- fluidPage(
  fluidRow(
    column(4,
         img(src = "SequencingEstimator.png", height = 250, width = 240, style="float:right"),
         titlePanel("Sequencing Estimator"),
         h4("Rhett M. Rautsaw"),
         h5(a(href="https://github.com/reptilerhett/SequencingEstimator","github.com/reptilerhett/SequencingEstimator")),
         helpText("This application is designed to help you choose the right sequencing platform and estimate costs."),
         helpText("The app currently provides space for Whole-Genome Sequencing (WGS), RNA-seq, Whole-Genome Bisulfite Sequencing (WGBS), ATAC-Seq, and ddRAD-Seq. However, given desired coverage or number of reads (in million), the fields could be used any type of sequencing."), 
         helpText("Most information was obtained from:",tags$a(href="https://www.illumina.com/systems/sequencing-platforms.html", "Illumina")),
         helpText("The recommended number of reads and coverage are default input for all options")
         ),
    
    column(2,
           # Input: Genome Size
           numericInput(inputId = "genome_size",
                        label = "Genome Size (Gb):",
                        value = 1.6),
            helpText("Not required for RNA-Seq, ATAC-Seq, ddRAD-Seq, or any sequencing type based on 
                     sequencing a specific number or reads rather than coverage."),
           hr(),
           # Input: Read Length
           selectInput(inputId = "read_length_choice", label = "Read Length:", 
                       choices = levels(sequencers$read_length), selected="2 x 150"),
           helpText("Read length is not applicable for PacBio sequencing and will be included regardless."),
           helpText("Illumina sequencers will be filtered by which are applicable with the chosen read length.")),
    
    column(2,
           # Input: Whole Genome Sequencing (WGS)
           numericInput(inputId = "wgs_quantity",
                        label = "WGS Samples:",
                        value = 0),
           numericInput(inputId = "wgs_coverage",
                        label = "WGS Coverage:",
                        value = 30),
           hr(),
           # Input: Transcriptome Sequencing (RNA-Seq)
           numericInput(inputId = "rna_quantity",
                        label = "RNA-Seq Samples:",
                        value = 0),
           numericInput(inputId = "rna_reads",
                        label = "RNA-Seq reads per sample (in millions):",
                        value = 20)),
    
    column(2,
           # Input: Whole Genome Bisulfite Sequencing (WGBS)
           numericInput(inputId = "wgbs_quantity",
                        label = "WGBS Samples:",
                        value = 0),
           numericInput(inputId = "wgbs_coverage",
                        label = "WGBS Coverage:",
                        value = 20),
           hr(),
           # Input: ATAC-Seq
           numericInput(inputId = "atac_quantity",
                        label = "ATAC-Seq Samples:",
                        value = 0),
           numericInput(inputId = "atac_reads",
                        label = "ATAC-Seq reads per sample (in millions):",
                        value = 60)),
    
    column(2,
           # Input: ddRAD-Seq
           numericInput(inputId = "ddrad_quantity",
                        label = "ddRAD-Seq Samples:",
                        value = 0),
           numericInput(inputId = "ddrad_reads",
                        label = "ddRAD-Seq reads per sample (in millions):",
                        value = 5),
           hr(),
           helpText("If you have updated information on sequencing platforms or costs, \
                     you can download an example of the necessary data, update it, and re-upload it here for your needs."),
           downloadButton("downloadData", "Download"),
           helpText(""),
           fileInput("file1", "Choose CSV File",
                     multiple = TRUE,
                     accept = c("text/csv",
                                "text/comma-separated-values,text/plain",
                                ".csv"))
           )
    ),
    
  hr(),
  
  fluidRow(
    column(2,
           h2("Total:"),
           h4(textOutput("txtOutput1")),
           h4(textOutput("txtOutput2")),
           h4(textOutput("txtOutput3"))
           ),
    column(10,
           dataTableOutput("dataTableOutput1"))
  ),
  
  fluidRow(
    column(6,plotOutput("pltOutput1")),
    column(6,plotOutput("pltOutput2"))
  ),
  
  fluidRow(
    column(12,plotOutput("pltOutput3"))
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {
  
  data<-reactive({
    if(is.null(input$file1)){
      data<-sequencers2 %>% filter(is.na(read_length) | read_length==input$read_length_choice) 
      data
    }else{
      req(input$file1)
      
      sequencers <- read.csv(input$file1$datapath)
      
      sequencers2 <- sequencers %>% gather("min_max","output_Gbp",c("min_output_Gbp","max_output_Gbp")) %>%
        mutate(min_max = gsub("_output_Gbp", "", min_max)) %>%
        mutate(output_bp=output_Gbp*(10^9)) %>% 
        mutate(output_reads=ifelse(paired==1,
                                   ((output_bp/length/2)*(10^-6)),
                                   ((output_bp/length)*(10^-6)))) %>%
        mutate(cost_per_lane=cost/num_lanes) %>%
        mutate(output_Gbp_per_lane=output_Gbp/num_lanes) %>%
        mutate(output_bp_per_lane=output_bp/num_lanes) %>%
        mutate(output_reads_per_lane=output_reads/num_lanes)
      
      data<-sequencers2 %>% filter(is.na(read_length) | read_length==input$read_length_choice) 
      data
    }
  })
  
  read_var<-reactive({
    read_var<-as.numeric(data() %>% select(length) %>% filter(!is.na(length)) %>% distinct())
    read_var
  })
  
  total_bp<-reactive({
    wgs_total_bp<-input$wgs_quantity*(input$genome_size*(10^9)*input$wgs_coverage)
    rna_total_bp<-input$rna_quantity*(input$rna_reads*(10^6)*read_var())
    wgbs_total_bp<-input$wgbs_quantity*(input$genome_size*(10^9)*input$wgbs_coverage)
    atac_total_bp<-input$atac_quantity*(input$atac_reads*(10^6)*read_var())
    ddrad_total_bp<-input$ddrad_quantity*(input$ddrad_reads*(10^6)*read_var())
    
    total_bp<-sum(wgs_total_bp,
                  rna_total_bp,
                  wgbs_total_bp,
                  atac_total_bp,
                  ddrad_total_bp)
    total_bp
  })
  
  total_Gbp<-reactive({
    total_Gbp<-total_bp()*(10^-9)
    total_Gbp
  })
  
  total_reads<-reactive({
    total_reads<-(total_bp()/read_var())*(10^-6)
    total_reads
  })
  
  data2<-reactive({
    data2<-data() %>% 
      mutate(Percent_filled = (total_Gbp()/output_Gbp)*100) %>% 
      mutate(Percent_filled_lane = (total_Gbp()/output_Gbp_per_lane)*100) %>%
      mutate(Percent_filled_cost = (cost*(Percent_filled/100)))
    data2
  })
  
  output$txtOutput1 = renderText({
    paste0(format(round(total_Gbp(),2),scientific = F, big.mark=","), " Gbp")
  })
  
  output$txtOutput2 = renderText({
    paste0(format(round(total_bp(),2),scientific = F, big.mark=","), " bp")
  })
  
  output$txtOutput3 = renderText({
    paste0(format(round(total_reads(),2),scientific = F, big.mark=","), " M reads")
  })
  
  output$dataTableOutput1 = DT::renderDataTable({
    data3<-dcast(setDT(data2()), sequencer~min_max, value.var = c(names(data2()))) %>% 
      select(sequencer, num_lanes_min, output_Gbp_min, output_Gbp_max, 
             output_Gbp_per_lane_min, output_Gbp_per_lane_max,
             cost_min, cost_per_lane_min, 
             Percent_filled_max,Percent_filled_min, Percent_filled_lane_max,Percent_filled_lane_min,
             Percent_filled_cost_max,Percent_filled_cost_min) %>% 
      mutate(overfilled=ifelse(Percent_filled_min>100, "yes","no")) %>%
      arrange(output_Gbp_max)
    
    datatable(data3, rownames = FALSE, container = sketch, 
              options = list(
                columnDefs = list(
                  list(targets = "_all", className = "dt-center")
                ),
                pageLength = 20
              )) %>%
      formatStyle(c(1,2,4,6,7,8,10,12,14,15), `border-right` = "solid 2px") %>%
      formatRound(c(9:14), digits=0) %>%
      formatCurrency(c(7,8,13,14), digits=0) %>%
      formatStyle(
        'overfilled',
        target = 'row',
        backgroundColor = styleEqual(c("yes"), c("gray"))
      )
  })
  
  output$pltOutput1 = renderPlot({
    total_Gbp<-round(total_Gbp(),2)
    total_reads<-round(total_reads(),2)
    
    A<-ggbarplot(data2(),"sequencer","output_Gbp",add="mean_range",fill="sequencer",
                 xlab="",ylab="Gbp",title="Flow Cell Output (Gbp)",subtitle="log10 scale",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=1)+
      yscale("log10")+coord_flip()+
      geom_hline(yintercept=total_Gbp, size=2)+
      geom_text(aes(10,total_Gbp,label = total_Gbp, hjust = 1.5))+
      theme_pubclean()+rremove("legend")
    
    C<-ggbarplot(data2(),"sequencer","output_reads",add="mean_range",fill="sequencer",
                 xlab="",ylab="Reads (millions)",title="Flow Cell Output (M reads)",subtitle="log10 scale",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
      yscale("log10")+coord_flip()+
      geom_hline(yintercept=total_reads, size=2)+
      geom_text(aes(10,total_reads,label = total_reads, hjust = 1.5))+
      theme_pubclean()+rremove("legend")

    A+C
  })
  
  output$pltOutput2 = renderPlot({
    total_Gbp<-round(total_Gbp(),2)
    total_reads<-round(total_reads(),2)
    
    B<-ggbarplot(data2(),"sequencer","output_Gbp_per_lane",add="mean_range",fill="sequencer",
                 xlab="",ylab="Gbp",title = "Lane Output (Gbp)",subtitle="log10 scale",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=1)+
      yscale("log10")+coord_flip()+
      geom_hline(yintercept=total_Gbp, size=2)+
      geom_text(aes(10,total_Gbp,label = total_Gbp, hjust = 1.5))+
      theme_pubclean()+rremove("legend")
    
    D<-ggbarplot(data2(),"sequencer","output_reads_per_lane",add="mean_range",fill="sequencer",
                 xlab="",ylab="Reads (millions)",title = "Lane Output (M reads)", subtitle="log10 scale",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
      yscale("log10")+coord_flip()+
      geom_hline(yintercept=total_reads, size=2)+
      geom_text(aes(10,total_reads,label = total_reads, hjust = 1.5))+
      theme_pubclean()+rremove("legend")
    
    B+D
  })
  
  output$pltOutput3 = renderPlot({
    E<-ggbarplot(data2(),"sequencer","Percent_filled",add="mean_range",fill="sequencer",
                 xlab="",ylab="Percent Filled",title = "Flow Cell Filled",subtitle=">100% is overfilled",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = -0.5,lab.nb.digits=2)+
      coord_flip(ylim=c(0,100))+theme_pubclean()+rremove("legend")
    
    F<-ggbarplot(data2(),"sequencer","Percent_filled_lane",add="mean_range",fill="sequencer",
                 xlab="",ylab="Percent Filled",title = "Lane Filled",subtitle=">100% is overfilled",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = -0.5, lab.nb.digits=2)+
      coord_flip(ylim=c(0,100))+theme_pubclean()+rremove("legend")
    
    G<-ggbarplot(data2(),"sequencer","Percent_filled_cost",add="mean_range",fill="sequencer",
                 xlab="",ylab="Cost (USD)",title = "Cost", subtitle="",
                 sort.val="desc",sort.by.groups = FALSE,
                 label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
      yscale("log10")+coord_flip()+
      theme_pubclean()+rremove("legend")
    
    E+F+G
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("sequencers.csv", sep = "")
    },
    content = function(file) {
      write.csv(sequencers, file, row.names = FALSE)
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)