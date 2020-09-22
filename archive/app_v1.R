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
    
    # Sidebar layout with input and output definitions 
    sidebarLayout(
        
        # Sidebar panel for inputs 
        sidebarPanel(width=3,
            img(src = "SequencingEstimator.png", height = 125, width = 120, style="float:right"),
            titlePanel("Sequencing Estimator"),
            
            helpText("This application is designed to help you choose the right sequencing platform and estimate costs."),
            helpText("The app currently provides space for WGS, RNA-seq, WGBS, ATAC-Seq, and ddRAD-Seq. 
                     However, given desired coverage or number of reads (in million), the fields could be used any type of sequencing."), 
            helpText("Most information was obtained from:",tags$a(href="https://www.illumina.com/systems/sequencing-platforms.html", "Illumina")),
            helpText("The recommended number of reads and coverage are default input for all options"),
            helpText("_______________________________________________________________"),
           
            helpText("If you have updated information on sequencing platforms or costs, \
                     you can download an example of the necessary data, update it, and re-upload it here for your needs."),
            downloadButton("downloadData", "Download"),
            helpText(""),
            fileInput("file1", "Choose CSV File",
                      multiple = TRUE,
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),

            # Input: Genome Size
            numericInput(inputId = "genome_size",
                         label = "Genome Size (Gb):",
                         value = 1.6),
            helpText("Not required for RNA-Seq, ATAC-Seq, ddRAD-Seq, or any sequencing type based on 
                     sequencing a specific number or reads rather than coverage."),
            helpText("_______________________________________________________________"),
            
            # Input: Read Length
            selectInput(inputId = "read_length_choice", label = "Read Length:", 
                        choices = levels(sequencers$read_length), selected="2 x 150"),
            helpText("Read length is not applicable for PacBio sequencing and will be included regardless."),
            helpText("Illumina sequencers will be filtered by which are applicable with the chosen read length."),
            helpText("_______________________________________________________________"),
            
            # Input: Whole Genome Sequencing (WGS)
            numericInput(inputId = "wgs_quantity",
                         label = "Number of samples for Whole Genome Sequencing (WGS):",
                         value = 0),
            numericInput(inputId = "wgs_coverage",
                         label = "Coverage for Whole Genome Sequencing (WGS):",
                         value = 30),
            helpText("_______________________________________________________________"),
            
            # Input: Transcriptome Sequencing (RNA-Seq)
            numericInput(inputId = "rna_quantity",
                         label = "Number of samples for RNA-Seq:",
                         value = 0),
            numericInput(inputId = "rna_reads",
                         label = "Number of reads (in millions) per sample for RNA-Seq:",
                         value = 20),
            helpText("_______________________________________________________________"),
            
            # Input: Whole Genome Bisulfite Sequencing (WGBS)
            numericInput(inputId = "wgbs_quantity",
                        label = "Number of samples for Whole Genome Bisulfite Sequencing (WGBS):",
                        value = 0),
            numericInput(inputId = "wgbs_coverage",
                        label = "Coverage for Whole Genome Bisulfite Sequencing (WGBS):",
                        value = 20),
            helpText("_______________________________________________________________"),
            
            # Input: ATAC-Seq
            numericInput(inputId = "atac_quantity",
                        label = "Number of samples for ATAC-Seq:",
                        value = 0),
            numericInput(inputId = "atac_reads",
                        label = "Number of reads (in millions) per sample for ATAC-Seq:",
                        value = 60),
            helpText("_______________________________________________________________"),
            
            # Input: ddRAD-Seq
            numericInput(inputId = "ddrad_quantity",
                        label = "Number of samples for ddRAD-Seq:",
                        value = 0),
            numericInput(inputId = "ddrad_reads",
                        label = "Number of reads (in millions) per sample for ddRAD-Seq:",
                        value = 5),
        ),
        
        # Main panel for displaying outputs 
        mainPanel(
            h2("Total Amount of Data Required:"),
            h4(textOutput("txtOutput1")),
            h4(textOutput("txtOutput2")),
            h4(textOutput("txtOutput3")),
            h3(""),
            dataTableOutput("dataTableOutput1"),
            plotOutput("pltOutput1"),
            plotOutput("pltOutput2"),
            plotOutput("pltOutput3"),
            plotOutput("pltOutput4"),
        )
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
        paste0(format(total_Gbp(),scientific = F, big.mark=","), " Gbp")
    })
    
    output$txtOutput2 = renderText({
        paste0(format(total_bp(),scientific = F, big.mark=","), " bp")
    })
    
    output$txtOutput3 = renderText({
        paste0(format(total_reads(),scientific = F, big.mark=","), " M reads")
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
        A<-ggbarplot(data2(),"sequencer","output_Gbp",add="mean_range",fill="sequencer",
                     xlab="",ylab="Gbp",title="Gbp Output by Sequencer",subtitle="on a log10 scale",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=1)+
            yscale("log10")+coord_flip()+
            geom_hline(yintercept=total_Gbp(), size=2)+
            geom_text(aes(10,total_Gbp(),label = total_Gbp(), hjust = 1.5))+
            theme_pubclean()+rremove("legend")
        
        B<-ggbarplot(data2(),"sequencer","output_Gbp_per_lane",add="mean_range",fill="sequencer",
                     xlab="",ylab="Gbp",title = "Gbp Output per Lane",subtitle="on a log10 scale",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=1)+
            yscale("log10")+coord_flip()+
            geom_hline(yintercept=total_Gbp(), size=2)+
            geom_text(aes(10,total_Gbp(),label = total_Gbp(), hjust = 1.5))+
            theme_pubclean()+rremove("legend")
        
        A+B
    })
     
    output$pltOutput2 = renderPlot({
        C<-ggbarplot(data2(),"sequencer","output_reads",add="mean_range",fill="sequencer",
                     xlab="",ylab="Reads (millions)",title="Millions of Reads Output by Sequencer",subtitle="on a log10 scale",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
            yscale("log10")+coord_flip()+
            geom_hline(yintercept=total_reads(), size=2)+
            geom_text(aes(10,total_reads(),label = total_reads(), hjust = 1.5))+
            theme_pubclean()+rremove("legend")
        
        D<-ggbarplot(data2(),"sequencer","output_reads_per_lane",add="mean_range",fill="sequencer",
                     xlab="",ylab="Reads (millions)",title = "Millions of Reads Output per Lane", subtitle="on a log10 scale",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
            yscale("log10")+coord_flip()+
            geom_hline(yintercept=total_reads(), size=2)+
            geom_text(aes(10,total_reads(),label = total_reads(), hjust = 1.5))+
            theme_pubclean()+rremove("legend")
        
        C+D
    })

    output$pltOutput3 = renderPlot({
        E<-ggbarplot(data2(),"sequencer","Percent_filled",add="mean_range",fill="sequencer",
                     xlab="",ylab="Percent Filled",title = "Amount of the Sequencer Filled",subtitle=">100% is overfilled",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = -0.5,lab.nb.digits=2)+
            coord_flip(ylim=c(0,100))+theme_pubclean()+rremove("legend")
        
        F<-ggbarplot(data2(),"sequencer","Percent_filled_lane",add="mean_range",fill="sequencer",
                     xlab="",ylab="Percent Filled",title = "Amount of the Sequencer Filled per Lane",subtitle=">100% is overfilled",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = -0.5, lab.nb.digits=2)+
            coord_flip(ylim=c(0,100))+theme_pubclean()+rremove("legend")
        
        E+F
    })

    output$pltOutput4 = renderPlot({
        G<-ggbarplot(data2(),"sequencer","cost",add="mean_range",fill="sequencer",
                     xlab="",ylab="Cost (USD)",title = "Total Cost per Flow Cell", subtitle="static graphic representation of flow cell cost",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
            yscale("log10")+coord_flip()+
            theme_pubclean()+rremove("legend")
        
        H<-ggbarplot(data2(),"sequencer","Percent_filled_cost",add="mean_range",fill="sequencer",
                     xlab="",ylab="Cost (USD)",title = "Cost given Percent Filled", subtitle="",
                     sort.val="desc",sort.by.groups = FALSE,
                     label=T,lab.hjust = 1.5, lab.vjust = 0.5,lab.nb.digits=2)+
            yscale("log10")+coord_flip()+
            theme_pubclean()+rremove("legend")
        
        G+H
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