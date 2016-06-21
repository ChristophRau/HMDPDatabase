# A Graphical User Interface for querying a genetics SQL Database using Shiny in R
# Version: 0.8
# Last Modified: 5/15/16
# 
# The following is an implementation of a GUI using the Shiny package in Rstudio.  Shiny programs have two scripts associated with them.  This script, ui.R, controls the appearance of
# The GUI and provides inputs to and displays outputs from Server.R which contains the actual functions.

#for details on how this page's layout works, please see http://shiny.rstudio.com/tutorial/ and http://shiny.rstudio.com/reference/shiny/latest/
shinyUI(fluidPage(
  titlePanel("Welcome to the HMDP Database Shiny Server v0.7"),
  tabsetPanel(
    tabPanel("Start Here/Login",
             h6("Welcome to the first iteration of the searchable HMDP Database.  Please click a relevant tab to begin."),
             textInput("Password","Enter password for full access"),
             actionButton("Password_Go","Login"),
             textOutput("PassOK")),
    tabPanel("Visualize GWAS Result",
             sidebarLayout(
               sidebarPanel(
      selectInput("DataViz_DataType",label=h3("Select a type of data"),choices = c("Clinical","Expression","Metabolite","Protein")), #a selectInput is a dropdown menu
			htmlOutput("DataViz_StudyUI"), #an 'htmlOutput' is actually a way to create dynamic inputs.  In this case, Server.R is taking the selection from above and creating a new
			                               #selectInput populated with all the studies which have that type of data
      htmlOutput("DataViz_FinalTableSelectUI"), #then this one is allowing for fine tuning of the selection (typically selecting which gender of mice to examine)
      htmlOutput("DataViz_PhenotypeUI"), #and finally this one gives you a list of all possible phenotypes that can be used.
			actionButton("DataViz_Calculate","Create Manhattan Plot"),  #This is a button which, when clicked, tells Server.R to start calculating.
			
      selectInput("DataViz_Chromosome",label="Which Chromosome?",choices=c("All",c(1:19),"X"),selected="All"), #another select input
			numericInput("DataViz_Lower_Bound","Lower Bound (In MB)",1,min=0), #a numeric input which will take any number
			numericInput("DataViz_Upper_Bound","Upper Bound (In MB)",999,min=0)
                 ),
               mainPanel(downloadButton('DataViz_Download', 'Download These Results') , #creates a download button which takes a created file from Server.R
                         plotOutput('DataViz_Manhattan'), #creates a plot
                         h5("At distances of less than 10Mb, the UCSC Genome Browser will Appear Below."),
                         htmlOutput("DV_GenomeBrowser"))) #once again an htmlOutput, but in this case it really is an output, namely a visualization of the UCSC genome browser
    ),
    tabPanel("Create Beeswarm Plot",
             sidebarLayout(
               sidebarPanel(
                 selectInput("Beeswarm_DataType",label=h3("Select a type of data"),choices = c("Clinical","Expression","Metabolite","Protein")),
                 htmlOutput("Beeswarm_StudyUI"),
                 htmlOutput("Beeswarm_FinalTableSelectUI"),
                 htmlOutput("Beeswarm_PhenotypeUI"),
                 textInput("Beeswarm_rsID","Enter your SNP of choice",value=""), #will take any string as an input
                 actionButton("Beeswarm_Calculate","Create Plot") ),
               mainPanel(h6("Each metabolite and protein dataset is entered into the database differently.  Until these data are standardized, 
                            functionality of the Beeswarm plots are not guaranteed! "),
                 plotOutput('BS_Plot'))              
               )             
             ),
    tabPanel("Visualize Values Across Strains and Tissues",
             sidebarLayout(
              sidebarPanel(
                selectInput("VVAS_DataType",label=h3("Select a type of data"),choices = c("Phenotype","Gene")),
                htmlOutput("VVAS_StudyUI"),
                htmlOutput("VVAS_SelectExperimentsUI"),
                checkboxGroupInput("VVAS_SelectStrains",label="Strain Groups",choices=c("Inbred","AxB","BxA","BxD","BxH","CxB"),selected=c("Inbred","AxB","BxA","BxD","BxH","CxB")),
                actionButton("VVAS_Calculate","Create Plot") ),
              mainPanel(downloadButton('VVAS_Download', 'Download These Results'),
                        plotOutput('VVAS_Plot')
                        #,textOutput("TEST_Checkbox")
                        )              
             ) ),
    tabPanel("Nonsynnonymous SNPs",
             sidebarLayout(
               sidebarPanel(
                 textInput("NonSynnon_Gene","Enter Gene (SYMBOL FOR NOW)"),
                 actionButton("NonSynnon_Calculate","Run")
                 ),
               mainPanel(htmlOutput("NonSynnon_Result")))),
    tabPanel("cis-eQTLs",
             sidebarLayout(
               sidebarPanel(
                 htmlOutput("ciseQTL_PhenotypeUI"),
                 numericInput("ciseQTL_Window","Size of cis-eQTL window in MB",min=0,value=2),
                 actionButton("ciseQTL_Calculate","Create Table")
                 ),
               mainPanel(dataTableOutput('ciseQTL_Table'),
                         downloadButton('ciseQTL_Download', 'Download These Results')))
               ),
    tabPanel("Gene/Phenotype Correlations",sidebarLayout(
      sidebarPanel(
        textInput("FC_Input","Enter your gene or phenotype name"),
        selectInput("FC_Type","What sort of data is this?",c("clinical","expression","metabolite","protein")),
        checkboxGroupInput("FC_Cor_Types","What should we correlate against?",c("clinical","expression","metabolite","protein","microbiota")),
        htmlOutput("FC_SelectExperimentsUI"),
        numericInput("FC_threshold","P-value threshold",min=0,value=.0000042,max=1),
        actionButton("FC_Calculate","Create Table")
      ),
      mainPanel(h6("Don't use HighFatHypothalmusMale for now... data formated differently, will break your output.  Correlations to transripts will take some time to run.  Keep this in mind."),
                dataTableOutput("FC_Output"),
                downloadButton('FC_Download', 'Download These Results'))
    )),
    tabPanel("Overlapping Loci",
             sidebarLayout(
               sidebarPanel(
                 selectInput("Overlap_window_or_rsID","Please select to begin",c("window","rsID"))
                ,htmlOutput("Overlap_chr"),
                htmlOutput("Overlap_LB"),
                htmlOutput("Overlap_additional"),
                
                numericInput("Overlap_threshold","P-value threshold",min=0,value=.0000042,max=1),
                checkboxInput("Overlap_includeGenes","Include eQTLs?",value=FALSE),
                actionButton("Overlap_Calculate","Create Table")
               ),
               mainPanel(dataTableOutput('Overlap_Table'),
                         downloadButton('Overlap_Download', 'Download These Results'))
               )),
    tabPanel("Generate LD Plot",
             sidebarLayout(
               sidebarPanel(
                 selectInput("LD_window_or_rsID","Please select to begin",c("window","rsID"))
                 ,htmlOutput("LD_chr"),
                 htmlOutput("LD_LB"),
                 htmlOutput("LD_additional"),
                 htmlOutput("LD_MAFCutoff"),
                 actionButton("LD_Calculate","Calculate!")
               ),
               mainPanel(h6("heatmap is now zoomable and hover-overable.  Click and drag a box around the area you want to examine.  Click without dragging to zoom out.  WARNING:  Memory intensive."),
                         textOutput("LD_rsIDOut"),d3heatmapOutput("LD_windowOut"))
             )),
    tabPanel("Gene Name Conversions",
             sidebarLayout(
               sidebarPanel(
                 textInput("Lookup_One","Please enter a gene name or probesetID"),
                 fileInput("Lookup_Batch","Or upload a file for batch conversion"),
                 actionButton("Lookup_Button", "Convert!")
                 ),
               mainPanel(
                 dataTableOutput("Lookup_Table")
                 )
               )),
    tabPanel("More Tools To Come!",h3("Soon...")),
    tabPanel("Bugs/Suggestions",
             textInput("Suggestion_Name","Name"),
             textInput("Suggestion_Report","Suggestion/Bug"),
             tags$style(type='text/css', "#Suggestion_Report { height: 300px; width: 600px; }"),
             actionButton("Suggestion_Button","Suggest!"),
             textOutput("Suggestion_Text"),
             h3("Planned Changes:"),
             h4("Make it Faster (Especially correlations)"),
             h4("Eliminate Bugs"),
             h4("Make it Look Nice")
             )
    )))